"""
Static User Equilibrium (UE) traffic assignment for the DLSim engine.

Solves the classic Beckmann formulation of the traffic assignment problem

    min  sum_a  integral_0^{x_a} t_a(w) dw
    s.t. path flows >= 0 reproduce the OD demand, link flows are the sum of
         path flows crossing the link (Wardrop's first principle at optimum)

over the routable graph built from GMNS node/link files, using BPR volume
delay functions t_a(v) = fftt_a * (1 + alpha * (v / cap_a)^beta).

Two solvers are provided:

  * solve_frank_wolfe        - link-based Frank-Wolfe: all-or-nothing towards
                               the shortest-path tree, exact line search by
                               bisection on the Beckmann derivative, relative
                               gap convergence. Path probabilities per OD are
                               tracked alongside the link flows so individual
                               simulation agents can be assigned UE-consistent
                               paths afterwards.
  * solve_gradient_projection - path-based gradient projection (Jayakrishnan):
                               keeps an explicit path set per OD and shifts
                               flow from costlier paths onto the shortest path
                               scaled by the inverse second derivative.

Both return an AssignmentResult; DLSim.py consumes it to seed the mesoscopic
simulation with equilibrium routes.
"""

import heapq
import math
from collections import defaultdict

INF = float('inf')


class AssignmentLink:
    """Minimal link view used by the solvers (decoupled from DLSim.Link)."""

    __slots__ = ('seq', 'from_no', 'to_no', 'fftt', 'capacity', 'alpha', 'beta')

    def __init__(self, seq, from_no, to_no, fftt, capacity, alpha=0.15, beta=4.0):
        self.seq = seq
        self.from_no = from_no
        self.to_no = to_no
        # free-flow travel time in minutes
        self.fftt = max(1e-6, float(fftt))
        # capacity in vehicles per hour over all lanes
        self.capacity = max(1e-4, float(capacity))
        self.alpha = float(alpha)
        self.beta = float(beta)


class AssignmentResult:

    def __init__(self, link_flows, link_times, od_path_flows, od_demand,
                 iteration_log, unreachable_od_list):
        # link_flows[a]: UE flow on link a (veh/h), link_times[a]: minutes
        self.link_flows = link_flows
        self.link_times = link_times
        # {(o_no, d_no): {path(tuple of link seqs): flow share in [0, 1]}}
        self.od_path_flows = od_path_flows
        # {(o_no, d_no): demand in veh/h}
        self.od_demand = od_demand
        # list of dict rows: iteration, step_size, objective, tstt, sptt, rel_gap
        self.iteration_log = iteration_log
        self.unreachable_od_list = unreachable_od_list

    @property
    def relative_gap(self):
        return self.iteration_log[-1]['rel_gap'] if self.iteration_log else INF


class UEAssignment:

    def __init__(self, num_nodes, links, od_demand):
        """
        num_nodes: node count, nodes addressed by internal sequence number
        links:     list of AssignmentLink
        od_demand: {(o_node_no, d_node_no): demand in veh/h}
        """
        self.num_nodes = num_nodes
        self.links = links
        self.num_links = len(links)
        self.outgoing = [[] for _ in range(num_nodes)]
        for link in links:
            self.outgoing[link.from_no].append(link)
        # group destinations by origin so one shortest-path tree serves them all
        self.od_by_origin = defaultdict(dict)
        for (o, d), demand in od_demand.items():
            if o != d and demand > 0:
                self.od_by_origin[o][d] = self.od_by_origin[o].get(d, 0.0) + demand
        self.od_demand = {(o, d): dem
                          for o, dests in self.od_by_origin.items()
                          for d, dem in dests.items()}

    # ------------------------------------------------------------------
    # BPR volume delay function and friends
    # ------------------------------------------------------------------
    def link_time(self, link, flow):
        return link.fftt * (1.0 + link.alpha * (flow / link.capacity) ** link.beta)

    def link_time_derivative(self, link, flow):
        if link.beta <= 0:
            return 0.0
        return (link.fftt * link.alpha * link.beta / link.capacity
                * (flow / link.capacity) ** (link.beta - 1.0))

    def link_times(self, flows):
        return [self.link_time(l, flows[l.seq]) for l in self.links]

    def beckmann_objective(self, flows):
        """sum_a integral_0^x t_a(w) dw, closed form for BPR."""
        total = 0.0
        for l in self.links:
            x = flows[l.seq]
            total += l.fftt * (x + l.alpha * l.capacity / (l.beta + 1.0)
                               * (x / l.capacity) ** (l.beta + 1.0))
        return total

    # ------------------------------------------------------------------
    # Shortest paths / all-or-nothing
    # ------------------------------------------------------------------
    def _shortest_path_tree(self, origin, times):
        """Dijkstra one-to-all. Returns (label costs, predecessor link seqs)."""
        dist = [INF] * self.num_nodes
        pred_link = [-1] * self.num_nodes
        dist[origin] = 0.0
        heap = [(0.0, origin)]
        while heap:
            d, u = heapq.heappop(heap)
            if d > dist[u]:
                continue
            for link in self.outgoing[u]:
                nd = d + times[link.seq]
                if nd < dist[link.to_no]:
                    dist[link.to_no] = nd
                    pred_link[link.to_no] = link.seq
                    heapq.heappush(heap, (nd, link.to_no))
        return dist, pred_link

    def _backtrack(self, pred_link, destination):
        """Path to destination as a tuple of link seq nos (origin -> dest)."""
        path = []
        node = destination
        while pred_link[node] >= 0:
            link = self.links[pred_link[node]]
            path.append(link.seq)
            node = link.from_no
        path.reverse()
        return tuple(path)

    def all_or_nothing(self, times):
        """Assign every OD's full demand to its current shortest path.

        Returns (aon link flows, {od: shortest path}, unreachable od list).
        """
        flows = [0.0] * self.num_links
        od_paths = {}
        unreachable = []
        for origin, dests in self.od_by_origin.items():
            dist, pred_link = self._shortest_path_tree(origin, times)
            for dest, demand in dests.items():
                if dist[dest] == INF:
                    unreachable.append((origin, dest))
                    continue
                path = self._backtrack(pred_link, dest)
                od_paths[(origin, dest)] = path
                for seq in path:
                    flows[seq] += demand
        return flows, od_paths, unreachable

    # ------------------------------------------------------------------
    # Frank-Wolfe
    # ------------------------------------------------------------------
    def _line_search(self, x, y, tol=1e-10, max_bisections=50):
        """Optimal step towards the all-or-nothing solution.

        Minimizes the Beckmann objective along x + lam * (y - x) by bisecting
        on its derivative g(lam) = sum_a (y_a - x_a) * t_a(x + lam * (y - x)),
        which is nondecreasing in lam by convexity.
        """
        def g(lam):
            total = 0.0
            for l in self.links:
                seq = l.seq
                diff = y[seq] - x[seq]
                if diff != 0.0:
                    total += diff * self.link_time(l, x[seq] + lam * diff)
            return total

        if g(1.0) <= 0.0:
            return 1.0
        lo, hi = 0.0, 1.0
        for _ in range(max_bisections):
            mid = 0.5 * (lo + hi)
            if g(mid) <= 0.0:
                lo = mid
            else:
                hi = mid
            if hi - lo < tol:
                break
        return 0.5 * (lo + hi)

    def solve_frank_wolfe(self, max_iterations=100, rel_gap_tolerance=1e-4,
                          verbose=True):
        iteration_log = []
        times = self.link_times([0.0] * self.num_links)
        x, od_paths, unreachable = self.all_or_nothing(times)
        # path_flows[od][path] = share of the OD demand carried by the path;
        # the convex combination step contracts every share by (1 - lam) and
        # credits lam to the fresh all-or-nothing path, so shares always
        # reproduce the current link flow vector exactly.
        path_flows = {od: {path: 1.0} for od, path in od_paths.items()}

        for iteration in range(1, max_iterations + 1):
            times = self.link_times(x)
            y, aon_paths, _ = self.all_or_nothing(times)
            tstt = sum(times[a] * x[a] for a in range(self.num_links))
            sptt = sum(times[a] * y[a] for a in range(self.num_links))
            rel_gap = (tstt - sptt) / max(tstt, 1e-12)

            converged = rel_gap < rel_gap_tolerance
            step = 0.0
            if not converged:
                step = self._line_search(x, y)
                for a in range(self.num_links):
                    x[a] += step * (y[a] - x[a])
                for od, aon_path in aon_paths.items():
                    shares = path_flows.setdefault(od, {})
                    for path in list(shares):
                        shares[path] *= (1.0 - step)
                        if shares[path] < 1e-9:
                            del shares[path]
                    shares[aon_path] = shares.get(aon_path, 0.0) + step

            iteration_log.append({
                'iteration': iteration,
                'step_size': step,
                'objective': self.beckmann_objective(x),
                'tstt': tstt,
                'sptt': sptt,
                'rel_gap': rel_gap,
            })
            if verbose:
                print('FW iter %3d  rel_gap=%.6f  step=%.4f  obj=%.2f'
                      % (iteration, rel_gap, step, iteration_log[-1]['objective']))
            if converged:
                break

        self._normalize_path_flows(path_flows)
        return AssignmentResult(x, self.link_times(x), path_flows,
                                dict(self.od_demand), iteration_log,
                                unreachable)

    # ------------------------------------------------------------------
    # Gradient projection (path-based)
    # ------------------------------------------------------------------
    def solve_gradient_projection(self, max_iterations=50,
                                  rel_gap_tolerance=1e-4, step_scale=1.0,
                                  verbose=True):
        times = self.link_times([0.0] * self.num_links)
        _, od_paths, unreachable = self.all_or_nothing(times)
        # explicit path flows in demand units per OD
        path_flows = {}
        for od, path in od_paths.items():
            demand = self.od_by_origin[od[0]][od[1]]
            path_flows[od] = {path: demand}

        iteration_log = []
        x = self._link_flows_from_paths(path_flows)
        for iteration in range(1, max_iterations + 1):
            # snapshot costs drive the shortest-path trees and the gap
            # measure; the flow shifts below read the live x (Gauss-Seidel)
            # so each OD reacts to the shifts already made this iteration,
            # which keeps the projected Newton step stable under congestion
            times = self.link_times(x)
            tstt = sum(times[a] * x[a] for a in range(self.num_links))
            sptt = 0.0
            for origin, dests in self.od_by_origin.items():
                dist, pred_link = self._shortest_path_tree(origin, times)
                for dest, demand in dests.items():
                    if dist[dest] == INF:
                        continue
                    sptt += dist[dest] * demand
                    od = (origin, dest)
                    shortest = self._backtrack(pred_link, dest)
                    paths = path_flows[od]
                    paths.setdefault(shortest, 0.0)
                    sp_links = set(shortest)
                    for path in list(paths):
                        if path == shortest:
                            continue
                        sp_cost = sum(self.link_time(self.links[a], x[a])
                                      for a in shortest)
                        cost = sum(self.link_time(self.links[a], x[a])
                                   for a in path)
                        # second derivative over the symmetric difference of
                        # the two paths (shared links cancel in the shift)
                        second = sum(self.link_time_derivative(self.links[a], x[a])
                                     for a in path if a not in sp_links)
                        second += sum(self.link_time_derivative(self.links[a], x[a])
                                      for a in shortest if a not in path)
                        second = max(second, 1e-9)
                        shift = min(paths[path],
                                    step_scale * (cost - sp_cost) / second)
                        if shift > 0.0:
                            paths[path] -= shift
                            paths[shortest] += shift
                            for a in path:
                                x[a] -= shift
                            for a in shortest:
                                x[a] += shift
                        if paths[path] < 1e-9:
                            for a in path:
                                x[a] -= paths[path]
                            paths[shortest] += paths[path]
                            del paths[path]

            rel_gap = (tstt - sptt) / max(tstt, 1e-12)
            iteration_log.append({
                'iteration': iteration,
                'step_size': step_scale,
                'objective': self.beckmann_objective(x),
                'tstt': tstt,
                'sptt': sptt,
                'rel_gap': rel_gap,
            })
            if verbose:
                print('GP iter %3d  rel_gap=%.6f  obj=%.2f'
                      % (iteration, rel_gap, iteration_log[-1]['objective']))
            if rel_gap < rel_gap_tolerance:
                break

        x = self._link_flows_from_paths(path_flows)
        shares = {}
        for od, paths in path_flows.items():
            demand = self.od_by_origin[od[0]][od[1]]
            shares[od] = {p: f / demand for p, f in paths.items() if f > 0.0}
        self._normalize_path_flows(shares)
        return AssignmentResult(x, self.link_times(x), shares,
                                dict(self.od_demand), iteration_log,
                                unreachable)

    # ------------------------------------------------------------------
    # helpers
    # ------------------------------------------------------------------
    def _link_flows_from_paths(self, path_flows):
        flows = [0.0] * self.num_links
        for paths in path_flows.values():
            for path, flow in paths.items():
                for seq in path:
                    flows[seq] += flow
        return flows

    @staticmethod
    def _normalize_path_flows(path_flows):
        for paths in path_flows.values():
            total = sum(paths.values())
            if total > 0.0:
                for path in paths:
                    paths[path] /= total


def apportion_counts(shares, count):
    """Split `count` integer agents across paths by largest-remainder rounding
    so the realized path counts match the UE path shares as closely as an
    integer split allows. Returns {path: agent count}."""
    if not shares or count <= 0:
        return {}
    quotas = [(path, count * share) for path, share in shares.items()]
    counts = {path: int(math.floor(q)) for path, q in quotas}
    remainder = count - sum(counts.values())
    for path, q in sorted(quotas, key=lambda item: item[1] - math.floor(item[1]),
                          reverse=True)[:remainder]:
        counts[path] += 1
    return {path: n for path, n in counts.items() if n > 0}
