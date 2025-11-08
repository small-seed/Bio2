class DirectedGraph:
    def __init__(self, adj_dict=None):
        """Initialize the graph with an optional adjacency dictionary; defaults to empty."""
        self.adjacency = adj_dict if adj_dict is not None else {}

    def display_adjacency(self):
        """Display the graph's adjacency list."""
        for node in self.adjacency:
            print(f"{node} -> {self.adjacency[node]}")

    # Basic graph information
    def nodes(self):
        """Retrieve a list of all nodes in the graph."""
        return list(self.adjacency.keys())

    def edges(self):
        """Retrieve all edges as a list of (source, target) tuples."""
        edge_list = []
        for source in self.adjacency:
            for target in self.adjacency[source]:
                edge_list.append((source, target))
        return edge_list

    def graph_size(self):
        """Return the graph size: (number of nodes, number of edges)."""
        return len(self.nodes()), len(self.edges())

    # Node and edge addition
    def insert_node(self, node):
        """Insert a node if it doesn't already exist."""
        if node not in self.adjacency:
            self.adjacency[node] = []

    def insert_edge(self, source, target):
        """Insert an edge from source to target, adding nodes if necessary."""
        self.insert_node(source)
        self.insert_node(target)
        if target not in self.adjacency[source]:
            self.adjacency[source].append(target)

    # Neighbor relations
    def outgoing_neighbors(self, node):
        """Get the list of nodes reachable directly from the given node."""
        return self.adjacency[node].copy()

    def incoming_neighbors(self, node):
        """Get the list of nodes that point directly to the given node."""
        predecessors = []
        for key in self.adjacency:
            if node in self.adjacency[key]:
                predecessors.append(key)
        return predecessors

    def neighbors(self, node):
        """Get all adjacent nodes (incoming and outgoing, without duplicates)."""
        out = self.outgoing_neighbors(node)
        in_ = self.incoming_neighbors(node)
        all_neighbors = in_[:]
        for neigh in out:
            if neigh not in all_neighbors:
                all_neighbors.append(neigh)
        return all_neighbors

    # Degree calculations
    def out_degree(self, node):
        """Compute the out-degree of the node."""
        return len(self.adjacency[node])

    def in_degree(self, node):
        """Compute the in-degree of the node."""
        return len(self.incoming_neighbors(node))

    def total_degree(self, node):
        """Compute the total degree (in + out)."""
        return len(self.neighbors(node))

    def compute_all_degrees(self, degree_kind="total"):
        """Calculate degrees for all nodes: 'in', 'out', or 'total'."""
        degrees = {}
        for node in self.adjacency:
            if degree_kind == "out" or degree_kind == "total":
                degrees[node] = len(self.adjacency[node])
            else:
                degrees[node] = 0
        if degree_kind == "in" or degree_kind == "total":
            for node in self.adjacency:
                for neigh in self.adjacency[node]:
                    degrees[neigh] = degrees.get(neigh, 0) + 1
        return degrees

    def top_degree_nodes(self, degrees=None, kind="total", limit=10):
        """Return the top 'limit' nodes by degree."""
        if degrees is None:
            degrees = self.compute_all_degrees(kind)
        sorted_deg = sorted(degrees.items(), key=lambda item: item[1], reverse=True)
        return [item[0] for item in sorted_deg[:limit]]

    # Statistical metrics on degrees
    def average_degree(self, kind="total"):
        """Compute the average degree across all nodes."""
        degrees = self.compute_all_degrees(kind)
        return sum(degrees.values()) / len(degrees)

    def degree_distribution(self, kind="total"):
        """Compute the probability distribution of degrees."""
        degrees = self.compute_all_degrees(kind)
        dist = {}
        for node in degrees:
            deg = degrees[node]
            dist[deg] = dist.get(deg, 0) + 1
        total_nodes = len(degrees)
        for deg in dist:
            dist[deg] /= total_nodes
        return dist

    # Traversal methods
    def bfs_reachable(self, start):
        """Perform BFS to find all reachable nodes from start (excluding start)."""
        queue = [start]
        visited = []
        while queue:
            current = queue.pop(0)
            if current != start:
                visited.append(current)
            for neigh in self.adjacency[current]:
                if neigh not in visited and neigh not in queue:
                    queue.append(neigh)
        return visited

    def dfs_reachable(self, start):
        """Perform DFS to find all reachable nodes from start (excluding start)."""
        stack = [start]
        visited = []
        while stack:
            current = stack.pop(0)
            if current != start:
                visited.append(current)
            insert_pos = 0
            for neigh in self.adjacency[current]:
                if neigh not in visited and neigh not in stack:
                    stack.insert(insert_pos, neigh)
                    insert_pos += 1
        return visited

    def path_length(self, start, end):
        """Find the shortest path length from start to end using BFS."""
        if start == end:
            return 0
        queue = [(start, 0)]
        visited = {start}
        while queue:
            current, length = queue.pop(0)
            for neigh in self.adjacency[current]:
                if neigh == end:
                    return length + 1
                if neigh not in visited:
                    queue.append((neigh, length + 1))
                    visited.add(neigh)
        return None

    def shortest_route(self, start, end):
        """Find the shortest path from start to end as a list of nodes."""
        if start == end:
            return [start]
        queue = [(start, [])]
        visited = {start}
        while queue:
            current, path = queue.pop(0)
            for neigh in self.adjacency[current]:
                if neigh == end:
                    return path + [current, neigh]
                if neigh not in visited:
                    queue.append((neigh, path + [current]))
                    visited.add(neigh)
        return None

    def reachable_distances(self, start):
        """Return list of (node, distance) for all reachable nodes from start (excluding start)."""
        result = []
        queue = [(start, 0)]
        visited = []
        while queue:
            current, dist = queue.pop(0)
            if current != start:
                result.append((current, dist))
            for neigh in self.adjacency[current]:
                if not check_tuple_presence(queue + result, neigh):
                    queue.append((neigh, dist + 1))
        return result

    # Global metrics
    def average_path_lengths(self):
        """Compute average shortest path length and reachability ratio."""
        total_dist = 0
        reachable_count = 0
        node_count = len(self.nodes())
        for node in self.nodes():
            dists = self.reachable_distances(node)
            total_dist += sum(d for _, d in dists)
            reachable_count += len(dists)
        avg_dist = total_dist / reachable_count if reachable_count > 0 else 0
        reach_ratio = reachable_count / ((node_count - 1) * node_count) if node_count > 0 else 0
        return avg_dist, reach_ratio

    def closeness_measure(self, node):
        """Compute closeness centrality for the node."""
        dists = self.reachable_distances(node)
        if not dists:
            return 0.0
        total_dist = sum(d for _, d in dists)
        reachable = len(dists)
        return reachable / total_dist if total_dist > 0 else 0.0

    def top_closeness_nodes(self, limit=10):
        """Return top 'limit' nodes by closeness centrality."""
        closeness = {node: self.closeness_measure(node) for node in self.nodes()}
        sorted_close = sorted(closeness.items(), key=lambda item: item[1], reverse=True)
        return [item[0] for item in sorted_close[:limit]]

    def betweenness_measure(self, node):
        """Compute betweenness centrality for the node."""
        total_paths = 0
        paths_through = 0
        nodes = self.nodes()
        for s in nodes:
            for t in nodes:
                if s != t and s != node and t != node:
                    path = self.shortest_route(s, t)
                    if path is not None:
                        total_paths += 1
                        if node in path:
                            paths_through += 1
        return paths_through / total_paths if total_paths > 0 else 0.0

    # Cycle detection
    def has_self_cycle_from(self, node):
        """Check if there's a cycle starting and ending at the node."""
        queue = [node]
        visited = {node}
        while queue:
            current = queue.pop(0)
            for neigh in self.adjacency[current]:
                if neigh == node:
                    return True
                if neigh not in visited:
                    queue.append(neigh)
                    visited.add(neigh)
        return False

    def contains_cycle(self):
        """Check if the graph contains any cycle."""
        for node in self.nodes():
            if self.has_self_cycle_from(node):
                return True
        return False

    # Clustering analysis
    def local_clustering(self, node):
        """Compute the local clustering coefficient for the node."""
        adj_nodes = self.neighbors(node)
        adj_count = len(adj_nodes)
        if adj_count <= 1:
            return 0.0
        connections = 0
        for i in range(adj_count):
            for j in range(i + 1, adj_count):
                node_i, node_j = adj_nodes[i], adj_nodes[j]
                if node_j in self.adjacency[node_i] or node_i in self.adjacency[node_j]:
                    connections += 1
        possible = adj_count * (adj_count - 1) // 2
        return connections / possible if possible > 0 else 0.0

    def all_clustering(self):
        """Compute clustering coefficients for all nodes."""
        return {node: self.local_clustering(node) for node in self.nodes()}

    def average_clustering(self):
        """Compute the average clustering coefficient."""
        clusters = self.all_clustering()
        return sum(clusters.values()) / len(clusters) if clusters else 0.0

    def clustering_by_degree(self, kind="total"):
        """Compute average clustering coefficient per degree value."""
        degrees = self.compute_all_degrees(kind)
        clusters = self.all_clustering()
        deg_groups = {}
        for node, deg in degrees.items():
            deg_groups.setdefault(deg, []).append(node)
        avg_per_deg = {}
        for deg, nodes in deg_groups.items():
            avg_cluster = sum(clusters[n] for n in nodes) / len(nodes)
            avg_per_deg[deg] = avg_cluster
        return avg_per_deg


def check_tuple_presence(tuples_list, value):
    """Helper: Check if a value is the first element in any tuple in the list."""
    for t_val, _ in tuples_list:
        if t_val == value:
            return True
    return False


if __name__ == "__main__":
    g = DirectedGraph()
    g.insert_node(1)
    g.insert_node(2)
    g.insert_node(3)
    g.insert_node(4)
    g.insert_edge(1, 2)
    g.insert_edge(2, 3)
    g.insert_edge(3, 2)
    g.insert_edge(3, 4)
    g.insert_edge(4, 2)
    g.display_adjacency()
    print(g.graph_size())

    print(g.outgoing_neighbors(2))
    print(g.incoming_neighbors(2))
    print(g.neighbors(2))

    print(g.in_degree(2))
    print(g.out_degree(2))
    print(g.total_degree(2))

    print(g.compute_all_degrees("total"))
    print(g.compute_all_degrees("in"))
    print(g.compute_all_degrees("out"))

    g2 = DirectedGraph({1: [2, 3, 4], 2: [5, 6], 3: [6, 8], 4: [8], 5: [7], 6: [], 7: [], 8: []})
    print(g2.bfs_reachable(1))
    print(g2.dfs_reachable(1))

    print(g2.path_length(1, 7))
    print(g2.shortest_route(1, 7))
    print(g2.path_length(1, 8))
    print(g2.shortest_route(1, 8))
    print(g2.path_length(6, 1))
    print(g2.shortest_route(6, 1))

    print(g2.reachable_distances(1))

    print(g.contains_cycle())
    print(g2.contains_cycle())

    print(g.average_degree())
    print(g.degree_distribution())
    print(g.average_path_lengths())
    print(g.local_clustering(1))
    print(g.local_clustering(2))
