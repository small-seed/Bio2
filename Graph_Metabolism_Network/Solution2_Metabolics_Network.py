from solution_directed_graph import DirectedGraph  # Assuming the rewritten MyGraph is now DirectedGraph in solution_directed_graph.py

class MetabolicNetwork(DirectedGraph):
    def __init__(self, network_kind="metabolite-reaction", reverse_split=False):
        """Initialize the metabolic network with specified type and optional reversible split."""
        super().__init__(None)
        self.network_kind = network_kind
        self.node_categories = {}
        if network_kind == "metabolite-reaction":
            self.node_categories["metabolite"] = []
            self.node_categories["reaction"] = []
        self.reverse_split = reverse_split

    def insert_node_category(self, node, category):
        """Add a node with its category if not already present."""
        self.insert_node(node)
        if category not in self.node_categories:
            self.node_categories[category] = []
        if node not in self.node_categories[category]:
            self.node_categories[category].append(node)

    def nodes_by_category(self, category):
        """Retrieve nodes of a specific category."""
        return self.node_categories.get(category, [])

    def import_from_file(self, file_path):
        """Load the network from a file and convert based on network kind."""
        with open(file_path, 'r') as file_handle:
            temp_mr = MetabolicNetwork("metabolite-reaction")
            for line in file_handle:
                line = line.strip()
                if not line:
                    continue
                if ":" in line:
                    parts = line.split(":", 1)
                    reaction_id = parts[0].strip()
                    temp_mr.insert_node_category(reaction_id, "reaction")
                    reaction_side = parts[1].strip()
                else:
                    raise ValueError(f"Invalid line format: {line}")
                
                if "<=>" in reaction_side:
                    left_side, right_side = reaction_side.split("<=>")
                    left_mets = [met.strip() for met in left_side.split("+")]
                    for met in left_mets:
                        if met not in temp_mr.adjacency:
                            temp_mr.insert_node_category(met, "metabolite")
                        if self.reverse_split:
                            backward_reaction = reaction_id + "_b"
                            temp_mr.insert_node_category(backward_reaction, "reaction")
                            temp_mr.insert_edge(met, reaction_id)
                            temp_mr.insert_edge(backward_reaction, met)
                        else:
                            temp_mr.insert_edge(met, reaction_id)
                            temp_mr.insert_edge(reaction_id, met)
                    
                    right_mets = [met.strip() for met in right_side.split("+")]
                    for met in right_mets:
                        if met not in temp_mr.adjacency:
                            temp_mr.insert_node_category(met, "metabolite")
                        if self.reverse_split:
                            temp_mr.insert_edge(met, reaction_id + "_b")
                            temp_mr.insert_edge(reaction_id, met)
                        else:
                            temp_mr.insert_edge(met, reaction_id)
                            temp_mr.insert_edge(reaction_id, met)
                
                elif "=>" in reaction_side:
                    left_side, right_side = reaction_side.split("=>")
                    left_mets = [met.strip() for met in left_side.split("+")]
                    for met in left_mets:
                        if met not in temp_mr.adjacency:
                            temp_mr.insert_node_category(met, "metabolite")
                        temp_mr.insert_edge(met, reaction_id)
                    
                    right_mets = [met.strip() for met in right_side.split("+")]
                    for met in right_mets:
                        if met not in temp_mr.adjacency:
                            temp_mr.insert_node_category(met, "metabolite")
                        temp_mr.insert_edge(reaction_id, met)
                
                else:
                    raise ValueError(f"Invalid reaction format: {line}")
        
        if self.network_kind == "metabolite-reaction":
            self.adjacency = temp_mr.adjacency
            self.node_categories = temp_mr.node_categories
        elif self.network_kind == "metabolite-metabolite":
            self.transform_to_metabolite_network(temp_mr)
        elif self.network_kind == "reaction-reaction":
            self.transform_to_reaction_network(temp_mr)
        else:
            self.adjacency = {}

    def transform_to_metabolite_network(self, mr_network):
        """Convert metabolite-reaction network to metabolite-metabolite network."""
        for met in mr_network.node_categories["metabolite"]:
            self.insert_node(met)
            successors = mr_network.outgoing_neighbors(met)
            for succ in successors:
                succ_successors = mr_network.outgoing_neighbors(succ)
                for next_met in succ_successors:
                    if met != next_met and next_met not in self.outgoing_neighbors(met):
                        self.insert_edge(met, next_met)

    def transform_to_reaction_network(self, mr_network):
        """Convert metabolite-reaction network to reaction-reaction network."""
        for reac in mr_network.node_categories["reaction"]:
            self.insert_node(reac)
            successors = mr_network.outgoing_neighbors(reac)
            for succ in successors:
                succ_successors = mr_network.outgoing_neighbors(succ)
                for next_reac in succ_successors:
                    if reac != next_reac and next_reac not in self.outgoing_neighbors(reac):
                        self.insert_edge(reac, next_reac)

    # Metabolic potential assessment
    def enabled_reactions(self, available_mets):
        """Identify reactions that can be activated by available metabolites."""
        if self.network_kind != "metabolite-reaction" or not self.reverse_split:
            return None
        result = []
        for reac in self.node_categories['reaction']:
            predecessors = set(self.incoming_neighbors(reac))
            if predecessors and predecessors.issubset(set(available_mets)):
                result.append(reac)
        return result

    def generated_mets(self, enabled_reacs):
        """Find metabolites produced by enabled reactions."""
        result = []
        for reac in enabled_reacs:
            successors = self.outgoing_neighbors(reac)
            for succ in successors:
                if succ not in result:
                    result.append(succ)
        return result

    def complete_generated_mets(self, starting_mets):
        """Iteratively compute all producible metabolites from starting ones."""
        metabolites = list(starting_mets)
        changed = True
        while changed:
            changed = False
            reactions = self.enabled_reactions(metabolites)
            new_mets = self.generated_mets(reactions)
            for new_met in new_mets:
                if new_met not in metabolites:
                    metabolites.append(new_met)
                    changed = True
        return metabolites

    # Chapter 14 Exercise 1
    def terminal_mets(self):
        """Identify metabolites with no outgoing reactions."""
        result = []
        for node in self.nodes():
            if node.startswith("M") and self.outgoing_neighbors(node) == [] and self.incoming_neighbors(node):
                result.append(node)
        return result

    # Chapter 14 Exercise 2
    def minimal_path_to_product(self, starting_mets, target_met):
        """Find the shortest reaction path to produce the target metabolite."""
        if target_met in starting_mets:
            return []
        paths = {met: [] for met in starting_mets}
        reactions = self.enabled_reactions(starting_mets)
        changed = True
        while changed:
            changed = False
            new_reactions = []
            for reac in reactions:
                succs = self.outgoing_neighbors(reac)
                preds = self.incoming_neighbors(reac)
                for succ in succs:
                    if succ not in paths:
                        prev_paths = []
                        for pred in preds:
                            for path in paths.get(pred, []):
                                if path not in prev_paths:
                                    prev_paths.append(path)
                        new_path = prev_paths + [reac] if prev_paths else [reac]
                        paths[succ] = new_path
                        if succ == target_met:
                            return new_path
                        changed = True
                        for pred in preds:
                            if pred not in starting_mets:
                                new_reactions.append(reac)
            if changed:
                reactions = self.enabled_reactions(list(paths.keys()))
        return None


def validation_test1():
    """Test basic construction of a small metabolic network."""
    net = MetabolicNetwork("metabolite-reaction")
    net.insert_node_category("R1", "reaction")
    net.insert_node_category("R2", "reaction")
    net.insert_node_category("R3", "reaction")
    net.insert_node_category("M1", "metabolite")
    net.insert_node_category("M2", "metabolite")
    net.insert_node_category("M3", "metabolite")
    net.insert_node_category("M4", "metabolite")
    net.insert_node_category("M5", "metabolite")
    net.insert_node_category("M6", "metabolite")
    net.insert_edge("M1", "R1")
    net.insert_edge("M2", "R1")
    net.insert_edge("R1", "M3")
    net.insert_edge("R1", "M4")
    net.insert_edge("M4", "R2")
    net.insert_edge("M6", "R2")
    net.insert_edge("R2", "M3")
    net.insert_edge("M4", "R3")
    net.insert_edge("M5", "R3")
    net.insert_edge("R3", "M6")
    net.insert_edge("R3", "M4")
    net.insert_edge("R3", "M5")
    net.insert_edge("M6", "R3")
    net.display_adjacency()
    print("Reactions: ", net.nodes_by_category("reaction"))
    print("Metabolites: ", net.nodes_by_category("metabolite"))
    # ex 1
    print(net.terminal_mets())


def validation_test2():
    """Test loading and conversion of example network file."""
    print("metabolite-reaction network:")
    mr_net = MetabolicNetwork("metabolite-reaction")
    mr_net.import_from_file("example-net.txt")
    mr_net.display_adjacency()
    print("Reactions: ", mr_net.nodes_by_category("reaction"))
    print("Metabolites: ", mr_net.nodes_by_category("metabolite"))
    print()
    print("metabolite-metabolite network:")
    mm_net = MetabolicNetwork("metabolite-metabolite")
    mm_net.import_from_file("example-net.txt")
    mm_net.display_adjacency()
    print()
    print("reaction-reaction network:")
    rr_net = MetabolicNetwork("reaction-reaction")
    rr_net.import_from_file("example-net.txt")
    rr_net.display_adjacency()
    print()
    print("metabolite-reaction network (splitting reversible):")
    mr_split_net = MetabolicNetwork("metabolite-reaction", True)
    mr_split_net.import_from_file("example-net.txt")
    mr_split_net.display_adjacency()
    print()
    print("reaction-reaction network (splitting reversible):")
    rr_split_net = MetabolicNetwork("reaction-reaction", True)
    rr_split_net.import_from_file("example-net.txt")
    rr_split_net.display_adjacency()
    print()
    print(mm_net.average_degree("out"))
    print(mm_net.degree_distribution("out"))
    print(mm_net.average_path_lengths())
    print(mr_net.average_path_lengths())
    print(mm_net.all_clustering())
    print(mm_net.average_clustering())
    print(mm_net.clustering_by_degree())
    print(mm_net.top_degree_nodes(limit=3))
    print(mm_net.top_closeness_nodes(limit=3))
    print(mm_net.betweenness_measure("M5"))


def validation_test3():
    """Test analysis on E. coli network (requires ecoli.txt file)."""
    print("metabolite-reaction network:")
    ec_mr = MetabolicNetwork("metabolite-reaction")
    ec_mr.import_from_file("ecoli.txt")
    print("Size:", ec_mr.graph_size())
    print("Mean degree: ", ec_m
