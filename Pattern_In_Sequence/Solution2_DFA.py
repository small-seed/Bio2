class StringMatchingAutomaton:    
    def __init__(self, char_set, target_pattern):
        self.char_set = char_set
        self.target_pattern = target_pattern
        self.num_states = len(target_pattern) + 1
        self.transitions = {}
        self._construct_transitions(target_pattern)
    
    def _compute_prefix_overlap(self, candidate_str, full_pattern):
        max_overlap_len = min(len(candidate_str), len(full_pattern))
        for overlap_len in range(max_overlap_len, 0, -1):
            if candidate_str[-overlap_len:] == full_pattern[:overlap_len]:
                return overlap_len
        return 0
    
    def _construct_transitions(self, target_pattern):
        for state in range(self.num_states):
            for char in self.char_set:
                candidate = target_pattern[:state] + char
                next_state = self._compute_prefix_overlap(candidate, target_pattern)
                self.transitions[(state, char)] = next_state
    
    def display_structure(self):
        print(f"States: {self.num_states}")
        print(f"Alphabet: {self.char_set}")
        print("Transition table:")
        for (state, char), next_state in sorted(self.transitions.items()):
            print(f"{state}, {char} -> {next_state}")
    
    def get_transition(self, current_state, symbol):
        return self.transitions.get((current_state, symbol))
    
    def trace_states(self, input_sequence):
        states_visited = [0]
        current_state = 0
        for char in input_sequence:
            current_state = self.get_transition(current_state, char) or 0
            states_visited.append(current_state)
        return states_visited
    
    def locate_matches(self, input_text):
        matches = []
        current_state = 0
        pattern_length = len(self.target_pattern)
        for index, char in enumerate(input_text):
            current_state = self.get_transition(current_state, char) or 0
            if current_state == self.num_states - 1:
                start_pos = index - pattern_length + 1
                matches.append(start_pos)
        return matches

def run_demonstration():
    automaton = StringMatchingAutomaton("ACGT", "ACA")
    automaton.display_structure()
    sequence = "CACATGACATG"
    print(automaton.trace_states(sequence))
    print(automaton.locate_matches(sequence))

if __name__ == "__main__":
    run_demonstration()
