class HiddenMarkovModel:
    def __init__(self, start_probs, emit_probs, trans_probs):
        """Initialize the HMM with start probabilities, emission probabilities, and transition probabilities.
        Emission and transition probabilities are dictionaries of dictionaries. States and symbols are derived from the probabilities.
        """
        self.start_probs = start_probs
        self.emit_probs = emit_probs
        self.trans_probs = trans_probs
        self.states = list(emit_probs.keys())
        self.symbols = list(next(iter(emit_probs.values())).keys())

    def get_start_prob(self, state):
        """Retrieve the start probability for the given state."""
        return self.start_probs.get(state, 0.0)

    def get_emit_prob(self, state, symbol):
        """Retrieve the emission probability for the given state and symbol."""
        return self.emit_probs.get(state, {}).get(symbol, 0.0)

    def get_trans_prob(self, from_state, to_state):
        """Retrieve the transition probability from from_state to to_state."""
        return self.trans_probs.get(from_state, {}).get(to_state, 0.0)

    def set_start_prob(self, state, prob):
        """Update the start probability for the given state."""
        if state in self.states:
            self.start_probs[state] = prob

    def set_emit_prob(self, state, symbol, prob):
        """Update the emission probability for the given state and symbol."""
        if state in self.states and symbol in self.symbols:
            self.emit_probs[state][symbol] = prob

    def set_trans_prob(self, from_state, to_state, prob):
        """Update the transition probability from from_state to to_state."""
        if from_state in self.states and to_state in self.states:
            self.trans_probs[from_state][to_state] = prob

    def joint_probability(self, obs_seq, state_seq):
        """Compute the joint probability of the observation sequence and state path."""
        seq_length = len(obs_seq)
        if seq_length == 0:
            return None
        if seq_length != len(state_seq):
            print("Observation sequence and state path have different lengths!")
            return None
        joint = self.get_start_prob(state_seq[0]) * self.get_emit_prob(state_seq[0], obs_seq[0])
        for idx in range(1, seq_length):
            joint *= self.get_trans_prob(state_seq[idx - 1], state_seq[idx]) * self.get_emit_prob(state_seq[idx], obs_seq[idx])
        return joint

    def forward(self, obs_seq):
        """Compute the forward probabilities for the observation sequence."""
        seq_length = len(obs_seq)
        if seq_length == 0:
            return []
        fwd_probs = [{} for _ in range(seq_length)]
        for state in self.states:
            fwd_probs[0][state] = self.get_start_prob(state) * self.get_emit_prob(state, obs_seq[0])
        for t in range(1, seq_length):
            for state in self.states:
                total = 0.0
                for prev_state in self.states:
                    total += fwd_probs[t - 1][prev_state] * self.get_trans_prob(prev_state, state)
                fwd_probs[t][state] = total * self.get_emit_prob(state, obs_seq[t])
        return fwd_probs

    def backward(self, obs_seq):
        """Compute the backward probabilities for the observation sequence."""
        seq_length = len(obs_seq)
        if seq_length == 0:
            return []
        bwd_rev = []
        current_bwd = {state: 1.0 for state in self.states}
        for t in range(seq_length - 1, -1, -1):
            new_bwd = {}
            for state in self.states:
                total = 0.0
                for next_state in self.states:
                    total += current_bwd[next_state] * self.get_trans_prob(state, next_state) * self.get_emit_prob(next_state, obs_seq[t])
                new_bwd[state] = total
            bwd_rev.append(new_bwd)
            current_bwd = new_bwd
        return bwd_rev[::-1]

    def viterbi(self, obs_seq):
        """Apply the Viterbi algorithm to find the most likely state path for the observation sequence."""
        seq_length = len(obs_seq)
        if seq_length == 0:
            return []
        vit_probs = {}
        backpointers = {}
        paths = {}
        for state in self.states:
            vit_probs[state] = self.get_start_prob(state) * self.get_emit_prob(state, obs_seq[0])
            backpointers[state] = None
            paths[state] = [state]
        for t in range(1, seq_length):
            new_vit = {}
            new_back = {}
            new_paths = {}
            for state in self.states:
                candidates = []
                for prev in self.states:
                    cand_prob = vit_probs[prev] * self.get_trans_prob(prev, state)
                    candidates.append((cand_prob, prev))
                max_cand_prob, best_prev = max(candidates)
                new_vit[state] = max_cand_prob * self.get_emit_prob(state, obs_seq[t])
                new_back[state] = best_prev
                new_paths[state] = paths[best_prev] + [state]
            vit_probs = new_vit
            backpointers = new_back
            paths = new_paths
        best_final_prob = max(vit_probs.values())
        best_final_state = max(vit_probs, key=vit_probs.get)
        return best_final_prob, paths[best_final_state]

    def baum_welch(self, obs_seq):
        """Perform one iteration of the Baum-Welch algorithm to update model parameters based on the observation sequence."""
        print(obs_seq)
        seq_length = len(obs_seq)
        if seq_length == 0:
            return
        alpha = self.forward(obs_seq)
        beta = self.backward(obs_seq)
        # Expectation step
        gamma = [{} for _ in range(seq_length)]
        xi = [{} for _ in range(seq_length - 1)]
        for t in range(seq_length):
            sum_ab = 0.0
            for state in self.states:
                gamma[t][state] = alpha[t][state] * beta[t][state]
                sum_ab += gamma[t][state]
            if sum_ab > 0:
                for state in self.states:
                    gamma[t][state] /= sum_ab
            if t == 0:
                for state in self.states:
                    self.set_start_prob(state, gamma[t][state])
            if t == seq_length - 1:
                continue
            # Compute xi[t]
            sum_xi = 0.0
            for from_state in self.states:
                xi[t][from_state] = {}
                for to_state in self.states:
                    x_prob = (alpha[t][from_state] * self.get_trans_prob(from_state, to_state) *
                              self.get_emit_prob(to_state, obs_seq[t + 1]) * beta[t + 1][to_state])
                    xi[t][from_state][to_state] = x_prob
                    sum_xi += x_prob
            if sum_xi > 0:
                for from_state in self.states:
                    for to_state in self.states:
                        xi[t][from_state][to_state] /= sum_xi
        # Maximization step
        # Update emission probabilities
        for state in self.states:
            denom = sum(gamma[t][state] for t in range(seq_length))
            for symbol in self.symbols:
                num = sum(gamma[t][state] for t in range(seq_length) if obs_seq[t] == symbol)
                if denom > 0:
                    self.set_emit_prob(state, symbol, num / denom)
                else:
                    self.set_emit_prob(state, symbol, 0.0)
        # Update transition probabilities
        for from_state in self.states:
            denom = sum(gamma[t][from_state] for t in range(seq_length - 1))
            for to_state in self.states:
                num = sum(xi[t][from_state][to_state] for t in range(seq_length - 1))
                if denom > 0:
                    self.set_trans_prob(from_state, to_state, num / denom)
                else:
                    self.set_trans_prob(from_state, to_state, 0.0)


def test():
    start_probs = {"5": 0.8, "M": 0.15, "3": 0.05}
    emit_probs = {"5": {"A": 0.20, "C": 0.30, "G": 0.30, "T": 0.20},
                  "M": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                  "3": {"A": 0.35, "C": 0.15, "G": 0.15, "T": 0.35}}
    trans_probs = {"5": {"5": 0.8, "M": 0.2, "3": 0.0}, "M": {"5": 0.0, "M": 0.9, "3": 0.1},
                   "3": {"5": 0.0, "M": 0.0, "3": 1.0}}
    model = HiddenMarkovModel(start_probs, emit_probs, trans_probs)
    # Joint probability
    obs = "ATGCAATGCGCATGCTAAAA"
    states = "555555MMMMMM33333333"
    joint_val = model.joint_probability(obs, states)
    print("\nJoint probability of " + obs + " and " + states + " : " + str(joint_val))
    # Forward
    obs = "ATGTGTGCACGCACCGTGCGACGCGTCGCGGAAGCTGTTATA"
    alpha_vals = model.forward(obs)
    fwd_prob = sum(alpha_vals[-1].values())
    print("\nForward probability of " + obs + " : " + str(fwd_prob))
    # Backward
    print("Calculate backward probabilities")
    beta_vals = model.backward(obs)
    # Viterbi
    obs = "ACAATGCCGTCTCCGCGACGCCTTTAATTAT"
    (vit_prob, vit_path) = model.viterbi(obs)
    print("\nRunning Viterbi algorithm for seq " + obs)
    print("Optimal probability: " + str(vit_prob))
    print("Optimal statepath: " + " ".join(vit_path))
    # Baum-Welch learning
    print("\n\nProbabilities before learning")
    print("Emission probabilities")
    print(model.emit_probs)
    print("Transition probabilities")
    print(model.trans_probs)
    print("Optimizing probabilities with Baum-Welch: ")
    obs = "AGGGACGCTAAGCTCGCGCGAGCGACGCCATTATAGCGTAGCTTTTTAT"
    print("Input sequence: " + obs)
    model.baum_welch(obs)
    obs = "ATGTGGCGCGCGGAAGCTGTTATA"
    print("Input sequence: " + obs)
    model.baum_welch(obs)
    obs = "AATCGCGAGCGGCCCGCGAAGCTGTTTTTTAATA"
    print("Input sequence: " + obs)
    model.baum_welch(obs)
    obs = "ATGATGCGCTCGATGCTATCGCGCCGCGCGCGAGCGGCCCGCGAAGCTGTTTTAGTTAATAATGATATTGTA"
    print("Input sequence: " + obs)
    model.baum_welch(obs)
    print("Probabilities after learning")
    print("Emission probabilities")
    print(model.emit_probs)
    print("Transition probabilities")
    print(model.trans_probs)

test()
