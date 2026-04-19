import numpy as np

class HiddenMarkovModel:
    """
    A Hidden Markov Model (HMM) with Viterbi algorithm implementation

    Attributes:
        states (list): The set of possible hidden states in the model.
        initial_probs (dict): Initial state probabilities mapping each state to its probability of being the starting state.
        transition_probs (dict): Transition probabilities mapping each state to a dictionary of probabilities for moving to another state.
        emission_probs (dict): Emission probabilities mapping each state to a dictionary of probabilities for observing each possible output.
    """
    def __init__(self, initial_probs, transition_probs, emission_probs, seed=None):
        self.states = []
        self.states = list(initial_probs.keys())
        self.initial_probs = initial_probs
        self.transition_probs = transition_probs
        self.emission_probs = emission_probs
        self.seed = seed
        self.nucleotide_map = {"A": 0, "C": 1, "G": 2, "T": 3}

    def get_transition_probs(self, state):
        """
        Returns the transition probabilities for a given state.
        :param state: The state whose transition probabilities should be returned.
        :return:
            dict: A dictionary mapping each state to its probability.
        """
        return self.transition_probs[state]

    def get_emission_probs(self, state):
        """
        Returns the emission probabilities for a given state.
        :param state: The state whose emission probabilities should be returned.
        :return:
            dict: A dictionary mapping each state to its probability.
        """
        return self.emission_probs[state]


    def viterbi_algorithm(self, observations):
        """
        Runs the Viterbi algorithm to find the optimal hidden state path for each observation.
        :param observations (list): A list of observations. If one observation sequence is provided, it is wrapped in a list
        :return:
            list: A list of optimal hidden state paths
        """
        if type(observations) != list:
            observations = [observations]

        #print(self.states)
        optimal_path = []
        for observation in observations:
            viterbi_matrix, traceback_matrix = self.build_viterbi_traceback_matrix(observation)
            #print(f"Observation: {observation}")
            #print(f"Viterbi matrix:\n{viterbi_matrix}")
            #print(f"Traceback matrix:\n{traceback_matrix}")
            optimal_path_index = self.viterbi_traceback(viterbi_matrix, traceback_matrix)
            #print(f"Optimal path index:\n{optimal_path_index}")
            path = self.convert_index_to_states(optimal_path_index)
            #print(f"Optimal path:\n{path}\n")
            optimal_path.append(path)


        return optimal_path


    def build_viterbi_traceback_matrix(self, observations):
        """
        Builds the Viterbi and traceback matrices for a given observation sequence.
        :param observations: A single observation
        :return
            tuple:
                - prob_matrix: Matrix of the most likely probabilities at each state and position.
                - traceback_matrix: Matrix storing the previous state indices for reconstructing the optimal path.
        """
        #Initialize viterbi and traceback matrix
        prob_matrix = np.zeros((len(self.states), len(observations)), dtype = float)
        traceback_matrix = np.zeros((len(self.states), len(observations)), dtype = int)

        ####### Iteration ########
        for i, observation in enumerate(observations):
            for j, state in enumerate(self.states):
                # Get initial and emission probabilities for current state
                state_init_probs = self.initial_probs[state]
                state_emit_probs = self.get_emission_probs(state)

                # If we are looking at the first observation
                if i == 0:
                    # Calculate the initial probability per state
                    state_prob = state_init_probs * state_emit_probs[observation] # Come back to you
                    # Use natural log to prevent numerical underflow
                    prob_matrix[j][i] = np.log(state_prob)

                # Else if we are looking at the second observation onward
                else:
                    # Calculate the possible probabilities based on transitioning for all states
                    possible_probs = [np.exp(prob_matrix[k][i-1]) * self.get_transition_probs(prev_state)[state] * state_emit_probs[observation] for k, prev_state in enumerate(self.states)]

                    # Get the max probability and add the log value to the viterbi matrix
                    max_prob = max(possible_probs)
                    prob_matrix[j][i] = np.log(max_prob)

                    # Get the index of the maximum probabilities and add that to the traceback matrix
                    # Theoretically, if there is a tie (although rare as we are dealing with floating point numbers), argmax would choose the first index making it reproducible
                    previous_coords = np.argmax(possible_probs)
                    traceback_matrix[j][i] = previous_coords

        return prob_matrix, traceback_matrix

    def viterbi_traceback(self, viterbi_matrix, traceback_matrix):
        """
        Traces back through the Viterbi and traceback matrices to recover the optimal hidden state path.

        :param viterbi_matrix: Matrix of the most likely probabilities at each state and position.
        traceback_matrix: Matrix storing the previous state indices for reconstructing the optimal path.

        :return:
            list: The optimal hidden state path as a sequence of indices
        """
        # Identify the final state with the highest probability
        end_state = np.argmax(viterbi_matrix[:,-1])

        # Add state to list
        predictions = [int(end_state)]

        # Traceback through the matrix starting at the end of the matrix
        for i in range(len(traceback_matrix[end_state])-1, 0, -1):
            # Append the state index to the predictions
            end_state = traceback_matrix[end_state][i]
            predictions.append(int(end_state))

        # Reverse the predictions so it starts at the beginning of matrix
        return predictions[::-1]

    def convert_index_to_states(self, predictions):
        """
        Converts a list of predicted indices into a list of state names.
        :param predictions (list): The list of indices that make up the optimal path
        :return:
            list: The list of state names in the optimal path
        """
        states_final = []
        for index in predictions:
            states_final.append(self.states[index])
        return states_final



class ForwardBackward(HiddenMarkovModel):
    def forward_matrix(self, observations):
        """
        The algorithm computes the forward probability matrix. Each entry represents the probability of having emitted the observed symbols up to the particular position and state.

        :param observations: list of observed sequences
        :return:
        prob_matrix: np.ndarray - Forward probability matrix of shape (num_states x len(observations)) in log space
        """
        prob_matrix = np.zeros((len(self.states), len(observations)), dtype=float)

        ####### Iteration ########
        for i, observation in enumerate(observations):
            for j, state in enumerate(self.states):
                # Get initial and emission probabilities for current state
                state_init_probs = self.initial_probs[state]
                state_emit_probs = self.get_emission_probs(state)

                # calculate probabilities
                # Remember log(x*y) = log(x) + log(y)
                # first column
                if i == 0:
                    state_prob = np.log(state_init_probs) + np.log(state_emit_probs[observation])
                    prob_matrix[j][i] = state_prob

                else:
                    joint_probs = [prob_matrix[k][i - 1] +
                                   np.log(self.get_transition_probs(prev_state)[state]) +
                                   np.log(state_emit_probs[observation])
                                   for k, prev_state in enumerate(self.states)]

                    sum_prob = np.logaddexp.reduce(joint_probs)
                    prob_matrix[j][i] = sum_prob

        return prob_matrix

    def backward_matrix(self, observations):
        """
        The algorithm computes the backward probability matrix. Each entry represents the probability of having emitted the observed symbols from the current positon to the end of the sequence given a specific state.

        :param observations: list of observed sequences
        :return:
        prob_matrix: np.ndarray - Forward probability matrix of shape (num_states x len(observations)) in log space
        """
        prob_matrix = np.zeros((len(self.states), len(observations)), dtype=float)

        for i in range(len(observations) - 2, -1, -1):
            for j, state in enumerate(self.states):
                joint_probs = [prob_matrix[k][i + 1] +
                               np.log(self.get_transition_probs(state)[next_state]) +
                               np.log(self.get_emission_probs(next_state)[observations[i + 1]])
                               for k, next_state in enumerate(self.states)]

                prob_matrix[j][i] = np.logaddexp.reduce(joint_probs)

        return prob_matrix

    def sequence_probability(self, prob_matrix):
        """
        Computes the total probability of the observation sequence from either the forward probability matrix.

        :param prob_matrix:  np.ndarray - Forward probability matrix of shape (num_states x len(observations))
        :return:
            float - Log probability of the observation sequence
        """
        return np.logaddexp.reduce(prob_matrix[:, -1])

    def forward_backward(self, observations):
        """
        The algorithm computes the forward-backward posterior matrix by combining the forward and backward matrices, normalized by the total probability of the observation sequence.

        :param observations: list of observed sequences
        :return:
        posterior_matrix: np.ndarray - posterior probability matrix of shape (num_states x len(observations)) in log space
        """
        self.posterior_matrix = np.zeros((len(self.states), len(observations)), dtype=float)

        forward_matrix = self.forward_matrix(observations)
        backward_matrix = self.backward_matrix(observations)

        final_col_prob = self.sequence_probability(forward_matrix)

        for i in range(len(observations)):
            for j in range(len(self.states)):
                self.posterior_matrix[j][i] = forward_matrix[j][i] + backward_matrix[j][i]

        return forward_matrix, backward_matrix, self.posterior_matrix

    def posterior_decode(self, observations):
        """
        Decodes the most probable state at each position independently using posterior decoding
        :param observations: list of observed sequences
        :return:
        path: list - The most probable state at each position
        """
        self.forward_backward(observations)
        path = []
        for i in range(len(observations)):
            best_state = np.argmax(self.posterior_matrix[:, i])
            path.append(self.states[best_state])

        return path

    def get_posterior_at_index(self, position):
        """
        Returns probabilities of a particular position (index) for all hidden states
        :param position: int index of desired position to check probabilities
        :return:
        list of probabilties

        """
        return self.posterior_matrix[:, position]
