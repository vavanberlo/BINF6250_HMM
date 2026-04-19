# Introduction
The project focuses on implementing the Baum–Welch algorithm to train a Hidden Markov Model (HMM) using observed DNA sequences. Building on previously developed Forward and Backward algorithm from Week 2, the goal was to construct a Expectation–Maximization (EM) pipeline that estimates initial, transition, and emission probabilities directly from the input data.

# Pseudocode
```
Baum Welch class:
  BW Initialization function:

    Create random number generator from seed
    Using Dirchlet distribution:
    Create random initial probabilities
    Create random transition probabilities
    Create random emission probabilities
    Add pseudocounts to each of the values
    Convert to log-space

    Return initial, transition and emission probabilities

(Tweaking the Forward-Backward algorithm slightly)
  Foward Matrix function:

    Initialize the matrix and fill the positions with -inf
    
    Looping through each observation and its index
    # Base case
    For each state:
    Initialize initial probabilities for each current state by multiplying initial probs by emission probs

    # Recursion
    For each next observation and index:
      Look at all possible previous states,
      Combine the probability of being in each previous state with the probability of transitioning to the current state
      Sum all these possibilities in log-space
      Multiply by emission probability for current symbol in log-space
      Store all the values in the matrix  

    Return the probability matrix

  Backward Matrix function:

    Reverse all observations
    Initialise the matrix and fill the positions with zeros

    Iterating through the reversed observations and their indices
    # Base case
    If i == 0
      Skip since it is already ste to 0

    # Recursion
    Else:
    For each state, look at all possible next states
    For each next state, combine the transition probability
    Sum over all next states using logaddexp.reduce along axis=1
    Add emission probability for the next symbol
    Combine the backward probability from next position
    Store the value in the current column of prob_matrix

  Flip matrix left to right to maintain the original orfer
  Return the probability matrix

 Sequence Probability Function:

    Take the last column of the forward matrix
    For all states, sum the probabilities using logaddexp.reduce
    Return the total log-probability of the observation sequence

 Forward-Backward Function:

    Initialize the posterior matrix withe zeros in all positions
    Call the forward matrix, backward matrix and sequence probability of the observation sequence

    For each position:
      For each state:
          Combine the forward and backward probability at that state and position
          Normalise by subtracting the sequence probability
          Store in posterior matrix
        
    Return the forward, backward and posterior matrices

 Expectation Function:

    For each sequence in observations:
      Call the forward, backward and posterior matrix
      Call the transition posterior matrix

      Accumulation:
      Use the posterior matrix to accumulate emission and initial state counts
      Use the transition posterior matrix to accumulate transition counts
      Use the forward matrix to accumulate total sequence probability
  
  Transition Posterior Matrix Function:

    Call the sequence probability
    Initialise the transition posterior matrix and fill the positions with -inf

    For each position n from 0 to T-2:
      Add forward probability at current state and current position
      Add transition probability from current state to next state
      Add emission probability of the next observed symbol from next state
      Add backward probability of next state at next position
      Normalize by subtracting the total sequence probability

      Collapse transition posterior into a 2D array

    Return the transition posterior matrix 

  Accumulate Emission Counts Function:
    For each position and symbol
      Obtain the symbol index from the nucleotide map
    For each state, add the posterior probability and emission count for that state and symbol

  Accumulate Transition Counts Function:
    For each pair of states
       Add the transition posterior value to the transition count for that pair

  Accumulate Initial Counts Function:
    For the first column of the posterior matrix
      Add the value to initial state count

  Accumulate Total Sequence Probability Function:
    Take the probability of the current sequence
    Add it to the running total probability
  
  Normalization Function:
  (in log space)
    Add all the initial counts to get the total
    Subtract the total count for each state's initial count to normalize

    Add all emission counts 
    Add pseudocount to emission counts
    Subtract total emission count from each value to normalize

    Add each row of transition counts 
    Add pseudocount to transition counts 
    Subtract total transition count from each value to normalize

    Return the normalized values

  Maximization Function:
    Call the normalization function

    Take the normalized initial probabilities and store them as the new initial probabilities
    Take the normalized emission probabilities and store them as the new emission probabilities
    Take the normalized transition probabilities and store them as the new transition probabilities

  Baum Welch Function:
     Set previous total sequence probability as None
      For i from 0 to number of sweeps
        Run expectation and maximization

        For every 10 sweeps
            Print the current and previous sequence probability to check for convergence

            If the difference is smaller
              Quit training

            If not, set previous sequence probability as total sequence probability

        Reset all accumulated counts and total sequence probability to -inf for the next iteration  
      If convergence is never reached after the maximum number of sweeps, report

  Model Probability Function:
    Initialize a dictionary for initial probabilities
    For each state at i
      Take the log value and convert to raw probability using exp
      Print

    Initialize a dictionary for transition probabilities
    Loop through every state with its index
      For each state create a nested empty dictionary
      Take the log value and convert to raw probability using exp
      Print

    Initialize a dictionary for emission probabilities
    Reverse nucleotide map so that the index maps to symbol
    For each state
      Create a nested empty dictionary
      For each symbol
      Take the log value and convert to raw probability using exp
      Print

  Return the dictionaries  

```

# Successes
We developed a full-implementation of Baum-Welch algorithm while deepening our understanding of its underlying mechanics, particularly the expectation and normalization step. Along the way, we gained experience with Numpy vectorization, made thoughtful design choices regarding data structures, and were able to incorporate peer's feedback in our design to improve our overall implementation. We were successful in translating complex mathematical concepts into working code. After finishing the main algorithm, we successfully validated our model using generated training data and comparing to BW estimates.

# Struggles
We spent a lot of time trying to understand the algorithm. Once we finally understood it, we had to refactor some of the old functions to work in the new class. We wanted to go back and update some of the older functions to work more properly but just didn't have time.

# Personal Reflections
## Group Leader
Victoria Van Berlo:
This project was definitely the hardest conceptually. We spent a lot of meetings just trying to figure out what to do. It seemed like every time I thought I understood something, at the next stage, half of my assumptions about the previous stage were incorrect. The time we spent trying to understand the algorithm did eventually help me to have a clear idea of flow and functions required, so that was at least productive time spent. Understanding when we should be using log space and how to get the vectorization working were real sticking points. When we ran the final test on known-probabilities sequence data and got BW estimates close to the original probabilities, it was a really great feeling. 

## Other member
Aaronie Jersha Jenyfred: 
I found this week's implementation to be the hardest of the HMM series. I took a very long time to get the conceptual understanding right. And even while pseudocoding, the boundaries between the E step, the accumulation functions, and the normalization kept blurring together. Vectorization of the forward and backward matrices was super hard to understand. Getting the log space and raw probability calculations right was a task. Also, small mistakes in indexing or normalization often led to incorrect probability distributions, which took time to fix. However, this week's pipeline led to an understanding about model training and performance, along with the influential parameters.  

Chantera Lazard:
This week's implementation was very challenging as we took a long time conceptualizing Baum-Welch (~4-5 meetings), specifically the Expectation step. Also as we made an earlier decision to work with numpy arrays, we chose to learn vectorization along the way so as to not overload our functions with for loops, but with not having used numypy previously, it was also hard to comprehend, specifically understanding broadcasting rules and using `newaxis` to reshape arrays for certain arithmetic operations. Nonetheless, we were able to comprehend the main goal of Baum-Welch, gaining valuable experience in probabalistic modeling. 


# Generative AI Appendix
We used Claude Sonnet 4.6 to ask conceptual questions about Baum-Welch, starting with questions we had concerning external resources we found such as Wikipedia articles, and to validate our own understanding.
