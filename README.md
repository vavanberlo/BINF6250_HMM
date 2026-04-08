# BINF6250_HMM
# Introduction
This project implements the Forward-Backward class of Hidden Markov Model. The Forward-Backward algorithm calculates the probability of each observed state at each position's probability of belonging to a hidden state. The Forward-Backward algorithm consists of 3 parts: the forward algorithm, which calculates probabilities from left to right, the backward algorithm, which calculates probabilities from right to left, and the combined forward-backward algorithm which uses the matrices created by the forward and backward algorithms to calculate the probability of a position in the observed sequence being assigned a particular hidden state. Our algorithm also calculates the probabilities for the entire matrix, ie, all possible observation position and hidden states.

# Pseudocode
Pseudocode for the major functions in our notebook

Note: This is the logic and does not reflect our exact implementation

**Foward_matrix**

```
Initialize numpy array num_states x len(observations)

# Base Case
for each state:
    Initialize initial probabilities for each current state by multiplying initial probs by emission probs
        Add initial probabilities in log space

# Recursion    
Iterate through each observation
    Iterate through each state    
    Calculate joint probabilities in log space
    Add joing probabilities in log space
    Add to probability matrix

Return probability matrix

```

**Backward_matrix**

```

Initialize numpy array num_states x len(observations)

# Base Case
for each state:
Initialize initial probabilities for each current state to log(1) = 0 (handled through numpy array init)

# Recursion    
Iterate through each observation starting from len(observations)-2 to 0 (backwards)
    Iterate through each state    
    Calculate joint probabilities in log space
    Add joint probabilities in log space
    Add to probability matrix
    
Return probability matrix 
        
```

**Forward-Backward**
```
Initialize numpy array num_states x len(observations)
Calculate the forward matrix
Calculate the backward matrix
Calculate the sequence probability


Iterate through each observation 
    Iterate through each state    
        Calculate the posterior probability by adding the backward and forward probabilities and subtracting 
        the sequence probability
        Store in posterior matrix
        
Return posterior matrix 
```

**Posterior-Decode**
```
Initialize empty list for path
Iterate through each oservation:
  Get the row index for the maximum posterior probability for each column
  Append to path

Return path
```


# Successes
Successful implementation of the Forward-Backward algorithm.
Our implementation from Week 1 really came in handy to extend the pipeline.
Maintained consistency in the class-based structuring. 
Successfully updated program to make calculations in log-space to avoid underflow.

# Struggles
Implementing the Backward algorithm was a bit tricky, especially in handling indexing correctly while reversing the observation sequence. 
Debugging mismatches in probabilities across positions took some time. 
Handling log-space computations and ensuring the correct use of logaddexp also took a while. 

# Personal Reflections
## Victoria Van Berlo
Overall, this algorithm was conceptually a bit more confusing for me than Viterbi, but balanced by requiring less coding from the ground up, since we re-used the class object framework. A few times I thought I had it, but then we had to go back and re-do something because it wasn't quite how we thought the first time around, which seems to be a running theme with me through many of these projects!

## Aaronie Jersha Jenyfred
I was a bit more comfortable this week since we had the class structure figured out by then. The algorithm was comparatively less complicated and implementing Viterbi first really eased the understanding. I also developed a more clearer and practical understanding of what expansion and reusability of class structures mean. Aligning the backward matrix and debugging indexing mismatches also required careful attention.

## Chantera Lazard 
The algorithm was easy to implement as building the probability matrix was similar to the same concept seen in Viterbi Matrix. However, understanding the difference between the forward and backward matrices and their necessity was pretty hard to wrap my head around, especially given Viterbi already provides state decoding. I am now understanding that it computes the posterior probability of being in each state at every position independently which will be useful for the last algorithm, Baum-Welch. 

# Generative AI Appendix
Generative AI was not used in this project.
