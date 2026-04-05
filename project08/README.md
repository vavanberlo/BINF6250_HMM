# Introduction
Hidden Markov Models (HMMs) contain hidden states that we are trying to infer from observed data. This is useful in bioinformatics, because the observed data is what we can directly measure, like sequenced DNA. On the other hand, the hidden states represent the underlying biological context we are trying to uncover or infer. We will use the Viterbi algorithm, which utilizes dynamic programming to find the optimal paath of hidden stages given a sequence of observations. 

# Pseudocode
Put pseudocode in this box:

```

# Viterbi:

probability matrix, traceback matrix = len(states) * len(observation)

for each i, observation
    for each j, state
        if observation [0]:
            prob = current state initial prob * current emission prob 
            prob_matrix[j][i] = log(prob)
        else # calculate all probs possible for this box
            possible_probs =
                for all previous states: [ prev state prob * previous state to current state transition prob * current state emission prob ]
            max_prob = max(possible_probs)
            prob_matrix[j][i] = log(max_prob)

            #fill traceback
            prev_coords = argmax(possible_probs)
            traceback_matrix[j][i] = prev_coords

# Traceback:

end_state = argmax(prob_matrix[:,-1]
predictions = end_state
for i in len(observation) backwards:
    predictions.append (trace_matrix[end_state][i])
reverse predictions

for observation in predictions:
    states_final.append(states[observation])
    
```

# Successes
We had a solid game plan going in and did some solid peer programming. We kept in touch over the week and met several times and tried to make sure all our group members understood code as we went. We were successful in implementing a class-based structure depite limited prior experience in object-oritented programming. Having the members teach Viterbil algorithm to the group also proved effective, as it deeped our collective understanding and directly informed our implemenmtation strategy.

# Struggles
Although successfull in implementing the class-based structure, we did spend a lot of time figuring out how to implement it. We were also under the impression, the viterbi algorithm should be implemented outside the class, but late on Tuesday, we realized we needed to move the function inside the class object which required time to update rather last minute. Lastly, we received a tip about handling ties in ptential max probabilities. Although they may be rare, they can still happen. We did not have time to implement it during this algorithm, but as this is an ongoing project, we plan to add during our next meeting.

# Personal Reflections
## Victoria Van Berlo
This project required an approach I'm not used to. This project required us to keep the class object implementation in the back of our minds, so that the functions would suit the class setup and the class methods and variables would suit the functions. It was also an interesting challenge to create such a generalized implementation to suit any number of states and observations that were fed into it.

## Aaronie Jersha Jenyfred
Building a pipeline from scratch was new and made me understand the number of factors that are put into consideration. Though the algorithm was conceptually manageable, implementation of class structures and overall OOPs programming application was the tricky part for me. Debugging indexing errors and implementing log conversions(math domain errors) without disturbing the flow of the algorithm took a while. Developing a code that supports reproducibility and scalability was also a new experience for me.

## Chantera Lazard 
Initially, we were unsure where to start as we were building from scratch and that was overwhelming, but in re-visiting the lecture, we were able to find an entry point in planning out our implementation. Writing the algorithm seemed straightforward, but implementing a class exposed gaps in my understanding of OOP. I am looking forward to improving my knowledge and implementation as we continue to build on this project. 

# Generative AI Appendix
Generative AI was not used in this project.
