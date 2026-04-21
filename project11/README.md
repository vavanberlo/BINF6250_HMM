# Introduction
Profile HMM is a class of Hidden Markov Model designed to identify protein sequence similarities. A Profile HMM takes in a set of multiple sequence alignments and identifies regions of consensus and variability. Compared to a normal HMM model which has overall emission and transition probabilities for all positions of the sequence (except for the initial position), a Profile HMM has emission probabilities for each individual position. A Profile HMM is also directed and acyclical, meaning transitions can move only from left to right and have specific rules about which states may transition to other states. 


# Pseudocode

## Assign State to Columns
```
Input: List of Sequences
Return: List of State Types

# Initialization
* Empty List for State Type
* match_case = 0
* sequences_by_positions = Transpose sequences by positions using zip(*list)

# Iteration
Iterate through each position:
    initialize a gap count to 0
    for index in list of residues of position:
        Update count if position is gap ('-')
    Calculate percentage of gaps by dividing by number of sequences
    if gap_count <= 0.5:
        This is a match case
        append M{match_case+1} to state type list
        Add 1 to match_case
    else if gap_count > 0.5:
        This is an insertion case
        append I{match_case} to state type list

    return state_type_list, match_case
```

## Label Each Sequence
```
Input: List of sequences, state_types
Return: List of Lists where each list is the labels of the residue for each position
# Initialization
* A list to hold labeled sequences

# Iteration
for each sequence:
    match_case = 0
    sequence_label = Initialize "Begin" to hold each sequence by label
    for i, residue in enumerate(sequence):
        if residue is '-' and state_type[i] starts with 'M':
            # This represents a deletion
            Append to sequence_label D{match_case+1}
            match_case+= 1
        elif res is in alphabet and state_type[i] starts with 'M':
            # This is a match case
            Append to sequence_label M{match_case +1}
            match_case +=1
        elif res is in alphabet and state_type[i] starts with "I":
            # This is an insertion
            Append to sequence_label I{match_case}
    sequence_label.append("End")
    Append sequence_label to labeled sequences

```

## Estimate Initial Emissions Count
Rules:
* Match Cases: Position-specific distributions (Count residues at each position for each sequence in msa)
  * This is a consensus column
* Insertion Cases: Compute global aa freq from all training sequences and use background distribution for all $I_i$ states
  * This is extra residues in between consensus columns
* Deletion Cases: No emission
  * This is silence between consensus columns
* Pseudocounts are applied to avoid zeroes

```
Input
* sequences_by_positions
* state_type

* Return an emissions probability dictionary

# Initialization
* emit_probs = {}
* pseudocount = 0.01

# Iteration
for position in len(state_types):
    Initialize a counts dictionary with amino acids labels as keys and pseudocounts as values
    # Match Cases
    if state_type[pos].startswith("M"):
        Get position residues = sequences_by_positions[pos]
        for res in position residues:
            if res not in ['-', "Begin", "End"]:
                counts dictionary[res] +=1
        Calculate total across values in count dictionary
        emit_probs[state_type[pos]] = {amino acids: count dictionary[amino acids] / total for specific amino acid in alphabet}

    # Insertion case
    if state_type_pos.startswith("I"):
        Initialize a background counts dictionary with amino acids labels as keys and pseudocounts as values
        Iterate through each sequence
            for each residue in sequence
                if res not in ['-', "Begin", "End"]:
                    Add count to amino acid in background count dictionary
        Calculate total across values in background count dictionary
        emit_probs[state_type[pos]] = {amino acids: background count dictionary[amino acids] / total for specific amino acid in alphabet}
```

## Estimate Initial Transitions Count
**Pre Step: Build list of states (once we figure out state types (Step 1))**
```
*Input: Consensus length
*Output: List of states
states = ["Begin", "I0"]
for index, range(1, Consensus length + 1):
  states += [f"M{i}", f"I{i}", f"D{i}"]
states.append("End")
```

**Step A: Define Transition Rules**
```
* Input: State, Next State (Example: M1, M2 or M1, D2), Consensus Length
* Returns Boolean which confirms if current state can transition to next_state

if state == "Begin":
    return True if next_state is in ["M1", "I0", "D1"]
if state == "End":
    return False since pHMM are directed graphs flowing from left to right

n = Get the number for the state  (Example: If M3, n=3)

if state starts with "M" # Match
    if n == consensus length:
            return next_state in [I[[n], "End"]
    else return True if next_state is in [M[n+1], I[[n], D[n+1]]
if state starts with "I" # Insertion
    if n == consensus length:
        return next_state in [I[[n], "End"]
    else return True if next_state is in [I[n], M[n+1]]
if state starts with "D" # Deletion
    if n == length:
            return next_state in ["End"]
    return True if next_state is in [M[n+1], D[n+1]
    
Else return False
```

**Step B: Build Transition Skeleton**
```
Input: states, pseudocounts = 0.01
Return: Dictionary transition_counts skeleton with tuples defining valid transitions (prev_state, next_state) as keys and pseudocounts as values 

transition_counts = {}
for prev_state in states:
    for next_state in states:
        if Transition_Rules(prev_state, next_state) returns True as valid transition states
            transition_counts(state, next_state) = pseudocount
return  transition_counts   
```

**Step C: Compute Transition Counts**
```
Input 
* Consensus Length, labeled sequences (msa)
Return
* Dictionary transition_counts with tuples defining valid transitions (prev_state, next_state) as keys and actual counts as values
 
states = build_states(consensus length)
Initialize transition counts with transition skeleton
for seqs in labeled sequences:
  Index through labeled sequences:
    transition_counts[(seq[i], seq[i+1] += 1
```

**Normalize to Transition Probabilities**
```
Input
* labeled sequences (msa)
* states
Return
* Dictionary of dictionaries transition probs where outer key == prev_state and inner keys are outgoing states, values are transition probs of prev_state -> outgoing state

Get Trans counts dictionary inputting labeled sequences in function
trans_probs = {}

for state in states
  if state != 'End'
    outgoing = {k[1]: v for k, v in trans_counts.items() if k[0] == state}
  if outgoing: 
    total = Calculate the total across values for specific state
    for next_state, count in outgoing.items()
      trans_probs[state][next_state] = count / total 
```

## Viterbi
```
input: sequence, trans, emit, states
return: optimal path
matrix = {}
traceback = {}

# in log space
matrix["Begin"] = log(1)
# handle I0 start edge case
matrix["I0"] = matrix["Begin"] + log(transition["Begin"]["I0"]) + log(emission["I0"][sequence[0]])
traceback["I0"] = "Begin"

all other states matrix[state] = -inf

for each i, residue in enumerate(sequence):
    for each currstate in (Mi, Di, Ii-1):
        for each valid prevstate: # using transition rules
            if currstate is M or I:
                score = matrix[prevstate] + log(transition[prevstate][currstate] + log(emit[currstate][residue]))
            elif currstate is D: #no emission
                score = matrix[prevstate] + log(transition[prev][curr]
            if score > matrix[currstate]: #update when it's better, default -inf
                update matrix[currstate] = score
                update traceback[currstate] = prevstate

path = trace back thru traceback dict End to Begin
reverse path
return path
```

## Forward, Backward, Forward-Backward
```
forward():
input: sequence, trans, emit, states
return: probmatrix, totalprob

forward = {}
forward["Begin"] = log(1)
all other states forward[state] = -inf

for each currstate in sequence
    if M or I:
        forward[currstate] = logaddexp over all prevstate:
            forward[prevstate] + log(trans[prevstate][currstate] + log(emit[currstate][residue]))
    elif D: # no emission
        forward[currstate] = logaddexp over all prevstate:
            forward[prevstate] + log(trans[prevstate][currstate]

totalprob = logaddexp over all states:
    forward[currstate]
return forward, totalprob
        
backward():
input: sequence, trans, emit, states
return: probmatrix, totalprob
    
backward = {}
backward["End"] = log(1)
all other states backward[state] = -inf

for each currstate in sequence end to 1:
    if M or I:
        reverse[currstate] = logaddexp over all nextstate:
            reverse[nextstate] + log(trans[currstate][nextstate] + log(emit[currstate][residue]))
    elif D: # no emission
        reverse[currstate] = logaddexp over all nextstate:
            forward[currstate] + log(trans[currstate][nextstate]

totalprob = logaddexp over all states:
    reverse[currstate]

return reverse, totalprob

forward_backward():
input: seq, trans, emit, states
return: posterior, path
    
forward, totalf = forward()
backward, totalb = backward()
check if totalf and totalb sufficiently close
total_prob = totalf+totalb/2

posterior = {}

for each i in len(sequence):
    for each state (M, I, D):
        posterior[state][i] = forward[state][i] + backward[state][i] - total_prob # log space

# find best path
path = []
for each i in len(sequence):
    path.append argmax of M, I, D for posterior[state][i]

return posterior, path
```

# Successes
Our team met many times over the week and each time we moved forward a bit more in understanding. Each time we met we focused on a specific part of the algorithm which helped us focus in on specific areas we needed more work on.
Using the state diagram helped guide us in understanding possible transitions and training probabilities.

One strength we found was that we really understood the transition rules which were an integral part of the looping.

In our implementation, we chose to focus on a specific state diagram and implement that, rather than trying to reconcile several models into one. Because of that, we chose to include I_O and I_i+1 as possible states. 
[Reference diagram](https://www.ebi.ac.uk/training/online/courses/pfam-creating-protein-families/what-are-profile-hidden-markov-models-hmms/)


We also decided to focus more on the pseudocode to really understand the algorithm logic. We did implement some functions to help us test our pseudocode, but we focused on understanding. 

# Struggles
We struggled with our conceptual understanding. We think we found several typos in the resources which set us back time wise as we spent time trying to figure those out.
We struggled to nail down the relationship between the pHMM and the viterbi, forward, and backward algorithms which made it more difficult to implement changes in those functions.
We struggled to keep track of which functions/loops needed to use the match_case length and which needed the sequence length.
It was really hard to write pseudocode without it accidentally turning into code!

# Personal Reflections
## Group Leader
Victoria Van Berlo -
This was conceptually the hardest algorithm by far. The lack of clear guidelines for the project required a lot of research even outside the given resources, and I had to do a lot of back and forth. It seemed like some of the given resources conflicted with other resources' guidelines for how pHMMs work and so that was a considerable amount of confusion.

## Other member
Aaronie Jersha Jenyfred: The algorithm was conceptually hard to understand. I couldn't comprehend the transition from HMM to profile HMM. I also was uncertain about how viterbi and posterior decoding translates to the states assigned for profile HMM. Focusing on the pseudocode helped me to get the concepts right. A major point of confusion for me was separating model position vs sequence position, since in profile HMMs they don’t always move together. I also realized that small implementation details like indexing and state transitions can completely break the logic, which made me pay more attention to the structure of the algorithm.


Chantera Lazard-
I am grateful to my teammates for meeting up multiple times so we can reason through profileHMM. I think our goal was to write just pseudocode for this algorithm as sometimes it is more challenging to think through an algorithm compared to implementing it. I think this was one of our major successes out of the HMM suite where we just conceptualized algorithmically. We still wrote code as sometimes it is easier to write code before writing pseudocode, but we would write our ideas down in plain English or draw out ideas using *drawio* and then code so we can provide a more suitable pseudocode. profileHMM was easier for me to conceptualize because our resources did a good job in framing out the differences between a standard HMM vs a profileHMM. I did struggle with how viterbi and forward algorithms related to the profileHMM, but we reasoned that they were the algorithms used to test our trained data (estimated transition and emission probabilities). 

# Generative AI Appendix
Claude was consulted for a considerable amount of the conceptual understanding of profile HMM rules.
