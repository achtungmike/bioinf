
# coding: utf-8

# ### Overview
# 
# This program will complete overlap alignment on 2 sequences using a dynamic programming method.
# It will produce two different alignments:
# 1. The suffix gaps on S1 and prefix gaps on S2 are penalty free.
#     * Aligning the S1 [blue] suffix to the S2 [orange] prefix 
#     
#     <img src="sc_1.png">
# 2. The prefix gaps on S1 and the suffix gaps on S2 are penalty free.
#     * Aligning the S1 [blue] prefix to the S2 [orange] suffix 
#     
#     <img src="sc_2.png">
# 
# 
# This program will utilize the following scoring parameters:
#     1. Match Premium = 1
#     2. Mismatch Penalty = -1
#     3. Gap Penalty = -1
# 
# 
# It will be tested first using the following 2 sequences:
# S1 = ATCTGTTCTT
# S2 = CTTTATGTT
# 
# I will also test more broadly (these results will be included in the report).
# 
# Include in directory:
#     1. Input files
#     2. PDF Report:
#         a. Overall Description
#         b. DP table for S1 and S2 (pathway table)
#         c. Resulting alignment.

# In[28]:


# Packages Required

# Numpy will be used for our DP table.
# tabulate will be used to pretty print the DP tables.
import numpy as np
import tabulate


# In[29]:


# Asks user for two sequences to do overlap alignment on.

# Input:
#   None
# Output:
#   Two sequences represented as strings in a tuple.
# Notes:
#   No limit on the alphabet for the sequences is currently implemented.

def get_seqs():
    seq1 = input("Please enter your first sequence for alignment.\n")
    seq2 = input("Now, please enter your second seqeuence for alignment.\n")

    # For testing:
    #seq1 = "ATCTGTTCTT"
    #seq2 = "CTTTATGTT"
    
    return (seq1, seq2)


# In[30]:


# Generates the DP tables, including a scoring table and a direction table

# Input:
#   2 Sequences represented as strings
# Output:
#   Two numpy 2D arrays, one for scoring (tab) and another for direction (dir)
# Notes:
#  

def sg_build_tab(seq1, seq2):

    s1 = list(seq1)
    s1_len = len(seq1)
    s2 = list(seq2)
    s2_len = len(seq2)


    # S1 goes across top, S2 down side.

    # Create our DP table:
    tab = np.zeros((s2_len, s1_len), np.object)
    tab.fill("_")
    # Create a table to store directions
    dir = np.zeros((s2_len, s1_len), np.object)
    dir.fill("_")

    #  Fill our first row and col with 0's
    for x in range(0, s2_len):
        for y in range(0,s1_len):
            if (x == 0):
                tab[x,y] = 0
            if (y == 0):
                # Origin
                if x == 0:
                    tab[x,y] = 0
                # Below origin (take score from above and sub 1 for new gap)
                else:
                    tab[x,y] = tab[x-1,y] - 1


    # Loop over our 2D array, first select a row then go through each column
    # in that row. (skipping row 1 and column 1 which stay 0)
    for x in range(1, s2_len):
        for y in range(1, s1_len):

            # Look at diag
            diag = tab[x-1,y-1]
            # If we match
            if (s1[y] == s2[x]):
                diag += 1
            else:
                diag -= 1

            # (gap) Coming from horizontal
            hor = tab[x,y-1] -1
            vert = tab[x-1,y] -1

            # Set the score in the table
            tab[x,y] = max(diag, hor, vert)

            # Determine direction, default to horizontal gap
            if (tab[x,y] == hor):
                dir[x,y] = '>'
            elif (tab[x,y] == vert):
                dir[x,y] = 'V'
            else:
                dir[x,y] = '*'

    return (tab,dir)


# In[43]:


# Aligns 2 sequences given the scoring and direction DP arrays.

# Input:
#   2 Sequences represented as strings
#   A numpy 2D array holding the scoring DP table
#   A numpy 2D array holding the direction information.
# Output:
#   Two strings representing the aligned sequences.
# Notes:
#  

def sg_align(seq1, seq2, tab, dir):

    al_seq1 = ""
    al_seq2 = ""
    
    # First thing we need to do: find the max in the last col.
    col = tab[:,(tab.shape[1] -1)]
    x_max = max(col)
    x_max_ind = np.where(col == x_max)
    x_max_ind = x_max_ind[0][-1]


    
    # Need to loop backwards through x and y.
    x = tab.shape[0] - 1
    y = tab.shape[1] - 1



    # Starting from last column

    # Loop backwards till we hit first row or col
    while x > 0 and y > 0:

        #print("Seq1: ", al_seq1)
        #print("Seq2: ", al_seq2)

        # Prepend Gap on Sequence 1 (and line to seq 2)
        if x > x_max_ind:
            al_seq2 = seq2[x] + al_seq2
            al_seq1 = "-" + al_seq1
            # decrement x
            x -= 1
            continue

        # If we should move diag (no gaps)
        if dir[x, y] == '*':
            al_seq1 = seq1[y] + al_seq1
            al_seq2 = seq2[x] + al_seq2

            # decrement both x and y
            x -= 1
            y -= 1

        # Horizontal gap
        elif dir[x, y] == '>':
            al_seq1 = seq1[y] + al_seq1
            al_seq2 = '-' + al_seq2
            y -= 1

        # Vertical gap
        else:
            al_seq1 = '-' + al_seq1
            al_seq2 = seq2[x] + al_seq2
            x -= 1


    # Dump the rest of the chars into seq and prepend with gaps for other
    if x > 0:
        while x > 0:
            al_seq1 = '-' + al_seq1
            al_seq2 = seq2[x] + al_seq2
            x -= 1

    elif y > 0:
        while y > 0:
            al_seq1 = seq1[y] + al_seq1
            al_seq2 = '-' + al_seq2
            y -= 1

    return (al_seq1, al_seq2, x_max)



# In[32]:


# Helper function to print out dynamic programming tables.

# Input:
#   2 Sequences represented as strings
#   A numpy 2D array holding the scoring DP table
#   A numpy 2D array holding the direction information.
# Output:
#   None - outputs directly to console.
# Notes:
#  

def pretty_print(seq1, seq2, tab, dir):

    s1 = list(seq1)
    s2 = list(seq2)

    print(tabulate.tabulate(tab, headers=s1, showindex=s2, tablefmt="grid"))
    print("--------------------------------------")
    print(tabulate.tabulate(dir, headers=s1, showindex=s2, tablefmt="grid"))


# In[44]:


# Get our sequences from the user:
seqs = get_seqs()
seq1 = seqs[0]
seq2 = seqs[1]

# Prepend a hyphen to each as a place holder for first row/col
seq1 = "-" + seq1
seq2 = "-" + seq2


# Print out the unaligned sequences
print("=========================")
print("Sequence 1: ", seq1[1:])
print("Sequence 2: ", seq2[1:])
print("=========================")


# Build our scoring and direction arrays.
build = sg_build_tab(seq1, seq2)

tab = build[0]
dir = build[1]

# Pretty print the determined tables
#pretty_print(seq1, seq2, tab, dir)

# Align the sequences
al = sg_align(seq1, seq2, tab, dir)
align_seq1 = al[0]
align_seq2 = al[1]
score = al[2]

# Print aligned sequences
print("Alignment Score: ", score)
print("=========================")
print(align_seq1)
print(align_seq2)
print("=========================")

# Run this all again in the opposite direction:
build = sg_build_tab(seq2, seq1)

tab = build[0]
dir = build[1]

# Pretty print the determined tables
#pretty_print(seq2, seq1, tab, dir)

al = sg_align(seq2, seq1, tab, dir)
align_seq1 = al[1]
align_seq2 = al[0]
score = al[2]
print("Alignment Score: ", score)
print("=========================")
print(align_seq1)
print(align_seq2)
print("=========================")


# In[47]:


# Get our sequences from the user:
seqs = get_seqs()
seq1 = seqs[0]
seq2 = seqs[1]

# Prepend a hyphen to each as a place holder for first row/col
seq1 = "-" + seq1
seq2 = "-" + seq2


# Print out the unaligned sequences
print("=========================")
print("Sequence 1: ", seq1[1:])
print("Sequence 2: ", seq2[1:])
print("=========================")


# Build our scoring and direction arrays.
build = sg_build_tab(seq1, seq2)

tab = build[0]
dir = build[1]

# Pretty print the determined tables
#pretty_print(seq1, seq2, tab, dir)

# Align the sequences
al = sg_align(seq1, seq2, tab, dir)
align_seq1 = al[0]
align_seq2 = al[1]
score = al[2]

# Print aligned sequences
print("Alignment Score: ", score)
print("=========================")
print(align_seq1)
print(align_seq2)
print("=========================")

# Run this all again in the opposite direction:
build = sg_build_tab(seq2, seq1)

tab = build[0]
dir = build[1]

# Pretty print the determined tables
#pretty_print(seq2, seq1, tab, dir)

al = sg_align(seq2, seq1, tab, dir)
align_seq1 = al[1]
align_seq2 = al[0]
score = al[2]
print("Alignment Score: ", score)
print("=========================")
print(align_seq1)
print(align_seq2)
print("=========================")


# In[48]:


# Get our sequences from the user:
seqs = get_seqs()
seq1 = seqs[0]
seq2 = seqs[1]

# Prepend a hyphen to each as a place holder for first row/col
seq1 = "-" + seq1
seq2 = "-" + seq2


# Print out the unaligned sequences
print("=========================")
print("Sequence 1: ", seq1[1:])
print("Sequence 2: ", seq2[1:])
print("=========================")


# Build our scoring and direction arrays.
build = sg_build_tab(seq1, seq2)

tab = build[0]
dir = build[1]

# Pretty print the determined tables
#pretty_print(seq1, seq2, tab, dir)

# Align the sequences
al = sg_align(seq1, seq2, tab, dir)
align_seq1 = al[0]
align_seq2 = al[1]
score = al[2]

# Print aligned sequences
print("Alignment Score: ", score)
print("=========================")
print(align_seq1)
print(align_seq2)
print("=========================")

# Run this all again in the opposite direction:
build = sg_build_tab(seq2, seq1)

tab = build[0]
dir = build[1]

# Pretty print the determined tables
#pretty_print(seq2, seq1, tab, dir)

al = sg_align(seq2, seq1, tab, dir)
align_seq1 = al[1]
align_seq2 = al[0]
score = al[2]
print("Alignment Score: ", score)
print("=========================")
print(align_seq1)
print(align_seq2)
print("=========================")

