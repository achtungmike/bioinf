
# coding: utf-8

# Using 3 protein families (globins, C2H2 zinc fingers, cytochromes P450) build a training set consisting of 10 representative proteins with known structures in each family. 
# 
# 1. Using cross-validation, train and validate 
#     a k-NN predictor 
#         that takes amino acid sequence as input 
#             and assigns one of the 3 secondary structure states to each amino acid 
#                 residue: helix, strand, coil. 
#                     Use a sliding window (of some length) 
#                         to define feature vectors for each amino acid residues, 
#                             compare a simple binary vector 
#                                 to Blosum62 feature vector

# In[68]:


import pandas as pd
import numpy as np
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import KFold
from sklearn.model_selection import cross_val_score
import matplotlib.pyplot as plt
from collections import defaultdict


# In[69]:


# read secondary struct file
sec = open("data/ss.txt","r").readlines()


# In[70]:


# blosum62 dict
blos = {'A': [4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0],
        'R': [-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3],
        'N': [-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3],
        'D': [-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3, ],
        'C': [0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, ],
        'Q': [-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2, ],
        'E': [-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2, ],
        'G': [0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3, ],
        'H': [-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3, ],
        'I': [-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3, ],
        'L': [-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1, ],
        'K': [-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2, ],
        'M': [-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1, ],
        'F': [-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1, ],
        'P': [-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2, ],
        'S': [1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2, ],
        'T': [0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0, ],
        'W': [-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3, ],
        'Y': [-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1, ],
        'V': [0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4, ]}


# In[71]:


# Pull list of PDBs from Excel file
# Input:
#     sheet_name: Sheet name that the list of PDBs are in that you want
# Ouput:
#     ret: List of PDBs found in that Excel file on specified sheet
def get_pdb():
    pdbs = pd.read_excel("data/pdbs.xlsx",header=None)
    ret = list(pdbs[0])
    return(ret)


# In[72]:


# Pull the full seq or secondary struct from secondary structure file
# Input:
#     ind: index within ss_file for seq or sec struct you want to pull
#     ss_file: a ss.txt from RSCB PDB read into python via readlines().
# Ouput:
#     txt: Full sequence or secondary structure found at index given in ss.txt
def pull_full_text(ind, ss_file):
    #ind is the index in the file where we start from

    txt = ""

    # skip the first line which contains the seq marker >
    ind += 1
    # Get current line and strip off new line
    tmp = ss_file[ind]
    tmp = tmp.rstrip('\n')
    # Until we hit another Seq tag (>)
    while ">" not in tmp:
        # We ran off the end of the file, break
        if ind == len(ss_file) - 1:
            break
        # append string here to txt removing the new line
        txt += tmp

        # Update index and tmp for next loop
        ind += 1
        tmp = ss_file[ind]
        tmp = tmp.rstrip('\n')
    return(txt)


# In[73]:


# Translate secondary structure into more simple form use in this project of
# Coil, Helix, or Extended (beta sheet)
# Input:
#     ss_str: secondary structure str of more complex form (used by RSCB PDB)
# Ouput:
#     ret: simplified secondary structure list with just c's, h's, and e's.
def trans_ss(ss_str):
    ret = []
    for x in range(len(ss_str)):
        if ss_str[x] in [' ', 'S','T']:
            ret.append('c')
        elif ss_str[x] in ['G', 'I', 'H']:
            ret.append('h')
        elif ss_str[x] in ['B', 'E']:
            ret.append('e')

    return(ret)


# In[74]:


# Get index of our sequence and our secondary structure for this PDB
# Input:
#     pdb: PDB identifier with chain in the form of 2AGV:A
#     sec: a ss.txt from RSCB PDB read into python via readlines().
# Ouput:
#     ret: index of sequnce id tag and secondary structure tag in the ss.txt file input.
def get_ind(pdb, sec):
    for x in range(len(sec)):
        # We have a match on sequence
        if(pdb + ":sequence" in sec[x]):
            seq_ind = x
        # We have a match secondary structure
        elif(pdb + ":secstr" in sec[x]):
            ss_ind = x

    return(seq_ind, ss_ind)


# In[75]:


# Helper function to actually generate the binary vector for a given residue
# Input:
#     res: amino acid residue identified by capital single letter code.
# Ouput:
#     ret: list matching aa defined in this function where all residues that
#          do not match input residue become a 0 and the one that does match 
#       becomes a 1.  Order is maintained.
def binarize(res):
    aa = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M',
          'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    ret = []
    for x in aa:
        if x == res:
            ret.append(1)
        else:
            ret.append(0)

    return(ret)


# In[76]:


# Create feature vector representation for sequences
# Input:
#     seq: string representing full AA sequence of interest
#     window: Sliding window of interest (should be an odd number).
#     method: Can use binary vector or blosum62 to represent residues
#             method = 0: binary vector representation
#             method = 1: blosum62 vector representation
# Ouput:
#     seq_vec: list of length(sec) with a representation of type == method
#               for each residue in input residue, order is maintained.
def vectorize(seq, window, method):
    # Stream (up or down stream) window // 2 (we want int division)
    stream = window // 2

    # Method of feature definition, currently supported binary
    # or blosum62
    # 0 = binary
    # 1 = blosum62


    # We need 2 arrays
    # curr_res_vec: temp one that will hold each residue as we
    #               work on it. [1 row, 20 AA * window columns]
    # seq_vec: we will append the curr_res_vec one onto after
    #               we are done proc it. [1 row, 20 AA * window columns]
    curr_res_vec = np.empty((0, 20 * window), int)
    seq_vec = np.empty((len(seq), 20 * window), int)

    for idx, val in enumerate(seq):
        # empty out the array
        curr_res_vec = np.empty((0, 20 * window), int)

        # Get list of other residues inside window of our current residue being considered
        winds = [idx]
        # Below
        winds += range(idx - stream, idx, 1)
        # Above
        winds += range(idx + 1, idx + stream + 1, 1)

        # We need to loop through all these guys in our window
        for x in winds:
            # check if x is outside of our sequence to pad the neighborhood if it is!
            if x >= len(seq) or x < 0:
                # Pad all 0's
                curr_res_vec = np.append(curr_res_vec, np.zeros(20))
            else:
                res = seq[x]

                if method == 0:
                    # create binary vector using aa vector
                    curr_res_vec = np.append(curr_res_vec, binarize(res))
                elif method == 1:
                    curr_res_vec = np.append(curr_res_vec, blos[res])
                else:
                    raise ValueError('method needs to be 0 for binary or 1 for blosum62, not %s as you entered!' % str(method))
        # Moving onto next residue so append our results for this one to
        # our seq_vec
        seq_vec[idx] =  curr_res_vec

    return(seq_vec)


# In[77]:


# SO the features have been generated I need to extract all those features and
# put them in one data structure while the corresponding secondary
# structure label in a
def prep_for_knn(dat):

    labs = list()
    vects = list()

    # Loop through each PDB in keep (ALL PDBS)
    for index, row in dat.iterrows():
        # Loop through each index in seq
        for key, val in enumerate(row['seq']):
            # Append the secondary structure (label) and vector to dat
            labs.append(row['sec'][key])
            vects.append(row['vect'][key])

    # Collapse the list
    vects = np.vstack(vects)
    
    return(vects, labs)


# In[85]:


# Setup KNN and CV testing grounds
# Need to test each metric and each number of neighbors
# Remove: chebyshev
def knn(vects, labs, fn):
    metrics = ['euclidean']#, 'manhattan', 'hamming']
    nebs = list(range(3, 15, 2))
    results = defaultdict(list)
    for metric in metrics:
        fname = fn + ".png"
        print("Testing metric: ", metric)
        for neb in nebs:
            print("Testing neb: ", neb)
            nn = KNeighborsClassifier(n_neighbors=neb, metric=metric)
            cv = KFold(n_splits=10, shuffle=True)
            scoring = cross_val_score(estimator = nn, X = vects, y = labs,
                                      scoring = 'accuracy', cv = cv)
            results[metric].append(scoring.mean())
        plt.plot(nebs,results[metric])
        plt.scatter(nebs,results[metric])
    plt.title(fn)
    plt.legend(metrics, loc=9, bbox_to_anchor=(0.5, -0.2), ncol=3)
    plt.xticks(list(range(3, 13, 2)))
    plt.xlabel("KNN K-Value")
    plt.ylabel("Mean Accuracy")
    plt.savefig(fname,bbox_inches='tight')
    plt.close()
        
    return(results)


# In[78]:


# Build all the PDBs in our internal representation.

pdb_df = []

pdb_list = get_pdb()
for pdb in pdb_list:
        inds = get_ind(pdb,sec)
        seq_ind = inds[0]
        ss_ind = inds[1]
        seq = pull_full_text(seq_ind, sec)
        ss = pull_full_text(ss_ind, sec)
        ss = trans_ss(ss)
        pdb_df.append([pdb,seq,ss])


# # Binary Feature Vector

# In[88]:


# Command and Control (main)
# Create DF to hold out work:
# First do Neuro:
#   Get the PDB codes from the Excel doc
#   Find these PDBs in ss.txt and get their indices
#   Pull the full sequence and full secondary struct
#   Translate the secondary structure
#   Store all of this in our DF
#   Vectorize our sequence

# Now ready to run our CV and k-NN

windows = [5,7,9,11,13,15,17,19, 21]
bin_results = []

for win in windows:
    bin_keep = pd.DataFrame(columns=['name','seq','sec', 'vect'])
    for ind, row in enumerate(pdb_df):
        pdb = row[0]
        seq = row[1]
        ss = row[2]
        vec = vectorize(seq, win, 0)
        bin_keep = bin_keep.append({'name': pdb, 'seq': seq,'sec': ss, 'vect': vec}, pdb)

    prep = prep_for_knn(bin_keep)
    vects = prep[0]
    labs = prep[1]

    fn = "bin_win" + str(win)
    tmp_res = knn(vects,labs, fn)
    print(win)
    print(tmp_res)


# # Blosum62 Feature Vectors

# In[91]:


# Command and Control (main)
# Create DF to hold out work:
# First do Neuro:
#   Get the PDB codes from the Excel doc
#   Find these PDBs in ss.txt and get their indices
#   Pull the full sequence and full secondary struct
#   Translate the secondary structure
#   Store all of this in our DF
#   Vectorize our sequence

# Now ready to run our CV and k-NN

windows = [5,7,9,11,13,15,17,19, 21]
blos_results = []


for win in windows:
    blos_keep = pd.DataFrame(columns=['name','seq','sec', 'vect'])
    for ind, row in enumerate(pdb_df):
        pdb = row[0]
        seq = row[1]
        ss = row[2]
        vec = vectorize(seq, win, 1)
        blos_keep = blos_keep.append({'name': pdb, 'seq': seq,'sec': ss, 'vect': vec}, pdb)

    prep = prep_for_knn(blos_keep)
    vects = prep[0]
    labs = prep[1]


    fn = "blos_win" + str(win)
    tmp_res = knn(vects,labs, fn)
    print(win)
    print(tmp_res)


