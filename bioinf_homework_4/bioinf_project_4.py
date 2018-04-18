
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

# In[33]:





import pandas as pd
import numpy as np
from sklearn.neighbors import NearestNeighbors




# read secondary struct file
sec = open("data/ss.txt","r").readlines()

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


# In[34]:


# Pull list of PDBs from Excel file
# Input:
#     Sheet name that the list of PDBs are in that you want
# Ouput:
#     List of PDBs found in that Excel sheet
def get_pdb(sheet_name):
    pdbs = pd.read_excel("data/pdbs.xlsx",sheet_name=sheet_name,header=None)
    ret = list(pdbs[0])
    return(ret)


# In[55]:
#
#
# def find_pdb(pdb, sec_file):
#     seq_ind = pdb + ":sequence" in sec_file
#     ss_ind = pdb + ":secstr" in sec_file
#     return((seq_ind, ss_ind))


# In[35]:


# Pull the full seq or secondary struct from file
def pull_full_text(ind, ss_file):
    #ind is the index in the file where we start from

    txt = ""

    # skip the first line which contains the seq marker >
    ind += 1
    # Get current line and strip off new line
    tmp = ss_file[ind]
    tmp = tmp.strip('\n')
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
        tmp = tmp.strip('\n')
    return(txt)


# In[36]:


# Translate secondary structure into more simple form of
# Coil, Helix, or Extended (beta sheet)
def trans_ss(ss_str):
    ret = ""
    for x in range(len(ss_str)):
        if ss_str[x] in [' ', 'S','T']:
            ret += 'c'
        elif ss_str[x] in ['G', 'I', 'H']:
            ret += 'h'
        elif ss_str[x] in ['B']:
            ret += 'e'

    return(ret)


# In[38]:


#Get index of our sequence and our secondary structure for this PDB
def get_ind(pdb, sec):
    for x in range(len(sec)):
        # We have a match on sequence
        if(pdb + ":sequence" in sec[x]):
            seq_ind = x
        # We have a match secondary structure
        elif(pdb + ":secstr" in sec[x]):
            ss_ind = x

    return(seq_ind, ss_ind)




# Helper function to actually translate the residue into a binary vector
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


# Create binary vector representation for sequences
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
    seq_vec = np.empty((0, 20 * window), int)

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
        seq_vec = np.vstack((seq_vec, curr_res_vec))

    return(seq_vec)






def main():
    keep = pd.DataFrame(columns=['name','seq','sec', 'vect'])

    tmp = get_pdb("neuro")

    for pdb in tmp:
        inds = get_ind(pdb,sec)
        seq_ind = inds[0]
        ss_ind = inds[1]
        seq = pull_full_text(seq_ind, sec)
        ss = pull_full_text(ss_ind, sec)
        ss = trans_ss(ss)

        vec = vectorize(seq, 5, 1)
        keep = keep.append({'name': pdb, 'seq': seq,'sec': ss, 'vect': vec}, pdb)






def run_knn(dat, labs):
    from sklearn.model_selection import KFold
    from sklearn.model_selection import cross_val_score

    metrics = ['euclidean', 'manhattan', 'chebyshev', 'canberra', 'braycurtis']
    nebs = list(range(3, 15, 2))

    # Need to test each metric and each number of neighbors
    for metric in metrics:
        print("Testing metric: ", metric)
        for neb in nebs:
            print("Testing neb: ", neb)
            nn = NearestNeighbors(n_neighbors=neb, metric=metric)
            cv = KFold(n_splits=10, shuffle=True)
            scoring = cross_val_score(estimator = nn, X = dat, y = labs,
                                      scoring = 'accuracy', cv = cv)








if __name__ == "__main__":
    main()
