
# coding: utf-8

# # Bowtie Prep
# Bowtie(1.2.2) was downloaded.
# Built indexes using a lower offrate value which will increase the size of index on disk and in memory but should help with the alignment when using the -a command.
#     python c:/bowtie/bowtie-build --offrate 3 ..\data\ames\ames.fasta  ..\data\ames\ames    
#     python c:/bowtie/bowtie-build --offrate 3 ..\data\sterne\sterne.fasta  ..\data\sterne\sterne 
#     
#     
# # Generate tandem repeat library

# In[3]:


import itertools as it

# Define our library of chars and build our tandem repeats strings
lib = ('A', 'T', 'C','G')


# Repeat Length (n)  ==> Number of possible tandem repeats (len(lib)^n = 4^n)
# 6  ==>        4,096
# 7  ==>       16,384
# 8  ==>       65,536
# 9  ==>      262,144
# 10 ==>    1,048,576

# This could be shrunk into 1 big for loop to output 1 list of all lenghts
# since bowtie outputs the matched query, you could collapse all the outputs into a map
# match1 -> index1, index2, index3
# match2 -> index4, index5
# and then you just need to take the len(match1) to determine the jump distance to determine if 
# index1 and 2 are close enough to be considered tandem repeats

# generate tandem repeat library
six = []
for i in it.product(lib, repeat=6):
    six.append(''.join(map(str, i)))

sev = []
for i in it.product(lib, repeat=7):
    sev.append(''.join(map(str, i)))

eight = []
for i in it.product(lib, repeat=8):
    eight.append(''.join(map(str, i)))

nine = []
for i in it.product(lib, repeat=9):
    nine.append(''.join(map(str, i)))
    
ten = []
for i in it.product(lib, repeat=10):
    ten.append(''.join(map(str, i)))


# Function to output results to individual files
def out(filename, dat):
    f = open(filename, "w")
    for i in dat:
        f.write("\'"+i+"\'\n")
    f.close()

out("data/lib/six", six)
out("data/lib/sev", sev)
out("data/lib/eight", eight)
out("data/lib/nine", nine)
out("data/lib/ten", ten)


# # Searching for tandem repeats
# We have our tandem repeat library and our bowtie indices of interest.
# We are only interested in tandem repeats of length 3+ (that is 3 or more repeats of the same library entry)
# 
# bowtie options:
# --suppress 1,2,3,6,7,8 [Leave only the columns we are intrested in, query and loc]
# -t [time to process]
# -p 4 [use 4 threads]
# -a   [find all matches]
# -v 0 [no mismatches allowed]
# -r input file [our file that contains the queries]
# 
# #### Ames
# python c:\bowtie\bowtie ames/ames -r lib/six -a -v 0 -t -p 4 --suppress 1,2,3,6,7,8 > ames/out/six.out 2> ames/out/six.err  
# Time loading forward index: 00:00:00
#     Time for 0-mismatch search: 00:00:13
#     #reads processed: 4096
#     #reads with at least one reported alignment: 4096 (100.00%)
#     #reads that failed to align: 0 (0.00%)
#     Reported 14139938 alignments
#     Time searching: 00:00:13
#     Overall time: 00:00:13
# 
# python c:\bowtie\bowtie ames/ames -r lib/sev -a -v 0 -t -p 4 --suppress 1,2,3,6,7,8 > ames/out/sev.out 2> ames/out/sev.err  
#     Time loading forward index: 00:00:00
#     Time for 0-mismatch search: 00:00:13
#     #reads processed: 16384
#     #reads with at least one reported alignment: 16384 (100.00%)
#     #reads that failed to align: 0 (0.00%)
#     Reported 14167641 alignments
#     Time searching: 00:00:13
#     Overall time: 00:00:13
# 
# python c:\bowtie\bowtie ames/ames -r lib/eight -a -v 0 -t -p 4 --suppress 1,2,3,6,7,8 > ames/out/eight.out 2> ames/out/eight.err  
#     Time loading forward index: 00:00:00
#     Time for 0-mismatch search: 00:00:13
#     #reads processed: 65536
#     #reads with at least one reported alignment: 65530 (99.99%)
#     #reads that failed to align: 6 (0.01%)
#     Reported 14177830 alignments
#     Time searching: 00:00:13
#     Overall time: 00:00:13
# 
# python c:\bowtie\bowtie ames/ames -r lib/nine -a -v 0 -t -p 4 --suppress 1,2,3,6,7,8 > ames/out/nine.out 2> ames/out/nine.err  
#     Time loading forward index: 00:00:00
#     Time for 0-mismatch search: 00:00:15
#     #reads processed: 262144
#     #reads with at least one reported alignment: 260224 (99.27%)
#     #reads that failed to align: 1920 (0.73%)
#     Reported 14181366 alignments
#     Time searching: 00:00:15
#     Overall time: 00:00:15
# 
# python c:\bowtie\bowtie ames/ames -r lib/ten -a -v 0 -t -p 4 --suppress 1,2,3,6,7,8 > ames/out/ten.out 2> ames/out/ten.err  
#     Time loading forward index: 00:00:00
#     Time for 0-mismatch search: 00:00:18
#     #reads processed: 1048576
#     #reads with at least one reported alignment: 950147 (90.61%)
#     #reads that failed to align: 98429 (9.39%)
#     Reported 14181700 alignments
#     Time searching: 00:00:18
#     Overall time: 00:00:18
# 
# 
# #### Sterne
# python c:\bowtie\bowtie sterne/sterne -r lib/six -a -v 0 -t -p 4 --suppress 1,2,3,6,7,8 > sterne/out/six.out 2> sterne/out/six.err  
#     Time loading forward index: 00:00:00
#     Time for 0-mismatch search: 00:00:13
#     #reads processed: 4096
#     #reads with at least one reported alignment: 4096 (100.00%)
#     #reads that failed to align: 0 (0.00%)
#     Reported 14143170 alignments
#     Time searching: 00:00:13
#     Overall time: 00:00:13
# 
# python c:\bowtie\bowtie sterne/sterne -r lib/sev -a -v 0 -t -p 4 --suppress 1,2,3,6,7,8 > sterne/out/sev.out 2> sterne/out/sev.err  
#     Time loading forward index: 00:00:00
#     Time for 0-mismatch search: 00:00:13
#     #reads processed: 16384
#     #reads with at least one reported alignment: 16384 (100.00%)
#     #reads that failed to align: 0 (0.00%)
#     Reported 14170873 alignments
#     Time searching: 00:00:13
#     Overall time: 00:00:13
# 
# 
# python c:\bowtie\bowtie sterne/sterne -r lib/eight -a -v 0 -t -p 4 --suppress 1,2,3,6,7,8 > sterne/out/eight.out 2> sterne/out/eight.err  
#     Time loading forward index: 00:00:00
#     Time for 0-mismatch search: 00:00:13
#     #reads processed: 65536
#     #reads with at least one reported alignment: 65530 (99.99%)
#     #reads that failed to align: 6 (0.01%)
#     Reported 14181071 alignments
#     Time searching: 00:00:13
#     Overall time: 00:00:13
# 
# python c:\bowtie\bowtie sterne/sterne -r lib/nine -a -v 0 -t -p 4 --suppress 1,2,3,6,7,8 > sterne/out/nine.out 2> sterne/out/nine.err  
#     Time loading forward index: 00:00:00
#     Time for 0-mismatch search: 00:00:14
#     #reads processed: 262144
#     #reads with at least one reported alignment: 260230 (99.27%)
#     #reads that failed to align: 1914 (0.73%)
#     Reported 14184603 alignments
#     Time searching: 00:00:14
#     Overall time: 00:00:14
# 
# 
# python c:\bowtie\bowtie sterne/sterne -r lib/ten -a -v 0 -t -p 4 --suppress 1,2,3,6,7,8 > sterne/out/ten.out 2> sterne/out/ten.err  
#     Time loading forward index: 00:00:00
#     Time for 0-mismatch search: 00:00:14
#     #reads processed: 1048576
#     #reads with at least one reported alignment: 950217 (90.62%)
#     #reads that failed to align: 98359 (9.38%)
#     Reported 14184944 alignments
#     Time searching: 00:00:14
#     Overall time: 00:00:14

# # Output
# 1st column is index of where the match starts.
# 2nd column is the matching query sequence.
# 
# ### Example:
# 280045	AACAA
# 290018	AACAA
# 9256	AACAA
# 741380	AACAA
# 145436	AACAA
# 29050	AACAA
# 532426	AACAA
# 82371	AACAA
# 266746	AACAA
# 

# In[4]:


# This function takes a file path (please have r appended to the front for paths on windows)
# and generates a repeat dictionary.  This repeat dictionary has sequences as the key
# and an ordered list with no dupes as the value.  This ordered list contains each location
# in the genome that this sequence was found in asc order.

def create_rep_dict(filepath):
    file = filepath
    # dictionary that will use repeat sequence as key, and keep location values
    tan_dict = dict()

    with open(file, "r") as input:
        line = input.read().splitlines()
        # Loop through each line in file and deal with it
        for i in line:
            val = i.split()
            seq = val[1]
            loc = int(val[0])
            if seq in tan_dict:
               # append the new loc to the existing array for this seq.
                tan_dict[seq].append(loc)
            else:
               # create a new array for this seq
                tan_dict[seq] = [loc]

        # Ugly probably slow hack, removes all seqs with less than 3 locs since they couldn't be a tandem repeat
        # also sorts and uniques the list of locs for seqs
        for k in list(tan_dict.keys()):
            if(len(tan_dict[k]) > 2):
                tan_dict[k] = list(sorted(set(tan_dict[k])))
            else:
                del tan_dict[k]

    return tan_dict


# In[5]:


# This function takes the a prepared tandem repeat dictionary.
# This dictionary needs to have the key set as a sequence, with the 
# value for each key being an ordered list with no dupes of the locations
# within the genome the seq was found.

# It determines the number max tandem repeat seen in the genome for
# a particular sequence, and returns a dict with the seq as they key
# and the number of repeats as the value.

# Note that this function only returns the max number of repeats so if your
# seq has to separate tandem repeats, 1 of 3 repeats and another of 7, it will 
# only return a value of 7.

def finds_reps(tan_rep_dict):
    REPEAT_LEN = len(list(tan_rep_dict.keys())[-1])

    tan_dict = tan_rep_dict

    # Loop through all sequences
    for k in list(tan_dict.keys()):
        # assign the list of locs for this seq to locs
        locs = tan_dict[k]
        rep_len = [0]
        count = 1
        max_count = 1
        prev_loc = locs[0]

        for loc in locs:
            dist = loc - prev_loc

            # Too short to be a repeat
            if dist < REPEAT_LEN:
                continue
            # We have a repeat
            if dist == REPEAT_LEN:
                count += 1

            if dist > REPEAT_LEN:
                # Store this count if its greater than our longest count seen so far and reset
                max_count = max(max_count, count)
                count = 1

            prev_loc = loc


        # no tandem repeats for this sequence so drop it from the dict
        # otherwise if there is replace the locs with the max number of repeats seen.
        if max_count < 3:
            del tan_dict[k]
        else:
            tan_dict[k] = max_count

    return tan_dict




# In[6]:


# Handling output.  Output files need to be read in for processing.

# dict to hold our identified ames tandem repeats
ames_reps = dict()

# go through each file 6-10 adding tandem repeats to ames_reps dict
# 6
reps_dict = create_rep_dict(r'C:\Users\Mike\Desktop\Mike\UC\Spring_2018\BMIN7099_Bioinformatics\hw\homework_1\data\ames\out\six.out')
ames_reps.update(finds_reps(reps_dict))

# 7
reps_dict = create_rep_dict(r'C:\Users\Mike\Desktop\Mike\UC\Spring_2018\BMIN7099_Bioinformatics\hw\homework_1\data\ames\out\sev.out')
ames_reps.update(finds_reps(reps_dict))

# 8
reps_dict = create_rep_dict(r'C:\Users\Mike\Desktop\Mike\UC\Spring_2018\BMIN7099_Bioinformatics\hw\homework_1\data\ames\out\eight.out')
ames_reps.update(finds_reps(reps_dict))

# 9
reps_dict = create_rep_dict(r'C:\Users\Mike\Desktop\Mike\UC\Spring_2018\BMIN7099_Bioinformatics\hw\homework_1\data\ames\out\nine.out')
ames_reps.update(finds_reps(reps_dict))

# 10
reps_dict = create_rep_dict(r'C:\Users\Mike\Desktop\Mike\UC\Spring_2018\BMIN7099_Bioinformatics\hw\homework_1\data\ames\out\ten.out')
ames_reps.update(finds_reps(reps_dict))



# dict to hold our identified sterne tandem repeats
sterne_reps = dict()

# go through each file 6-10 adding tandem repeats to sterne_reps dict
# 6
reps_dict = create_rep_dict(r'C:\Users\Mike\Desktop\Mike\UC\Spring_2018\BMIN7099_Bioinformatics\hw\homework_1\data\sterne\out\six.out')
sterne_reps.update(finds_reps(reps_dict))

# 7
reps_dict = create_rep_dict(r'C:\Users\Mike\Desktop\Mike\UC\Spring_2018\BMIN7099_Bioinformatics\hw\homework_1\data\sterne\out\sev.out')
sterne_reps.update(finds_reps(reps_dict))

# 8
reps_dict = create_rep_dict(r'C:\Users\Mike\Desktop\Mike\UC\Spring_2018\BMIN7099_Bioinformatics\hw\homework_1\data\sterne\out\eight.out')
sterne_reps.update(finds_reps(reps_dict))

# 9
reps_dict = create_rep_dict(r'C:\Users\Mike\Desktop\Mike\UC\Spring_2018\BMIN7099_Bioinformatics\hw\homework_1\data\sterne\out\nine.out')
sterne_reps.update(finds_reps(reps_dict))

# 10
reps_dict = create_rep_dict(r'C:\Users\Mike\Desktop\Mike\UC\Spring_2018\BMIN7099_Bioinformatics\hw\homework_1\data\sterne\out\ten.out')
sterne_reps.update(finds_reps(reps_dict))





# In[7]:


import pandas as pd
import tabulate as tab

ames_set = set(ames_reps.keys())
sterne_set = set(sterne_reps.keys())

inter = ames_set & sterne_set
print("Number of tandem repeats in Ames: ", len(ames_set))
print("Number of tandem repeats in Sterne: ", len(sterne_set))
print("Number of tandem repeats in existing in both: ", len(inter))



# In[ ]:


inter_reps = dict()

for i in inter:
    # Loop through and save any where the number seen in ames and stern differ
    if ames_reps[i] != sterne_reps[i]:
        inter_reps[i] = (ames_reps[i], sterne_reps[i])
        print(len(i))

print("Number of tandem repeats seen at different lengths: ", len(inter_reps))



# In[18]:


import pandas as pd
import tabulate as tab
# Generate sets of the keys so we can intersect

ames_set = set(ames_reps.keys())
sterne_set = set(sterne_reps.keys())

print(ames_reps)

print(sterne_reps)
inter_reps = dict()

for i in inter:
    # Loop through and save any where the number seen in ames and stern differ
    if ames_reps[i] != sterne_reps[i]:
        inter_reps[i] = (ames_reps[i], sterne_reps[i])

df = pd.DataFrame(inter_reps)
df = df.transpose()
df['Difference'] = df[0] - df[1]     
total = df['Difference'].sum()
head = ['Tandem Repeat', 'Ames Len', 'Sterne Len', 'Difference']

print(tab.tabulate(df, headers=head))
print("Total (+ means ames is longer, - means Sterne is longer)", total)

