import numpy as np
import tabulate

# Test Sequences:
# S1: ATCTGTTCTT
# S2: CTTTATGTT

# Build our Matrix in a Pandas Dataframe.

def get_seqs():
    #seq1 = input("Please enter your first sequence for alignment.\n")
    #seq2 = input("Now, please enter your second seqeuence for alignment.\n")

    # For testing:
    seq1 = "ATCTGTTCTT"
    #seq1 = "ATCTGTTCT"
    #seq2 = "TATGT"
    seq2 = "CTTTATGTT"
    return (seq1, seq2)


# Build the DP table
def build_tab(seq1, seq2):

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
            if (x == 0) or (y == 0):
                tab[x,y] = 0

    #pretty_print(seq1, seq2, tab, dir)

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

            #pretty_print(seq1, seq2, tab, dir)
            #print("================================================")
    return (tab,dir)


def align(seq1, seq2, tab, dir):

    al_seq1 = ""
    al_seq2 = ""

    # First thing we need to do: find the max in the last col.
    col = tab[:,(tab.shape[1] -1)]
    x_max = max(col)
    x_max_ind = np.where(col == x_max)
    x_max_ind = x_max_ind[0][-1]

    row = tab[(tab.shape[0] -1),:]
    y_max = max(row)
    y_max_ind = np.where(row == y_max)
    y_max_ind = y_max_ind[0][-1]

    # Need to loop backwards through x and y.
    x = tab.shape[0] - 1
    y = tab.shape[1] - 1


    # Starting from last row
    if y_max > x_max:
        # Loop backwards till we hit first row or col
        while x > 0 and y > 0:

            print("Seq1: ", al_seq1)
            print("Seq2: ", al_seq2)
            print("===========================")

            # Prepend Gap on Sequence 2
            if y > y_max_ind:
                al_seq1 = seq1[y] + al_seq1
                al_seq2 = "-" + al_seq2
                # decrement y
                y -= 1
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

        print(al_seq1)
        print(al_seq2)


    # Starting from last column
    else:
        # Loop backwards till we hit first row or col
        while x > 0 and y > 0:

            print("Seq1: ", al_seq1)
            print("Seq2: ", al_seq2)

            # Prepend Gap on Sequence 2
            if x > x_max_ind:
                al_seq1 = seq1[y] + al_seq1
                al_seq2 = "-" + al_seq2
                # decrement y
                y -= 1
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

        print(al_seq1)
        print(al_seq2)


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

    return (al_seq1, al_seq2)

def pretty_print(seq1, seq2, tab, dir):

    s1 = list(seq1)
    s2 = list(seq2)

    print(tabulate.tabulate(tab, headers=s1, showindex=s2, tablefmt="grid"))
    print("--------------------------------------")
    print(tabulate.tabulate(dir, headers=s1, showindex=s2, tablefmt="grid"))



def main():

    seqs = get_seqs()
    seq1 = seqs[0]
    seq2 = seqs[1]

    seq1 = "-" + seq1
    seq2 = "-" + seq2

    print("Sequence 1: ", seq1)
    print("Sequence 2: ", seq2)



    build = build_tab(seq1, seq2)

    tab = build[0]
    dir = build[1]

    # Pretty print
    # pretty_print(seq1, seq2, tab, dir)


    al = align(seq1, seq2, tab, dir)
    align_seq1 = al[0]
    align_seq2 = al[1]

    print("=========================")
    print(align_seq1)
    print(align_seq2)
    print("=========================")

    # Run other way
    build = build_tab(seq2, seq1)

    tab = build[0]
    dir = build[1]


    al = align(seq2, seq1, tab, dir)
    align_seq1 = al[0]
    align_seq2 = al[1]
    print("=========================")
    print(align_seq1)
    print(align_seq2)
    print("=========================")

if __name__ == "__main__":
    main()