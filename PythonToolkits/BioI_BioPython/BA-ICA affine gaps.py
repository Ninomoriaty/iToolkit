# Perform Smith Waterman Alignment in Python (from first principles)

# contains rows lists each of length cols initially set to 0
# index as my_matrix[1][2] my_matrix[R][C]
def create_matrix(sequence1, sequence2, D, E):
    rows=len(sequence2)
    cols=len(sequence1)
    # Original create_matrix section
    def x(i, j):
        if i > 0 and j == 0:
            return MIN
        elif j > 0:
            return D + (E * j)
        else:
            return 0

    def y(i, j):
        if j > 0 and i == 0:
            return MIN
        elif i > 0:
            return D + (E * i)
        else:
            return 0

    def m(i, j):
        if j == 0 and i == 0:
            return 0
        elif j == 0 or i == 0:
            return MIN
        else:
            return 0

    # initialise values and matrices
    m_matrix = [[m(col, row) for col in range(cols+1)] for row in range(rows+1)]
    x_matrix = [[x(col, row) for col in range(cols+1)] for row in range(rows+1)]
    y_matrix = [[y(col, row) for col in range(cols+1)] for row in range(rows+1)]

    return m_matrix, x_matrix, y_matrix

def match(seq1, seq2, i, j):
    if seq2[i-1] == seq1[j-1]:
        return seqmatch
    else:
        return seqmismatch

# x is row index, y is column index
# follows[r][c]
def calc_score(matrix, x, y):
    # print("seq1:",sequence1[y- 1]," seq2: "+sequence2[x - 1],"x:",x," y:",y)

    sc = seqmatch if sequence1[y - 1] == sequence2[x - 1] else seqmismatch

    base_score = matrix[x - 1][y - 1] + sc
    insert_score = matrix[x - 1][y] + seqgap
    delete_score = matrix[x][y - 1] + seqgap
    v = max(0, base_score, insert_score, delete_score)
    return v


# makes a single traceback step
def traceback(mymatrix, maxv):
    x = maxv[0]
    y = maxv[-1]
    val = mymatrix[x][y]

    # oops bug...now fixed
    sc = seqmatch if sequence1[y - 2] == sequence2[x - 2] else seqmismatch
    base_score = mymatrix[x - 1][y - 1] + sc
    if base_score == val:
        return [x - 1, y - 1]

    insert_score = mymatrix[x - 1][y] + seqgap
    if insert_score == val:
        return [x - 1, y]
    else:
        return [x, y - 1]


# builds the initial scoring matrix used for traceback
def build_matrix(m_matrix, x_matrix, y_matrix, sequence1, sequence2, D, E):
    rows=len(sequence2)
    cols=len(sequence1)
    for j in range(1, cols+1):
        for i in range(1, rows+1):
            x_matrix[i][j] = max((D + E + m_matrix[i][j-1]), (E + x_matrix[i][j-1]), (D + E + y_matrix[i][j-1]))
            y_matrix[i][j] = max((D + E + m_matrix[i-1][j]), (D + E + x_matrix[i-1][j]), (E + y_matrix[i-1][j]))
            m_matrix[i][j] = max(match(sequence1, sequence2, i, j) + m_matrix[i-1][j-1], x_matrix[i][j], y_matrix[i][j])

    return m_matrix, x_matrix, y_matrix


# gets the max value from the built matrix
def get_max(mymatrix):
    max = mymatrix[0][0]
    mrow = 0
    mcol = 0

    rows = len(mymatrix)
    cols = len(mymatrix[0])

    for i in range(1, rows):
        for j in range(1, cols):
            if mymatrix[i][j] > max:
                max = mymatrix[i][j]
                mrow = i
                mcol = j
    print("max score: ", max)
    return [mrow, mcol]


# print out the best scoring path from the SW matrix
def print_matrix(mymatrix):
    rows = len(mymatrix)
    cols = len(mymatrix[0])
    s1 = "  " + sequence1
    s2 = " " + sequence2

    print("Dimensions: r= %2d , c= %2d" % (rows, cols))

    for a in s1:
        print(a, end="")
        print(" \t", end="")
    print("\n", end="")

    for i in range(0, rows):
        print(s2[i], end="")
        print(" \t", end="")
        for j in range(0, cols):
            print("%02d\t" % (mymatrix[i][j]), end="")
        print("\n", end="")


# print out the traceback of the best scoring alignment
def print_traceback(seq1, seq2, m_matrix, x_matrix, y_matrix):
   #traverse the matrix to find the traceback elements
   #if more than one path just pick one
    topstring=""
    midstring=""
    bottomstring=""

    i = len(seq2)
    j = len(seq1)
    while (i > 0 or j > 0):
        if (i > 0 and j > 0 and m_matrix[i][j] == m_matrix[i-1][j-1] + match(seq1, seq2, i, j)):
            topstring += seq1[j-1]
            bottomstring += seq2[i-1]
            midstring += "|"
            i -= 1; j -= 1
        elif (i>0 and m_matrix[i][j] == y_matrix[i][j]):
            topstring += '-'
            bottomstring += seq2[i-1]
            midstring += " "
            i -= 1
        elif (j>0 and m_matrix[i][j] == x_matrix[i][j]):
            topstring += seq1[j-1]
            bottomstring += '-'
            midstring += " "
            j -= 1

    topstringr = ''.join([topstring[j] for j in range(-1, -(len(topstring)+1), -1)])
    midstringr = ''.join([midstring[j] for j in range(-1, -(len(midstring)+1), -1)])
    bottomstringr = ''.join([bottomstring[j] for j in range(-1, -(len(bottomstring)+1), -1)])

    print(topstringr)
    print(midstringr)
    print(bottomstringr)


# build the SW alignment...
def perform_smith_waterman(D=-5, E=-0.2):
    #values for weights
    global seqmatch
    global seqmismatch
    global seqgap
    global sequence1
    global sequence2
    global MIN


    #note these are not the exact weights used in the original SW paper
    seqmatch =1
    seqmismatch=-1
    seqgap=-1

    #input sequences
    #sequence1="AGTGATAAACTAGTAATTTTT"
    #sequence2="TTGGGGGTAAACAGGGG"

    sequence1 = "GTGTATT"
    sequence2 = "GTGTTATT"

    print("Sequence1: "+sequence1)
    print("Sequence2: "+sequence2)
    MIN = max(len(sequence1) * E + D, len(sequence2) * E + D)

    m_matrix, x_matrix, y_matrix=create_matrix(sequence1, sequence2, D, E)
    m_matrix, x_matrix, y_matrix=build_matrix(m_matrix, x_matrix, y_matrix, sequence1, sequence2, D, E)
    for mymatrix in m_matrix, x_matrix, y_matrix:
        print_matrix(mymatrix)

    print_traceback(sequence1, sequence2, m_matrix, x_matrix, y_matrix)

##this calls the SW algorithm when the script loads
perform_smith_waterman()
