import numpy as np
import sys
#This program assumes that the sequences to be inserted are found on the second line (1st index)
#of the file that they come from. The example files given all follow this format
#program can be run on windows commandline with nw_algorithm.py [query file] [ref file] [output file]

def systemInput():#systemInput checks if the function was called on commandline with a input file, if not it defaults to pre-selected file
    defaultFile = "pfizer_mrna_truncated_example.fna"
    try:
        if sys.argv[1] is None:
            return defaultFile
        else:
            return sys.argv[1]
    except IndexError as e:
        return defaultFile
def systemInput2():#checks second commandline argument
    defaultFile = "sars_spike_protein_truncated_example.fna"
    try:
        if sys.argv[2] is None:
            return defaultFile
        else:
            return sys.argv[2]
    except IndexError as e:
        return defaultFile

def fileToOutput():#fileToOutput checks the commandline for an outputfile, if none was chosen then it will default to a preselected output file
    defaultOutputFile = 'Example_output.fna'
    try:
        if sys.argv[3] is None:
            return defaultOutputFile
        else:
            return sys.argv[3]
    except IndexError as e:
        return defaultOutputFile

#this will create a static 11 gaps, optionally the first line can be uncommented
#and the second line commented, this will make the gap length represent how many gaps
#were actually created
def create_gaps(total_length,sub_length):
    #i = total_length - sub_length
    i = 10
    gaps = ""
    while i >= 0:
        gaps += "_"
        i -= 1
    return gaps
def create_spaces(total_length,sub_length):
    #i = total_length - sub_length
    i = 10
    gaps = ""
    while i >= 0:
        gaps += " "
        i -= 1
    return gaps
    
#Needleman-Wunsch algorithm, gap, mismatch, and match parameters optional.
#if you want Init/terminal gap penalty of 0 pass in False for InitPenalty parameter
#it will default to counting all gap penalties if not specified
def dynamic_alignment(query,ref,InitPenalty = True,gap = -2, mismatch = -1, match = 1):
    #lambda function to check if diagonal movement is a match or mismatch
    mis_or_match = lambda query,ref,i,j: mismatch if ref[i] != query[j] else match

    #sets first row and column values to 0 if InitPenalty = false
    adjusted_gap = 0
    #if InitPenalty is true then gaps are counted normally
    if InitPenalty:
        adjusted_gap = gap
    #if not true then start and end gaps are counted as 0
    else:
        adjusted_gap = 0
    #create a matrix with length ref +1 and query + 1. initializes to 0
    matrix = np.zeros((len(ref)+1,len(query)+1))
    
    #adds the gap penalty to the first row and column
    for i in range(len(ref)+1):
        matrix[i][0] = i* adjusted_gap
    for j in range(len(query)+1):
        matrix[0][j] = j * adjusted_gap


    #solving for each row, begin at index 1
    for i in range(1, len(ref)+1):
        for j in range(1, len(query)+1):
            #find max between diagonal, delete, or insert
            matrix[i][j] = max(
                            matrix[i-1][j-1] + mis_or_match(query,ref,i-1,j-1),
                            matrix[i-1][j] + gap,
                            matrix[i][j-1] + gap)
    #solution matrix has been fully filled out now

    #print matrix for debugging
    #for line in matrix:
    #    print(line)

    #TRACEBACK
    i = len(ref)
    j = len(query)
    #strings for output and the relationship between them (x _ |)
    alignment1 = relationship = alignment2 =""

    ##debugging variables
    match_count =0
    mismatch_count = 0
    gap_count = 0
    #print("diag, up, left")
    ##end debugging variables

    #If init gap penalties arent counted
    if InitPenalty == False:
        #appendRef will later help determine which alignment gets gaps appended to the front or back
        appendRef = True
        #tracks the largest value in row and column respectively
        maxInRow = matrix[i][j]
        maxInColumn = matrix[i][j]
        row_i = i
        col_j = j
        #finds starting location of traceback by locating index of the largest value in final column or row
        for n in range(0, len(query)+1):
            if maxInRow < matrix[len(ref)][n]:
                maxInRow = matrix[len(ref)][n]
                col_j = n
                
        for m in range(0,len(ref)+1):
            if maxInColumn < matrix[m][len(query)]:
                maxInColumn = matrix[m][len(query)]
                row_i = m
        #the score is the largest value found
        score = max(maxInRow,maxInColumn)
        #this if else block creates a string of gaps and spaces to append to the sequences based on starting location
        if maxInRow > maxInColumn:
            gap_string = create_gaps(len(ref),row_i)
            space_string = create_spaces(len(ref),row_i)
            #i = len(ref) selects the final row, j selects the column location
            i = len(ref)
            j = col_j
        else:
            appendRef = False
            gap_string = create_gaps(len(query),col_j)
            space_string = create_spaces(len(query),col_j)
            #j = len(query) selects the final column, i selects the row location
            i = row_i
            j = len(query)

    #If init gap penalties are counted - Normal NW Algorithm
    else:
        #score is bottom right value of matrix
        score = matrix[i][j]

    #while loop begins traceback based on i and j location. If init gap penalties arent counted the start location is varied.
    #if init gap penalties are counted begin from bottom right
    while i > 0 and j > 0:
        diagonal = matrix[i-1][j-1]
        up = matrix[i-1][j]
        left = matrix[i][j-1]
        #print(diagonal,up,left, end = " ")
        #finds the highest value (best) of the three options and chooses it
        getPath = max(diagonal,up,left)
        getDiag = max(diagonal,up+gap,left+gap)
        #this if block checks if the diagonal element is a match, if it is it bypasses the getPath check because matches should automatically take priority.
        #checks the diagonal against the other two options adjusted for their gap penalty
        if ref[i-1] == query[j-1] and getDiag == diagonal:
            #if its a match then relationship is a |
            alignment1 = ref[i-1] + alignment1
            relationship = "|" + relationship
            alignment2 = ref[i-1] + alignment2
            match_count += 1
            #decrement both
            i -= 1
            j -= 1
        #if getPath is the highest value option we already know its a mismatch due to prior check.
        elif getPath == diagonal:
            #print("diagonal")
            #if not a match then relationship is an x
            alignment1 = ref[i-1] + alignment1
            relationship = "x" + relationship
            alignment2 = query[j-1] + alignment2
            mismatch_count += 1
            #decrement both    
            i -=1
            j -=1
        #if getPath is up we insert a gap in alignment 2
        elif getPath == up:
            #print("up")
            alignment1 = ref[i-1] + alignment1
            relationship = " " + relationship
            alignment2 = "_" + alignment2
            # decrement row pointer
            gap_count += 1
            i -= 1
        #if getPath is left we insert a gap in alignment 1 
        elif getPath == left:
            #print("left")
            alignment1 = "_" + alignment1
            relationship = " " + relationship
            alignment2 = query[j-1] + alignment2
            gap_count += 1
            # decrement column pointer
            j -= 1

   #if we dont count start/terminal gaps then we have to adjust the alignments for extra gaps
    if InitPenalty == False:
        relationship = space_string + relationship + space_string
        if appendRef:
            alignment2 = gap_string + alignment2
            alignment1 = space_string + alignment1+gap_string
        else:
            alignment1 = gap_string + alignment1
            alignment2 = space_string + alignment2+gap_string

    #print("")
    #lines below would count the mismatches, matches, and gaps and then add them together to check if they matched the score. this indicates whether traceback was successful.
    #print("matches:",match_count,"mismatches:",mismatch_count,"gaps:", gap_count)
    #print("computed score:", (match_count*match)+(mismatch_count*mismatch)+
    #      (gap_count*gap))

    #print("expected score is:", score)
    #print the alignments and their relationship
    #print(alignment1)
    #print(relationship)
    #print(alignment2)

    
    #converts score from numpy.float64 to int so it can be written in file
    file_score = str(score)
    f = open(fileToOutput(), 'w+', newline ='')
    # writing the data into the file
    f.write(file_score)
    f.write("\n")
    f.write(systemInput2())
    f.write("\n")
    f.write(alignment1)
    f.write("\n")
    f.write(relationship)
    f.write("\n")
    f.write(alignment2)
    f.write("\n")
    f.write(systemInput())
    f.close()

##MAIN
if __name__ == '__main__':
    #uses systemInput functions to read in files and generate arguments
    with open(systemInput(), "r+") as query_input:
        #strips lines of \n and places them in an array
        Lines = [line.strip() for line in query_input.readlines()]
    with open(systemInput2(), "r+") as ref_input:
        #strips lines of \n and places them in an array
        Lines2 = [line.strip() for line in ref_input.readlines()]
    #gets query sequence from second line of first file
    query = Lines[1]
    #gets ref sequence from second line of second file
    ref = Lines2[1]
    
    #query = "atgct"
    #ref = "agct"
    
    #change the third parameter to False for no init/teminal gap penalty
    #change the third parameter to True for init/terminal penalties to count
    dynamic_alignment(query, ref, True)
    
