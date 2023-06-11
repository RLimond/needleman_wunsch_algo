import csv
import sys
#
#
#THIS PROGRAM CONVERTS CODON SEQUENCES TO AMINO ACIDS STARTING FROM THE FIRST START CODON ENDING WITH THE FIRST STOP CODON
#
#
def systemInput():#systemInput checks if the function was called on commandline with a input file, if not it defaults to pre-selected file
    defaultFile = "to_be_converted_to_amino.fna"
    try:
        if sys.argv[1] is None:
            return defaultFile
        else:
            return sys.argv[1]
    except IndexError as e:
        return defaultFile

def fileToOutput():#fileToOutput checks the commandline for an outputfile, if none was chosen then it will default to a preselected output file
    defaultOutputFile = 'sars_spike_protein.aa'
    try:
        if sys.argv[2] is None:
            return defaultOutputFile
        else:
            return sys.argv[2]
    except IndexError as e:
        return defaultOutputFile
def getSubStrings(DNA, position = 0):
    return [DNA[i:i+3] for i in range(position, len(DNA) - 2, 3)]
DNA_Codons = {
    "GCT": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "TGT": "Cys", "TGC": "Cys",
    "GAT": "Asp", "GAC": "Asp",
    "GAA": "Glu", "GAG": "Glu",
    "TTT": "Phe", "TTC": "Phe",
    "GGT": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",
    "CAT": "His", "CAC": "His",
    "ATA": "Ile", "ATT": "Ile", "ATC": "Ile",
    "AAA": "Lys", "AAG": "Lys",
    "TTA": "Leu", "TTG": "Leu", "CTT": "Leu", "CTC": "Leu", "CTA": "Leu", "CTG": "Leu",
    "ATG": "Met",
    "AAT": "Asn", "AAC": "Asn",
    "CCT": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "CAA": "Gln", "CAG": "Gln",
    "CGT": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg", "AGA": "Arg", "AGG": "Arg",
    "TCT": "Ser", "TCC": "Ser", "TCA": "Ser", "TCG": "Ser", "AGT": "Ser", "AGC": "Ser",
    "ACT": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "GTT": "Val", "GTC": "Val", "GTA": "Val", "GTG": "Val",
    "TGG": "Trp",
    "TAT": "Tyr", "TAC": "Tyr",
    "TAA": "Stp", "TAG": "Stp", "TGA": "Stp"
    }

with open(systemInput(), "r+") as writeF:
    # Reading from a file
    Lines = [line.strip() for line in writeF.readlines()]#strips lines of \n and places them in an array
    Codons = []#stores individual codons ('AAA')
    Codon_Lines = []#stores the sequence(s)
    i = 0
    #This implementation is a bit lazy. Should be fixed at a later date.
    for line in Lines:
        if (i%2 == 0): #the genome/genes are found on odd lines only in the .fna (even list index). If list index is not even we ignore the line
            Codon_Lines += [line]
            i = i + 1
        else:
            i = i + 1
    Codon_Lines.pop(0) #currently each seperate gene is in its own element in a list
    Codon_Lines[0:len(Codon_Lines)] = [''.join(Codon_Lines[0:len(Codon_Lines)])] #separate codons have been merged into single list element
    Codons = getSubStrings(Codon_Lines[0])
    #Codons.sort()
    print(Codons)

    tracker = 0
    for element in Codons:
        if element == 'ATG':
            break;
        else:
            tracker += 1
    stop_tracker = 0
    for elm in Codons:
        if elm == 'TAA' or elm == 'TAG' or elm == 'TGA':
            break;
        else:
            stop_tracker += 1
    Amino = ""
    for codon in range(tracker, stop_tracker+1):
        Amino = Amino + (DNA_Codons[Codons[codon]])
    #print(Amino)

    f = open(fileToOutput(), 'w+', newline ='')
    f.write(fileToOutput())
    f.write(" converted to amino")
    f.write("\n")
    f.write(Amino)
    f.close()
    
        
        
    
