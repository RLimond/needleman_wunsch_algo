# Background
The Needleman-Wunsch algorithm is a classic dynamic programming algorithm used in bioinformatics for sequence alignment. It was developed by Saul Needleman and Christian Wunsch in 1970 and has since become a fundamental tool in computational biology.

Sequence alignment is the process of arranging two or more biological sequences (such as DNA, RNA, or protein sequences) to identify regions of similarity or homology. The Needleman-Wunsch algorithm specifically focuses on global sequence alignment, which means it attempts to align the entire length of two sequences.

The algorithm uses a dynamic programming approach to calculate an optimal alignment score and to determine the actual alignment of the sequences. It constructs a matrix, often referred to as the scoring matrix or similarity matrix, to evaluate the alignment possibilities between individual characters or elements of the sequences.
The scoring matrix assigns scores to different match, mismatch, and gap penalties. A match score represents a positive score when two elements match, while a mismatch score represents a penalty when two elements do not match. Gap penalties are applied when introducing gaps (insertions or deletions) into the alignment. The algorithm aims to find the alignment that maximizes the overall score based on these penalties.

# Additional Information and Running the program
This folder includes two runable programs, nw_algorithm.py and convert_to_amino.py\
nw_algorithm.py is the needleman-wunsch algorithm used in bioinformatics to align protein or nucleotide sequences\
it can be run in windows command line as python nw_algorithm.py [query file] [ref file] [outputfile]\
Some things to note:\
the program expects that the sequences within query file and ref file are all on the second line of the file\
this is how they are all structured within the folder, so they should all run properly.\
Additionally, if you want to run nw_algorithm.py with init/terminal gaps not counted then the third parameter of 'dynamic_alignment(query,ref, InitPenalty)' should be set to False.\
This is all explained within the program itself as well\
the place to change this would be in main at the bottom of the program, where dynamic alignment is called.\

Originally, this paramter has been set to True so that init/terminal gap penalties ARE counted.\
