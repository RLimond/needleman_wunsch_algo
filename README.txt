This folder includes two runable programs, nw_algorithm.py and convert_to_amino.py
nw_algorithm.py is the needleman-wunsch algorithm we were tasked to implement in question 1.
it can be run in windows command line as nw_algorithm.py [query file] [ref file] [outputfile]
Some things to note:
the program expects that the sequences within query file and ref file are all on the second line of the file
this is how they are all structured within the folder, so they should all run properly.
Additionally, if you want to run nw_algorithm.py with init/terminal gaps not counted
the third parameter of 'dynamic_alignment(query,ref, InitPenalty)' should be set to False. This is all explained within the program itself as well
the place to change this would be in main at the bottom of the program, where dynamic alignment is called.

Originally, this paramter has been set to True so that init/terminal gap penalties ARE counted.