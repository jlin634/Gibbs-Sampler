# Gibbs-Sampler

## Program Purpose:  
This program detects local motifs found within a set of multiple sequences

## Directions on Compilation of the Program
compilation should follow the following format:    
    python3 gibbs.py < [insert FASTA file] > output
    
Use the following compilation instructions instead, if compiling in powershell: \n
Get-Content [insert FASTA file] | python3 gibbs.py > output
    
In this case, we are piping a FASTA-formatted file into stdin and sending
all output to some output.txt file. Sending to a separate output file is 
an optional component. If the client chooses not to send to an output text
file, then all output will just print to terminal. The default values are 
as follows: motif_len = 6, lang = [’A’, ’C’, ’G’, ’T’], 
lang_prob = 0.25, and seed = 20. No command line interaction has been 
implented.
