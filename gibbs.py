# Gibb's Sampler
# author: Jeffrey Lin
# Date: 2/24/2023
#
# Purpose: Detects local motifs found within a set of multiple sequences

import numpy as np 
import sys as sys
import math

# parse_FASTA_file
# purpose: Finds and returns a list of sequences from the FASTA file 
# input: FASTA file containing the sequences
# output: returns the two sequences as strings x and y respectively
# Effects: none

def parse_FASTA_file(file):
    Sequence_list = []
    index = -1

    #reads in sequences into strings
    for line in file:
        if ">" in line:
            Sequence_list.append("")
            index += 1
            continue
        else:
            Sequence_list[index] += line.rstrip()

    return Sequence_list

class SimpleGibbsSampler:

    # Constructor
    # Purpose: initializes SimpleGibbsSampler class
    # Input: list of seqs, length of motif, the language, language probability,
    #        seed for random generator 
    # Output: none
    # Effects: creates instance of SimpleGibbsSample class

    def __init__(self, seqs, motif_len, lang, lang_prob, seed):
        self.seqs = seqs
        self.motif_len = motif_len
        self.lang = lang
        self.lang_prob = lang_prob
        self.seed = seed

        # sets up randoms state generator
        self.random_state = np.random.RandomState(seed)

        self.init_positions = self.pick_init_positions()

        # list holding best position for each withheld sequence 
        self.withheld_list = [None for x in range(len(self.seqs)) ]

    # run_sampler
    # Purpose: Executes the SimpleGibbsSampler program
    # Input: none
    # Output: none
    # Effects: None

    def run_sampler(self):
        
        # sets starting values so that simulation can begin
        Rem_seq = 0
        prev_rem = -1
        iteration = 1
        best_score_pos = -1
                
        while (self.__continue(best_score_pos, Rem_seq)):
            best_score_prev = best_score_pos
            
            motifs, Rem_seq = self.__find_motifs(self.init_positions, prev_rem)
            count_matrix = self.build_pseudo_count_matrix(motifs)
            pssm = self.build_pssm(count_matrix)
            
            scores_list = self.score_seq_windows(self.seqs[Rem_seq], pssm)
            best_score_pos = self.__best_scoring_pos(scores_list, Rem_seq)
            self.init_positions[Rem_seq] = best_score_pos

            msa = self.get_msa()


            print("Iteration: " + str(iteration))
            print("S* index: " + str(Rem_seq))
            print("S*: " + self.seqs[Rem_seq])
            print("pssm:")
            print(pssm)
            self.__print_msa(msa)

            prev_rem = Rem_seq
            iteration += 1
        
        return msa

    # pick_init_positions
    # Purpose: Determines initial position of motifs in sequences
    # Input: none
    # Output: none
    # Effects: None
    def pick_init_positions(self):
        init_pos_list = []

        for i in range(0, len(self.seqs)):
            upper_limit = len(self.seqs[i]) - self.motif_len
            init_pos = self.random_state.randint(0, upper_limit )

            init_pos_list.append(init_pos) 

        return init_pos_list
    
    # find motifs
    # Purpose: Determines motifs corresponding to the starting positions
    # Input: list of initial positons. previous removed sequence
    # Output: returns list of strings representing each of the motifs
    # Effects: None

    def __find_motifs(self, init_pos_list, prev_removed):
        # determines new removed sequence for iteration
        Removed_sequence = self.random_state.randint(0, len(self.seqs))
        
        # Reselects removed sequence if it is the same as prev iteration
        if (Removed_sequence == prev_removed):
            while (Removed_sequence == prev_removed):
                Removed_sequence = self.random_state.randint(0, len(self.seqs))

        motifs = []

        # Finds the subsequences of length motif_len from the non-removed seqs
        for i in range(0, len(init_pos_list)):
            if (i == Removed_sequence):
                continue
            else:
                start = init_pos_list[i]
                stop = start + self.motif_len
                motif = self.seqs[i][start:stop]

                motifs.append(motif)
        
        return motifs, Removed_sequence 



    # build_pseudo_count_matrix
    # Purpose: builds pseudo count matrix
    # Input: list of strings representing all the motifs
    # Output: returns the pseudo count matrix
    # Effects: creates a pseudo count matrix
    
    def build_pseudo_count_matrix(self, motifs):
        dimY = len(self.lang)
        dimX = self.motif_len

        count_matrix = np.full([dimY, dimX], 1, dtype = float)

        for i in range(0, self.motif_len):
            self.build_pseudo_count_helper(count_matrix, i, motifs)
        
        return count_matrix

    
    # build_pseudo_count_helper
    # Purpose: Fills in column of count matrix according to freq of each base
    # Input: count_matrix, index, motifs
    # Output: none
    # Effects: Fills pseudo count matrix

    def build_pseudo_count_helper(self, count_matrix, index, motifs):
        # lists of the frequency for each char in the language
        lang_count = [0 for x in range(len(self.lang))]

        # Determines counts of each char in language for an index of the motifs
        for i in range(0, len(self.lang)):
            for j in range(0, len(motifs)):
                if self.lang[i] == motifs[j][index]:
                    lang_count[i] += 1

        # fills the count_matrix according to counts above
        for k in range(0, len(self.lang)):
            count_matrix[k][index] += lang_count[k]
    
    # build_pssm
    # Purpose: creates the pssm matrix
    # Input: count matrix
    # Output: returns pssm matrix
    # Effects: creates pssm matrix
    
    def build_pssm(self, count_matrix):
        dimY = len(self.lang)
        dimX = self.motif_len

        pssm = np.empty([dimY, dimX], dtype = float)

        # determines sum of a column in count matrix
        col_sum = np.sum(count_matrix, axis = 0)
        col_sum = col_sum[0]

        # background probability
        backgr = self.lang_prob

        for i in range(0, len(self.lang)):
            for j in range(0, self.motif_len):
                pssm[i][j] = math.log2((count_matrix[i][j] / col_sum) / backgr)
        
        return pssm

    # score_seq_windows
    # Purpose: scores every window of the withheld sequence
    # Input: seqs representing the removed sequence and the pssm
    # Output: returns a list corresponding to scores of each window
    # Effects: creates scores list

    def score_seq_windows(self, seqs, pssm):
        scores_list = []

        for i in range(len(seqs) - self.motif_len + 1):
            scores_list.append(self.__score_seq_windows_helper(seqs, pssm, i))

        return scores_list

    # scores_seq_windows_helper
    # Purpose: calculates scores of the current window
    # Input: seqs representing the removed sequence, pssm, and index for 
    #        current window
    # Output: returns score of the window
    # Effects: none

    def __score_seq_windows_helper(self, seqs, pssm, index):
        log_odds_score = 0

        for i in range(self.motif_len):
            if (seqs[index + i] == 'A'):
                log_odds_score += pssm[0][i]
            elif (seqs[index + i] == 'C'):
                log_odds_score += pssm[1][i]
            elif (seqs[index + i] == 'G'):
                log_odds_score += pssm[2][i]
            elif (seqs[index + i] == 'T'): 
                log_odds_score += pssm[3][i]
            
        return log_odds_score

    # get_msa
    # Purpose: determines the current msa
    # Input: none
    # Output: returns the msa
    # Effects: none

    def get_msa(self):
        msa = []

        for i in range(0, len(self.init_positions)):
            start = self.init_positions[i]
            stop = start + self.motif_len
            motif = self.seqs[i][start:stop]
            msa.append(motif)

        return msa
    
    # get_msa
    # Purpose: prints the current MSA
    # Input: list representing the msa
    # Output: none
    # Effects: none

    def __print_msa(self, msa):
        print("Current MSA: ")

        for i in range(len(msa)):
            print(msa[i])
        
        print("\n")

    # best_scoring_pos
    # Purpose: returns index of highest scoring position in withheld sequence
    # Input: score_list, removed_sequence
    # Output: index of highest scoring position  in withheld sequence
    # Effects: none
    def __best_scoring_pos(self, score_list, removed_sequence):

        best_score = np.argmax(score_list)
        
        return best_score

    # continue
    # Purpose: Determines whether the program should continue running based on
    #          whether any of the highest scoring positions in a seq converge
    # Input: the current best scoring position and the index of the removed 
    #        sequence
    # Output: boolean for whether to continue or not
    # Effects: determines whether program will quit
    def __continue(self, best_scoring_pos, removed_seq):
        if (best_scoring_pos != self.withheld_list[removed_seq]):
            self.withheld_list[removed_seq] = best_scoring_pos
            return True
        elif (best_scoring_pos == self.withheld_list[removed_seq]):
            return False
    

def main():
    lang = ['A', 'C', 'G', 'T']

    seqs = parse_FASTA_file(sys.stdin)

    sgs = SimpleGibbsSampler(seqs, 6, lang, 0.25, 39) 
    msa = sgs.run_sampler()

if __name__ == "__main__":
    main()
