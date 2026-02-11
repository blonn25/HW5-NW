# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    # TODO Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    sub_matrix_file = "./substitution_matrices/BLOSUM62.mat"
    gap_open = -10
    gap_extend = -1
    nw_algo = NeedlemanWunsch(sub_matrix_file, gap_open, gap_extend)
    alignments = []
    for seq, header in [(gg_seq, gg_header), (mm_seq, mm_header), (br_seq, br_header), (tt_seq, tt_header)]:
        
        # get the species frm the header
        species = header.split('OS=')[1].split('OX=')[0].strip()
        
        # perform the alignment and add species, score, and aligned sequences to alignments list
        alignment_score, seqA_align, seqB_align = nw_algo.align(hs_seq, seq)
        alignments.append([species, alignment_score, seqA_align, seqB_align])

    # sort the alignments by the score in descending order (higher score is more similar) and print the species and score
    alignments.sort(key=lambda x: x[1], reverse=True)
    for species, _, _, _ in alignments:
        print(species)

    # TODO print all of the alignment score between each species BRD2 and human BRD2
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    for _, alignment_score, _, _ in alignments:
        print(alignment_score)
    

if __name__ == "__main__":
    main()
