# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    # get seqs
    seq1, _ = read_fasta("./data/test_seq1.fa") # MYQR
    seq2, _ = read_fasta("./data/test_seq2.fa") # MQR

    # initialize the algorithm and
    sub_matrix_file = "./substitution_matrices/BLOSUM62.mat"
    gap_open = -10
    gap_extend = -1
    nw_algo = NeedlemanWunsch(sub_matrix_file, gap_open, gap_extend)

    # perform the alignment
    nw_algo.align(seq1, seq2)

    # initialize matrices solved by hand for comparison
    #                                      M      Q      R
    my_align_matrix = np.array([[  0.0, -11.0, -12.0, -13.0],
                                [-11.0,   5.0,  -6.0,  -7.0],  # M
                                [-12.0,  -6.0,   4.0,  -7.0],  # Y
                                [-13.0,  -7.0,  -1.0,   5.0],  # Q
                                [-14.0,  -8.0,  -6.0,   4.0]]) # R
    #                                        M        Q        R
    my_gapA_matrix = np.array([[-np.inf, -np.inf, -np.inf, -np.inf],
                               [-np.inf,   -22.0,   -23.0,   -24.0],  # M
                               [-np.inf,    -6.0,   -17.0,   -18.0],  # Y
                               [-np.inf,    -7.0,    -7.0,   -18.0],  # Q
                               [-np.inf,    -8.0,    -8.0,    -6.0]]) # R
    #                                        M        Q        R
    my_gapB_matrix = np.array([[-np.inf, -np.inf, -np.inf, -np.inf],
                               [-np.inf,   -22.0,    -6.0,    -7.0],  # M
                               [-np.inf,   -23.0,   -17.0,    -7.0],  # Y
                               [-np.inf,   -24.0,   -18.0,   -12.0],  # Q
                               [-np.inf,   -25.0,   -19.0,   -17.0]]) # R

    # assert that the alignment matrices are correctly filled out
    assert nw_algo._align_matrix == pytest.approx(my_align_matrix)
    assert nw_algo._gapA_matrix == pytest.approx(my_gapA_matrix)
    assert nw_algo._gapB_matrix == pytest.approx(my_gapB_matrix)

def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    # get seqs
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")

    # initialize the algorithm
    sub_matrix_file = "./substitution_matrices/BLOSUM62.mat"
    gap_open = -10
    gap_extend = -1
    nw_algo = NeedlemanWunsch(sub_matrix_file, gap_open, gap_extend)

    # perform the alignment
    alignment_score, seq3_align, seq4_align = nw_algo.align(seq3, seq4)

    # assert that the score and aligned sequences produced by the backtrace are correct
    assert alignment_score == 17
    assert seq3_align == "MAVHQLIRRP"
    assert seq4_align == "M---QLIRHP"




