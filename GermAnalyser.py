import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from Bio import Align, SeqIO
from Bio.Align import substitution_matrices
from Bio.Seq import Seq

### Initialize Parameters for Pairwise Alignment same as Needleman-Wunsch 
aligner = Align.PairwiseAligner()
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")

### Define Functions

def mutation_count(alignment):
    """Takes Biopython Alignment Object in a pandas Dataframe, splits into strings and counts mismatches"""
    count = 0
    wt_seq = format(alignment[0]).split("\n")[0]
    gl_seq = format(alignment[0]).split("\n")[2]
    for element, value in enumerate(wt_seq):
        if value != gl_seq[element]:
            count += 1
    return count


def mutation_identifier(alignment):
    """Takes Biopython Alignment Object in a pandas Dataframe, splits into strings and returns all mismatches"""
    mismatches = []
    wt_seq = format(alignment[0]).split("\n")[0]
    gl_seq = format(alignment[0]).split("\n")[2]
    for element, value in enumerate(wt_seq):
        if value != gl_seq[element]:
            mismatches.append("{0}: {1} -> {2}".format((element+1), gl_seq[element], value))
    if len(mismatches) == 0:
        pass
    else:
        return mismatches


def pd_aligner(df, seq1, seq2, name):
    """Function that creates alignment object in pandas dataframe using Biopython Align
        Input: Pandas Dataframe, Columnname Seq1, Columnname Seq2, Columnname of alignment
        Output: Alignment object stored in a new dataframe column"""
    for index in df.index:
        df.loc[index, name] = aligner.align(df.loc[index, seq1], df.loc[index, seq2])
    return df



def new_mutation_count(alignment):
    """Takes Biopython Alignment Object in a pandas Dataframe, splits into strings and counts mismatches
        Function version 1.1"""
    count = 0 
    wt_seq = format(alignment).split("\n")[0]
    gl_seq = format(alignment).split("\n")[2]
    for element, value in enumerate(wt_seq):
        if value != gl_seq[element]:
            count += 1
    return count


def pd_aligner(df, seq1, seq2, name):
    """Function that creates alignment object in pandas dataframe using Biopython Align
        Input: Pandas Dataframe, Columnname Seq1, Columnname Seq2, Columnname of alignment
        Output: Alignment object stored in a new dataframe column"""
    number_of_alignments = []
    for index in df.index:
        alignment = aligner.align(df.loc[index, seq1], df.loc[index, seq2])
        df.loc[index, name] = alignment[0]
    return df