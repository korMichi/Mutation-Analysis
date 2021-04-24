import pandas as pd
import os
import subprocess
from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline


def write_alignment(file_name, sequences, column_sequence_name):
    """Opens file with the group name and writes all corresponding sequences into it."""
    f = open("{0}.fasta".format(file_name), "w")
    for n, i in enumerate(sequences[column_sequence_name]):
        fasta_string = ">" + str(sequences[column_sequence_name].index[n]) + "\n" + str(sequences[column_sequence_name][n]) + "\n"
        f.write(fasta_string)
    f.close()


def get_v_gen_groups(df, column_v_gen_name, column_sequence_name):
    """Identifies groups of sequences (e.g. v-gen) depending on column name and passes dataframe of groups to write_alignment()
       Inputs: dataframe, name_of_columns_with_group, name_of_column_with_sequence"""
    v_gen_list = df[column_v_gen_name].unique()
    for v_gen in v_gen_list:
         write_alignment(v_gen, df.loc[df[column_v_gen_name] == v_gen], column_sequence_name)


def muscle_alignment(path):
    """Performs MUSCLE alignment using the command line tool and writes to output file.
       CAVE: filenames can not include special characters such as (*, /, &)"""
    for files in os.listdir(path):
        if "fasta" in files:
            file_path = path + "/{0}".format(files)
            output_file = str(files)
            muscle_cline = MuscleCommandline(input=file_path, out=output_file)
            child = subprocess.Popen(str(muscle_cline),
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              universal_newlines=True,
                              shell=True)
            child.wait()
            with open(output_file) as align_handle:
                align = AlignIO.read(align_handle, "fasta")


if __name__ == "__main__":
    seq_df = pd.read_excel("", index_col=0)
    get_v_gen_groups(seq_df, "", "")
    path = ""
    muscle_alignment(path)
