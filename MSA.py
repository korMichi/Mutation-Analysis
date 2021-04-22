import pandas as pd
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


if __name__ == "__main__":
    seq_df = pd.read_excel("", index_col=0)
    get_v_gen_groups(seq_df, "", "")
