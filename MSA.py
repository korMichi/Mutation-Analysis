import pandas as pd
import os
import subprocess
from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline


def write_sequences(file_name, sequences, column_sequence_name):
    """Opens file with the group name and writes all corresponding sequences into it."""
    f = open("{0}.fasta".format(file_name), "w")
    for n, i in enumerate(sequences[column_sequence_name]):
        fasta_string = ">" + str(sequences[column_sequence_name].index[n]) + "\n" + str(sequences[column_sequence_name][n]) + "\n"
        f.write(fasta_string)
    f.close()


def get_v_gen_groups(df, germline_genes, column_v_gen_name, column_sequence_name, switch="forgot"):
    v_gen_list = df[column_v_gen_name].unique()
    reduced_df = df[[column_v_gen_name, column_sequence_name]].copy().rename(columns={"V-Gene": "gene", column_sequence_name: "sequence"})
    for v_gen in v_gen_list:
        if switch == "yes":
            for_groups = df.loc[df["V-Gen-Group"] == v_gen, "V-Gene"].unique() # must be commented out if germline genes not included
            result_frame = pd.concat([reduced_df.loc[df[column_v_gen_name] == v_gen], germline_genes.loc[germline_genes["gene"] == (for_groups[0] + "*01")]]) #for_groups[0] if with germline gene or str(v_gen) if without
            write_sequences(v_gen, result_frame, "sequence") 
        elif switch == "no":
            result_frame = reduced_df.loc[df[column_v_gen_name] == v_gen]
            write_sequences(v_gen, result_frame, "sequence") 
        elif switch == "forgot":
            print("Please add if you would like to include the V-Gen!")
            break


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
#    seq_df = seq_df.loc[~seq_df.index.str.contains("HC")] # if only heavy chain or light chain is used in the grouped analysis
    germline_genes = pd.read_excel("Germline_Genes_AA.xlsx") # Read in database file of aminoacid sequences from germline genes
    germline_genes.index = germline_genes["gene"]

    get_v_gen_groups(seq_df, germline_genes, "", "") # Include if you would like to include the germline genes as a sequence into the file
    path_to_files = ""
    muscle_alignment(path_to_files)
