import os
import csv
import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
def fasta_creator(data_input, fasta_output_path):
    """
    This function will create a FASTA file from the input CSV file. For this it will use the **Biopyton** module.

    :param path_input: pandas data frame
    :type path_input: pandas data frame

    :param fasta_output_path: Path to where we want to save the FASTA file.
    :type fasta_output_path: string

    :return: All the data from the CSV in a FASTA format.
    :rtype: FASTA File
    """
    matrix = []
    for index, (_, sequence) in enumerate(data_input.iterrows()):
        # index += 1 # To start the index in 1
        rec = SeqRecord(Seq(sequence.loc["sseq"]),  # In the 5 position is the seq
                        id="Seq_" + str(index),
                        description="Leishmania infantum"  # argument maybe?
                        )
        matrix.append(rec)
    SeqIO.write(matrix, fasta_output_path, "fasta")

def columns_to_numeric(data_input, columns_to_convert):
    """
    This function will convert the columns of a pandas DataFrame to numeric values.

    :param data_input: pandas data frame
    :type data_input: pandas data frame

    :return: pandas data frame with numeric values
    :rtype: pandas data frame
    """
    for column in columns_to_convert:
        data_input[column] = pd.to_numeric(data_input[column], errors='coerce')
    return data_input