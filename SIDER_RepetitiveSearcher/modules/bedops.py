# Modules needed
import numpy as np
import pandas as pd
import subprocess
import os

from modules.files_manager import columns_to_numeric

def get_data_sequence(data, strand, genome_fasta):
    """
    This function gets the sequence of the data from the fasta file. It will keep the Chromosome ID, start coordinate, end coordinate and strand.

    :param data: A pandas data frame with the data read of the BED files.
    :type data: pandas.core.frame.DataFrame

    :param strand: The strand of the sequence. It can be "plus" or "minus".
    :type strand: str

    :param genome_fasta: Path to the whole genome sequence in FASTA format.
    :type genome_fasta: string
    """
    sequences = []
    for _, row in data.iterrows():
        sseqid = row["sseqid"]
        start = row["sstart"]
        end = row["send"]
        cmd = [
            "blastdbcmd",
            "-db", genome_fasta,
            "-entry", sseqid,
            "-range", f"{start}-{end}",
            "-strand", strand,
            "-outfmt", "%s"
        ]

        sequence = subprocess.check_output(cmd, universal_newlines=True).replace('\n', '')

        sequences.append({
            "sseqid": sseqid,
            "sstart": start,
            "send": end,
            "sstrand": strand,
            "sseq": sequence
        })

    sequences_df = pd.DataFrame(sequences)
    return sequences_df

def bedops_main(data_input, genome_fasta, writing_path_input):
    """
    This function will implement BEDOPS to filter duplicates and overlaps in a CSV file.

    :param path_input: Path to the CSV file we want to filter data.
    :type path_input: string

    :param genome_fasta: Path to the whole genome sequence in FASTA format.
    :type genome_fasta: string
    """
    # -----------------------------------------------------------------------------
    # 1) Filter and sort data
    # -----------------------------------------------------------------------------
    columns_ids = data_input.columns  # gets the columns names
    df_plus = data_input[data_input["sstrand"] == "plus"]  # filters the "+" strand
    df_minus = data_input[data_input["sstrand"] == "minus"]  # filters the "-" strand

    # Sort the data by the start coordinate
    df_plus = df_plus.sort_values(by=["sseqid", "sstart"])  # sorts the "+" strand by the start coordinate
    df_minus = df_minus.sort_values(by=["sseqid", "sstart"])  # sorts the "-" strand by the start coordinate

    # -----------------------------------------------------------------------------
    # 2) BEDOPS files creation:
    # -----------------------------------------------------------------------------
    #  row[1] == Chromosome ID, row[10] == Start coordinate, row[11] == End coordinate
    plus_path = os.path.join(writing_path_input, os.path.basename(writing_path_input) + "_plus.bed")
    minus_path = os.path.join(writing_path_input, os.path.basename(writing_path_input) + "_minus.bed")
    
    df_plus[["sseqid", "sstart", "send"]].to_csv(plus_path, sep="\t", header=False, index=False)  # creates a BED file for the "+" strand
    # In the minus strand its important to change the order of the coordinates, because "bedops" reads them like "plus" strand.
    # If not, it will not merge them.
    df_minus[["sseqid", "send", "sstart"]].to_csv(minus_path, sep="\t", header=False, index=False)  # creates a BED file for the "-" strand

    # -----------------------------------------------------------------------------
    # 3) BEDOPS function call with subprocess.
    # -----------------------------------------------------------------------------
    # BEDOPS call to "plus.bed" and "minus.bed" files
    # Using subprocess to call BEDOPS.
    # Using subprocess.check_output() to get the output from the command.
    # shell=True to have the input of .check_output() as a string.
    # universal_newlines=True to have the EoL character as "\n" (Unix-like systems).
    # The output will be a variable of strings.
    df_plus_bedops = subprocess.check_output(f"bedops --merge {plus_path}", shell=True, universal_newlines=True)  # merges the "+" strand BED file
    df_minus_bedops = subprocess.check_output(f"bedops --merge {minus_path}", shell=True, universal_newlines=True)  # merges the "-" strand BED file

     # Now let's transform then into Data Frames
    df_plus_bedops = pd.DataFrame([x.split("\t") for x in df_plus_bedops.split("\n") if x],
                                    columns=["sseqid", "sstart", "send"])  # transforms the "+" strand BEDOPS output into a Data Frame
    df_minus_bedops = pd.DataFrame([x.split("\t") for x in df_minus_bedops.split("\n") if x],
                                    columns=["sseqid", "sstart", "send"])  # transforms the "-" strand BEDOPS output into a Data Frame
    # -----------------------------------------------------------------------------
    # 4) Call `blastdbcmd` to get the sequences with the function get_data_sequence()
    # -----------------------------------------------------------------------------
    if df_plus_bedops.empty:  # In case the original data is empty, the code needs to keep going
        df_plus_bedops_wseq = pd.DataFrame(columns=columns_ids)  # creates an empty Data Frame with 5 columns
    else:  # If the original data is not empty, tit uses get_data_sequence
        df_plus_bedops_wseq = get_data_sequence(df_plus_bedops, "plus", genome_fasta)

    # The same for the minus strand:
    if df_minus_bedops.empty:
        df_minus_bedops_wseq = pd.DataFrame(columns=columns_ids)
    else:   
        df_minus_bedops_wseq = get_data_sequence(df_minus_bedops, "minus", genome_fasta)

 
    # Let's reorderthe `df_minus_bedps_wseq` data frame:
    df_minus_bedops_wseq[["sstart", "send"]] = df_minus_bedops_wseq[["send", "sstart"]].copy()  # swap only values
    # -----------------------------------------------------------------------------
    # 5) Processing data
    # -----------------------------------------------------------------------------
    # Join both data frames
    all_data = pd.concat([df_plus_bedops_wseq, df_minus_bedops_wseq], ignore_index=True)  # joins both Data Frames

    # Adding sequence length to the DataFrame:
    new_column = [len(x) for x in all_data.loc[:,"sseq"]]  # creates a list with the length of each sequence
    all_data.insert(1, "length", new_column)  # inserts the new column with the sequence length. Column index are shifted.


    # -----------------------------------------------------------------------------
    # 6) Correctly modeling the output Data Frame to 15 columns and output as CSV file.
    # -----------------------------------------------------------------------------
    new_data = pd.DataFrame(index=range(all_data.shape[0]), columns=columns_ids)  # creates a new Data Frame with 15 columns. The rows depends on the .shape[0]

    new_data.loc[:,["sseqid", "length", "sstart", "send", "sstrand", "sseq"]] = all_data.loc[:,["sseqid", "length", "sstart", "send", "sstrand", "sseq"]].copy()
    new_data = columns_to_numeric(new_data, ["pident", "length", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen"])
    
    return new_data  # returns the new Data Frame


def bedops_coincidence(last_df, old_df, folder_path, strand, genome_fasta):
    """
    Will tell the elements from old_df that are in last_df.
    last_df is the last data frame that we have.
    old_df is the first input data frame.
    """
    last_df = last_df.sort_values(by=["sseqid", "sstart"])  # Sort the data frame by the start coordinate
    old_df = old_df.sort_values(by=["sseqid", "sstart"])  # Sort the data frame by the start coordinate
    last_length = last_df.shape[0]
    old_length = old_df.shape[0]

    last_df_path = os.path.join(f"{folder_path}_1.1_last.bed")
    old_df_path = os.path.join(f"{folder_path}_1.2_old.bed")

    last_df[["sseqid", "sstart","send"]].to_csv(last_df_path, sep="\t", header=False, index=False)
    old_df[["sseqid", "sstart","send"]].to_csv(old_df_path, sep="\t", header=False, index=False)

    # -----------------------------------------------------------------------------
    # Call BEDOPS to find coincidences
    cmd = f"bedops --element-of 10 {last_df_path} {old_df_path}"
    Last_in_Old = subprocess.check_output(cmd, shell=True, universal_newlines=True)
    Last_in_Old = pd.DataFrame([x.split("\t") for x in Last_in_Old.split("\n") if x],
                             columns=["sseqid", "sstart", "send"])
    Last_in_Old = columns_to_numeric(Last_in_Old, ["sstart", "send"])
    print(f"\t\t\t- Last in old: {Last_in_Old.shape[0]}/{last_length} - {Last_in_Old.shape[0]/last_length*100:.2f}%")

    cmd = f"bedops --element-of 10 {old_df_path} {last_df_path}"
    Old_in_Last = subprocess.check_output(cmd, shell=True, universal_newlines=True)
    Old_in_Last = pd.DataFrame([x.split("\t") for x in Old_in_Last.split("\n") if x],
                             columns=["sseqid", "sstart", "send"])
    Old_in_Last = columns_to_numeric(Old_in_Last, ["sstart", "send"])
    print(f"\t\t\t- Old in last: {Old_in_Last.shape[0]}/{old_length} - {Old_in_Last.shape[0]/old_length*100:.2f}%")
    # -----------------------------------------------------------------------------
    # Let's merge
    Old_in_Last_path = os.path.join(f"{folder_path}_2.1_Old_in_Last.bed")
    Last_in_Old_path = os.path.join(f"{folder_path}_2.2_Last_in_Old.bed")

    Old_in_Last[["sseqid", "sstart", "send"]].to_csv(Old_in_Last_path, sep="\t", header=False, index=False)
    Last_in_Old[["sseqid", "sstart", "send"]].to_csv(Last_in_Old_path, sep="\t", header=False, index=False)

    cmd = f"bedops --merge {Old_in_Last_path} {Last_in_Old_path}"
    merged_LastOld = subprocess.check_output(cmd, shell=True, universal_newlines=True)
    merged_LastOld = pd.DataFrame([x.split("\t") for x in merged_LastOld.split("\n") if x],
                                   columns=["sseqid", "sstart", "send"])
    merged_LastOld = columns_to_numeric(merged_LastOld, ["sstart", "send"])

    # Now recapture the elements with the genome
    coincidence_data = get_data_sequence(merged_LastOld, strand, genome_fasta)
    # -----------------------------------------------------------------------------
    # -----------------------------------------------------------------------------
    # Now let's check the elements that are not in df2 (the first input). They would be the new elements.
    cmd = f"bedops --not-element-of 10 {last_df_path} {old_df_path}"
    Last_notin_Old = subprocess.check_output(cmd, shell=True, universal_newlines=True)
    Last_notin_Old = pd.DataFrame([x.split("\t") for x in Last_notin_Old.split("\n") if x],
                             columns=["sseqid", "sstart", "send"])
    Last_notin_Old = columns_to_numeric(Last_notin_Old, ["sstart", "send"])
    print(f"\t\t\t- Last NOT in old: {Last_notin_Old.shape[0]}/{last_length} - {Last_notin_Old.shape[0]/last_length*100:.2f}%")

    if not Last_notin_Old.empty:  # If the data frame is not empty
        new_data = get_data_sequence(Last_notin_Old, strand, genome_fasta)
    else:  # If the data frame is empty
        new_data = pd.DataFrame()
    # -----------------------------------------------------------------------------
    # Now check the elements in Old that are not in Last
    cmd = f"bedops --not-element-of 10 {old_df_path} {last_df_path}"
    Old_notin_Last = subprocess.check_output(cmd, shell=True, universal_newlines=True)
    Old_notin_Last = pd.DataFrame([x.split("\t") for x in Old_notin_Last.split("\n") if x],
                             columns=["sseqid", "sstart", "send"])
    Old_notin_Last = columns_to_numeric(Old_notin_Last, ["sstart", "send"])
    print(f"\t\t\t- Old NOT in last: {Old_notin_Last.shape[0]}/{old_length} - {Old_notin_Last.shape[0]/old_length*100:.2f}%")

    if not Old_notin_Last.empty:  # If the data frame is not empty
        old_data_exclusive = get_data_sequence(Old_notin_Last, strand, genome_fasta)
    else:  # If the data frame is empty
        old_data_exclusive = pd.DataFrame()

    # -----------------------------------------------------------------------------
    return coincidence_data, new_data, old_data_exclusive

def bedops_stopping(new_data, old_data, folder_path, query_cover):
    new_data_length = new_data.shape[0]
    old_data_length = old_data.shape[0]

    # -----------------------------------------------------------------------------
    # prepare paths
    new_data_path = os.path.join(folder_path, "new_data.bed")
    old_data_path = os.path.join(folder_path, "old_data.bed")
    # -----------------------------------------------------------------------------
    # Save the data in BED format
    new_data[["sseqid", "sstart", "send"]].to_csv(new_data_path, sep="\t", header=False, index=False)
    old_data[["sseqid", "sstart", "send"]].to_csv(old_data_path, sep="\t", header=False, index=False)
    # -----------------------------------------------------------------------------
    # Call BEDOPS to find coincidences
    cmd = f"bedops --element-of {query_cover}% {new_data_path} {old_data_path}"
    new_in_old = subprocess.check_output(cmd, shell=True, universal_newlines=True)
    new_in_old = pd.DataFrame([x.split("\t") for x in new_in_old.split("\n") if x],
                              columns=["sseqid", "sstart", "send"])
    new_in_old = columns_to_numeric(new_in_old, ["sstart", "send"])
    print(f"\t\t- QC {query_cover}% vs old: {new_in_old.shape[0]}/{new_data_length} - {new_in_old.shape[0]/new_data_length*100:.2f}%")

    if new_in_old.shape[0]/new_data_length*100 >= query_cover:
        return True
    else:
        return False
    # -----------------------------------------------------------------------------

