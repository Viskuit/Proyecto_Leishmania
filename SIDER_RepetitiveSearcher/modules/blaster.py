import os
import pandas as pd
import subprocess
import time
import shutil
from pathlib import Path
from datetime import datetime

from modules.aesthetics import boxymcboxface  # Some aesthetics function
from modules.identifiers import genome_specific_chromosome_main
from modules.filters import global_filters_main
from modules.files_manager import columns_to_numeric
from modules.compare import compare_main

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


def blastn_dic(path_input, path_output):
    """
    Creation af a BLAST database of our whole genome. It uses the BLAST :sup:`R` command line, see BLAST
    `Command Line Application User Manual`_ for more information.


    The generation of the properly database will be placed in the directory where ``path_input`` is.
    It is recommended to use a dedicated folder to this FASTA file so the database is written next to it.

    :param path_input: path to a FASTA file.
    :type path_input: string

    :param path_output: path to the output folder where the BLAST database will be created.
    :type path_output: string

    :return: a BLAST database.
    :rtype: Muitiples files (**.nhr**, **.nin**, **.nog**, **.nsd**, **.nsi** and **.nsq** extensions)
    """

    # Remember is "path.input.dic_path" for "argparse".
    try:
        # "parse_seqids" is used to keep the sequence ID in the output.
        cmd = f"makeblastdb -in {path_input} -dbtype nucl -parse_seqids -out {path_output}"
        subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except Exception:
        print("\nError: Blast Dictionary couldn't be created")

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


def blastn_blaster(query_path, dict_path, perc_identity):
    """
    This module calls for `blastn` in the BLAST :sup:`R` command line.
    See `Command Line Application User Manual`_ for more information.

    :param query_path: Path to our FASTA query used in BLASTn against our database.
    :type query_path: string

    :param dict_path: Path to the FASTA file ``query_path`` will be launched to. In the same directory should be the BLAST data base created with :func:`~blastn_dic`.
    :type dict_path: string

    :param outfile_path: Path where the results will be saved. Name the whole file, not only the path.
    :type outfile_path: string

    :param perc_identity: Percent of sequence identity which we want to make the BLASTn. **Important**.
    :type perc_identity: string

    :return: A CSV file with the results of the BLASTn.
    :rtype: CSV file


    .. list-table:: Arguments used in ``blastn``
       :header-rows: 1

       * - Argument
         - Type
         - Description

       * - ``-word_size``
         - integer
         - | 28 by default.
           | BLAST search starts with finding a perfect sequence match of length given by this parameter.
           | See `world_size info`_ for more information.

       * - ``-query``
         - string
         - Query file name.

       * - ``-db``
         - string
         - BLAST database name.

       * - ``-out``
         - string
         - Output file name.

       * - ``-perc_identity``
         - integer
         - Percent identity cutoff.

       * - ``-outfmt``
         - string
         - | Alignment view options. Written inside in a nested string.
           | Explained in :ref:`Table <outfmt_Table>`
           | See `Command Line Application User Manual`_ and `Metagenomics BLASTn manual`_ for help.


    ..
        # This is to mark the next table

    .. _outfmt_Table:


    .. list-table:: Arguments used in ``--outfmt`` inside ``blastn``
       :header-rows: 1


       * - Argument
         - Description

       * - ``10``
         - Comma-separated values.

       * - ``qseqid``
         - | Query Seq-id.
           | ID of our query sequence.

       * - ``sseqid``
         - | Subject Seq-id.
           | ID of our subject sequence.

       * - ``pident``
         - Percentage of identical matches.

       * - ``length``
         - Alignment length.

       * - ``qlen``
         - Query sequence length.

       * - ``slen``
         - Subject sequence length.

       * - ``mismatch``
         - Number of mismatches.

       * - ``gapopen``
         - Number of gap openings.

       * - ``qstart``
         - Start of alignment in query.

       * - ``qend``
         - End of alignment in query.

       * - ``sstart``
         - Start of alignment in subject.

       * - ``send``
         - End of alignment in subject.

       * - ``evalue``
         - Expect value.

       * - ``bitscore``
         - Bit score.

       * - ``sstrand``
         - Subject strand.

       * - ``sseq``
         - Aligned part of subject sequence.
    """

    cmd = "blastn -word_size 15 -query " \
        + query_path + " -db " \
        + dict_path \
        + " -perc_identity " + str(perc_identity) \
        + " -outfmt '10 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen sstrand sseq'"
    data = subprocess.check_output(cmd, shell=True, universal_newlines=True)  # Important the E value
    data = pd.DataFrame([x.split(",") for x in data.split("\n") if x])
    data.columns = ["qseqid", "sseqid", "pident", "length", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "sstrand", "sseq"]
    return data
    

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


def repetitive_blaster(data_input, genome_fasta, folder_path, numbering, start_time, identity_1, tic_start, coincidence_data=None):
    """
    This function will iterate till a number of ``maximun_runs`` defined.

    :param genome_fasta: Path to our whole genome sequence in .fasta. It should already be a BLASTn dictionary.
    :type genome_fasta: string

    :param path_input: pandas DataFrame containing data from the first BLASTn. It's the output from :func:`~blastn_blaster`.
    :type path_input: pandas Data Frame

    :param folder_path: Path to a folder where the results will be placed. Subfolder will be created with the cromosome names.
    :type folder_path: string

    :param numbering: Indicates the number showed for the first program run. If it's 0, then, the first run will be name 0.
    :type numbering: integer

    :param maximun_runs: Indicates the last iterative execution of function.
    :type maximun_runs: integer
    """

    # Call the aesthetics function RUN identifier.
    boxymcboxface("RUN " + str(numbering))
    tic_main = time.perf_counter()  # Start the timer

    # -----------------------------------------------------------------------------
    tic = time.perf_counter()
    # First let's order the data by "sseqid", "sstrand", "sstart".
    data_ordered = data_input.sort_values(by=["sseqid", "sstrand", "sstart"])

    # Now let's group the data by "sseqid". We'll have a pandas groupby object.
    data_grouped = data_ordered.groupby("sseqid")
    toc = time.perf_counter()
    print("")
    print(f"1. Initial data:\n",
          f"\t- Data row length: {data_input.shape[0]}\n",
          f"\t- Execution time: {toc - tic:0.2f} seconds")
    # -----------------------------------------------------------------------------
    # Now let's call  `genome_specific_chromosome_main` for each chromosome_ID in the data using the groupby object.
    terminal_width = shutil.get_terminal_size().columns  # Get the terminal width
    
    print("")
    print(f"2. Individual searching and cleaning:")
    whole_group = pd.DataFrame()  # This will be the final data frame for each chromosome
    stop_dic = {}  # This will be the stop data for each chromosome
    stop_bedops_dic = {}  # This will be the stop data for each chromosome using BEDOPS
    for _, (chromosome, group) in enumerate(data_grouped):
        tic = time.perf_counter()
        now_time = datetime.now()
        formatted_now_time = now_time.strftime("%Y %B %d at %H:%M")
        print("")
        print(f"{" "*7}{"-"*74}")
        print(f"\t- {chromosome}:") 
        start_time_text = f"Program started: {start_time}"
        end_time_text = f"Program time now: {formatted_now_time}"
        RUN_text = f"RUN {numbering}"
        print(f"{RUN_text:>{terminal_width}}")
        print(f"{start_time_text:>{terminal_width}}")
        print(f"{end_time_text:>{terminal_width}}")
        
        data = genome_specific_chromosome_main(data_input=group,
                                               chromosome_ID=chromosome,
                                               main_folder_path=folder_path,
                                               genome_fasta=genome_fasta,
                                               identity_1=identity_1,
                                               run_phase=numbering,
                                               coincidence_data=coincidence_data)
        toc = time.perf_counter()
        print("")
        print(f"\t\t- Data row length: {data.shape[0]}\n",
              f"\t\t- Execution time: {toc - tic:0.2f} seconds")
        whole_group = pd.concat([whole_group, data])
    print(f"{" "*7}{"-"*74}")
    print("")
    print(f"\t- Blast data row length: {whole_group.shape[0]}\n",
          f"\t- Execution time: {toc - tic:0.2f} seconds\n")
    # -----------------------------------------------------------------------------   
    tic = time.perf_counter()
    print("")
    print(f"4. Global filtering:")
    whole_group_filtered = global_filters_main(data_input=whole_group,
                                               genome_fasta=genome_fasta,
                                               writing_path=folder_path)
    toc = time.perf_counter()
    print("")
    print(f"\t- Data row length: {whole_group_filtered.shape[0]}\n",
          f"\t- Execution time: {toc - tic:0.2f} seconds")

    # -----------------------------------------------------------------------------
    ## Save the RUN
    RUNS_folder = os.path.join(folder_path, "RUNS")  # Creates the folder for the RUNS
    os.makedirs(RUNS_folder, exist_ok=True)  # Creates the folder for the RUNS
    RUN_saver_path = os.path.join(RUNS_folder, "run_" + str(numbering) + ".csv")  # Path to save the RUN
    whole_group_filtered.to_csv(RUN_saver_path, sep=",", header=True, index=False)  # Saves the RUN
    # -----------------------------------------------------------------------------
    # Compare part
    # Prepare folders and path
    comparison_folder = os.path.join(folder_path, "comparison")
    os.makedirs(comparison_folder, exist_ok=True)

    print("")
    print(f"5. Comparison vs Old Run:")
    if coincidence_data is not None:
        print("")
        print(f"\t- Last Run data:\n",
              f"\t\t- Coincidence data row length: {coincidence_data.shape[0]}\n",
              f"\t\t- New data row length: {data_input.shape[0]}")
        # This part is important to campere with the last run whole data "whoel_group_fildered", and not only the "new_data" subset.
        data_input = pd.concat([coincidence_data, data_input], ignore_index=True)
        data_input.sort_values(by=["sseqid", "sstrand", "sstart"], inplace=True)  # Sort the data frame by the start coordinate
        print(f"\t\t- Total data row length: {data_input.shape[0]}")
    else:  # when coincidence_data == None
        print(f"\t- Last Run data:\n",
                f"\t\t- First run row length: {data_input.shape[0]}")
        
    tic = time.perf_counter()
    print("")
    print(f"\t- Results in this RUN:")
    coincidence_data, new_data, old_data_exclusive = compare_main(last_df=whole_group_filtered,
                                                                  old_df=data_input,
                                                                  folder_path=comparison_folder,
                                                                  genome_fasta=genome_fasta)
    toc = time.perf_counter()
    print("")
    print(f"\t\t- Coincidence data row length: {coincidence_data.shape[0]}\n",
          f"\t\t- New data row length: {new_data.shape[0]}\n",
          f"\t\t- Old data row length: {old_data_exclusive.shape[0]}")    
    
    # Join coincidence_data with old_data_exclusive
    if not coincidence_data.empty and not old_data_exclusive.empty:
        coincidence_data = pd.concat([coincidence_data, old_data_exclusive], ignore_index=True)
        coincidence_data.sort_values(by=["sseqid", "sstrand", "sstart"], inplace=True)
        print(f"\t\t- Coincidence data + Old data: {coincidence_data.shape[0]}")
    else:
        pass
    print(f"\t\t- Execution time: {toc - tic:0.2f} seconds")
    # -----------------------------------------------------------------------------
    # Stopping part
    stopping_folder = os.path.join(folder_path, "stopping")
    os.makedirs(stopping_folder, exist_ok=True)

    if new_data.shape[0] != 0:
        print("")
        print(f"6. Stopping criteria:")
        print(f"\t- New data found. Continuing the program.")
        print("")
    else:
        print("")
        print(f"6. Stopping criteria:")
        print(f"\t- No new data found.")
        print(f"\t- Checking similarities with the last run")
        


    # -----------------------------------------------------------------------------
    toc_main = time.perf_counter()
    print("")
    print(f"RUN {numbering} finished:\n",
        f"\t- Execution time: {toc_main - tic_main:0.2f} seconds")
    # -----------------------------------------------------------------------------
    numbering += 1  # Increase the numbering
    repetitive_blaster(data_input=new_data,
                       genome_fasta=genome_fasta,
                       folder_path=folder_path,
                       numbering=numbering,
                       start_time=start_time,
                       identity_1=identity_1,
                       tic_start=tic_start,
                       coincidence_data=coincidence_data)
                        