import os
import csv

from modules.files_manager import folder_creator, csv_creator, fasta_creator, csv_mixer
from modules.blaster import blastn_dic, blastn_blaster
from modules.seq_modifier import specific_sequence_1000nt, specific_sequence_corrected
from modules.filters import filter_by_column, global_filters_main
from modules.subfamilies_finder import subfamily_sorter


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


def specific_sequence_extractor(path_input, chromosome_ID, main_folder_path):
    """
    First it reads a .csv file and takes the rows which include in the second column (the one with the chromosome IDs), the defined chromosome. 
    Then, it creates a .csv file  with all the selectes rows in a subfolder in a specified directory using :func:`~modules.files_manager.folder_creator` and :func:`~modules.files_manager.csv_creator`.

    .. admonition:: Example of use

       If ``chromosome_ID = LinJ.02``, it will take all rows in the **.csv** where the second column contains ``LinJ.02``

    :param path_input: Path where the .csv file we'll use to read and filter data is.
    :type path_input: string

    :param chromosome_ID: Chromosome identifier, e.g., *LinJ.07*.
    :type chromosome_ID: string

    :param main_folder_path: Path where the results will be placed. It will create a subfolder with the chromosome_ID name.
    :type main_folder_path: string

    :return: A .csv file with the selected "chromosome_ID".
    :rtype: csv file
    """

    chr_x_seqs = []
    with open(path_input, "r") as main_file:
        reader = csv.reader(main_file, delimiter=",")
        for row in reader:
            if chromosome_ID in row[1]:
                chr_x_seqs.append(row)

    folder_path = main_folder_path + "/" + chromosome_ID
    folder_creator(folder_path)

    writing_path_input = main_folder_path + "/" + chromosome_ID + "/" + chromosome_ID + ".csv"

    csv_creator(writing_path_input, chr_x_seqs)

    return (folder_path, writing_path_input)  # Important. We return the new paths.


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


def genome_specific_chromosome_main(path_input, chromosome_ID, main_folder_path, genome_fasta, naming_short, max_diff):
    """
    """

    new_directories = specific_sequence_extractor(path_input, chromosome_ID, main_folder_path)  # We get a .csv with the specified chromosome_ID.
    # folder_path = new_directories[0]  # Chromosome's directory, i.e., folder_path from return (folder_path, writing_path_input)
    last_output = new_directories[1]  # Chromosome's .csv file inside Chromosome's directory, i.e., writing_path_input from return (folder_path_writing_path_input).

    # -----------------------------------------------------------------------------
    nucleotides1000_directory = specific_sequence_1000nt(last_output, chromosome_ID, main_folder_path)  # Extend sequence to 1000 nt.

    fasta_creator_output = main_folder_path + "/" + chromosome_ID + "/" + chromosome_ID + "_1000nt.fasta"
    fasta_creator(nucleotides1000_directory, fasta_creator_output)

    blastn_dic(fasta_creator_output)

    blaster_output = main_folder_path + "/" + chromosome_ID + "/" + chromosome_ID + "_1000nt_Blaster.csv"
    blastn_blaster(fasta_creator_output,
                   fasta_creator_output,
                   blaster_output,
                   "85")

    filter_by_column(blaster_output,
                     "length",
                     100,
                     blaster_output)

    # -----------------------------------------------------------------------------
    corrected_sequences = specific_sequence_corrected(blaster_output, nucleotides1000_directory, main_folder_path, chromosome_ID)

    # -----------------------------------------------------------------------------
    subfamilies_file_path_writing = main_folder_path + "/" + chromosome_ID + "/" + chromosome_ID + "_Subfamily.csv"
    subfamily_sorter(blaster_output, corrected_sequences, subfamilies_file_path_writing)

    # -----------------------------------------------------------------------------
    second_fasta_creator_output = main_folder_path + "/" + chromosome_ID + "/" + chromosome_ID + "_Corrected.fasta"
    fasta_creator(corrected_sequences, second_fasta_creator_output)

    second_blaster_output = main_folder_path + "/" + chromosome_ID + "/" + chromosome_ID + "_BLAST_MAIN.csv"
    blastn_blaster(second_fasta_creator_output,
                   genome_fasta,
                   second_blaster_output,
                   "60")

    # -----------------------------------------------------------------------------
    global_filters_main(second_blaster_output,
                        second_blaster_output,
                        genome_fasta,
                        naming_short,
                        max_diff)

    # -----------------------------------------------------------------------------
    csv_mixer_output = main_folder_path + "/" + "MIXER.csv"

    if os.path.isfile(csv_mixer_output) is False:  # #Cuando no existe, se crea
        csv_mixer(path_input, second_blaster_output, csv_mixer_output)  # Para mezclar
    else:  # Si existe ya el archivo porque ha sido creado, se cambia el path_input por csv_mixer_output
        csv_mixer(csv_mixer_output, second_blaster_output, csv_mixer_output)

# genome_specific_chromosome_main(path_input, chromosome_ID, main_folder_path, genome_fasta, naming_short, max_diff)

    # Arg 0: STRING. Directorio del archivo en formato CSV de donde leeremos y filtraremos los datos
    # Arg 1: STRING. Identificacion del cromosoma, e.g., "LinJ.07"
    # Arg 2: STRING. Directorio de la carpeta en donde se disponen los resultados del programa
    # Arg 4: STRING. Directorio del archivo en formato fasta al que queremos leer la cantidad de cromosomas, es el fasta FASTA del genoma entero
    # Arg 5: STRING. Etiqueta para leer de identificacion y numeracion de cada cromosoma en el archivo CSV. Depende del propio archivo CSV. En el caso de L. infantum es "LinJ"
    # Arg 6: INT. Numeracion con la que le indicamos el maximo valor de proximidad para las diferentes secuancias cuando tienen que ser agrupadas. MUY IMPORTANTE