import pdb  # In case of debbuging
import csv
import subprocess

# from modules.filters import chromosome_filter  # Don't call --> ciruclar import
from modules.files_manager import csv_creator


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


def genome_solap_location_filter(chromosome_rows):  # Todo STRING menos chromosome_rows
    """
    For the DNA ("minus" or "plus" strand) we get the *start* and *end* coordinates.

    :param chromosome_rows: Given by :func:`~modules.overlap.genome_solap_location_grouping`. It's a Python list with all the rows from a CSV of one specific chromosome.
    :type chromosome_rows: python list

    :return: 4 lists with the start and end coordinates for "minus" and "plus" DNA strands.
    :rtype: python list
    """
    pos_plus_start = []  # start position for "plus" strand
    pos_plus_end = []  # end position for "plus" strand
    pos_minus_start = []  # start position for "minus" strand
    pos_minus_end = []  # end position for "minus" strand

    for row in chromosome_rows:
        if "plus" in row[14]:  # row[14] is "minus" or "plus".
            pos_plus_start.append(int(row[10]))
            pos_plus_end.append(int(row[11]))
        else:  # For "minus" strand
            pos_minus_start.append(int(row[10]))
            pos_minus_end.append(int(row[11]))

    return (pos_plus_start, pos_plus_end, pos_minus_start, pos_minus_end)


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


def genome_solap_location_grouping(chromosome_rows, DNA_sense, max_diff):
    """
    This function will analyze each position of coordinates (start and end) for each DNA strand. It will group them depending on nearness.

    :param chromosome_rows: Given by :func:`~modules.overlap.genome_solap_minmax`. It's a Python list with all the rows from a CSV of one specific chromosome.
    :type chromosome_rows: python list

    :param DNA_sense: Maximun proxomity value for the different sequences when they have to be grouped. **Important**.
    :type DNA_sense: integer

    :param max_diff:
    :type max_diff:

    :return: An Array 3D with the positions for "starts" and ends" in "groups" like ``[[[group1 starts][group2 starts][...]]   [[group1 ends][group2 ends][...]]]``
    :rtype: Array 3D
    """
    # Numbers: 0 for Plus:Start, 1 for Plus:End, 2 for Minus:Start, 3 for Minus:End
    # -----------------------------------------------------------------------------
    # Here it will adjust the matrix depending on the DNA strand.
    # -----------------------------------------------------------------------------
    position_list = genome_solap_location_filter(chromosome_rows)
    if DNA_sense == "plus":
        position_list_main = [position_list[0], position_list[1]]
    elif DNA_sense == "minus":
        position_list_main = [position_list[2], position_list[3]]

    matrix1 = []
    matrix2 = []

    for main_list in position_list_main:  # position_list_main = [[Starts][Ends]]. First it will chose [Starts] and then [Ends].
        if main_list == position_list_main[0]:  # If the list == Start, it goes to matrix1
            matrix = matrix1
        elif main_list == position_list_main[1]:  # If the list == End, it goes to matrix2
            matrix = matrix2

        # -----------------------------------------------------------------------------
        # This code is hard to understand. It will group the coordinates depending on nearness.
        # For example, the list [10, 12, 30, 34, 60, 62] will be grouped in [[12, 12] [30, 34] [60, 62]]
        # -----------------------------------------------------------------------------
        for position in main_list:  # "position" is an INT. A specific coordinate.
            main_statement = False  # If it stays "False", it will create a new group with `matrix.append([position])`
            for group in matrix:  # In the first iteration, matrix is empty
                for member in group:  # "member" is inside "group" which is inside "matrix" python list.
                    if abs(member - position) <= max_diff:
                        group.append(position)  # If it's near, we'll append it to the "group"
                        main_statement = True
                        break  # It breaks the code, exit this loot "for member..." and continue with "if not main_statement"
            if not main_statement:  # If "group" its empty at first, it will go here first. The number will be added like ([]) inside matrix, creating a "group"
                matrix.append([position])  # If not found, We create a List inside a list [[x]], which I called before as "group".

    # After "for main_list..." we'll have coordinates for "Starts" and "Ends". We simply join them.
    matrix_main = [matrix1, matrix2]

    return (matrix_main)


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


def genome_solap_minmax(chromosome_rows, max_diff):  # Todo STRING menos max_diff
    """
    After doing the groupings with :func:`~modules.overlap.genome_solap_location_grouping` for each DNA strand. We simply get the minimun or maximun value depending if it's a groupings of "starts" or "ends", but mainly on the *strand* (because "plus" and "minus" strands have like inverted coordinates).

    :param chromosome_rows: Given by :func:`~modules.overlap.genome_solap_main`. It's a Python list with all the rows from a CSV of one specific chromosome.
    :type chromosome_rows: python list

    :param max_diff: Maximun proxomity value for the different sequences when they have to be grouped. **Important**.
    :type max_diff: integer

    :return: 4 list with the minimum or maximun values depending if it's a groupings of "starts" or "ends", but mainly on the *strand*
    :rtype: 4 python lists.
    """
    # We get the coordinates groupings in an array 3D for each strand.
    plus = genome_solap_location_grouping(chromosome_rows, "plus", max_diff)
    minus = genome_solap_location_grouping(chromosome_rows, "minus", max_diff)

    # For each group (see "genome_solap_location_grouping"):
    plus_min = [min(x) for x in plus[0]]  # "plus" strand "starts", we get a list with the minimum value of each group
    plus_max = [max(x) for x in plus[1]]  # "plus" strand "ends", we get a list with the maximum value of each group
    minus_max = [max(x) for x in minus[0]]  # "minus" strand "starts", we get  a list with the maximum value of each group
    minus_min = [min(x) for x in minus[1]]  # "minus" strand "ends", we get a list with the minimum value of each group

    return (plus_min, plus_max, minus_max, minus_min)


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


def genome_solap_by_pairs(rows_to_filter, genome_fasta):
    """
    Partials overlaps will be filtered and joined.

    :param rows_to_filter: A Array 3D given by :func:`~modules.overlap.genome_solap_main` with all the supossed rows with overlaps.
    :type rows_to_filter: array 3D

    :return: a 3D list with all the real coordinates for overlaps
    :rtype: a Python 3D list
    """
    # -----------------------------------------------------------------------------
    # IMPORTANT. It organize the rows by their coordinates by the DNA strand. So the row <i> coordinate is the nearest to the row <i+1> coordinate.
    # -----------------------------------------------------------------------------
    plus_list = [x for x in rows_to_filter if "plus" in x[14]]  # We select the "plus" side of the list
    plus_list = sorted(plus_list, key=lambda x: x[10], reverse=True)  # We order them. The order is with "strings" not with numbers. The result is the same.
    minus_list = [x for x in rows_to_filter if "minus" in x[14]]  # We do the same with minus strand
    minus_list = sorted(minus_list, key=lambda x: x[10], reverse=True)

    rows_to_filter = plus_list + minus_list  # We update "rows_to_filter" with the new ones.
    # -----------------------------------------------------------------------------

    rows_final = []
    for first, second in zip(*[iter(rows_to_filter)] * 2):  # IMPORTANT. It uses iter to select the rows by pairs. So it selects at the same time row <i> and row <i+1>, which coordinates are the nearest.
        two_sequence_rec = []
        two_sequence_rec.append(first)
        two_sequence_rec.append(second)
        # At the end, "two_sequence_rec" will have only 2 rows/sequences

        sequence_start = []
        sequence_end = []
        homology1 = []
        e_value1 = []
        bit_score1 = []

        for sequence in two_sequence_rec:
            sequence_start.append(sequence[10])
            sequence_end.append(sequence[11])
            homology1.append(float(sequence[2]))  # This one, and the ones after are here to maintain the data, but they are not neccesary
            e_value1.append(float(sequence[12]))
            bit_score1.append(float(sequence[13]))

        homology2 = str(round((homology1[0] + homology1[1]) / 2, 3))
        e_value2 = (e_value1[0] + e_value1[1]) / 2
        e_value2 = str("{:.2e}".format(e_value2))
        bit_score2 = str(round((bit_score1[0] + bit_score1[1]) / 2, 1))

        if "plus" in first[14] and abs(int(first[10]) - int(second[10])) <= 1000:  # This number is important. It wll only continue if they are near (just in case).
            min_start = min(sequence_start)
            max_end = max(sequence_end)
            seq_length = str(int(max_end) - int(min_start) + 1)

            seq = subprocess.check_output("blastdbcmd -db " + genome_fasta + " -entry "
                                          + first[1] + " -range " + min_start + "-" + max_end
                                          + " -strand plus -outfmt %s",
                                          shell=True,
                                          universal_newlines=True)  # subprocess is really important
            seq = seq.strip()  # Remove EoL characters

            new_row = [first[0], first[1], homology2, seq_length, first[4], first[5], "", "", "", "", str(min_start), str(max_end), e_value2, bit_score2, first[14], seq]

            rows_final.append(new_row)

        elif "minus" in first[14] and abs(int(first[10]) - int(second[10])) <= 1000:
            max_start = max(sequence_start)
            min_end = min(sequence_end)
            seq_length = str(int(max_start) - int(min_end) + 1)

            seq = subprocess.check_output("blastdbcmd -db " + genome_fasta + " -entry "
                                          + first[1] + " -range " + min_end + "-" + max_start
                                          + " -strand minus -outfmt %s",
                                          shell=True,
                                          universal_newlines=True)  # subprocess is really important
            seq = seq.strip()  # Remove EoL characters

            new_row = [first[0], first[1], homology2, seq_length, first[4], first[5], "", "", "", "", str(max_start), str(min_end), e_value2, bit_score2, first[14], seq]

            rows_final.append(new_row)

    return (rows_final)

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


def genome_solap_main(genome_fasta, naming_short, path_input, max_diff, writing_path_input):  # Todo STRING menos max_diff
    """
    Will call a variety of functions with the purpose of filtering overlaps in the data:

    1. First it will get all the coordinates from ``path_input`` with :func:`~modules.overlap.genome_solap_location_filter`.
    2. The output will be given to :func:`~modules.overlap.genome_solap_location_grouping`, that will group them in a 3D list by nearness.
    3. The output will be given to :func:`~modules.overlap.genome_solap_minmax` that will get depending on the DNA strand, the maximum or minimum value for each group.
    4. The ouput will be given to :func:`~modules.overlap.genome_solap_main` that will make two filters:
    
       - The first one, will filters all sequences which have the minimum value as well as the maximum (depending again on the DNA strand).
       - The second one will get all the sequences which have only one minimum or one maximum. In the end, since it will search with the coordinates given by :func:`~modules.overlap.genome_solap_minmax`, the number of results will be even (one with the minimum and one with the maximum).
       
    5. The ourput of the second filter of the 4º step will be given to `func:`~modules.overlap.genome_solap_by_pairs` that will first order the inputs, so the when it analyze them by pairs, the pair are in fact, the nearest ones. Then, it will join them to form one sequence which will be added to a final 3D list.
    6. That 3D list will be used to make a CSV file.

    It uses :func:`~modules.filters.chromosome_filter`

    :param genome_fasta: Path to our whole genome sequence in FASTA format.
    :type genome_fasta: string

    :param naming_short: Label needed to read the ID of each cromosome in the .csv file. In the case of **L. infantum** for example, would be *LinJ* since the .csv file IDs are *LinJ.XX*.
    :type naming_short: string

    :param path_input: Path to the CSV file we want to filter data. It's the output file created by :func:`~modules.blaster.blastn_blaster` and given here by :func:`modules.filters.global_filters_main`. It's called "_BLAST_MAIN.csv".
    :type path_input: string

    :param max_diff: Maximun proxomity value for the different sequences when they have to be grouped. **Important**.
    :type max_diff: intenger

    :param writing_path_input: Path where the CSV file will be saved.
    :type writing_path_input: string
    """
    from modules.filters import chromosome_filter  # Delayed import --> to break the ciruclar import. Need to be at the start of function.

    print("\n", "=" * 50, "\nFiltering overlaps proceeding:\n", "=" * 50, sep="")

    genome_solap_main_matrix = []
    chromosome_number = chromosome_filter(genome_fasta, naming_short)  # It gets the names for the all the chromosome's IDs in the fasta file, e.g., "LinJ.01", "LinJ.02", etc.

    # -----------------------------------------------------------------------------
    # 1) It'll extract all the rows from the CSV file that has the same chromosome number, e.g., "LinJ.01".
    # -----------------------------------------------------------------------------
    for chromosome in chromosome_number:
        solap_main_matrix = []
        chromosome_rows = []  # We get the rows from the CSV from a "chromosome"
        with open(path_input, "r") as main_file:  # We read the CSV "_BLAST_MAIN.csv"
            reader = csv.reader(main_file, delimiter=",")
            for row in reader:
                if chromosome in row[1]:  # chromosome filter
                    chromosome_rows.append(row)

    # -----------------------------------------------------------------------------
    # 2.1) For each DNA strand and start or end position --> it will group them by nearness.
    # 2.2) From those groups --> it will get ONLY the minimum and maximum values depending on DNA strand.
    # 2.3) It will save them in an 3D Array, called `minmax`.
    # -----------------------------------------------------------------------------
        minmax = genome_solap_minmax(chromosome_rows, max_diff)  # Here we get minimun and maximums in an 3D Array.
        # [0] --> plus_min | [1] --> plus_max | [2] --> pinus_max | [3] --> pinus_min
        plus_start_matrix = []
        plus_end_matrix = []
        minus_start_matrix = []
        minus_end_matrix = []

        # -----------------------------------------------------------------------------
        # 3) For the cases where there are OVERLAPS, but one seq's got the MAX and MIN values. Both of them.
        # -----------------------------------------------------------------------------
        with open(path_input, "r") as main_file:
            reader = csv.reader(main_file, delimiter=",")
            for row in reader:
                if chromosome in row[1]:  # chromosome ID filter
                    if "plus" in row[14]:
                        if int(row[10]) in minmax[0] and int(row[11]) in minmax[1]:  # It makes sure the seq's got the MIN and MAX values.
                            solap_main_matrix.append(row)  # It saves the whole row.
                            plus_start_matrix.append(int(row[10]))  # It saves the start coordinate, to not repeat in the next iterations.
                            plus_end_matrix.append(int(row[11]))  # It saves the end coordinate, to not repeat in the next iterations.
                    elif "minus" in row[14]:
                        if int(row[10]) in minmax[2] and int(row[11]) in minmax[3]:
                            solap_main_matrix.append(row)
                            minus_start_matrix.append(int(row[10]))
                            minus_end_matrix.append(int(row[11]))
        # -----------------------------------------------------------------------------
        # 4.1) For the cases where there are OVERLAPS, but one seq's got the MAX or MIN values. Only one of them.
        # 4.2) Probably another seq will have the other value. It needs to find it and join them.
        # 4.3) It will get 4 ARRAYS (2 for plus strand and 2 for minus strand):
                    # 4.3.1) Seqs plus with only Start coordinate found.
                    # 4.3.2) Seqs plus with only End coordinate found.
                    # 4.3.3) Seqs minus with only Start coordinate found.
                    # 4.3.4) Seqs minus with only End coordinate found.
        # -----------------------------------------------------------------------------
        solap_segments = []  # Here are all the small overlaps segments
        solap_segments_plus_start = []
        solap_segments_plus_end = []
        solap_segments_minus_start = []
        solap_segments_minus_end = []
        with open(path_input, "r") as main_file:
            reader = csv.reader(main_file, delimiter=",")
            for row in reader:
                if chromosome in row[1]:
                    if "plus" in row[14]:
                        if int(row[10]) not in plus_start_matrix and int(row[10]) in minmax[0]:
                            if int(row[10]) not in solap_segments_plus_start:
                                solap_segments.append(row)  # Here we save the whole row
                                solap_segments_plus_start.append(int(row[10]))  # Here we save only the coordinate

                        if int(row[11]) not in plus_end_matrix and int(row[11]) in minmax[1]:
                            if int(row[11]) not in solap_segments_plus_end:
                                solap_segments.append(row)
                                solap_segments_plus_end.append(int(row[11]))

                    elif "minus" in row[14]:
                        if int(row[10]) not in minus_start_matrix and int(row[10]) in minmax[2]:
                            if int(row[10]) not in solap_segments_minus_start:
                                solap_segments.append(row)
                                solap_segments_minus_start.append(int(row[10]))

                        if int(row[11]) not in minus_end_matrix and int(row[11]) in minmax[3]:
                            if int(row[11]) not in solap_segments_minus_end:
                                solap_segments.append(row)
                                solap_segments_minus_end.append(int(row[11]))

        # -----------------------------------------------------------------------------
        # 5) Here, `genome_solap_by_pairs` it will join the small segments of overlaps.
        # -----------------------------------------------------------------------------
        solap_by_pairs_definitive = genome_solap_by_pairs(solap_segments, genome_fasta)  # Very IMPORTANT.
        solap_main_matrix += solap_by_pairs_definitive

        genome_solap_main_matrix += solap_main_matrix

    csv_creator(writing_path_input, genome_solap_main_matrix)
