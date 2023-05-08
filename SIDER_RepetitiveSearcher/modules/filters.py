import csv

from modules.files_manager import csv_creator


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


def filter_by_column(path_input, column, size_filter, writing_path_input):
    """
    This function will filter a CSV data depending on ``length`` (if we want to firlter by sequence length) or ``percent`` (if we want to filter by identity percent).

    :param path_input: Path to the .csv file we want to filter data.
    :type path_input: string

    :param column: Can be ``length`` (if we want to firlter by sequence length) or ``percent`` (if we want to filter by identity percent)
    :type column: string

    :param size_filter: Number to filter dependint of the **column** argument.
    :type size_filter: integer

    :param writing_path_input: Path to the CSV file this function will create and save
    :type writing_path_input: string

    :return: A CSV file with the dalta filtered depending on the **column** and **size_filter** argumetns.
    :rtype: CSV file
    """

    if column == "length":
        column = 3
    elif column == "percent":
        column = 2

    matrix_filter_by_column = []
    with open(path_input, "r") as main_file:
        reader = csv.reader(main_file, delimiter=",")  # Recordar que antes al poner outfmt 10, ahora estan separados por comas.
        for row in reader:
            if column == 3:
                if 1000 >= int(row[column]) >= size_filter:  # ##1000 por el problema de duplicaciones que me dijo Requena
                    matrix_filter_by_column.append(row)
            elif column == 2:
                if float(row[column]) >= size_filter:  # Necesario para pasar de STRING a FLOAT
                    matrix_filter_by_column.append(row)

    csv_creator(writing_path_input, matrix_filter_by_column)

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------abs

# 1)Este junta todos los anteriores de la primera prueba, cada uno escribe un CSV y lo va sobreescribiendo constantemente hasta el final.


def Global_Filters_Main(path_input, writing_path_input, genome_fasta, naming_short, max_diff):

    column = "length"
    size_filter = 100
    filter_by_column(path_input, column, size_filter, writing_path_input)

    path_input = writing_path_input  # Así le decimos que el archivo de entrada es el de salida del anterior, y que en el mismo, escriba los nuevos datos
    dash_filter(path_input, writing_path_input)
    
    genome_duplicate_filter(genome_fasta, naming_short, path_input, max_diff, writing_path_input)

    Genome_Solap_Main(genome_fasta, naming_short, path_input, max_diff, writing_path_input)
    # En este ultima ya imprime los resultados finales

# Global_Filters_Main(path_input, writing_path_input, genome_fasta, naming_short, max_diff):

    # Arg 0: STRING. Directorio del archivo en formato CSV de donde leeremos y filtraremos los datos
    # Arg 1: STRING. Directorio del archivo en formato CSV en donde guardaremos los resultados del filtrado, Recordar la extension .csv
    # Arg 2: STRING. Directorio del archivo en formato fasta al que queremos leer la cantidad de cromosomas, es el fasta FASTA del genoma entero
    # Arg 3: STRING. Etiqueta para leer de identificacion y numeracion de cada cromosoma en el archivo CSV. Depende del propio archivo CSV. En el caso de L. infantum es "LinJ"
    # Arg 4: INT. Numeracion con la que le indicamos el maximo valor de proximidad para las diferentes secuancias cuando tienen que ser agrupadas. MUY IMPORTANTE