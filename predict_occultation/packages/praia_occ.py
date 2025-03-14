import logging
import pathlib
import numpy as np
from typing import Optional
import subprocess
import os
from datetime import datetime
import shutil
from utils import count_lines

def get_best_projected_search_radius(
        object_ephemeris:pathlib.Path, 
        object_diameter: Optional[float]):
    """
    Calculate the best projected search circle based on object ephemeris.

    Parameters:
        object_ephemeris (str): The path to the object ephemeris file.

    Returns:
        float: The size of the projected search circle in arcseconds.

    Notes:
        The object ephemeris file should contain distance data in a specific format.
        The function calculates the best distance from the ephemeris and computes the projected search circle size.
    """
    earth_radius = 6371  # km
    body_radius_compensation = (
        0 if object_diameter is None else object_diameter / 2
    )  # km
    distances = []
    with open(object_ephemeris, "r") as file:
        for i, line in enumerate(file, start=3):
            distances.append(line[120:160])
    distances = np.array(distances[3:], dtype=float)
    best_distance = distances.min()
    projected_search_diameter = (
        2 * (earth_radius + body_radius_compensation * 2) / best_distance
    )
    projected_search_diameter *= 3600 * 180 / np.pi  # converts to arcsec
    projected_search_radius = projected_search_diameter / 2
    return np.around(projected_search_radius, 4)

def get_minimum_outreach_radius(projected_search_radius):
    # The minimum outreach radius is X % greater than the projected search radius.
    # All stars outside this region will be discounted to accelarate computations
    return projected_search_radius * 1

def get_position_from_occ_table(data_array, index_list):
    """Function to extract right ascension or declination from a array
    The index of specific columns (RA or Dec) is defined inside of index_list
    Args:
        data_array ([type]): [description]
        index_list ([type]): [description]

    Returns:
        [type]: [description]
    """
    return [" ".join(pos) for pos in data_array[:, index_list]]

def praia_occ_input_file(
    star_catalog_filepath: pathlib.Path,
    ephemeris_filepath: pathlib.Path,
    projected_search_circle: float,
    minimum_outreach_radius: float,
    cwd: str,
):
    
    # Atenção! Nome do Arquivo de input está HARDCODED!
    filename = "praia_occ_star_search_12.dat"

    # TODO Melhorar o nome destes arquivos!
    stars_catalog_mini_filename = "g4_micro_catalog_JOHNSTON_2018"
    stars_catalog_xy_filename = "g4_occ_catalog_JOHNSTON_2018"
    stars_parameters_of_occultation_filename = "g4_occ_data_JOHNSTON_2018"
    stars_parameters_of_occultation_plot_filename = "g4_occ_data_JOHNSTON_2018_table" # Este é o arquivo de resultados com as occultaçoes

    input_template =("""{stellar_catalog}| input  file : final xy catalog with all stars from global reduction
{object_ephemeris}| input  file : ephemeris of occulting body
1                                                 | ephemeris format: 1 - NAIF ; 2 - Horizons general format; 3 -  Horizons Pluto and satellites format
DE435/JPL                                         | Ephemeris label
{stars_catalog_mini}| output file : mini-catalog with only (RA,DEC)s of all stars from the input catalog
{stars_catalog_xy}| output file : candidate stars xy catalog
{stars_parameters_of_occultation}| output file : star parameters of occultation (minimum geocentric distance,t_occ,t_initial,t_final, etc)
{stars_parameters_of_occultation_plot}| output file : star parameters of occultation (Bruno Sicardy data plot format)
{projected_search_circle}| radius (arcsec) of projected search circle (~ projected Earth figure + occultation body)
{minimum_outreach_radius}| minimum outreach radius (arcsec) for fast elimination of farway ephemeris points
12d0  12d0                                        | exclusion range of Local Solar Time (day light): min max (hours)
2005                                              | to (years - see below)
+00.00000d0                                       | bofra
+0.00000d0                                        | aofra	 Linear (RA,DEC) ephemeris drift: offra = aofra * t + bofra
+00.0000d0                                        | bofde					  offde = aofde * t + bofde
+0.00000d0                                        | aofde					  offra, offde in (mas), t in years
***********************************************************************************************************************************************************""").format(
        stellar_catalog=star_catalog_filepath.name.ljust(50),
        object_ephemeris=ephemeris_filepath.name.ljust(50),
        stars_catalog_mini=stars_catalog_mini_filename.ljust(50),
        stars_catalog_xy=stars_catalog_xy_filename.ljust(50),
        stars_parameters_of_occultation=stars_parameters_of_occultation_filename.ljust(50),
        stars_parameters_of_occultation_plot=stars_parameters_of_occultation_plot_filename.ljust(50),
        projected_search_circle=f"{projected_search_circle:2.7f}".ljust(50),
        minimum_outreach_radius=f"{minimum_outreach_radius:2.7f}".ljust(50),
    )

    filepath = pathlib.Path(cwd, filename)
    with open(filepath, "w") as f:
        f.write(input_template)

    if not filepath.exists():
        raise FileNotFoundError(f"PRAIA Input template file not found.")

    return filepath

def run_praia_occ(
    input_filepath: pathlib.Path,
    cwd: str,
    logger: logging.Logger,
):

    logger.info("Running Praia Occultation")

    log = pathlib.Path(cwd, "praia_star_search.log")

    # Diretório onde o script está sendo executado. 
    original_cwd = os.getcwd()
    logger.info(f"Original Execution CWD: [{original_cwd}]")

    try:
        logger.info(f"Changing to CWD: [{cwd}]")
        os.chdir(cwd)
        logger.info(f"Current Execution CWD: [{os.getcwd()}]")

        # Utiliza apenas o filename do input para executar o PRAIA
        # TODO: quando o praia for alterado para receber filepaths completos
        # basta utilizar o filepath.
        input = input_filepath.name

        with open(log, "w") as fp:
            p = subprocess.Popen(
                "PRAIA_occ_star_search_12 < " + input,
                stdin=subprocess.PIPE,
                shell=True,
                stdout=fp,
            )
            p.communicate()

    except Exception as e:
        msg=f"Error in PRAIA OCC: {e}"
        logger.error(msg)
        raise Exception(msg)
    finally:
        # OBRIGATÓRIO voltar para o diretório original da execução.
        os.chdir(original_cwd)
        logger.info(f"Returning to original CWD: [{original_cwd}]")

def fix_praia_occ_table(filepath: pathlib.Path):

    inoutFile = open(filepath, "r+b")
    contents = inoutFile.readlines()

    contents[4] = b" G: G magnitude from Gaia\n"
    contents[5] = contents[5][:41] + b"cluded)\n"
    contents[6] = b" G" + contents[6][2:]
    contents[17] = contents[17][:27] + b"\n"
    contents[26] = contents[26][:6] + b"only Gaia DR1 stars are used\n"
    contents[27] = contents[27][:-1] + b" (not applicable here)\n"
    contents[35] = contents[35][:34] + b"10)\n"
    contents[36] = contents[36][:41] + b"\n"
    contents[37] = contents[37][:36] + b"/yr); (0 when not provided by Gaia DR1)\n"
    contents[39] = contents[39][:115] + b"G" + contents[39][116:]

    for i in range(41, len(contents)):
        contents[i] = contents[i][:169] + b"-- -" + contents[i][173:]

    inoutFile.seek(0)  # go at the begining of the read/write file
    inoutFile.truncate()  # clean the file (delete all content)
    inoutFile.writelines(contents)  # write the new content in the blank file
    inoutFile.close()

def ascii_to_csv(ascii_filepath: pathlib.Path, csv_filepath: pathlib.Path):
    """Function to convert data from ascii table (generate by PRAIA OCC) to csv file

    Args:
        inputFile ([type]): [description]
        outputFile ([type]): [description]
    """
    data = np.loadtxt(ascii_filepath, skiprows=41, dtype=str, ndmin=2)

    nRows, nCols = data.shape

    # To avoid 60 in seconds (provided by PRAIA occ),
    date = []
    for d in data[:, range(6)]:
        if d[5] == "60.":
            d[4] = int(d[4]) + 1
            d[5] = "00."
        date.append(datetime.strptime(" ".join(d), "%d %m %Y %H %M %S."))

    # use this definition when seconds = 0..59
    # date = [datetime.strptime(' '.join(d), "%d %m %Y %H %M %S.") for d in data[:,range(6)]]

    dateAndPositions = []
    dateAndPositions.append(date)

    # Extracting positions of stars and objects and save it in a array
    for i in range(6, 17, 3):
        dateAndPositions.append(get_position_from_occ_table(data, [i, i + 1, i + 2]))

    dateAndPositions = np.array(dateAndPositions)
    dateAndPositions = dateAndPositions.T

    # Extracting others parameters (C/A, P/A, etc.)
    otherParameters = data[:, range(18, nCols)]

    newData = np.concatenate((dateAndPositions, otherParameters), 1)

    # Defining the column's names
    colNames = (
        "occultation_date;ra_star_candidate;dec_star_candidate;ra_object;"
        "dec_object;ca;pa;vel;delta;g;j;h;k;long;loc_t;"
        "off_ra;off_de;pm;ct;f;e_ra;e_de;pmra;pmde"
    )

    np.savetxt(csv_filepath, newData, fmt="%s", header=colNames, delimiter=";")

def search_candidates(
    star_catalog_filepath: pathlib.Path,
    ephemeris_filepath: pathlib.Path,
    object_diameter: float,
    cwd: str,
    logger: logging.Logger,
):

    logger.info("Calculating best projected search circle")
    projected_search_circle = get_best_projected_search_radius(
        ephemeris_filepath, object_diameter
    )
    logger.debug(f"Projected search circle: [{projected_search_circle}]")

    logger.info("Calculating minimum outreach radius")
    minimum_outreach_radius = get_minimum_outreach_radius(
        projected_search_circle
    )
    logger.debug(f"Minimum outreach radius: [{minimum_outreach_radius}]")

    logger.info("Preparing Praia Occultation Input file")
    input_filepath =  praia_occ_input_file(
        star_catalog_filepath=star_catalog_filepath, 
        ephemeris_filepath=ephemeris_filepath, 
        projected_search_circle=projected_search_circle, 
        minimum_outreach_radius=minimum_outreach_radius, 
        cwd=cwd)
    logger.info("Praia Occultation Input file generated successfully.")
    logger.debug(f"PRAIA OCC .dat file: [{input_filepath}]")


    # Running Praia Occultation
    run_praia_occ(
        input_filepath=input_filepath, 
        cwd=cwd, 
        logger=logger)

    # OBS: nome do arquivo de output está hardcoded no input file
    praia_occ_table = pathlib.Path(cwd, "g4_occ_data_JOHNSTON_2018_table")

    if not praia_occ_table.exists():
        raise FileNotFoundError(f"PRAIA Occultation table not found. {praia_occ_table}")

    # Após a criação do arquivo de occultações pelo praia é necessário corrigir o formato
    logger.info("Fixing PRAIA Occultation table")
    fix_praia_occ_table(filepath=praia_occ_table)

    logger.debug(f"PRAIA Occultation table: [{praia_occ_table}]")

    # Converter o arquivo de occultações de ascii para csv
    logger.info("Converting PRAIA Occultation table to CSV")

    occultation_table = pathlib.Path(cwd, "occultation_table.csv")
    ascii_to_csv(
        ascii_filepath=praia_occ_table,
        csv_filepath=occultation_table,
    )
    if not occultation_table.exists():
        raise FileNotFoundError(f"Occultation table not found. {occultation_table}")

    logger.debug(f"Occultation table: [{occultation_table}]")

    # COPY file to keep orinal praia csv for debug
    praia_occultation_table = pathlib.Path(cwd, "praia_occultation_table.csv")
    shutil.copy(occultation_table, praia_occultation_table)

    # Count number of occultation events for debug
    count = count_lines(occultation_table) -1   # -1 to exclude header
    logger.debug(f"Number of occultation events: [{count}]")

    return occultation_table