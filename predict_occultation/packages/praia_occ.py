import logging
from generate_dates import generate_dates_file
from pathlib import Path

def run_praia_occ(
    name: str,
    start_date: str,
    end_date: str,
    step: int,
    maximum_visual_magnitude: float,
    leap_seconds: str,
    planetary_ephemeris: str,
    cwd: str,
    logger: logging.Logger,
):

    logger.info("Preparing Praia Occultation")
    logger.debug("Name: %s", name)
    logger.debug("Start Date: %s", start_date)
    logger.debug("End Date: %s", end_date)
    logger.debug("Step: %s", step)
    logger.debug("Leap Seconds: %s", leap_seconds)
    logger.debug("Planetary Ephemeris: %s", planetary_ephemeris)
    logger.debug("Maximum Visual Magnitude: %s", maximum_visual_magnitude)


    # # If the dates.txt file already exists in the specified directory, it will use the existing file.
    # dates_filepath = Path(cwd, "dates.txt")
    # if dates_filepath.exists():
    #     logger.info("Using previously created date file.")
    # else:
    #     dates_filepath = generate_dates_file(start_date, end_date, step, cwd, logger)



    print(dates_filepath)

    # Gera o arquivo de efemérides para o período de interesse