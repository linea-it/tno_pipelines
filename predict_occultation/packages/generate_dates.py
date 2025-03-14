import os
import subprocess
import logging
from pathlib import Path

def generate_dates_file(start_date:str, end_date:str, ephemeris_step:int, cwd:str, logger:logging.Logger) -> Path:
    """
    Generates a dates file within the specified directory.

    This function generates a file named 'dates.txt' containing a list of dates
    from `start_date` to `end_date` with a step interval of `step` in seconds.

    Args:
        start_date (str): The start date in the format 'YYYY-MM-DD'.
        end_date (str): The final date in the format 'YYYY-MM-DD'.
        ephemeris_step (int): The step interval in seconds between each date.
        cwd (str): The current working directory where the file will be generated.
        logger (logging.Logger): A logger instance to log information and errors.

    Returns:
        Path: The path to the generated or existing 'dates.txt' file.

    Raises:
        Exception: If there is an error during the file generation process or if the file is not generated.
    """

    logger.info("Generating dates file")

    filepath = Path(cwd, "dates.txt")

    # Diretório onde o script está sendo executado. 
    original_cwd = os.getcwd()
    logger.info(f"Original Execution CWD: [{original_cwd}]")

    try:
        logger.info(f"Changing to CWD: [{cwd}]")
        os.chdir(cwd)
        logger.info(f"Current Execution CWD: [{os.getcwd()}]")

        # IMPORTANT: Precisa do arquivo leap_seconds (naif0012.tls) no diretório de execução
        logger.info(f"Generating dates file: start:[{start_date}] end: [{end_date}] step: [{ephemeris_step}]")
        with open(filepath, "w") as fp:
            parameters = [start_date, end_date, ephemeris_step]
            strParameters = "\n".join(map(str, parameters))
            p = subprocess.Popen(
                "geradata", stdin=subprocess.PIPE, stdout=fp, shell=True
            )
            p.communicate(str.encode(strParameters))

    except Exception as e:
        msg="Error generating dates file: %s", e
        logger.error(msg)
        raise Exception(msg)
    finally:
        # OBRIGATÓRIO voltar para o diretório original da execução.
        os.chdir(original_cwd)
        logger.info(f"Returning to original CWD: [{original_cwd}]")

    if filepath.exists():
        logger.info(f"Dates file created: [{filepath}]")
        return filepath
    else:
        msg="Date file not generated. [%s]" % filepath
        logger.error(msg)
        raise Exception(msg)