

import pathlib
from datetime import datetime


def check_asteroid_ephemeris(
        bsp_filepath: pathlib.Path,
        bsp_period_start: str,
        bsp_period_end: str,
        required_period_start: str,
        required_period_end: str,
) -> bool:
    """
        Checks if the asteroid ephemeris file exists and covers the required period.

        Args:
            bsp_filepath (pathlib.Path): Path to the asteroid ephemeris file (.bsp).
            bsp_period_start (str): Start date of the ephemeris data in 'YYYY-MM-DD' format.
            bsp_period_end (str): End date of the ephemeris data in 'YYYY-MM-DD' format.
            required_period_start (str): Start date of the required period in 'YYYY-MM-DD' format.
            required_period_end (str): End date of the required period in 'YYYY-MM-DD' format.

        Returns:
            bool: True if the file exists and covers the required period.

        Raises:
            FileNotFoundError: If the ephemeris file does not exist.
            ValueError: If the ephemeris data does not cover the required period.
    """
    if not bsp_filepath.exists():
        raise FileNotFoundError(f"Asteroid ephemeris file '{bsp_filepath.name}' not found in '{bsp_filepath.parent}'.")

    bsp_start = datetime.strptime(bsp_period_start, "%Y-%m-%d").date()
    bsp_end = datetime.strptime(bsp_period_end, "%Y-%m-%d").date()

    required_start = datetime.strptime(required_period_start, "%Y-%m-%d").date()
    required_end = datetime.strptime(required_period_end, "%Y-%m-%d").date()

    if not (
        bsp_start < required_start
        and bsp_end > required_end
    ):
        raise ValueError(
            f"Asteroid ephemeris data is NOT available for the period {required_period_start} to {required_period_end}. "
            f"Available period is from {bsp_start} to {bsp_end}."
        )

    return True


def check_uncertainty_file(uncertainty_filepath: pathlib.Path) -> bool:
    """
    Checks whether the specified uncertainty file exists.

    Args:
        uncertainty_filepath (pathlib.Path): The path to the uncertainty file.

    Returns:
        bool: True if the file exists and passes all checks.

    Raises:
        FileNotFoundError: If the uncertainty file does not exist at the specified path.

    Note:
        Additional content checks can be added as needed.
    """
    if not uncertainty_filepath.exists():
        raise FileNotFoundError(f"Uncertainty file '{uncertainty_filepath.name}' not found in '{uncertainty_filepath.parent}'.")

    # TODO: Add more checks on the content of the uncertainty file if needed

    return True
