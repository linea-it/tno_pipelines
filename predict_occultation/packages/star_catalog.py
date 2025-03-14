import logging
import pathlib
from typing import Dict, Optional, List
from dao.gaia import GaiaDao
def generate_star_catalog(
    ra: List,
    dec: List,
    angular_diameter: List,
    maximum_visual_magnitude: int, 
    catalog_name: str,
    catalog_display_name: str,
    catalog_schema: str,
    catalog_tablename: str,
    catalog_ra_property: str,
    catalog_dec_property: str,
    cwd: str, 
    logger: logging.Logger
) -> pathlib.Path:

    logger.info("Connecting to Catalog Database.")

    dao = GaiaDao(
        name = catalog_name,
        display_name=catalog_display_name,
        schema=catalog_schema,
        tablename=catalog_tablename,
        ra_property=catalog_ra_property,
        dec_property=catalog_dec_property,
        logger=logger
    )

    # Quering catalog
    df_catalog = dao.catalog_by_polygons(
        ra=ra,
        dec=dec,
        angular_diameter=angular_diameter,
        max_mag=maximum_visual_magnitude,
    )

    # Cria um arquivo no formato especifico do praia_occ
    logger.info("Writing star catalog to files (.cat, .csv).")
    praia_star_catalog_filepath = dao.write_gaia_catalog(
        df_catalog.to_dict("records"), pathlib.Path(cwd, "star_catalog.cat")) 
    logger.debug(f"Star Catalog .cat file: {praia_star_catalog_filepath}")

    # Cria um arquivo csv com os dados do catalogo
    star_catalog_filepath = dao.gaia_catalog_to_csv(df_catalog, pathlib.Path(cwd, "star_catalog.csv"))
    logger.debug(f"Star Catalog .csv file: {star_catalog_filepath}")
    
    logger.info("Star Catalog generated successfully.")

    return praia_star_catalog_filepath, star_catalog_filepath
