from typing import Dict
import pathlib
import logging
import pandas as pd
from datetime import datetime, timezone
from library import (
    ra_hms_to_deg, 
    dec_hms_to_deg, 
    get_mag_ra_dec_uncertainties_interpolator,
    asteroid_visual_magnitude,
    compute_magnitude_drop,
    get_apparent_diameter,
    get_event_duration,
    get_moon_and_sun_separation,
    get_moon_illuminated_fraction,
    get_instant_uncertainty,
    get_closest_approach_uncertainty,
    generate_hash
)
from occviz import occultation_path_coeff
import json
import numpy as np
from astropy.time import Time

def prepare_occultation_table_dataframe(
    occultation_table: pathlib.Path, 
    object_data: Dict,
):
    # Le o arquivo occultation table e cria um dataframe
    df = pd.read_csv(
        occultation_table,
        delimiter=";",
        header=None,
        skiprows=1,
        names=[
            "occultation_date",
            "ra_star_candidate",
            "dec_star_candidate",
            "ra_object",
            "dec_object",
            "ca",
            "pa",
            "vel",
            "delta",
            "g",
            "j",
            "h",
            "k",
            "long",
            "loc_t",
            "off_ra",
            "off_de",
            "pm",
            "ct",
            "f",
            "e_ra",
            "e_de",
            "pmra",
            "pmde",
        ],
    )

    # Adiciona as colunas de coordenadas de target e star convertidas para degrees.
    df["ra_target_deg"] = df["ra_object"].apply(ra_hms_to_deg)
    df["dec_target_deg"] = df["dec_object"].apply(dec_hms_to_deg)
    df["ra_star_deg"] = df["ra_star_candidate"].apply(ra_hms_to_deg)
    df["dec_star_deg"] = df["dec_star_candidate"].apply(dec_hms_to_deg)

    # Remover valores como -- ou -
    df["ct"] = df["ct"].str.replace("--", "")
    df["f"] = df["f"].str.replace("-", "")

    # Altera o nome das colunas
    df = df.rename(
        columns={
            "occultation_date": "date_time",
            "ra_object": "ra_target",
            "dec_object": "dec_target",
            "ca": "closest_approach",
            "pa": "position_angle",
            "vel": "velocity",
            "g": "g_star",
            "j": "j_star",
            "h": "h_star",
            "k": "k_star",
            "off_de": "off_dec",
            "pm": "proper_motion",
            "f": "multiplicity_flag",
            "e_de": "e_dec",
            "pmde": "pmdec",
        }
    )

    # Correcao de valores nao validos
    # Fix https://github.com/linea-it/tno_pipelines/issues/10.
    df.loc[df["j_star"] == 50, "j_star"] = None
    df.loc[df["h_star"] == 50, "h_star"] = None
    df.loc[df["k_star"] == 50, "k_star"] = None

    # -------------------------------------------------
    # MPC asteroid data used for prediction
    # -------------------------------------------------
    ast_data_columns = [
        "name",
        "number",
        "base_dynclass",
        "dynclass",
        "albedo",
        "albedo_err_max",
        "albedo_err_min",
        "alias",
        "aphelion",
        "arg_perihelion",
        "astorb_dynbaseclass",
        "astorb_dynsubclass",
        "density",
        "density_err_max",
        "density_err_min",
        "diameter",
        "diameter_err_max",
        "diameter_err_min",
        "epoch",
        "eccentricity",
        "g",
        "h",
        "inclination",
        "last_obs_included",
        "long_asc_node",
        "mass",
        "mass_err_max",
        "mass_err_min",
        "mean_anomaly",
        "mean_daily_motion",
        "mpc_critical_list",
        "perihelion",
        "pha_flag",
        "principal_designation",
        "rms",
        "semimajor_axis",
    ]
    for column in ast_data_columns:
        df[column] = object_data.get(column)

    # Correcao de valores nao validos para H (magnitude absoluta do asteroide)
    if object_data["h"] is not None and object_data["h"] > 99:
        df["h"] = None

    # -------------------------------------------------
    # Provenance Fields
    # Adiciona algumas informacoes de Proveniencia a cada evento de predicao
    # -------------------------------------------------
    df["job_id"] = object_data["job_id"]
    df["catalog"] = object_data["predict_params"]["star_catalog"]["display_name"]
    df["predict_step"] = object_data["predict_params"]["ephemeris_step"]
    df["bsp_source"] = object_data["bsp_jpl"]["source"]
    df["bsp_planetary"] = object_data["predict_params"]["planetary_ephemeris"]["display_name"]
    df["leap_seconds"] = object_data["predict_params"]["leap_seconds"]["display_name"] 
    df["nima"] = False # NIMA nao esta sendo utilizado
    df["created_at"] = datetime.now(tz=timezone.utc)

    # -------------------------------------------------
    # Coeff paths e outras grandezas
    # -------------------------------------------------
    df["have_path_coeff"] = False
    df["occ_path_coeff"] = {}
    columns = [
        "gaia_source_id", 
        "g_star",
        "apparent_magnitude",
        "magnitude_drop",
        "apparent_diameter",
        "event_duration",
        "instant_uncertainty",
        "closest_approach_uncertainty",
        "closest_approach_uncertainty_km",
        "moon_separation",
        "sun_elongation",
        "moon_illuminated_fraction",
        "e_ra_target",
        "e_dec_target",
        "e_ra",
        "e_dec",
        "occ_path_max_longitude",
        "occ_path_min_longitude",
        "occ_path_is_nightside",
        "occ_path_max_latitude",
        "occ_path_min_latitude",
    ]
    for column in columns:
        df[column] = None

    # ------------------------------------------------------
    # Colunas que aparentemente não esto sendo preenchidas
    # ------------------------------------------------------
    columns_for_future = [
        "g_mag_vel_corrected",
        "rp_mag_vel_corrected",
        "h_mag_vel_corrected",
        "ra_star_with_pm",
        "dec_star_with_pm",
        "ra_star_to_date",
        "dec_star_to_date",
        "ra_target_apparent",
        "dec_target_apparent",
        "ephemeris_version",
        "probability_of_centrality",
    ]
    for column in columns_for_future:
        df[column] = None

    # Hash id
    df["hash_id"] = None

    # Altera a ordem das colunas para coincidir com a da tabela
    df = df.reindex(
        columns=[
            "name",
            "number",
            "date_time",
            "gaia_source_id",
            "ra_star_candidate",
            "dec_star_candidate",
            "ra_target",
            "dec_target",
            "closest_approach",
            "position_angle",
            "velocity",
            "delta",
            "g",
            "j_star",
            "h",
            "k_star",
            "long",
            "loc_t",
            "off_ra",
            "off_dec",
            "proper_motion",
            "ct",
            "multiplicity_flag",
            "e_ra",
            "e_dec",
            "pmra",
            "pmdec",
            "ra_star_deg",
            "dec_star_deg",
            "ra_target_deg",
            "dec_target_deg",
            "created_at",
            "apparent_diameter",
            "aphelion",
            "apparent_magnitude",
            "dec_star_to_date",
            "dec_star_with_pm",
            "dec_target_apparent",
            "diameter",
            "e_dec_target",
            "e_ra_target",
            "ephemeris_version",
            "g_mag_vel_corrected",
            "h_mag_vel_corrected",
            "inclination",
            "instant_uncertainty",
            "magnitude_drop",
            "perihelion",
            "ra_star_to_date",
            "ra_star_with_pm",
            "ra_target_apparent",
            "rp_mag_vel_corrected",
            "semimajor_axis",
            "have_path_coeff",
            "occ_path_max_longitude",
            "occ_path_min_longitude",
            "occ_path_coeff",
            "occ_path_is_nightside",
            "occ_path_max_latitude",
            "occ_path_min_latitude",
            "base_dynclass",
            "bsp_planetary",
            "bsp_source",
            "catalog",
            "dynclass",
            "job_id",
            "leap_seconds",
            "nima",
            "predict_step",
            "albedo",
            "albedo_err_max",
            "albedo_err_min",
            "alias",
            "aphelion",
            "arg_perihelion",
            "astorb_dynbaseclass",
            "astorb_dynsubclass",
            "density",
            "density_err_max",
            "density_err_min",
            "diameter_err_max",
            "diameter_err_min",
            "epoch",
            "eccentricity",
            "last_obs_included",
            "long_asc_node",
            "mass",
            "mass_err_max",
            "mass_err_min",
            "mean_anomaly",
            "mean_daily_motion",
            "mpc_critical_list",
            "perihelion_dist",
            "pha_flag",
            "principal_designation",
            "rms",
            "g_star",
            "h_star",
            "event_duration",
            "moon_separation",
            "sun_elongation",
            "closest_approach_uncertainty",
            "moon_illuminated_fraction",
            "probability_of_centrality",
            "hash_id",
            "closest_approach_uncertainty_km",
        ]
    )
    df.to_csv(occultation_table, index=False, sep=";")


def execute_path_coeff(
    object_data: Dict,
    occultation_table: pathlib.Path,
    star_catalog: pathlib.Path,
    object_ephemeris: pathlib.Path,
    mag_and_uncert: pathlib.Path,
    planetary_ephemeris: pathlib.Path,
    leap_seconds: pathlib.Path,    
    cwd: str,
    logger: logging.Logger,
):

    # Le o arquivo de ocultações
    df = pd.read_csv(
        occultation_table,
        delimiter=";")

    # Le o arquivo de incertezas
    logger.info("Reading the uncertainties file")
    has_uncertainties = False
    mag_cs, ra_cs, dec_cs = None, None, None
    if mag_and_uncert.exists():
        # Arquivo com incertezas nao foi criado
        with open(mag_and_uncert, "r") as f:
            mag_and_uncertainties = json.load(f)
        # Interpolador para as incertezas
        mag_cs, ra_cs, dec_cs = get_mag_ra_dec_uncertainties_interpolator(
            **mag_and_uncertainties
        )

    if ra_cs and dec_cs:
        has_uncertainties = True
        logger.debug(f"Has uncertainties: [{has_uncertainties}]")
    else:
        logger.warning(f"Uncertainties file does not exist. {str(mag_and_uncert)}")
   
    logger.info("Reading the star catalog")
    # Le o source_id, ra, dec e g magnitude do catalogo gaia csv
    df_gaia_csv = pd.read_csv(
        star_catalog,
        usecols=(0, 1, 2, 3, 4, 13),
        delimiter=";",
    )

    logger.info("Iterate over rows to calculate the occultation path")
    coeff_paths = []
    # Para cada Ocultacao e necessario calcular o occultation path.
    for row in df.to_dict(orient="records"):

        new_row = {
            "hash_id": None,
            "gaia_source_id": None,
            "gaia_g_mag": None,
            "apparent_magnitude": None,
            "magnitude_drop": None,
            "apparent_diameter": None,
            "event_duration": None,
            "instant_uncertainty": None,
            "closest_approach_uncertainty": None,
            "closest_approach_uncertainty_km": None,
            "moon_separation": None,
            "sun_elongation": None,
            "moon_illuminated_fraction": None,
            "e_ra_target": None,
            "e_dec_target": None,
            "e_ra_star": None,
            "e_dec_star": None,
            "have_path_coeff": False,
            "occ_path_min_longitude": None,
            "occ_path_max_longitude": None,
            "occ_path_min_latitude": None,
            "occ_path_max_latitude": None,
            "occ_path_is_nightside": None,
            "occ_path_coeff": {},
        }

        # Obtem a magnitude da estrela do catalogo gaia (a magnitude retornada
        # pelo praia esta corrigida por padrão por velocidade da sombra, por isso é substituída)
        minimum_distance = np.sqrt(
            (df_gaia_csv["ra"] - row["ra_star_deg"]) ** 2
            + (df_gaia_csv["dec"] - row["dec_star_deg"]) ** 2
        )
        star_index = np.argmin(minimum_distance)
        source_id, gaia_g_mag = (
            df_gaia_csv["source_id"][star_index],
            df_gaia_csv["phot_g_mean_mag"][star_index],
        )
        new_row.update({"gaia_source_id": source_id, "gaia_g_mag": gaia_g_mag})

        e_ra_star, e_dec_star = (
            df_gaia_csv["ra_error"][star_index],
            df_gaia_csv["dec_error"][star_index],
        )
        new_row.update({"e_ra_star": e_ra_star, "e_dec_star": e_dec_star})

        # ------------------------------------------------------------------------
        # Calcula a magnitude visual do asteroide no instante da ocultação
        # Alerta: os arquivos bsp estão na memoria global por alguma razão,
        #         convém analisar esse comportamente no futuro.
        # ------------------------------------------------------------------------
        ast_vis_mag = None
        if mag_cs:
            # usa intepolacao a partir de dados do JPL
            # mais preciso em casos como plutao etc..
            datetime_object = datetime.fromisoformat(row["date_time"])
            # Convert to Julian Date using astropy
            time_obj = Time(datetime_object)
            julian_date = time_obj.jd
            ast_vis_mag = mag_cs(julian_date)  #
        else:
            # tenta calcular a partir de h e g do bsp e só funciona para asteroides em geral
            if (
                object_data["h"] < 99
            ):  # some objects have h defined as 99.99 when unknown in the asteroid table inherited from MPC
                ast_vis_mag = asteroid_visual_magnitude(
                    asteroid_bsp = object_ephemeris,
                    naif_tls=leap_seconds,
                    planetary_bsp=planetary_ephemeris,
                    instant=datetime.strptime(row["date_time"], "%Y-%m-%d %H:%M:%S"),
                    h=object_data["h"],
                    g=object_data["g"],
                    spice_global=True,
                )

        # Calcula a queda em magnitude durante a ocultacao
        magnitude_drop = compute_magnitude_drop(ast_vis_mag, gaia_g_mag)

        # Calcula o diametro apararente o diametro em km existe
        apparent_diameter = get_apparent_diameter(
            object_data["diameter"], row["delta"]
        )

        # Calcula a duração do evento se o diametro existir
        event_duration = get_event_duration(object_data["diameter"], row["velocity"])

        # Calcula a separação angular da lua e do sol do objeto no instante da ocultação
        moon_separation, sun_elongation = get_moon_and_sun_separation(
            row["ra_target_deg"], row["dec_target_deg"], row["date_time"]
        )

        # Calcula a fração ilumnada da lua
        moon_illuminated_fraction = get_moon_illuminated_fraction(row["date_time"])

        # Obtem o valor da incerteza do objeto no instante da ocultação
        # e calcula a incerteza no instante central
        e_ra_target, e_dec_target = None, None
        closest_approach_uncertainty = None
        closest_approach_uncertainty_km = None
        instant_uncertainty = None
        if has_uncertainties:
            datetime_object = datetime.fromisoformat(row["date_time"])
            # Convert to Julian Date using astropy
            time_obj = Time(datetime_object)
            julian_date = time_obj.jd
            # Errors are divided by 3 because the interpolation is done using 3-sigma values
            e_ra_target = ra_cs(julian_date) / 3.0
            e_dec_target = dec_cs(julian_date) / 3.0

            instant_uncertainty = get_instant_uncertainty(
                row["position_angle"],
                row["delta"],
                row["velocity"],
                e_ra_target,
                e_dec_target,
                e_ra_star=df_gaia_csv["ra_error"][star_index] / 1000,
                e_dec_star=df_gaia_csv["dec_error"][star_index] / 1000,
            )

            closest_approach_uncertainty = get_closest_approach_uncertainty(
                row["position_angle"],
                e_ra_target,
                e_dec_target,
                e_ra_star=df_gaia_csv["ra_error"][star_index] / 1000,
                e_dec_star=df_gaia_csv["dec_error"][star_index] / 1000,
            )

            # TODO: se o angulo for menor que 0 ou maior que 90 graus a incerteza
            # tende ao infinito e por razoes praticas estao sendo defininas como None.
            # Verificar se essa é a melhor solução para o caso
            rad_angle = np.deg2rad(closest_approach_uncertainty / 3600) / 2

            if rad_angle > 0 and rad_angle < np.pi / 2:
                closest_approach_uncertainty_km = (
                    (row["delta"] * 149597870.7) * 2 * abs(np.tan(rad_angle))
                )
            else:
                closest_approach_uncertainty_km = None

        new_row.update(
            {
                "apparent_magnitude": ast_vis_mag,
                "magnitude_drop": magnitude_drop,
                "apparent_diameter": apparent_diameter,
                "event_duration": event_duration,
                "instant_uncertainty": instant_uncertainty,
                "closest_approach_uncertainty": closest_approach_uncertainty,
                "closest_approach_uncertainty_km": closest_approach_uncertainty_km,
                "moon_separation": moon_separation,
                "sun_elongation": sun_elongation,
                "moon_illuminated_fraction": moon_illuminated_fraction,
                "e_ra_target": e_ra_target,
                "e_dec_target": e_dec_target,
            }
        )

        occ_coeff = occultation_path_coeff(
            date_time=datetime.strptime(row["date_time"], "%Y-%m-%d %H:%M:%S")
            .replace(tzinfo=timezone.utc)
            .isoformat(),
            ra_star_candidate=row["ra_star_candidate"],
            dec_star_candidate=row["dec_star_candidate"],
            closest_approach=row["closest_approach"],
            position_angle=row["position_angle"],
            velocity=row["velocity"],
            delta_distance=row["delta"],
            offset_ra=row["off_ra"],
            offset_dec=row["off_dec"],
            closest_approach_error=closest_approach_uncertainty_km,
            object_diameter=row.get("diameter", None),
            object_diameter_error=row.get("diameter_err_max", None),
        )

        if (
            len(occ_coeff["coeff_latitude"]) > 0
            and len(occ_coeff["coeff_longitude"]) > 0
        ):
            new_row.update(
                {
                    "have_path_coeff": True,
                    "occ_path_min_longitude": (
                        float(occ_coeff["min_longitude"])
                        if occ_coeff["min_longitude"] != None
                        else None
                    ),
                    "occ_path_max_longitude": (
                        float(occ_coeff["max_longitude"])
                        if occ_coeff["max_longitude"] != None
                        else None
                    ),
                    "occ_path_min_latitude": (
                        float(occ_coeff["min_latitude"])
                        if occ_coeff["min_latitude"] != None
                        else None
                    ),
                    "occ_path_max_latitude": (
                        float(occ_coeff["max_latitude"])
                        if occ_coeff["max_latitude"] != None
                        else None
                    ),
                    "occ_path_is_nightside": bool(occ_coeff["nightside"]),
                    "occ_path_coeff": json.dumps(occ_coeff),
                }
            )

        coeff_paths.append(new_row)

    if len(coeff_paths) > 0:
        df_coeff = pd.DataFrame.from_dict(coeff_paths)
        df["gaia_source_id"] = df_coeff["gaia_source_id"]
        df["g_star"] = df_coeff["gaia_g_mag"]
        df["apparent_magnitude"] = df_coeff["apparent_magnitude"]
        df["magnitude_drop"] = df_coeff["magnitude_drop"]
        df["apparent_diameter"] = df_coeff["apparent_diameter"]
        df["event_duration"] = df_coeff["event_duration"]
        df["instant_uncertainty"] = df_coeff["instant_uncertainty"]
        df["closest_approach_uncertainty"] = df_coeff[
            "closest_approach_uncertainty"
        ]
        df["closest_approach_uncertainty_km"] = df_coeff[
            "closest_approach_uncertainty_km"
        ]
        df["moon_separation"] = df_coeff["moon_separation"]
        df["sun_elongation"] = df_coeff["sun_elongation"]
        df["moon_illuminated_fraction"] = df_coeff["moon_illuminated_fraction"]
        df["e_ra_target"] = df_coeff["e_ra_target"]
        df["e_dec_target"] = df_coeff["e_dec_target"]
        df["e_ra"] = df_coeff["e_ra_star"]
        df["e_dec"] = df_coeff["e_dec_star"]

        df["have_path_coeff"] = df_coeff["have_path_coeff"]
        df["occ_path_max_longitude"] = df_coeff["occ_path_max_longitude"]
        df["occ_path_min_longitude"] = df_coeff["occ_path_min_longitude"]
        df["occ_path_coeff"] = df_coeff["occ_path_coeff"]
        df["occ_path_is_nightside"] = df_coeff["occ_path_is_nightside"]
        df["occ_path_max_latitude"] = df_coeff["occ_path_max_latitude"]
        df["occ_path_min_latitude"] = df_coeff["occ_path_min_latitude"]

        del df_coeff


    # Converter as strings date_time do instante da ocultação em objetos datetime utc
    df["date_time"] = pd.to_datetime(df["date_time"], utc=True)

    #  Identificar as estrelas do gaia catalog que estao sendo utilizadas

    logger.info("Calculating unique HashIDs.")
    # Gerar o hash id
    # Gera um hash unico para cada evento de predicao
    df["hash_id"] = df.apply(
        lambda x: generate_hash(
            x["name"],
            x["gaia_source_id"],
            x["date_time"],
            x["ra_star_candidate"],
            x["dec_star_candidate"],
        ),
        axis=1,
    )

    logger.info("Writing the occultation table with path coefficients.")

    df.to_csv(occultation_table, index=False, sep=";")
    del df