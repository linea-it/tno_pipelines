import math

import spiceypy as spice
import pathlib
import logging

def findIDSPK(n, key):
    # OBS: Esta função está lendo o bsp do objeto direto da biblioteca spice.
    loc = 2  # order setting bsp files (1=DEXXX.bsp,2=Ast.bsp)
    m, header, flag = spice.dafec(loc, n)
    spk = ""
    for row in header:
        if row[: len(key)] == key:
            spk = row[len(key) :].strip()
    return spk

def angle(v1, v2):
    """
    Compute the angle between two vector (return the angle in degrees)
    """
    rad = math.acos(dotproduct(v1, v2) / (norm(v1) * norm(v2)))
    return math.degrees(rad)


def dotproduct(v1, v2):
    return sum((a * b) for a, b in zip(v1, v2))


def norm(v):
    return math.sqrt(dotproduct(v, v))

def HMS2deg(ra="", dec=""):
    RA, DEC, rs, ds = "", "", 1, 1
    if dec:
        D, M, S = [float(i) for i in dec.split()]
        if str(D)[0] == "-":
            ds, D = -1, abs(D)
        deg = D + (M / 60) + (S / 3600)
        DEC = deg * ds

    if ra:
        H, M, S = [float(i) for i in ra.split()]
        if str(H)[0] == "-":
            rs, H = -1, abs(H)
        deg = (H * 15) + (M / 4) + (S / 240)
        RA = deg * rs

    if ra and dec:
        return [RA, DEC]
    else:
        return RA or DEC

def ra2HMS(rarad=""):
    radeg = math.degrees(rarad)
    raH = int(radeg / 15.0)
    raM = int((radeg / 15.0 - raH) * 60)
    raS = 60 * ((radeg / 15.0 - raH) * 60 - raM)
    RA = "{:02d} {:02d} {:07.4f}".format(raH, raM, raS)
    return RA


def dec2DMS(decrad=""):
    decdeg = math.degrees(decrad)
    ds = "+"
    if decdeg < 0:
        ds, decdeg = "-", abs(decdeg)
    deg = int(decdeg)
    decM = abs(int((decdeg - deg) * 60))
    decS = 60 * (abs((decdeg - deg) * 60) - decM)
    DEC = "{}{:02d} {:02d} {:06.3f}".format(ds, deg, decM, decS)
    return DEC

def write_ephemeris_ascii_file(
        filepath: pathlib.Path, 
        dates: list, 
        ra: list, 
        dec: list, 
        distance: list, 
        elongation: list):
    
    n = len(dates)
    with open(filepath, "w") as outFile:
        outFile.write(
            "\n\n     Data Cal. UTC" + " ".ljust(51) + "R.A.__(ICRF//J2000.0)__DEC"
        )
        outFile.write(" ".ljust(43) + "DIST (km)" + " ".ljust(24) + "S-O-A\n")
        for i in range(n):
            outFile.write(dates[i] + " ".ljust(44) + ra[i] + "  " + dec[i] + " ".ljust(35))
            outFile.write("{:.16E}".format(distance[i]) + " ".ljust(17))
            outFile.write("{:.4f}".format(elongation[i]) + "\n")


def generate_ephemeris_file(
        dates_filepath: pathlib.Path, 
        object_ephemeris: pathlib.Path, 
        planetary_ephemeris: pathlib.Path, 
        leap_seconds: pathlib.Path, 
        cwd: str,
        logger: logging.Logger
    ):
    logger.info("Generating ephemeris file")

    logger.debug(f"Dates filepath: [{dates_filepath}]")
    logger.debug(f"Object Ephemeris: [{object_ephemeris}]")
    logger.debug(f"Planetary Ephemeris: [{planetary_ephemeris}]")
    logger.debug(f"Leap Seconds: [{leap_seconds}]")

    # Path for output file
    eph_filepath = pathlib.Path(cwd, "ephemeris.eph")

    logger.info("Load the planetary ephemeris, the leap second and asteroid bsp")
    # IMPORTANT: The order of the files is important
    spice.furnsh(str(planetary_ephemeris))
    spice.furnsh(str(leap_seconds))
    spice.furnsh(str(object_ephemeris))

    logger.info("Finding the IDSPK in the header of the bsp file")
    # Values specific for extract all comments of header from bsp files (JPL, NIMA)
    source = {"NIMA": (45, "ASTEROID_SPK_ID ="), "JPL": (74, "Target SPK ID   :")}
    n, key = source["NIMA"]
    idspk = findIDSPK(n, key)

    # TODO: Plutao nao segue as regras de bsp dos demais asteroides. Essa é uma solução temporária
    if object_ephemeris.name == "Pluto.bsp":
        idspk = "999"

    if idspk == "":
        n, key = source["JPL"]
        idspk = findIDSPK(n, key)

    logger.debug(f"IDSPK: [{idspk}]")

    # Read the file with dates
    logger.info("Reading the file with dates")
    with open(dates_filepath, "r") as inFile:
        dates = inFile.read().splitlines()

    n = len(dates)
    logger.debug(f"Number of dates: [{n}]")

    # Convert dates from utc to et format
    datesET = [spice.utc2et(utc) for utc in dates]

    # Compute geocentric positions (x,y,z) for each date with light time correction
    logger.info("Compute geocentric positions")
    rAst, ltAst = spice.spkpos(idspk, datesET, "J2000", "LT", "EARTH")
    rSun, ltSun = spice.spkpos("SUN", datesET, "J2000", "NONE", "EARTH")

    elongation = [angle(rAst[i], rSun[i]) for i in range(n)]

    data = [spice.recrad(xyz) for xyz in rAst]
    distance, rarad, decrad = zip(*data)

    # TODO: Revisar se este arquivo radec.txt está sendo utilizado
    # # ================= for graphics =================
    # radec_filepath = pathlib.Path(cwd, "radec.txt")    
    # radecFile = open(radec_filepath, "w")
    # for row in data:
    #     radecFile.write(str(row[1]) + ";" + str(row[2]) + "\n")
    # radecFile.close()

    # # Altera permissão do arquivo para escrita do grupo
    # # ================================================

    logger.info("Convert cartesian to angular coordinates")
    ra = [ra2HMS(alpha) for alpha in rarad]
    dec = [dec2DMS(delta) for delta in decrad]

    logger.info("Writing the ephemeris ascii file")
    write_ephemeris_ascii_file(eph_filepath, dates, ra, dec, distance, elongation)

    if eph_filepath.exists():
        logger.info("Ephemeris file generated successfully")
        logger.debug(f"Ephemeris filepath: [{eph_filepath}]")
        return eph_filepath
    else:
        raise Exception("Ephemeris file not generated.")