from pathlib import Path

import urllib.request

import os
import shutil
from zipfile import ZipFile
import tempfile

try:
    DATADIR = Path(__file__).parent
except NameError:
    DATADIR = Path("./")

"""
Bioclimatic variables:

BIO1 = Annual Mean Temperature
BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
BIO3 = Isothermality (BIO2/BIO7) (* 100)
BIO4 = Temperature Seasonality (standard deviation *100)
BIO5 = Max Temperature of Warmest Month
BIO6 = Min Temperature of Coldest Month
BIO7 = Temperature Annual Range (BIO5-BIO6)
BIO8 = Mean Temperature of Wettest Quarter
BIO9 = Mean Temperature of Driest Quarter
BIO10 = Mean Temperature of Warmest Quarter
BIO11 = Mean Temperature of Coldest Quarter
BIO12 = Annual Precipitation
BIO13 = Precipitation of Wettest Month
BIO14 = Precipitation of Driest Month
BIO15 = Precipitation Seasonality (Coefficient of Variation)
BIO16 = Precipitation of Wettest Quarter
BIO17 = Precipitation of Driest Quarter
BIO18 = Precipitation of Warmest Quarter
BIO19 = Precipitation of Coldest Quarter
"""

urls = [
    # Bioclimatic variables, averages 1970–2000
    "http://biogeo.ucdavis.edu/data/worldclim/v2.0/tif/base/wc2.0_30s_bio.zip",
    # Average monthly precipitation
    "http://biogeo.ucdavis.edu/data/worldclim/v2.0/tif/base/wc2.0_30s_prec.zip",

    # GLOBE
    "https://www.ngdc.noaa.gov/mgg/topo/DATATILES/elev/all10g.zip",
]

for url in urls:
    # Download the file from `url` and save it locally under `file_name`:
    with urllib.request.urlopen(url) as response:
        print("Received {:}, downloading…".format(url))
        handle, file_name = tempfile.mkstemp(dir=DATADIR, suffix=".zip")
        with os.fdopen(handle, "wb") as temp_file:
            shutil.copyfileobj(response, temp_file)
        print("Done.")
        print("Extracting…")
        with ZipFile(file_name, 'r') as zip_file:
            zip_file.extractall()
        print("Done.")
