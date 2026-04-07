from astroquery.mast import Catalogs

ps1id = 153851661137291254  # your PS1ID here

result = Catalogs.query_criteria(catalog="Panstarrs", objID=ps1id)

ra = result["raMean"][0]
dec = result["decMean"][0]

print(ra, dec)

import pandas as pd
import numpy as np
import os
from urllib.request import urlretrieve

from astroquery.mast import Catalogs
from astroquery.irsa import Irsa
from astropy.coordinates import SkyCoord
import astropy.units as u
import requests
from urllib.parse import urlencode

def query_ztf_metadata(ra, dec, radius_arcsec=2):

    radius_deg = radius_arcsec / 3600

    url = f"https://irsa.ipac.caltech.edu/ibe/search/ztf/products/sci?" + "CIRCLE,{ra},{dec},{radius_deg}"
    r = requests.get(url)
    print(r.text[:500])  # preview the first 500 characters
    df = pd.read_csv(url)
    df.to_csv("ASTRO501/seminar/ztf_cutouts_153851661137291254/request.csv")
    return df

def build_image_url(row):
    base = "https://irsa.ipac.caltech.edu/ibe/data/ztf/products/sci"
    return f"{base}/{row['filefracday']}/{row['filename']}"

def get_ztf_cutouts_from_zubercal(
    csv_path,
    output_dir="ASTRO501/seminar/ztf_cutouts_153851661137291254/",
    ps1id=153851661137291254,
    radius_arcsec=2,
    mjd_tolerance=0.01,
    max_images_per_point=1
):
    """
    Download ZTF image cutouts for a Zubercal light curve.

    Parameters
    ----------
    csv_path : str
        Path to Zubercal CSV file
    output_dir : str
        Where to save images
    radius_arcsec : float
        Search radius for ZTF images
    mjd_tolerance : float
        Time matching tolerance (days)
    max_images_per_point : int
        Limit number of images per MJD

    Returns
    -------
    None (downloads FITS files)
    """

    # --- Load light curve ---
    df = pd.read_csv(csv_path)

    if len(df) == 0:
        print("Empty CSV.")
        return

    # --- Resolve PS1ID → RA/Dec ---
    #ps1id = df["#PS1ID"].iloc[0]

    print(f"Resolving PS1ID {ps1id}...")

    result = Catalogs.query_criteria(catalog="Panstarrs", objID=ps1id)

    if len(result) == 0:
        print("Could not resolve PS1ID.")
        return

    ra = result["raMean"][0]
    dec = result["decMean"][0]

    print(f"RA, Dec = {ra}, {dec}")

    # --- Query ZTF images ---
    coord = SkyCoord(ra, dec, unit="deg")


    from astroquery.irsa import Irsa

    ra_deg = coord.ra.deg
    dec_deg = coord.dec.deg
    radius_deg = radius_arcsec / 3600

    print("Querying ZTF image archive (IBE)...")

    images = query_ztf_metadata(ra, dec, radius_arcsec)

    if len(images) == 0:
        print("No images found.")
        return

    print(f"Found {len(images)} candidate images")
    print(images.columns)
    matches = images[np.abs(images["mjd"] - mjd) < mjd_tolerance]
    # --- Prepare output ---
    os.makedirs(output_dir, exist_ok=True)

    zubercal_mjds = df["MJD"].values

    total_downloaded = 0

    # --- Match each MJD ---
    for i, mjd in enumerate(zubercal_mjds):

        # Find matching images in time
        mask = np.abs(images["mjd"] - mjd) < mjd_tolerance
        matches = images[mask]

        if len(matches) == 0:
            continue

        # Limit number of images
        matches = matches[:max_images_per_point]

        for j, row in enumerate(matches):

            # Construct download URL
            # (this field usually exists in IRSA tables)
            if "access_url" in row.colnames:
                url = row["access_url"]
            else:
                # fallback (less reliable)
                continue

            filename = f"{output_dir}/mjd_{mjd:.5f}_{i}_{j}.fits"

            try:
                urlretrieve(url, filename)
                total_downloaded += 1
                print(f"Downloaded: {filename}")
            except Exception as e:
                print(f"Failed: {e}")

    print(f"\nDone. Downloaded {total_downloaded} images.")
