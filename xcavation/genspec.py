
#-----------------------------------------------------------------------#
# xcavation.genspec v0.3.0
# By Hunter Brooks, at UToledo, Toledo: Jan. 13, 2026
#
# Purpose: Main API Function for SphereX Data Retrieval and Photometry
#-----------------------------------------------------------------------#



# Import Internal Modules
# ------------------------------------------------------ #
from xcavation.motion import proper_motion, time
from xcavation.aperture import spherex_aperature_phot
from xcavation.psf import spherex_psf_phot
from xcavation.quality import finder_chart, spectra_plot
# ------------------------------------------------------ #



# Import Data Management
# ------------------------------------------------------ #
import numpy as np
import pandas as pd
import astropy.units as u
from astroquery.irsa import Irsa
# ------------------------------------------------------ #



# Import WCS, Photometry, and Plotting
# ------------------------------------------------------ #
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
# ------------------------------------------------------ #



# Import Multithreading and Watching Packages
# ------------------------------------------------------ #
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor
# ------------------------------------------------------ #




# Multi-Thread the SphereX Query Function
# ------------------------------------------------------ #
def genspec(ra, dec, style, pmra=None, pmdec=None, mjd = None, verification = False, threads = 8):

    # ----- Proper Motion Propagation ----- #
    if pmra is not None and pmdec is not None:
      time_passed = time(mjd) # years since CW observation
      print(f"\nTime Passed Since Observation: {round(time_passed, 6)} years")
      ra_deg, dec_deg = proper_motion(ra, dec, pmra, pmdec, time_passed) # Adjust for Proper Motion
      print(f"Adjusted Coordinates: RA = {round(ra_deg, 6)} deg, Dec = {round(dec_deg, 6)} deg") # Print Adjusted Coordinates
    else:
      ra_deg = ra
      dec_deg = dec
    # ------------------------------------- #



    # ----- SphereX Query ----- #
    coord = SkyCoord(ra_deg, dec_deg, unit='deg') # Convert to astropy Units
    results = Irsa.query_sia(pos=(coord, 10 * u.arcsec), collection='spherex_qr2') # Query SphereX QR2 (Will Need Updates)
    urls = results["access_url"] # Gets the Image URLs
    print(f'\n{len(urls)} Total Images Found in SphereX')
    # ------------------------- #



    # ----- Multi-Thread Wavelength and Aperature Calc. ----- #
    output = []
    with ThreadPoolExecutor(max_workers=threads) as executor:

        print('\nPerforming Photometric Calculations')

        # Multi-Thread w/ tqdm
        if style == 'aperture':
            futures = [executor.submit(spherex_aperature_phot, url, coord) for url in urls]
        elif style == 'psf':
            futures = [executor.submit(spherex_psf_phot, url, coord) for url in urls]
        else:
            return np.nan, np.nan, np.nan

        for future in tqdm(futures):
            result = future.result()

            output.append({
                "wavelength": result["wavelength"],
                "delta_lambda": result["delta_lambda"],
                "flux": result["flux"],
                "flux_err": result["flux_err"],
                "flux_cutout": result["flux_cutout"],
                "aperture": result["aperture"],
                "annulus": result["annulus"],
                "x_loc": result["x_loc"],
                "y_loc": result["y_loc"],
                "ap_radius": result['ap_radius'],
                "inner_annulus": result['inner_annulus'],
                "outer_annulus": result['outer_annulus'],
                "url": result["url"]
            })
    # ------------------------------------------------------- #



    # ----- Finder Chart and Spectral Plots ----- #
    if verification is True:
      finder_chart(output)
      spectra_plot(output)
    # # ------------------------------------------- #



    # Returns The Output Dictionary
    return pd.DataFrame(output)
# ------------------------------------------------------ #
