

# WORK IN PROGRESS FILE - DO NOT SUGGEST ADDITIONS
# USE APERTURE.PY INSTEAD FOR NOW

#-----------------------------------------------------------------------#
# xcavation.psf v0.3.0
# By Hunter Brooks, at UToledo, Toledo: Jan. 13, 2026
#
# Purpose: Perform PSF Photometry on SphereX Data
#-----------------------------------------------------------------------#



# Import Data Management
# ------------------------------------------------------ #
import numpy as np
from astropy.io import fits
# ------------------------------------------------------ #



# Import WCS, Photometry, and Plotting
# ------------------------------------------------------ #
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from astropy.nddata.utils import Cutout2D
from photutils.detection import DAOStarFinder
from photutils.psf import PSFPhotometry, CircularGaussianSigmaPRF
# ------------------------------------------------------ #



# Import API Mangament Tools
# ------------------------------------------------------ #
import requests
from io import BytesIO
# ------------------------------------------------------ #



# Wavelength Resolving Power Table
# ------------------------------------------------------ #
def resolving_table(wave):
  if 0.75 <= wave < 2.42:
    R = 41
  elif 2.42 <= wave < 3.82:
    R = 35
  elif 3.82 <= wave < 4.42:
    R = 110
  elif 4.42 <= wave <= 5:
    R = 130
  return wave/R
# ------------------------------------------------------ #



# Perform SphereX PSF Photometry
# ------------------------------------------------------ #
def spherex_psf_phot(url, coord, pixel_scale_deg=0.001708333):
  # ----- Download FITS entirely in RAM ----- #
  resp = requests.get(url)
  resp.raise_for_status()
  file_bytes = BytesIO(resp.content)
  # ----------------------------------------- #



  # Opens File
  with fits.open(url, lazy_load=True) as hdul:



    # ---------------- WCS ---------------- #
    # Make sure the SIP warning message doesnt show
    spectral_header = hdul["IMAGE"].header.copy()
    for key in list(spectral_header.keys()):
        if key.startswith(("A_", "B_", "AP_", "BP_")):
            del spectral_header[key]

    # Obtain Header Data
    celestial_wcs = WCS(hdul[1].header.copy()) # Normal WCS
    spectral_wcs = WCS(spectral_header, fobj=hdul, key="W") # Wavelength WCS
    flux_image = hdul[1].data.astype(float) - hdul[4].data.astype(float) # IMAGE Image - ZODI Image
    variance_image = hdul['VARIANCE'].data.astype(float) # VARIANCE Image
    # ------------------------------------- #



    # ---------------- Cutouts ---------------- #
    # Cutout of Flux, Variance and WCS Using Inputted Radius
    cutout_flux = Cutout2D(flux_image, position=coord, size=20, wcs=celestial_wcs)
    cutout_var = Cutout2D(variance_image, position=coord, size=20, wcs=celestial_wcs)
    cutout_wcs = cutout_flux.wcs

    # Convert Flux and Variance Cutouts from MJyr/sr to uJy
    pixel_area_sr = (pixel_scale_deg * np.pi / 180)**2
    flux_ujy = cutout_flux.data * 1e6 * pixel_area_sr * 1e6
    variance_ujy = cutout_var.data * (1e6 * pixel_area_sr * 1e6)**2
    uncertainty_ujy = np.sqrt(variance_ujy)
    # ----------------------------------------- #



    # ---------------- PSF Photometry ---------------- #
    data = flux_ujy - np.nanmedian(flux_ujy)
    error = uncertainty_ujy

    # PSF parameters
    fwhm = hdul[1].header['PSF_FWHM']/6.15
    sigma = fwhm / 2.355
    bkg_rms = np.nanstd(data)
    finder = DAOStarFinder(5.0 * bkg_rms, fwhm)

    # Pixel-integrated PSF (non-deprecated)
    psf_model = CircularGaussianSigmaPRF(flux=1.0, sigma=sigma)
    psf_model.sigma.fixed = False
    fit_shape = (11, 11)
    psfphot = PSFPhotometry( psf_model, fit_shape, finder=finder, aperture_radius=4)

    phot = psfphot(data, error=error)
    if phot is None:
      flux_point = np.nan
      error_point = np.nan
    elif phot is not None:
      best_idx = np.nanargmin(phot['reduced_chi2'])
      flux_point = phot['flux_fit'][best_idx]
      error_point = phot['flux_err'][best_idx]
    # ------------------------------------------------ #



    # ----- Wavelength ----- #
    # General Set-Up
    y, x = np.indices(hdul["IMAGE"].data.shape) # Obtain Detector Size
    wave, *_ = spectral_wcs.pixel_to_world(x, y) # Convert the x,y Pixels to Wavelength Using Wavelength WCS
    wave = np.asarray(wave)

    # Obtains Wavelength Average
    yslice, xslice = cutout_flux.slices_original
    wave_cutout = wave[yslice, xslice] # Cutout the Wavelength Cutout Used By Aperature (FIX LATER)
    flux_cutout = flux_image[yslice, xslice] # Cutout the Flux Cutout Used By Aperature (FIX LATER)
    mask = np.isfinite(wave_cutout) & np.isfinite(flux_cutout) & (flux_cutout > 0) # Enough Valid Points
    wave_point = np.average(wave_cutout[mask], weights=flux_cutout[mask]) # Average the Cutouts and Weigh by Flux
    # ---------------------- #

    return {
        "wavelength": wave_point,
        "flux": flux_point,
        "flux_err": error_point,
        "flux_cutout": data,
        # "aperture": ap_mask,
        # "annulus": an_mask,
        # "x_loc": x_query,
        # "y_loc": y_query,
        "url": url
    }
# ------------------------------------------------------ #
