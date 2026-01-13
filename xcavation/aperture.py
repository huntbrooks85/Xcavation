
#-----------------------------------------------------------------------#
# xcavation.aperture v0.2.1
# By Hunter Brooks, at UToledo, Toledo: Jan. 13, 2026
#
# Purpose: Perform Aperture Photometry on SphereX Data
#-----------------------------------------------------------------------#



# Import Data Management
# ------------------------------------------------------ #
import numpy as np
import astropy.units as u
from astropy.io import fits
# ------------------------------------------------------ #



# Import WCS, Photometry, and Plotting
# ------------------------------------------------------ #
from astropy.wcs import WCS
from astropy.stats import SigmaClip
from astropy.nddata.utils import Cutout2D
from photutils.aperture import CircularAperture, CircularAnnulus, aperture_photometry, ApertureStats
# ------------------------------------------------------ #



# Import API Mangament Tools
# ------------------------------------------------------ #
import requests
from io import BytesIO
# ------------------------------------------------------ #



# Wavelength Resolving Power Table
# ------------------------------------------------------ #
def resolving_table(wave):
  if wave < 2.42:
    R = 41
  elif 2.42 <= wave < 3.82:
    R = 35
  elif 3.82 <= wave < 4.42:
    R = 110
  elif 4.42 <= wave:
    R = 130

  return wave/R
# ------------------------------------------------------ #




# Perform SphereX Wavelength and Aperature Calculation
# ------------------------------------------------------ #
def spherex_aperature_phot(url, coord, pixel_scale_deg=0.001708333):

    # ----- Download FITS entirely in RAM ----- #
    resp = requests.get(url)
    resp.raise_for_status()
    file_bytes = BytesIO(resp.content)
    # ----------------------------------------- #



    # Opens File
    with fits.open(file_bytes, memmap=False, lazy_load_hdus=True, ignore_missing_simple=True, checksum=False) as hdul:


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
        cutout_flux = Cutout2D(flux_image, position=coord, size=35, wcs=celestial_wcs)
        if cutout_flux.shape[0] != 35 or cutout_flux.shape[1] != 35:
          return {
                    "wavelength": np.nan,
                    "delta_lambda": np.nan,
                    "flux": np.nan,
                    "flux_err": np.nan,
                    "flux_cutout": np.nan,
                    "aperture": np.nan,
                    "annulus": np.nan,
                    "x_loc": np.nan,
                    "y_loc": np.nan,
                    "ap_radius": np.nan,
                    "inner_annulus": np.nan,
                    "outer_annulus": np.nan,
                    "url": url
                }

        cutout_var = Cutout2D(variance_image, position=coord, size=35, wcs=celestial_wcs)
        cutout_wcs = cutout_flux.wcs

        # Convert Flux and Variance Cutouts from MJyr/sr to uJy
        pixel_area_sr = (pixel_scale_deg * np.pi / 180)**2
        flux_ujy = np.nan_to_num(cutout_flux.data * 1e6 * pixel_area_sr * 1e6, nan=0.0)
        variance_ujy = np.nan_to_num(cutout_var.data * (1e6 * pixel_area_sr * 1e6)**2, nan=0.0)
        # ----------------------------------------- #



        # ----- Aperture Photometry With Local Background ----- #
        # Aperature Set-up
        x_query, y_query = cutout_wcs.world_to_pixel(coord) # Puts the RA and DEC into x, y Coord in the Cutout Frame
        fwhm = hdul[1].header['PSF_FWHM']

        r_fwhm = 2.5
        ap_radius_pix = (r_fwhm * fwhm * u.arcsec).to(u.deg).value / pixel_scale_deg # Aperature Radius
        aperture = CircularAperture((x_query, y_query), r=ap_radius_pix) # Set-up Aperature Matrix

        # Annulus Set-up
        an_inner = (45 * u.arcsec).to(u.deg).value / pixel_scale_deg # Inner Radius
        an_outer = (90 * u.arcsec).to(u.deg).value / pixel_scale_deg # Outer Radius
        annulus_ap = CircularAnnulus((x_query, y_query), r_in=an_inner, r_out=an_outer)

        # Calculates Background Frame
        sigclip = SigmaClip(sigma=5.0, maxiters=5) # Calculates the Mean and Standard Deviation of the Background Field
        bkg_per_pix = (ApertureStats(flux_ujy, annulus_ap, sigma_clip=sigclip)).median # Average Annulus Background
        aperture_area = aperture.area_overlap(flux_ujy) # Aperature Area in Pixels
        total_bkg = bkg_per_pix * aperture_area # The Total Background Frame in Aperature Area

        # Calculates the Aperature Photometry
        phot_table = aperture_photometry(flux_ujy, aperture)
        flux_ap = phot_table['aperture_sum'][0] - total_bkg # Aperature - Total Background Frame in Aperature Area

        # Calculate Progated Errors
        phot_var_table = aperture_photometry(variance_ujy, aperture)
        var_flux_ap = phot_var_table['aperture_sum'][0]
        var_bkg_per_pix = (((ApertureStats(flux_ujy, annulus_ap, sigma_clip=sigclip)).std)**2)/(annulus_ap.area_overlap(variance_ujy))
        flux_err = np.sqrt(var_flux_ap + (aperture_area**2)*var_bkg_per_pix)
        # ----------------------------------------------------- #



        # ----- Wavelength ----- #
        # General Set-Up
        y, x = np.indices(hdul["IMAGE"].data.shape) # Obtain Detector Size
        wave, *_ = spectral_wcs.pixel_to_world(x, y) # Convert the x,y Pixels to Wavelength Using Wavelength WCS
        wave = np.asarray(wave)

        # Obtains Wavelength Average
        yslice, xslice = cutout_flux.slices_original
        mask = aperture.to_mask(method='exact')
        wave_cutout = wave[yslice, xslice] # Cutout the Wavelength Cutout Used By Aperature (FIX LATER)
        flux_cutout = flux_image[yslice, xslice] # Cutout the Flux Cutout Used By Aperature (FIX LATER)
        mask = np.isfinite(wave_cutout) & np.isfinite(flux_cutout) & (flux_cutout > 0) # Enough Valid Points
        wave_point = np.average(wave_cutout[mask], weights=flux_cutout[mask]) # Average the Cutouts and Weigh by Flux
        # ---------------------- #

    ap_mask_obj = aperture.to_mask(method='exact')
    an_mask_obj = annulus_ap.to_mask(method='exact')

    # Convert to image-sized arrays
    ap_mask = ap_mask_obj.to_image(cutout_flux.data.shape)
    an_mask = an_mask_obj.to_image(cutout_flux.data.shape)

    # Return Output Dictionary
    return {
        "wavelength": wave_point,
        "delta_lambda": resolving_table(wave_point),
        "flux": flux_ap,
        "flux_err": flux_err,
        "flux_cutout": flux_ujy,
        "aperture": ap_mask,
        "annulus": an_mask,
        "x_loc": x_query,
        "y_loc": y_query,
        "ap_radius": r_fwhm * fwhm,
        "inner_annulus": 45,
        "outer_annulus": 90,
        "url": url
    }
# ------------------------------------------------------ #
