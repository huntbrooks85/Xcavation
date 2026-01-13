
<p align="center">
    <a href="https://ibb.co/pvx3rgVs"><img src="/misc/xcavation.png" width="60%"></a> <br>
</p>

[![SPHEREx](https://img.shields.io/badge/SPHEREx-QR2-green)](https://irsa.ipac.caltech.edu/data/SPHEREx/docs/overview_qr.html)

[![PyPi](https://img.shields.io/pypi/v/xcavation)](https://pypi.org/project/xcavation/)

[![PyPi - Downloads](https://img.shields.io/pypi/dm/xcavation)]

<h1 align="center" id="title"> 🪏 <b> Xcavation </b> 🪏 </h1>

<div align="center">
  <p id="description"> <b>Xcavation</b> is a Python package designed for efficient retrieval, photometric extraction, and visualization of <a href="https://spherex.caltech.edu/">SPHEREx</a> survey data. Streamlined API for querying SPHEREx QR2 data products via IRSA, performing either aperture or PSF photometry, accounting for proper motion. Multi-threading is supported to significantly accelerate polar regions.  </p>
</div>

<div align="center">
  <pp><b> It is recommended that a user uses Google Colab, as it has high-speed internet for free. <a href="https://drive.google.com/file/d/1S0iLap2IoNdk4ErpdEyMyipBvJeaJ3ay/view?usp=sharing">Here</a> is an example. </b></pp> 
</div>

<div align="center">
  <h2>🛠️ Installation 🛠️</h2>
</div>

<div align="center">
<pp><b> pip Installation </b><pp>
</div>
<div align="center">
</div>

1. **Download Python:** Visit [here](https://www.python.org/downloads/) to install Python 
2. **Download pip:** Visit [here](https://pip.pypa.io/en/stable/installation/) to install pip
3. **Run Install Command:** Run the command in terminal:
   ```bash
   pip install xcavation

<div align="center">
  <h2>⚙️ Using Xcavation ⚙️</h2>
</div>

<div align="center">
  <p><b> How to Use genspec </b></p>
</div>
<div align="center">
</div>

1. After Xcavation is installed, verify the installation by running the following command: ```from xcavation.genspec import *```. If you encounter any issues during installation, please reach out to Hunter Brooks for assistance. 
2. Assign the relavent variables as described below. 
3. Execute the command: ```genspec(ra, dec, style)```. These are the minimum required parameters for Xcavation to run. You can include optional variables if needed.


<div align="center">
  <pp><b> Relavent Variables For genspec </b></pp> 
</div>

- **Required Variables:**
  - **ra:** Right Accension in Degrees: *float*:
     - *example:* ```131.123```

  - **dec:** Declination in Degrees: *float*:
     - *example:* ```-12.31254```

  - **style:** Aperature or PSF Photometry: *string*:
     - *example:* ```aperture``` or ```psf```

- **Optional Variables:**
  - **pmra:** Proper Motion in Right Accension (in arcsec/year): *float*
    - *example:* ```-0.981```, default=```None```
  - **pmdec:** Proper Motion in Declination (in arcsec/year): *float*
    - *example:* ```0.123```, default=```None```
  - **mjd:** Modified Julian Date of inputed R.A. and Decl. from Above: *float*
    - *example:* ```57170```, default=```None```
  - **verification:** Plots Q.A. and Spectrum if True: *boolean*
    - *example:* ```True```, default=```False```
  - **threads:** Number of Threads for Multi-Threading: *list*
    - *example:* ```2```, default=```8```

<div align="center">
  <h2>📞 Support & Development Team 📞</h2>
</div>

- **Mr. Hunter Brooks**
  - Email: hbrooks8@rockets.utoledo.edu

<div align="center">
  <h2>📖 Acknowledgments 📖</h2>
</div>

1. If you intend to publish any calculations done by fossyl, please reference Brooks et al. (in prep.).

2. Please reference the relavent SPHEREx publication.

