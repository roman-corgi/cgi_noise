# CGI_Noise

This package contains a function library and associated constants that describe the noise properties of the Roman Space Telescope Coronagraph Instrument (CGI) when conducting Direct Imaging (DI) or Spectroscopic observations.

The formulae and approach are documented in the SPIE JATIS paper:
> Nemati, B., Krist, J., Poberezhskiy, I., Kern, B., 2023. Analytical performance model and error budget for the Roman coronagraph instrument. Journal of Astronomical Telescopes, Instruments, and Systems 9. https://doi.org/10.1117/1.jatis.9.3.034007

## Requirements
- Python 3.12 or higher (required for compatibility with related Roman Space Telescope tools)
- Required dependencies (installed automatically): `numpy`, `pandas`, `pyyaml`, `xlwings`, `chardet`, `prettytable`
- Excel installation (required for `xlwings` to load Excel files)

## Installation
Clone the repository and install the package:
```bash
git clone https://github.com/roman-corgi/cgi_noise.git
cd cgi_noise
pip install .
```

Alternatively, install directly from GitHub:
```bash
pip install git+https://github.com/roman-corgi/cgi_noise.git
```

## Data Files
The package requires CSV and YAML data files located in a `data/` folder, included in the GitHub repository. After cloning, ensure the `data/` folder is in your working directory or specify its path in your scripts. The structure is:
```
data/
├── Scenarios/
│   └── OPT_IMG_NFB1_HLC.yml
├── Photometry/
│   └── *.csv
├── Calibration/
│   └── *.csv
├── Cstability/
│   └── *.csv
└── Spectra/
    └── *.csv
```

Note: The installed package includes example data files (e.g., `site-packages/cgi_noise/data/Scenarios/OPT_IMG_NFB1_HLC.yml`), but the full dataset is available in the repository’s `data/` folder.

## Usage
The `CGI_Noise` package is organized under the `cgi_noise` module, with all Python files (e.g., `cginoiselib.py`, `main.py`) in the `cgi_noise/` subfolder. It involves multiple steps to simulate exoplanet imaging, as demonstrated in `cgi_noise/main.py`. Below is a simplified example to calculate the time required to achieve a desired signal-to-noise ratio (SNR) for an exoplanet observation. See `cgi_noise/main.py` for a complete implementation.

```python
import os
from pathlib import Path
import yaml
from cgi_noise.cginoiselib import Target, getScenFileNames, loadCSVs, getSpectra
from cgi_noise.unitsConstants import mas

# Ensure data/ folder is accessible
data_path = Path.cwd() / "data"
if not data_path.exists():
    raise FileNotFoundError("Download the 'data/' folder from https://github.com/roman-corgi/cgi_noise")

# Step 1: Load scenario configuration
scenario_file = data_path / "Scenarios" / "OPT_IMG_NFB1_HLC.yml"
with open(scenario_file, "r") as f:
    config = yaml.safe_load(f)

# Step 2: Define instrument parameters
wavelength = config['instrument']['wavelength']
diameter = config['instrument']['Diam']
bandwidth = config['instrument']['bandwidth']

# Step 3: Define target star and planet
target = Target(
    v_mag=5.0,
    dist_pc=10.0,
    specType='g0v',
    phaseAng_deg=65.0,
    sma_AU=4.1536,
    radius_Rjup=5.6211,
    geomAlb_ag=0.44765,
    exoZodi=1,
)

# Step 4: Load CSV data files
filename_list = getScenFileNames(config)
cg_data, qe_data, det_data, stray_data, thpt_data, cal_data, cs_data = loadCSVs(filename_list)

# Step 5: Calculate planet working angle
lam_d = wavelength / diameter
sep_mas = target.phaseAng_to_sep(target.sma_AU, target.dist_pc, target.phaseAng_deg)
planet_wa = sep_mas * mas / lam_d  # Convert mas to lambda/D

# Step 6: Compute fluxes (simplified)
star_flux = getSpectra(target, wavelength, bandwidth)[2]
planet_flux = target.albedo_from_geomAlbedo(target.phaseAng_deg, target.geomAlb_ag) * star_flux

# Step 7: Compute time to desired SNR
# Requires additional setup (coronagraph parameters, noise calculations) from main.py
desired_snr = 5.0
# time_to_snr = compute_tsnr(...)  # See cgi_noise/main.py for full implementation

print(f"Separation: {sep_mas:.0f} mas")
print(f"Star Flux: {star_flux:.3e} ph/s/m^2")
```

**Note**: The full simulation involves additional steps (e.g., coronagraph setup, photon rates, noise calculations) detailed in `cgi_noise/main.py`. Copy or adapt `main.py` for your specific use case, ensuring the `data/` folder is correctly configured. All package functionality is accessed via the `cgi_noise` module (e.g., `cgi_noise.cginoiselib.Target`).

## Documentation
- [JATIS Paper](https://doi.org/10.1117/1.jatis.9.3.034007)
- API documentation: [link to Sphinx docs once hosted, e.g., GitHub Pages]

## Contributing
See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines on submitting issues or pull requests.

## License
BSD 3-Clause License, with restrictions on commercial use. See [LICENSE](LICENSE) for details. Contact Caltech’s Office of Technology Transfer for commercial inquiries.

6/2/2025