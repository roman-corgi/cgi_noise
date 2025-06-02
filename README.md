# CGI_Noise

This package, developed by the California Institute of Technology, provides a function library and constants to compute the time to achieve a desired signal-to-noise ratio (tSNR) for exoplanet imaging with the Roman Space Telescope Coronagraph Instrument (CGI). It supports multiple observation scenarios defined by YAML files in the `data/Scenarios/` folder.

The methodology is documented in the SPIE JATIS paper:
> Nemati, B., Krist, J., Poberezhskiy, I., Kern, B., 2023. Analytical performance model and error budget for the Roman coronagraph instrument. Journal of Astronomical Telescopes, Instruments, and Systems 9. https://doi.org/10.1117/1.jatis.9.3.034007

## Requirements
- Python 3.12 or higher
- Dependencies (installed via `setup.py`):
  - `numpy>=1.20`
  - `pandas>=1.0`
  - `pyyaml>=5.4`
  - `xlwings>=0.24`
  - `chardet>=4.0`
  - `prettytable>=2.0`
- Excel installation (required for `xlwings` to load Excel files)
- The `data/` folder from the repository (included when cloning)

## Setup
1. Clone the repository:
   ```bash
   git clone https://github.com/roman-corgi/cgi_noise.git
   cd cgi_noise
   ```
2. Install dependencies using `setup.py`:
   ```bash
   pip install -e .
   ```
   The `-e` flag allows you to modify `cgi_noise/main.py` without reinstalling.
3. Ensure the `data/` folder is in the project root (parallel to `cgi_noise/`), included in the cloned repository. The structure is:
   ```
   data/
   ├── Scenarios/
   │   └── OPT_IMG_NFB1_HLC.yml
   │   └── REFERENCE_IMG_NFB1_HLC.yml
   │   └── CON_IMG_WFB4_SPC.yml
   │   └── ...
   ├── Photometry/
   │   └── *.csv
   ├── Calibration/
   │   └── *.csv
   ├── Cstability/
   │   └── *.csv
   └── Spectra/
       └── *.csv
   ```
   If running scripts from another directory, download `data/` from https://github.com/roman-corgi/cgi_noise or set the `CGI_NOISE_DATA_PATH` environment variable:
   ```bash
   export CGI_NOISE_DATA_PATH=/path/to/data
   ```

## Usage
The package is organized in the `cgi_noise/` subfolder, with `main.py` providing a complete example to compute tSNR. To use the package:

1. Run `main.py` with a scenario:
   ```bash
   python cgi_noise/main.py --scenario OPT_IMG_NFB1_HLC.yml
   ```
   Replace `OPT_IMG_NFB1_HLC.yml` with any YAML file from `data/Scenarios/` (e.g., `REFERENCE_IMG_NFB1_HLC.yml`, `CON_IMG_WFB4_SPC.yml`).

2. Create a custom script by copying `cgi_noise/main.py` and modifying it. Below is a simplified example:

```python
import os
from pathlib import Path
import yaml
from cgi_noise import check_data_path
from cgi_noise.cginoiselib import Target, getScenFileNames, loadCSVs, getSpectra
from cgi_noise.constants import mas

# Ensure data/ folder is accessible
data_path = check_data_path()
scenario_file = data_path / "Scenarios" / "OPT_IMG_NFB1_HLC.yml"  # Change to another .yml

# Step 1: Load scenario configuration
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
filename_list = getScenFileNames(config, data_path)
cg_data, qe_data, det_data, stray_data, thpt_data, cal_data, cs_data = loadCSVs(filename_list)

# Step 5: Calculate planet working angle
lam_d = wavelength / diameter
sep_mas = target.phaseAng_to_sep(target.sma_AU, target.dist_pc, target.phaseAng_deg)
planet_wa = sep_mas * mas / lam_d

# Step 6: Compute fluxes (simplified)
star_flux = getSpectra(target, wavelength, bandwidth)[2]
planet_flux = target.albedo_from_geomAlbedo(target.phaseAng_deg, target.geomAlb_ag) * star_flux

# Step 7: Compute time to desired SNR
# See cgi_noise/main.py for full implementation
desired_snr = 5.0
# time_to_snr = compute_tsnr(...)

print(f"Separation: {sep_mas:.0f} mas")
print(f"Star Flux: {star_flux:.3e} ph/s/m^2")
```

**Note**: The full tSNR calculation requires additional steps (e.g., coronagraph setup, photon rates, noise calculations) in `cgi_noise/main.py`. Copy and modify `main.py` to customize scenarios or parameters. Ensure the `data/` folder contains the required CSV files for the selected YAML.

## Documentation
- [JATIS Paper](https://doi.org/10.1117/1.jatis.9.3.034007)
- Source code: [GitHub Repository](https://github.com/roman-corgi/cgi_noise)

## License
BSD 3-Clause License, with restrictions on commercial use. See [LICENSE](LICENSE) for details. Contact Caltech’s Office of Technology Transfer for commercial inquiries.

6/2/2025