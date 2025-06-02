# Roman CGI Noise Estimation Tool

# for Exoplanet Imaging

This package contains a function library and associated constants that describe the noise properties of the Roman Space Telescope Coronagraph Instrument (CGI) when conducting Direct Imaging (DI) or Spectroscopic observations.  

The formulae and approach are documented in the SPIE JATIS paper on the CGI error budget:

> Nemati, B., Krist, J., Poberezhskiy, I., Kern, B., 2023. Analytical performance model and error budget for the Roman coronagraph instrument. Journal of Astronomical Telescopes, Instruments, and Systems 9. https://doi.org/10.1117/1.jatis.9.3.034007

This Python package calculates the required integration time to achieve a target signal-to-noise ratio (SNR) in direct imaging of exoplanets. It models star and planet flux, coronagraph throughput, detector noise, and background noise (zodi, speckle, stray light).

## 📦 Features
- Modular architecture
- YAML-based scenario configuration
- No hardcoded paths
- Python 3.10+ required
- Ready for command-line or notebook use

## 📁 Folder Structure
```
cgi_noise/
├── Main.py                       # Entry point
├── setup.py                      # Installation script
├── requirements.txt              # Python dependencies
├── README.md                     # This file
└── cgi_noise/
    ├── __init__.py
    ├── main_core.py              # Core simulation pipeline
    ├── cginoiselib.py            # Physics and modeling logic
    ├── unitsConstants.py         # Physical units and constants
    ├── loadCSVrow.py             # CSV loader with comment support
    └── data/
        ├── Scenarios/            # YAML configuration files
        │   └── OPT_IMG_NFB1_HLC.yml
        ├── Photometry/
        ├── Calibration/
        ├── Cstability/
        └── Spectra/
```
Upon import, `cgi_noise` creates a new environment variable (`CGI_NOISE_DATA_DIR`) that points to the location of the data directory wherever it is installed. 

## 🚀 Quick Start
```bash
# Clone the repository
$ git clone https://github.com/roman-corgi/cgi_noise.git
$ cd cgi_noise

# Create virtual environment (optional)
$ python -m venv venv
$ source venv/bin/activate  # or venv\Scripts\activate on Windows

# Install in editable mode
$ pip install -e .

# Run the script
$ python Main.py
```

## 🧪 Customizing Targets and SNR
To use a custom observation scenario, SNR threshold, and target parameters, define them as a dictionary at the top of `main.py`:

```python
obs_params = {
    "scenario": "OPT_IMG_NFB1_HLC.yml",
    "target_params": {
        'v_mag': 5.0,
        'dist_pc': 10.0,
        'specType': 'g0v',
        'phaseAng_deg': 65,
        'sma_AU': 4.1536,
        'radius_Rjup': 5.6211,
        'geomAlb_ag': 0.44765,
        'exoZodi': 1,
    },
    "snr": 5.0
}

run_snr_scenario(obs_params)
```

This keeps all user-modifiable values in one place, improving readability and flexibility.

## 🔧 Requirements
- Python 3.10 or higher
- Dependencies listed in `requirements.txt`

