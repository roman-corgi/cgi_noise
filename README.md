# CGISNR Estimation Tool for Exoplanet Imaging

This package contains a function library and associated constants that describe the noise properties of the Roman Space Telescope Coronagraph Instrument (CGI) when conducting Direct Imaging (DI) or Spectroscopic observations.  

The formulae and approach are documented in the SPIE JATIS paper on the CGI error budget:

> Nemati, B., Krist, J., Poberezhskiy, I., Kern, B., 2023. Analytical performance model and error budget for the Roman coronagraph instrument. Journal of Astronomical Telescopes, Instruments, and Systems 9. https://doi.org/10.1117/1.jatis.9.3.034007

This Python package calculates the required integration time to achieve a target signal-to-noise ratio (SNR) in direct imaging of exoplanets. It models star and planet flux, coronagraph throughput, detector noise, and background noise (zodi, speckle, stray light).

## 📦 Features
- YAML-based scenario configuration
- No hardcoded paths
- Python 3.12+ required
- Ready for command-line or notebook use

## 📁 Folder Structure
```
cgi_noise/
├── main.py                       # Entry point
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

# Run the default scenario
$ python main.py

# Run a different config file
$ python main.py --config CON_SPEC_NFB3_SPC.yml
```

## 🧪 Customizing Targets and SNR
To use custom target parameters or SNR thresholds, edit `main.py` before the `run_pipeline()` follow this example:

```python
target_params = {
    'v_mag': 5.0,
    'dist_pc': 10.0,
    'specType': 'g0v',
    'phaseAng_deg': 65,
    'sma_AU': 4.1536,
    'radius_Rjup': 5.6211,
    'geomAlb_ag': 0.44765,
    'exoZodi': 1,
}
SNRdesired = 10.0  # instead of the default 5.0
run_pipeline(config, DATA_DIR, target_params, SNRdesired)
```

## 🔧 Requirements
- Python 3.12 or higher
- Dependencies listed in `requirements.txt`

## 
