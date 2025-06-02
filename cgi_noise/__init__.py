__version__ = "1.0.0"

import os
from pathlib import Path

def check_data_path():
    """Check if the data/ folder exists in the project root or CGI_NOISE_DATA_PATH."""
    # Check environment variable
    env_data_path = os.environ.get("CGI_NOISE_DATA_PATH")
    if env_data_path:
        data_path = Path(env_data_path)
        if data_path.exists():
            return data_path
    
    # Check project root (parallel to cgi_noise/)
    project_root = Path(__file__).parent.parent  # CGI_Noise/
    data_path = project_root / "data"
    if data_path.exists():
        return data_path
    
    raise FileNotFoundError(
        "The 'data/' folder is missing. Please download it from "
        "https://github.com/roman-corgi/cgi_noise and place it in the project root "
        "(parallel to cgi_noise/), or set the CGI_NOISE_DATA_PATH environment variable."
    )