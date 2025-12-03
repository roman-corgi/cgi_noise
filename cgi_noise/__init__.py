import os
import importlib.resources

name = "cgi_noise"
__version__ = "1.3.0"

# identify data directory and add to environment variables for this session
datapath = importlib.resources.files("cgi_noise").joinpath("data")
assert datapath.exists(), (
    "Could not identify cgi_noise datapath. Check that the "
    "cgi_noise installation was fully successful."
)
os.environ["CGI_NOISE_DATA_DIR"] = str(datapath)
