import os
from pathlib import Path
import cgi_noise
from cgi_noise.loadCSVrow import loadCSVrow
import cgi_noise.cginoiselib as fl

D = Path(os.environ["CGI_NOISE_DATA_DIR"]) / "Photometry"
for name in ["DET_CBE_250530", "DET_CBE_261116"]:
    det = loadCSVrow(str(D / f"{name}.csv"))
    r = det.df.iloc[0]
    n = fl.detector_noise_rates(det, 21, 1.0, 1, True)   # per pixel, 1 s frame
    print(f"{name}: EMGain={r.EMGain:.0f}  RN={r.ReadNoise_e} e-  "
          f"dark@21mo={n.dark*3600:.2f} e-/pix/hr  CIC@21mo={n.CIC:.5f} e-/pix/fr  "
          f"FWC img/ser={r.FWCimage:.0f}/{r.FWCserial:.0f}")