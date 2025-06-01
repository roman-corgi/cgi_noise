# CGI_Noise

This package contains a function library and associated constants that describe the noise properties of the Roman Space Telescope Coronagraph Instrument (CGI) when conduction Direct Imaging (DI) or Spectroscopic observations.  

The formulae and approach are documented in the SPIE JATIS paper on the CGI error budget:

> Nemati, B., Krist, J., Poberezhskiy, I., Kern, B., 2023. Analytical performance model and error budget for the Roman coronagraph instrument. Journal of Astronomical Telescopes, Instruments, and Systems 9. https://doi.org/10.1117/1.jatis.9.3.034007

The accompanying script (Main.py) details the steps needed to calculate the photon rates and the noise (variance) rates required for the computation of the integration time for a desired SNR on a given target. 

Roman CGI has a number of discrete coronagraph operating modes and for each of these an "optimistic" and "conservative" configuration file is provided in the **scenarios** folder. 

 

5/31/2025
