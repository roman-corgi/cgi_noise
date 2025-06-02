from setuptools import setup, find_packages

setup(
    name="cgi_noise",
    version="1.0.0",
    description="Noise modeling for the Roman Space Telescope Coronagraph Instrument",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="California Institute of Technology",
    author_email="your.email@caltech.edu",  # Replace with appropriate contact
    url="https://github.com/your-org/cgi_noise",  # Replace with your actual GitHub repo URL
    packages=find_packages(),
    include_package_data=True,
    package_data={
        "cgi_noise": [
            "data/Scenarios/*.yml",
            "data/Photometry/*.csv",
            "data/Calibration/*.csv",
            "data/Cstability/*.csv",
            "data/Spectra/*.csv",
        ],
    },
    install_requires=[
        "numpy>=1.20",
        "pandas>=1.0",
        "pyyaml>=5.4",
        "xlwings>=0.24",
        "chardet>=4.0",
        "prettytable>=2.0",
    ],
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
        "Topic :: Scientific/Engineering :: Astronomy",
    ],
    python_requires=">=3.12",
    entry_points={
        "console_scripts": [
            "cgi-noise=cgi_noise.main:main",
        ],
    },
)