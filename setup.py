from setuptools import setup, find_packages

setup(
    name="cgi_noise",
    version="0.1.0",
    description="SNR estimation tool for exoplanet direct imaging",
    author="roman-corgi",
    url="https://github.com/roman-corgi/cgi_noise",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "numpy",
        "pandas",
        "PyYAML",
        "prettytable",
        "chardet"
    ],
    python_requires='>=3.12',
    entry_points={
        'console_scripts': [
            'run_cgi_snr = main:run_snr_scenario'
        ]
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
