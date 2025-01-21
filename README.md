# SIBPF_2.0
The current code identifies steep-spectrum radio sources in an image. It then cross-matches these sources with the TGSS (TIFR GMRT Sky Survey) and NVSS (National Radio Astronomy Observatory's VLA Sky Survey) catalogs to determine their spectral index. The script calculates spectral indices for matched sources and identifies potential pulsar candidates. The current version has been tailored for radio images from GMRT at 320 MHz. 
The new version of the code is designed to be more flexible, with improvements allowing it to handle a wider range of data, frequency bands, and source characteristics, making it adaptable for future use across various surveys and radio frequencies. 

## Requirements
To run the script, make sure you have the following prerequisites installed:

- Python 
- Required Python packages: pandas, astropy, numpy, matplotlib, scipy, psrqpy, astroquery, bdsf
  
## Download Data

Before running the script, you need to download the following data files:

1. Download TGSS Catalog (TGSSADR1_7sigma_catalog.tsv):
   [Download TGSS Catalog](http://tgssadr.strw.leidenuniv.nl/catalogs/TGSSADR1_7sigma_catalog.tsv)

2. Download Spectral Index Catalog (spidxcat_v1.1b.fits):
   [Details about the catalog](https://tgssadr.strw.leidenuniv.nl/doku.php?id=spidx)
   [Download Spectral Index Catalog](http://tgssadr.strw.leidenuniv.nl/spidx/spidxcat_v1.1b.fits)

## Usage

1. Edit the configuration JSON file (e.g., `config.json`), The code will create a directory called {region} eg. ~/dir/{G033.0-5.0} where the files generated from the code are saved.  

3. Execute the script using the command:

   ```
   python sibpf.py --c config.json
   ```

## Functionality

The script provides the following functionalities:

### Crossmatching

- Matches the sources within a specified circular region in the provided FITS image (`image.fits`) obtained using AIPS or CASA with the TGSS and NVSS catalogs.
- Generates plots showing matched sources on an RA-Dec plot, with marker sizes proportional to resolution.

### Spectral Analysis

- Calculates spectral indices for the matched sources using TGSS and NVSS flux data.
- Produces histogram of the spectral index of the matched sources. 
- Identifies potential pulsar candidates based on a spectral index threshold and these candidates are saved to a separate file. 

### Image Inverter(Experimental) 

- The Image inverter code is included which allows the user to find artifacts that may have been detected by pybdsf as true sources
- The code creates an inverted image and then runs pybdsf on the image to locate any strong artifact usually around strong sources.
- This can be done before the crossmatching to check if any of the artifacts are confused as true sources and are crossmatched because of their proximity to the true source. 

### Variability (experimental)
- An experimental variability code is being worked on to check the variability of the pulsar candidates. Pulsars show a variability with time and this can be used to constrain the candidates further. 

### Output

The script generates the following outputs:

- Matched source information in CSV format (`{your_region}_spidx.csv`). This includes a catalog-like format generated using pybdsf along with the spectral index values from the spidx catalog from the TGSS survey, please check the example file for details. 
- Summary plots of spectral index distribution and SPIDX in the image field (`spectral_index_hist.png` and `spectral_index_hist_spidx.png`).
- The matched sources are plotted in (`matched_sources_{region}.png`).
- Pulsar candidates list in CSV format (`Pulsar_candidates_{your_region}.csv`).


# Acknowledgments
This code was developed as part of my Master's project. I thank Dr. Mayuresh Surnis for his invaluable insights and guidance during the development process. I also thank him for generously sharing data from the 24_051 GMRT cycle, which greatly contributed to this work. 

## TGSS Catalog

I acknowledge the use of data from the TGSS (TIFR GMRT Sky Survey) catalog in this project. The TGSS survey has provided invaluable radio astronomy data that has contributed to the analysis. For more information about the TGSS catalog and data access.

## PyBDSF

I acknowledge the developers of PyBDSF (Python Blob Detection and Source Finder), a powerful tool that has been instrumental in our source extraction and analysis. PyBDSF's functionality and capabilities have greatly facilitated our radio astronomy research. For more information about PyBDSF and its features, please refer to their GitHub repository: [PyBDSF GitHub Repository](https://github.com/lofar-astron/PyBDSF).

I acknowledge the use of various other libraries and tools that have contributed to the success of this project.

## Citation
If you use this code or the results obtained from it in your research or work, kindly consider citing relevant publications or resources associated with the TGSS catalog, PyBDSF, and NVSS catalog including any other relevant tool that has contributed to your analysis

---

# Note 
Make sure you have properly referenced the TGSS  and the Spidx.fits catalogs and provided accurate paths to the required files in your configuration. The files are present for download


