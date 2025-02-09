import os
import warnings
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier
import astropy.units as u
from scipy.stats import linregress
from psrqpy import QueryATNF

warnings.simplefilter(action='ignore')




def read_fits_image(file_path):
    """Open and read the header from a FITS file."""
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File {file_path} not found.")
    with fits.open(file_path) as hdul:
        header = hdul[0].header
    return header

def get_pixel_scale(header):
    """Calculate the pixel scale from the FITS header."""
    cdelt1 = header.get('CDELT1',0)
    cdelt2 = header.get('CDELT2',0)
    return np.sqrt(cdelt1**2 +cdelt2**2)

def get_image_radius(header,pixel_scale):
    """Calculate the radius of the image in degrees based on pixel size."""
    nx = header.get('NAXIS1',0)
    ny = header.get('NAXIS2',0)
    return (np.sqrt((nx/2)**2 +(ny/2)**2) *pixel_scale)/2

def create_source_catalog(file_path_image,dir,region):
    """Run PyBDSF to create source catalog from a FITS image."""
    import bdsf  # Assumed external dependency
    #save_file = os.path.join(save_dir, f"{region}.sav") #this was actually not required ! 
    img = bdsf.process_image(file_path_image, quiet=True)
    img.write_catalog(outfile=os.path.join(dir,region, f"{region}.pybdsm.gaul"), format='csv',catalog_type='gaul',clobber=True)
    img.export_image(outfile=os.path.join(dir,region,f"{region}.pybdsm_gaus_model.fits"), img_format='fits', img_type='gaus_model',clobber=True)

def load_observed_sources(file_path_catalog):
    """Load the observed sources catalog."""
    data = pd.read_csv(file_path_catalog, index_col=None, skiprows=5)
    observed_sources = pd.DataFrame(data)
    observed_sources.set_index("# Gaus_id", inplace=True)
    return observed_sources

def match_sources(spix,observed_sources,theta_obs):
    """Match SPIDX sources with observed sources."""
    for i, source in spix.iterrows():
        ra_1, dec_1 = source['RA'],source['DEC']
        min_theta = np.inf
        matched_source_index = None

        for j, obs_source in observed_sources.iterrows():
            ra_2, dec_2 = obs_source[' RA'],obs_source[' DEC']
            theta = calculate_angular_separation(ra_1, dec_1,ra_2,dec_2)

            if abs(theta) <=theta_obs and theta <min_theta:
                min_theta = theta
                matched_source_index = j

        if matched_source_index is not None:
            update_matched_sources(observed_sources, matched_source_index,source,min_theta)
    
    # Remove unmatched sources
    observed_sources.dropna(subset=['Spidx', 'E_Spidx'], inplace=True)
    observed_sources.reset_index(drop=True, inplace=True)
    return observed_sources

def calculate_angular_separation(ra_1, dec_1,ra_2, dec_2):
    """Calculate the angular separation between two RA, Dec coordinates."""
    return np.degrees(np.arccos(np.sin(np.radians(dec_1)) *np.sin(np.radians(dec_2)) +np.cos(np.radians(dec_1))* np.cos(np.radians(dec_2)) *np.cos(np.radians(ra_1 - ra_2))))

def update_matched_sources(observed_sources, matched_index, source, min_theta):
    """Update the matched sources with spectral index and other info."""
    observed_sources.at[matched_index,'RA_SPIDX'] = source['RA']
    observed_sources.at[matched_index, 'DEC_SPIDX'] = source['DEC']
    observed_sources.at[matched_index,'sep_arcsecods'] = min_theta * 3600
    observed_sources.at[matched_index,'Total_flux_NVSS'] = source['Total_flux_NVSS']
    observed_sources.at[matched_index,'E_Total_flux_NVSS'] = source['E_Total_flux_NVSS']
    observed_sources.at[matched_index,'Total_flux_TGSS'] = source['Total_flux_TGSS']
    observed_sources.at[matched_index,'E_Total_flux_TGSS'] = source['E_Total_flux_TGSS']
    observed_sources.at[matched_index,'Spidx'] =source['Spidx']
    observed_sources.at[matched_index,'E_Spidx'] =source['E_Spidx']

def plot_matched_sources(region, dir_path,ref_ra, ref_dec, matched_sources, observed_sources,NVSS, TGSS, theta_obs):
    """Plot the matched sources and save the plot as an image."""
    resolution_nvss = 45
    resolution_tgss = 25
    resolution_observed = theta_obs *3600

    #scaling it to the original resolution 
    marker_size_nvss = 5 *resolution_nvss 
    marker_size_tgss = 5 *resolution_tgss
    marker_size_gmrt = 5* resolution_observed

    fig, ax = plt.subplots(figsize=(15, 8))
    ax.scatter(NVSS['RAJ2000'], NVSS['DEJ2000'], color="black", label='NVSS', s=marker_size_nvss, alpha=0.2)
    ax.scatter(TGSS['RA'], TGSS['DEC'], color="black", label='TGSS', s=marker_size_tgss, alpha=0.3)
    ax.scatter(observed_sources[' RA'], observed_sources[' DEC'], color="black", label='Observed_sources', s=marker_size_gmrt, alpha=0.4)
    ax.scatter(matched_sources[' RA'], matched_sources[' DEC'], color="none", edgecolor="red", label='Matched_sources', s=marker_size_gmrt)

    plt.gca().invert_xaxis()
    ax.set_xlabel('RA (deg)')
    ax.set_ylabel('Dec (deg)')
    ax.set_title(f'Matched sources_{region}')
    ax.legend()
    ax.grid(True)

    file_loc = os.path.join(dir_path, region, 'matched_sources_{region}.png')
    plt.savefig(file_loc)
    #plt.show()

def primary_beam_correction(sep,oflux):
    a, b, c, d = -3.397, 47.192, -30.931, 7.803  #GMRT at 320 MHz 
    x = sep * 0.325

    corrfac = 1 + x**2 *(a / 1e3) + x**4 *(b / 1e7) + x**6*(c / 1e10) + x**8*(d / 1e13)
    mflux = oflux*corrfac
    return mflux

def compute_spectral_index(matched_sources,FREQUENCIES):
    """Compute the spectral index for each matched source."""
    spectral_index = matched_sources.copy()
    for i, source in matched_sources.iterrows():
        fluxes = np.array([source['Total_flux_TGSS'], source[' Total_flux'], source['Total_flux_NVSS']])
        log_of_flux = np.log10(fluxes)
        slope, intercept, r_value, p_value, std_err = linregress(np.log10(FREQUENCIES), log_of_flux) #this needs an update ! (WIP)
        spectral_index.at[i, 'Sr.no'] = i
        spectral_index.at[i, 'Spectral_index'] = slope
        spectral_index.at[i, 'E_Spectral_index'] = std_err/abs(slope)
    return spectral_index


def apply_primary_beam_correction(data, ref_ra, ref_dec):
    """
    Apply primary beam correction to the Total flux values in the observed sources dataframes.

    """
    data['Pbcorr_Sep(arcmin)'] = data.apply(lambda row: calculate_angular_separation(ref_ra, ref_dec, row[' RA'], row[' DEC']) * 60,axis=1)

    data[' Total_Flux'] = data.apply(lambda row: primary_beam_correction(row['Pbcorr_Sep(arcmin)'], row[' Total_flux']),axis=1)

    return data


def sibpf(dir, region, file_path_image, file_path_spidx, spec_index_constr, file_path_TGSS, show_matches=True, get_spectral_index=True, get_candidates=True, get_pulsars=True,do_pbcorr=True):
    """
    Main function to process FITS image, match sources, and compute spectral index.
    Parameters:
        dir (str): Directory to save files.
        region (str): Region name for file naming.
        file_path_image (str): Path to the FITS image file.
        file_path_spidx (str): Path to the spectral index file.
        file_path_TGSS (str): Path to the TGSS data file.
        show_matches (bool): Whether to display matched sources.
        get_spectral_index (bool): Whether to calculate spectral index.
    """
    header = read_fits_image(file_path_image)
    freq_obs = header.get('CRVAL3',0)

    FREQUENCIES = np.array([1.50e08, freq_obs , 1.4e09])  # in Hz    
    pixel_scale = get_pixel_scale(header)
    radius_deg = get_image_radius(header, pixel_scale)
    ref_ra = header.get('CRVAL1', 0)
    ref_dec = header.get('CRVAL2', 0)

    create_source_catalog(file_path_image, dir, region)
    observed_sources = load_observed_sources(os.path.join(dir, region, f"{region}.pybdsm.gaul"))

    if do_pbcorr == True :

        print('Primary beam correction is being applied /n Caution! : The correction is only valid for GSB band 3 at GMRT')

        observed_sources = apply_primary_beam_correction(observed_sources, ref_ra, ref_dec)

    dat = Table.read(file_path_spidx, format='fits')
    spix_main = dat.to_pandas()
    spix = spix_main[(spix_main['S_code'] == b'S') | (spix_main['S_code'] == b'M')]
    distances = np.sqrt((spix['RA'] - ref_ra)**2 + (spix['DEC'] - ref_dec)**2)
    spix['dist'] = distances
    spix = spix[spix['dist'] <= radius_deg]

    theta_obs = 0.5*header.get('BMAJ', 0)  # assuming theta_obs comes from the FITS header
    matched_sources = match_sources(spix, observed_sources, theta_obs)

    if show_matches:
        main_TGSS = pd.read_csv(file_path_TGSS,sep='\t')
        distances = np.sqrt((main_TGSS['RA'] - ref_ra)**2 + (main_TGSS['DEC'] - (ref_dec))**2)
        mask = distances <= radius_deg
        tgss = main_TGSS[mask]
        tgss = tgss.reset_index(drop=True)

        Vizier.ROW_LIMIT = -1
        print(f'The image centre is : {ref_ra},{ref_dec}, The radius is {radius_deg}')
        result = Vizier.query_region(coord.SkyCoord(ra=ref_ra, dec=ref_dec,unit=(u.deg, u.deg),frame='icrs'),radius= radius_deg*u.degree,catalog=["NVSS"])
        table = result[0]
        NVSS = table.to_pandas()

        coords = SkyCoord(NVSS['RAJ2000'], NVSS['DEJ2000'], unit=(u.hourangle, u.deg))
        NVSS['RAJ2000'] = coords.ra.degree
        NVSS['DEJ2000'] = coords.dec.degree

        plot_matched_sources(region, dir, ref_ra, ref_dec, matched_sources, observed_sources, NVSS, tgss, theta_obs)

    if get_spectral_index:
        spectral_index = compute_spectral_index(matched_sources,FREQUENCIES)
        file_name_spidx = os.path.join(dir, region, f'{region}_spidx.csv')
        spectral_index.to_csv(file_name_spidx)

    

    if get_candidates == True:
        print(f"'The spectral index constraint for selecting candidates is set to {spec_index_constr}: {spec_index_constr})")

    	#print(f: {spec_index_constr}')
        # TO DO :
        # How to deal with large error bars spectral index values?
        # Give a warning if the error bar is unusually high for a candidate.
        # If other error bars are low, what would the high error bar indicate ? > Should we try a different fit if not linear? 
        # Warn the user if the majority of the spectral index values have an unusually large error bar, which currently can be diagnosed using the candidate image 

        pulsar_candidates = spectral_index[spectral_index['Spectral_index'] < spec_index_constr]
        
        filename = str(dir)+'/'+str(region)+'/Pulsar_candidates_'+str(region)+'.csv'
        pulsar_candidates.to_csv(filename, index=False)
        fig, ax = plt.subplots()
        #plotting the candidates
        plt.scatter(spectral_index['Sr.no'],spectral_index['Spectral_index'],marker = '.')
        plt.scatter(pulsar_candidates['Sr.no'],pulsar_candidates['Spectral_index'],s=50, linewidth=0.5, alpha=0.7,color='r',label = 'Selected_candidates')
        plt.errorbar(spectral_index['Sr.no'],spectral_index['Spectral_index'],yerr =spectral_index['E_Spectral_index'],alpha=0.7,ecolor='black',capsize=2,ls = 'none')
        plt.axhline(y=spec_index_constr, color='b', linestyle='--', label=f'Spectral index = {spec_index_constr}')
        plt.axhspan(ymin=spectral_index['Spectral_index'].min(), ymax=-0.9, alpha=0.3, color='g', label=f'Spectral index < {spec_index_constr}')
        plt.legend()
        plt.grid()

        plt.title('Pulsar_candidates')
        plt.ylabel('Spectral Index')
        plt.xlabel('Sr.no')
        loc = str(dir)+'/'+str(region)+'/Pulsar_candidates.png'
        plt.savefig(loc)
        plt.show()


        ax_1 = spectral_index.hist(column='Spectral_index', bins='auto', grid=True, figsize=(12,8), color='#86bf91')

        plt.axvline(x={spec_index_constr}, color='b', linestyle='--', label=f'Spectral index = {spec_index_constr}')


        plt.axvspan(xmin=spectral_index['Spectral_index'].min(), xmax=spec_index_constr, alpha=0.3, color='r', label=f'Spectral index < {spec_index_constr}')


        plt.title('Spectral Index Distribution')
        plt.xlabel('Spectral Index')
        plt.ylabel('Counts')


        plt.legend()

        loc = str(dir)+'/'+str(region)+'/spectral_index_hist.png'
        plt.savefig(loc)
        plt.show


        ax_2 = spectral_index.hist(column='Spidx', bins='auto', grid=True, figsize=(12,8), color='#86bf91') #smooth green 

        plt.axvline(x={spec_index_constr}, color='b', linestyle='--', label=f'Spectral index = {spec_index_constr}') #updated for the new pulsar list (This needs to be updated to get the user to select the contraint !) 


        plt.axvspan(xmin=spectral_index['Spidx'].min(), xmax={spec_index_constr}, alpha=0.3, color='r', label=f'Spectral index < {spec_index_constr}')


        plt.title('Spectral Index Distribution of SPIDX')
        plt.xlabel('Spectral Index')
        plt.ylabel('Counts')


        plt.legend()

        loc = str(dir)+'/'+str(region)+'/spectral_index_hist_spidx.png'
        plt.savefig(loc)
        plt.show


    return spectral_index, pulsar_candidates


