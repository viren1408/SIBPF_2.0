import argparse
import json
import os
import logging
from datetime import datetime
from crossxspec import sibpf  # Import the updated version of the crossmatching function

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def setup_output_dir(base_dir, region):
    """Create output directory structure based on region and timestamp."""
    #timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = os.path.join(base_dir, region)
    os.makedirs(output_dir, exist_ok=True)
    logging.info(f"Output directory created: {output_dir}")
    return output_dir

def save_to_file(data, filename):
    """Save data to a file."""
    with open(filename, 'w') as f:
        f.write(data)
    logging.info(f"Data saved to {filename}")

def main():
    # Create a parser
    parser = argparse.ArgumentParser(description="Crossmatching and spectral analysis tool")

    # Add an argument for the configuration file
    parser.add_argument('--config', required=True, help="Path to configuration JSON file")
    parser.add_argument('--verbose', action='store_true', help="Increase output verbosity")
    

    args = parser.parse_args()


    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)


    with open(args.config, 'r') as config_file:
        config = json.load(config_file)

    # Setup output directory
    output_dir = setup_output_dir(config['dir'], config['region'])

    try:
        # Call the crossmatching function
        spectral_index, pulsar_candidates = sibpf(
            dir=config['dir'],
            region=config['region'],
            spec_index_constr= config['spec_index_constr'],
            file_path_image=config['image'],
            file_path_spidx=config['spidx'],
            file_path_TGSS=config['tgss'],
            show_matches=config['show_matches'],
            get_spectral_index=config['get_spectral_index'],
            get_candidates=config['get_candidates'],
            get_pulsars=config['get_pulsars']
        )


        spectral_index_file = os.path.join(output_dir, f'{config["region"]}_spectral_index.txt')
        pulsar_candidates_file = os.path.join(output_dir, f'{config["region"]}_pulsar_candidates.txt')

        save_to_file(str(spectral_index), spectral_index_file)
        save_to_file(str(pulsar_candidates), pulsar_candidates_file)

        # Print results if needed (optional)
        print("Spectral Index and Pulsar Candidates have been saved to:")
        print(spectral_index_file)
        print(pulsar_candidates_file)

    except Exception as e:
        logging.error(f"An error occurred during processing: {e}")
        raise

if __name__ == "__main__":
    main()
