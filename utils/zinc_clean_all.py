import os
import glob
import sys
import logging
from tqdm import tqdm
import zinc_clean
from zinc_clean import prepare_zinc

# Set up logging
log_file = "zinc_clean.log"  # Set your log file name here
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s', handlers=[
    logging.FileHandler(log_file),  # Log to a file
    logging.StreamHandler(sys.stdout)  # Also print to console
])

# Redirect stdout to logging
class LoggerWriter:
    def __init__(self, logger):
        
        self.logger = logger

    def write(self, message):
        if message.strip():  # To avoid logging empty messages
            self.logger.info(message.strip())

    def flush(self):
        pass  # No need to do anything here

if __name__ == "__main__":
    sys.stdout = LoggerWriter(logging.getLogger())
    pattern = os.path.join("./zinc_drug_like/*/*", '*.smi')
    smi_files = glob.glob(pattern, recursive=True)
    print(f"Number of files: {len(smi_files)}")

    # Filter out empty files and files with "clean" in the name
    filtered_smi_files = [f for f in smi_files if os.path.getsize(f) > 0 and "clean" not in os.path.basename(f)] 
    print(f"Number of non-empty files: {len(filtered_smi_files)}")
    del smi_files

    for f in tqdm(filtered_smi_files, desc="Processing files", unit="file"):
        bn = os.path.basename(f)
        dn = os.path.dirname(f)

        out_name = bn.replace(".smi", "_clean.smi")
        if not os.path.exists(dn + "/" + out_name): # check if the file with the "clean" name already exists
            #print(f)
            df = prepare_zinc(f) # clean the file
            df.to_csv(dn + "/" + out_name, sep=" ", index=False) # save the cleaned file by adding "_clean" to the file name
            print("Created: " + dn + "/" + out_name)
        else:
            print(f'{dn + "/" + out_name} already exists')
    logging.info(f"Logging informationsaved to {log_file}")
    sys.exit(0)
