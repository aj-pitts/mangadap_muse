import os
import re
from tqdm import tqdm
import time
import argparse

def remove_older_files(directory, dry = False):
    # Regex pattern to match the filenames
    file_pattern = re.compile(fr"{directory}-binid-(\d+)-(\d+)-samples-run-(\d+)\.fits")

    # Store files grouped by "run"
    files_by_run = {}

    # Iterate through files in the directory
    for filename in os.listdir(directory):
        print(filename)
        match = file_pattern.match(filename)
        if match:
            print("Match True")
            # Extract run number from filename
            run = int(match.group(3))

            # Get the full file path
            file_path = os.path.join(directory, filename)

            # Get the file's creation time (or last modification time)
            creation_time = os.path.getctime(file_path)

            # Group files by "run" and store (filename, creation time)
            if run not in files_by_run:
                files_by_run[run] = []
            files_by_run[run].append((filename, creation_time))

    if dry:
        files_kept = []
    # Process each group
    for run, files in files_by_run.items():
        if len(files) == 1:  # Skip if only one file
            print(f"Only one file found for run {run}: {files[0][0]} (no deletion)")
        else:
            # Sort files by creation time (oldest first)
            files.sort(key=lambda x: x[1])

            # Keep the most recent file (last in the sorted list)
            most_recent_file = files[-1][0]
            if dry:
                files_kept.append(most_recent_file)
            # Delete all files except the most recent
            for file, _ in files[:-1]:
                if dry: 
                    print(f"Would remove {file_path_to_delete}")
                else:
                    file_path_to_delete = os.path.join(directory, file)
                    os.remove(file_path_to_delete)
                    print(f"Deleted: {file}")
    if dry:
        print(f"Files Kept: {files_kept}")

# Specify the directory containing the files

def main():
    parser = argparse.ArgumentParser(description="Remove duplicate MCMC files based on creation date.")
    parser.add_argument("gal_directory", type=str, nargs="?", default=None, help="Galaxy directory defined as <GALNAME>-<BIN_METHOD> containing the files to check (optional)")
    parser.add_argument("-dr","--dry-run", action="store_true", help="If set, will only print files to be deleted, without actually deleting them")
    args = parser.parse_args()


    dry = args.dry_run
    input_dir = args.gal_directory
    

    datadir = '/data2/muse/mcmc_outputs'

    if input_dir is None:
        subdir_list = os.listdir(datadir)
        for dir in tqdm(subdir_list, desc="Cleaning MCMC outputs."):
            datapath = os.path.join(dir, "BETA-CORR")
            remove_older_files(datapath, dry=dry)
    else:
        dir = os.path.join(datadir, input_dir, "BETA-CORR")
        remove_older_files(dir, dry=dry)

if __name__ == "__main__":
    main()