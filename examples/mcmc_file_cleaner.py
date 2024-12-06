import os
import re
from tqdm import tqdm
import time
import argparse

def remove_older_files(directory, galstring, dry=False, time_window=48*3600):
    # Regex pattern to match the filenames

    file_pattern = re.compile(fr"{galstring}-binid-(\d+)-(\d+)-samples-run-(\d+)\.fits")

    files = os.listdir(directory)

    # Dictionary to store the files grouped by run number
    run_files = {}

    # Group files by run number based on the filename pattern
    for file in files:
        match = file_pattern.match(file)
        if match:
            run_number = match.group(3)  # Extract the run number
            if run_number not in run_files:
                run_files[run_number] = []
            file_path = os.path.join(directory, file)
            run_files[run_number].append(file_path)

    # List to store the files that would be kept (most recent file for each run group)
    files_to_keep = []

    # Track the most recent run (based on the modification time of files in the group)
    for run_number, files_in_run in run_files.items():
        # Sort files by modification time (most recent first)
        files_in_run.sort(key=lambda f: os.path.getmtime(f), reverse=True)

        # Most recent file based on modification time within the run
        most_recent_file = files_in_run[0]
        most_recent_time = os.path.getmtime(most_recent_file)

        # Identify the time window (within most recent file's time +/- window)
        time_threshold_start = most_recent_time - time_window
        time_threshold_end = most_recent_time + time_window

        # Check if each file in the run falls within the time window
        for file_to_check in files_in_run:
            file_time = os.path.getmtime(file_to_check)
            if time_threshold_start <= file_time <= time_threshold_end:
                # If it's within the time window, keep it
                files_to_keep.append(file_to_check)
            else:
                # If it's outside the time window, delete it
                if dry:
                    mod_time = os.path.getmtime(file_to_check)
                    formatted_time = time.ctime(mod_time)  # Format the modification time
                    print(f"Would delete: {file_to_check} (Last modified: {formatted_time})")
                else:
                    os.remove(file_to_check)  # Actual delete

    # If dry_run is enabled, also print the list of files to be kept
    if dry:
        print("\nFiles that would be kept:")
        # Sort files to keep by run number (extracted from filenames)
        files_to_keep.sort(key=lambda f: re.search(r"run-(\d+)", f).group(1))  # Sort by run number
        for file_to_keep in files_to_keep:
            mod_time = os.path.getmtime(file_to_keep)
            formatted_time = time.ctime(mod_time)  # Format the modification time
            print(f"{file_to_keep} (Last modified: {formatted_time})")


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
        for dir in tqdm(subdir_list, desc=f"Cleaning MCMC outputs: {dir}."):
            datapath = os.path.join(datadir, dir, "BETA-CORR")
            remove_older_files(datapath, os.path.basename(dir), dry=dry)
    else:
        dir = os.path.join(datadir, input_dir, "BETA-CORR")
        remove_older_files(dir, input_dir, dry=dry)

if __name__ == "__main__":
    main()