import os
import re
from tqdm import tqdm
import time
import argparse
import glob

def remove_older_files(directory, galstring, dry=False):
    # Regex pattern to match the filenames
    print(directory)

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

    # Step 1: Identify all files with "run-0" in the filename
    run_0_files = glob.glob(os.path.join(directory, '*run-0*.fits'))

    if not run_0_files:
        print("No files with 'run-0' found. Cannot determine the most recent run.")
        return

    # Step 2: Get the most recently modified file with "run-0" in its name
    most_recent_file = max(run_0_files, key=os.path.getmtime)
    most_recent_time = os.path.getmtime(most_recent_file)

    # Step 3: Set a time window (e.g., +/- 48 hours) for identifying recent files
    time_threshold = 48 * 3600  # 48 hours in seconds
    lower_bound = most_recent_time - time_threshold
    upper_bound = most_recent_time + time_threshold

    print(f"Most recent run identified by {most_recent_file}. Time window: "
          f"{time.ctime(lower_bound)} to {time.ctime(upper_bound)}")

    # Step 4: Iterate through all the files in the directory and delete outdated ones
    files_to_keep = []

    for file in run_0_files:
        file_mod_time = os.path.getmtime(file)
        if lower_bound <= file_mod_time <= upper_bound:
            files_to_keep.append(file)

    # Track the most recent run (based on the run numbers)
    most_recent_run = max(run_files.keys(), key=lambda run: max(os.path.getctime(f) for f in run_files[run]))

    # Go through each run group and remove duplicates
    for run_number, files_in_run in run_files.items():
        if len(files_in_run) > 1:
            # Sort files by modification time (oldest first)
            files_in_run.sort(key=lambda f: os.path.getmtime(f))

            # If it's not the most recent run, delete all but the most recently modified file
            if run_number != most_recent_run:
                # Print the files that would be deleted (dry run)
                for file_to_delete in files_in_run[:-1]:  # Exclude the most recently modified
                    if dry:
                        print(f"Would delete: {file_to_delete}")
                    else:
                        os.remove(file_to_delete)  # Actual delete
            else:
                # Add the most recent file of the most recent run to the list of files to keep
                most_recent_file = files_in_run[-1]
                files_to_keep.append(most_recent_file)

        else:  # If there's only one file in the run group, keep it
            files_to_keep.append(files_in_run[0])

    # If dry_run is enabled, also print the list of files to be kept
    if dry:
        print("\nFiles that would be kept:")
        for file_to_keep in files_to_keep:
            print(file_to_keep)


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
            datapath = os.path.join(datadir, dir, "BETA-CORR")
            remove_older_files(datapath, os.path.basename(dir), dry=dry)
    else:
        dir = os.path.join(datadir, input_dir, "BETA-CORR")
        remove_older_files(dir, input_dir, dry=dry)

if __name__ == "__main__":
    main()