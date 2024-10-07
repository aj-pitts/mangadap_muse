#### UNDER CONSTRUCTION

import os
import argparse
from glob import glob
import warnings

def get_args():

    parser = argparse.ArgumentParser(description="A script to create an equivalent width map of ISM Na I for beta-corrected DAP outputs.")
    
    parser.add_argument('galname',type=str,help='Input galaxy name.')
    parser.add_argument('bin_method',type=str,help='Input DAP patial binning method.')
    parser.add_argument('--air', action='store_true',default=False)
    
    return parser.parse_args()

def main(args):
    dap_output_dir = "/data2/muse/dap_outputs/"
    data_dir = os.path.join(dap_output_dir, f"{args.galname}-{args.bin_method}")

    if args.air:
        plot_inspect_dir = "/Users/andrew/Repo/nai_analysis/plot_inspect"
    else:
        plot_inspect_dir = "/Users/apitts4030/Repo/NaI_Analysis/plot_inspect"
    
    local_galdir = os.path.join(plot_inspect_dir, f"{args.galname}-{args.bin_method}")


    qa_dirs = glob(os.path.join(data_dir, "**", "qa", "*.png"), recursive=True)
    qa_dir = None
    if len(qa_dirs) > 1:
        for dir in qa_dirs:
            if "NO-CORR" in dir:
                qa_dir = dir
                print(f"Found qa plot directory: {qa_dir}")
    elif len(qa_dirs) == 0:
        warnings.warn(f"No directory qa/ found within {data_dir}")
    else:
        qa_dir = qa_dirs[0]
        print(f"Found qa plot: {qa_dir}")



    beta_plot_dirs = glob(os.path.join(data_dir, "**", "beta_plots"), recursive=True)
    beta_plot_dir = None
    if len(beta_plot_dirs) > 1:
        warnings.warn(f"More than one beta_plots/ found within {data_dir}\nDefaulting to {beta_plot_dirs[0]}")
        beta_plot_dir = beta_plot_dirs[0]
    
    elif len(beta_plot_dirs) == 0:
        warnings.warn(f"No beta_plots/ found within {data_dir}")

    else:
        beta_plot_dir = beta_plot_dirs[0]
        print(f"Found beta plot directory: {beta_plot_dir}")


    ewmap_dirs = glob(os.path.join(data_dir, "**", "*EW_map*.png"), recursive=True)
    ewmap_dir = None
    if len(ewmap_dirs) == 0:
        warnings.warn(f"No EW map plots found in {data_dir}")
    else:
        for dir in ewmap_dirs:
            if os.path.exists(dir):
                ewmap_dir = dir
                print(f"Found equivalent width plot: {ewmap_dir}")


    pw = input(f"Enter inrainbows password for convenience or press enter to skip.")

    print("\n")
    print("SCP commands:")
    if len(pw)>0:
        print(pw)
    print("\n")
    if qa_dir:
        print(f"scp apitts@inrainbows:{qa_dir} {os.path.join(local_galdir, 'qa/')}")
        print("\n")
    if beta_plot_dir:
        subdirs = os.listdir(beta_plot_dir)
        print(f"scp -r apitts@inrainbows:{beta_plot_dir}/{subdirs[0]} {os.path.join(local_galdir, 'beta_plots/')}")
        print(f"scp -r apitts@inrainbows:{beta_plot_dir}/{subdirs[1]} {os.path.join(local_galdir, 'beta_plots/')}")
        print("\n")
    if ewmap_dir:
        print(f"scp apitts@inrainbows:{ewmap_dir} {os.path.join(local_galdir, 'ew_map/')}")

if __name__ == "__main__":
    args = get_args()
    main(args)
