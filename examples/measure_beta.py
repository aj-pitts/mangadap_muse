import os
import glob
import argparse
import warnings

from mangadap.util.parser import DefaultConfig
from beta_corr import CubeData
import beta_corr_old

def get_args():
    """
    Command line input arguments.

    Returns:
            :class:`argparse.Namespace`: Converts argument
            strings to objects and assigns them as attributes to the class `Namespace`.
    """
    parser = argparse.ArgumentParser(description='File to run the MaNGA DAP on a MUSE cube.')

    parser.add_argument('galname', type=str,help='input galaxy name.')

    parser.add_argument('bin_key', type=str, help='input DAP spatial binning method.')

    parser.add_argument('-V', '--verbose', action='store_true', default=False, help='Keyword to print information about the bin modes and beta medians being fit in each iteration.')

    parser.add_argument('-o', '--old', action='store_true', default=False, help="Keyword to create beta measurements and plot results using `beta_corr_old.py`. Does NOT write beta tables.")

    parser.add_argument('--std', action='store_true', default=False, help="Keyword to include the standard devation of the beta distributions as the uncertainty on each beta median.")

    return parser.parse_args()

def main(args):
    # fit_one_cube_muse.py directory path
    file_dir = os.path.dirname(os.path.abspath(__file__))

    # main cube directory path
    #main_cube_dir = os.path.join(os.path.dirname(file_dir), 'MUSE_cubes')
    main_cube_dir = "/data2/muse/muse_cubes/"

    # cube directory path
    cube_dir = os.path.join(main_cube_dir, args.galname)
    if not os.path.isdir(cube_dir):
        raise ValueError(f'{cube_dir} is not a directory within /data2/')
    
    # check if there is only one config file in the cube directory
    # & input configuration file path
    ini_fils = glob.glob(f"{cube_dir}/*.ini")
    if len(ini_fils) == 0:
        raise ValueError(f"No configuration file found in {cube_dir}")
    
    if len(ini_fils) > 1:
        warnings.warn(f'Multiple .ini files within {cube_dir}')
        for fil in ini_fils:
            if args.galname in fil:
                config_fil = fil
            else:
                config_fil = ini_fils[0]
    
        print(f"Defaulting configuration to {config_fil}")
    
    else:
        config_fil = ini_fils[0]
        print(f"Found configuration file {config_fil}")

    # get parameter values from config file
    cfg = DefaultConfig(config_fil, interpolate=True)
    plate = cfg.getint('plate',default=None)
    ifu = cfg.getint('ifu',default=None)

    if args.old:
        cube = beta_corr_old.CubeData(galname=args.galname,bin_key=args.bin_key,plate=plate,ifu=ifu)
        cube.get_beta()
        cube.mk_hist_plots(show_plots=False, save_plots=True)
        cube.mk_violin_plots(show_plots=False, save_plots=True)
        return

    # initialize object and write beta correction files
    cube = CubeData(galname=args.galname,bin_key=args.bin_key,plate=plate,ifu=ifu)
    # get beta values
    cube.get_beta()
    # save histogram plots
    cube.mk_hist_plots(show_plots=False, save_plots=True)
    # write files to mangadap/data/beta_tables directory under the galaxy name
    cube.create_beta_tables(verbose=args.verbose, std=args.std)
    # save violin plots
    cube.mk_violin_plots(show_plots=False, save_plots=True)


if __name__ == "__main__":
    args = get_args()
    main(args)

