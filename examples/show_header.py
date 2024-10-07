import numpy as np
from astropy.io import fits
import argparse
import os
from glob import glob
import logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)


def get_args():
    parser = argparse.ArgumentParser(description="Display the header information from the Reduced MUSE datacube.")
    parser.add_argument("galname", type=str, help="Input galaxy name.")

    return parser.parse_args()

def main(args):
    scriptdir = os.path.abspath(os.path.dirname(__file__))
    repodir = os.path.dirname(scriptdir)

    cubedir = os.path.join(repodir, "MUSE_cubes")
    galdir = os.path.join(cubedir,args.galname)

    fil = os.path.join(galdir, f"{args.galname}.fits")
    if not os.path.exists(fil):
        logging.warning(f"{args.galname}.fits not found in {galdir}\nSearching for .fits file")

        fils = glob(os.path.join(galdir,"*.fits"))
        if len(fils)<1:
            raise NameError(f"No .fits files found for {args.galaname} in {galdir}")
        if len(fils)>1:
            logging.warning(f"Multiple .fits files found in {galdir}")
        
        fil = fils[0]
    
    header = fits.getheader(fil,ext=0)

    print(header)

if __name__ == "__main__":
    args = get_args()
    main(args)