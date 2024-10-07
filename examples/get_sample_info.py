import numpy as np
from astropy.io import fits
from astropy.time import Time
import argparse
import os
from glob import glob
import logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)

def main():
    scriptdir = os.path.abspath(os.path.dirname(__file__))
    repodir = os.path.dirname(scriptdir)

    cubedir = os.path.join(repodir, "MUSE_cubes")
    cubesubdir = os.path.join(cubedir,"datacubes")

    galdirs = glob(os.path.join(cubedir, "NGC*"))
    galdirs2 = glob(os.path.join(cubesubdir, "NGC*"))
    sampledirs = np.concatenate((galdirs,galdirs2))
    
    sky_res = []
    spec_res = []
    exp_times = []
    dates = []
    ao = {}
    count=0
    logging.info("Acquiring sample info")
    for dir in sampledirs:
        fitsfil = glob(os.path.join(dir, "*.fits"))
        header = fits.getheader(fitsfil[0])

        sky_res.append(header['SKY_RES'])
        spec_res.append(header['SPEC_RES'])
        exp_times.append(header['TEXPTIME'])
        dates.append(header['DATE-OBS'][:10])

        ins_mode = header['HIERARCH ESO INS MODE']
        if 'NOAO' not in ins_mode:
            galname = header['OBJECT']
            ao[galname] = ins_mode

        count+=1
    logging.info(f"Done. Information for {count} galaxies acquired")

    print(f"SKY RES: {np.min(sky_res)} | {np.max(sky_res)} | {np.mean(sky_res)} (arcsec)")
    print(f"SPEC Res: {np.min(spec_res)} | {np.max(spec_res)} | {np.mean(spec_res)}")
    print(f"EXP TIME: {np.min(exp_times)} | {np.max(exp_times)} | {np.median(exp_times)}")
    print(f"OBSDATE: {Time(dates).min()} | {Time(dates).max()}")

    if len(ao)>0:
        for key in ao.keys():
            print(f"{key}: {ao[key]}")


if __name__ == "__main__":
    main()
