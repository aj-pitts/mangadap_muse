#!/usr/bin/env python3

import os
import glob
import argparse

from IPython import embed

class DAP_MUSE:
    """
    Main wrapper class for running the DAP on a MUSE cube.

    """
    def __init__(self,galname=None,plate=None,ifudesign=None,bin_key=None,
                 beta_corr=False,beta_dir=None):

        self.galname = galname
        self.plate = plate
        self.ifudesign = ifudesign
        self.bin_key = bin_key
        self.beta_corr = beta_corr
        self.beta_dir = beta_dir

    def fit_one_cube_muse(self,config_file, directory_path=None,analysis_path=None,
                          analysis_plan=None):
        r"""
           Method to execute the MaNGA DAP on a single MUSE cube.

           This function is designed to be called once per datacube. The
           :class:`mangadap.par.analysisplan.AnalysisPlanSet` instance sets
           the types of analyses to perform on this observation. Each
           analysis plan results in a MAPS and model LOGCUBE file.

           The procedure is as follows.  For each plan in ``analysis_plan``:

               - Determine basic assessments of the MUSE_cubes, including the S/N
                 and spaxel coordinates. See
                 :class:`mangadap.proc.reductionassessments.ReductionAssessment`.
               - Bin the spectra. See
                 :class:`mangadap.proc.spatiallybinnedspectra.SpatiallyBinnedSpectra`.
               - Fit the stellar continuum for stellar kinematics. See
                 :class:`mangadap.proc.stellarcontinuummodel.StellarContinuumModel`.
               - Measure the emission-line moments. See
                 :class:`mangadap.proc.emissionlinemoments.EmissionLineMoments`.
               - Fit parameterized line profiles to the emission lines. See
                 :class:`mangadap.proc.emissionlinemodel.EmissionLineModel`.
               - Subtract the fitted emission-line models and measure the
                 spectral indices. See
                 :class:`mangadap.proc.spectralindices.SpectralIndices`.
               - Construct the primary output files based on the plan
                 results. See :func:`mangadap.dapfits.construct_maps_file`
                 and :func:`mangadap.dapfits.construct_cube_file`.

           Verbose levels (still under development):
               0. Nothing but errors.
               1. Basic status updates. E.g., start and end of each block,
                  minor progress within blocks.
               2. Block-level and sub-function updates.
               3. Same as above, with figures.

           Args:
                config_file (str):
                    Configuration file containing the initial galaxy inputs.
               directory_path (str):
                    Path to directory containing the MUSE cube file.
               analysis_path (str):
                   Top-level directory for the DAP output MUSE_cubes; default is
                   defined by :func:`mangadap.config.defaults.dap_analysis_path`.

           Returns:
               manga_dap (:func:`mangadap.survey.manga_dap.manga_dap`):
                    Wrapper function that runs the DAP entirely.
        """

        # mask the first pixel in the cube
        cube = MUSEDataCube.from_config(config_file)
        cube.mask[0, 0, :] = True
        # input correlation correction parameters
        cube.beta_corr = self.beta_corr
        cube.beta_dir = self.beta_dir

        # Run it!
        return manga_dap(cube, analysis_plan, verbose=2, directory_path=directory_path,
                         analysis_path=analysis_path)

    def run_MUSE_cube(self,config_fil=None,directory_path=None,analysis_plan=None):

        # fit_one_cube_muse_test.py directory path
        file_dir = os.path.dirname(os.path.abspath(__file__))

        # main output directory
        output_root_dir = os.path.join(os.path.dirname(file_dir), 'outputs')
        # make the directory if there isn't one already
        if not os.path.isdir(output_root_dir):
            os.makedirs(output_root_dir)

        # galaxy output directory name
        # beta subdirectory in mangadap.MUSE_cubes.beta_tables should have the same name as gal_name
        gal_name = self.galname
        bin_method = analysis_plan['bin_key'][0]
        output_gal_dir = os.path.join(output_root_dir, gal_name + '-' + bin_method)
        if not os.path.isdir(output_gal_dir):
            os.makedirs(output_gal_dir)

        # output galaxy subdirectory
        if not self.beta_corr:
            output_gal_sub_dir = os.path.join(output_gal_dir, 'NO-CORR')
        else:
            output_gal_sub_dir = os.path.join(output_gal_dir, 'BETA-CORR')

        # fit a MUSE cube!
        self.fit_one_cube_muse(config_fil, directory_path=directory_path, analysis_path=output_gal_sub_dir,
                               analysis_plan=analysis_plan)
        # create output plot from ppxf fitting results
        ppxffit_qa_plot(self.plate, self.ifudesign, analysis_plan ,drpver=None, redux_path=directory_path,
                        dapver=None, analysis_path=output_gal_sub_dir, tpl_flux_renorm=None)

# -----------------------------------------------------------------------------
def get_args():
    parser = argparse.ArgumentParser(description='Wrapper File to run the MaNGA DAP on a MUSE cube.')

    parser.add_argument('-bc', '--beta_corr', action='store_true', default=False,
                        help='Flag to specify the DAP to perform a correlation correction on the MUSE cube.')

    parser.add_argument('-no-bc', '--no-beta_corr', dest='beta_corr',action='store_false',
                        help='Flag to specify the DAP to not perform a correlation correction on the MUSE cube.')

    parser.add_argument('--beta_dir', type=str, default=None,
                        help='Beta table data directory name for the correlation correction process.')
    return parser.parse_args()

# -----------------------------------------------------------------------------
def main(args):

    # fit_one_cube_muse_test.py directory path
    file_dir = os.path.dirname(os.path.abspath(__file__))

    # test cube directory path
    cube_dir = os.path.join(os.path.dirname(file_dir), 'MUSE_cubes','test_cube')
    # input configuration file path
    config_fil = glob.glob(f"{cube_dir}/*.ini")[0]

    # get plate and ifu parameter values from config file
    cfg = DefaultConfig(config_fil, interpolate=True)
    plate = cfg.getint('plate',default=None)
    ifu = cfg.getint('ifu',default=None)
    galname = cfg.get('galname',default=None)

    if plate is None or ifu is None:
        raise ValueError('Configuration file must define the plate and IFU.')
    if galname is None:
        raise ValueError('Configuration file muse define the galaxy name.')

    # Define how you want to analyze the MUSE_cubes
    plan = AnalysisPlanSet([AnalysisPlan(
                                        # Data reduction quality assessment keyword.
                                        # - Defines how to assess the reduced MUSE_cubes
                                        # - 'SNRG' keyword assesses quality of each spaxel by its S/N ratio.
                                        drpqa_key='SNRG',
                                        # Overwrite existing MUSE_cubes-quality assessment reference files if True.
                                        # - if False, the DAP searches for these references files to skip
                                        #   the process if files exist already from a previous execution.
                                        drpqa_clobber=True,
                                        # Spatial binning method keyword.
                                        # - Places spaxels in square bins with a certain
                                        #   bin size length (in units of arcsec).
                                        # - Available SQUARE method bin size options are
                                        #   SQUARE2.0, SQUARE1.0 and SQUARE0.6
                                        bin_key='SQUARE2.0',
                                        # Overwrite any existing spatial binning reference files.
                                        # Same as dqpa_clobber
                                        bin_clobber=True,
                                        # Stellar-continuum fitting method keyword.
                                        # - Defines which stellar continuum model to use.
                                        # -'MILESHC-NOISM' uses a new version of the MILES Hierarchical
                                        #   Cluster Library containing a significant
                                        #   reduction in the model overestimation of stellar absorption.
                                        continuum_key='MILESHC-NOISM',
                                        # Overwrite any existing stellar-continuum fitting reference files
                                        continuum_clobber=True,
                                        # Emission-line moments measurement method keyword
                                        # -Measures the moments of the observed emission lines.
                                        # 'EMOMMPL11' uses the ELBMPL9 emission-line bandpass database.
                                        #  These Line wavelengths are "Ritz" wavelengths from NIST:
                                        #  http://physics.nist.gov/PhysRefData/ASD/Html/help.html
                                        elmom_key='EMOMMPL11',
                                        # Overwrite any existing emission-line moments reference files
                                        elmom_clobber=True,
                                        # Emission-line modeling method keyword
                                        # - Used to define which emission-line moment database to use
                                        # - This process conducts the emission-line profile fitting.
                                        # - 'EFITMPL11-ISMMASK-HCNOISM' places an ISM mask on the ELPMPL11
                                        #    library and uses the MaStar Hierarchically Clustered Library v2
                                        #    for the continuum templates.
                                        elfit_key='EFITMPL11-ISMMASK-HCNOISM',
                                        # Overwrite any existing emission-line modeling reference files
                                        elfit_clobber=True,
                                        # Spectral-index measurement method keyword
                                        # - Defines which spectral indices library to use.
                                        # - 'INDEX' uses the bandhead indices from:
                                        #     D4000   - Bruzual (1983, ApJ, 273, 105)
                                        #     Dn4000  - Balogh et al. (1999, ApJ, 527, 54)
                                        #     TiOCvD  - Conroy & van Dokkum (2012, ApJ, 747, 69)
                                        spindex_key='INDXEN',
                                        # Overwrite any existing spectral-index reference files
                                        spindex_clobber=True)])

    # argparse check
    if not args.beta_corr:
        if args.beta_dir is not None:
            file_name = os.path.basename(os.path.abspath(__file__))
            raise ValueError('Must enter flag -bc or --beta_corr on the command line '
                             'if providing a beta table directory! \n'
                            f'{" "*12}Type python {file_name} -h on the command line for '
                             'more information.')
    else:
        if args.beta_dir is None:
            file_name = os.path.basename(os.path.abspath(__file__))
            raise ValueError('Must provide a beta table directory name to the command line! '
                             f'{" "*12}Type python {file_name} -h on the command line for '
                             'more information.')

    # instantiate object with input parameters
    muse_obj = DAP_MUSE(galname=galname,plate=plate,
                        ifudesign=ifu, beta_corr=args.beta_corr,beta_dir=args.beta_dir)

    # fit MUSE cube!
    muse_obj.run_MUSE_cube(config_fil=config_fil,directory_path=cube_dir,analysis_plan=plan)

#-----------------------------------------------------------------------------
if __name__ == '__main__':
    # get arguments from command line
    args = get_args()

    # import the MaNGA DAP modules
    from mangadap.datacube import MUSEDataCube
    from mangadap.survey.manga_dap import manga_dap
    from mangadap.par.analysisplan import AnalysisPlan, AnalysisPlanSet
    from mangadap.scripts.ppxffit_qa import ppxffit_qa_plot
    from mangadap.util.parser import DefaultConfig

    # pass args to main class
    main(args)