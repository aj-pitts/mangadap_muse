#!/usr/bin/env python3

import os
from mangadap.datacube import MUSEDataCube
from mangadap.survey.manga_dap import manga_dap
from mangadap.par.analysisplan import AnalysisPlan, AnalysisPlanSet
from mangadap.scripts.ppxffit_qa import ppxffit_qa_plot
from IPython import embed

#-----------------------------------------------------------------------------
def fit_one_cube_muse(config_file, analysis_plan, directory_path=None, analysis_path=None):
    r"""
       Custom wrapper function to execute the MaNGA DAP on MUSE cubes.

       This function is designed to be called once per datacube. The
       :class:`mangadap.par.analysisplan.AnalysisPlanSet` instance sets
       the types of analyses to perform on this observation. Each
       analysis plan results in a MAPS and model LOGCUBE file.

       The procedure is as follows.  For each plan in ``analysis_plan``:

           - Determine basic assessments of the data, including the S/N
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

           analysis_plan (:class:`mangadap.par.analysisplan.AnalysisPlanSet`):
               Object with the analysis plans to implement.

           directory_path (str):
                Path to directory containing the MUSE cube file.

           analysis_path (str):
               Top-level directory for the DAP output data; default is
               defined by :func:`mangadap.config.defaults.dap_analysis_path`.

       Returns:
           manga_dap (:func:`mangadap.survey.manga_dap.manga_dap`):
                Wrapper function that runs the DAP entirely.
    """
    # mask the first pixel in the cube
    cube = MUSEDataCube.from_config(config_file)
    cube.mask[0,0,:] = True

    # Run it!
    return manga_dap(cube, analysis_plan, verbose=2, directory_path=directory_path,
                     analysis_path=analysis_path)
#-----------------------------------------------------------------------------


if __name__ == '__main__':

    # file directory path
    file_dir = os.path.dirname(__file__)

    # MUSE cube directory path
    directory_path = os.path.join(file_dir,'data/test_cube_data')

    # Need to make up plate and ifu design numbers
    plate = 100000
    ifudesign = 1

    # file path to the input configuration file
    config_fil = os.path.join(directory_path, 'mangadap-100000-1-LINCUBE.ini')
    # configuration file that implements correlation ratio correction
    config_fil_beta = os.path.join(directory_path, 'mangadap-100000-1-LINCUBE-BETA.ini')

    # Define how you want to analyze the data [include more descriptive definitions]
    plan = AnalysisPlanSet([ AnalysisPlan(drpqa_key='SNRG', # Data reduction quality assessment
                                          # Overwrite existing data-quality assessment reference files
                                          drpqa_clobber=True,
                                          # Spatial binning method keyword
                                          bin_key='SQUARE2.0', # SQUARE bin size options are 0.6, 1.0 and 2.0
                                          # Overwrite any existing spatial binning reference files
                                          bin_clobber=True,
                                          # Stellar-continuum fitting method keyword
                                          continuum_key='MILESHC-NOISM',
                                          # Overwrite any existing stellar-continuum fitting reference files
                                          continuum_clobber=True,
                                          # Emission-line moments measurement method keyword
                                          elmom_key='EMOMMPL11',
                                          # Overwrite any existing emission-line moments reference files
                                          elmom_clobber=True,
                                          # Emission-line modeling method keyword
                                          elfit_key= 'EFITMPL11-ISMMASK-HCNOISM',
                                          # Overwrite any existing emission-line modeling reference files
                                          elfit_clobber=True,
                                          # Spectral-index measurement method keyword
                                          spindex_key='INDXEN',
                                          # Overwrite any existing spectral-index reference files
                                          spindex_clobber=True) ])

    # main output directory
    output_root_dir = os.path.join(file_dir, 'outputs')
    # make the directory if there isn't one already
    if not os.path.isdir(output_root_dir):
        os.makedirs(output_root_dir)

    # output galaxy directory
    gal_name = 'test'
    bin_method = plan['bin_key'][0]
    output_gal_dir = os.path.join(output_root_dir,gal_name+'-'+bin_method)
    if not os.path.isdir(output_gal_dir):
        os.makedirs(output_gal_dir)

    # directories for the corrected and non-corrected MUSE cubes
    no_corr_dir = os.path.join(output_gal_dir,'NO-CORR')
    corr_dir = os.path.join(output_gal_dir,'BETA-CORR')

    # initial run of the DAP
    fit_one_cube_muse(config_fil, plan, directory_path=directory_path, analysis_path=no_corr_dir)
    # quality check
    ppxffit_qa_plot(plate, ifudesign, plan, drpver=None, redux_path=directory_path, dapver=None,
                    analysis_path=no_corr_dir,tpl_flux_renorm=None)

    # final run
    fit_one_cube_muse(config_fil_beta,plan,directory_path=directory_path,analysis_path=corr_dir)

    ppxffit_qa_plot(plate, ifudesign, plan, drpver=None, redux_path=directory_path, dapver=None,
                    analysis_path=corr_dir,tpl_flux_renorm=None)