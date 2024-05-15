#!/usr/bin/env python3

import os
import numpy

import astropy.constants
from astropy.io import fits

from mangadap.datacube import MUSEDataCube
from mangadap.survey.manga_dap import manga_dap
from mangadap.par.analysisplan import AnalysisPlan, AnalysisPlanSet
from mangadap.proc.reductionassessments import available_reduction_assessments
from mangadap.scripts.spotcheck_dap_maps import spotcheck_images
from mangadap.scripts.ppxffit_qa import ppxffit_qa_plot
from mangadap.scripts.fit_residuals_muse import fit_residuals_muse
from mangadap.scripts.manga_dap_inspector import manga_dap_inspector
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
           analysis_plan (:class:`mangadap.par.analysisplan.AnalysisPlanSet`):
               Object with the analysis plans to implement.

           directory_path (:obj:`str`, optional):
               Direct path to directory containing the DRP output file;
               default is defined by
               :func:`mangadap.config.defaults.drp_directory_path`

           analysis_path (:obj:`str`, optional):
               Top-level directory for the DAP output data; default is
               defined by
               :func:`mangadap.config.defaults.dap_analysis_path`.

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

    # directory path for test cube NGC0000
    directory_path = os.path.join(os.getcwd(),'data/test_cube_data/')

    # Need to make up plate and ifu design numbers
    plate = 100000
    ifudesign = 1

    # make sure the directory path is correct within the .ini file
    config_fil = os.path.join(directory_path,'mangadap-100000-1-LINCUBE.ini')

    # Define how you want to analyze the data
    plan = AnalysisPlanSet([ AnalysisPlan(drpqa_key='SNRG',
                                          drpqa_clobber=False,
                                          # square spatial binning size options are 0.6, 1.0 and 2.0
                                          bin_key='SQUARE2.0', #'HYB10'
                                          bin_clobber=False,
                                          continuum_key='MILESHC-NOISM',  #'MILESHCMPL11',
                                          continuum_clobber=False,
                                          elmom_key='EMOMMPL11',
                                          elmom_clobber=False,
                                          elfit_key= 'EFITMPL11-ISMMASK-HCNOISM', #'EFITMPL11SSP', #'EFITMPL9DB',
                                          elfit_clobber=False,
                                          spindex_key='INDXEN',
                                          spindex_clobber=False) ])

    # when re-running this function make sure to delete the output folders in the analysis path
    # in order for the DAP to rerun properly.
    fit_one_cube_muse(config_fil, plan ,directory_path=directory_path, analysis_path='./output2.0_test_no_corr')
    #fit_one_cube_muse(config_file, directory_path=directory_path, analysis_path='./output2.0_test_err_corr')
    #fit_one_cube_muse(config_file, directory_path=directory_path, analysis_path='./output2.0_NGC4030')

    # daptype='SQUARE2.0-MILESHC-MASTARSSP'
    # #daptype = 'SQUARE0.6-MILESHC-MASTARSSP'
    # dapver='2.2.3dev'
    # drpver ='v3_0_1'
    # drpver = '3.1.3dev'
    # #spotcheck_images('./output2.0_test_merge', daptype, plate, ifudesign, ofile=None, drpver=None, dapver=None)

    # redux path is the directory path containing the original unprocessed MUSE cube
    redux_path = directory_path
    ppxffit_qa_plot(plate, ifudesign, plan, drpver=None, redux_path=redux_path, dapver=None, analysis_path='./output2.0_test_no_corr',
                  tpl_flux_renorm=None)

    #fit_residuals_muse(dapver, './output2.0_NGC4030', daptype, plate, ifudesign)
    #manga_dap_inspector(maps_file, model_file, ext=None, masked_spectra=True)