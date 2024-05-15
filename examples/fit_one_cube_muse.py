#!/usr/bin/env python3

import os
import numpy

import astropy.constants
from astropy.io import fits

from mangadap.datacube import MaNGADataCube
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