# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Implements a wrapper class for pPXF.

*License*:
    Copyright (c) 2015, SDSS-IV/MaNGA Pipeline Group
        Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/proc/ppxffit.py

*Imports and python version compliance*:
    ::

        from __future__ import division
        from __future__ import print_function
        from __future__ import absolute_import
        from __future__ import unicode_literals

        import sys
        import warnings
        if sys.version > '3':
            long = int

        import numpy
        import scipy.interpolate
        import scipy.signal
        import astropy.constants

        from ..par.parset import ParSet
        from ..util.bitmask import BitMask
        from ..util.fileio import init_record_array
        from ..util.instrument import spectrum_velocity_scale, resample_vector
        from ..contrib.ppxf import ppxf
        from .spatiallybinnedspectra import SpatiallyBinnedSpectra
        from .templatelibrary import TemplateLibrary
        from .pixelmask import PixelMask
        from .spectralfitting import StellarKinematicsFit

*Class usage examples*:
        Add examples

*Revision history*:
    | **26 Apr 2016**: Moved from spectralfitting.py to its own file by
        K. Westfall (KBW)
    | **05 Jul 2016**: (KBW) V6.0.0 of pPXF does not use the oversample
        keyword given a better solution; see Cappellari (in prep).  This
        keyword was therefore removed from the parameter set.
    | **06 Jul 2016**: (KBW) Use v6.0.0 pPXF functions to compute models
        using new LOSVD kernel functionality.
    | **10 Oct 2016**: (KBW) Fixed error in calculation of velocity
        offset between template and object spectra to account for
        different size pixels.
    | **31 Oct 2016**: (KBW) Allow the spectral resolution to be a
        vector per spaxel or a vector per set of input spectra in
        :func:`PPXFFit.fit`.
    | **01 Nov 2016**: (KBW) Added new iteration method that does not do
        a first fit to the global spectrum but does include the
        rejection iteration.  Allows users to treat each spectrum
        provided to :func:`PPXFFit.fit` individually.  Fixed goodpixel
        mask when not first fitting the global spectrum.
    | **02 Nov 2016**: (KBW) Added ability to limit which templates are
        fit to each spectrum in :func:`PPXFFit.fit` using usetpl kwarg.
    | **17 Feb 2017**: (KBW) Included filtering options.  Changed to use
        ppxf v6.0.4; mpfit object no longer returned, only mpfit status.

.. todo::
    - Include input to :func:`PPXFFit.fit` that allows the user to
      specify a subset of templates to use for each object spectrum.
    - Additional iteration modes?
        - Allow the fit to continue to iterate until the velocity is not
          up against one of the +/- 2000 km/s limits?  Could be useful
          for poor redshift guesses.
    - Keep track of which templates were non-zero in the fit to the
      global spectrum, if it's fit

.. _astropy.io.fits.hdu.hdulist.HDUList: http://docs.astropy.org/en/v1.0.2/io/fits/api/hdulists.html
.. _glob.glob: https://docs.python.org/3.4/library/glob.html
.. _configparser.ConfigParser: https://docs.python.org/3/library/configparser.html#configparser.ConfigParser

"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
import warnings
if sys.version > '3':
    long = int

import logging

import time

import numpy
from scipy import interpolate
import astropy.constants

from ..par.parset import ParSet
from ..util.bitmask import BitMask
from ..util.pixelmask import PixelMask, SpectralPixelMask
from ..util.fileio import init_record_array
from ..util.filter import BoxcarFilter
from ..util.log import log_output
from ..util.instrument import spectrum_velocity_scale, resample_vector, spectral_resolution
from ..util.constants import constants
from ..contrib.ppxf import ppxf, _templates_rfft, _losvd_rfft
from ..contrib import ppxf_util
from .spatiallybinnedspectra import SpatiallyBinnedSpectra
from .templatelibrary import TemplateLibrary
from .spectralfitting import StellarKinematicsFit
# from .spectralstack import SpectralStack
from .util import residual_growth, optimal_scale

from matplotlib import pyplot

__author__ = 'Kyle B. Westfall'
# Add strict versioning
# from distutils.version import StrictVersion


class PPXFFitPar(ParSet):
    r"""

    Define a parameter set used by the pPXF fitting method.

    .. todo::
        The overlap between this and
        :class:`mangadap.proc.stellarcontinuummodel.StellarContinuumModelDef`
        is not well designed.

    Args:
        template_library_key (str): Keyword of the library to fit.  See
            :func:`mangadap.proc.templatelibrary.available_template_libraries`.

        template_library
            (:class:`mangadap.proc.templatelibrary.TemplateLibrary`):
            Object with the spectra in the template library that have
            been prepared for analysis of the data.

        guess_redshift (array-like): Initial guess for the redshift
            (:math:`cz`) of each binned spectrum.

        guess_dispersion (array-like): Initial guess for the velocity
            dispersion for each binned spectrum.

        iteration_mode (str): (**Optional**) Iteration mode to use; see
            :func:`PPXFFit.iteration_modes`.

        match_resolution (bool): (**Optional**) Match the spectral
            resolution of the template to that of the galaxy data.  This
            is used only when constructing the template library.
            Default is True.
        
        velscale_ratio (int): (**Optional**) The **integer** ratio
            between the velocity scale of the pixel in the galaxy data
            to that of the template data.  This is used only when
            constructing the template library.  Default is None, which
            is the same as assuming that the velocity scales are
            identical.
        
        minimum_snr (float): (**Optional**) Minimum S/N ratio to include
            in the fitting.

        pixelmask (:class:`mangadap.proc.pixelmask.PixelMask`):
            (**Optional**) Pixel mask to include during the fitting.
        
        bias, degree, mdegree, moments: (**Optional**) See
            :class:`mangadap.contrib.ppxf.ppxf` documentation.

    """
    def __init__(self, template_library_key, template_library, guess_redshift, guess_dispersion,
                 iteration_mode='global_template', reject_boxcar=None, filter_boxcar=None,
                 filter_iterations=None, match_resolution=True, velscale_ratio=None,
                 minimum_snr=None, pixelmask=None, bias=None, degree=None, mdegree=None,
                 filt_degree=None, filt_mdegree=None, moments=None):
    
        arr = [ numpy.ndarray, list ]                   # sky, p0, lam, component
        arr_in_fl = [ numpy.ndarray, list, int, float ] # component, reg_dim
        in_fl = [ int, float ]                          # Reddening, bias, regul

        _def = self._keyword_defaults()
        
        iter_opt = PPXFFit.iteration_modes()
        moment_opt = [ 2, 4, 6 ]

        pars =     [ 'template_library_key', 'template_library', 'guess_redshift',
                     'guess_dispersion', 'iteration_mode', 'reject_boxcar', 'filter_boxcar',
                     'filter_iterations', 'match_resolution', 'velscale_ratio', 'minimum_snr',
                     'pixelmask', 'bias', 'degree', 'mdegree', 'filt_degree', 'filt_mdegree',
                     'moments' ]
        values =   [ template_library_key, template_library, guess_redshift, guess_dispersion,
                     iteration_mode, reject_boxcar, filter_boxcar, filter_iterations,
                     match_resolution, velscale_ratio, minimum_snr, pixelmask, bias, degree,
                     mdegree, filt_degree, filt_mdegree, moments ]
        options =  [ None, None, None, None, iter_opt, None, None, None, None, None, None, None,
                     None, None, None, None, None, moment_opt ]
        defaults = [ None, None, None, None, 'global_template', None, None, 0, True, None, None,
                     None, _def['bias'], _def['degree'], _def['mdegree'], _def['filt_degree'],
                     _def['filt_mdegree'], _def['moments'] ]
        dtypes =   [ str, TemplateLibrary, arr_in_fl, arr_in_fl, str, int, int, int, bool, int,
                     in_fl, PixelMask, in_fl, int, int, int, int, int ]

        ParSet.__init__(self, pars, values=values, defaults=defaults, options=options,
                        dtypes=dtypes)
        self._check()


    @staticmethod
    def _keyword_defaults():
        """
        Return the keyword defaults.  Pulled from
        :class:`mangadap.contrib.ppxf.ppxf`.
        """
#        return { 'bias':None, 'clean':False, 'degree':4, 'mdegree':0, 'moments':2, 'regul':0,
#                 'reddening':None, 'component':0, 'reg_dim':None }
        return { 'bias':None, 'degree':8, 'mdegree':0, 'filt_degree':8, 'filt_mdegree':0,
                 'moments':2 }


    def _check(self):
        """
        Perform some preliminary checks on the values of the parameters.
        """
        if self['filter_iterations'] > 0 and self['iteration_mode'] != 'fit_reject_filter':
            warnings.warn('Only the \'fit_reject_filter\' iteration mode includes the filtering '
                          'iterations.')
#        if self['filter_iterations'] > 0 and not self['mdegree'] > 0:
#            raise ValueError('If filtering, the multiplicative polynomial must have a non-zero'
#                             ' order (mdegree > 0)')
        
#        if self['reddening'] is not None:
#            if self['mdegree'] > 0:
#                warnings.warn('Cannot both fit multiplicative polynomial and reddening.' \
#                              'Ignoring mdegree.')
#                self['mdegree'] = 0
        # Other checks (and the one above) done within pPXF


        pars =     [ 'template_library_key', 'template_library', 'guess_redshift',
                     'guess_dispersion', 'iteration_mode', 'reject_boxcar', 'filter_boxcar',
                     'filter_iterations', 'match_resolution', 'velscale_ratio', 'minimum_snr',
                     'pixelmask', 'bias', 'degree', 'mdegree', 'filt_degree', 'filt_mdegree',
                     'moments' ]

    def toheader(self, hdr):
        hdr['PPXFTPLK'] = (self['template_library_key'], 'Template library key used with pPXF')
        hdr['PPXFMODE'] = (self['iteration_mode'], 'pPXF iteration mode')
        hdr['PPXFBIAS'] = (str(self['bias']) if self['bias'] is None else self['bias'],
                            'pPXF bias value')
        hdr['PPXFAO'] = (self['degree'], 'Additive order in pPXF')
        hdr['PPXFMO'] = (self['mdegree'], 'Multiplicative order in pPXF')
        hdr['PPXFFAO'] = (self['filt_degree'], 'Additive order for filtered spectra')
        hdr['PPXFFMO'] = (self['filt_mdegree'], 'Multiplicative order for filtered spectra')
        hdr['PPXFMOM'] = (self['moments'], 'Number of fitted LOSVD moments in pPXF')
        if self['reject_boxcar'] is not None:
            hdr['PPXFRBOX'] = (self['reject_boxcar'], 'pPXF rejection boxcar')
        if self['filter_boxcar'] is not None:
            hdr['PPXFFBOX'] = (self['filter_boxcar'], 'pPXF filtering boxcar')
        if self['filter_tterations'] is not None:
            hdr['PPXFFILT'] = (self['filter_iteration'], 'pPXF number of filtering iterations')
        return hdr


    def fromheader(self, hdr):
        self['template_library_key'] = hdr['PPXFTPLK']
        self['bias'] = eval(hdr['PPXFBIAS'])
        self['degree'] = hdr['PPXFAO']
        self['mdegree'] = hdr['PPXFMO']
        self['filt_degree'] = hdr['PPXFFAO']
        self['filt_mdegree'] = hdr['PPXFFMO']
        self['moments'] = hdr['PPXFMOM']

        try:
            self['reject_boxcar'] = hdr['PPXFRBOX']
        except KeyError as e:
            warnings.warn('Input header does not specify rejection boxcar.')
            self['reject_boxcar'] = None

        try:
            self['filter_boxcar'] = hdr['PPXFFBOX']
        except KeyError as e:
            warnings.warn('Input header does not specify filtering boxcar.')
            self['filter_boxcar'] = None

        try:
            self['filter_iteration'] = hdr['PPXFFILT']
        except KeyError as e:
            warnings.warn('Input header does not specify number of filtering iterations.')
            self['filter_iteration'] = None
        


class PPXFFitResult(object):
    """
    A basic utility to save the critical parts of the pPXF model.
    """
    def __init__(self, degree, mdegree, start, end, tpl_to_use, ppxf_fit, ntpl):
        self.start = start
        self.end = end
        self.npixtot = end-start
        self.tpl_to_use = tpl_to_use
        self.status = None if ppxf_fit is None else ppxf_fit.status
        self.gpm = None if ppxf_fit is None else ppxf_fit.goodpixels.copy()
        self.bestfit = None if ppxf_fit is None else ppxf_fit.bestfit.copy()
        self.tplwgt = None if ppxf_fit is None else ppxf_fit.weights[0:ntpl].copy()
        self.addcoef = None if ppxf_fit is None or degree < 0 else ppxf_fit.polyweights.copy()
        self.multcoef = None if ppxf_fit is None or mdegree <= 0 else ppxf_fit.mpolyweights.copy()
        self.kin = None if ppxf_fit is None else ppxf_fit.sol.copy()
        self.kinerr = None if ppxf_fit is None else ppxf_fit.error.copy()
        self.robust_rchi2 = None if ppxf_fit is None else ppxf_fit.chi2

    def empty_fit(self):
        return self.status is None

    def reached_maxiter(self):
        return self.status == 5

    def fit_failed(self):
        return self.empty_fit() or self.status <= 0


class PPXFFit(StellarKinematicsFit):
    """
    Use pPXF to measure the stellar kinematics.  Although it can also
    fit the composition and emission lines, for now just force it to be
    a :class:`StellarKinematicsFit` objec.

    Attributes:
        bitmask (BitMask): Bitmask to use for masking spectral pixels
            and parameter values.  For
            :func:`fit_SpatiallyBinnedSpectra`, must have bit for
            'LOW_SNR'.  For :func:`fit` must have bits for 'TPL_PIXELS',
            'TRUNCATED', 'PPXF_REJECT', 'LARGE_CHI2', 'LARGE_RESID',
            'INSUFFICIENT_DATA', 'FIT_FAILED', 'NEAR_BOUND'.
    """
    def __init__(self, bitmask, par=None):
        StellarKinematicsFit.__init__(self, 'ppxf', bitmask, par=par)
        # Specify some mask bit names
        self.snr_flag = 'LOW_SNR'
        self.rng_flag = 'OUTSIDE_RANGE'
        self.tpl_flag = 'TPL_PIXELS'
        self.trunc_flag = 'TRUNCATED'
        self.rej_flag = 'PPXF_REJECT'
        # Logging and terminal output
        self.loggers = None
        self.quiet = False

        # Imposed fitting boundaries
        self.velocity_limits = None
        self.sigma_limits = None
        self.gh_limits = None

        # Fitting parameters
        self.tpl_wave = None
        self.tpl_sres = None
        self.tpl_flux = None
        self.ntpl = None
        self.obj_wave = None
        self.obj_sres = None
        self.obj_flux = None
        self.obj_ferr = None
        self.input_obj_mask = None
        self.obj_to_fit = None
        self.nobj = None
        self.velscale = None
        self.velscale_ratio = None
        self.matched_resolution = None
        self.guess_kin = None
        self.spectrum_start = None
        self.spectrum_end = None
        self.dof = None
        self.filt_dof = None
        self.base_velocity = None

        # Fitting options
        self.iteration_mode = None
        self.reject_boxcar = None
        self.filter_boxcar = None
        self.filter_iterations = None
        self.bias = None
        self.degree = None
        self.mdegree = None
        self.filt_degree = None
        self.filt_mdegree = None
        self.moments = None


    @staticmethod
    def iteration_modes():
        r"""
        
        Possible iteration methods:

            ``none``: Fit all bins with all templates with a single call
            to pPXF.

            ``fit_reject_filter``: Perform the following procedure:
                - Fit each spectrum
                - for n iterations:
                    - Reject outliers
                    - Filter the object and template spectra
                    - Fit the filtered spectra
                - Fit the unfiltered spectra with the kinematics fixed
                  to result of the final filtered fit

            ``no_global_wrej``: Do not fit the global spectrum
            first, but include a rejection iteration.  All templates are
            fit in each step.

            ``global_template``:  Fit the global spectrum with all
            templates and include a single rejection iteration.  The
            pixel mask for this fit is the base mask for all fits to the
            individual bins.  A single rejection iteration is done for
            each bin.  **Only the global template is used when fitting
            each bin.**

            ``nonzero_templates``:  Fit the global spectrum with
            all templates and include a single rejection iteration.  The
            pixel mask for this fit is the base mask for all fits to the
            individual bins.  A single rejection iteration is done for
            each bin.  **Only the templates with non-zero weights are
            used when fitting each bin.**

            ``all_templates``:  Fit the global spectrum with all
            templates and include a single rejection iteration.  The
            pixel mask for this fit is the base mask for all fits to the
            individual bins.  A single rejection iteration is done for
            each bin.  **All templates are used when fitting each bin.**

            ``global_template_plus_emlines``:  Same as
            ``global_template``, but includes an iteration where the
            emission-lines are simultaneously fit in an effort to fix
            extrapolation errors beneath the each line.

            ``nonzero_templates_plus_emlines``:  Same as
            ``nonzero_templates``, but includes an iteration where the
            emission-lines are simultaneously fit in an effort to fix
            extrapolation errors beneath the each line.

            ``all_templates_plus_emlines``:  Same as ``all_templates``,
            but includes an iteration where the emission-lines are
            simultaneously fit in an effort to fix extrapolation errors
            beneath the each line.

        Returns:
            list: List of allowed options.
        """
        return [ 'none',
                 'fit_reject_filter',
                 'no_global_wrej',
                 'global_template',
                 'nonzero_templates',
                 'all_templates',
                 'global_template_plus_emlines',
                 'nonzero_templates_plus_emlines',
                 'all_templates_plus_emlines'
               ]


    def _mode_uses_global_spectrum(self):
        return self.iteration_mode in [ 'global_template', 'global_template_plus_emlines',
                                        'nonzero_templates', 'nonzero_templates_plus_emlines',
                                        'all_templates', 'all_templates_plus_emlines' ]

    def _mode_uses_global_template(self):
        return self.iteration_mode in [ 'global_template', 'global_template_plus_emlines' ]

            
    def _mode_uses_nonzero_templates(self):
        return self.iteration_mode in [ 'nonzero_templates', 'nonzero_templates_plus_emlines' ]


    def _mode_uses_all_templates(self):
        return self.iteration_mode in [ 'none', 'fit_reject_filter', 'no_global_wrej',
                                        'all_templates', 'all_templates_plus_emlines' ]

    def _mode_includes_emission_line_fit(self):
        return self.iteration_mode in [ 'global_template_plus_emlines',
                                        'nonzero_templates_plus_emlines',
                                        'all_templates_plus_emlines' ]

    def _mode_includes_rejection(self):
        return self.iteration_mode != 'none'

    def _mode_uses_filter(self):
        return self.iteration_mode == 'fit_reject_filter'


    def _check_mode(self, iteration_mode, reject_boxcar, filter_boxcar, filter_iterations, mdegree):
        if iteration_mode not in self.iteration_modes():
            raise ValueError('Do not understand iteration mode \'{0}\''.format(iteration_mode))
        self.iteration_mode = iteration_mode
        self.filter_iterations = 0 if self.iteration_mode != 'fit_reject_filter' \
                                        or filter_iterations is None else filter_iterations
        if self.iteration_mode == 'fit_reject_filter' and filter_boxcar is None:
            warnings.warn('Must provide boxcar for filtering iterations.  Using default (100).')
            self.filter_boxcar = 100
        else:
            self.filter_boxcar = filter_boxcar
        self.reject_boxcar = reject_boxcar
#        if self.filter_iterations > 0 and not mdegree > 0:
#            raise ValueError('If filtering, the multiplicative polynomial must have a non-zero'
#                             ' order (mdegree > 0)')


    def _check_template_usage_flags(self, usetpl):
        if usetpl is not None and usetpl.shape != (self.nobj, self.ntpl) \
                and usetpl.shape != (self.ntpl,):
            raise ValueError('Provided template selection object does not have the correct shape!')
        if usetpl is None:
            self.usetpl = numpy.ones((self.nobj,self.ntpl), dtype=numpy.bool)
        else:
            self.usetpl = usetpl.astype(bool)
            if self.usetpl.shape == (self.ntpl,):
                self.usetpl = numpy.array([self.usetpl]*self.nobj)


    def _check_input_kinematics(self, guess_redshift, guess_dispersion):
        # Get the input guess kinematics
        _guess_redshift = numpy.atleast_1d(guess_redshift)
        if len(_guess_redshift) != 1 and len(_guess_redshift) != self.nobj:
            raise ValueError('Must provide a single redshift or one per object spectrum.')
        if len(_guess_redshift) == 1:
            _guess_redshift = numpy.full(self.nobj, _guess_redshift[0], dtype=numpy.float)
        _guess_dispersion = numpy.atleast_1d(guess_dispersion)
        if len(_guess_dispersion) != 1 and len(_guess_dispersion) != self.nobj:
            raise ValueError('Must provide a single dispersion or one per object spectrum.')
        if len(_guess_dispersion) == 1:
            _guess_dispersion = numpy.full(self.nobj, _guess_dispersion[0], dtype=numpy.float)

        # Set the input redshifts
        self.input_cz = _guess_redshift*astropy.constants.c.to('km/s').value
        # Set the input pPXF, pixel-based kinematics
        # (see convert_velocity)
        self.guess_kin = numpy.array([
                            numpy.log(_guess_redshift+1)*astropy.constants.c.to('km/s').value,
                            _guess_dispersion]).T


    def _check_resolution_match(self, matched_resolution):
        # Confirm there is enough information to handle an unmatched
        # resolution
        if not matched_resolution and (self.obj_sres is None or self.tpl_sres is None):
            raise ValueError('If the spectral resolution is not matched between the template and '
                             'the object data, you must provide the spectral resolution for both.')
        self.matched_resolution = matched_resolution


    def _check_wavelength_range(self, waverange):
        self.waverange = numpy.array([[self.obj_wave[0]-1, self.obj_wave[-1]+1]]) \
                                if waverange is None else numpy.atleast_2d(waverange)
        if self.waverange.shape[0] == 1:
            self.waverange = numpy.array([self.waverange[0,:]]*self.nobj)
        if self.waverange.shape != (self.nobj,2):
            raise ValueError('Input wavelength range array does not have the correct shape.')


    def _initialize_output(self, mask):
        # Initialize the output arrays
        #  Model flux:
        model_flux = numpy.zeros(self.obj_flux.shape, dtype=numpy.float)
        #  Model mask:
        model_mask = numpy.zeros(self.obj_flux.shape, dtype=self.bitmask.minimum_dtype())
        if mask is not None:
            # Include the mask, if provided
            if isinstance(mask, numpy.ndarray):
                if mask is not None and mask.shape != self.obj_flux.shape:
                    raise ValueError('Shape of object mask array must match its flux array.')
                model_mask = mask.astype(self.bitmask.minimum_dtype())
            if isinstance(mask, SpectralPixelMask):
                model_mask = mask.bits(self.bitmask, self.obj_wave, nspec=self.nobj,
                                       velocity_offsets=self.input_cz)

        # Record array with model parameters
        # TODO: Fits using emission lines have hardwired polynomial
        # degrees; Need to fix this!
        model_par = init_record_array(self.nobj,
                        self._per_stellar_kinematics_dtype(self.ntpl, self.degree+1,
                                                           max(self.mdegree,0), self.moments,
                                                           self.bitmask.minimum_dtype())
                            if self.iteration_mode not in [ 'global_template_plus_emlines',
                                                            'nonzero_templates_plus_emlines',
                                                            'all_templates_plus_emlines' ] else
                        self._per_stellar_kinematics_dtype(self.ntpl, 0, 8, self.moments,
                                                           self.bitmask.minimum_dtype()))
        # Set the bins; here the ID and index are identical
        model_par['BINID'] = numpy.arange(self.nobj)
        model_par['BINID_INDEX'] = numpy.arange(self.nobj)

        return model_flux, model_mask, model_par


    def _initialize_pixels_to_fit(self, model_mask, ensemble=True, max_velocity_range=400.,
                                  alias_window=None):

        err = numpy.empty(self.nobj, dtype=object)  # Empty error messages

        # Assess the regions that need to be masked during fitting
        if ensemble:
            print('ensemble is true!')
            self.waverange = numpy.array([numpy.amax(self.waverange[:,0]),
                                          numpy.amin(self.waverange[:,1])])
            try:
                fit_indx, waverange_mask, npix_mask, alias_mask \
                        = self.fitting_mask(self.tpl_wave, self.obj_wave, self.velscale,
                                            velscale_ratio=self.velscale_ratio,
                                            waverange=self.waverange,
                                            velocity_offset=self.input_cz,
                                            max_velocity_range=max_velocity_range,
                                            alias_window=alias_window, loggers=self.loggers,
                                            quiet=self.quiet)
            except ValueError as e:
                return model_mask, numpy.array([e]*self.nobj)

            fit_indx = numpy.array([ fit_indx ]*self.nobj)
            waverange_mask = numpy.array([ waverange_mask ]*self.nobj)
            npix_mask = numpy.array([ npix_mask ]*self.nobj)
            alias_mask = numpy.array([ alias_mask ]*self.nobj)
        else:
            fit_indx = numpy.zeros(self.obj_flux.shape, dtype=bool)
            waverange_mask = numpy.zeros(self.obj_flux.shape, dtype=bool)
            npix_mask = numpy.zeros(self.obj_flux.shape, dtype=bool)
            alias_mask = numpy.zeros(self.obj_flux.shape, dtype=bool)
            for i in range(self.nobj):
                try:
                    fit_indx[i,:], waverange_mask[i,:], npix_mask[i,:], alias_mask[i,:] \
                            = self.fitting_mask(self.tpl_wave, self.obj_wave, self.velscale,
                                                velscale_ratio=self.velscale_ratio,
                                                waverange=self.waverange[i,:],
                                                velocity_offset=self.input_cz[i],
                                                max_velocity_range=max_velocity_range,
                                                alias_window=alias_window, loggers=self.loggers,
                                                quiet=True)#self.quiet)

                except ValueError as e:
                    model_par['MASK'][i] = self.bitmask.turn_on(model_par['MASK'], 'NO_FIT')
                    err[i] = e

        # Add these regions to the mask
        model_mask[waverange_mask] = self.bitmask.turn_on(model_mask[waverange_mask], self.rng_flag)
        model_mask[npix_mask] = self.bitmask.turn_on(model_mask[npix_mask], self.tpl_flag)
        model_mask[alias_mask] = self.bitmask.turn_on(model_mask[alias_mask], self.trunc_flag)

        # Make sure that the errors are valid
        indx = ~(self.obj_ferr.data > 0) | ~(numpy.isfinite(self.obj_ferr.data))
        if numpy.sum(indx) > 0:
            model_mask[indx] = self.bitmask.turn_on(model_mask[indx], 'INVALID_ERROR')
        # To avoid having pPXF throw an exception, make sure that all
        # of these bad errors are set to unity
        self.obj_ferr[indx] = 1.0

        # Update the internal mask of the data
        self.obj_flux[model_mask > 0] = numpy.ma.masked
        self.obj_ferr[model_mask > 0] = numpy.ma.masked

        # Determine the starting and ending pixels
        pix = numpy.ma.MaskedArray(numpy.array([numpy.arange(self.obj_flux.shape[1])]*self.nobj),
                                   mask = ~fit_indx)
        self.spectrum_start = numpy.ma.amin(pix, axis=1)
        self.spectrum_end = numpy.ma.amax(pix, axis=1)+1

        return model_mask, err


    def _run_fit_iteration(self, obj_flux, obj_ferr, start, end, base_velocity, tpl_flux,
                           obj_to_fit=None, tpl_to_use=None, plot=False, fixed_kin=None,
                           degree=None, mdegree=None, dof=None):
        r"""
        Fit all the object spectra in obj_flux.

        Args:
            obj_flux (array): Size is :math:`N_{\rm spec}\times N_{\rm
                chan}`, object spectra
            obj_ferr (array): Size is :math:`N_{\rm spec}\times N_{\rm
                chan}`, object errors
            start (array): Size is :math:`N_{\rm spec}`, starting pixel
                for each spectrum
            end (array): Size is :math:`N_{\rm spec}`, ending pixel (+1)
                for each spectrum
            base_velocity (array): Size is :math:`N_{\rm spec}`, base
                velocity offset between each object spectrum and the
                template spectra.
            tpl_flux (array): Size is :math:`N_{\rm tpl}\times N_{\rm
                tpl chan}`, template spectra
            obj_to_fit (array): (**Optional**) Size is :math:`N_{\rm
                spec}`, boolean flag to fit object spectrum
            tpl_to_use (array): (**Optional**) Size is :math:`N_{\rm
                spec}\times N_{\rm tpl}`, boolean flag to use a template
                for the fit to each object spectrum
            plot (bool): (**Optional**) Produce the default ppxf fit
                plot.

        Returns:
            numpy.ndarray : Array with :math:`N_{\rm spec}` instances of
            :class:`PPXFFitResult`.
        """

        # Get the list of templates to use
        nspec = obj_flux.shape[0]
        _obj_to_fit = numpy.ones(nspec, dtype=bool) if obj_to_fit is None else obj_to_fit
        ntpl = tpl_flux.shape[0]
        _tpl_to_use = numpy.ones((nspec,ntpl), dtype=bool) if tpl_to_use is None else tpl_to_use

        input_kin = self.guess_kin if fixed_kin is None else fixed_kin
        moments = self.moments if fixed_kin is None else -self.moments
        degree = self.degree if degree is None else degree
        mdegree = self.mdegree if mdegree is None else mdegree
        dof = self.dof if dof is None else dof

        linear = fixed_kin is not None and mdegree < 1

        # Create the object to hold all the fits
        result = numpy.empty(nspec, dtype=object)

        # Fit each spectrum
        for i in range(nspec):
            print('Running pPXF fit on spectrum: {0}/{1}'.format(i+1,nspec), end='\r')
            # Meant to ignore this spectrum
            if not _obj_to_fit[i]:
                result[i] = None
                continue

            # Get the pixels to fit for this spectrum
            gpm = numpy.where(~(obj_flux.mask[i,start[i]:end[i]]))[0]

            # Check if there is sufficient data for the fit
            ntpl_to_use = numpy.sum(_tpl_to_use[i,:])
            if len(gpm) < dof+ntpl_to_use:
                if not self.quiet:
                    warnings.warn('Insufficient data points ({0}) to fit spectrum {1}'
                                  '(dof={2}).'.format(len(gpm), i+1, dof+ntpl_to_use))
                result[i] = PPXFFitResult(degree, mdegree, start[i], end[i], _tpl_to_use[i,:],
                                          None, ntpl)
                continue

            # Run ppxf
            if plot:
                pyplot.clf()
            result[i] = PPXFFitResult(degree, mdegree, start[i], end[i], _tpl_to_use[i,:],
                            ppxf(tpl_flux[tpl_to_use[i,:],:].T, obj_flux.data[i,start[i]:end[i]],
                                 obj_ferr.data[i,start[i]:end[i]], self.velscale,
                                 input_kin[i,:], velscale_ratio=self.velscale_ratio,
                                 goodpixels=gpm, bias=self.bias, degree=degree, mdegree=mdegree,
                                 moments=moments, vsyst=-base_velocity[i], quiet=(not plot),
                                 plot=plot, linear=linear), ntpl)
#            print('Fitted kinematics: ', result[i].kin)
            if result[i].reached_maxiter() and not self.quiet:
                warnings.warn('pPXF optimizer reached maximum number of iterations for spectrum '
                              '{0}.'.format(i+1))
            if plot:
                pyplot.show()

        print('Running pPXF fit on spectrum: {0}/{1}'.format(nspec,nspec))
        return result


    def _fit_global_spectrum(self, obj_to_include=None, plot=False):
        """
        Fit the global spectrum.  This:
            - Sets the base-level good pixel mask for the fits to the individual
              bins
            - Gives the template weights for the global template

        .. todo::
            - Only include spectra above a given S/N in global spectrum?
            - Allow for a number of iterations as input.
        """

        _obj_to_include = numpy.ones(self.nobj, dtype=bool) \
                            if obj_to_include is None else obj_to_include
        if len(_obj_to_include) != self.nobj:
            raise ValueError('Incorrect number of object flags.')

        # Create the global spectrum and its error
        # TODO: Does not include covariance!
        global_spectrum = numpy.ma.sum(self.obj_flux[_obj_to_include,:], axis=0).reshape(1,-1)
        global_spectrum_err = numpy.ma.sqrt(numpy.ma.sum(
                                                numpy.square(self.obj_ferr[_obj_to_include,:]),
                                                         axis=0)).reshape(1,-1)
        global_spectrum_err[numpy.ma.getmaskarray(global_spectrum)] = 1.0   # To avoid pPXF error

        # TODO: Because of how it's used, setting start, end, and
        # base_vel this way will mess things up later in the fit()
        # function UNLESS all the spectra have the same start and end;
        # ie., when fitting the global spectrum, the object spectra
        # provided to fit() must be treated as an ensemble.

        # Set the fitting region and base velocity offset
        start = numpy.array([numpy.amax(self.spectrum_start)]).reshape(1,1)
        base_vel = numpy.array([self.base_velocity[numpy.argmax(self.spectrum_start)]]).reshape(1,1)
        end = numpy.array([numpy.amin(self.spectrum_end)]).reshape(1,1)

        # Use any template that is request for any of the individual
        # spectra
        usetpl = numpy.any(self.usetpl, axis=0).reshape(1,-1)

        # Fit the spectrum
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'First fit to global spectrum.')
        result = self._run_fit_iteration(global_spectrum, global_spectrum_err, start, end,
                                         base_vel, self.tpl_flux, tpl_to_use=usetpl, plot=plot)

        # Return if pPXF failed
        if result[0].fit_failed():
            if not self.quiet:
                log_output(self.loggers, 1, logging.INFO, 'pPXF failed to fit global spectrum!')
            return None, None

        # Perform a single fit rejection
        global_spectrum = PPXFFit.reject_model_outliers(global_spectrum, result, rescale=True)

        # refit the spectrum
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'Fit to global spectrum after rejection.')
        return self._run_fit_iteration(global_spectrum, global_spectrum_err, start, end,
                                       base_vel, self.tpl_flux, tpl_to_use=usetpl, plot=plot)[0]


    def _fill_ppxf_par(self, kin, no_shift=True):

        # Moments for each kinematic component
        ncomp = 1
        moments = numpy.atleast_1d(self.moments)
        _moments = numpy.full(ncomp, numpy.absolute(moments), dtype=int) if moments.size == 1 \
                              else numpy.absolute(moments)

        # Construct the LOSVD parameter vector
        vj = numpy.append(0, numpy.cumsum(_moments)[:-1])
        par = kin.copy()
        par[0:2] /= self.velscale
        if no_shift:
            par[0] = 0.0

        return par, _moments, vj


    def _get_losvd_kernels(self, result, no_shift=True):
        
        nspec = len(result)
        losvd_kernel_rfft = numpy.empty(nspec, dtype=object)
        for i in range(nspec):
            if result[i].fit_failed():
                continue
            par, _moments, vj = self._fill_ppxf_par(result[i].kin, no_shift=no_shift)
            losvd_kernel_rfft[i] = _losvd_rfft(par, 1, _moments, vj, self.tpl_npad, 1,
                                        0.0 if no_shift else -self.base_velocity[i]/self.velscale,
                                               self.velscale_ratio, 0.0)[:,0,0]
        return losvd_kernel_rfft
            

    def _matched_mask_filter(self, bf, obj_mask, obj_to_fit, tpl_flux, tpl_to_use, result):
        """
        tpl_rfft and tpl_npad must exist!
        """
        # Get the LOSVD kernels for each fit
        losvd_kernel_rfft = self._get_losvd_kernels(result)
        npix_tpl_resampled = self.npix_tpl // self.velscale_ratio

        # Instantiate the output arrays
        ntpl_per_obj = numpy.sum(tpl_to_use, axis=1)
        tpl_flux_filt = numpy.ma.MaskedArray(numpy.ma.zeros((numpy.sum(ntpl_per_obj),
                                                             self.npix_tpl), dtype=float))
        tpl_to_use_filt = numpy.zeros((self.nobj, numpy.sum(ntpl_per_obj)), dtype=bool)

        # For each spectrum:
        for i in range(losvd_kernel_rfft.shape[0]):
            print('Masking and smoothing templates for object spectrum: {0}/{1}'.format(
                    i+1, losvd_kernel_rfft.shape[0]), end='\r')
#            t = time.clock()

            # Get all the templates convolved by the LOSVD for this fit
            cnvlv_tpl_flux = numpy.fft.irfft(self.tpl_rfft[tpl_to_use[i],:]
                                                    * losvd_kernel_rfft[i][None,:],
                                             axis=1)[:,:self.npix_tpl]
            if self.velscale_ratio > 1:
                cnvlv_tpl_flux = numpy.mean(cnvlv_tpl_flux.reshape(ntpl_per_obj[i], -1,
                                                                   self.velscale_ratio), axis=2)
#            print('fft: time: {0} seconds'.format(time.clock() - t))

            # Get the object-spectrum mask shifted to the template frame
            _obj_mask = numpy.array([obj_mask[i,:]]*ntpl_per_obj[i], dtype=bool)
            shift = numpy.floor((result[i].kin[0]-self.base_velocity[i])/self.velscale).astype(int)
            if shift != 0:
                _obj_mask = numpy.roll(_obj_mask, shift)
                if shift > 0:
                    _obj_mask[:,:shift] = True
                else:
                    _obj_mask[:,shift:] = True
           
            # Apply the mask to the convolved spectra
            _cnvlv_tpl_flux = numpy.ma.MaskedArray(cnvlv_tpl_flux,
                                                   mask=_obj_mask[:,:npix_tpl_resampled])

            # Smooth the template spectra using the same smoothing
            # function as used for the object data
#            t2 = time.clock()
            sm_cnvlv_tpl_flux = bf.smooth(_cnvlv_tpl_flux)
#            print('smooth: time: {0} seconds'.format(time.clock() - t2))

            # Interpolate the smoothing function to the original pixels
            # of the template spectra
#            t2 = time.clock()
            pixcoo = numpy.arange(self.npix_tpl*self.ntpl)
            interpolator = interpolate.interp1d(numpy.mean(pixcoo.reshape(-1,self.velscale_ratio),
                                                           axis=1), sm_cnvlv_tpl_flux.ravel(),
                                                fill_value='extrapolate', assume_sorted=True)
            sm_tpl_flux = interpolator(pixcoo).reshape(self.ntpl, -1)
#            print('interpolate: time: {0} seconds'.format(time.clock() - t2))

            # Set the filtered templates for this object spectrum
            tpls = numpy.sum(ntpl_per_obj[:i])
            tple = numpy.sum(ntpl_per_obj[:i+1])
            tpl_to_use_filt[i,tpls:tple] = True
            tpl_flux_filt[tpls:tple,:] = numpy.ma.divide(self.tpl_flux[tpl_to_use[i],:],
                                                         sm_tpl_flux)
#            print('obj: {0}, total time: {1} seconds'.format(i+1, time.clock() - t))
        
        print('Masking and smoothing templates for object spectrum:              DONE')
        return tpl_flux_filt, tpl_to_use_filt


    def _fit_all_spectra(self, templates, tpl_to_use, plot=False, plot_file_root=None):
        """
        Fit all spectra provided.

        - Get an initial fit
        - Reject
        - Mask and smooth templates and objects
        - Fit ratio
        - Mask and smooth
        - Fit ratio
        - Fit unmasked with fixed kinematics to get models (with
          emission lines?)

        """
        #---------------------------------------------------------------
        # Fit the spectra
        print('Objects to fit: {0}/{1}'.format(numpy.sum(self.obj_to_fit), len(self.obj_to_fit)))
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'Fit to all spectra.')
        result = self._run_fit_iteration(self.obj_flux, self.obj_ferr, self.spectrum_start,
                                         self.spectrum_end, self.base_velocity, templates,
                                         obj_to_fit=self.obj_to_fit, tpl_to_use=tpl_to_use,
                                         plot=plot)
        if not self._mode_includes_rejection():
            # Only a single fit so return
            return result

        #---------------------------------------------------------------
        # Prep for rejection and filtering

        # Copy the input as to not overwrite the input masks
        obj_flux = self.obj_flux.copy()
        obj_ferr = self.obj_ferr.copy()
        obj_to_fit = self.obj_to_fit.copy()

        # Save which were not fit successfully
        obj_to_fit &= numpy.invert(numpy.array([ r is None or r.fit_failed() for r in result ]))
        print('Objects to fit: {0}/{1}'.format(numpy.sum(obj_to_fit), len(obj_to_fit)))

        #---------------------------------------------------------------
        # Reject model outliers
        obj_flux = PPXFFit.reject_model_outliers(obj_flux, result, rescale=False,
                                                 local_sigma=True, boxcar=self.reject_boxcar)
        # Copy the new mask to the errors
        obj_ferr[numpy.ma.getmaskarray(obj_flux)] = numpy.ma.masked

        if self.filter_iterations == 0:
            return self._run_fit_iteration(obj_flux, obj_ferr, self.spectrum_start,
                                           self.spectrum_end, self.base_velocity, templates,
                                           obj_to_fit=obj_to_fit, tpl_to_use=tpl_to_use,
                                           plot=plot)

        #---------------------------------------------------------------
        # Iteratively filter, fit and reject outliers
        print('Filtering. Boxcar size is {0}'.format(self.filter_boxcar))
        bf = BoxcarFilter(self.filter_boxcar)
        for i in range(self.filter_iterations):

            # Get the filtered object spectra and errors
            sm_obj_flux = bf.smooth(obj_flux)
            obj_flux_filt = numpy.ma.divide(obj_flux, sm_obj_flux)
            obj_ferr_filt = numpy.ma.absolute(numpy.ma.divide(obj_ferr, sm_obj_flux))
            obj_ferr_filt[numpy.ma.getmaskarray(obj_ferr_filt)] = 1.0

#            pyplot.imshow(obj_flux, origin='lower', interpolation='nearest', aspect='auto')
#            pyplot.show()
#            pyplot.imshow(numpy.ma.log10(obj_flux_filt), origin='lower', interpolation='nearest',
#                          aspect='auto')
#            pyplot.show()
        
            # Get the filtered template spectra
            tpl_flux_filt, tpl_to_use_filt \
                        = self._matched_mask_filter(bf, numpy.ma.getmaskarray(obj_flux_filt),
                                                    obj_to_fit, templates, tpl_to_use,
                                                    result)

#            pyplot.imshow(templates, origin='lower', interpolation='nearest', aspect='auto')
#            pyplot.show()
#            pyplot.imshow(numpy.ma.log10(tpl_flux_filt), origin='lower', interpolation='nearest',
#                          aspect='auto')
#            pyplot.show()
       
            # Fit the filtered spectra        
            if numpy.sum(numpy.ma.getmaskarray(tpl_flux_filt)) > 0:
                warnings.warn('There are masked template pixels!')
            result = self._run_fit_iteration(obj_flux_filt, obj_ferr_filt, self.spectrum_start,
                                             self.spectrum_end, self.base_velocity,
                                             tpl_flux_filt.data, obj_to_fit=obj_to_fit,
                                             tpl_to_use=tpl_to_use_filt, plot=plot,
                                             degree=self.filt_degree, mdegree=self.filt_mdegree,
                                             dof=self.filt_dof)

            if i == self.filter_iterations - 1:
                break

            # Reject outliers
            obj_flux_filt = PPXFFit.reject_model_outliers(obj_flux_filt, result, rescale=False,
                                                          local_sigma=True,
                                                          boxcar=self.reject_boxcar)
            # Copy the new mask to the unfiltered spectra
            obj_flux[numpy.ma.getmaskarray(obj_flux_filt)] = numpy.ma.masked
            obj_ferr[numpy.ma.getmaskarray(obj_flux)] = numpy.ma.masked

        # Save the best-fit kinematics and errors
        best_fit_kin = numpy.array([ r.kin for r in result ])
        best_fit_kinerr = numpy.array([ r.kinerr for r in result ])

        # Refit the unfiltered spectra with the fixed kinematics and
        # with a multiplicative polynomial
        result = self._run_fit_iteration(obj_flux, obj_ferr, self.spectrum_start,
                                         self.spectrum_end, self.base_velocity, templates,
                                         obj_to_fit=obj_to_fit, tpl_to_use=tpl_to_use, plot=plot,
                                         fixed_kin=best_fit_kin)

        # Replace the kinematic errors
        for i in range(len(result)):
            result[i].kinerr = best_fit_kinerr[i,:]

        return result


#    def _rescale_to_fit_emission_lines(self, ppxf_fit, i, global_weights, gpm, dvtol):
#
#        # - updated mask with emission-lines unmasked
#        # - Keep stellar kinematics fixed
#
#        # Construct a set of Gaussian emission line templates
#        z = (numpy.exp(ppxf_fit.sol[0]/astropy.constants.c.to('km/s').value)-1.0)
#        obj_rest_waverange = numpy.array([self.obj_wave[self.spectrum_start[i]],
#                                          self.obj_wave[self.spectrum_end[i]-1] ]) / (1+z)
#        
#        # TODO: Just adopting FWHM = 2.5 for now
#        gas_templates, line_names, line_wave = \
#                ppxf_util.emission_lines(numpy.log(self.tpl_wave), obj_rest_waverange, 2.5)
##        print('GAS')
##        print( gas_templates.shape )
##        print( numpy.sum(gas_templates, axis=0) )
##        for ii in range(gas_templates.shape[1]):
##            pyplot.plot(numpy.arange(gas_templates.shape[0]), gas_templates[:,ii])
##        pyplot.show()
#
#        # Construct the best-fitting stellar template based on the first
#        # fit
#
#        # TODO: Need to check that this works with the change in how the
#        # templates are selected for individual spectra
#        if self.iteration_mode == 'global_template_plus_emlines':
#            global_weights *= ppxf_fit.weights[0]
##        elif self.iteration_mode == 'nonzero_templates_plus_emlines':
##            global_weights[tpl_to_fit[i,:]] = ppxf_fit.weights[0:ntpl_to_fit]
#            #templates.shape[0]]
#        else:
#            global_weights[tpl_to_fit[i,:]] = ppxf_fit.weights[0:ntpl_to_fit]
##            global_weights = ppxf_fit.weights[0:self.ntpl]
#        global_template = numpy.dot(global_weights, self.tpl_flux).reshape(1,-1)
#
##        pyplot.plot(self.tpl_wave, global_template[0,:])
##        pyplot.show()
#
#        # Include emission-line templates with others
#        _templates = numpy.row_stack([global_template, gas_templates.T])
##        for ii in range(_templates.shape[0]):
##            pyplot.plot(numpy.arange(_templates.shape[1]), _templates[ii,:])
##        pyplot.show()
#
#        # Set the component values
#        _component = [0]*1 + [1]*gas_templates.shape[1]
#        # do not fit stars but use previous solution
#        _moments = [-self.moments, 2]
#        # initialize the gas kinematics based on the stars
#        _guess_kin = [ppxf_fit.sol[0:2], [ppxf_fit.sol[0], 50]]
##        _guess_kin = [ppxf_fit.sol[0:2], ppxf_fit.sol[0:2]]
#        # unmask the emission lines
#        _gpm = numpy.union1d(gpm, numpy.where(self.bitmask.flagged(
#                                        model_mask[i,self.spectrum_start[i]:self.spectrum_end[i]],
#                                        flag='EML_REGION'))[0])
#        # Fit with emission lines, using only multiplicative polynomials
#        _ppxf_fit = ppxf(_templates.T,
#                         self.obj_flux.data[i,self.spectrum_start[i]:self.spectrum_end[i]],
#                         self.obj_ferr.data[i,self.spectrum_start[i]:self.spectrum_end[i]],
#                         self.velscale, _guess_kin, velscale_ratio=self.velscale_ratio,
#                         goodpixels=_gpm, bias=self.bias, clean=self.clean, degree=-1,
#                         mdegree=8, moments=_moments, vsyst=-self.base_velocity[i],
#                         quiet=(not plot), plot=plot, component=_component)
#        if plot:
#            pyplot.show()
#
#        # Fit failed
#        if not ppxf_fit.status > 0:
#            model_mask[i,:], model_par['MASK'][i] = self._set_and_report_failed_status(
#                            model_mask[i,:], model_par['MASK'][i],
#                           'Emission-line iteration pPXF status for spectrum {0}; '
#                           'nothing saved.'.format(i+1))
#            continue
#
#        # Save the necessary new data to the old fit (without
#        # emission lines)
#        ppxf_fit.weights = global_weights * _ppxf_fit.weights[0]
#        ppxf_fit.polyweights = None
#        ppxf_fit.mpolyweights = _ppxf_fit.mpolyweights
#
##        old_bestfit = ppxf_fit.bestfit
#        ppxf_fit.bestfit = self.reconstruct_model(self.tpl_wave, self.tpl_flux,
#                                self.obj_wave, ppxf_fit.sol, ppxf_fit.weights[0:self.ntpl],
#                                self.velscale, polyweights=ppxf_fit.polyweights,
#                                mpolyweights=ppxf_fit.mpolyweights,
#                                start=self.spectrum_start[i], end=self.spectrum_end[i],
#                                velscale_ratio=self.velscale_ratio, dvtol=dvtol,
#                                revert_velocity=False)[self.spectrum_start[i]:self.spectrum_end[i]]
##        pyplot.plot(self.obj_wave[self.spectrum_start:self.spectrum_end], ppxf_fit.bestfit)
##        pyplot.plot(self.obj_wave[self.spectrum_start:self.spectrum_end], _ppxf_fit.bestfit)
##        pyplot.plot(self.obj_wave[self.spectrum_start:self.spectrum_end], old_bestfit)
##        pyplot.show()
#        return ppxf_fit


    def _dispersion_correction(self, obj_sres, gpm):
        """
        Calculate the dispersion corrections as the quadrature
        difference between the spectral resolution of the template and
        object spectra.
        """
        cnst = constants()
        fwhm_inst_obj = astropy.constants.c.to('km/s').value/obj_sres[gpm]
        fwhm_inst_tpl = astropy.constants.c.to('km/s').value/self.tpl_sres(self.obj_wave[gpm])
        return numpy.sqrt(numpy.mean(numpy.square(fwhm_inst_obj)
                                     - numpy.square(fwhm_inst_tpl)))/cnst.sig2fwhm


    def _is_near_bounds(self, result, guess_velocity, tol_frac=1e-2):
        """
        Check if the fitted kinematics are near the imposed limits.
        
        The definition of "near" is that the velocity and higher moments
        cannot be closer than the provided fraction of the total width
        to the boundary.  For the velocity dispersion, the fraction is
        done in log space.
        """

        near_bounds = numpy.zeros(self.moments, dtype=numpy.bool)
        near_lower_sigma_bound = False

        # Velocity
        _velocity_limits = self.velocity_limits + guess_velocity
        Dv = numpy.diff(_velocity_limits)[0]
        tol = Dv*tol_frac

        near_bounds[0] = result.kin[0]-_velocity_limits[0] < tol \
                            or _velocity_limits[1]-result.kin[0] < tol

        # Velocity dispersion
        Ds = numpy.diff(numpy.log10(self.sigma_limits))[0]
        tol = Ds*tol_frac

        near_lower_sigma_bound \
                = numpy.log10(result.kin[1]) - numpy.log10(self.sigma_limits[0]) < tol
        near_bounds[1] = near_lower_sigma_bound \
                        or numpy.log10(self.sigma_limits[1]) - numpy.log10(result.kin[1]) < tol

        if self.moments == 2:
            return near_bounds, near_lower_sigma_bound

        # H3 and H4
        Dh = numpy.diff(self.gh_limits)
        tol = Dh*tol_frac
        near_bounds[2] = result.kin[2] - self.gh_limits[0] < tol \
                                or self.gh_limits[1] - result.kin[2] < tol
        near_bounds[3] = result.kin[3] - self.gh_limits[0] < tol \
                                or self.gh_limits[1] - result.kin[3] < tol
        
        if self.moments == 4:
            return near_bounds, near_lower_sigma_bound

        # H5 and H6
        near_bounds[4] = result.kin[4] - self.gh_limits[0] < tol \
                                or self.gh_limits[1] - result.kin[4] < tol
        near_bounds[5] = result.kin[5] - self.gh_limits[0] < tol \
                                or self.gh_limits[1] - result.kin[5] < tol

        return near_bounds, near_lower_sigma_bound


#    def _set_and_report_failed_status(self, model_mask, model_par, message):
#        if not self.quiet:
#            log_output(self.loggers, 1, logging.INFO, message)
#        return self.bitmask.turn_on(model_mask, 'FIT_FAILED'), \
#               self.bitmask.turn_on(model_par, 'FIT_FAILED')


    def _validate_kinematics(self, model_mask, model_par):
        """
        Validate the returned kinematics.

        Checks:
            - corrected velocity dispersion must be in the range 50-400
              km/s
        """
        sigcor = numpy.square(model_par['KIN'][:,1]) - numpy.square(model_par['SIGMACORR'])
        indx = (sigcor < 2500.) | (sigcor > 1.6e5)
        if numpy.sum(indx) == 0:
            return
        model_mask[indx,:] = self.bitmask.turn_on(model_mask[indx,:], 'BAD_SIGMA')
        model_par['MASK'][indx] = self.bitmask.turn_on(model_par['MASK'][indx], 'BAD_SIGMA')


    def _save_results(self, global_fit_result, result, model_flux, model_mask, model_par):

        #---------------------------------------------------------------
        # Get the model spectra
        model_flux = PPXFFit.compile_model_flux(self.obj_flux, result)
        # Calculate the model residuals; both should be masked where the
        # data were not fit
        residual = self.obj_flux - model_flux
        fractional_residual = numpy.ma.divide(self.obj_flux - model_flux, model_flux)
        # Get the chi-square for each spectrum
        model_par['CHI2'] = numpy.sum(numpy.square(residual/self.obj_ferr), axis=1)
        # Get the (fractional) residual RMS for each spectrum
        model_par['RMS'] = numpy.sqrt(numpy.ma.mean(numpy.square(residual), axis=1))
        model_par['FRMS'] = numpy.sqrt(numpy.ma.mean(numpy.square(fractional_residual), axis=1))

        # Flag the pixels that were not used
        model_mask[numpy.ma.getmaskarray(self.obj_flux)] \
                        = self.bitmask.turn_on(model_mask[numpy.ma.getmaskarray(self.obj_flux)],
                                               flag='DIDNOTUSE')

        #---------------------------------------------------------------
        # Need to iterate over each spectrum
        for i in range(self.nobj):

            #-----------------------------------------------------------
            # Set output flags
            # No fit was performed
            if result[i] is None:
                model_mask[i,:] = self.bitmask.turn_on(model_mask[i,:], 'NO_FIT')
                continue

            # No fit attempted because of insufficient data
            if result[i].empty_fit():
                model_mask[i,:] = self.bitmask.turn_on(model_mask[i,:], 'NO_FIT')
                model_par['MASK'][i] = self.bitmask.turn_on(model_par['MASK'][i],
                                                            'INSUFFICIENT_DATA')
                continue

            # Fit attempted but failed
            if result[i].fit_failed():
                model_mask[i,:] = self.bitmask.turn_on(model_mask[i,:], 'FIT_FAILED')
                model_par['MASK'][i] = self.bitmask.turn_on(model_par['MASK'][i], 'FIT_FAILED')

            # Fit successful but hit maximum iterations.
            if result[i].reached_maxiter():
                model_par['MASK'][i] = self.bitmask.turn_on(model_par['MASK'][i], 'MAXITER')

            # Test if the kinematics are near the imposed boundaries.
            # The best fitting model is still provided, but masked.
            near_bound, near_lower_sigma_bound = self._is_near_bounds(result[i],
                                                                      self.guess_kin[i,0])

            # If the velocity dispersion has hit the lower limit, ONLY
            # flag the value as having a MIN_SIGMA.
            if near_lower_sigma_bound:
                near_bound[1] = False
                model_par['MASK'][i] = self.bitmask.turn_on(model_par['MASK'][i], 'MIN_SIGMA')
            # Otherwise, flag both the model and parameter set
            if numpy.any(near_bound):
                if not self.quiet:
                    warnings.warn('Returned parameters for spectrum {0} too close to '
                                  'bounds.'.format(i+1))
                model_par['MASK'][i] = self.bitmask.turn_on(model_par['MASK'][i], 'NEAR_BOUND')
                model_mask[i,:] = self.bitmask.turn_on(model_mask[i,:], 'NEAR_BOUND')

            # Mask rejected pixels
            original_gpm = numpy.where(numpy.invert(
                               numpy.ma.getmaskarray(self.obj_flux)[i,self.spectrum_start[i]
                                                                      :self.spectrum_end[i]]))[0]
            rejected_pixels = list(set(original_gpm) - set(result[i].gpm))
            if len(rejected_pixels) > 0:
                model_mask[i,self.spectrum_start[i]:self.spectrum_end[i]][rejected_pixels] \
                        = self.bitmask.turn_on(model_mask[i,self.spectrum_start[i]:
                                                            self.spectrum_end[i]][rejected_pixels],
                                               flag=self.rej_flag)

            #-----------------------------------------------------------
            # Save the model parameters and figures of merit
            # Number of fitted pixels
            model_par['NPIXFIT'][i] = len(result[i].gpm)
            # Template weights
            if self._mode_uses_global_spectrum():
                model_par['TPLWGT'][i] = result[i].tplwgt[0]*global_fit_result.tplwgt
            else:
                model_par['TPLWGT'][i][result[i].tpl_to_use] = result[i].tplwgt
            # Templates used
            model_par['USETPL'][i] = result[i].tpl_to_use
            # Additive polynomial coefficients
            if self.degree >= 0 and result[i].addcoef is not None:
                model_par['ADDCOEF'][i] = result[i].addcoef
            if self.mdegree > 0 and result[i].multcoef is not None:
                model_par['MULTCOEF'][i] = result[i].multcoef
            # Input kinematics
            model_par['KININP'][i] = self.guess_kin[i,:]
            # Best-fit kinematics
            model_par['KIN'][i] = result[i].kin
            # Kinematic errors
            model_par['KINERR'][i] = result[i].kinerr
            # Chi-square
            resid = self.obj_flux[i,self.spectrum_start[i]:self.spectrum_end[i]] - result[i].bestfit
            model_par['RCHI2'][i] = model_par['CHI2'][i] \
                                        / (model_par['NPIXFIT'][i] 
                                            - self.dof - numpy.sum(model_par['TPLWGT'][i] > 0))
            model_par['ROBUST_RCHI2'][i] = result[i].robust_rchi2

            # Get growth statistics for the residuals
            model_par['ABSRESID'][i] = residual_growth((residual[i,:]).compressed(),
                                                       [0.68, 0.95, 0.99])
            model_par['FABSRESID'][i] = residual_growth(fractional_residual[i,:].compressed(),
                                                        [0.68, 0.95, 0.99])

            # Calculate the dispersion correction if necessary
            if not self.matched_resolution:
                model_par['SIGMACORR'][i] = self._dispersion_correction(self.obj_sres[i],
                                                                        result[i].gpm)

        #---------------------------------------------------------------
        # Test if kinematics are reliable
        self._validate_kinematics(model_mask, model_par)

        #---------------------------------------------------------------
        # Convert the velocities from pixel units to redshift
        model_par['KININP'][:,0], _ = PPXFFit.convert_velocity(model_par['KININP'][:,0],
                                                               numpy.zeros(self.nobj))
        model_par['KIN'][:,0], model_par['KINERR'][:,0] \
                = PPXFFit.convert_velocity(model_par['KIN'][:,0], model_par['KINERR'][:,0])

        #---------------------------------------------------------------
        return model_flux, model_mask, model_par


    def fit_SpatiallyBinnedSpectra(self, binned_spectra, par=None, loggers=None, quiet=False):
        """

        This is a basic interface that is geared for the DAP that
        interacts with the rest of the, more general, parts of the
        class.

        This should not declare anything to self!

        """
        # Parameter must be provided
        if par is None:
            raise ValueError('Required parameters for PPXFFit have not been defined.')
        # Check the parameters
        _def = PPXFFitPar._keyword_defaults()
#        if par['regul'] != _def['regul'] or par['reddening'] != _def['reddening'] \
#                or par['component'] != _def['component'] or par['reg_dim'] != _def['reg_dim']:
#            raise NotImplementedError('Cannot use regul, reddening, component, or regul_dim yet.')

        # SpatiallyBinnedSpectra object always needed
        if binned_spectra is None or not isinstance(binned_spectra, SpatiallyBinnedSpectra):
            raise TypeError('Must provide a SpatiallyBinnedSpectra object for fitting.')
        if binned_spectra.hdu is None:
            raise ValueError('Provided SpatiallyBinnedSpectra object is undefined!')

        # TemplateLibrary object always needed
        if par['template_library'] is None \
                or not isinstance(par['template_library'], TemplateLibrary):
            raise TypeError('Must provide a TemplateLibrary object for fitting.')
        if par['template_library'].hdu is None:
            raise ValueError('Provided TemplateLibrary object is undefined!')

        # Select the spectra that meet the selection criteria
        # TODO: Link this to the StellarContinuumModel._bins_to_fit()
        # function...
        good_spec = binned_spectra.above_snr_limit(par['minimum_snr'])

        # Get the object data
        obj_wave = binned_spectra['WAVE'].data.copy()
        obj_flux = binned_spectra.copy_to_masked_array(flag=binned_spectra.do_not_fit_flags())
        obj_ferr = numpy.ma.power(binned_spectra.copy_to_masked_array(ext='IVAR',
                                                    flag=binned_spectra.do_not_fit_flags()) , -0.5)
        obj_sres = binned_spectra.copy_to_array(ext='SPECRES')

        # Warn the user that only a single spectral resolution is used
        # for the templates
        if not self.quiet:
            warnings.warn('Adopting mean spectral resolution of all templates!')
        tpl_sres = numpy.mean(par['template_library']['SPECRES'].data, axis=0).ravel()

        # Perform the fit
        model_wave, model_flux, model_mask, model_par \
                = self.fit(par['template_library']['WAVE'].data.copy(),
                           par['template_library']['FLUX'].data.copy(),
                           binned_spectra['WAVE'].data.copy(), obj_flux[good_spec,:],
                           obj_ferr[good_spec,:], par['guess_redshift'][good_spec],
                           par['guess_dispersion'][good_spec], iteration_mode=par['iteration_mode'],
                           velscale_ratio=par['velscale_ratio'], mask=par['pixelmask'],
                           matched_resolution=par['match_resolution'],
                           tpl_sres=tpl_sres, obj_sres=obj_sres[good_spec,:],
                           waverange=par['pixelmask'].waverange, bias=par['bias'],
                           degree=par['degree'], mdegree=par['mdegree'],
                           filt_degree=par['filt_degree'], filt_mdegree=par['filt_mdegree'],
                           moments=par['moments'], loggers=loggers, quiet=quiet, dvtol=1e-9)

        # DEBUG
        assert numpy.sum(model_wave - binned_spectra['WAVE'].data) == 0, \
                    'Incorrect wavelength range'

        # Save the the bin ID numbers indices based on the spectra
        # selected to be fit
        model_par['BINID'] = binned_spectra['BINS'].data['BINID'][good_spec]
        model_par['BINID_INDEX'] = numpy.arange(binned_spectra.nbins)[good_spec]

        # Only return model and model parameters for the *fitted*
        # spectra
        return model_wave, model_flux, model_mask, model_par

    
    def fit(self, tpl_wave, tpl_flux, obj_wave, obj_flux, obj_ferr, guess_redshift,
            guess_dispersion, iteration_mode='global_template', reject_boxcar=100, 
            filter_boxcar=100, filter_iterations=0, ensemble=True, velscale_ratio=None, mask=None,
            usetpl=None, matched_resolution=True, tpl_sres=None, obj_sres=None, waverange=None,
            bias=None, degree=4, mdegree=0, filt_degree=4, filt_mdegree=0, moments=2, loggers=None,
            quiet=False, max_velocity_range=400., alias_window=None, dvtol=1e-10, plot=False,
            plot_file_root=None):
        """
        Wrapper for pPXF with some additional convenience functions.
        Limited implementation at the moment.

        Args:
            tpl_wave (numpy.ndarray): 1D vector of template wavelengths
                at rest in angstroms.
            tpl_flux (numpy.ndarray): N-templates x N-wavelengths array
                of template spectra to fit.
            obj_wave (numpy.ndarray): 1D vector of object wavelengths in
                angstroms.  Does NOT need to be same as the template
                wavelengths.
            obj_flux (numpy.ndarray): N-spec x N-wavelengths array
                of object spectra to fit.
            obj_ferr (numpy.ndarray): N-spec x N-wavelengths array
                with the errors in the object spectra.
            guess_redshift (float or numpy.ndarray): Single or
                spectrum-specific redshift used to set the initial guess
                kinematics.
            guess_dispersion (float or numpy.ndarray): Single or
                spectrum-specific velocity dispersion used to set the
                initial guess kinematics.
            iteration_mode (str): (**Optional**) Iteration sequence to
                perform.  See :func:`iteration_modes`.
            reject_boxcar (int): (**Optional**) Size of the boxcar to
                use during the rejection iteration.  Default is 100.  If
                None, rejection uses the entire residual spectrum.
            filter_boxcar (int): (**Optional**) Size of the boxcar to
                use when filtering the spectra.  Default is 100.  Cannot
                be None.
            filter_iterations (int): (**Optional**) Number of filtering
                iterations for the 'fit_reject_filter' iteration mode.
            ensemble (bool): (**Optional**) Treat the list of input
                spectra as an ensemble.  Currently, this only affects
                how the spectra are masked.  Default is to treat them as
                an ensemble.  When not treated as an ensemble, each
                spectrum is masked individually according to the input
                wavelength range and velocity offsets.  *It does not
                make sense to set the iteration_mode to something that
                will include a fit to the global spectrum if you're not
                treating the list of object spectra as an ensemble.*
            velscale_ratio (int): (**Optional**) Ratio of velocity scale
                per pixel in the object spectra to that for the template
                spectra.  Default is that they should be identical.
            mask (numpy.ndarray): (**Optional**) A
                baseline pixel mask to use during the fitting  Other
                pixels may be masked via the convenience functions, but
                these pixels will always be masked.
            matched_resolution (bool): (**Optional**)  Flag that the
                object and template spectra have identical spectral
                resolution.  Default is True.
            tpl_sres (numpy.ndarray): (**Optional**) One-dimensional
                vector with the spectral resolution (:math:`R =
                \lambda/\Delta\lambda`) at each wavelength of the
                template spectra.  Default is the resolution is not
                provided and assumed to be same as the object
                resolution.
            obj_sres (numpy.ndarray): (**Optional**) One- or Two-dimensional
                array with the spectral resolution (:math:`R =
                \lambda/\Delta\lambda`) at each wavelength for (each of)
                the object spectra.  Default is the resolution is not
                provided and assumed to be same as the template
                resolution.
            waverange (array-like): (**Optional**) Lower and upper
                wavelength limits to *include* in the fit.  This can be
                a two-element vector to apply the same limits to all
                spectra, or a N-spec x 2 array with wavelength ranges
                for each spectrum to be fit.  Default is to use as much
                of the spectrum as possible.
            bias (float): (**Optional**) Defaults to 0.0. From the pPXF
                documentation: This parameter biases the (h3, h4, ...)
                measurements towards zero (Gaussian LOSVD) unless their
                inclusion significantly decreases the error in the fit.
                Set this to BIAS=0.0 not to bias the fit: the solution
                (including [V, sigma]) will be noisier in that case. The
                default BIAS should provide acceptable results in most
                cases, but it would be safe to test it with Monte Carlo
                simulations. This keyword precisely corresponds to the
                parameter \lambda in the Cappellari & Emsellem (2004)
                paper. Note that the penalty depends on the *relative*
                change of the fit residuals, so it is insensitive to
                proper scaling of the NOISE vector. A nonzero BIAS can
                be safely used even without a reliable NOISE spectrum,
                or with equal weighting for all pixels.
            degree (int): (**Optional**) Default is 4.  From the pPXF
                documentation: degree of the *additive* Legendre
                polynomial used to correct the template continuum shape
                during the fit (default: 4).  Set DEGREE = -1 not to
                include any additive polynomial.
            mdegree (int): (**Optional**) Default is 0.  From the pPXF
                documentation: degree of the *multiplicative* Legendre
                polynomial (with mean of 1) used to correct the
                continuum shape during the fit (default: 0). The zero
                degree multiplicative polynomial is always included in
                the fit as it corresponds to the weights assigned to the
                templates.  Note that the computation time is longer
                with multiplicative polynomials than with the same
                number of additive polynomials.

                .. note::
                
                    **IMPORTANT**: Multiplicative polynomials cannot be
                    used when the REDDENING keyword is set.

            filt_degree (int): (**Optional**) The order of the additive
                polynomial to use when fitting the filtered spectra.
            filt_mdegree (int): (**Optional**) The order of the
                multiplicative polynomial to use when fitting the
                filtered spectra.
            moments (int): (**Optional**) Default is 2.  From the pPXF
                documentation: Order of the Gauss-Hermite moments to
                fit. Set this keyword to 4 to fit [h3, h4] and to 6 to
                fit [h3, h4, h5, h6]. Note that in all cases the G-H
                moments are fitted (non-linearly) *together* with [V,
                sigma].

                    - If MOMENTS=2 or MOMENTS is not set then only [V,
                      sigma] are fitted and the other parameters are
                      returned as zero.
                    - If MOMENTS is negative then the kinematics of the
                      given COMPONENT are kept fixed to the input
                      values.
                    - EXAMPLE: We want to keep fixed component 0, which
                      has an LOSVD described by [V, sigma, h3, h4] and
                      is modelled with 100 spectral templates; At the
                      same time we fit [V, sigma] for COMPONENT=1, which
                      is described by 5 templates (this situation may
                      arise when fitting stellar templates with
                      pre-determined stellar kinematics, while fitting
                      the gas emission).  We should give in input to
                      ppxf() the following parameters: component =
                      [0]*100 + [1]*5   # --> [0, 0, ..., 0, 1, 1, 1, 1,
                      1] moments = [-4, 2] start = [[V, sigma, h3, h4],
                      [V, sigma]]

            loggers (list): (**Optional**) List of `logging.Logger`_ objects
                to log progress; ignored if quiet=True.  Logging is done
                using :func:`mangadap.util.log.log_output`.  Default is
                no logging.
            quiet (bool): (**Optional**) Suppress all terminal and
                logging output.  Default is False.
            max_velocity_range (float): (**Optional**) Maximum range
                (+/-) expected for the fitted velocities in km/s.
                Default is 400 km/s.
            alias_window (float) : (**Optional**) The window to mask to
                avoid aliasing near the edges of the spectral range in
                km/s.  Default is six times *max_velocity_range*.
            dvtol (float): (**Optional**) The velocity scale of the
                template spectra and object spectrum must be smaller
                than this tolerance.  Default is 1e-10.
            plot (bool): (**Optional**) Show the automatically generated
                pPXF fit plots.
                
        Returns:
            numpy.ndarray: Returns 4 objects:

                1. The wavelengths of the best fitting model spectra.
                Nominally the same as the wavelengths of the input
                object spectra (*obj_wave*).

                2. The fluxes of the best-fitting model spectra.

                3. A mask for the best-fitting models spectra, following
                from the internal bitmask.

                4. A record array with the fitted model parameters; see
                :class:`spectralfitting.StellarKinematicsFit._per_stellar_kinematics_dtype`.

        Raises:
            ValueError: Raised if the input arrays are not of the
                correct shape or if the pixel scale of the template and
                object spectra is greater than the specified tolerance.
        """
        #---------------------------------------------------------------
        # Initialize the reporting
        if loggers is not None:
            self.loggers = loggers
        self.quiet = quiet

        #---------------------------------------------------------------
        # Check the input; the following operations are necessarily
        # sequential in some cases because of attributes assigned along
        # the way

        # - Mode
        self._check_mode(iteration_mode, reject_boxcar, filter_boxcar, filter_iterations, mdegree)
        # - Templates
        self.tpl_wave, self.tpl_flux, self.tpl_sres \
                = PPXFFit.check_templates(tpl_wave, tpl_flux, tpl_sres=tpl_sres)
        self.ntpl, self.npix_tpl = self.tpl_flux.shape
        # - Objects
        self.obj_wave, self.obj_flux, self.obj_ferr, self.obj_sres \
                = PPXFFit.check_objects(obj_wave, obj_flux, obj_ferr=obj_ferr, obj_sres=obj_sres)
        self.nobj, self.npix_obj = self.obj_flux.shape
        self.input_obj_mask = numpy.ma.getmaskarray(self.obj_flux).copy()
        self.obj_to_fit = numpy.any(numpy.invert(self.input_obj_mask), axis=1)
        print('Objects to fit: {0}/{1}'.format(numpy.sum(self.obj_to_fit), self.nobj))
        # - Template usage (needs to know the number of object spectra)
        self._check_template_usage_flags(usetpl)
        # - Input kinematics
        self._check_input_kinematics(guess_redshift, guess_dispersion)
        # - Input pixel scale
        self.velscale, self.velscale_ratio \
                = PPXFFit.check_pixel_scale(self.tpl_wave, self.obj_wave,
                                            velscale_ratio=velscale_ratio, dvtol=dvtol)
        # Force the number of template pixels to be an integer number of
        # object pixels
        if self.velscale_ratio > 1:
            self.npix_tpl -= self.npix_tpl % self.velscale_ratio
            self.tpl_flux = self.tpl_flux[:,:self.npix_tpl]
            self.tpl_wave = self.tpl_wave[:self.npix_tpl]
#            self.tpl_sres = self.tpl_sres[:self.npix_tpl]

        # - Spectral resolution
        self._check_resolution_match(matched_resolution)
        # - Selected wavelength range: always has shape (self.nobj,2)
        self._check_wavelength_range(waverange)

        #---------------------------------------------------------------
        # Get the real FFT of the templates if filtering
        if self._mode_uses_filter():
            _tpl_rfft, self.tpl_npad = _templates_rfft(self.tpl_flux.T)
            self.tpl_rfft = _tpl_rfft.T
        else:
            self.tpl_rfft = None
            self.tpl_npad = None

        #---------------------------------------------------------------
        # Make sure that any mode that fits the global spectrum treats
        # the individual spectra as part of an ensemble.  This step is
        # important for setting the spectral range and the masking that
        # is assumed throughout the rest of the function!
        if not ensemble and self._mode_uses_global_spectrum():
            warnings.warn('When fitting the global spectrum, the spectra MUST be treated as an '
                          'ensemble.')
            _ensemble = True
        else:
            _ensemble = ensemble

        #---------------------------------------------------------------
        # Save the basic pPXF parameters
        self.velocity_limits, self.sigma_limits, self.gh_limits \
                    = PPXFFit.losvd_limits(self.velscale)
        self.bias = bias
        self.degree = degree
        self.mdegree = mdegree
        self.filt_degree = filt_degree
        self.filt_mdegree = filt_mdegree
        self.moments = moments
        self.dof = self.moments + max(self.mdegree, 0)
        if self.degree >= 0:
            self.dof += self.degree+1
        self.filt_dof = self.moments + max(self.filt_mdegree, 0)
        if self.filt_degree >= 0:
            self.filt_dof += self.filt_degree+1

        #---------------------------------------------------------------
        # Report the input checks/results
        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'Number of templates: {0}'.format(self.ntpl))
            log_output(self.loggers, 1, logging.INFO, 'Number of object spectra: {0}'.format(
                                                                                        self.nobj))
            log_output(self.loggers, 1, logging.INFO, 'Pixel scale: {0} km/s'.format(self.velscale))
            log_output(self.loggers, 1, logging.INFO, 'Pixel scale ratio: {0}'.format(
                                                                            self.velscale_ratio))
            log_output(self.loggers, 1, logging.INFO, 'Dispersion limits: {0} - {1}'.format(
                                                                            *self.sigma_limits))
            log_output(self.loggers, 1, logging.INFO, 'Iteration mode: {0}'.format(
                                                                            self.iteration_mode))
            log_output(self.loggers, 2, logging.INFO, 'Model degrees of freedom: {0}'.format(
                                                                            self.dof+self.ntpl))
            if self.iteration_mode == 'fit_reject_filter':
                log_output(self.loggers, 2, logging.INFO,
                           'Model degrees of freedom for filtered spectra: {0}'.format(
                           self.filt_dof+self.ntpl))

        #---------------------------------------------------------------
        # Initialize the output data
        model_flux, model_mask, model_par = self._initialize_output(mask)

        #---------------------------------------------------------------
        # Initialize the mask and the spectral range to fit
        model_mask, err = self._initialize_pixels_to_fit(model_mask, ensemble=_ensemble,
                                                         max_velocity_range=max_velocity_range,
                                                         alias_window=alias_window)
        ended_in_error = numpy.array([e is not None for e in err])
        if numpy.any(ended_in_error):
            if not self.quiet:
                warnings.warn('Masking failures in some/all spectra.  Errors are: {0}'.format(
                                numpy.array([(i,e) for i,e in enumerate(err)])[ended_in_error]))
            model_par['MASK'][ended_in_error] \
                    = self.bitmask.turn_on(model_par['MASK'][ended_in_error], 'NO_FIT')
        if numpy.all(ended_in_error):
            return self.obj_wave, model_flux, model_mask, model_par

        # Save the pixel statistics
        model_par['BEGPIX'] = self.spectrum_start
        model_par['ENDPIX'] = self.spectrum_end
        model_par['NPIXTOT'] = self.spectrum_end - self.spectrum_start

        #---------------------------------------------------------------
        # Get the input pixel shift between the object and template
        # wavelength vectors; interpretted by pPXF as a base velocity
        # shift between the two
        self.base_velocity = numpy.array([self.ppxf_tpl_obj_voff(self.tpl_wave, self.obj_wave[s:e],
                                                                 self.velscale,
                                                                 velscale_ratio=self.velscale_ratio)
                                                for s,e in zip(self.spectrum_start,
                                                               self.spectrum_end)])

        #---------------------------------------------------------------
        # Fit the global spectrum if requested by the iteration mode
        global_fit_result = None
        if self._mode_uses_global_spectrum():
            if not self.quiet:
                log_output(self.loggers, 1, logging.INFO, 'Fitting global spectrum.')
            global_fit_result = self._fit_global_spectrum(obj_to_include=self.obj_to_fit, plot=plot)
            if global_fit_result.fit_failed():
                # pPXF failed!  
                model_mask[:] = self.bitmask.turn_on(model_mask[:], 'FIT_FAILED')
                model_par['MASK'][:] = self.bitmask.turn_on(model_par['MASK'][:], 'FIT_FAILED')
                return self.obj_wave, model_flux, model_mask, model_par

        #---------------------------------------------------------------
        # Initialize the template set according to the iteration mode
        if self._mode_uses_global_template():
            templates = numpy.dot(global_fit_result.weights, self.tpl_flux).reshape(1,-1)
            tpl_to_use = numpy.ones(1, dtype=numpy.bool)
        elif self._mode_uses_nonzero_templates():
            templates = self.tpl_flux
            tpl_to_use = numpy.array([global_weights0 > 0]*self.nobj) & self.usetpl
        elif self._mode_uses_all_templates():
            templates = self.tpl_flux
            tpl_to_use = self.usetpl
        else:
            raise ValueError('Unrecognized iteration mode: {0}'.format(self.iteration_mode))

#        #---------------------------------------------------------------
#        # Print heading of fitted kinematics for each bin
#        if not self.quiet:
#            log_output(self.loggers, 2, logging.INFO, '{0:>5}'.format('INDX')
#                       + ' {:>7} {:>7}'.format('SWAVE', 'EWAVE')
#                       + (' {:>9}'*self.moments).format(*(['KIN{0}'.format(i+1)
#                                                                for i in range(self.moments)]))
#                       + ' {0:>5} {1:>9} {2:>5} {3:>4}'.format('NZTPL', 'CHI2', 'RCHI2', 'STAT'))

        # Fit all spectra
        result = self._fit_all_spectra(templates, tpl_to_use, plot=plot,
                                       plot_file_root=plot_file_root)

        #---------------------------------------------------------------
        # Refit with emission lines
        if self._mode_includes_emission_line_fit():
            raise NotImplementedError('Iterations including emission lines have been deprecated '
                                      'for the time-being.')
#            result = self._rescale_to_fit_emission_lines(result, templates, tpl_to_use, plot=plot)

        #---------------------------------------------------------------
        # Save the results
        _,_,model_par = self._initialize_output(mask)
        model_flux, model_mask, model_par \
                = self._save_results(global_fit_result, result, model_flux, model_mask, model_par)

#        pyplot.imshow(numpy.ma.log10(self.obj_flux), origin='lower',
#                      interpolation='nearest', aspect='auto')
#        pyplot.show()
#
#        pyplot.imshow(numpy.ma.log10(model_flux), origin='lower',
#                      interpolation='nearest', aspect='auto')
#        pyplot.show()
#
#        print(numpy.sum(numpy.ma.getmaskarray(self.obj_flux)), numpy.sum(self.input_obj_mask))
#        for i in range(self.nobj):
#            pyplot.plot(self.obj_wave, self.obj_flux[i,:])
#            pyplot.plot(self.obj_wave, model_flux[i,:])
#            pyplot.show()

        if not self.quiet:
            log_output(self.loggers, 1, logging.INFO, 'pPXF finished')

        return self.obj_wave, model_flux, model_mask, model_par


    @staticmethod
    def obj_tpl_pixelmatch(velscale, tpl_wave, velscale_ratio=None, dvtol=1e-10):
        """
        Confirm that the pixel scale of the template and object data are
        the same within some tolerance, accounting for and input ratio
        of the two.
        """
        _velscale_ratio = 1 if velscale_ratio is None else velscale_ratio
#        print(numpy.absolute(velscale - spectrum_velocity_scale(tpl_wave)*_velscale_ratio), dvtol)
        return numpy.absolute(velscale - spectrum_velocity_scale(tpl_wave)*_velscale_ratio) < dvtol


    @staticmethod
    def fitting_mask(tpl_wave, obj_wave, velscale, velscale_ratio=None, waverange=None,
                     velocity_offset=None, max_velocity_range=400., alias_window=None,
                     loggers=None, quiet=False):
        """
        Return a list of pixels in the object spectrum to be fit using pPXF.
 
        Be clear between velocity (ppxf) vs. redshift (cz) !

        The limits applied to the fitted pixels are:
    
            - Apply the provided wavelength range limit (*waverange*).
            - pPXF will only allow a fit when the number of template
              pixels is the same as or exceeds the number of pixels in
              the object spectrum.  The first step toward limiting
              template spectra that are too long is to truncate the blue
              and red edges that likely won't be used given the provided
              velocity offsets (*velocity_offset*) and the expected
              velocity range (*max_velocity_range*).
            - Remove leading and trailing pixels that will cause alias
              problems during the convolution with the LOSVD
              (*alias_window*).
    
        Args:
            obj_wave (array): Wavelength vector of the object spectrum to be
                fit.
            tpl_wave (array): Wavelength vector of the template library to
                fit to the object spectrum.
            velscale (float): Velocity scale of the pixel.
            velscale_ratio (int): (**Optional**) Ratio of the object
                velscale to the template velscale.  Default is 1 (i.e.
                the two have the same pixel scale).
            waverange (array): (**Optional**) Lower and upper wavelength
                limits to *include* in the analysis.  The array can
                either define a single wavelength range -- shape is (2,)
                -- or a set of wavelength ranges -- shape is (n,2).
                Default is to apply no wavelength range limitation.
            velocity_offset (array): (**Optional**) Vector with the
                velocity offset (expected or actual) between the
                template and the object spectrum in km/s.  Used to
                estimate which wavelengths can be removed from the
                template.  This can be a single offset or a set of
                offsets.  If both waverange and velocity_offset are 2D
                arrays, the number of wavelength ranges and offsets must
                be the same.  Default is that there is no velocity
                offset.
            max_velocity_range (float): (**Optional**) Maximum range
                (+/-) expected for the fitted velocities in km/s.
                Default is 400 km/s.
            alias_window (float) : (**Optional**) The window to mask to
                avoid aliasing near the edges of the spectral range in
                km/s.  Default is six times *max_velocity_range*.
            loggers (list): (**Optional**) List of `logging.Logger`_
                objects to log progress; ignored if quiet=True.  Logging
                is done using :func:`mangadap.util.log.log_output`.
                Default is no logging.
            quiet (bool): (**Optional**) Suppress all terminal and
                logging output.  Default is False.
    
        Returns:
            numpy.ndarray: Four boolean vectors are returned:

                - flags for pixels to include in the fit
                - flags for pixels that were excluded because they were
                  outside the designated wavelength range
                - flags for pixels that were excluded to ensure the
                  proper length of the template spectrum wrt the object
                  spectrum.
                - flags for pixels that were truncated to avoid
                  convolution aliasing

        Raises:
            ValueError: Raised if (1) no pixels are valid for the fit,
                (2) the template and object do not have overlapping
                spectral regions given the expected velocity offset
                between the two, or (3) the turncation to deal with
                aliasing removes all remaining pixels.

        """
        # Check the input
        _velocity_offset = numpy.zeros(1, dtype=float) if velocity_offset is None \
                        else numpy.atleast_1d(velocity_offset)
        _alias_window = 6*max_velocity_range if alias_window is None else alias_window

        # 1. Apply the wavelength range limit, if provided
        now=len(obj_wave)                               # Number of object wavelengths
        if not quiet:
            log_output(loggers, 1, logging.INFO,
                                            'Original number of object pixels: {0}'.format(now))
        if waverange is not None and len(waverange) == 2:
            fit_indx = numpy.logical_and(obj_wave > waverange[0], obj_wave < waverange[1])
            if numpy.sum(fit_indx) == 0:
                raise ValueError('Selected wavelength range for analysis contains no pixels!')
        else:
            fit_indx = numpy.full(now, True, dtype=numpy.bool)
        waverange_mask = ~fit_indx

        # Minimum and maximum redshift about primary offsets
        c = astropy.constants.c.to('km/s').value
        z_min = (numpy.amin(_velocity_offset) - max_velocity_range)/c
        z_max = (numpy.amax(_velocity_offset) + max_velocity_range)/c
        if not quiet:
            log_output(loggers, 1, logging.INFO,
                       'Minimum/Maximum redshift: {0}/{1}'.format(z_min, z_max))

        # 2. If the number of template pixels is not >= number of fitted galaxy pixels,
        #    further limit the blue and red edges of the galaxy spectra.
        #    **This must account for the relative pixel scale.**
        now=numpy.sum(fit_indx)                 # Number of good object pixels
        _velscale_ratio = 1 if velscale_ratio is None else velscale_ratio
        ntw=len(tpl_wave)//_velscale_ratio       # Number of template pixels
        if not quiet:
            log_output(loggers, 1, logging.INFO,
                       'After selecting the wavelength range to analyze: {0}'.format(now))
            log_output(loggers, 1, logging.INFO,
                'Number of template pixels (in units of the galaxy pixel scale): {0}'.format(ntw))
        if ntw < now:
            # Indices of wavelengths redward of the redshifted template
            # TODO: Change this to the rigorous calculation of the pPXF
            # velocity: see 
            indx = obj_wave > tpl_wave[0]*(1. + z_min)

            if numpy.sum(indx) == 0:
                raise ValueError('No overlapping wavelengths between galaxy and template!')
            if not quiet:
                log_output(loggers, 1, logging.INFO,
                           'Initial wavelength of template: {0}'.format(tpl_wave[0]))
                log_output(loggers, 1, logging.INFO,
                           'Initial wavelength of redshifted template: {0}'.format(
                                                                    tpl_wave[0]*(1. + z_min)))
                log_output(loggers, 1, logging.INFO,
                           'Initial wavelength of spectrum: {0}'.format(obj_wave[0]))
                log_output(loggers, 1, logging.INFO,
                           'Pixels blueward of redshifted template: {0}'.format(
                                                                    len(obj_wave)-numpy.sum(indx)))
            # Merge with current index
            fit_indx &= indx
            if numpy.sum(fit_indx) == 0:
                raise ValueError('No overlapping wavelengths between galaxy and template!')
            now_= numpy.sum(fit_indx)
            if not quiet:
                log_output(loggers, 1, logging.INFO,
                           'After merging with the specified wavelength range: {0}'.format(now_))
        else:
            now_ = now
    
        # New number of good object pixels
        if ntw < now_:
            fit_indx[ numpy.where(fit_indx)[0][ntw:] ] = False      # Truncate the red edge as well
            if not quiet:
                log_output(loggers, 1, logging.INFO, 
                           'Limit to at most the number of template pixels: {0} {1}'.format(
                                                                        numpy.sum(fit_indx), ntw))
        npix_mask = ~(fit_indx) & ~(waverange_mask)

        # 3. Limit wavelength range to avoid aliasing problems in the template convolution
        nalias = int(numpy.floor(_alias_window/velscale))            # Number of pixels to mask
        # Mask to the range that should be unaffected by alias errors
        waverange_tpl = [ tpl_wave[nalias*_velscale_ratio]*(1+z_max),
                          tpl_wave[(ntw-nalias)*_velscale_ratio-1]*(1+z_min) ]
        if not quiet:
            log_output(loggers, 1, logging.INFO,
                       'Mask to these wavelengths to avoid convolution aliasing: {0} - {1}'.format(
                                                            waverange_tpl[0], waverange_tpl[1]))
        indx = numpy.logical_and(obj_wave > waverange_tpl[0], obj_wave < waverange_tpl[1])

        # Merge with current index
        fit_indx &= indx
        if numpy.sum(fit_indx) == 0:
            raise ValueError('No wavelengths available in this range!')

        if not quiet:
            log_output(loggers, 1, logging.INFO,
                       'Final wavelength range to fit: {0} {1}'.format(obj_wave[fit_indx][0],
                                                                       obj_wave[fit_indx][-1]))
        alias_mask = ~(fit_indx) & ~(npix_mask)

        return fit_indx, waverange_mask, npix_mask, alias_mask


    @staticmethod
    def ppxf_tpl_obj_voff(tpl_wave, obj_wave, velscale, velscale_ratio=None):
        """
        Determine the pseudo offset in velocity between the template and
        object spectra, just due to the difference in the starting
        wavelengths.
    
        This calculation is independent of the base of the logarithm used
        the sampling of the spectra.

        Assumes wavelengths are logarithmically binned.
    
        Args:
            tpl_wave (numpy.ndarray): Wavelength vector for the template
                library to fit to the object spectrum.
            obj_wave (numpy.ndarray): Wavelength vector for the object
                spectrum to be fit.
            velscale (float): Velocity step per pixel in km/s for the
                **object** spectrum.
            velscale_ratio (int): (**Optional**) The **integer** ratio
                between the velocity scale of the pixel in the galaxy
                data to that of the template data.  This is used only
                when constructing the template library.  Default is
                None, which is the same as assuming that the velocity
                scales are identical.
    
        Returns:
            float: Velocity offset in km/s between the initial wavelengths
            of the template and object spectra.

        .. todo::
            - Implement a check that calculates the velocity ratio directly?
        """
        dlogl = numpy.log(obj_wave[0])-numpy.log(tpl_wave[0]) if velscale_ratio is None \
                    else numpy.log(obj_wave[0])-numpy.mean(numpy.log(tpl_wave[0:velscale_ratio]))
        return dlogl*velscale / numpy.diff(numpy.log(obj_wave[0:2]))[0]


    @staticmethod
    def check_templates(tpl_wave, tpl_flux, tpl_sres=None):
        r"""
        Check that the input template data is valid for use with
        pPXFFit.  Static so that it can also be called by
        :func:`construct_models`.
        """
        if len(tpl_wave.shape) != 1:
            raise ValueError('Input template wavelengths must be a vector; all spectra should '
                             'have the same wavelength solution.')
        if isinstance(tpl_flux, numpy.ma.MaskedArray):
            raise TypeError('Template spectra cannot be masked arrays!')
        if tpl_flux.shape[1] != len(tpl_wave):
            raise ValueError('The template spectra fluxes must have the same length as the '
                             'wavelength array.')
        if tpl_sres is not None and tpl_sres.shape != tpl_wave.shape:
            raise ValueError('Provided template resolution vector does not have the correct shape.')
        _tpl_sres = None if tpl_sres is None \
                            else spectral_resolution(tpl_wave, tpl_sres, log10=True)
        # TODO: Allow spectral resolution to be spectrum dependent?
        return tpl_wave, tpl_flux, _tpl_sres


    @staticmethod
    def check_objects(obj_wave, obj_flux, obj_ferr=None, obj_sres=None):
        """
        Check that the input object data is valid for use with pPXFFit.
        Static so that it can also be called by
        :func:`construct_models`.
        """
        if len(obj_wave.shape) != 1:
            raise ValueError('Input object wavelengths must be a vector; all spectra should '
                             'have the same wavelength solution.')
        if obj_flux.shape[1] != len(obj_wave):
            raise ValueError('The object spectra fluxes must have the same length as the '
                             'wavelength array.')
        _obj_flux = obj_flux if isinstance(obj_flux, numpy.ma.MaskedArray) \
                                else numpy.ma.MaskedArray(obj_flux)
        if obj_ferr is not None and obj_ferr.shape != _obj_flux.shape:
            raise ValueError('The shape of any provided error array must match the flux array.')
        _obj_ferr = obj_ferr if isinstance(obj_ferr, numpy.ma.MaskedArray) \
                                else numpy.ma.MaskedArray(obj_ferr)
        if obj_sres is not None and obj_sres.shape != obj_wave.shape \
                and obj_sres.shape != _obj_flux.shape:
            raise ValueError('Provided object resolution vector does not have the correct shape.')
        _obj_sres = None if obj_sres is None else obj_sres.copy()
        if _obj_sres is not None and _obj_sres.shape != _obj_flux.shape:
            _obj_sres = numpy.array([_obj_sres]*_obj_flux.shape[0])
        return obj_wave, _obj_flux, _obj_ferr, _obj_sres


    @staticmethod
    def check_pixel_scale(tpl_wave, obj_wave, velscale_ratio=None, dvtol=1e-10):
        """
        Confirm that the pixel scale of the template and object spectra
        are identical within a certain tolerance, accounting for an
        input pixel-scale ratio.  Returns the velocity scale of the
        object spectra and the velocity scale ratio wrt the template
        spectra.

        Static so that it can be called by :func:`construct_models`.
        """
        # Get the pixel scale
        velscale = spectrum_velocity_scale(obj_wave)
        _velscale_ratio = 1 if velscale_ratio is None else velscale_ratio
        if not PPXFFit.obj_tpl_pixelmatch(velscale, tpl_wave, velscale_ratio=_velscale_ratio,
                                          dvtol=dvtol):
            raise ValueError('Pixel scale of the object and template spectra must be identical.')
        return velscale, _velscale_ratio


    @staticmethod
    def losvd_limits(velscale):
        r"""
        Return the limits on the LOSVD parameters used by pPXF.

            - Velocity limits are :math:`\pm 2000` km/s
            - Velocity-disperison limits are from 1/10 pixels to 1000 km/s
            - Limits of the higher orders moments are from -0.3 to 0.3
        """
        velocity_limits = numpy.array([-2e3, 2e3])
        sigma_limits = numpy.array([0.1*velscale, 1e3])
        gh_limits = numpy.array([-0.3, 0.3])
        return velocity_limits, sigma_limits, gh_limits


    @staticmethod
    def reject_model_outliers(obj_flux, ppxf_result, rescale=False, local_sigma=False, boxcar=None,
                              nsigma=3.0, niter=9):
        if boxcar is None and local_sigma:
            raise ValueError('For local sigma determination, must provide boxcar size.')
        model_flux = PPXFFit.compile_model_flux(obj_flux, ppxf_result, rescale=rescale)

        if local_sigma:
            print('Rejecting. Boxcar size is {0}'.format(boxcar))
            residual = obj_flux-model_flux # This should be masked where the data were not fit
            bf = BoxcarFilter(boxcar, lo=nsigma, hi=nsigma, niter=niter, y=residual,
                              local_sigma=True)
            obj_flux[bf.output_mask] = numpy.ma.masked
            return obj_flux

        for i in range(niter):
            residual = obj_flux-model_flux  # This should be masked where the data were not fit
            sigma = numpy.array([numpy.ma.std(residual, axis=1)]*residual.shape[1]).T
            indx = numpy.absolute(residual) > nsigma*sigma
            old_mask = numpy.ma.getmaskarray(obj_flux).copy()
            obj_flux[indx] = numpy.ma.masked
            if numpy.sum(numpy.ma.getmaskarray(obj_flux)) == numpy.sum(old_mask):
                break
        return obj_flux


    @staticmethod
    def compile_model_flux(obj_flux, ppxf_result, rescale=False):
        """
        Return the model flux, masked where the data were not fit.
        """
        model_flux = numpy.ma.MaskedArray(numpy.zeros(obj_flux.shape, dtype=float),
                                          mask=numpy.ones(obj_flux.shape, dtype=bool))
        for i in range(obj_flux.shape[0]):
            if ppxf_result[i] is None or ppxf_result[i].fit_failed():
                continue
            s = ppxf_result[i].start
            e = ppxf_result[i].end
            g = ppxf_result[i].gpm
            scale = optimal_scale(obj_flux[i,s:e][g], ppxf_result[i].bestfit[g]) if rescale else 1.
            model_flux[i,s:e] = scale*ppxf_result[i].bestfit
        return model_flux


    @staticmethod
    def convert_velocity(v, verr):
        r"""
        Convert kinematics from pPXF from pixel shifts to redshifts.
        pPXF determines the velocity offset by making the approximation
        that every pixel (logarithmically binned in wavelength) is a
        constant change in velocity.  Over large velocity shifts, this
        approximation can become poor.  An e-mail summary from Michele
        Cappellari:

        The velocity scale per pixel is input as
        
        .. math::
            
            \delta v = c \delta\ln\lambda = c (\ln\lambda_1 -
            \ln\lambda_0)

        The velocites output by pPXF are:

        .. math::

            V = \delta V N_{\rm shift} = c \ln(\lambda_{N_{\rm
            shift}}/\lambda_0)

        which implies that the relation between PPXF output velocity and
        redshift is

        .. math::
            
            1 + z = exp(V/c),

        which reduces z~vel/c in the low-redshift limit.  This function
        converts the :math:`V` values provided by pPXF to :math:`cz`
        velocities.

        .. note::

            **IMPORTANT**: After this conversion, you must revert the
            velocities back to the "pixel-based velocities" (using
            :func:`_revert_velocity`) before using the velocities to
            reconstruct the pPXF fitted model.
        """
        c=astropy.constants.c.to('km/s').value
        return (numpy.exp(v/c)-1.0)*c, verr*numpy.absolute(numpy.exp(v/c))


    @staticmethod
    def revert_velocity(v,verr):
        """
        Revert the velocity back to the "pixelized" velocity returned by
        pPXF.

        Order matters here.  The error computation is NOT a true error
        propagation; it's just the inverse of the "convert" operation
        """
        c=astropy.constants.c.to('km/s').value
        _v = c*numpy.log(v/c+1.0)
        return _v, verr/numpy.absolute(numpy.exp(_v/c))


#    @staticmethod
#    def shift_via_resample(wave, flux, redshift, inLog, outRange, outpix, outLog):
#        inRange = wave[[0,-1]] * (1.0 + redshift)
#        new_wave, new_flux = resample_vector(flux, xRange=inRange, inLog=inLog, newRange=outRange,
#                                             newpix=outpix, newLog=outLog, flat=False)
#        return new_flux


    @staticmethod
    def reconstruct_model(tpl_wave, templates, obj_wave, kin, weights, velscale, polyweights=None,
                          mpolyweights=None, start=None, end=None, redshift_only=False,
                          sigma_corr=0.0, velscale_ratio=None, dvtol=1e-10, revert_velocity=True):
        """
        Construct a pPXF model spectrum based on a set of input spectra and parameters.
        """
        # Make sure the pixel scales match
        _velscale_ratio = 1 if velscale_ratio is None else velscale_ratio
        if not PPXFFit.obj_tpl_pixelmatch(velscale, tpl_wave, velscale_ratio=_velscale_ratio,
                                          dvtol=dvtol):
            raise ValueError('Pixel scale of the object and template spectra must be identical.')

        # Moments for each kinematic component
        ncomp = 1
        moments = numpy.atleast_1d(kin.size)
        _moments = numpy.full(ncomp, numpy.absolute(moments), dtype=int) if moments.size == 1 \
                            else numpy.absolute(moments)

        # Get the redshift to apply
        redshift = kin[0]/astropy.constants.c.to('km/s').value

        # Check that the corrected sigma is defined
        corrected_sigma = numpy.square(kin[1]) - numpy.square(sigma_corr)
        if not redshift_only and not corrected_sigma > 0:
            warnings.warn('Corrected sigma is 0 or not defined.  Redshifting template only.')
            _redshift_only = True
        else:
            _redshift_only = redshift_only

        # Start and end pixel in the object spectrum to fit
        _start = 0 if start is None else start
        _end = obj_wave.size if end is None else end

        # Construct the composite template
        composite_template = numpy.dot(weights, templates)

#        pyplot.step(tpl_wave, composite_template, where='mid', linestyle='-', lw=0.5,
#                    color='k')
#        pyplot.show()

        # Construct the output models
        model = numpy.ma.zeros(obj_wave.size, dtype=numpy.float)
        if _redshift_only:
            # Resample the redshifted template to the wavelength grid of
            # the binned spectra
            inRange = tpl_wave[[0,-1]] * (1.0 + redshift)
            _, model[:] = resample_vector(composite_template, xRange=inRange, inLog=True,
                                                 newRange=obj_wave[[0,-1]], newpix=obj_wave.size,
                                                 newLog=True, flat=False)
            model[:_start] = 0.0
            model[:_start] = numpy.ma.masked
            model[_end:] = 0.0
            model[_end:] = numpy.ma.masked
        else:
            # Perform the same operations as pPXF v6.0.0

            # Get the offset velocity just due to the difference in the
            # initial wavelength of the template and object data
            vsyst = -PPXFFit.ppxf_tpl_obj_voff(tpl_wave, obj_wave[_start:_end], velscale,
                                               velscale_ratio=_velscale_ratio)
            # Account for a modulus in the number of object pixels in
            # the template spectra
            if _velscale_ratio > 1:
                npix_temp = composite_template.size - composite_template.size % _velscale_ratio
                _composite_template = composite_template[:npix_temp].reshape(-1,1)
                npix_temp //= _velscale_ratio
            else:
                _composite_template = composite_template
                npix_temp = composite_template.size

            # Get the FFT of the composite template
            ctmp_rfft, npad = _templates_rfft(_composite_template)

            # Construct the LOSVD parameter vector
            vj = numpy.append(0, numpy.cumsum(_moments)[:-1])
            par = kin.copy()
            # Convert the velocity to pixel units
            if revert_velocity:
                par[0], verr = PPXFFit.revert_velocity(par[0], 1.0)
            # Convert the velocity dispersion to ignore the
            # resolution difference
            par[1] = numpy.sqrt(numpy.square(par[1]) - numpy.square(sigma_corr))
           
            # Convert the kinematics to pixel units
            par[0:2] /= velscale

            # Construct the FFT of the LOSVD kernel
            kern_rfft = _losvd_rfft(par, 1, _moments, vj, npad, 1, vsyst/velscale, _velscale_ratio,
                                    0.0)

            # Construct the model spectrum
            _model = numpy.fft.irfft(ctmp_rfft[:,0] * kern_rfft[:,0,0]).ravel()
            if _velscale_ratio > 1:
                _model = numpy.mean(_model[:npix_temp*_velscale_ratio].reshape(-1,_velscale_ratio),
                                    axis=1)
            
            # Copy the model to the output vector
            model[_start:_end] = _model[:_end-_start]

        # Account for the polynomials
        x = numpy.linspace(-1, 1, _end-_start)
        if mpolyweights is not None:
            model[_start:_end] *= numpy.polynomial.legendre.legval(x,numpy.append(1.0,mpolyweights))
        if polyweights is not None:
            model[_start:_end] += numpy.polynomial.legendre.legval(x,polyweights)

        return model


    @staticmethod
    def construct_models(tpl_wave, tpl_flux, obj_wave, obj_flux, model_par, select=None,
                         velscale_ratio=None, redshift_only=False, dvtol=1e-10):
        """
        Construct models using the provided set of model parameters.
        This is a wrapper for :func:`PPXFFit.reconstruct_model`.

        Allows for a replacement template library that must have the
        same shape as :attr:`tpl_flux`.

        If redshift_only is true, ignore the LOSVD and simply shift the
        model to the correct redshift.
        """
        # Check the input spectra
        _tpl_wave, _tpl_flux, _ = PPXFFit.check_templates(tpl_wave, tpl_flux)
        ntpl = _tpl_flux.shape[0]
        _obj_wave, _obj_flux, _, _ = PPXFFit.check_objects(obj_wave, obj_flux)
        nobj = _obj_flux.shape[0]

        # Check the spectral sampling
        _velscale, _velscale_ratio \
                = PPXFFit.check_pixel_scale(_tpl_wave, _obj_wave, velscale_ratio=velscale_ratio,
                                             dvtol=dvtol)

        # Check the shape of the input model parameter database
        if model_par['BINID'].size != nobj:
            raise ValueError('Incorrect number of model-parameter sets.')
        if model_par['TPLWGT'].shape[1] != ntpl:
            raise ValueError('The number of weights does not match the number of templates.')

        # Only produce selected models
        skip = numpy.zeros(nobj, dtype=bool) if select is None else numpy.invert(select)

        # Instantiate the output model array
        models = numpy.ma.zeros(_obj_flux.shape, dtype=numpy.float)

        # Construct the model for each (selected) object spectrum
        for i in range(nobj):
            if skip[i]:
                models[i,:] = numpy.ma.masked
                continue

            models[i,:] = PPXFFit.reconstruct_model(_tpl_wave, _tpl_flux, _obj_wave,
                                                    model_par['KIN'][i,:], model_par['TPLWGT'][i],
                                                    _velscale, polyweights=model_par['ADDCOEF'][i],
                                                    mpolyweights=model_par['MULTCOEF'][i],
                                                    start=model_par['BEGPIX'][i],
                                                    end=model_par['ENDPIX'][i],
                                                    redshift_only=redshift_only,
                                                    velscale_ratio=_velscale_ratio, dvtol=dvtol)
        return models

