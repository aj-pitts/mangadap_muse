
import numpy as np
import numpy.ma as ma
import os
from scipy import stats
from scipy.optimize import curve_fit
from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

import logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)

from mangadap.config import defaults

plt.rcParams['figure.figsize'] = (10,6)
plt.rc('axes', labelsize = 25)
plt.rc('axes', titlesize = 18)
#plt.rc('axes', titleweight = 'bold')
plt.rc('axes', lw = 2)
plt.rc('xtick',labelsize=20)
plt.rc('ytick',labelsize=20)
plt.rcParams['font.family'] = 'stixgeneral'
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['xtick.major.size'] = 10
plt.rcParams['ytick.major.size'] = 10
plt.rcParams['xtick.minor.size'] = 8
plt.rcParams['ytick.minor.size'] = 8


class CubeData:

    def __init__(self, galname=None, bin_key=None, plate=None, ifu=None):

        # mangadap_muse root directory path
        mangadap_muse_dir = os.path.dirname(os.path.dirname(defaults.dap_data_root()))
        # output directory path
        output_root_dir = os.path.join(mangadap_muse_dir, 'outputs')

        self.galname = galname
        self.bin_key = bin_key

        output_gal_dir = os.path.join(output_root_dir, f"{self.galname}-{self.bin_key}")
        if not os.path.isdir(output_gal_dir):
            raise ValueError(f'{output_gal_dir} is not a directory within {output_root_dir}.')

        # output beta plot directory
        output_beta_dir = os.path.join(output_gal_dir, 'beta_plots')
        if not os.path.isdir(output_beta_dir):
            os.makedirs(output_beta_dir)
        self.output_beta_dir = output_beta_dir

        # non-corrected MUSE cube directory
        output_gal_sub_dir = os.path.join(output_gal_dir, 'NO-CORR')
        # key methdos from analysis plan
        analysisplan_methods = 'MILESHC-MASTARHC2-NOISM'
        # cube directory
        cube_dir = os.path.join(output_gal_sub_dir, f"{bin_key}-{analysisplan_methods}", str(plate), str(ifu))
        # paths to the LOGCUBE and MAPS files
        cube_file_path = os.path.join(cube_dir,
                                      f"manga-{plate}-{ifu}-LOGCUBE-{bin_key}-{analysisplan_methods}.fits")
        maps_file_path = os.path.join(cube_dir,
                                      f"manga-{plate}-{ifu}-MAPS-{bin_key}-{analysisplan_methods}.fits")
        cube = fits.open(cube_file_path)
        maps = fits.open(maps_file_path)
        # bin ID has multiple layers of the same bin id map so use first one
        self.binid_map = maps['BINID'].data[0]
        # obtain wavelength, flux, model and error array from MUSE cube
        self.wave = cube['WAVE'].data
        self.flux = cube['FLUX'].data
        self.model = cube['MODEL'].data
        self.error = np.sqrt(1 / cube['IVAR'].data)

    def SN_dict(self):
        # create a dictionary containing that stores the number of spaxels in each bin ID
        # and beta values for each S/N bin
        SN_dict = {}
        SN_channel_keys = ['SN_range', 'bin_size', 'beta_bin']

        SN_lims = [[0, 50], [50, 75], [75, 100], [100, 1000]]

        for i in range(len(SN_lims)):
            channel = f'SN_CHANNEL{i}'
            SN_dict[channel] = {}

            for key in SN_channel_keys:
                if key == 'SN_range':
                    SN_dict[channel][key] = SN_lims[i]
                else:
                    SN_dict[channel][key] = None

        return SN_dict

    def N_spx_dict(self):
        # create a dictionary containing the different N_spx bins
        N_spx_dict = {}
        N_spx_channel_keys = ['N_spx_range', 'SN_dict']

        N_spx_lims = [[0, 5.5], [5.5, 22], [22, 52], [52, 122]]

        for i in range(len(N_spx_lims)):
            channel = f'N_spx_CHANNEL{i}'
            N_spx_dict[channel] = {}

            for key in N_spx_channel_keys:
                if key == 'N_spx_range':
                    N_spx_dict[channel][key] = N_spx_lims[i]
                # add S/N dictionary
                if key == 'SN_dict':
                    N_spx_dict[channel][key] = self.SN_dict()

        return N_spx_dict

    def beta_all_dict(self):
        # create a dictionary containing the different N_spx bins
        beta_all_dict = {}
        beta_all_channel_keys = ['SN_range', 'bin_size_all', 'beta_all', 'bin_modes', 'beta_medians']

        SN_lims = [[0, 50], [50, 75], [75, 100], [100, 1000]]

        for i in range(len(SN_lims)):
            channel = f'SN_CHANNEL{i}'
            beta_all_dict[channel] = {}

            for key in beta_all_channel_keys:
                if key == 'SN_range':
                    beta_all_dict[channel][key] = SN_lims[i]

                elif (key == 'bin_size_all') or (key == 'beta_all'):
                    beta_all_dict[channel][key] = []

                else:
                    beta_all_dict[channel][key] = []

        return beta_all_dict

    def wv_dict(self):
        # create a dictionary for each wavelength channel
        wv_dict = {}
        wv_channel_keys = ['wv_range', 'bin_masks', 'N_spx', 'SN', 'beta', 'N_spx_dict', 'beta_all_dict']

        # wavelength ranges
        wv_lims = [[4751.42, 5212], [5212, 5672], [5672, 6132], [6132, 6592],
                   [6592, 7052], [7052, 7513], [7513, 7973], [7973, 8433], [8433, 8893], [8893, 9353.44]]

        for i in range(len(wv_lims)):
            channel = f'WAVE_CHANNEL{i}'
            wv_dict[channel] = {}

            for key in wv_channel_keys:
                if key == 'wv_range':
                    wv_dict[channel][key] = wv_lims[i]
                # add N_spx dictionary
                elif key == 'N_spx_dict':
                    wv_dict[channel][key] = self.N_spx_dict()
                # add Beta_all dictionary
                elif key == 'beta_all_dict':
                    wv_dict[channel][key] = self.beta_all_dict()
                else:
                    wv_dict[channel][key] = None

        return wv_dict

    def filter_data(self, array, conditional_array, lim):
        # filter out the input array based on the conditional array and input limits
        return array[(conditional_array >= lim[0]) & (conditional_array <= lim[1])]

    def apply_mask(self, wave_bin, flux_bin, model_bin, error_bin):
        # mask for large unaccounted emission or absorption lines if present

        # calculate the standard deviation in the fit residuals (rN)
        # https://quantifyinghealth.com/residual-standard-deviation-error/ for reference
        rN_1sigma = np.sqrt(np.mean((flux_bin - model_bin) ** 2))
        # calculate rN per pixel
        rN = np.sqrt((flux_bin - model_bin) ** 2)

        # perform a rudimentary sigma clip and
        # mask indices where the residual noise exceeds 10 sigma
        wave_mask = ma.masked_where(rN / rN_1sigma >= 10, wave_bin)
        flux_mask = ma.masked_where(rN / rN_1sigma >= 10, flux_bin)
        model_mask = ma.masked_where(rN / rN_1sigma >= 10, model_bin)
        error_mask = ma.masked_where(rN / rN_1sigma >= 10, error_bin)

        return wave_mask, flux_mask, model_mask, error_mask

    def calc_beta(self, wave, flux, model, error, binid_map):
        # calculate number of spaxels, S/N and correlation ratio value for each bin ID

        # total number of bins
        n_bins = np.max(binid_map) + 1

        bin_masks = []
        bin_sizes = np.zeros(n_bins)
        sn = np.zeros(n_bins)
        betas = np.zeros(n_bins)

        for i in range(n_bins):
            id_mask = binid_map == i
            # get bin sizes, S/N and beta value for each bin ID
            n_spx = len(flux[:, id_mask][0, :])
            # get flux, stellar continuum model and variance array from bin ID
            # print(flux)
            flux_bin = flux[:, id_mask][:, 0]
            model_bin = model[:, id_mask][:, 0]
            error_bin = error[:, id_mask][:, 0]

            # mask values
            wave_mask, flux_mask, model_mask, error_mask = self.apply_mask(wave, flux_bin, model_bin,
                                                                           error_bin)

            # calculate the S/N and correlation ratio with the masked data
            snr = np.median(flux_mask / error_mask)
            # calculate the standard deviation in the fit residuals (rN)
            rN_1sigma = np.sqrt(np.mean((flux_mask - model_mask) ** 2))
            # calculate the correlation ratio (between residual noise
            # and propagated noise (beta = rN/N)
            beta = (rN_1sigma / np.median(error_mask))

            # don't store correlation ratio if its non-finite and less than 0
            if (np.isfinite(beta) == False) or (beta < 0):
                continue
            else:
                bin_masks.append(flux_mask.mask)
                bin_sizes[i] = n_spx
                sn[i] = snr
                betas[i] = beta

        return bin_masks, bin_sizes, sn, betas

    def get_beta(self):
        # dictionary for each wavelength channel
        self.wv_dict = self.wv_dict()

        for wv_key in self.wv_dict.keys():
            # filter data
            # wave_filt, flux_filt, model_filt, err_filt = self.filter_data(self.wv_dict[key]['wv_range'])
            wave_filt = self.filter_data(self.wave, self.wave, self.wv_dict[wv_key]['wv_range'])
            flux_filt = self.filter_data(self.flux, self.wave, self.wv_dict[wv_key]['wv_range'])
            model_filt = self.filter_data(self.model, self.wave, self.wv_dict[wv_key]['wv_range'])
            err_filt = self.filter_data(self.error, self.wave, self.wv_dict[wv_key]['wv_range'])

            # get number of spaxels in each bin, S/N and beta value for each bin ID
            bin_masks, bin_sizes, sn, betas = self.calc_beta(wave_filt, flux_filt,
                                                             model_filt, err_filt, self.binid_map)
            # save the mask from each bin ID
            self.wv_dict[wv_key]['bin_masks'] = bin_masks

            # loop through and place each spaxel bin size and beta value into their respective
            # N_spx and S/N bin
            N_spx_dict = self.wv_dict[wv_key]['N_spx_dict']
            beta_all_dict = self.wv_dict[wv_key]['beta_all_dict']

            for N_spx_key in N_spx_dict.keys():
                # filter bin size, S/N and beta values according to number of spaxel (N_spx) range
                bin_sizes_N_spx_filt = self.filter_data(bin_sizes,
                                                        bin_sizes, N_spx_dict[N_spx_key]['N_spx_range'])
                sn_N_spx_filt = self.filter_data(sn,
                                                 bin_sizes, N_spx_dict[N_spx_key]['N_spx_range'])
                beta_N_spx_filt = self.filter_data(betas,
                                                   bin_sizes, N_spx_dict[N_spx_key]['N_spx_range'])

                SN_dict = N_spx_dict[N_spx_key]['SN_dict']
                for SN_key in SN_dict.keys():
                    # filter bin size and beta values according to S/N range
                    bin_size_SN_filt = self.filter_data(bin_sizes_N_spx_filt,
                                                        sn_N_spx_filt, SN_dict[SN_key]['SN_range'])
                    beta_SN_filt = self.filter_data(beta_N_spx_filt,
                                                    sn_N_spx_filt, SN_dict[SN_key]['SN_range'])
                    # save the data
                    SN_dict[SN_key]['bin_size'] = bin_size_SN_filt
                    SN_dict[SN_key]['beta_bin'] = beta_SN_filt

                N_spx_dict[N_spx_key]['SN_dict'] = SN_dict

            for all_key in beta_all_dict.keys():
                for N_spx_key in N_spx_dict.keys():
                    bin_size = N_spx_dict[N_spx_key]['SN_dict'][all_key]['bin_size']
                    beta_bin = N_spx_dict[N_spx_key]['SN_dict'][all_key]['beta_bin']

                    beta_all_dict[all_key]['bin_size_all'].append(bin_size)
                    beta_all_dict[all_key]['beta_all'].append(beta_bin)

            for all_key in beta_all_dict.keys():
                bin_size_all = beta_all_dict[all_key]['bin_size_all']
                beta_all = beta_all_dict[all_key]['beta_all']

                for bin_array, beta_array in zip(bin_size_all, beta_all):
                    if len(bin_array) == 0:
                        beta_all_dict[all_key]['bin_modes'].append(-999)
                        beta_all_dict[all_key]['beta_medians'].append(-999)

                    else:
                        beta_all_dict[all_key]['bin_modes'].append(stats.mode(bin_array)[0])
                        beta_all_dict[all_key]['beta_medians'].append(np.median(beta_array))

            self.wv_dict[wv_key]['N_spx_dict'] = N_spx_dict
            self.wv_dict[wv_key]['beta_all_dict'] = beta_all_dict

    # Sarzi+2018 quadratic function
    def beta_func_quad(self, N_spx, a, b):
        return 1 + a * (np.log10(N_spx)) ** b

    def create_beta_tables(self):

        # dictionary for each wavelength channel
        wv_dict = self.wv_dict

        for wv_key in wv_dict.keys():
            # dictionary for each wavelength channel
            wv_range = wv_dict[wv_key]['wv_range']
            beta_all_dict = wv_dict[wv_key]['beta_all_dict']

            # list to store median beta values to write into .dat file
            beta_table_data = []

            # list to store for values for fitting Sarzi+2018 relation
            bin_fit = []
            beta_fit = []

            # loop through each S/N channel within Beta_all dicitonary
            for all_key in beta_all_dict.keys():
                bin_modes = beta_all_dict[all_key]['bin_modes']
                beta_medians = beta_all_dict[all_key]['beta_medians']

                # don't keep last beta median value because the correlation
                # correction does not do well for N_spx > 100
                beta_table_data.append(beta_medians[:3])

                for bin_mode, beta_median in zip(bin_modes, beta_medians):

                    if bin_mode == -999:
                        continue
                    else:
                        bin_fit.append(bin_mode)
                        beta_fit.append(beta_median)

            beta_func_quad = self.beta_func_quad
            # get parameters
            # fit data points usings Sarzi+2018 relationshiop
            popt_sarzi, pcov_sarzi = curve_fit(beta_func_quad, np.array(bin_fit),np.array(beta_fit))

            # create an Astropy data table containting each median beta histogram distribution
            # for a given S/N and N_spx size
            beta_table = Table([beta_table_data[0], beta_table_data[1], beta_table_data[2],
                                [popt_sarzi[0], -999, -999], [popt_sarzi[1], -999, -999]],
                               names=('S_N_0-50', 'S_N_50-75', 'S_N_75-100', 'param_fit_a', 'param_fit_b'))

            self.wv_dict[wv_key]['beta_table'] = beta_table

            # write the beta table into a .dat file
            beta_tables_dir = os.path.join(defaults.dap_data_root(), 'beta_tables')
            beta_tables_gal_dir = os.path.join(beta_tables_dir, f'{self.galname}-{self.bin_key}')
            if not os.path.isdir(beta_tables_gal_dir):
                os.makedirs(beta_tables_gal_dir)

            # create the file name
            file_name = f'beta_corr_{wv_range[0]}_{wv_range[1]}.dat'
            beta_table_path = os.path.join(beta_tables_gal_dir, file_name)
            beta_table.write(beta_table_path, format='ascii', overwrite=True)
            logging.info(f'Writing {file_name} to {beta_tables_gal_dir}')

    def mk_residual_plot(self, wv_channel, bin_id, zoom=None):
        # plot the flux, model and error in main plot
        # lower plot showcases the ratio between the residual-noise-per-pixel and
        # 1 standard deviation of the residual noise over the entire wavelength channel (rN/1sigma_rN)

        # wavelength channel dictionary
        wv_dict = self.wv_dict[f'WAVE_CHANNEL{wv_channel}']
        # mask array
        mask_bin = wv_dict['bin_masks'][bin_id]

        wave_bin = self.wave
        # get flux, stellar continuum model and variance array from bin ID
        id_mask = self.binid_map == bin_id
        flux_bin = self.flux[:, id_mask][:, 0]
        model_bin = self.model[:, id_mask][:, 0]
        error_bin = self.error[:, id_mask][:, 0]

        # filter out data and apply mask from bin ID
        wave = ma.masked_array(self.filter_data(wave_bin, wave_bin, wv_dict['wv_range']), mask=mask_bin)
        flux = ma.masked_array(self.filter_data(flux_bin, wave_bin, wv_dict['wv_range']), mask=mask_bin)
        model = ma.masked_array(self.filter_data(model_bin, wave_bin, wv_dict['wv_range']), mask=mask_bin)
        err = ma.masked_array(self.filter_data(error_bin, wave_bin, wv_dict['wv_range']), mask=mask_bin)

        # 1 sigma residual noise before and afer being masked
        rN_1sigma_unmasked = np.sqrt(np.mean((flux.data - model.data) ** 2))
        rN_1sigma_masked = np.sqrt(np.mean((flux - model) ** 2))

        # show good and bad beta values per pixel before being masked
        beta_good_unmasked = rN_1sigma_unmasked / err
        beta_bad_unmasked = rN_1sigma_unmasked / err.data

        # show beta values after being masked
        beta_good_masked = rN_1sigma_masked / err.data

        # make plots
        fig = plt.figure(figsize=(10, 6))

        # create an extra masked plot if there are masked pixels
        if flux.mask.any() == True:
            gs = fig.add_gridspec(3, hspace=0, height_ratios=(3, 1, 1))
            axs = gs.subplots(sharex=True)
        else:
            gs = fig.add_gridspec(2, hspace=0, height_ratios=(3, 1))
            axs = gs.subplots(sharex=True)

        axs[0].set_title(f'Wavelength Channel {wv_channel} ({wv_dict["wv_range"][0]}-'
                         f'{wv_dict["wv_range"][1]} $\mathrm{{\AA}}$)', fontsize=22)

        axs[0].step(wave.data, flux.data, c='k', zorder=2, label='Data')
        axs[0].step(wave.data, model.data, c='tab:purple', zorder=3, label='Model')
        axs[0].step(wave.data, err.data, c='tab:green', zorder=1, label='Error')
        axs[0].set_ylabel(r'Flux ($10^{-20} \mathrm{erg}$ $\mathrm{cm}^{-2}$ $\mathrm{s}^{-1}$)', fontsize=19)
        axs[0].set_ylim(-0.01, np.median(flux) * 1.25)

        if flux.mask.any() == True:

            # manually define a patch for the masked region
            ax0_handles, labels = axs[0].get_legend_handles_labels()
            mask_patch = mpatches.Patch(color='cadetblue', label='Masked Region', alpha=0.3)
            # append manual patch
            ax0_handles.append(mask_patch)
            axs[0].legend(handles=ax0_handles, fontsize='large', loc='center left')

            # plot beta values before being maksed
            axs[1].step(wave, beta_good_unmasked,
                        label=r'$\beta$ Before Masking', c='grey', zorder=2)
            # plot bad beta values
            axs[1].step(wave.data, beta_bad_unmasked,
                        label=r'Masked Pixels', c='r', zorder=1)

            axs[1].set_ylabel(r'$\beta$', fontsize=20)
            axs[1].set_ylim(-0.5, np.median(beta_good_unmasked) * 2.2)
            axs[1].legend(fontsize='large', frameon=False, loc='upper left', ncol=2)

            # plot beta values after being masked
            axs[2].step(wave.data, beta_good_masked,
                        label=r'$\beta$ After Masking', c='grey', zorder=1)

            axs[2].set_xlabel(r'Wavelength ($\mathrm{\AA}$)', fontsize=20)
            axs[2].set_ylabel(r'$\beta$', fontsize=20)
            axs[2].set_ylim(-0.5, np.median(beta_good_unmasked) * 2.2)
            axs[2].legend(fontsize='large', frameon=False, loc='upper left')

            # highlight masked regions
            pix_size = wave.data[1] - wave.data[0]
            for masked_pix in wave.data[wave.mask]:
                axs[0].axvspan(masked_pix - (pix_size / 1), masked_pix, alpha=0.3,
                               color='cadetblue', zorder=4)
                axs[1].axvspan(masked_pix - (pix_size / 1), masked_pix, alpha=0.3,
                               color='cadetblue', zorder=4)
                axs[2].axvspan(masked_pix - (pix_size / 1), masked_pix, alpha=0.3,
                               color='cadetblue', zorder=4)

            if zoom is not None:
                axs[0].set_xlim(zoom[0], zoom[1])
                axs[1].set_xlim(zoom[0], zoom[1])
                axs[2].set_xlim(zoom[0], zoom[1])

        else:
            # only create one beta plot if theres no masking
            logging.info('No masked pixels!')

            # create legend without masked region patch
            axs[0].legend(fontsize='large', loc='center left')

            # plot beta values before being maksed
            axs[1].step(wave, beta_good_masked, c='grey', zorder=1)

            axs[1].set_xlabel(r'Wavelength ($\mathrm{\AA}$)', fontsize=20)
            axs[1].set_ylabel(r'$\beta$', fontsize=20)
            axs[1].set_ylim(-0.5, np.median(beta_good_masked) * 2.2)

            # zoom must be a list with two wavelegnth values in it
            if zoom is not None:
                axs[0].set_xlim(zoom[0], zoom[1])
                axs[1].set_xlim(zoom[0], zoom[1])

    def mk_hist_plot(self, wv_channel, show_plot=True, save_plot=False):
        # wavelength channel dictionary
        wv_dict = self.wv_dict[f'WAVE_CHANNEL{wv_channel}']

        beta_bins = np.arange(0, 15, 0.2)
        median_colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']

        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 10))
        fig.tight_layout(pad=4.0)

        N_spx_dict = wv_dict['N_spx_dict']

        for N_spx_key, ax in zip(N_spx_dict.keys(), fig.get_axes()):
            N_spx_range = N_spx_dict[N_spx_key]['N_spx_range']
            SN_dict = N_spx_dict[N_spx_key]['SN_dict']

            for i, SN_key in enumerate(SN_dict.keys()):
                SN_range = SN_dict[SN_key]['SN_range']
                beta_SN_bin = SN_dict[SN_key]['beta_bin']

                if SN_key == 'SN_CHANNEL3':
                    label = f'{SN_range[0]} < SN'
                else:
                    label = f'{SN_range[0]} < SN < {SN_range[1]}'

                ax.hist(beta_SN_bin, bins=beta_bins, label=label,
                        histtype='step', stacked=False, alpha=0.9, fill=False)
                ax.tick_params(axis='both', which='both', direction='in', top=True, right=True, length=5)
                ax.set_xticks(np.arange(1, 10))

                if np.isnan(np.median(beta_SN_bin)) == True:
                    continue
                else:
                    ax.axvline(np.median(beta_SN_bin), ls='--', c=median_colors[i], lw=1)

            if N_spx_key == 'N_spx_CHANNEL0':
                ax.set_title(f'$N_{{ \mathrm{{spax}} }}$ < {N_spx_range[1]}', fontsize=22)
                ax.legend(fontsize=16)

            elif N_spx_key == 'N_spx_CHANNEL3':
                ax.set_title(f'$N_{{ \mathrm{{spax}} }}$ > {N_spx_range[1]}', fontsize=22)

            else:
                ax.set_title(f'{N_spx_range[0]} < $N_{{\mathrm{{spax}}}}$ < {N_spx_range[1]}', fontsize=22)

            ax.set_xlim(0, 10)
            ax.set_xlabel(r'$\beta$ (rN/N)', fontsize=22)
            ax.set_ylabel('Counts')

        fig.suptitle(f' $\\beta$ Histograms for Wavelength Channel {wv_channel} ({wv_dict["wv_range"][0]}-'
                     f'{wv_dict["wv_range"][1]} $\mathrm{{\AA}}$)', fontsize=25, y=1.03)

        # save the figure if save_plot is True
        if save_plot:
            # create historgram directory if it doesn't exist already
            beta_hist_dir = os.path.join(self.output_beta_dir, 'beta_histogram_plots')
            if not os.path.isdir(beta_hist_dir):
                os.makedirs(beta_hist_dir)

            file_name = f'beta_hist_{wv_dict["wv_range"][0]}_{wv_dict["wv_range"][1]}.png'
            # save the figure image
            fig.savefig(f'{beta_hist_dir}/{file_name}', bbox_inches='tight')
            logging.info(f'Saving {file_name} to {beta_hist_dir}')

        # don't show the plot if show_plot is False
        if not show_plot:
            plt.close()
            logging.info('Not displaying plot!')

    def mk_violin_plot(self, wv_channel, show_plot=True, save_plot=False):
        # create violin plots for each wavelength channel

        # dictionary containing all the bin sizes and beta values for each
        # wavelength channel
        beta_all_dict = self.wv_dict[f'WAVE_CHANNEL{wv_channel}']['beta_all_dict']

        # slightly offset each violin plot within the same S/N channel
        x_offset = [-3, -1, 1, 3]
        colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']
        # get all the bin modes for each S/N channel
        bin_mode_all = []

        plt.figure(figsize=(10, 6))
        # legend patches
        color_patches = []
        color_patches.append(plt.Line2D([], [], color='k', ls="--", label='Sarzi+2018'))

        # loop through each S/N channel and create a violin plot
        # for each non-zero array using the beta distribution and bin size mode
        for i, SN_key in enumerate(beta_all_dict.keys()):
            SN_range = beta_all_dict[SN_key]['SN_range']
            bin_modes = np.array(beta_all_dict[SN_key]['bin_modes'])
            beta_all = beta_all_dict[SN_key]['beta_all']

            # check to see if there are any beta distributions for the given
            # S/N range
            if np.all(bin_modes == -999):
                logging.info(f'No beta distributions for S/N range between {SN_range[0]}-{SN_range[1]}')
                continue
            else:
                for bin_mode, betas in zip(bin_modes, beta_all):
                    # if array empty, skip
                    if bin_mode == -999:
                        continue
                    else:
                        bin_mode_all.append(bin_mode)
                        violinplots = plt.violinplot(betas, positions=[bin_mode + x_offset[i]],
                                                     showmedians=True, widths=7)

                        # configure each individiual violin plot colors
                        for partname in ('cbars', 'cmins', 'cmaxes', 'cmedians', 'bodies'):
                            violin = violinplots[partname]
                            if partname == 'bodies':
                                for vp in violin:
                                    vp.set_color(colors[i])
                                    vp.set_facecolor(colors[i])
                                    vp.set_edgecolor(colors[i])
                                    vp.set_alpha(0.2)
                                    vp.set_label('uhh')

                            else:
                                violin.set_edgecolor(colors[i])
                                violin.set_linewidth(2)

                        if SN_key == 'SN_CHANNEL3':
                            SN_str = f'{SN_range[0]} < SN'
                        else:
                            SN_str = f'{SN_range[0]} < SN < {SN_range[1]}'

                        # add to legend handle
                        if bin_modes[bin_modes > -999][-1] == bin_mode:
                            c_patch = mpatches.Patch(color=colors[i], label=SN_str)
                            color_patches.append(c_patch)

        # Sarzi+2018 fit parameters
        a = self.wv_dict[f'WAVE_CHANNEL{wv_channel}']['beta_table']['param_fit_a'][0]
        b = self.wv_dict[f'WAVE_CHANNEL{wv_channel}']['beta_table']['param_fit_b'][0]

        # show Sarzi relation against violin plots
        N_spx_range = np.arange(0, bin_mode_all[-1] + x_offset[-1] + 1)
        plt.plot(N_spx_range, self.beta_func_quad(N_spx_range, a, b), c='k', ls='--', label='Sarzi+2018')

        plt.tick_params(axis='both', which='both', direction='in', top=True, right=True, length=7)
        plt.ylabel(r'$\beta$ (rN/N)')
        plt.xlabel('Mode $N_{spx}$')
        # {x:.2f}
        plt.legend(handles=color_patches, fontsize='xx-large', loc='upper center', frameon=False)

        plt.figtext(x=0.15, y=0.50, s='Parameter Fits\n', fontsize=16)
        text_str = f'a = {a:.2f}, b = {b:.2f}'
        plt.figtext(x=0.15, y=0.49, s=text_str, fontsize=16)

        wv_range = self.wv_dict[f'WAVE_CHANNEL{wv_channel}']["wv_range"]
        plt.title(f'Wavelength Channel {wv_channel} ({wv_range[0]}-'
                  f'{wv_range[1]} $\mathrm{{\AA}}$)', fontsize=22)

        # save the figure if save_plot is True
        if save_plot:
            # create historgram directory if it doesn't exist already
            beta_violin_dir = os.path.join(self.output_beta_dir, 'beta_violin_plots')
            if not os.path.isdir(beta_violin_dir):
                os.makedirs(beta_violin_dir)

            file_name = f'beta_violin_{wv_range[0]}_{wv_range[1]}.png'
            # save the figure image
            plt.savefig(f'{beta_violin_dir}/{file_name}', bbox_inches='tight')
            logging.info(f'Saving {file_name} to {beta_violin_dir}')

        # don't show the plot if show_plot is False
        if not show_plot:
            plt.close()

    def mk_hist_plots(self, show_plots=False, save_plots=True):

        # create histogram plots for all the wavelength channels
        for i in range(len(self.wv_dict.keys())):
            self.mk_hist_plot(wv_channel=i, show_plot=show_plots, save_plot=save_plots)

    def mk_violin_plots(self, show_plots=False, save_plots=True):

        # create violin plots for all the wavelength channels
        for i in range(len(self.wv_dict.keys())):
            self.mk_violin_plot(wv_channel=i, show_plot=show_plots, save_plot=save_plots)