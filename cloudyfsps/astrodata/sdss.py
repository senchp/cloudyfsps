#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)
__all__ = ["plot_bpt"]
import numpy as np
import pkg_resources
from .kewley import (NII_OIII_agn_lim, NII_OIII_sf_lim,
                     OI_OIII_agn_lim, SII_OIII_agn_lim)
import matplotlib.pyplot as plt
import matplotlib.colors as mpl_colors
from matplotlib import cm as cmx

def load_spec():
    linefile = pkg_resources.resource_filename(__name__, "data/sdss_data_ls.npz")
    data = np.load(linefile)
    i, = np.where((data['lineindex_cln'] == 4) | (data['lineindex_cln'] == 5))
    outdata = dict()
    for key, val in data.items():
        outdata[key] = val
    def logify(a,b):
        if a is None or b is None:
            to_return = None
        else:
            np.seterr(all="ignore")
            to_return = np.log10(a/b)
            np.seterr(all=None)
        return to_return
    outdata['log_OIII_OII'] = logify(data['strength_OIII'][i], data['strength_OII'][i])
    outdata['log_NII_OII'] = logify(data['strength_NII'][i], data['strength_OII'][i])
    outdata['log_OIII_Hb'] = logify(data['strength_OIII'][i], data['strength_Hb'][i])
    outdata['log_OIIIb_Hb'] = logify(data['strength_OIIIb'][i], data['strength_Hb'][i])
    outdata['log_NII_Ha'] = logify(data['strength_NII'][i], data['strength_Ha'][i])
    outdata['log_NIIb_Ha'] = logify(data['strength_NIIb'][i], data['strength_Ha'][i])
    outdata['log_SII_Ha'] = logify(data['strength_SII'][i], data['strength_Ha'][i])
    outdata['log_OI_Ha'] = logify(data['strength_OI'][i], data['strength_Ha'][i])
    outdata['log_OIa_Ha'] = logify(data['strength_OIa'][i], data['strength_Ha'][i])
    outdata['log_OII_Ha'] = logify(data['strength_OII'][i], data['strength_Ha'][i])
    outdata['log_OIII_OII'] = logify(data['strength_OIII'][i], data['strength_OII'][i])
    outdata['HaHb'] = data['strength_Ha'][i]/data['strength_Hb'][i]
    outdata['R23'] = np.log10((data['strength_OII'][i] + data['strength_OIII'][i])/data['strength_Hb'][i])
    return outdata

def get_line_ratio(data, line_ratio, **kwargs):
    both_OIII = kwargs.get('both_OIII', False)
    yratio = 'log_OIIIb_Hb'
    xratio = 'log_NIIb_Ha'
    if line_ratio == 'OII': # this produces NII/OII by OIII/OII plot
        yratio = 'log_OIII_OII'
        xratio = 'log_NII_OII'
    elif line_ratio == 'R23':
        xratio = 'R23'
        yratio = 'log_OIII_OII'
    else:
        if line_ratio == 'OI':
            xratio = 'log_OIa_Ha'
            yratio = 'log_OIIIb_Hb'
        else:
            xratio = 'log_{}_Ha'.format(line_ratio)
            if (line_ratio[-1] == 'b' or line_ratio[-1] == 'a'):
                yratio = 'log_OIIIb_Hb'
            else:
                yratio = 'log_OIII_Hb'
            if both_OIII:
                yratio = 'log_OIII_Hb'
    return xratio, yratio

def plot_bpt(var_label, ax=None, color_code=False, line_ratio='NIIb', **kwargs):
    '''
    sdss.plot_bpt(True)
    SDSS data generated with astroML.fetch_corrected_sdss_spectra()
    '''
    assert line_ratio in ['NII','NIIb','SII','OI', 'OIa', 'OII', 'R23']
    if var_label:
        lab = kwargs.get('lab', 'SDSS')
    else:
        lab = '__nolegend__'
    
    data = load_spec()
    lineindex_cln = 'lineindex_cln'
    il = np.where((data['lineindex_cln'] == 4) | (data['lineindex_cln'] == 5))
    xratio, yratio = get_line_ratio(data, line_ratio, **kwargs)
    
    if ax is None:
        plt.figure()
        ax = plt.gca()
    if color_code:
        color_by = kwargs.get('color_by', 'bpt')
        if color_by == 'bpt':
            ax.scatter(data[xratio], data[yratio],
                       c=data[lineindex_cln][il], s=9, lw=0,
                       label=lab)
        elif color_by == 'HaHb':
            gi, = np.where(data[color_by] <= 15.)
            sM = retColors(data[color_by][gi], cname='gist_heat')
            for g in gi:
                if g == gi[0]:
                    plab = lab
                else:
                    plab = '__nolegend__'
                ax.plot(data[xratio][g], data[yratio][g],
                        color=sM.to_rgba(data[color_by][g]),
                        marker='.', markersize=6, label=plab)
                fig = plt.gcf()
            cb = fig.colorbar(sM)
            cb.set_label(r'$H \alpha / H\beta$')
    else:
        ax.plot(data[xratio], data[yratio], 'o',
                markersize=2.0, color='k', alpha=0.5, label=lab)
    if line_ratio[0] == 'N':
        NII_OIII_agn_lim(ax=ax)
        NII_OIII_sf_lim(ax=ax)
        ax.set_xlim(-2.0, 1.0)
        ax.set_ylim(-1.2, 1.5)
    if line_ratio[0] == 'S':
        SII_OIII_agn_lim(ax=ax)
        ax.set_xlim(-2.0, 0.3)
        ax.set_ylim(-1.2, 1.5)
    if (line_ratio == 'OI' or line_ratio == 'OIa'):
        OI_OIII_agn_lim(ax=ax)
        ax.set_xlim(-2.0, 0.0)
        ax.set_ylim(-1.2, 1.5)
    if (line_ratio == 'OII'):
        ax.set_ylim(-2.0, 1.0)
        ax.set_xlim(-1.3, 1.3)
    return

def retColors(vals, cname='CMRmap', minv=0.05, maxv=0.8, cmap=None,
               set_bad_vals=False, return_cNorm=False, logNorm=False):
    '''
    sM = get_colors(arr, cname='jet', minv=0.0, maxv=1.0)
    sM = get_colors(arr, cmap=cubehelix.cmap())
    '''
    if cmap is None:
        cmap = plt.get_cmap(cname)
    new_cmap = mpl_colors.LinearSegmentedColormap.from_list('trunc({0}, {1:.2f}, {2:.2f})'.format(cmap.name, minv, maxv), cmap(np.linspace(minv, maxv, 100)))
    if set_bad_vals:
        new_cmap.set_bad('white', alpha=1.0)
    if logNorm:
        cNorm = mpl_colors.LogNorm(vmin=vals.min(), vmax=vals.max())
    else:
        cNorm = mpl_colors.Normalize(vmin=vals.min(), vmax=vals.max())
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=new_cmap)
    if return_cNorm:
        return scalarMap, cNorm
    else:
        scalarMap.set_array(vals)
        return scalarMap
