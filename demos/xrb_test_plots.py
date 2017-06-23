#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)
from itertools import product

import numpy as np
import matplotlib.pyplot as plt

from cloudyfsps import outObj as ob

dir_ = '/home/psenchyna/cloudyfsps/mist_xrb/output_mist_ssp_xrb/'
mod_prefix='ZAUX'
mods = ob.allmods(dir_, mod_prefix) #read_out=True for dust

plt.style.use('publications')

plt.figure()

for age, logZ in product(mods.age_vals, mods.logZ_vals):
    logU = -3   # fixed; ran -4, -3, -2
    w = np.where( (np.abs(mods.logU-logU)<0.1) &
                    (np.abs(mods.age-age)<10) &
                    (np.abs(mods.logZ-logZ)<0.01) )
    print(w[0].size, mods.xrbnorm_vals.size)
    print(mods.xrbnorm[w])

    logheii_hb = [np.log10(mods.mods[i].HeII / mods.mods[i].Hb)
                    for i in w[0]]
    lognii_ha = [np.log10(mods.mods[i].NIIb / mods.mods[i].Ha)
                    for i in w[0]]

    print(logheii_hb)
    print(np.std(logheii_hb))

    plt.plot(lognii_ha, logheii_hb, 
            color='blue', ls='--',
            marker='.')

# plt.xlim(-3, 1)
# plt.ylim(-3, 1)
plt.ylabel(r'$\log$ HeII 4686 / H$\beta$')
plt.xlabel(r'$\log$ NII 6584 / H$\alpha$')

plt.savefig('./heii_bpt_test.pdf', dpi=400,
        bbox_inches='tight')
