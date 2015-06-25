# import os
import numpy as np
import astropy.units as u
import astropy.constants as c
import astropy.io as io
from astropy.table import Column

names = ['freq', 'freq_err', 'logI', 'df', 'El_cm',
         'gu', 'tag', 'qncode', 'qn', 'specname']


def load_spec(tag):
    """"""
    tb = io.ascii.read(
        'specdata/{}.txt'.format(tag),
        col_starts=[0, 14, 21, 29, 31, 41, 44, 51, 55, 81],
        col_ends=[13, 20, 28, 30, 40, 43, 50, 54, 80, 100],
        format='fixed_width_no_header',
        names=names)
    Eu_cm = Column(name='Eu_cm',
                   data=tb['El_cm']*(1/u.cm)+(tb['freq']*u.MHz/c.c).to(1/u.cm))
    tb.add_column(Eu_cm)
    return tb


def compute_qpart(table, T):
    """"""
    q = 1
    for row in table:
        q += row['gu']*np.exp(-row['Eu_cm']*(1/u.cm)*c.h*c.c/(c.k_B*T*u.K))
        # print row['gu'], np.exp(-row['Eu_cm']*(1/u.cm)*c.h*c.c/(c.k_B*T*u.K))
    return q



if __name__ == '__main__':
    tag = 28503
    T = 9.375
    tb = load_spec(tag)
    tb.pprint()
    qval = compute_qpart(tb,T)
    compute_qpart_approx(tag,T)
    print qval
