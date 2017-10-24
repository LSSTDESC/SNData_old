import numpy as np
import pandas as pd

from astropy.io import fits
from astropy.table import Table

def read_fitres(fname):
    with open(fname) as f:
        data = f.read()
    lines = list(l for l in data.split('NVAR')[1].split('\n') if len(l.strip()) > 0)
    data = list(l.split()[1:] for l in lines[2:])
    cols = lines[1].split()[1:]
    numCols = int(lines[0].split(':')[1].strip())
    # validate
    assert numCols == len(cols)
    x = np.unique(map(len, data))
    assert len(x) == 1
    assert x[0] == numCols
    return pd.DataFrame(data, columns=cols).convert_objects(convert_numeric=True).set_index('CID')


def fname(i):
    return "RH_LSST_SNMIX_WFD_FITOPT000_SPLIT{0:03d}.FITRES.TEXT".format(i)

data = list(read_fitres(fname(i)) for i in range(1, 51))
exit()
hdulist_DDF = fits.open('LSST_Ia_HEAD.FITS')
head = Table(hdulist[1].data).to_pandas().convert_objects(convert_numeric=True)

hdulist = fits.open('LSST_Ia_PHOT.FITS')
phot = Table(hdulist[1].data).to_pandas().convert_objects(convert_numeric=True)


a = list((i, x) for (i, x) in head.loc[0:10].iterrows())


def maxsnr(df):
    arr = np.zeros(shape=(len(df), 6), dtype=np.float)
    for i, vals in df[['PTROBS_MIN', 'PTROBS_MAX']].iterrows():
        imin = vals['PTROBS_MIN'] - 1
        imax = vals['PTROBS_MAX'] - 1
        arr[i] = phot.loc[imin-1:imax-1].query('FLT == r').SNR.max()
    return arr


fig, ax = plt.subplots()
ax.hist(head.query('SNRMAX > 20.').REDSHIFT_FINAL, bins=np.arange(0., 1.4, 0.05), histtype='step')
ax.grid(True)
fig.savefig('/Users/rbiswas/Desktop/ddf.png')

