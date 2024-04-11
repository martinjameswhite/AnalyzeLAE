import numpy as np
import os
from   astropy.table import Table


def get_bb_brns(bb, photdir=None):
    assert bb in ["HSC", "LS"]
    if photdir is None:
        photdir = "/global/cfs/cdirs/desi/users/raichoor/laelbg/odin/phot/"
    if bb == "HSC":
        bands = ["g", "r", "i", "z"]
        keys = ["forced_nexp_{}".format(band) for band in bands]
        keys += ["forced_nexp_{}2".format(band) for band in ["r", "i"]]
        fn = os.path.join(photdir, "ODIN_N419_tractor_HSC_forced_all.fits.gz")
    if bb == "LS":
        bands = ["g", "r", "z"]
        keys = ["forced_nexp_{}".format(band) for band in bands]
        fn = os.path.join(photdir, "ODIN_N419_tractor_DR10_forced_all.fits.gz")
    d   = Table.read(fn, format='fits')
    sel = np.ones(len(d), dtype=bool)
    for band in bands:
        key = "forced_nexp_{}".format(band)
        if (bb == "HSC") & (band in ["r", "i"]):
            key2 = "forced_nexp_{}2".format(band)
            sel &= ((d[key] > 0) & (d[key] < 999999)) | ((d[key2] > 0) & (d[key2] < 999999))
        else:
            sel &= (d[key] > 0) & (d[key] < 999999)
    brns = np.unique(d["brickname"][sel])
    return(brns)
