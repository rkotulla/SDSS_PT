#!/usr/bin/env python

import os
import sys
import pyfits
import numpy


import sdss2fits


def reduce_sdss(fn,
                overscan=True,
                trim=True,
                caldir=None,
                subtract_bias=None,
                correct_flat=None,

                bias_hdu=None,
                flat_hdus=None,):

    if (caldir is not None):
        bias_fn = "%s/masterbias.fits" % (caldir)
        if (os.path.isfile(bias_fn) and
            subtract_bias != False and
            bias_hdu is None):
            bias_hdu = pyfits.open(bias_fn)

        filterlist = {}
        _flat_hdus = {}
        for filter in ['u', 'g', 'r', 'i', 'z']:
            _fn = "%s/masterflat_%s.fits" % (caldir, filter)
            print os.path.abspath(_fn)
            filterlist[filter] = _fn
            if (os.path.isfile(_fn) and correct_flat != False and
                flat_hdus is None):
                _flat_hdus[filter] = pyfits.open(_fn)
        if (correct_flat != False and flat_hdus is None):
            flat_hdus = _flat_hdus
            correct_flat = True

    print flat_hdus

    hdulist = sdss2fits.open_sdss_fits(fn)
    # hdulist.info()

    filtername = hdulist[0].header['FILTER']

    #
    # Apply overscan correction
    #
    linewise_overscan = False

    rawdata = hdulist[0].data
    overscan_block1 = rawdata[:, 0:40]
    overscan_block2 = rawdata[:, 2089:2128]

    overscan_block1 = rawdata[:, 0:20]
    overscan_block2 = rawdata[:, 2109:2128]

    if (linewise_overscan):
        overscan1 = numpy.median(overscan_block1, axis=1).reshape((-1,1))
        overscan2 = numpy.median(overscan_block2, axis=1).reshape((-1,1))
    else:
        overscan1 = numpy.median(overscan_block1)
        overscan2 = numpy.median(overscan_block2)

    # print overscan1, overscan2

    left_read = rawdata[:, 41:1064]
    left_read -= overscan1
    right_read = rawdata[:, 1064:2087]
    right_read -= overscan2

    #
    # trim off overscan region
    #
    data = rawdata[22:, 41:2087]

    #
    # subtract off bias if valid data
    #
    if (bias_hdu is not None):
        print "subtracting bias"
        data -= bias_hdu[0].data

    if (flat_hdus is not None and filtername in flat_hdus):
        print "correcting flat-field"
        data /= flat_hdus[filtername][0].data

    #
    # write results
    #
    hdulist[0].data = data

    return hdulist


if __name__ == "__main__":

    fn = sys.argv[1]

    caldir = sys.argv[2]

    hdulist = reduce_sdss(fn, caldir=caldir)

    object = hdulist[0].header['NAME']
    filtername = hdulist[0].header['FILTER']

    out_fn = "%s_%s_%s.fits" % (fn[:-4], object, filtername) #+'.red.fits'
    print "Writing results to %s" % (out_fn)
    hdulist.writeto(out_fn, clobber=True)

    os.system("ds9 %s &" % (out_fn))

