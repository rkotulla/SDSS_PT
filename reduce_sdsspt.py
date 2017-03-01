#!/usr/bin/env python

import os
import sys
import pyfits
import numpy
from optparse import OptionParser
import tempfile

import sdss2fits

import config
sys.path.append(config.qr_dir)

import dev_ccmatch
from podi_definitions import SXcolumn

def reduce_sdss(fn,
                overscan=True,
                trim=True,
                caldir=None,
                subtract_bias=None,
                correct_flat=None,
                fixwcs=False,

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
            # print os.path.abspath(_fn)
            filterlist[filter] = _fn
            if (os.path.isfile(_fn) and correct_flat != False and
                flat_hdus is None):
                _flat_hdus[filter] = pyfits.open(_fn)
        if (correct_flat != False and flat_hdus is None):
            flat_hdus = _flat_hdus
            correct_flat = True

    # print flat_hdus

    hdulist = sdss2fits.open_sdss_fits(fn)
    # hdulist.info()

    filtername = hdulist[0].header['FILTER']

    #
    # Apply overscan correction
    #
    linewise_overscan = False

    rawdata = hdulist['SCI'].data.astype(numpy.float)

    overscan_block1 = rawdata[:, 0:40]
    overscan_block2 = rawdata[:, 2089:2128]

    overscan_block1 = rawdata[:, 0:20] #.astype(numpy.float)
    overscan_block2 = rawdata[:, 2109:2128] #.astype(numpy.float)

    if (linewise_overscan):
        overscan1 = numpy.median(overscan_block1, axis=1).reshape((-1,1))
        overscan2 = numpy.median(overscan_block2, axis=1).reshape((-1,1))
    else:
        overscan1 = numpy.median(overscan_block1)
        overscan2 = numpy.median(overscan_block2)

    # print overscan1, overscan2

    left_read = rawdata[:, 41:1064] #.astype(numpy.float)
    left_read -= overscan1
    right_read = rawdata[:, 1064:2087] #.astype(numpy.float)
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
        data -= bias_hdu['SCI'].data

    if (flat_hdus is not None and filtername in flat_hdus):
        print "correcting flat-field"
        data /= flat_hdus[filtername]['SCI'].data

    #
    # write results
    #
    hdulist[1].data = data

    if (fixwcs):
        # write current frame to temporary file
        tmpfile = tempfile.NamedTemporaryFile(
            suffix=".fits",
            delete=True,
        )
        print tmpfile.name
        hdulist.writeto(tmpfile)

        # Run sextractor to get source catalog
        catfile, catfilename = tempfile.mkstemp(suffix=".cat")
        basedir,_ = os.path.split(os.path.abspath(__file__))
        sex_config = "%s/config/wcsfix.sex" % (config.qr_dir)
        sex_param = "%s/config/wcsfix.sexparam" % (config.qr_dir)
        sex_cmd = "sex -c %s -PARAMETERS_NAME %s -CATALOG_NAME %s %s" % (
            sex_config, sex_param, catfilename, tmpfile.name
        )
        print sex_cmd
        os.system(sex_cmd)

        # load catalog
        source_catalog = numpy.loadtxt(catfilename)
        source_catalog[:, SXcolumn['ota']] = 0
        print source_catalog.shape

        hdulist[0].header['FILTER'] = "odi_%s" % (hdulist[0].header['FILTER'])
        hdulist[1].header['OTA'] = 0
        ccmatch_results = dev_ccmatch.ccmatch(
            source_catalog=source_catalog,
            reference_catalog=None,
            input_hdu=hdulist,
            mode='otashear',
            max_pointing_error=15,
            use_ota_coord_grid=False,
        )
        print ccmatch_results




    return hdulist


if __name__ == "__main__":

    parser = OptionParser()
    parser.add_option("", "--show", dest="show",
                       action="store_true", default=False)
    parser.add_option("", "--cals", dest="caldir",
                       default=None, type=str)
    parser.add_option("", "--fixwcs", dest="fixwcs",
                       action="store_true", default=False)
    (options, cmdline_args) = parser.parse_args()


    show_list = []

    for fn in cmdline_args[1:]:
        hdulist = reduce_sdss(fn,
                              caldir=options.caldir,
                              fixwcs=options.fixwcs,)

        object = hdulist[0].header['NAME']
        filtername = hdulist[0].header['FILTER']

        out_fn = "%s_%s_%s.fits" % (fn[:-4], object, filtername) #+'.red.fits'
        print "Writing results to %s" % (out_fn)
        hdulist.writeto(out_fn, clobber=True)

        # os.system("ds9 %s &" % (out_fn))

        if (options.show):
            show_list.append(out_fn)

    if (options.show):
        os.system("ds9 %s &" % (" ".join(show_list)))
