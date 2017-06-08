#!/usr/bin/env python

import os
import sys
import pyfits
import numpy
from optparse import OptionParser
import tempfile
import logging
import glob
import math

import sdss2fits
import makemastercals
import reduce_sdsspt
import warnings
import config
sys.path.append(config.qr_dir)

import dev_ccmatch
from podi_definitions import SXcolumn
from podi_collectcells import apply_wcs_distortion
import podi_logging

ticker = ['1\b', '2\b', '3\b', '4\b', '5\b', '6\b', '7\b', '8\b', '9\b', '='] #"123456789="
if __name__ == "__main__":

    logsetup = {}
    podi_logging.setup_logging(logsetup)

    parser = OptionParser()
    # parser.add_option("", "--show", dest="show",
    #                    action="store_true", default=False)
    # parser.add_option("", "--cals", dest="caldir",
    #                    default=None, type=str)
    parser.add_option("", "--reusecals", dest="reuse_cals",
                       action="store_true", default=False)
    parser.add_option("", "--outdir", dest="out_dir",
                      default=None, type=str)
    parser.add_option("", "--object", dest="object",
                      default=None, type=str)
    (options, cmdline_args) = parser.parse_args()


    for night_dir in cmdline_args:

        logger = logging.getLogger(night_dir)

        #
        # Create all master calibration frames
        #
        filelist = glob.glob("%s/*.fit" % (night_dir))
        filelist.extend(glob.glob("%s/*.fit.gz" % (night_dir)))
        logger.debug("Full input file list:\n%s" % ("\n".join(filelist)))

        #
        # open all frames and decide which frames are calibrations
        # and which frames are science-frames
        #
        cal_list = []
        sci_list = []
        logger.info("Sorting %s frames into SCI and CALs" % (len(filelist)))
        n = int(math.ceil(float(len(filelist))/10.))
        sys.stdout.write(" "*n+"x\r")
        sys.stdout.flush()

        for i,in_file in enumerate(filelist):

            sys.stdout.write(ticker[i%10])
            sys.stdout.flush()

            logger.debug("Checking out %s" % (in_file))
            hdulist = sdss2fits.open_sdss_fits(in_file)
            # with warnings.catch_warnings():
            #     warnings.simplefilter("ignore")
            #     hdulist = pyfits.open(in_file)

            flavor = hdulist[0].header['FLAVOR'].lower()
            hdulist.close()

            logger.debug("%(FLAVOR)s %(EXPTIME)d %(FILTER)s" % hdulist[0].header)
            if (flavor in ['ignore']):
                continue

            elif (flavor in ['bias',
                             'dome', 'flat']):
                cal_list.append(in_file)

            else:
                object_name = hdulist[0].header['OBJECT']
                if (options.object is not None):
                    if (not object_name.startswith(options.object)):
                        logger.debug("%s, OBJECT %s --> not selected (%s)" % (
                            in_file, object_name, options.object))
                        continue
                    else:
                        logger.debug("%s, OBJECT %s --> valid object (%s)" % (
                            in_file, object_name, options.object))
                sci_list.append(in_file)
        print
        logger.info("Found %d SCIs and %d CALs" % (
            len(sci_list), len(cal_list)))

        #
        # Now create all master-calibrations (bias & flats)
        #
        cals_dir = os.path.join(night_dir, 'cals/')
        if (os.path.isdir(cals_dir) and options.reuse_cals):
            logger.info("Re-using existing calibrations")
        else:
            if (not os.path.isdir(cals_dir)):
                os.makedirs(cals_dir)
            logger.info("Creating mastercals, writing to %s" % (cals_dir))
            makemastercals.make_mastercals_from_filelist(
                filelist=cal_list,
                cals_dir=cals_dir)

        #
        # And now go and reduce each of the science frames
        #
        for sci_frame in sci_list:

            logger.info("Working on %s" % (sci_frame))

            _, bn_ext = os.path.split(os.path.abspath(sci_frame))
            bn = os.path.splitext(bn_ext)[0]
            if (options.out_dir is not None):
                out_base = os.path.join(options.out_dir,bn)
            else:
                out_base = os.path.join(night_dir,bn)

            try:
                hdulist, extras = reduce_sdsspt.reduce_sdss(
                    fn=sci_frame,
                    caldir=cals_dir,
                    fixwcs=True,
                    photcalib=True,
                    out_basename=out_base,
                )
            except:
                podi_logging.log_exception()
                continue

            object = hdulist[0].header['NAME']
            filtername = hdulist[0].header['FILTER']

            _, bn_ext = os.path.split(os.path.abspath(sci_frame))
            bn = os.path.splitext(bn_ext)[0]
            out_base_fn = "%s_%s_%s.fits" % (bn, object, filtername)
            if (options.out_dir is not None):
                out_fn = os.path.join(options.out_dir,out_base_fn)
            else:
                out_fn = os.path.join(night_dir,out_base_fn)

            logger.debug("Writing results to %s" % (out_fn))
            if (os.path.isfile(out_fn)):
                os.remove(out_fn)
            hdulist.writeto(out_fn, clobber=True)

    podi_logging.shutdown_logging(logsetup)
