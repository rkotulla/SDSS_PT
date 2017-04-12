#!/usr/bin/env python

import os
import sys
import pyfits
import numpy
from optparse import OptionParser
import tempfile
import logging
import scipy.stats

import sdss2fits

import config
sys.path.append(config.qr_dir)

import dev_ccmatch
from podi_definitions import SXcolumn
from podi_collectcells import apply_wcs_distortion
import podi_logging
import podi_sitesetup as qr_sitesetup
import podi_photcalib
import podi_commandline
import podi_fitskybackground
import podi_diagnosticplots


def reduce_sdss(fn,
                overscan=True,
                trim=True,
                caldir=None,
                subtract_bias=None,
                correct_flat=None,
                fixwcs=False,
                photcalib=False,

                bias_hdu=None,
                flat_hdus=None,
                out_basename=None,
                ):

    _,bn = os.path.split(fn)
    logger = logging.getLogger("ReduceSDSS(%s)" % (bn[:-4]))

    if (caldir is not None):
        logger.debug("Checking master-cal inventory")
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
        logger.debug("subtracting bias")
        data -= bias_hdu['SCI'].data

    if (flat_hdus is not None and filtername in flat_hdus):
        logger.debug("correcting flat-field")
        data /= flat_hdus[filtername]['SCI'].data

    #
    # write results
    #
    hdulist[1].data = data


    #
    # Import the WCS solution from the pre-canned distortion model
    #
    wcsmodel = True #False
    if (wcsmodel):
        basedir, _ = os.path.split(os.path.abspath(__file__))
        wcsmodel = "%s/wcs/wcs.fits" % (basedir)
        apply_wcs_distortion(wcsmodel, hdulist['SCI'], binning=1,
                             skip_keywords=['CD1_1', 'CD1_2', 'CD2_1', 'CD2_2'])

    filtername = hdulist[0].header['FILTER']
    exptime = hdulist[0].header['EXPTIME']
    object = hdulist[0].header['NAME']

    qr_options = podi_commandline.set_default_options()
    qr_options['otalevelplots'] = False
    ccmatch_results = None
    if (fixwcs):
        # write current frame to temporary file
        tmpfile = tempfile.NamedTemporaryFile(
            suffix=".fits",
            delete=True,
        )
        logger.debug("Writing temp file for sextractor: %s" % (tmpfile.name))
        print("Writing temp file for sextractor: %s" % (tmpfile.name))
        hdulist.info()
        hdulist.writeto(tmpfile.name)

        # Run sextractor to get source catalog
        catfile, catfilename = tempfile.mkstemp(suffix=".cat")
        sex_config = "%s/config/wcsfix.sex" % (config.qr_dir)
        sex_param = "%s/config/wcsfix.sexparam" % (config.qr_dir)
        sex_cmd = """
        sex -c %s
        -PARAMETERS_NAME %s
        -CATALOG_NAME %s
        -PHOT_APERTURES 2,3,4,5,6,8,10,12
        %s""" % (
            sex_config, sex_param, catfilename, tmpfile.name
        )
#        -PHOT_APERTURES 4,6,8,10,12,16,20,24
        logger.info("Running sextractor: %s" % (sex_cmd))
        os.system(" ".join(sex_cmd.split()))

        # load catalog
        source_catalog = numpy.loadtxt(catfilename)
        source_catalog[:, SXcolumn['ota']] = 0
        logger.info("Found %d sources" % (source_catalog.shape[0]))

        hdulist[0].header['FILTER'] = "odi_%s" % (filtername)
        hdulist[1].header['OTA'] = 0
        ccmatch_results = dev_ccmatch.ccmatch(
            source_catalog=source_catalog,
            reference_catalog=None,
            input_hdu=hdulist,
            mode='otashear',
            max_pointing_error=[5,15,30],
            use_ota_coord_grid=False,
            catalog_order=qr_sitesetup.wcscalib_order,
            #mag_limit=19.,
            #mag_limit_filter=filtername,
        )
        # print ccmatch_results

        #
        # Save the calibrated source catalog
        #
        title_info = hdulist[0].header.copy()

        wcs_matched_cat = ccmatch_results['matched_src+2mass']
        global_source_cat = ccmatch_results['calibrated_src_cat']
        wcs_quality = dev_ccmatch.global_wcs_quality(wcs_matched_cat, hdulist)

        # Create the WCS scatter plot
        plotfilename = "%s.wcs1" % (out_basename)
        podi_diagnosticplots.wcsdiag_scatter(matched_radec_odi=wcs_matched_cat[:,0:2],
                                             matched_radec_2mass=wcs_matched_cat[:,-2:],
                                             matched_ota=wcs_matched_cat[:,SXcolumn['ota']],
                                             matched_odierror=wcs_matched_cat[:, SXcolumn['mag_err_auto']],
                                             filename=plotfilename,
                                             options=qr_options,
                                             ota_wcs_stats=wcs_quality,
                                             also_plot_singleOTAs=False,
                                             title_info=title_info)

        # Create the WCS shift plot
        plotfilename = "%s.wcs2" % (out_basename)
        podi_diagnosticplots.wcsdiag_shift(matched_radec_odi=wcs_matched_cat[:,0:2],
                                           matched_radec_2mass=wcs_matched_cat[:,-2:],
                                           matched_ota=wcs_matched_cat[:,SXcolumn['ota']],
                                           filename=plotfilename, #outputfile[:-5]+".wcs2",
                                           options=qr_options,
                                           ota_wcs_stats=wcs_quality,
                                           ota_outlines=None,
                                           also_plot_singleOTAs=False,
                                           title_info=title_info)

        hdulist[0].header['FILTER'] = filtername
        del hdulist[1].header['OTA']

        numpy.savetxt(out_basename+".astrmref.cat", ccmatch_results['2mass-catalog'])
        numpy.savetxt(out_basename+".astrmref.cat.raw", ccmatch_results['astrom_ref_cat'])
        print "Used catalogs:\n%s" % "\n".join(ccmatch_results['catalog_filenames'])

    if (fixwcs and photcalib):
        photcalib_details = {}
        titlestring = "SDSS-PT: %s %s %dsec" % (object, filtername, exptime)
        zeropoint_median, zeropoint_std, odi_sdss_matched, zeropoint_exptime = \
            podi_photcalib.photcalib(
                 source_cat=ccmatch_results['calibrated_src_cat'],
                 output_filename=out_basename+".xxxx" if out_basename is not None else "testtesttest.xxxx",
                 filtername='odi_%s' % (filtername),
                 exptime=exptime,
                 diagplots=True,
                 plottitle=titlestring,
                 otalist=None,
                 options=qr_options,
                 detailed_return=photcalib_details,
                photcalib_odi_aperture=8.0 #'auto',
            )
        numpy.savetxt("sdss_matched.cat", odi_sdss_matched)
        # print photcalib_details

    #
    # Sample the background intensity
    #
    skyregions = podi_fitskybackground.sample_background(
        data=hdulist['SCI'].data,
        wcs=None, starcat=None,
    )
    skyregions = numpy.array(skyregions)
    # print skyregions
    skystats = scipy.stats.scoreatpercentile(skyregions[:, 4], [16,50,84])
    # print skystats


    logger.info("done!")

    return hdulist, ccmatch_results


if __name__ == "__main__":

    logsetup = {}
    podi_logging.setup_logging(logsetup)

    parser = OptionParser()
    parser.add_option("", "--show", dest="show",
                       action="store_true", default=False)
    parser.add_option("", "--cals", dest="caldir",
                       default=None, type=str)
    parser.add_option("", "--fixwcs", dest="fixwcs",
                       action="store_true", default=False)
    parser.add_option("", "--photcalib", dest="photcalib",
                       action="store_true", default=False)
    parser.add_option("", "--outdir", dest="outdir",
                       default=None, type=str)
    (options, cmdline_args) = parser.parse_args()


    show_list = []

    for fn in cmdline_args:

        _basedir, _bn = os.path.split(os.path.abspath(fn))
        bn = os.path.splitext(_bn)[0]
        print _basedir, _bn
        if (options.outdir is not None):
            out_fn = "%s/%s" % (options.outdir, bn)
        else:
            out_basename = "%s/%s" % (_basedir, bn)

        hdulist, wcscalib = reduce_sdss(
            fn,
            caldir=options.caldir,
            fixwcs=options.fixwcs,
            photcalib=options.photcalib,
            out_basename=out_basename,
        )

        object = hdulist[0].header['NAME']
        filtername = hdulist[0].header['FILTER']

        if (options.outdir is not None):
            _,bn = os.path.split(fn)
            out_fn = "%s/%s_%s_%s" % (
                options.outdir, bn[:-4], object, filtername) #+'.red.fits'
        else:
            out_fn = "%s_%s_%s" % (fn[:-4], object, filtername) #+'.red.fits'
        print "Writing results to %s" % (out_fn)
        if (os.path.isfile(out_fn)):
            os.remove(out_fn)
        hdulist.writeto(out_fn, clobber=True)

        # os.system("ds9 %s &" % (out_fn))

        if (options.show):
            show_list.append(out_fn)

    if (options.show):
        os.system("ds9 %s &" % (" ".join(show_list)))


    podi_logging.shutdown_logging(logsetup)
