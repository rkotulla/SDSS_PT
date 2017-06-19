#!/usr/bin/env python


import os
import sys
import pyfits
import numpy
import logging

import glob
import warnings

import sdss2fits
import reduce_sdsspt

import config
sys.path.append(config.qr_dir)

import podi_imcombine

def make_master_bias(filelist):

    datablocks = []
    for fn in filelist:

        valid_hdu, extras = reduce_sdsspt.reduce_sdss(
            fn,
            overscan=True,
            trim=True,
            subtract_bias=False,
            correct_flat=False,
        )
        if (valid_hdu is None):
            continue

        # valid_hdu.info()
        datablocks.append(valid_hdu['SCI'].data)

    # print datablocks

    masterbias = podi_imcombine.imcombine_data(
        datas=datablocks,
        operation='sigmaclipmean'
    )

    masterbias_hdu = pyfits.HDUList([
        pyfits.PrimaryHDU(),
        pyfits.ImageHDU(
            data=masterbias,
            name='SCI',
    )])

    return masterbias_hdu



def make_master_flat(filelist, bias_hdu, write_norm_flat=False,
                     good_flux=None):

    if (good_flux is None):
        good_flux=[10000,45000]

    logger = logging.getLogger("MakeMasterFlat")
    datablocks = []
    for fn in filelist:

        valid_hdu, extras = reduce_sdsspt.reduce_sdss(
            fn=fn,
            overscan=True, trim=True,
            subtract_bias=True,
            bias_hdu=bias_hdu,
            correct_flat=False,
        )
        if (valid_hdu is None):
            continue

        # normalize flatfield using the central 50%
        flatraw = valid_hdu['SCI'].data

        norm_area = flatraw[
            int(0.25*flatraw.shape[0]):int(0.75*flatraw.shape[0]),
            int(0.25*flatraw.shape[1]):int(0.75*flatraw.shape[1])]
        norm_intensity = numpy.median(norm_area)
        logger.debug("%s - mean flux: %.1f" % (fn, norm_intensity))

        if (norm_intensity < good_flux[0]):
            logger.warning("Excluding flat %s from mastercal, flux (%7.1f) to LOW" % (
                fn, norm_intensity
            ))
            continue
        elif (norm_intensity > good_flux[1]):
            logger.warning("Excluding flat %s from mastercal, flux (%7.1f) to HIGH" % (
                fn, norm_intensity
            ))
            continue

        flat_norm = flatraw / norm_intensity
        flat_norm[flat_norm < 0.1] = numpy.NaN

        if (write_norm_flat):
            pyfits.PrimaryHDU(data=flat_norm).writeto(fn[:-4]+".norm.fits", clobber=True)
        datablocks.append(flat_norm)

    masterflat = podi_imcombine.imcombine_data(
        datas=datablocks,
        operation='sigmaclipmean'
    )
    masterflat_hdu = pyfits.HDUList([
        pyfits.PrimaryHDU(),
        pyfits.ImageHDU(
            data=masterflat,
            name='SCI',
    )])

    return masterflat_hdu




def make_mastercals_from_filelist(filelist, cals_dir):

    #
    # Select all bias frames from file list
    #
    bias_list = []
    flat_list = {'u': [], 'g': [], 'r': [], 'i': [], 'z': []}

    for filename in filelist:

        hdulist = sdss2fits.open_sdss_fits(filename)
        if (hdulist is None):
            continue
        # with warnings.catch_warnings():
        #     warnings.simplefilter("ignore")
        #     hdulist = pyfits.open(filename)

        hdr = hdulist[0].header
        filtername = hdr['FILTER']
        obstype = hdr['FLAVOR']

        if (obstype.lower() == 'bias'):
            bias_list.append(filename)

        elif (obstype.lower() == 'flat'):
            flat_list[filtername].append(filename)


    if (len(bias_list) > 0):
        print "\n\nBIAS:\n--%s" % ("\n--".join(bias_list))
        bias_hdu = make_master_bias(bias_list)
        bias_out = "%s/masterbias.fits" % (cals_dir)
        bias_hdu.writeto(bias_out, clobber=True)


    for filtername in flat_list:

        filelist = flat_list[filtername]
        print "\n\nFLATS for %s:\n--%s" % (filtername, "\n--".join(filelist))

        if (len(filelist) <= 0):
            print "No files found!"
            continue

        masterflat_hdu = make_master_flat(filelist, bias_hdu=bias_hdu)

        flat_out = "%s/masterflat_%s.fits" % (cals_dir, filtername)
        print "Writing master-flat to %s" % (flat_out)
        masterflat_hdu.writeto(flat_out, clobber=True)



if __name__ == "__main__":

    dirname = sys.argv[1]
    filelist = glob.glob("%s/*.fit" % (dirname))
    print filelist

    cals_dir = sys.argv[2]

    make_mastercals_from_filelist(filelist, cals_dir)