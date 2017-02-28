#!/usr/bin/env python

import os
import sys
import pyfits

import warnings
import tempfile


scratch_dir = "/tmp"

def convert_sdss(in_file, out_file):

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        hdulist = pyfits.open(in_file)

    hdulist[0].header['SIMPLE'] = True
    hdulist[0].header['BITPIX'] = 16
    hdulist[0].header['NAXIS'] = 2
    del hdulist[0].header['SDSS']
    del hdulist[0].header['UNSIGNED']

    pixelscale = 1.16 / 3600.
    hdulist[0].header['CRVAL1'] = hdulist[0].header['RADEG']
    hdulist[0].header['CRVAL2'] = hdulist[0].header['DECDEG']
    hdulist[0].header['CTYPE1'] = 'RA---TAN'
    hdulist[0].header['CTYPE2'] = 'DEC--TAN'
    hdulist[0].header['CD1_1'] = pixelscale
    hdulist[0].header['CD2_2'] = pixelscale
    hdulist[0].header['CRPIX1'] = hdulist[0].header['NAXIS1']/2
    hdulist[0].header['CRPIX2'] = hdulist[0].header['NAXIS2']/2

    hdulist[0].header['OBJECT'] = \
        "%s -- %s" % (hdulist[0].header['NAME'], hdulist[0].header['FILTER'])

    #return hdulist
    hdulist.writeto(out_file, clobber=True)


def open_sdss_fits(filename):

    tmp_file = tempfile.TemporaryFile(
        suffix=".fits", prefix="sdsspt_",
    )
    #print tmp_file
    hdu_raw = convert_sdss(filename, tmp_file)
    #hdu_raw.writeto(tmp_file) #, mode='spool')

    tmp_file.seek(0)

    #print "opening file"
    hdulist = pyfits.open(tmp_file)

    return hdulist


if __name__ =="__main__":


    infile = sys.argv[1]

    for infile in sys.argv[1:]:
        outfile = infile[:-4]+".fits"
        if (os.path.isfile(outfile)):
            continue
        print "%s ==> %s" % (infile, outfile)
        convert_sdss(infile, outfile)
