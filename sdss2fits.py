#!/usr/bin/env python

import os
import sys
import pyfits

import warnings
import tempfile
import numpy
import gzip

import logging

from optparse import OptionParser

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

    #return hdulist
    hdulist.writeto(out_file, clobber=True)




def open_sdss_fits(filename):

    logger = logging.getLogger("OpenSDSS")

    tmp_file = tempfile.TemporaryFile(
        suffix=".fits", prefix="sdsspt_",
    )
    if (filename.endswith(".gz")):
        unpack_fn = tempfile.NamedTemporaryFile(
            suffix=".fits", prefix="gunzip_sdss"
        )
        with gzip.open(filename,"r") as gz:
            unpack_fn.write(gz.read())
        hdu_raw = convert_sdss(unpack_fn.name, tmp_file)

    else:
        #print tmp_file
        hdu_raw = convert_sdss(filename, tmp_file)
        #hdu_raw.writeto(tmp_file) #, mode='spool')

    tmp_file.seek(0)
    hdulist = pyfits.open(tmp_file)

    #
    # Now we have a very basic FITS-compatible image.
    # separate image data from primary header to make things more compatible
    # with the quickreduce workflow.
    #

    #
    # Fix the data to account for it being UNSIGNED rather than signed int
    #
    data = hdulist[0].data.copy().astype(numpy.int32)
    # print data.dtype
    # data[data<0] *= -1
    data[data<0] += 2**16
    #data = data.astype(numpy.uint16)
    #data = data.astype(numpy.uint16)
    #data.dtype = numpy.uint16

    primhdu = pyfits.PrimaryHDU(header=hdulist[0].header)
    imagehdu = pyfits.ImageHDU(
        data=data,
        name="SCI"
    )

    pixelscale = 1.161 / 3600.
    if ('CROTA2' in hdulist[0].header):
        crota = hdulist[0].header['CROTA2']
    else:
        crota = 0
    logger.debug("setting CROTA paramter to %d" % (crota))
    imagehdu.header['CRVAL1'] = hdulist[0].header['RADEG']
    imagehdu.header['CRVAL2'] = hdulist[0].header['DECDEG']
    imagehdu.header['CTYPE1'] = 'RA---TAN'
    imagehdu.header['CTYPE2'] = 'DEC--TAN'
    imagehdu.header['CD1_1'] = pixelscale*numpy.cos(numpy.radians(crota))
    imagehdu.header['CD2_2'] = pixelscale*numpy.cos(numpy.radians(crota))
    imagehdu.header['CD1_2'] = -pixelscale*numpy.sin(numpy.radians(crota))
    imagehdu.header['CD2_1'] = pixelscale*numpy.sin(numpy.radians(crota))
    imagehdu.header['CRPIX1'] = hdulist[0].header['NAXIS1']/2
    imagehdu.header['CRPIX2'] = hdulist[0].header['NAXIS2']/2

    imagehdu.header['OBJECT'] = \
        "%s -- %s" % (hdulist[0].header['NAME'], hdulist[0].header['FILTER'])

    #
    # Delete WCS-related keywords from primary header -
    # these go in the ImageHDU header
    #
    for hdrkey in [
        'CRPIX1', 'CRPIX2',
        'CRVAL1', 'CRVAL2',
        'CDELT1', 'CDELT2',
        'CROTA2',
        'CTYPE1', 'CTYPE2',]:
        if hdrkey in primhdu.header:
            del primhdu.header[hdrkey]


    # data = hdulist[0].data
    # print data.shape, data.dtype
    #
    # data.dtype=numpy.int16
    # ny = hdulist[0].header['NAXIS2']
    # nx = hdulist[0].header['NAXIS1']
    # data = data.reshape((ny,nx))
    # print data.shape, data.dtype
    # imagehdu = pyfits.ImageHDU(data=hdulist[0].data)

    # hdulist[0].header['SIMPLE'] = True
    # imagehdu.header['BITPIX'] = 16
    # imagehdu.header['NAXIS'] = 2
    # del hdulist[0].header['SDSS']
    # del hdulist[0].header['UNSIGNED']



    # pixelscale = 1.16 / 3600.
    #
    # hdulist[0].header['NAXIS'] = 0
    # del hdulist[0].header['NAXIS1']
    # del hdulist[0].header['NAXIS2']

    #return hdulist
    # primhdu = pyfits.PrimaryHDU(header=hdulist[0].header)
    out_hdulist = pyfits.HDUList([primhdu, imagehdu])
    # hdulist.writeto(out_file, clobber=True)
    #out_hdulist.writeto(out_file, clobber=True)



    #print "opening file"
    #hdulist = pyfits.open(tmp_file)

    return out_hdulist


if __name__ =="__main__":


    parser = OptionParser()
    parser.add_option("", "--show", dest="show",
                       action="store_true", default=False)
    (options, cmdline_args) = parser.parse_args()


    show_list = []
    for infile in cmdline_args[0:]:

        outfile = infile[:-4]+".fits"
        if (os.path.isfile(outfile)):
            continue
        print "%s ==> %s" % (infile, outfile)
        # convert_sdss(infile, outfile)

        hdu = open_sdss_fits(infile)
        hdu.writeto(outfile, clobber=True)

        if (options.show):
            show_list.append(outfile)

    if (options.show):
        os.system("ds9 %s &" % (" ".join(show_list)))
