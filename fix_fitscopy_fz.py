#!/usr/bin/env python

# fitscopy can workwith fz compressed cubes but seems to produce
# strange output files.
# This program attempts to fix such files.
#
# Stephen Bourke
# Onsala Space Observatory / Chalmers
# November 2018

from astropy.io import fits
import argparse


def fix_fitscopy_fz(infits, outfits):
    hdu_in = fits.open(infits)
    hdu_out = fits.PrimaryHDU()
    hdu_out.header = hdu_in[1].header
    hdu_out.header.rename_keyword("XTENSION", "SIMPLE")
    hdu_out.header.set("SIMPLE", True, "file does conform to FITS standard")
    for keyword in ("ZQUANTIZ", "ZDITHER0", "CHECKSUM", "DATASUM"):
        if keyword in hdu_out.header:
            del hdu_out.header[keyword]
    
    hdu_out.data = hdu_in[1].data
    hdu_out.writeto(outfits)


def main():
    parser = argparse.ArgumentParser(description='Fix FITS fz cube after running fitscopy')
    parser.add_argument('incube', type=str, help="Input FITS cube")
    parser.add_argument('outcube', type=str, help="Output FITS cube")
    args = parser.parse_args()
    fix_fitscopy_fz(args.incube, args.outcube)
    

if __name__ == "__main__":
    main()
