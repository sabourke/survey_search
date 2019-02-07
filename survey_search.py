#!/usr/bin/env python

# Locate LoTSS DR2 data based on coords or common name
# Stephen Bourke
# Onsala Space Observatory / Chalmers
# November 2018

from __future__ import print_function
import numpy as np
from astropy.coordinates import Angle, SkyCoord
import astropy.units as u
from astroquery.simbad import Simbad
from astropy.wcs import WCS
from astropy.io import fits
import sys
import argparse
import os.path
import json
try:
    # Python3
    from urllib.request import urlopen
except ImportError:
    # Python2
    from urllib2 import urlopen
import subprocess


pointings_url = "http://skogul.oso.chalmers.se/Lofar_Tier1_Survey/status/pointings_db.json"


def field_coords_rad():
    """Return arrays of 'field names', 'pointing centre', 'status'"""
    pointings = json.load(urlopen(pointings_url))
    ids = [p[0] for p in pointings]
    ra_dec = [p[1:3] for p in pointings]
    status = [p[3] for p in pointings]
    return np.array(ids), np.radians(np.array(ra_dec)), np.array(status)


def simbad_object_skycoord(target):
    """Return a SkyCoord object for a target"""
    cs = Simbad()
    cs.ROW_LIMIT = 1
    cs.TIMEOUT = 5

    result = cs.query_object(target)
    try:
        assert len(result) == 1
    except:
        raise LookupError("Simbad Object Query failed")
    result = result[0]

    return SkyCoord(result['RA'], result['DEC'], unit=(u.hourangle, u.degree))


def cos_distance_rad(p0, p1):
    """Return the distance between args. Inputs can be np.arrays"""
    r0, d0 = p0
    r1, d1 = p1
    cos_d = np.sin(d0) * np.sin(d1) +  np.cos(d0) * np.cos(d1) * np.cos(r0 - r1)
    return cos_d


def get_fields(target_skycoord):
    """Return sorted list of fields that contain the target in the ~ FWHM beam"""
    ids, coords, status = field_coords_rad()
    p0 = (target_skycoord.ra.radian, target_skycoord.dec.radian)
    d = np.degrees(np.arccos(cos_distance_rad(p0, coords.T)))
    d_sort = np.argsort(d)
    d_sort = d_sort[d[d_sort] < 1.75]
    return ids[d_sort], coords[d_sort], d[d_sort], status[d_sort]


def survey_wcs():
    """Returns a WCS centred on 0,0 with the geometry of a typcial LoTSS field"""
    header = fits.Header()
    header["NAXIS"] = 2
    header["NAXIS1"] = 6075
    header["CTYPE1"] = "RA---SIN"
    header["CRVAL1"] = 0.0
    header["CRPIX1"] = 3038.0
    header["CDELT1"] = -0.00125
    header["NAXIS2"] = 6075
    header["CTYPE2"] = "DEC--SIN"
    header["CRVAL2"] = 0.0
    header["CRPIX2"] = 3038.0
    header["CDELT2"] = 0.00125
    return WCS(header)


def download_commands(field, outname="crop", cutout=""):
    url_i = "https://lofar-surveys.org/downloads/DR2/fields/{field}/image_full_low_m.int.restored.fits".format(field=field)
    url_qu = "https://lofar-surveys.org/downloads/DR2/fields/{field}/image_full_low_QU.cube.dirty.corr.fits.fz".format(field=field)
    outname_i = "{field}_low_I_corr.fits".format(field=field)
    outname_qu = "{field}_low_QU_corr.fits.fz".format(field=field)
    commands = []
    commands.append("wget -O {} {}".format(outname_i, url_i))
    commands.append("wget -O {} {}".format(outname_qu, url_qu))
    if cutout != "":
        outname_i_crop = "{name}_low_I_corr.fits".format(name=outname)
        outname_qu_crop = "{name}_low_QU_corr.fits.fz".format(name=outname)
        commands.append("fitscopy {}{} {}".format(outname_i, cutout, outname_i_crop))
        commands.append("fitscopy {}{} {}".format(outname_qu, cutout, outname_qu_crop))
        commands.append("fix_fitscopy_fz.py {} {}".format(outname_qu_crop, os.path.splitext(outname_qu_crop)[0]))
    return commands


def run_commands(commands):
    for line in commands:
        subprocess.check_call(line, shell=True)

def main():
    parser = argparse.ArgumentParser(description='Get field and pixel position of a target.')
    parser.add_argument('-s', '--simbad', action='store_true', default=False, help='target is a simbad')
    parser.add_argument('-c', '--cutout', metavar='angle', type=str, default="", help='cutout size. eg: "30arcmin"')
    parser.add_argument('-d', '--download', action='store_true', default=False, help='run commands to retrieve data')
    parser.add_argument('target', type=str, help="target position or name")
    args = parser.parse_args()

    if args.simbad:
        target_sc = simbad_object_skycoord(args.target)
    else:
        target_sc = SkyCoord(args.target, unit=(u.hourangle, u.degree))

    id_list, coord_list, distance_list, status_list = get_fields(target_sc)
    for i in range(len(id_list)):
        field_name = id_list[i]
        coords = coord_list[i]
        distance = distance_list[i]
        status = status_list[i]

        field_sc = SkyCoord(coords[0]*u.radian, coords[1]*u.radian)

        if not args.cutout:
            target_offset = target_sc.transform_to(field_sc.skyoffset_frame())
            pixels = SkyCoord(target_offset.lon, target_offset.lat).to_pixel(survey_wcs())
            pixels = "[{}:{}]".format(*[x.round().astype(int) for x in pixels])
            print("# {:>11} {:>11}, Distance: {:3.2f}deg, Status: {}".format(field_name, pixels, distance, status))
            commands = download_commands(field_name, "{}_{}".format(args.target, i))
        else:
            cutout = Angle(args.cutout)
            pos = SkyCoord(cutout/2, -cutout/2, frame=target_sc.skyoffset_frame())
            pos = pos.transform_to(field_sc.skyoffset_frame())
            pixels = SkyCoord(pos.lon, pos.lat).to_pixel(survey_wcs())
            ra0, dec0 = [x.round().astype(int) for x in pixels]

            pos = SkyCoord(-cutout/2, cutout/2, frame=target_sc.skyoffset_frame())
            pos = pos.transform_to(field_sc.skyoffset_frame())
            pixels = SkyCoord(pos.lon, pos.lat).to_pixel(survey_wcs())
            ra1, dec1 = [x.round().astype(int) for x in pixels]

            region = "[{}:{},{}:{}]".format(ra0, ra1, dec0, dec1)
            print("# {:>11} {:>21}, Distance: {:3.2f}deg, Status: {}".format(field_name, region, distance, status))
            commands = download_commands(field_name, "{}_{}".format(args.target, i), region)

        if args.download:
            run_commands(commands)
        else:
            [print(c) for c in commands]
 

if __name__ == "__main__":
    main()
