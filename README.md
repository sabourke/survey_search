# `survey_search`

This is a tool to get the status of a target in the LOFAR Survey (LoTSS).
I.E. what survey field the target is in and if it is observed / processed.
It can attempt to download the data if it exists. It can also generate
cutouts of a user defined size centred on the target.

## Installing

Dependancies: `git python-astropy python-astroquery wget libcfitsio-bin`
```
sudo apt install git python-astropy python-astroquery wget libcfitsio-bin
```

Copy the `survey_search.py` and `fix_fitscopy_fz.py` programs to a
directory in your path.
```
mkdir -p ~/.local/bin
cp survey_search.py fix_fitscopy_fz.py ~/.local/bin
```

To setup passwordless downloading of survey data create a `.netrc` file
in your home directory and add the following line:
```
machine lofar-surveys.org login surveys password survey_password
```
replacing `survey_password` with the actual password. You will need to
do this if you want to use the automated download feature of this program.

## Usage
```
usage: survey_search.py [-h] [-s] [-c angle] target

Get field and pixel position of a target.

positional arguments:
  target                target position or name

optional arguments:
  -h, --help            show this help message and exit
  -s, --simbad          target is a simbad
  -c angle, --cutout angle
                        cutout size. eg: "30arcmin"
  -d, --download        run commands to retrieve data

```

The `fitscopy` utility can work with compressed cubes and make
cutouts avoiding decompression of the full cube. The output cubes
seem to be non-standard. The `fix_fitscopy_fz.py` program is
provided to fix these output cubes.

### Example:
```
$ survey_search.py 13:29:53+47:11:43
$ survey_search.py --cutout 30arcmin --download --simbad m51
```
