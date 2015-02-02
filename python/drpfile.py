from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import sys
if sys.version > '3':
    long = int

import os.path
from os import environ
from astropy.io import fits
import time
import numpy

__author__ = 'Kyle Westfall'

def parse_drp_file_name(name):
    """
    Parse the name of a DRP file and return the plate, ifudesign, and
    mode.
    """

    if (type(name) != str):
        raise Exception("File name must be a string.")
    if (str.find(name, 'manga-') == -1 or str.find(name, '-LOG') == -1 or 
        str.find(name, '.fits.gz') == -1):
        raise Exception("String does not look like a DRP fits file name.")

    plate_start = str.find(name, '-')+1
    plate_end = str.find(name,'-',plate_start)
    plate = long(name[plate_start:plate_end])

    ifudesign_end = str.find(name,'-',plate_end+1)
    ifudesign = long(name[plate_end+1:ifudesign_end])

    mode_start = str.find(name,'LOG')+3
    mode = name[mode_start:str.find(name,'.fits.gz')]

    return plate, ifudesign, mode


def arginp_to_list(inp, evaluate=False, quiet=True):
    """
    Convert an input argument to a list of separated values.

    If evaluate is True, evaluate each element in the list before
    returning it.
    """

    out = inp

    # Simply return a None value
    if out is None:
        return out

    # If the input value is a string, convert it to a list of strings
    if type(out) == str:
        out = out.replace("[", "")
        out = out.replace("]", "")
        out = out.replace(" ", "")
        out = out.split(',')

    # If the input is still not a list, make it one
    if type(out) != list:
        out = [out]

    # Evaluate each string, if requested
    if evaluate:
        n = len(out)
        for i in range(0,n):
            try:
                tmp = eval(out[i])
            except (NameError, TypeError):
                if not quiet:
                    print('Could not evaluate value {0}. Skipping.'.format(out[i]))
            else:
                out[i] = tmp

    return out


def drpfile_list(drpver, platelist, ifudesignlist, modelist, combinatorics=False):
    """
    Provided a list of plates, ifu designs, and modes, return a list of
    DRP files to be analyzed by the DAP.

    ARGUMENTS:
        drpver              - The DRP version, which MUST be a single
                              string value used for all DRP files.
        platelist           - List of plates to use.
        ifudesignlist       - List of IFU designs to use.
        modelist            - List of DRP output modes ('CUBE' or 'RSS')
        combinatorics       - Based on the input, create a list with all
                              possible combinations.

    If the number of elements in each list is the same, the matching is
    assumed to be finished unless combinatorics is True.  If the number
    of elements is not the same, or cominatorics is True, the matched
    list is all the combinations of the provided elements.

    REVISION HISTORY:
        20 Nov 2014: (KBW) Original implementation
    """

    if (type(drpver) != str):
        raise Exception('DRP drpver must be a string')

    if platelist is None or ifudesignlist is None or modelist is None:
        raise ValueError('Must provide platelist, ifudesignlist, and modelist!')

    platelist_ = arginp_to_list(platelist, evaluate=True)
    ifudesignlist_ = arginp_to_list(ifudesignlist, evaluate=True)
    modelist_ = arginp_to_list(modelist)

    n_plates = len(platelist_)
    n_designs = len(ifudesignlist_)
    n_modes = len(modelist_)

    # Perform the combinatorics
    if (n_plates != n_designs or n_plates != n_modes or combinatorics is True):

        # Force elements in platelist, ifudesignlist, and modelist
        # to be unique before performing combinatorics
        platelist_ = list(set(platelist_))
        ifudesignlist_ = list(set(ifudesignlist_))
        modelist_ = list(set(modelist_))

        n_plates = len(platelist_)
        n_designs = len(ifudesignlist_)
        n_modes = len(modelist_)

        # Each element of plate list is repeated n_designs*n_modes times
        for i in range(0,n_plates):
            for j in range (1,n_designs*n_modes):
                platelist_.insert(i*n_designs*n_modes, platelist_[i*n_designs*n_modes])

        # First repeat each element of ifudesign list n_modes times
        for i in range(0,n_designs):
            for j in range (1,n_modes):
                ifudesignlist_.insert(i*n_modes, ifudesignlist_[i*n_modes])
        # Then repeat the entire list n_plates times
        ifudesignlist__ = ifudesignlist_[:]
        for i in range(1,n_plates):
            ifudesignlist_ += ifudesignlist__

        # Mode list iterated most quickly, so repeat the full list
        # n_plates*n_designs times
        modelist__ = modelist_[:]
        for i in range(1,n_plates*n_designs):
            modelist_ += modelist__

        nn = n_plates*n_designs*n_modes
    else:
        nn = n_plates

    # Create and return the list of DRP files
    return [drpfile(drpver, platelist_[i], ifudesignlist_[i], modelist_[i]) for i in range(0,nn)]



class drpfile:
    """
    Object to hold the properties of a DRP-produced file to be analyzed by the DAP.

    REVISION HISTORY:
        20 Nov 2014: (KBW) Original Implementation
    """


    def __init__(self, drpver, plate, ifudesign, mode):
        """
        ARGUMENTS:
            drpver - DRP version to use
            plate - plate designation
            ifudesign - ifu design on the plate
            mode - mode of the produced file (CUBE or RSS)
        """

        # Set the attributes, forcing a known type
        self.drpver = str(drpver)
        self.plate = long(plate)
        self.ifudesign = long(ifudesign)
        self.mode = str(mode)

        self.manga_id = None
        self.object_ra = None
        self.object_dec = None


    def _gather_data(self):
        """
        Internal support routine used to gather header information from
        the file.
        """
        if (self.manga_id is not None and self.object_ra is not None and
            self.object_dec is not None):
            return

        inp = self.file_path()
        if (not os.path.exists(inp)):
            raise Exception('Cannot open file: {0}'.format(inp))

        # TODO: This is very slow!  memmap in fits.open() does not help.
        # A way to read header without setting up the whole hdulist
        # object?  hdulist = fits.open(inp,memmap=True)
        hdr = fits.getheader(inp,0)
        self.manga_id = hdr['MANGAID']
        self.object_ra = hdr['OBJRA']
        self.object_dec = hdr['OBJDEC']

#        hdulist = fits.open(inp)
#        self.manga_id = hdulist[0].header['MANGAID']
#        self.object_ra = hdulist[0].header['OBJRA']
#        self.object_dec = hdulist[0].header['OBJDEC']
#        hdulist.close()


    def file_name(self):
        """Return the name of the DRP file"""
        return 'manga-{0}-{1}-LOG{2}.fits.gz'.format(self.plate, self.ifudesign, self.mode)


    def file_path(self):
        """Return the full path to the DRP file"""
        return os.path.join(environ['MANGA_SPECTRO_REDUX'], self.drpver, str(self.plate), 'stack',
                            self.file_name())


    def object_data(self):
        """Return some selected object data from the fits file"""
        self._gather_data()
        return self.manga_id, self.object_ra, self.object_dec
    

    def object_coo(self):
        """Return the object coordinates read from the fits header"""
        self._gather_data()
        return self.object_ra, self.object_dec

    def created_today(self):
        """Return True if the file was created today."""

        # Get the current time
        current_time = time.localtime()

        # Last time the DRP file was created (or modified?)
        created_time = time.gmtime(os.path.getctime(self.file_path()))

        return (current_time.tm_year == created_time.tm_year and
                current_time.tm_mon == created_time.tm_mon and 
                current_time.tm_mday == created_time.tm_mday)



