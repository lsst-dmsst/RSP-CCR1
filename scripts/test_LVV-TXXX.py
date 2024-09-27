# Test script for LSST Data Management acceptance test campaign
# LVV-T1946:<https://jira.lsstcorp.org/secure/Tests.jspa#/testCase/LVV-T1946>
# LVV-T1946 tests verification element LVV-178, <https://jira.lsstcorp.org/browse/LVV-178>
# which verifies requirement DMS-REQ-0347-V-01: Measurements in catalogs
#
# The lsst environment must be set up to run this test,
# see <https://pipelines.lsst.io/install/setup.html> for more details

import os
import numpy as np
import astropy.units as u

from lsst.daf.butler import Butler

# Confirm the version of the Science Pipelines:
os.system('eups list -s | grep lsst_distrib')

# For RC2 data on the USDF:
config = '/repo/main'
collection = 'HSC/runs/RC2/w_2023_39/DM-40985'
butler = Butler(config, collections=collection)

def LVVT1946(instrument, tract, patch, band, skymap):
    """ Execute test case LVV-T1946

    Parameters
    ----------
    tract: `int`
        Input tract for test.
    patch: `int`
        Input patch for test.
    band: `str`
        Band on which to test schema.
    """

    dataIdCoadd = {'instrument':instrument, 'tract':tract, 'band':band,
                   'patch':patch, 'skymap':skymap}
    print('Input dataId: ', dataIdCoadd)

    try:
        forced_src = butler.get('deepCoadd_forced_src', dataId = dataIdCoadd)
    except:
        print('Invalid dataId. Please check and try again.')

    fsch = forced_src.schema
    forced_src_flag = check_schema(fsch)

    # If any of the units are not "count," print "FALSE."
    assert(np.all(forced_src_flag)), 'FALSE: not all instFlux entries have units of counts.'
    
    # If all of the flux units are "count," print "TRUE."
    print('\n All forced_src instFlux entries have units of counts: ', np.all(forced_src_flag))

    try:
        obj_table = butler.get('objectTable', dataId = dataIdCoadd)
    except:
        print('Invalid dataId. Please check and try again.')

    tractname = dataIdCoadd['tract']
    print('\n-------------------------------------------------')
    print(f'\nChecking flux values for visit {tractname}')
    check_flux_values(obj_table)


def check_flux_values(obj_table):
    """ Check that the flux values are within a range of reasonable
    values, and print the results to the screen

    Parameters
    ----------
    table: Source Table
        src table to check
    """

    bands = ['g', 'r', 'i', 'z', 'y']
    fluxes = ['psfFlux', 'free_psfFlux', 'gaapPsfFlux',
              'gaap0p5Flux', 'gaap0p7Flux', 'gaap1p0Flux',
              'gaap1p5Flux', 'gaap2p5Flux', 'gaap3p0Flux',
              'gaapOptimalFlux', 'kronFlux', 'calibFlux',
              'ap03Flux', 'ap06Flux', 'ap09Flux', 'ap12Flux',
              'ap17Flux', 'ap25Flux', 'ap35Flux', 'ap50Flux',
              'ap70Flux', 'cModelFlux', 'cModelFlux_inner',
              'free_cModelFlux', 'free_cModelFlux_inner']

    bright_maglim = -5.0*u.ABmag
    bright_fluxlim = bright_maglim.to(u.nJy)
    faint_fluxlim = -1e6*u.nJy

    for band in bands:
        # Flags to apply for all flux columns:
        okisprimary = obj_table['detect_isPrimary'] == 0
        # okpixflags = src_table[band+"_psfFlux_flag"] == 0
        okpixflags = obj_table[band+'_pixelFlags_saturatedCenter'] == 0
        okcentroidflags = obj_table['xy_flag'] == 0
        notskysource = obj_table['sky_object'] == 0

        print('\nflux column          % of bad flux values\n')
        for fluxcol_stub in fluxes:
            fluxcol = band+'_'+fluxcol_stub
            toofaint = (obj_table[fluxcol].values < faint_fluxlim.value)
            toobright = (obj_table[fluxcol].values > bright_fluxlim.value)
            if (fluxcol == 'calibFlux'):
                okfluxflag = obj_table['ap12Flux_flag'] == 0
            elif ('gaap' in fluxcol):
                okfluxflag = obj_table[band+'_gaapFlux_flag'] == 0
            elif ('cModel' in fluxcol):
                okfluxflag = obj_table[band+'_kronFlux_flag'] == 0
            else:
                okfluxflag = obj_table[fluxcol+'_flag'] == 0
            fluxflag = (toofaint | toobright) & okisprimary & okfluxflag & okpixflags & okcentroidflags & notskysource
            bad_percentage = (np.sum(fluxflag)/len(fluxflag))*100*u.percent
            print(f'{fluxcol:20} {bad_percentage:7.2}')


def check_schema(sch):
    """ Check the units on each flux entry in the schema

    Parameters
    ----------
    sch: Schema 
        forced_src schema to check

    Returns
    -------
    forced_src_flag: `list` of `bool`
        Boolean list containing "True" for elements with flux units in "count"
    """

    forced_src_flag = []
    for entry in sch:
        field = entry.getField().getName()
        if 'instFlux' in field and 'flag' not in field and 'Cov' not in field:
            flux_units = entry.getField().getUnits()
            forced_src_flag.append('count' in flux_units)
            print(f'{field} ..... {flux_units:20}')
    return forced_src_flag


if __name__ == "__main__":

    # Select an arbitrary source catalog from a deepCoadd:
    instrument = 'HSC'
    band = 'i'
    tract = 9615
    patch = 13
    skymap = 'hsc_rings_v1'

    # Run the test
    LVVT1946(instrument, tract, patch, band, skymap)
