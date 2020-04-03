# HACKED VERSION OF LUCA'S CODE. ALL ADDTIONAL BAD CODING PRACTICE CURTOSEY OF C. ANGUS. SORRY!
from __future__ import print_function
import sys
import numpy

##### IMAGE QUERY selection
from astropy.table import Table
import requests
from io import BytesIO
import warnings
warnings.filterwarnings("ignore")

def search(snname,ra, dec, width):
    tab=panstarrs_query_pos(ra, dec, width)
    return(tab)

def getimages(ra,dec,size=240,filters="grizy"):
    """Query ps1filenames.py service to get a list of images

    ra, dec = position in degrees
    size = image size in pixels (0.25 arcsec/pixel)
    filters = string with filters to include
    Returns a table with the results
    """
    service = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
    url = ("{service}?ra={ra}&dec={dec}&size={size}&format=fits"
           "&filters={filters}").format(**locals())
    table = Table.read(url, format='ascii')
    return table


def geturl(ra, dec, size=240, output_size=None, filters="grizy", format="jpg", color=False):
    """Get URL for images in the table

    ra, dec = position in degrees
    size = extracted image size in pixels (0.25 arcsec/pixel)
    output_size = output (display) image size in pixels (default = size).
                  output_size has no effect for fits format images.
    filters = string with filters to include
    format = data format (options are "jpg", "png" or "fits")
    color = if True, creates a color image (only for jpg or png format).
            Default is return a list of URLs for single-filter grayscale images.
    Returns a string with the URL
    """
    if color and format == "fits":
        raise ValueError("color images are available only for jpg or png formats")
    if format not in ("jpg","png","fits"):
        raise ValueError("format must be one of jpg, png, fits")
    table = getimages(ra,dec,size=size,filters=filters)
    url = ("https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?"
           "ra={ra}&dec={dec}&size={size}&format={format}").format(**locals())
    if output_size:
        url = url + "&output_size={}".format(output_size)
    # sort filters from red to blue
    flist = ["yzirg".find(x) for x in table['filter']]
    table = table[numpy.argsort(flist)]
    if color:
        if len(table) > 3:
            # pick 3 filters
            table = table[[0,len(table)//2,len(table)-1]]
        for i, param in enumerate(["red","green","blue"]):
            url = url + "&{}={}".format(param,table['filename'][i])
    else:
        urlbase = url + "&red="
        url = []
        for filename in table['filename']:
            url.append(urlbase+filename)
    return url

##### DATA QUERY SECTION
from astroquery.vizier import Vizier
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import Column
import astropy
import numpy as np
import pandas as pd

def panstarrs_query_pos(ra, dec, width):
    v = Vizier(columns=['RAJ2000','DEJ2000','gmag','e_gmag','rmag','e_rmag','imag','e_imag','zmag','e_zmag','ymag','e_ymag'])  #setting Vizier for specific columns
    v.ROW_LIMIT = -1                                                                                           #no limits on the table length
    table_pre = v.query_region(astropy.coordinates.ICRS(ra=ra*u.deg, dec=dec*u.deg), width=width*u.arcmin,catalog='pan-STARRS')      #the real query
    #selection of lines with determined in magnitudes in I, R, G, Z
    if len(table_pre)!=0:
        index=[]
        for i in range(0, len(table_pre[0]['RAJ2000'].data)):
            index.append(i)
        #rename the columns and add seperation coloumn
        table=table_pre[0][index]
        sep=SkyCoord(ra,dec,frame='fk5',unit='deg').separation(SkyCoord(table['RAJ2000'],table['DEJ2000'],frame='fk5',unit='deg'))*3600
#        sep = Column(np.cos(np.arccos(dec-table['DEJ2000']) + (np.cos(dec)*np.cos(table['DEJ2000'])*np.arccos(ra-table['RAJ2000']))))
#        (np.arccos(np.sin(np.radians(dec))*(np.sin(np.radians(table['DEJ2000'])) + np.cos(np.radians(dec))*np.cos(np.radians(table['DEJ2000']))*np.cos(np.radians(abs(ra-table['RAJ2000'])))))))
        table.rename_column('RAJ2000', 'raMean')
        table.rename_column('DEJ2000', 'decMean')
        table.add_column(sep, name='Seperation')
        table.rename_column('gmag', 'g')
        table.rename_column('e_gmag', 'g_err')
        table.rename_column('rmag', 'r')
        table.rename_column('e_rmag', 'r_err')
        table.rename_column('imag', 'i')
        table.rename_column('e_imag', 'i_err')
        table.rename_column('zmag', 'z')
        table.rename_column('e_zmag', 'z_err')
        table.rename_column('ymag', 'y')
        table.rename_column('e_ymag', 'y_err')
        return table
    else:
        return table_pre


##### FINAL SECTION - DATA & IMAGE
from astropy.io import ascii
from astropy.wcs import WCS
from astropy.io import fits
from astropy.visualization import PercentileInterval, AsinhStretch
#visualization libraries
import matplotlib.pyplot as plt
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.visualization import SqrtStretch
from photutils import SkyCircularAperture

def get_image_datatab(snname,ra, dec, width):
    fitsurl = geturl(ra, dec, size=int(width*240), filters="r", format="fits")
    fh = fits.open(fitsurl[0])
    fhead = fh[0].header
    wcs = WCS(fhead)
    fim = fh[0].data
    # replace NaN values with zero for display
    fim[numpy.isnan(fim)] = 0.0
    # set contrast to something reasonable
    transform = AsinhStretch() + PercentileInterval(85)
    bfim = transform(fim)
    #query PS1 catalog
    datatab = panstarrs_query_pos(ra, dec, width)
    #query the table with positions and aperture pre-defined
    positions = SkyCoord(ra=datatab['raMean'],dec=datatab['decMean'], unit='deg')
    aper = SkyCircularAperture(positions, 0.5*u.arcsec)
    pix_aperture = aper.to_pixel(wcs)
    targ_coord= SkyCoord(ra=ra,dec=dec, unit='deg').to_pixel(wcs)
    #plot
    fig = plt.figure()
    fig.add_subplot(111, projection=wcs)
    #norm = ImageNormalize(stretch=SqrtStretch())
    plt.imshow(bfim, cmap='Greys', origin='lower')
    pix_aperture.plot(color='red', lw=1.5, alpha=1)
    plt.plot(targ_coord[0],targ_coord[1],color='limegreen',marker='+',lw=2,alpha=1)
    plt.xlabel('RA')
    plt.ylabel('Dec')
    image = plt.savefig('%s_hosts.png'%(snname), dpi=1000)
    # fits.writeto('test_run.fits', fim, fhead, overwrite=True)
    ascii.write(datatab, '%s_hosts.csv'%(snname), format='csv', fast_writer=False)
    return datatab

#main code
if __name__ == "__main__":
   # main()
     get_image_datatab(snname,ra, dec, width)
