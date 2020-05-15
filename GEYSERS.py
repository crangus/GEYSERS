#Getting Exciting Young Supernova Experiment Recent Supernovae (GEYSERS)
#
# You need to run the SQL Explorer Script "DARK_YSE_Recent_Transients" and download the csv file this produces.
#

from __future__ import print_function
import sys
sys.path.insert(0, './ps1/')
import argparse
import numpy as np
import scipy as sp
import pandas as pd
from ps1 import PS1_query
import astropy
from astropy.time import TimeISO
from astropy.time import Time
from astropy.cosmology import FlatLambdaCDM 
import astropy.units as u

cosmo = FlatLambdaCDM(H0=70.0, Om0=0.3)
CGS_H=6.6260755e-27
CGS_C=2.99792458e10
CGS_K=1.380658e-16 


#weighting for grouped photometry
g_dawm = lambda x: np.average(x, weights=g_da.loc[x.index, "mag_err"])
r_dawm = lambda x: np.average(x, weights=r_da.loc[x.index, "mag_err"])
i_dawm = lambda x: np.average(x, weights=i_da.loc[x.index, "mag_err"])
z_dawm = lambda x: np.average(x, weights=z_da.loc[x.index, "mag_err"])
ztfg_dawm = lambda x: np.average(x, weights=ztfg_da.loc[x.index, "mag_err"])
ztfr_dawm = lambda x: np.average(x, weights=ztfr_da.loc[x.index, "mag_err"])
ztfi_dawm = lambda x: np.average(x, weights=ztfi_da.loc[x.index, "mag_err"])
ztfz_dawm = lambda x: np.average(x, weights=ztfz_da.loc[x.index, "mag_err"])

parser = argparse.ArgumentParser(description='Transient properties:')
parser.add_argument('file', type=str, help='input csv file',default=None)

#g colour args
parser.add_argument('--grmin', type=float, help='g-r colour (minimum)',default=None)
parser.add_argument('--grmax', type=float, help='g-r colour (maximum)',default=None)
parser.add_argument('--gimin', type=float, help='g-i colour (minimum)',default=None)
parser.add_argument('--gimax', type=float, help='g-i colour (maximum)',default=None)
parser.add_argument('--gzmin', type=float, help='g-z colour (minimum)',default=None)
parser.add_argument('--gzmax', type=float, help='g-z colour (maximum)',default=None)

# r colour args
parser.add_argument('--rimin', type=float, help='r-i colour (minimum)',default=None)
parser.add_argument('--rimax', type=float, help='r-i colour (maximum)',default=None)
parser.add_argument('--rzmin', type=float, help='r-z colour (minimum)',default=None)
parser.add_argument('--rzmax', type=float, help='r-z colour (maximum)',default=None)

#i colour args
parser.add_argument('--izmin', type=float, help='i-z colour (minimum)',default=None)
parser.add_argument('--izmax', type=float, help='i-z colour (maximum)',default=None)

#rise args
parser.add_argument('--risemin', type=float, help='rise rate (minimum). Use negative values for rising, positive for declining. Must also specify photometric band.',default=None)
parser.add_argument('--risemax', type=float, help='rise rate (maximum). Use negative values for rising, positive for declining. Must also specify photometric band.',default=None)
parser.add_argument('--riseband', type=str, help='rise band (g,r,i,z)',default=None)

#distance args
parser.add_argument('--photoz', type=float , help='photo-z for apparent host in YSE-PZ',default=None)
parser.add_argument('--photozrange', type=float , help='photo-z range for apparent host in YSE-PZ (defualt = +/- 0.05)',default=0.05)
parser.add_argument('--absmag', type=float , help='minimum peak absolute magnitude for transient based on photo-z (must include absmagband option)',default=None)
parser.add_argument('--absmagband', type=str , help=' photometric band for absolute magnitude limit for transient based on photo-z ("g"/"r"/"i"/"z")',default=None)
parser.add_argument('--absmagrange', type=float , help='range around absolute magnitude for transient based on photo-z (absmag and absmagband options, defualt = +/- 0.1)',default=0.1)
parser.add_argument('--absmaglowlim', type=str , help='treat absolute magnitude as lower limit ("True"/"False")',default='False')

#host args
parser.add_argument('--host', type=str , help='host environment ("host"/"nuclear"/"orphan"/"any")',default=None)
parser.add_argument('--hostmagdiff', type=float , help='transient current brightness above host environment (requires "host" and "hostband" options)',default=None)
parser.add_argument('--hostband', type=str , help='host environment band (required for "hostmagdiff" option)',default=None)

#photometry/followup args
parser.add_argument('--maglim', type=float, help='faintest objects of interest (all bands). Default is 20.0',default=20.0)
parser.add_argument('--telescope', type=float, help='Telescope for spectroscopic observations ("Asiago"/"NOT"/"VLT")',default=None)
parser.add_argument('--showall', type=str, help='Show good and partial candidate matches ("True"/"False")',default='False')

condits=parser.parse_args()

#read in photometry
File = condits.file
file=File#.decode(encoding='UTF-8',errors='strict')
g_r_min=condits.grmin
g_r_max=condits.grmax
g_i_min=condits.gimin
g_i_max=condits.gimax
g_z_min=condits.gzmin
g_z_max=condits.gzmax
r_i_min=condits.rimin
r_i_max=condits.rimax
r_z_min=condits.rzmin
r_z_max=condits.rzmax
i_z_min=condits.izmin
i_z_max=condits.izmax
risemin=condits.risemin
risemax=condits.risemax
riseband=condits.riseband
hostenv=condits.host
hostmagdiff=condits.hostmagdiff
hostband=condits.hostband
photoz=condits.photoz
photozrange=condits.photozrange
absmag=condits.absmag
absmagband=condits.absmagband
absmagrange=condits.absmagrange
absmaglowlim=condits.absmaglowlim
maglim=condits.maglim
telescope=condits.telescope
showall=condits.showall

Phot=pd.read_csv(file,sep=',',header=0)

#Convert date to MJD
Phot.observed_date=Time([x for x in Phot.observed_date.values]).mjd

#Drop bad mags, mjds
indexmjds = Phot[Phot['observed_date']==0].index
Phot.drop(indexmjds , inplace=True)
indexnas = Phot[Phot['mag_err'].isna()==True].index
Phot.drop(indexnas , inplace=True)
photometry=Phot
photometry.drop(['count_date'], axis=1,inplace=True)

#Drop rows with gaia or "unknown" photometric bands
indexGaia = photometry[photometry['filter'] == 'G' ].index
photometry.drop(indexGaia , inplace=True)
indexUnknown = photometry[photometry['filter'] == 'Unknown' ].index
photometry.drop(indexUnknown , inplace=True)
indexbadphot = photometry[photometry['data_quality_id'] == 1 ].index
photometry.drop(indexbadphot , inplace=True)
photometry.drop(['data_quality_id'], axis=1,inplace=True)

g_colours=np.array([g_r_min,g_r_max,g_i_min,g_i_max,g_z_min,g_z_max])
r_colours=np.array([r_i_min,r_i_max,r_z_min,r_z_max])
i_colours=np.array([i_z_min,i_z_max])
all_colours=np.array([g_r_min,g_r_max,g_i_min,g_i_max,g_z_min,g_z_max,r_i_min,r_i_max,r_z_min,r_z_max,i_z_min,i_z_max])
rises=np.array([risemin,risemax])
distances=np.array([absmag,photoz])
followup_resources=["NOT","VLT","Asiago"]

#output lists
col_key_objects=[]
col_flagged=[]
rise_key_objects=[]
rise_flagged=[]
host_key_objects=[]
host_flagged=[]
no_host_flag=[]
nuclear=[]
hostless=[]
redshift_range=[]
bright=[]
bright_range=[]
dist_key_objects=[]

#parameter checks
if len(rises[rises!=None])>0 and riseband==None:
    print()
    print('Error: You need to enter the rise/delcine photometric band (g/r/i/z)')
    print()
    exit()

if absmag!=None and absmagband==None:
    print()
    print('Error: You need to enter the photometric band for absolute magnitude (g/r/i/z)')
    print()
    exit()

if absmag!=None and absmagband==None:
    print()
    print('Error: You need to enter the photometric band for absolute magnitude (g/r/i/z)')
    print()
    exit()
    
if hostmagdiff!=None and hostenv==None:
    print()
    print('Error: You must include a host environment')
    print()
    exit()

if hostmagdiff!=None and hostband==None:
    print()
    print('Error: You must include a photometric band for the host')
    print()
    exit()

if telescope!=None and telescope not in followup_resources:
    print()
    print('Current follow-up resources incl.: NOT, VLT, Asiago')
    print()
    exit()

if telescope=='NOT':
    maglim=20.5
    hostmagdiff=1.0
    hostband='r'

if telescope=='Asiago':
    maglim=18.5
    hostmagdiff=1.0
    hostband='r'

if telescope=='VLT':
    maglim=21.5
    hostmagdiff=1.0
    hostband='r'

#Do some shit
candidates=photometry.groupby('name')

print("Ingesting %i recent candidates from YSE"%(len(candidates.groups.keys())))

#Colour cuts & rise
if len(all_colours[all_colours!= None])>0:
    print('Checking candidate photometry')
    for sn in candidates.groups.keys():
        g_dat=list()
        r_dat=list()
        i_dat=list()
        z_dat=list()
        tempdata=candidates.get_group(sn).reset_index()
    
        phot=tempdata.groupby('filter')

        if 'g' in phot.groups.keys():
            g_da=phot.get_group('g')
            tempg=g_da.groupby(np.round(g_da['observed_date']),as_index=False).agg({'ra':lambda r: np.median(r),'dec':lambda d: np.median(d),'observed_date':lambda z: np.median(z),'mag': g_dawm, 'mag_err': lambda   x: np.sqrt(np.sum(x**2))/x.count()})
            tempg.insert(0,'name',sn)
            tempg.insert(4,'filter','g')
            g_dat.append(tempg)
        if 'g-ZTF' in phot.groups.keys():
            ztfg_da=phot.get_group('g-ZTF')
            tempztfg=ztfg_da.groupby(np.round(ztfg_da['observed_date']),as_index=False).agg({'ra':lambda r: np.median(r),'dec':lambda d: np.median(d),'observed_date':lambda z: np.median(z),'mag': ztfg_dawm, 'mag_err': lambda   x: np.sqrt(np.sum(x**2))/x.count()})
            tempztfg.insert(0,'name',sn)
            tempztfg.insert(4,'filter','g')
            g_dat.append(tempztfg)
        if 'r' in phot.groups.keys():
            r_da=phot.get_group('r')
            tempr=r_da.groupby(np.round(r_da['observed_date']),as_index=False).agg({'ra':lambda r: np.median(r),'dec':lambda d: np.median(d),'observed_date':lambda z: np.median(z),'mag': lambda y: np.median(y), 'mag_err': lambda   x: np.sqrt(np.sum(x**2))/x.count()})
            tempr.insert(0,'name',sn)
            tempr.insert(4,'filter','r')
            r_dat.append(tempr)
        if 'r-ZTF' in phot.groups.keys():
            ztfr_da=phot.get_group('r-ZTF')
            tempztfr=ztfr_da.groupby(np.round(ztfr_da['observed_date']),as_index=False).agg({'ra':lambda r: np.median(r),'dec':lambda d: np.median(d),'observed_date':lambda z: np.median(z),'mag': lambda y: np.median(y), 'mag_err': lambda   x: np.sqrt(np.sum(x**2))/x.count()})
            tempztfr.insert(0,'name',sn)
            tempztfr.insert(4,'filter','r')
            r_dat.append(tempztfr)
        if 'i' in phot.groups.keys():
            i_da=phot.get_group('i')
            tempi=i_da.groupby(np.round(i_da['observed_date']),as_index=False).agg({'ra':lambda r: np.median(r),'dec':lambda d: np.median(d),'observed_date':lambda z: np.median(z),'mag': lambda y: np.median(y), 'mag_err': lambda   x: np.sqrt(np.sum(x**2))/x.count()})
            tempi.insert(0,'name',sn)
            tempi.insert(4,'filter','i')
            i_dat.append(tempi)
        if 'i-ZTF' in phot.groups.keys():
            ztfi_da=phot.get_group('i-ZTF')
            tempztfi=ztfi_da.groupby(np.round(ztfi_da['observed_date']),as_index=False).agg({'ra':lambda r: np.median(r),'dec':lambda d: np.median(d),'observed_date':lambda z: np.median(z),'mag': lambda y: np.median(y), 'mag_err': lambda   x: np.sqrt(np.sum(x**2))/x.count()})
            tempztfi.insert(0,'name',sn)
            tempztfi.insert(4,'filter','i')
            i_dat.append(tempztfi)
        if 'z' in phot.groups.keys():
            z_da=phot.get_group('z')
            tempz=z_da.groupby(np.round(z_da['observed_date']),as_index=False).agg({'ra':lambda r: np.median(r),'dec':lambda d: np.median(d),'observed_date':lambda z: np.median(z),'mag': lambda y: np.median(y), 'mag_err': lambda   x: np.sqrt(np.sum(x**2))/x.count()})
            tempz.insert(0,'name',sn)
            tempz.insert(4,'filter','z')
            z_dat.append(tempz)
        if 'z-ZTF' in phot.groups.keys():
            ztfz_da=phot.get_group('z-ZTF')
            tempztf=ztfz_da.groupby(np.round(ztfz_da['observed_date']),as_index=False).agg({'ra':lambda r: np.median(r),'dec':lambda d: np.median(d),'observed_date':lambda z: np.median(z),'mag': lambda y: np.median(y), 'mag_err': lambda   x: np.sqrt(np.sum(x**2))/x.count()})
            tempztfz.insert(0,'name',sn)
            tempztfz.insert(4,'filter','z')
            z_dat.append(tempztfz)
    
        if len(g_dat)>0:
        	g_data=pd.concat(g_dat, ignore_index=True)
        else: g_data=pd.DataFrame()

        if len(r_dat)>0:
        	r_data=pd.concat(r_dat, ignore_index=True)
        else: r_data=pd.DataFrame()

        if len(i_dat)>0:
        	i_data=pd.concat(i_dat, ignore_index=True)
        else: i_data=pd.DataFrame()

        if len(z_dat)>0:
        	z_data=pd.concat(z_dat, ignore_index=True)
        else: z_data=pd.DataFrame()

        if len(g_colours[g_colours!=None])>0:
        
            if np.array([g_r_max,g_r_min]).any !=None:
                if g_r_max==None:
                    g_r_max=np.inf
                if g_r_min==None:
                    g_r_min=-np.inf
                    
                if not g_data.empty:
                    if not r_data.empty:
                            if r_data.mag[0]<maglim:
                                g_data.sort_values(by=['observed_date'],ascending=False)
                                r_data.sort_values(by=['observed_date'],ascending=False)
                                g_minus_r=g_data.mag[0]-r_data.mag[0]
                                g_minus_r_date_diff=abs(g_data.observed_date[0]-r_data.observed_date[0])
                                if g_minus_r >g_r_min:
                                    if g_minus_r <g_r_max:
                                        if g_data.name[0] not in (col_key_objects):
                                            col_key_objects.append(g_data.name[0])
                                            if abs(g_data.observed_date[0]-r_data.observed_date[0])>=6:
                                                col_flagged.append(g_data.name[0])
        
            if np.array([g_i_max,g_i_min]).any !=None:
                if g_i_max==None:
                    g_i_max=np.inf
                if g_i_min==None:
                   g_i_min=-np.inf

                if not g_data.empty:
                    if not i_data.empty:
                        if i_data.mag[0]<maglim:
                            g_data.sort_values(by=['observed_date'],ascending=False)
                            i_data.sort_values(by=['observed_date'],ascending=False)
                            g_minus_i=g_data.mag[0]-i_data.mag[0]
                            g_minus_i_date_diff=abs(g_data.observed_date[0]-i_data.observed_date[0])
                            if g_minus_i >g_i_min:
                                if g_minus_i <g_i_max:
                                    if g_data.name[0] not in (col_key_objects):
                                        col_key_objects.append(g_data.name[0])
                                        if abs(g_data.observed_date[0]-i_data.observed_date[0])>=6:
                                            col_flagged.append(g_data.name[0])


                if np.array([g_z_max,g_z_min]).any !=None:
                    if g_z_max==None:
                        g_z_max=np.inf
                    if g_z_min==None:
                       g_z_min=-np.inf

                if not g_data.empty:
                    if not z_data.empty:
                        if z_data.mag[0]<maglim:
                            g_data.sort_values(by=['observed_date'],ascending=False)
                            z_data.sort_values(by=['observed_date'],ascending=False)
                            g_minus_z=g_data.mag[0]-z_data.mag[0]
                            g_minus_z_date_diff=abs(g_data.observed_date[0]-z_data.observed_date[0])
                            if g_minus_z >g_z_min:
                                if g_minus_z <g_z_max:
                                    if g_data.name[0] not in (col_key_objects):
                                        col_key_objects.append(g_data.name[0])
                                        if abs(g_data.observed_date[0]-z_data.observed_date[0])>=6:
                                            col_flagged.append(g_data.name[0])
        if len(r_colours[r_colours!=None]):
        
            if np.array([r_i_max,r_i_min]).any !=None:
                if r_i_max==None:
                    r_i_max=np.inf
                if r_i_min==None:
                    r_i_min=-np.inf
                    
                if not r_data.empty:
                    if not i_data.empty:
                            if r_data.mag[0]<maglim:
                                r_data.sort_values(by=['observed_date'],ascending=False)
                                i_data.sort_values(by=['observed_date'],ascending=False)
                                r_minus_i=r_data.mag[0]-i_data.mag[0]
                                r_minus_i_date_diff=abs(r_data.observed_date[0]-i_data.observed_date[0])
                                if r_minus_i >r_i_min:
                                    if r_minus_i <r_i_max:
                                        if r_data.name[0] not in (col_key_objects):
                                            col_key_objects.append(r_data.name[0])
                                            if abs(r_data.observed_date[0]-i_data.observed_date[0])>=6:
                                                col_flagged.append(r_data.name[0])
        
            if np.array([r_z_max,r_z_min]).any !=None:
                if r_z_max==None:
                    r_z_max=np.inf
                if r_z_min==None:
                   r_z_min=-np.inf

                if not r_data.empty:
                    if not z_data.empty:
                        if r_data.mag[0]<maglim:
                            r_data.sort_values(by=['observed_date'],ascending=False)
                            z_data.sort_values(by=['observed_date'],ascending=False)
                            r_minus_z=r_data.mag[0]-z_data.mag[0]
                            r_minus_z_date_diff=abs(r_data.observed_date[0]-z_data.observed_date[0])
                            if r_minus_z >r_z_min:
                                if r_minus_z <r_z_max:
                                    if r_data.name[0] not in (col_key_objects):
                                        col_key_objects.append(r_data.name[0])
                                        if abs(r_data.observed_date[0]-z_data.observed_date[0])>=6:
                                            col_flagged.append(r_data.name[0])


        if len(i_colours[i_colours!=None]):
        
            if any([i_z_max,i_z_min]) !=None:
                if i_z_max==None:
                    i_z_max=np.inf
                if i_z_min==None:
                    i_z_min=-np.inf
                    
                if not i_data.empty:
                    if not z_data.empty:
                            if i_data.mag[0]<maglim:
                                i_data.sort_values(by=['observed_date'],ascending=False)
                                z_data.sort_values(by=['observed_date'],ascending=False)
                                i_minus_z=i_data.mag[0]-z_data.mag[0]
                                i_minus_z_date_diff=abs(i_data.observed_date[0]-z_data.observed_date[0])
                                if i_minus_z >i_z_min:
                                    if i_minus_z <i_z_max:
                                        if i_data.name[0] not in (col_key_objects):
                                            col_key_objects.append(i_data.name[0])
                                            if abs(i_data.observed_date[0]-z_data.observed_date[0])>=6:
                                                col_flagged.append(i_data.name[0])                          

if len(rises[rises!= None])>0:
    print('Checking rise times')
    for sn in candidates.groups.keys():
        g_dat=list()
        r_dat=list()
        i_dat=list()
        z_dat=list()
        tempdata=candidates.get_group(sn).reset_index()
    
        phot=tempdata.groupby('filter')

        if 'g' in phot.groups.keys():
            g_da=phot.get_group('g')
            tempg=g_da.groupby(np.round(g_da['observed_date']),as_index=False).agg({'ra':lambda r: np.median(r),'dec':lambda d: np.median(d),'observed_date':lambda z: np.median(z),'mag': lambda y: np.median(y), 'mag_err': lambda   x: np.sqrt(np.sum(x**2))/x.count()})
            tempg.insert(0,'name',sn)
            tempg.insert(4,'filter','g')
            g_dat.append(tempg)
        if 'g-ZTF' in phot.groups.keys():
            ztfg_da=phot.get_group('g-ZTF')
            tempztfg=ztfg_da.groupby(np.round(ztfg_da['observed_date']),as_index=False).agg({'ra':lambda r: np.median(r),'dec':lambda d: np.median(d),'observed_date':lambda z: np.median(z),'mag': lambda y: np.median(y), 'mag_err': lambda   x: np.sqrt(np.sum(x**2))/x.count()})
            tempztfg.insert(0,'name',sn)
            tempztfg.insert(4,'filter','g')
            g_dat.append(tempztfg)
        if 'r' in phot.groups.keys():
            r_da=phot.get_group('r')
            tempr=r_da.groupby(np.round(r_da['observed_date']),as_index=False).agg({'ra':lambda r: np.median(r),'dec':lambda d: np.median(d),'observed_date':lambda z: np.median(z),'mag': lambda y: np.median(y), 'mag_err': lambda   x: np.sqrt(np.sum(x**2))/x.count()})
            tempr.insert(0,'name',sn)
            tempr.insert(4,'filter','r')
            r_dat.append(tempr)
        if 'r-ZTF' in phot.groups.keys():
            ztfr_da=phot.get_group('r-ZTF')
            tempztfr=ztfr_da.groupby(np.round(ztfr_da['observed_date']),as_index=False).agg({'ra':lambda r: np.median(r),'dec':lambda d: np.median(d),'observed_date':lambda z: np.median(z),'mag': lambda y: np.median(y), 'mag_err': lambda   x: np.sqrt(np.sum(x**2))/x.count()})
            tempztfr.insert(0,'name',sn)
            tempztfr.insert(4,'filter','r')
            r_dat.append(tempztfr)
        if 'i' in phot.groups.keys():
            i_da=phot.get_group('i')
            tempi=i_da.groupby(np.round(i_da['observed_date']),as_index=False).agg({'ra':lambda r: np.median(r),'dec':lambda d: np.median(d),'observed_date':lambda z: np.median(z),'mag': lambda y: np.median(y), 'mag_err': lambda   x: np.sqrt(np.sum(x**2))/x.count()})
            tempi.insert(0,'name',sn)
            tempi.insert(4,'filter','i')
            i_dat.append(tempi)
        if 'i-ZTF' in phot.groups.keys():
            ztfi_da=phot.get_group('i-ZTF')
            tempztfi=ztfi_da.groupby(np.round(ztfi_da['observed_date']),as_index=False).agg({'ra':lambda r: np.median(r),'dec':lambda d: np.median(d),'observed_date':lambda z: np.median(z),'mag': lambda y: np.median(y), 'mag_err': lambda   x: np.sqrt(np.sum(x**2))/x.count()})
            tempztfi.insert(0,'name',sn)
            tempztfi.insert(4,'filter','i')
            i_dat.append(tempztfi)
        if 'z' in phot.groups.keys():
            z_da=phot.get_group('z')
            tempz=z_da.groupby(np.round(z_da['observed_date']),as_index=False).agg({'ra':lambda r: np.median(r),'dec':lambda d: np.median(d),'observed_date':lambda z: np.median(z),'mag': lambda y: np.median(y), 'mag_err': lambda   x: np.sqrt(np.sum(x**2))/x.count()})
            tempz.insert(0,'name',sn)
            tempz.insert(4,'filter','z')
            z_dat.append(tempz)
        if 'z-ZTF' in phot.groups.keys():
            ztfz_da=phot.get_group('z-ZTF')
            tempztfz=ztfz_da.groupby(np.round(ztfz_da['observed_date']),as_index=False).agg({'ra':lambda r: np.median(r),'dec':lambda d: np.median(d),'observed_date':lambda z: np.median(z),'mag': lambda y: np.median(y), 'mag_err': lambda   x: np.sqrt(np.sum(x**2))/x.count()})
            tempztfz.insert(0,'name',sn)
            tempztfz.insert(4,'filter','z')
            z_dat.append(tempztfz)
    
        if len(g_dat)>0:
            g_data=pd.concat(g_dat, ignore_index=True)
        else: g_data=pd.DataFrame()

        if len(r_dat)>0:
            r_data=pd.concat(r_dat, ignore_index=True)
        else: r_data=pd.DataFrame()

        if len(i_dat)>0:
            i_data=pd.concat(i_dat, ignore_index=True)
        else: i_data=pd.DataFrame()

        if len(z_dat)>0:
            z_data=pd.concat(z_dat, ignore_index=True)
        else: z_data=pd.DataFrame()

        if risemin==None:
            risemin=np.inf
        if risemax==None:
            risemax=np.inf
        
        if riseband=='g':
            if not g_data.empty:
                if g_data.mag.count()>1:
                    riseg=g_data.sort_values(by=['observed_date'],ascending=False).reset_index()
                    if (np.round(riseg.observed_date[0])-np.round(riseg.observed_date[1]))>=1:
                        if (np.round(riseg.observed_date[0])-np.round(riseg.observed_date[1]))<=6:
                            g_rise=np.diff(riseg.mag)/np.diff(riseg.observed_date)
                            if np.isnan(g_rise[0]) == False:
                                if g_rise[0]>= risemin:
                                    if g_rise[0]<= risemax:
                                        rise_key_objects.append(sn)
                                        if len(g_rise)==1:
                                            rise_flagged.append(sn)

        if riseband=='r':
            if not r_data.empty:
                if r_data.mag.count()>1:
                    riser=r_data.sort_values(by=['observed_date'],ascending=False).reset_index()
                    if (np.round(riser.observed_date[0])-np.round(riser.observed_date[1]))>=1:
                        if (np.round(riser.observed_date[0])-np.round(riser.observed_date[1]))<=6:
                            r_rise=np.diff(riser.mag)/np.diff(riser.observed_date)
                            if np.isnan(r_rise[0]) == False:
                                if r_rise[0]>= risemin:
                                    if r_rise[0]<= risemax:
                                        rise_key_objects.append(sn)
                                        if len(r_rise)==1:
                                            rise_flagged.append(sn)

        if riseband=='i':
            if not i_data.empty:
                if i_data.mag.count()>1:
                    risei=g_data.sort_values(by=['observed_date'],ascending=False).reset_index()
                    if (np.round(risei.observed_date[0])-np.round(risei.observed_date[1]))>=1:
                        if (np.round(risei.observed_date[0])-np.round(risei.observed_date[1]))<=6:
                            i_rise=np.diff(risei.mag)/np.diff(risei.observed_date)
                            if np.isnan(i_rise[0]) == False:
                                if i_rise[0]>= risemin:
                                    if i_rise[0]<= risemax:
                                        rise_key_objects.append(sn)
                                        if len(i_rise)==1:
                                            rise_flagged.append(sn)

        if riseband=='z':
            if not z_data.empty:
                if z_data.mag.count()>1:
                    risez=z_data.sort_values(by=['observed_date'],ascending=False).reset_index()
                    if (np.round(risez.observed_date[0])-np.round(risez.observed_date[1]))>=1:
                        if (np.round(risez.observed_date[0])-np.round(risez.observed_date[1]))<=6:
                            z_rise=np.diff(risez.mag)/np.diff(risez.observed_date)
                            if np.isnan(z_rise[0]) == False:
                                if z_rise[0]>= risemin:
                                    if z_rise[0]<= risemax:
                                        rise_key_objects.append(sn)
                                        if len(z_rise)==1:
                                            rise_flagged.append(sn)

if hostenv != None:
    print('Checking host environments')
    for sn in candidates.groups.keys():       
        g_dat=list()
        r_dat=list()
        i_dat=list()
        z_dat=list()

        tempdata=candidates.get_group(sn).reset_index()
        phot=tempdata.groupby('filter')
  
        if hostmagdiff!=None:

            if 'g' in phot.groups.keys():
                g_da=phot.get_group('g')
                tempg=g_da.groupby(np.round(g_da['observed_date']),as_index=False).agg({'ra':lambda r: np.median(r),'dec':lambda d: np.median(d),'observed_date':lambda z: np.median(z),'mag': lambda y: np.median(y), 'mag_err': lambda   x: np.sqrt(np.sum(x**2))/x.count()})
                tempg.insert(0,'name',sn)
                tempg.insert(4,'filter','g')
                g_dat.append(tempg)
            if 'g-ZTF' in phot.groups.keys():
                ztfg_da=phot.get_group('g-ZTF')
                tempztfg=ztfg_da.groupby(np.round(ztfg_da['observed_date']),as_index=False).agg({'ra':lambda r: np.median(r),'dec':lambda d: np.median(d),'observed_date':lambda z: np.median(z),'mag': lambda y: np.median(y), 'mag_err': lambda   x: np.sqrt(np.sum(x**2))/x.count()})
                tempztfg.insert(0,'name',sn)
                tempztfg.insert(4,'filter','g')
                g_dat.append(tempztfg)
            if 'r' in phot.groups.keys():
                r_da=phot.get_group('r')
                tempr=r_da.groupby(np.round(r_da['observed_date']),as_index=False).agg({'ra':lambda r: np.median(r),'dec':lambda d: np.median(d),'observed_date':lambda z: np.median(z),'mag': lambda y: np.median(y), 'mag_err': lambda   x: np.sqrt(np.sum(x**2))/x.count()})
                tempr.insert(0,'name',sn)
                tempr.insert(4,'filter','r')
                r_dat.append(tempr)
            if 'r-ZTF' in phot.groups.keys():
                ztfr_da=phot.get_group('r-ZTF')
                tempztfr=r_da.groupby(np.round(ztfr_da['observed_date']),as_index=False).agg({'ra':lambda r: np.median(r),'dec':lambda d: np.median(d),'observed_date':lambda z: np.median(z),'mag': lambda y: np.median(y), 'mag_err': lambda   x: np.sqrt(np.sum(x**2))/x.count()})
                tempztfr.insert(0,'name',sn)
                tempztfr.insert(4,'filter','r')
                r_dat.append(tempztfr)
            if 'i' in phot.groups.keys():
                i_da=phot.get_group('i')
                tempi=r_da.groupby(np.round(i_da['observed_date']),as_index=False).agg({'ra':lambda r: np.median(r),'dec':lambda d: np.median(d),'observed_date':lambda z: np.median(z),'mag': lambda y: np.median(y), 'mag_err': lambda   x: np.sqrt(np.sum(x**2))/x.count()})
                tempi.insert(0,'name',sn)
                tempi.insert(4,'filter','i')
                i_dat.append(tempi)
            if 'i-ZTF' in phot.groups.keys():
                ztfi_da=phot.get_group('i-ZTF')
                tempztfi=r_da.groupby(np.round(ztfi_da['observed_date']),as_index=False).agg({'ra':lambda r: np.median(r),'dec':lambda d: np.median(d),'observed_date':lambda z: np.median(z),'mag': lambda y: np.median(y), 'mag_err': lambda   x: np.sqrt(np.sum(x**2))/x.count()})
                tempztfi.insert(0,'name',sn)
                tempztfi.insert(4,'filter','i')
                i_dat.append(tempztfi)
            if 'z' in phot.groups.keys():
                z_da=phot.get_group('z')
                tempz=r_da.groupby(np.round(z_da['observed_date']),as_index=False).agg({'ra':lambda r: np.median(r),'dec':lambda d: np.median(d),'observed_date':lambda z: np.median(z),'mag': lambda y: np.median(y), 'mag_err': lambda   x: np.sqrt(np.sum(x**2))/x.count()})
                tempz.insert(0,'name',sn)
                tempz.insert(4,'filter','z')
                z_dat.append(tempz)
            if 'z-ZTF' in phot.groups.keys():
                ztfz_da=phot.get_group('z-ZTF')
                tempztfz=r_da.groupby(np.round(ztfz_da['observed_date']),as_index=False).agg({'ra':lambda r: np.median(r),'dec':lambda d: np.median(d),'observed_date':lambda z: np.median(z),'mag': lambda y: np.median(y), 'mag_err': lambda   x: np.sqrt(np.sum(x**2))/x.count()})
                tempztfz.insert(0,'name',sn)
                tempztfz.insert(4,'filter','z')
                z_dat.append(tempztfz)
        
            if len(g_dat)>0:
                g_data=pd.concat(g_dat, ignore_index=True)
            else: g_data=pd.DataFrame()

            if len(r_dat)>0:
                r_data=pd.concat(r_dat, ignore_index=True)
            else: r_data=pd.DataFrame()

            if len(i_dat)>0:
                i_data=pd.concat(i_dat, ignore_index=True)
            else: i_data=pd.DataFrame()

            if len(z_dat)>0:
                z_data=pd.concat(z_dat, ignore_index=True)
            else: z_data=pd.DataFrame()
        
        successful=False
        try:
            nearby=PS1_query.search(sn,tempdata.ra.mean(),tempdata.dec.mean(),1).to_pandas().sort_values('Seperation').reset_index()
            successful=True
        except:
            continue
            
        if successful==False:
            no_host_flag.append(sn)
        if successful==True:
            if hostenv=='host':
                if nearby.iloc[0]['Seperation']>0.1:
                    if nearby.iloc[0]['Seperation']<=2:
                        if hostmagdiff !=None:
                            if hostband=='g':
                                maghostg=nearby.iloc[0]['g']
                                checkg=g_data.sort_values(by=['observed_date'],ascending=False).reset_index()
                                if maghostg-checkg.mag.values[0]>=hostmagdiff:
                                    host_key_objects.append(sn)
                                if maghostg-checkg.mag.values[0]<hostmagdiff:
                                    host_flagged.append(sn)
                            if hostband=='r':
                                maghostr=nearby.iloc[0]['r']
                                checkr=r_data.sort_values(by=['observed_date'],ascending=False).reset_index()
                                if maghostr-checkr.mag.values[0]>=hostmagdiff:
                                    host_key_objects.append(sn)
                                if maghostr-checkr.mag.values[0]<hostmagdiff:
                                    host_flagged.append(sn)
                            if hostband=='r':
                                maghosti=nearby.iloc[0]['i']
                                checki=i_data.sort_values(by=['observed_date'],ascending=False).reset_index()
                                if maghosti-checki.mag.values[0]>=hostmagdiff:
                                    host_key_objects.append(sn)
                                if maghosti-checki.mag.values[0]<hostmagdiff:
                                    host_flagged.append(sn)
                            if hostband=='z':
                                maghostz=nearby.iloc[0]['z']
                                checkz=z_data.sort_values(by=['observed_date'],ascending=False).reset_index()
                                if maghostz-checkz.mag.values[0]>=hostmagdiff:
                                    host_key_objects.append(sn)
                                if maghostz-checkz.mag.values[0]<hostmagdiff:
                                    host_flagged.append(sn)
                        if hostmagdiff ==None:
                            host_key_objects.append(sn)
            if hostenv=='nuclear':

                if nearby.iloc[0]['Seperation']<=0.1:
                    if hostmagdiff !=None:
                        if hostband=='g':
                            maghostg=nearby.iloc[0]['g']
                            checkg=g_data.sort_values(by=['observed_date'],ascending=False).reset_index()
                            if maghostg-checkg.mag.values[0]>=hostmagdiff:
                                nuclear.append(sn)
                            if maghostg-checkg.mag.values[0]<hostmagdiff:
                                nuclear_flagged.append(sn)
                        if hostband=='r':
                            maghostr=nearby.iloc[0]['r']
                            checkr=r_data.sort_values(by=['observed_date'],ascending=False).reset_index()
                            if maghostr-checkr.mag.values[0]>=hostmagdiff:
                                nuclear.append(sn)
                            if maghostr-checkr.mag.values[0]<hostmagdiff:
                                nuclear_flagged.append(sn)
                        if hostband=='r':
                            maghosti=nearby.iloc[0]['i']
                            checki=i_data.sort_values(by=['observed_date'],ascending=False).reset_index()
                            if maghosti-checki.mag.values[0]>=hostmagdiff:
                                nuclear.append(sn)
                            if maghosti-checki.mag.values[0]<hostmagdiff:
                                nuclear_flagged.append(sn)
                        if hostband=='z':
                            maghostz=nearby.iloc[0]['z']
                            checkz=z_data.sort_values(by=['observed_date'],ascending=False).reset_index()
                            if maghostz-checkz.mag.values[0]>=hostmagdiff:
                                nuclear.append(sn)
                            if maghostz-checkz.mag.values[0]<hostmagdiff:
                                nuclear_flagged.append(sn)
                    if hostmagdiff ==None:
                        nuclear.append(sn)
            if hostenv=='orphan':
                if nearby.iloc[0]['Seperation']>2:
                    hostless.append(sn)

if len(distances[distances!= None])>0:
    print('Checking photo-z and absolute mags')
    for sn in candidates.groups.keys():
        g_dat=list()
        r_dat=list()
        i_dat=list()
        z_dat=list()
        tempdata=candidates.get_group(sn).reset_index()
        tempdist=tempdata.photo_z.mean()

        if photoz!=None:
           if tempdist<photoz+photozrange:
                if  tempdist>photoz-photozrange:
                    redshift_range.append(sn)
                    dist_key_objects.append(sn)

    
        if absmag!=None:
            if tempdist != np.nan:
                phot=tempdata.groupby('filter')

                if 'g' in phot.groups.keys():
                    g_da=phot.get_group('g')
                    tempg=g_da.groupby(np.round(g_da['observed_date']),as_index=False).agg({'ra':lambda r: np.median(r),'dec':lambda d: np.median(d),'observed_date':lambda z: np.median(z),'mag': lambda y: np.median(y), 'mag_err': lambda   x: np.sqrt(np.sum(x**2))/x.count()})
                    tempg.insert(0,'name',sn)
                    tempg.insert(4,'filter','g')
                    g_dat.append(tempg)
                if 'g-ZTF' in phot.groups.keys():
                    ztfg_da=phot.get_group('g-ZTF')
                    tempztfg=ztfg_da.groupby(np.round(ztfg_da['observed_date']),as_index=False).agg({'ra':lambda r: np.median(r),'dec':lambda d: np.median(d),'observed_date':lambda z: np.median(z),'mag': lambda y: np.median(y), 'mag_err': lambda   x: np.sqrt(np.sum(x**2))/x.count()})
                    tempztfg.insert(0,'name',sn)
                    tempztfg.insert(4,'filter','g')
                    g_dat.append(tempztfg)
                if 'r' in phot.groups.keys():
                    r_da=phot.get_group('r')
                    tempr=r_da.groupby(np.round(r_da['observed_date']),as_index=False).agg({'ra':lambda r: np.median(r),'dec':lambda d: np.median(d),'observed_date':lambda z: np.median(z),'mag': lambda y: np.median(y), 'mag_err': lambda   x: np.sqrt(np.sum(x**2))/x.count()})
                    tempr.insert(0,'name',sn)
                    tempr.insert(4,'filter','r')
                    r_dat.append(tempr)
                if 'r-ZTF' in phot.groups.keys():
                    ztfr_da=phot.get_group('r-ZTF')
                    tempztfr=r_da.groupby(np.round(ztfr_da['observed_date']),as_index=False).agg({'ra':lambda r: np.median(r),'dec':lambda d: np.median(d),'observed_date':lambda z: np.median(z),'mag': lambda y: np.median(y), 'mag_err': lambda   x: np.sqrt(np.sum(x**2))/x.count()})
                    tempztfr.insert(0,'name',sn)
                    tempztfr.insert(4,'filter','r')
                    r_dat.append(tempztfr)
                if 'i' in phot.groups.keys():
                    i_da=phot.get_group('i')
                    tempi=r_da.groupby(np.round(i_da['observed_date']),as_index=False).agg({'ra':lambda r: np.median(r),'dec':lambda d: np.median(d),'observed_date':lambda z: np.median(z),'mag': lambda y: np.median(y), 'mag_err': lambda   x: np.sqrt(np.sum(x**2))/x.count()})
                    tempi.insert(0,'name',sn)
                    tempi.insert(4,'filter','i')
                    i_dat.append(tempi)
                if 'i-ZTF' in phot.groups.keys():
                    ztfi_da=phot.get_group('i-ZTF')
                    tempztfi=r_da.groupby(np.round(ztfi_da['observed_date']),as_index=False).agg({'ra':lambda r: np.median(r),'dec':lambda d: np.median(d),'observed_date':lambda z: np.median(z),'mag': lambda y: np.median(y), 'mag_err': lambda   x: np.sqrt(np.sum(x**2))/x.count()})
                    tempztfi.insert(0,'name',sn)
                    tempztfi.insert(4,'filter','i')
                    i_dat.append(tempztfi)
                if 'z' in phot.groups.keys():
                    z_da=phot.get_group('z')
                    tempz=r_da.groupby(np.round(z_da['observed_date']),as_index=False).agg({'ra':lambda r: np.median(r),'dec':lambda d: np.median(d),'observed_date':lambda z: np.median(z),'mag': lambda y: np.median(y), 'mag_err': lambda   x: np.sqrt(np.sum(x**2))/x.count()})
                    tempz.insert(0,'name',sn)
                    tempz.insert(4,'filter','z')
                    z_dat.append(tempz)
                if 'z-ZTF' in phot.groups.keys():
                    ztfz_da=phot.get_group('z-ZTF')
                    tempztfz=r_da.groupby(np.round(ztfz_da['observed_date']),as_index=False).agg({'ra':lambda r: np.median(r),'dec':lambda d: np.median(d),'observed_date':lambda z: np.median(z),'mag': lambda y: np.median(y), 'mag_err': lambda   x: np.sqrt(np.sum(x**2))/x.count()})
                    tempztfz.insert(0,'name',sn)
                    tempztfz.insert(4,'filter','z')
                    z_dat.append(tempztfz)
            
                if len(g_dat)>0:
                    g_data=pd.concat(g_dat, ignore_index=True)
                else: g_data=pd.DataFrame()

                if len(r_dat)>0:
                    r_data=pd.concat(r_dat, ignore_index=True)
                else: r_data=pd.DataFrame()

                if len(i_dat)>0:
                    i_data=pd.concat(i_dat, ignore_index=True)
                else: i_data=pd.DataFrame()

                if len(z_dat)>0:
                    z_data=pd.concat(z_dat, ignore_index=True)
                else: z_data=pd.DataFrame()


                if absmagband=='g':
                    if not g_data.empty:
                        maxg=np.min(g_data.mag.values)
                        lumdist=cosmo.luminosity_distance(tempdist).value
                        gabs=maxg-(5*np.log10((lumdist*(10**5))))
                        if absmaglowlim != 'True':
                            if gabs<=absmag-absmagrange:
                               if gabs>=absmag+absmagrange:
                                    bright_range.append(sn)
                                    dist_key_objects.append(sn) 
                        elif absmaglowlim =='True':
                            if gabs<=absmag:
                                bright.append(sn)
                                dist_key_objects.append(sn)

                if absmagband=='r':
                    if not r_data.empty:
                        maxr=np.min(r_data.mag.values)
                        lumdist=cosmo.luminosity_distance(tempdist).value
                        rabs=maxr-(5*np.log10((lumdist*(10**5))))
                        if absmaglowlim != 'True':
                            if rabs<=absmag-absmagrange:
                               if rabs>=absmag+absmagrange:
                                    bright_range.append(sn)
                                    dist_key_objects.append(sn) 
                        elif absmaglowlim =='True':
                            if rabs<=absmag:
                                bright.append(sn)
                                dist_key_objects.append(sn)

                if absmagband=='i':
                    if not i_data.empty:
                        maxi=np.min(i_data.mag.values)
                        lumdist=cosmo.luminosity_distance(tempdist).value
                        iabs=maxi-(5*np.log10((lumdist*(10**5))))
                        if absmaglowlim != 'True':
                            if iabs<=absmag-absmagrange:
                               if iabs>=absmag+absmagrange:
                                    bright_range.append(sn)
                                    dist_key_objects.append(sn) 
                        elif absmaglowlim =='True':
                            if iabs<=absmag:
                                bright.append(sn)
                                dist_key_objects.append(sn)

                if absmagband=='z':
                    if not z_data.empty:
                        maxz=np.min(z_data.mag.values)
                        lumdist=cosmo.luminosity_distance(tempdist).value
                        zabs=maxz-(5*np.log10((lumdist*(10**5))))
                        if absmaglowlim != 'True':
                            if zabs<=absmag-absmagrange:
                               if zabs>=absmag+absmagrange:
                                    bright_range.append(sn)
                                    dist_key_objects.append(sn) 
                        elif absmaglowlim =='True':
                            if zabs<=absmag:
                                bright.append(sn)
                                dist_key_objects.append(sn)
 
# #Ugly way to sort output lists  
if len(all_colours[all_colours!= None])>0:
    if len(rises[rises!= None])>0:
        if hostenv != None:
            if len(distances[distances!= None])>0:
                s1 = set(col_key_objects) 
                s2 = set(rise_key_objects) 
                s3 = set(host_key_objects)
                s4 = set(dist_key_objects)  
                set1 = s1.intersection(s2)
                next_set = set1.intersection(s3)
                result_set = next_set.intersection(s4)  
                gold = list(result_set) 
                ss1=s1.difference(result_set)
                ss2=s2.difference(result_set)
                ss3=s3.difference(result_set)
                ss4=s4.difference(result_set)
                silver=list(ss1)+list(ss2)+list(ss3)+list(ss4)

if len(all_colours[all_colours!= None])>0:
    if len(rises[rises!= None])>0:
        if hostenv != None:
            if len(distances[distances!= None])==0:
                s1 = set(col_key_objects) 
                s2 = set(rise_key_objects) 
                s3 = set(host_key_objects)
                set1 = s1.intersection(s2)
                result_set = set1.intersection(s3)  
                gold = list(result_set) 
                ss1=s1.difference(result_set)
                ss2=s2.difference(result_set)
                ss3=s3.difference(result_set)
                silver=list(ss1)+list(ss2)+list(ss3)

if len(all_colours[all_colours!= None])==0:
    if len(rises[rises!= None])>0:
        if hostenv != None:
            if len(distances[distances!= None])>0:
                s1 = set(dist_key_objects) 
                s2 = set(rise_key_objects) 
                s3 = set(host_key_objects)
                set1 = s1.intersection(s2)
                result_set = set1.intersection(s3)  
                gold = list(result_set) 
                ss1=s1.difference(result_set)
                ss2=s2.difference(result_set)
                ss3=s3.difference(result_set)
                silver=list(ss1)+list(ss2)+list(ss3)    

if len(all_colours[all_colours!= None])>0:
    if len(rises[rises!= None])==0:
        if hostenv != None:
            if len(distances[distances!= None])>0:
                s1 = set(col_key_objects) 
                s2 = set(dist_key_objects) 
                s3 = set(host_key_objects)
                set1 = s1.intersection(s2)
                result_set = set1.intersection(s3)  
                gold = list(result_set) 
                ss1=s1.difference(result_set)
                ss2=s2.difference(result_set)
                ss3=s3.difference(result_set)
                silver=list(ss1)+list(ss2)+list(ss3)    

if len(all_colours[all_colours!= None])>0:
    if len(rises[rises!= None])>0:
        if hostenv == None:
            if len(distances[distances!= None])>0:
                s1 = set(col_key_objects) 
                s2 = set(rise_key_objects) 
                s3 = set(dist_key_objects)
                set1 = s1.intersection(s2)
                result_set = set1.intersection(s3)  
                gold = list(result_set) 
                ss1=s1.difference(result_set)
                ss2=s2.difference(result_set)
                ss3=s3.difference(result_set)
                silver=list(ss1)+list(ss2)+list(ss3)    

if len(all_colours[all_colours!= None])>0:
    if len(rises[rises!= None])>0:
        if hostenv == None:
            if len(distances[distances!= None])==0:
                s1 = set(col_key_objects) 
                s2 = set(rise_key_objects) 
                result_set = s1.intersection(s2)
                gold = list(result_set)
                ss1=s1.difference(result_set)
                ss2=s2.difference(result_set)
                silver=list(ss1)+list(ss2)  

if len(all_colours[all_colours!= None])>0:
    if len(rises[rises!= None])==0:
        if hostenv != None:
            if len(distances[distances!= None])==0:
                s1 = set(col_key_objects) 
                s2 = set(host_key_objects) 
                result_set = s1.intersection(s2)
                gold = list(result_set)
                ss1=s1.difference(result_set)
                ss2=s2.difference(result_set)
                silver=list(ss1)+list(ss2) 

if len(all_colours[all_colours!= None])>0:
    if len(rises[rises!= None])==0:
        if hostenv == None:
            if len(distances[distances!= None])>0:
                s1 = set(col_key_objects) 
                s2 = set(dist_key_objects) 
                result_set = s1.intersection(s2)
                gold = list(result_set)
                ss1=s1.difference(result_set)
                ss2=s2.difference(result_set)
                silver=list(ss1)+list(ss2) 

if len(all_colours[all_colours!= None])==0:
    if len(rises[rises!= None])>0:
        if hostenv != None:
            if len(distances[distances!= None])==0:
                s1 = set(host_key_objects) 
                s2 = set(rise_key_objects) 
                result_set = s1.intersection(s2)
                gold = list(result_set)
                ss1=s1.difference(result_set)
                ss2=s2.difference(result_set)
                silver=list(ss1)+list(ss2) 

if len(all_colours[all_colours!= None])==0:
    if len(rises[rises!= None])>0:
        if hostenv == None:
            if len(distances[distances!= None])>0:
                s1 = set(dist_key_objects) 
                s2 = set(rise_key_objects) 
                result_set = s1.intersection(s2)
                gold = list(result_set)
                ss1=s1.difference(result_set)
                ss2=s2.difference(result_set)
                silver=list(ss1)+list(ss2) 

if len(all_colours[all_colours!= None])==0:
    if len(rises[rises!= None])==0:
        if hostenv != None:
            if len(distances[distances!= None])>0:
                s1 = set(dist_key_objects) 
                s2 = set(host_key_objects) 
                result_set = s1.intersection(s2)
                gold = list(result_set)
                ss1=s1.difference(result_set)
                ss2=s2.difference(result_set)
                silver=list(ss1)+list(ss2) 

if len(all_colours[all_colours!= None])==0:
    if len(rises[rises!= None])==0:
        if len(distances[distances!= None])==0:
            if hostenv != None:
                if hostenv=='host':
                    gold=host_key_objects
                    silver=host_flagged
                if hostenv=='orphan':
                    gold=hostless
                    silver=host_flagged
                if hostenv=='nuclear':
                    gold=nuclear
                    silver=host_flagged

if len(all_colours[all_colours!= None])>0:
    if len(rises[rises!= None])==0:
        if hostenv == None:
            if len(distances[distances!= None])==0:
                gold=col_key_objects
                silver=col_flagged

if len(all_colours[all_colours!= None])==0:
    if len(rises[rises!= None])>0:
        if hostenv == None:
            if len(distances[distances!= None])==0:
                gold=rise_key_objects
                silver=rise_flagged

if len(all_colours[all_colours!= None])==0:
    if len(rises[rises!= None])==0:
        if hostenv == None:
            if len(distances[distances!= None])>0:
                gold=dist_key_objects


print('Total objects meeting all selection criteria (%i/%i):'%(len(gold),(len(candidates.groups.keys()))))

for i in gold:
	print(i)

#print marginal candidate matches if requested
if showall=='True':
    print('Candidates which partially meet the criteria (%i):'%(len(silver)))
    for i in silver:
        print(i)

    if hostmagdiff!=None:
        print('Candidates which do not meet the host mag difference: (%i):'%(len(host_flagged)))
        for i in host_flagged:
            print(i)

    if len(all_colours[all_colours!= None])==0:
        print ('Flagged candidates with large gap between photometry: %i):'%(len(col_flagged)))
        for i in col_flagged:
            print(i)


out_file = open("GEYSERS_output.txt", "w")
out_file.write('Interesting candidates meeting all serach criteria:\n')
for obj in gold:
    out_file.write('https://ziggy.ucolick.org/yse/transient_detail/'+obj+'/\n')
out_file.write('\n')
out_file.write('Marginal match candidates:\n')
for obj in silver:
     out_file.write('https://ziggy.ucolick.org/yse/transient_detail/'+obj+'/\n')

