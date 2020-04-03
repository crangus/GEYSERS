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



#weighting for grouped photometry
def wavg(Flux,Err):
    d = Flux
    w = 1/Err
    try:
        return (d * w).sum() / w.sum()
    except ZeroDivisionError:
        return d.mean()

parser = argparse.ArgumentParser(description='Transient properties:')
#g colour args
parser.add_argument('file', type=str, help='input csv file',default=None)

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


#host args
parser.add_argument('--host', type=str , help='host environment ("host"/"nuclear"/"orphan")',default=None)

parser.add_argument('--maglim', type=float, help='faintest objects of interest (all bands). Default is 20.0',default=20.0)

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
maglim=condits.maglim



photometry=pd.read_csv(file,sep=',',header=0)

#Convert date to MJD
photometry.observed_date=Time([x for x in photometry.observed_date.values]).mjd

#Drop rows with gaia or "unknown" photometric bands
indexGaia = photometry[ photometry['filter'] == 'G' ].index
indexUnknown = photometry[ photometry['filter'] == 'Unknown' ].index
photometry.drop(['count_date'], axis=1,inplace=True)
photometry.drop(indexGaia , inplace=True)
photometry.drop(indexUnknown , inplace=True)

hosted=[]
hostless=[]
nuclear=[]
key_objects=[]
flagged=[]
hostflag=[]
rising=[]
rate=[]

g_colours=np.array([g_r_min,g_r_max,g_i_min,g_i_max,g_z_min,g_z_max])
r_colours=np.array([r_i_min,r_i_max,r_z_min,r_z_max])
i_colours=np.array([i_z_min,i_z_max])
all_colours=np.array([g_r_min,g_r_max,g_i_min,g_i_max,g_z_min,g_z_max,r_i_min,r_i_max,r_z_min,r_z_max,i_z_min,i_z_max])

rises=np.array([risemin,risemax])

if len(rises[rises!=None])>0 and riseband==None:
    print()
    print('Error: You need to enter the rise/delcine photometric band (g/r/i/z)')
    print()
    exit()


candidates=photometry.groupby('name')

print("Ingesting %i candidates detected within the last 20 days"%(len(candidates.groups.keys())))

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
                                        if g_data.name[0] not in (key_objects):
                                            key_objects.append(g_data.name[0])
                                            if abs(g_data.observed_date[0]-r_data.observed_date[0])>=6:
                                                flagged.append(g_data.name[0])
        
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
                                    if g_data.name[0] not in (key_objects):
                                        key_objects.append(g_data.name[0])
                                        if abs(g_data.observed_date[0]-i_data.observed_date[0])>=6:
                                            flagged.append(g_data.name[0])


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
                                    if g_data.name[0] not in (key_objects):
                                        key_objects.append(g_data.name[0])
                                        if abs(g_data.observed_date[0]-z_data.observed_date[0])>=6:
                                            flagged.append(g_data.name[0])
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
                                        if r_data.name[0] not in (key_objects):
                                            key_objects.append(r_data.name[0])
                                            if abs(r_data.observed_date[0]-i_data.observed_date[0])>=6:
                                                flagged.append(r_data.name[0])
        
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
                                    if r_data.name[0] not in (key_objects):
                                        key_objects.append(r_data.name[0])
                                        if abs(r_data.observed_date[0]-z_data.observed_date[0])>=6:
                                            flagged.append(r_data.name[0])


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
                                        if i_data.name[0] not in (key_objects):
                                            key_objects.append(i_data.name[0])
                                            if abs(i_data.observed_date[0]-z_data.observed_date[0])>=6:
                                                flagged.append(i_data.name[0])                          
if len(rises[rises!= None])==0:
    KEY_objects=key_objects
    FLAGGED=flagged        


if len(rises[rises!= None])>0:
    KEY_objects=[]
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


        if risemin==None:
            risemin=np.inf
        if risemax==None:
            risemax=np.inf
        
        if riseband=='g':
            if not g_data.empty:
                if g_data.mag.count()>1:
                    rise=g_data.sort_values(by=['observed_date'],ascending=False).reset_index()
                    if (np.round(rise.observed_date[0])-np.round(rise.observed_date[1]))>=1:
                        if (np.round(rise.observed_date[0])-np.round(rise.observed_date[1]))<=6:
                            g_rise=np.diff(rise.mag)/np.diff(rise.observed_date)
                            if np.isnan(g_rise[0]) == False:
                                if g_rise[0]>= risemin:
                                    if g_rise[0]<= risemax:
                                        KEY_objects.append(sn)
                                        if len(g_rise)==1:
                                            FLAGGED.append(sn)

        if riseband=='r':
            if not r_data.empty:
                if r_data.mag.count()>1:
                    rise=r_data.sort_values(by=['observed_date'],ascending=False).reset_index()
                    if (np.round(rise.observed_date[0])-np.round(rise.observed_date[1]))>=1:
                        if (np.round(rise.observed_date[0])-np.round(rise.observed_date[1]))<=6:
                            r_rise=np.diff(rise.mag)/np.diff(rise.observed_date)
                            if np.isnan(r_rise[0]) == False:
                                if r_rise[0]>= risemin:
                                    if r_rise[0]<= risemax:
                                        KEY_objects.append(sn)
                                        if len(r_rise)==1:
                                            FLAGGED.append(sn)

        if riseband=='i':
            if not i_data.empty:
                if i_data.mag.count()>1:
                    rise=g_data.sort_values(by=['observed_date'],ascending=False).reset_index()
                    if (np.round(rise.observed_date[0])-np.round(rise.observed_date[1]))>=1:
                        if (np.round(rise.observed_date[0])-np.round(rise.observed_date[1]))<=6:
                            i_rise=np.diff(rise.mag)/np.diff(rise.observed_date)
                            if np.isnan(i_rise[0]) == False:
                                if i_rise[0]>= risemin:
                                    if i_rise[0]<= risemax:
                                        KEY_objects.append(sn)
                                        if len(i_rise)==1:
                                            FLAGGED.append(sn)

        if riseband=='z':
            if not z_data.empty:
                if z_data.mag.count()>1:
                    rise=z_data.sort_values(by=['observed_date'],ascending=False).reset_index()
                    if (np.round(rise.observed_date[0])-np.round(rise.observed_date[1]))>=1:
                        if (np.round(rise.observed_date[0])-np.round(rise.observed_date[1]))<=6:
                            z_rise=np.diff(rise.mag)/np.diff(rise.observed_date)
                            if np.isnan(z_rise[0]) == False:
                                if z_rise[0]>= risemin:
                                    if z_rise[0]<= risemax:
                                        KEY_objects.append(sn)
                                        if len(z_rise)==1:
                                            FLAGGED.append(sn)

#Chicking host environments
if hostenv == None:
    promoted=KEY_objects
    final_FLAG=FLAGGED


if hostenv != None:
    final_FLAG=[]
    promoted=[]
    print('Checking host environments')
    if len(key_objects)==0:
        for sn in candidates.groups.keys():
            tempdata=candidates.get_group(sn).reset_index()
            successful=False
            try:
                nearby=PS1_query.search(sn,tempdata.ra.mean(),tempdata.dec.mean(),1).to_pandas().sort_values('Seperation').reset_index()
                successful=True
            except:
                continue
                
            if successful==False:
                hostflag.append(sn)
            if successful==True:
                if hostenv=='host':
                    if nearby.iloc[0]['Seperation']>0.1:
                        if nearby.iloc[0]['Seperation']<=2:
                           promoted.append(sn)
                    elif nearby.iloc[0]['Seperation']<=0.1:
                        nuclear.append(sn)
                    elif nearby.iloc[0]['Seperation']>2:
                        hostless.append(sn)
                elif hostenv=='nuclear':
                    if nearby.iloc[0]['Seperation']<0.1:
                        promoted.append(sn)

                elif hostenv=='orphan':
                    if nearby.iloc[0]['Seperation']>=2:
                        promoted.append(sn)
                    elif nearby.iloc[0]['Seperation']<2.0:
                        key_objects.remove(sn)


    if len(key_objects)>0:
        for sn in key_objects:
            tempdata=candidates.get_group(sn).reset_index()
            successful=False
            try:
                nearby=PS1_query.search(sn,tempdata.ra.mean(),tempdata.dec.mean(),1).to_pandas().sort_values('Seperation').reset_index()
                successful=True
            except:
                continue
                
            if successful==False:
                hostflag.append(sn)
            if successful==True:
                if hostenv=='host':
                    if nearby.iloc[0]['Seperation']>0.1:
                        if nearby.iloc[0]['Seperation']<=2:
                           promoted.append(sn)
                    elif nearby.iloc[0]['Seperation']<=0.1:
                        nuclear.append(sn)
                    elif nearby.iloc[0]['Seperation']>2:
                        hostless.append(sn)
                elif hostenv=='nuclear':
                    if nearby.iloc[0]['Seperation']<0.1:
                        promoted.append(sn)

                elif hostenv=='orphan':
                    if nearby.iloc[0]['Seperation']>=2:
                        promoted.append(sn)
                    elif nearby.iloc[0]['Seperation']<2.0:
                        key_objects.remove(sn)


gold=list(set(promoted) - set(final_FLAG))

print('Total objects meeting selection criteria (%i/%i):'%(len(promoted),(len(candidates.groups.keys()))))


print('Interesting objects (%i/%i):'%(len(gold),(len(promoted))))
for i in gold:
	print(i)

if len(final_FLAG)>0:
    print ('Flagged objects (large gap between photometry) (%i/%i):'%(len(final_FLAG),(len(promoted))))
    for i in flagged:
        print(i)


if hostenv=='host':
    print('Of which (%i/%i) are nuclear transients:'%(len(nuclear),(len(promoted))))
    for i in nuclear:
        print(i)

if hostenv!=None:
    print('Of which (%i/%i) had no available host information:'%(len(hostflag),(len(promoted))))
    for i in hostflag:
        print(i)




