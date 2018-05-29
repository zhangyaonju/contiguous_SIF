##### data prepare
##### generate daily data
####   use sif_inst/PAR_inst*PAR_daily
####

import multiprocessing
#import h5py
import os
import numpy as np
import fnmatch
import sys
from netCDF4 import Dataset

def getallfile (directory):
    fileList=[]
    for path, subdirs, files in os.walk(directory):
        for name in files:
            if ".nc4" == name[-4:] or '.nc' == name[-3:]: fileList.append(os.path.join(path,name))
    #print fileList
    fileList.sort()
    return fileList

def calculate_clear_inst(date):
    sza_file = fnmatch.filter(szafiles,'*'+str(date)[4:7]+'.nc')[0]
    print('sza_files:')
    print(sza_file)
    with Dataset(sza_file, 'r') as fh:
        try:
            cossza_lat = fh.variables["cossza"][:].T
            print(cossza_lat.dtype)
        except:
            print("failed to import:"+sza_file)
        cossza = np.repeat(cossza_lat,7200).reshape(3600,7200)
    #sza = np.arccos(cossza)
    rtoa = 1360.8*0.98*(1+0.033*np.cos(2*np.pi*(date%1000+2)/365))*cossza
    tb = a0 + a1*np.exp((-k)/cossza)
    td = 0.271-0.294*tb
    rtot=rtoa*(tb+td)
    return rtot

def get_bess_par(date):
    bessdate = []
    for i in range(4):
        bess_day = fnmatch.filter(bessfiles,'*'+str(date+i)+'.nc')
        bessdate +=bess_day
    if len(bessdate)==0:
        return None
    print('bess file: ')
    print(bessdate)
    bessrad = []
    for f in bessdate:
        with Dataset(f,mode='r') as fh:
            try:
                parday = fh.variables['surface_downwelling_photosynthetic_radiative_flux_in_air'][:,:]*1.0
                print('bess radiation shape: '+str(parday.shape))
                parday[parday<=-999]=np.nan
                bessrad +=[parday.T]
            except:
                print('cannot read bess file: '+f)
    bessrad = np.array(bessrad)
    print('bessstack'+str(bessrad.shape))
    avgrad = np.flipud(np.nanmean(bessrad,axis=0))
    print('radiaiton shape: '+str(avgrad.shape))
    return avgrad

def get_clear_sif(doy):
    sif_file = fnmatch.filter(siffiles,'*'+str(doy)+'.nc')[0]
    print('sif file: ')
    print(sif_file)
    with Dataset(sif_file,mode='r') as fh:
        try:
            sif = fh.variables['clear_inst_SIF'][:,:]
            print('sif size: '+str(sif.shape))
            ##sif[sif< -990] = np.nan
            sif[sif< -990] = 0
        except:
            print('cannot read sif file: '+f)
    return sif

def exportsif(predY,outfile):
    print('creating: '+outfile)
    f = Dataset(outfile, 'w', format='NETCDF4')
    lati = f.createDimension('lat', len(lat))
    loni = f.createDimension('lon', len(lon))
    latitudes = f.createVariable('lat', 'f4', 'lat')
    longitudes = f.createVariable('lon', 'f4', 'lon')

    latitudes.units = 'degrees_north'
    longitudes.units = 'degrees_east'
    latitudes.standard_name = 'latitude'
    longitudes.standard_name = 'longitude'
    latitudes.axis = 'Y'
    longitudes.axis = 'X'
    latitudes.long_name = 'latitude'
    longitudes.long_name = 'longitude'
    latitudes[:] = lat
    longitudes[:] = lon
    oco2_sif = f.createVariable('all_daily_SIF', 'f4', ('lat','lon'),zlib=True)
    oco2_sif.missing_value=-999.9
    oco2_sif [:] = predY.astype('f4').reshape(3600,7200)
    f.close()


#### read dem and calculate the a0, a1, and k
demfile = '/rigel/glab/users/zy2309/DATA/dem005.nc'
with Dataset(demfile, 'r') as fh:
    try:
        dem = fh.variables["dem"][:,:].T
        print('dem shape:'+str(dem.shape))
    except:
        print("failed to import:"+demfile)
dem[dem>2500]=2500
dem/=1000
a0 = 0.4237-0.00821*np.square(6-dem)
a1 = 0.5055+0.00595*np.square(6.5-dem)
k = 0.2711+0.01858*np.square(2.5-dem)

def calc_daily_SIF(doy):
    #get inst sif
    clear_sif = get_clear_sif(doy)
    #get daily PAR
    all_par = get_bess_par(doy)
    if all_par is None:
        return
    #get inst PAR
    clear_par = calculate_clear_inst(doy)
    clear_par[clear_par<25] = 25
    #clear_file = '/rigel/glab/users/zy2309/DATA/all_daily_SIF_4day_CMG/clear_day_PAR.'+str(doy)+'.nc'
    #exportsif(clear_par,clear_file)
    daily_sif = clear_sif*all_par/(clear_par*0.45)

    ##out file
    outfile = '/rigel/glab/users/zy2309/DATA/all_daily_SIF_4day_CMG/OCO2.SIF.all.daily.'+str(doy)+'.nc'
    exportsif(daily_sif,outfile)



lat = np.arange(-90 + 0.05 / 2., (90 + 1e-8) - 0.05 / 2., 0.05)
lon = np.arange(-180 + 0.05 / 2., (180 + 1e-8) - 0.05 / 2., 0.05)

bess_dir = '/rigel/glab/users/Datasets/BESS/BESSRadiation/realtimedata.snu.ac.kr/BESSRadiation/BESS_PAR_Daily/'
cos_sza_dir = '/rigel/glab/users/zy2309/DATA/cos_sza/'
sif_inst_dir = '/rigel/glab/users/zy2309/DATA/clear_inst_SIF_4day_CMG/'

bessfiles = getallfile(bess_dir)
szafiles = getallfile(cos_sza_dir)
siffiles = getallfile(sif_inst_dir)
print(bessfiles)
print(szafiles)
print(siffiles)


doy = np.repeat(np.arange(2000,2017),92)*1000+np.tile(np.arange(1,366,4),17)
#doy = np.repeat(np.arange(2003,2018),92)*1000+np.tile(np.arange(1,366,4),15)
pool = multiprocessing.Pool(8)
pool.map(calc_daily_SIF,doy)
#calc_daily_SIF(2000153)


