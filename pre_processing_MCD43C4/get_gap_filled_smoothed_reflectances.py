import multiprocessing
import os, sys, time
import numpy as np
import gdal
import time
import datetime
from netCDF4 import Dataset
import fnmatch
import pandas as pd
import math

#global root, outdir
root = '/rigel/glab/users/zy2309/DATA/'
outdir = '/rigel/glab/users/zy2309/DATA/MCD43C4_4day_GF/'
interval =4
obs = np.int(math.ceil(366.0/interval))

def getfile (directory):
    fileList=[]
    for path, subdirs, files in os.walk(directory):
        for name in files:
            if ".nc4" == name[-4:] or '.hdf' == name [-4:] or '.nc'==name[-3:]:
                fileList.append(os.path.join(path,name))
    #print fileList
    fileList.sort()
    return fileList


def getyeardoy (file):
    filetime = datetime.datetime.strptime(os.path.basename(file)[19:27],"%Y%m%d")
    yday = filetime.year*1000+(filetime - datetime.datetime(filetime.year,1,1)).days+1
    return yday


#use to smooth the data: 
def dataframeapply (df):
    df = pd.DataFrame(np.concatenate([df[(obs/2):obs, :], df, df[:obs/2, :]]))
    print df.shape
    df_smoothed = df.interpolate(axis=0)
    df_smooth = np.array(df_smoothed)
    df_select = df_smooth[(obs/2):(obs/2*3),:]
    return df_select


ncfiles = getfile(root+'/MCD43C4_4day/')
reffiles = getfile(root+'/MCD43C4_4day_ref/')

######
## main function for retrieval
band_names = ['b1','b2','b3','b4']#,'b5','b7']


dlat = np.arange(-90, 90 + 1e-8, 0.05)
dlon = np.arange(-180, 180 + 1e-8, 0.05)
lat = np.arange(-90 + 0.05 / 2., (90 + 1e-8) - 0.05 / 2., 0.05)
lon = np.arange(-180 + 0.05 / 2., (180 + 1e-8) - 0.05 / 2., 0.05)


def read_annual(year,bname):
    pat = '*.4day.'+str(year)+'*'
    yearfile = fnmatch.filter(ncfiles,pat)
    yearfile.sort()
    reffiles.sort()
    refdata = np.zeros((92,3600*7200))
    banddata = np.zeros((92,3600*7200))
    for i in range(0,len(yearfile)):
        with Dataset(yearfile[i],'r') as fin:
            #print (os.path.basename(yearfile[i])[20:23])
            cmd = 'banddata['+str((int(os.path.basename(yearfile[i])[20:23])-1)/4)+',:]=(fin.variables["'+bname+'"])[:,:].flatten()'
            print str(year)+'    '+cmd
            exec cmd
        with Dataset(reffiles[i],'r') as fin:
            cmd = 'refdata['+str((int(os.path.basename(yearfile[i])[20:23])-1)/4)+',:]=(fin.variables["'+bname+'"])[:,:].flatten()'   
            exec cmd
    banddata[banddata<=0] = np.nan
    refdata[refdata < -1000] = np.nan
    banddata -=refdata*1.0
    banddata = dataframeapply(banddata)+refdata
    banddata = np.where(np.isnan(banddata),-9999,banddata)
    return banddata

def process_year_reflectance(year):
    for doy in range(0,366,interval):
        outfile = outdir+'MCD43C4CMG.'+str(interval)+'day.'+str(year)+str(doy+1).zfill(3)+'.rec.nc'
        print 'creating: '+outfile
        f = Dataset(outfile, 'w', format='NETCDF4')

        lati = f.createDimension('lat', len(lat))
        loni = f.createDimension('lon', len(lon))
        latitudes = f.createVariable('lat', 'f8', 'lat')
        longitudes = f.createVariable('lon', 'f8', 'lon')

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
        f.close()
    for bname in band_names:
        vec_data=read_annual(year,bname)
        for doy in range (0,366,interval):
            outfile = outdir+'MCD43C4CMG.'+str(interval)+'day.'+str(year)+str(doy+1).zfill(3)+'.rec.nc'
            print 'exporting: '+outfile
            f = Dataset(outfile, 'a', format='NETCDF4')
            cmd = bname+'=f.createVariable("'+bname+'", "i2", ("lat","lon"),zlib=True)'
            exec cmd
            cmd = bname+'.missing_value=-9999'
            exec cmd
            cmd = bname+'[:]=np.round(vec_data['+str(doy/4)+',:].reshape(3600,7200))'
            exec cmd
            f.close()
        vec_data = None

pool = multiprocessing.Pool(2)

years = np.arange(2000,2006)
#years = np.arange(2006,2012)
#years = np.arange(2012,2018)

#process_year_reflectance(2012)
pool.map(process_year_reflectance, years)
