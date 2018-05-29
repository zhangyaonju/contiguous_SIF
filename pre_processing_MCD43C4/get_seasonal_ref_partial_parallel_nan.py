import multiprocessing
import os
from ctypes import *
import numpy as np
from netCDF4 import Dataset
import time
import pandas as pd
import fnmatch
#import h5py
import math
import datetime

def identifybad (timeseries):
    mask = np.concatenate(([False],np.isnan(timeseries),[False]))
    bad_ind = np.repeat(False,len(timeseries))
    if ~mask.any():
        return bad_ind
    else:
        idx = np.nonzero(mask[1:] != mask[:-1])[0]
        consective_na = idx[1::2] - idx[::2]
        start_na_id = np.nonzero(consective_na>=15)[0]
        for i in start_na_id*2:
            bad_ind[idx[i]:idx[i+1]]=True
        return bad_ind


def getallfile (directory):
    fileList=[]
    for path, subdirs, files in os.walk(directory):
        for name in files:
            if ".nc" == name[-3:] : fileList.append(os.path.join(path,name))
    #print fileList
    fileList.sort()
    return fileList

def getbandmedian (doy):
    pat = "*"+str(doy+1).zfill(3)+'.nc'
    doyfile = fnmatch.filter(ncfiles,pat)
    num_years = len(doyfile)
    print 'number of years to read in:'+str(num_years)+'   for DOY:'+str(doy+1)
    #initialize array for doy multiyear
    cmd = 'doy_data=np.zeros(('+str(num_years)+',3600*7200))'
    exec cmd
    print 'band created for:'+bname+'   DOY:'+str(doy) 
    #read in doy-year data
    for i in range(0,num_years):
        fin = Dataset(doyfile[i], 'r')
        cmd = 'doy_data['+str(i)+',:]=(fin.variables["'+bname+'"])[:,:].flatten()'
        print str(doy+1)+'   '+cmd
        exec cmd
    #get the median for each doy
    cmd = 'doy_data[doy_data<=0]=np.nan'
    exec cmd
    cmd = 'median_bands=np.nanmedian(doy_data,axis=0).flatten()'
    print str(doy+1)+'   '+cmd
    exec cmd
    return median_bands


def parallelize_dataframe(df):
    df_split = np.array_split(df, 8, axis=1)
    pool = multiprocessing.Pool(8)
    df = np.concatenate(pool.map(dataframeapply, df_split), axis=1)
    pool.close()
    pool.join()
    return df

#use to smooth the data: 
def dataframeapply (df):
    # test the reflectance with a moving filter   
    # add 4/27/2018 by yaozhang
    baddata = np.apply_along_axis(identifybad,0,df)
    df = pd.DataFrame(np.concatenate([df[int(obs/2):obs, :], df, df[:int(obs/2), :]]))
    print df.shape

    df_smoothed = df.interpolate(axis=0)
    df_smooth = np.array(df_smoothed)
    df_select = df_smooth[(obs/2):(obs/2*3),:]
    return np.where(baddata,-9999,df_select) 
        
global ncfiles, interval, indir, obs, band_name, outdir, bname
interval = 4
obs = np.int(math.ceil(366.0/interval))
#band_name = ['b1','b2','b3','b4']#,'b5','b7']
band_name = ['b5','b6','b7']

print obs
print 'start time:'+str(datetime.datetime.now())

root = '/rigel/glab/users/zy2309/'
indir = root+'/DATA/MCD43C4_'+str(interval)+'day/'
outdir = root+'/DATA/MCD43C4_'+str(interval)+'day_ref/'
if not os.path.exists(outdir):
    os.makedirs(outdir)

ncfiles = getallfile(indir)
#print ncfiles
#getref(ncfiles)
dlat = np.arange(-90, 90 + 1e-8, 0.05)
dlon = np.arange(-180, 180 + 1e-8, 0.05)
lat = np.arange(-90 + 0.05 / 2., (90 + 1e-8) - 0.05 / 2., 0.05)
lon = np.arange(-180 + 0.05 / 2., (180 + 1e-8) - 0.05 / 2., 0.05)

'''
for i in range(0,366,interval):
    outfile = outdir+'MCD43C4.'+str(interval)+'day.'+str(i+1).zfill(3)+'.ref.nc'   
    print 'creating:'+outfile
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
'''
for bname in band_name:
    i = np.arange(0,366,interval)

    pool = multiprocessing.Pool(8)
    medianband= pool.map(getbandmedian, i)
    pool.close()
    pool.join()
    print 'median_band shape:'+str(np.array(medianband).shape)
    med = np.array(medianband)
 
    print bname+" imported!!!!!!!!!"

    print "start interpolation for band:"+bname
    t_start=time.time()
    #vec_dat=parallelize_dataframe(med)
    vec_dat=dataframeapply(med)
    t_end = time.time()
    med = None
    medianband = None
    print 'time used for '+bname+':'+str(datetime.timedelta(seconds = t_end-t_start))
    
    for i in range(0,366,interval):
        outfile = outdir+'MCD43C4.'+str(interval)+'day.'+str(i+1).zfill(3)+'.ref.nc'   
        print 'exporting:'+outfile
        f = Dataset(outfile, 'a', format='NETCDF4')
        cmd = bname+'=f.createVariable("'+bname+'", "i2", ("lat","lon"),zlib=True)'
        exec cmd
        cmd = bname+'.missing_value=-9999'
        exec cmd
        cmd = bname+'[:]=np.round(vec_dat['+str(i/4)+',:].reshape(3600,7200))'
        exec cmd
        f.close()
    vec_data=None


