import multiprocessing
import os, sys, time
import numpy as np
import gdal
import time
import datetime
from netCDF4 import Dataset
import fnmatch
import math


def getfile(directory):
    fileList=[]
    for path, subdirs, files in os.walk(directory):
        for name in files:
            if ".nc4" == name[-4:] or '.hdf' == name [-4:]: 
                fileList.append(os.path.join(path,name))
    #print fileList
    fileList.sort()
    return fileList

def gethdffile (yday):
    pat = "*A"+str(yday)+"*"
    #print(pat)
    try:
        hdffile = fnmatch.filter(hdffiles,pat)[0]
        return hdffile
    except:
        return None

dict_name = {'b1': 'Nadir_Reflectance_Band1',
            'b2': 'Nadir_Reflectance_Band2',
            'b3': 'Nadir_Reflectance_Band3',
            'b4': 'Nadir_Reflectance_Band4'}
band_name = ['b1','b2','b3','b4']

def combinenc (startdate):
    dlat = np.arange(-90, 90 + 1e-8, 0.05)
    dlon = np.arange(-180, 180 + 1e-8, 0.05)
    lat = np.arange(-90 + 0.05 / 2., (90 + 1e-8) - 0.05 / 2., 0.05)
    lon = np.arange(-180 + 0.05 / 2., (180 + 1e-8) - 0.05 / 2., 0.05)
    
    dayfile = []
    for i in range(0,interval):
        dayfile.append(gethdffile(startdate+i))
    dayfile = [x for x in dayfile if x is not None]
    print dayfile
    if len(dayfile)==0:
        return
    outfile = outdir+'MCD43C4CMG.'+str(interval)+'day.'+str(startdate)+'.nc'
    #if file exist, remove first
    try:
        os.remove(outfile)
    except OSError:
        pass
    #create the new files
    f = Dataset(outfile, 'w', format='NETCDF4')
    for i in dict_name:
        cmd = 'f.' + i + '="' + dict_name[i] + '"'
        try:
            exec cmd
        except:
            print 'Error executing ' + cmd

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
    
    for bname in band_name:
        cmd = bname+'=f.createVariable("'+bname+'", "i2", ("lat","lon"),zlib=True)'
        exec cmd
        cmd = bname+'.missing_value= -9999'
        exec cmd
    
    for bname in band_name:
        cmd = 'all_'+bname +'=np.zeros((3600,7200,4))'
        exec cmd
        cmd = 'all_'+bname +'[:]=np.nan'
        exec cmd
    
    for ifile in range(len(dayfile)):
        print ifile
        bands = {}
        ds_mod09 = gdal.Open(dayfile[ifile])
        subdata = ds_mod09.GetSubDatasets()
        subdataname = [item[0] for item in subdata]
        for i in dict_name:
            dat_name = fnmatch.filter(subdataname,"*"+dict_name[i]+"*")[0]
            dataset = gdal.Open(dat_name)
            bands[i] = dataset.GetRasterBand(1).ReadAsArray()
            bands[i] = np.where(bands[i]>= 10000,np.nan,bands[i])
            print bands[i].shape
            cmd = 'print all_'+i+'.shape'
            exec cmd 
            cmd = 'all_'+i +'[:,:,ifile]= bands["'+i+'"][:,:]'
            print cmd
            exec cmd
    
    for bname in band_name:
       	cmd = bname +'[:,:] = np.flipud(np.round(np.nanmean(all_'+bname+',axis=2)))'
        exec cmd
    f.close()

def process_list(start_date, mp = True, count = 1):    
    if mp:
        pool = multiprocessing.Pool(processes=count)
        #print files
        pool.map(combinenc, start_date)

global root, outdir, hdffiles, interval
interval = 4
root = '/rigel/glab/users/zy2309/DATA/'
outdir = root+'/MCD43C4_'+str(interval)+'day/'
if not os.path.exists(outdir):
    os.makedirs(outdir)
hdffiles = getfile(root+'/MCD43C4CMG/')

#start_date = np.repeat(np.arange(2000,2009),math.ceil(366.0/interval))*1000+\
#            np.tile(np.arange(1,366,interval),9)
start_date = np.repeat(np.arange(2009,2018),math.ceil(366.0/interval))*1000+\
            np.tile(np.arange(1,366,interval),9)
process_list(start_date,mp=True,count=24)    

#combinenc(2016009)
