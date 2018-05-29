#import tensorflow as tf
import multiprocessing
#import h5py
import numpy as np
#import matplotlib.pyplot as plt
import random
import os
import time
from netCDF4 import Dataset

################################################################
# The SIF is reproduced using the 4 band refectacne and the solar zenith angle
root = '/rigel/glab/users/zy2309/'
datadir = root+'DATA/MCD43C4_4day_GF/'
szadir = root+'DATA/cos_sza/'
conver = root+'DATA/conversion_factor/'
outdir = root+'DATA/clear_inst_SIF_4day_CMG/'
modeldir = root+'PROJECT/SIF_OCO2_ML/trained_model/1layer_5band/'
# two files used for normalizaiton
mean_file = root+'PROJECT/SIF_OCO2_ML/habanero/training/mcd43_b7_sza_vpd_temp_sif_training_mean.nc'
std_file = root+'PROJECT/SIF_OCO2_ML/habanero/training/mcd43_b7_sza_vpd_temp_sif_training_std.nc'

inputX = ['b1','b2','b3','b4']
layer = ['1','out']
inputY = ['sif_757nm']

lat = np.arange(-90 + 0.05 / 2., (90 + 1e-8) - 0.05 / 2., 0.05)
lon = np.arange(-180 + 0.05 / 2., (180 + 1e-8) - 0.05 / 2., 0.05)

#####generate the prediction by importing the weights and biases from inspect_checkpoint files
def readfile(lay,ind):
    file = modeldir+'layer'+lay+'_'+ind+'_Variable'  # these files were generated using the inspect_checkpoints 
    para = np.loadtxt(file,delimiter=' ',dtype='f8',ndmin=2)
    return para

###### export the RSIF data as an nc file
def exportsif(predY,dailyY,outfile):
    print('creating: '+outfile)
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
    predY = np.where(np.isnan(predY),-9999,predY)
    inst_RSIF = f.createVariable('clear_inst_SIF', 'f4', ('lat','lon'),zlib=True)
    daily_RSIF = f.createVariable('clear_daily_SIF','f4',('lat','lon'),zlib=True)
    inst_RSIF.missing_value=-9999
    daily_RSIF.missing_value=-9999
    inst_RSIF[:] = predY.astype('f4').reshape(3600,7200)
    daily_RSIF[:] = dailyY.astype('f4').reshape(3600,7200)
    f.close()


######### this function read the doy reflectance, the calcualted SZA and the model for prediction. 
def predict_doy(doy):
    inputs = []
    data = []
    datstd = {}
    datmean = {}
    #get the output file
    outfile = outdir+'OCO2.SIF.clear.inst.'+str(doy)+'.nc'
    #get the sza file
    szafile = szadir+'cos_sza.4day.'+str(doy)[4:7]+'.nc'
    conversionfile = conver+'con_fac.4day.'+str(doy)[4:7]+'.nc'
    #get the reflectance file
    datfile = datadir+'MCD43C4CMG.4day.'+str(doy).zfill(3)+'.rec.nc'
    if not os.path.isfile(datfile):
        return
    with Dataset(std_file, mode='r') as fh:
        for k in inputX:
            try:
                datstd[k] = fh.variables[k][None,:]
            except:
                datstd[k] = np.array(fh.variables[k])[None]
    with Dataset(mean_file, mode='r') as fh:
        for k in inputX:
            try:
                datmean[k] = fh.variables[k][None,:]
            except:
                datmean[k] = np.array(fh.variables[k])[None]
    with Dataset(datfile, mode='r') as fh:
        for k in inputX:
            try:
                arr = fh.variables[k][:,:]*1.0
            except:
                arr = np.array(fh.variables[k][:])[None,:].T
            arr[arr<0.50] = np.nan 
            arr -=datmean[k]
            arr /=datstd[k]
            #print(arr.shape)
            inputs += [arr.flatten()]
    print(szafile)
    with Dataset(szafile, 'r') as fh:
        try:
            cossza_lat = fh.variables["cossza"][:].T
            print(cossza_lat.dtype)
        except:
            print("failed to import:"+szafile)
        cossza = np.repeat(cossza_lat,7200)
        print(len(cossza))
        inputs +=[cossza]
    with Dataset(conversionfile,'r') as fh:
        try:
            confac_lat = fh.variables["confac"][:].T
            print(confac_lat.dtype)
        except:
            print("failed to import:"+conversionfile)
        confac = np.array(np.repeat(confac_lat,7200))
        #print('conver_factor'+str(len(confac)))
        #print(confac)
    #print(input)
    num = 3600*7200
    net = np.array(inputs).T
    #with tf.Session() as sess:
    # build the models
    for i in layer:
        wei = readfile(i,'weights')
        #print(wei.dtype)
        bia = readfile(i,'biases').T*np.ones((num,1))
        #print(bia)
        #print(i)
        #print(net.shape)
        #print(wei.shape)
        #print(bia.shape)
        prenet = np.dot(net,wei)+bia
        if i == 'out':
            predY = prenet
        else:
            prenet[prenet<0]=0
            net = prenet
    #print(predY.T)
    dailyY = predY.T[:]*confac
    exportsif(predY,dailyY,outfile)


#main 
DOY = np.repeat(np.arange(2000,2009),92)*1000+np.tile(np.arange(1,366,4),9)
#DOY = np.repeat(np.arange(2009,2018),92)*1000+np.tile(np.arange(1,366,4),9)
pool = multiprocessing.Pool(24)
pool.map(predict_doy,DOY)
#predict_doy(DOY[1])



    
