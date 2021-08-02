#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""scrSimul
Function to perform SCR simulation over specified location on the 
Sentinel-1 coregistered stack

input = (stackDir, latitude [deg], longitude [deg], elevation [m], 
         cropSize [m], RCS_CR [dBm2], ovsFactor = 16)
output = Single-band GeoTIFF with SCR [dB] in WGS-84

"""

print('Geodetic Corner Reflector (In)SAR Toolbox (GECORIS) v.1.0')
print('Copyright (c) 2021 Richard Czikhardt, czikhardt.richard@gmail.com')
print('Dept. of Theoretical Geodesy, Slovak University of Technology')
print('-----------------------------------------------------------------')
print('License: GNU GPL v3+')
print('-----------------------------------------------------------------')

import sys
import numpy as np
from scipy import interpolate
try:
    import gdal
except:
    from osgeo import gdal
from gecoris import classes, radarUtils, s1Utils, geoUtils, ioUtils

def scrSimul(stackDir, outTiff, latitude, longitude, elevation, cropSize = 200, 
             RCS_CR = 30, ovsFactor = 16):
    
    # central point
    coords = classes.Coordinates(float(latitude)/180*np.pi,
                                 float(longitude)/180*np.pi,
                                 float(elevation))        
    # read stack data:
    # TODO: for now hardcoded S-1
    stack = classes.Stack('scr_stack',
                           stackDir,'Sentinel-1',None,stackType = 'coreg')
    stack.readData()
    
    # radarcode centrail point and crop:
    (Azimuth,Range) = radarUtils.radarcode(coords.plh,stack.masterMetadata)
    boundingBox = radarUtils.getCrop(Azimuth[0],Range[0],
                                     stack.masterMetadata,cropSize)
    # allocate time series array
    azSize = boundingBox[0][1]-boundingBox[0][0]+1
    rgSize = boundingBox[1][1]-boundingBox[1][0]+1
    beta0array = np.zeros((azSize*ovsFactor,rgSize*ovsFactor,
                           stack.nSLC))
    # loop through SLCs:
    idx = 0
    for file,acqDate in zip(stack.files,stack.acqDates):
        metadata = stack.metadata[acqDate]
        SLCderamp = s1Utils.readSLC(file,boundingBox,
                                    method = 'coreg',deramp = True)
        SLCovs = radarUtils.oversample(SLCderamp,ovsFactor)
        beta0array[:,:,idx] = np.power(np.abs(SLCovs),2)/(metadata['beta0']**2)
        idx += 1
    
    # temporal average beta0 in dB scale (closer to normal distrib.):
    meanBeta_dB = np.nanmean(10*np.log10(beta0array),axis=2)
    del beta0array
    # compute oversampled SCR:
    A_res = 10*np.log10(metadata['azimuthResolution']*metadata['rangeResolution'])
    SCR = RCS_CR - (meanBeta_dB + A_res);
    
    # load DEM heights of AOI
    hGrid = s1Utils.readHGT(file,boundingBox)
    # original radar grid:
    [azGrid_orig,rGrid_orig] = np.mgrid[boundingBox[0][0]:boundingBox[0][1]+1,
                                        boundingBox[1][0]:boundingBox[1][1]+1]
    # geocode orig. grid:
    xyz = radarUtils.radar2xyz(azGrid_orig.flatten(),rGrid_orig.flatten(),
                               hGrid.flatten(),stack.masterMetadata)
    plh = geoUtils.xyz2plh(xyz)
    # grid plh:
    latGrid_orig = plh[:,0].reshape(hGrid.shape)
    lonGrid_orig = plh[:,1].reshape(hGrid.shape)
    # up-sampled radar grid:
    azVec = np.linspace(boundingBox[0][0], boundingBox[0][1], num=SCR.shape[0])
    rVec = np.linspace(boundingBox[1][0], boundingBox[1][1], num=SCR.shape[1])    
    rGrid,azGrid = np.meshgrid(rVec,azVec)
    # up-sample plh using bi-linear interpolation:
    latGrid_orig = interpolate.griddata((azGrid_orig.flatten(), rGrid_orig.flatten()),
                                   latGrid_orig.flatten(),(azGrid, rGrid), 
                                   method='linear')  
    lonGrid_orig = interpolate.griddata((azGrid_orig.flatten(), rGrid_orig.flatten()),
                                   lonGrid_orig.flatten(),(azGrid, rGrid), 
                                   method='linear')    
    # create regular geo-grid for output raster:
    spacing = max(metadata['azimuthSpacing'],metadata['rangeSpacing'])/ovsFactor
    latSpacing = spacing/6371e3
    lonSpacing = spacing/6371e3*np.cos(np.mean(latGrid_orig))
    [latGrid,lonGrid] = np.mgrid[np.max(latGrid_orig):np.min(latGrid_orig):
                                 -latSpacing,
                                 np.min(lonGrid_orig):np.max(lonGrid_orig):
                                 lonSpacing]
    # interpolate SCR over regular grid:
    scrGrid = interpolate.griddata((latGrid_orig.flatten(), 
                                    lonGrid_orig.flatten()),
                                    SCR.flatten(),(latGrid, lonGrid), 
                                    method='nearest')   
    scrGrid = scrGrid.astype('float32') 
    # write geotiff using gdal:
    try:
        drv = gdal.GetDriverByName("GTiff")
        ds = drv.Create(outTiff, scrGrid.shape[1], scrGrid.shape[0], 
                        1, gdal.GDT_Float32)
        geotransform = (np.min(lonGrid)*180/np.pi, lonSpacing*180/np.pi, 0, 
                        np.max(latGrid)*180/np.pi, 0, -latSpacing*180/np.pi)
        sr = gdal.osr.SpatialReference()
        sr.ImportFromEPSG(4326)
        ds.SetGeoTransform(geotransform)
        ds.SetProjection(sr.ExportToWkt())
        ds.GetRasterBand(1).WriteArray(scrGrid)
        ds = None
    except:
        print('Cannot write GeoTIFF.')
    else:
        print('Simulated SCR GeoTIFF written.')
        
if __name__ == "__main__":
    # load parameters:
    # if len(sys.argv) == 6:
    #     scrSimul(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5])
    # elif len(sys.argv) == 7:
    #     scrSimul(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],
    #              sys.argv[6])
    # elif len(sys.argv) == 8:
    #     scrSimul(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],
    #              sys.argv[6],sys.arv[7])
    # elif len(sys.argv) == 9:
    #     scrSimul(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],
    #              sys.argv[6],sys.argv[7],sys.argv[8])
    # else:
    #     print('Not enough input arguments!')
    parms = ioUtils.loadJSON(sys.argv[1])
    scrSimul(parms['stackDir'],parms['outTiff'],parms['latitude'],
             parms['longitude'],parms['elevation'],parms['cropSize'],
             parms['RCS'],parms['oversamplingFactor'])
