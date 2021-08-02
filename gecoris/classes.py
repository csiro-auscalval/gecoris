# -*- coding: utf-8 -*-
"""GECORIS classes

Module containing definition of classes

Copyright (C) 2021 by R.Czikhardt

Email: czikhardt.richard@gmail.com
Last edit: 7.3.2021

This file is part of GECORIS - Geodetic Corner Reflector (In)SAR Toolbox.

    GECORIS is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    GECORIS is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with GECORIS. If not, see <https://www.gnu.org/licenses/>.
"""
###############################################################################

import os
import sys
import glob
from collections import OrderedDict
import datetime
import json
import numpy as np
import scipy.stats as ss
from scipy.spatial.distance import cdist
if sys.stdout.isatty(): # test if interactive shell or not
    import matplotlib
    matplotlib.use('Agg') # change plotting backend accordingly
import matplotlib.pyplot as plt
from matplotlib import dates
formatter = dates.DateFormatter('%m-%y')
# gecoris packages:
from gecoris import geoUtils, crUtils, radarUtils, s1Utils
from gecoris import plotUtils, atmoUtils, insarUtils
# constants:
speedOfLight = 299792458.0
# types:
stackType = np.dtype({'names':('acqDate','acqTime','status','RCS','SLC','azimuth','range',
                               'dAz','dR','dist','IFG'),
                      'formats':('U10','U10','?','f8','c16','f8','f8',
                                 'f8','f8','f8','c16')})
# search window parameters (do not change):
coregCropSizeApprox = 2
coregCropSizePrecise = 0.5
rawCropSizeApprox = 3
rawCropSizePrecise = 2


class Reflector:

    def __init__(self,reflId,reflType,ascending,descending,
                 installDate,startDate,endDate):
        self.id = reflId
        self.type = reflType
        self.ascending = ascending
        self.descending = descending
        self.installDate = installDate
        self.startDate = startDate
        self.endDate = endDate
        self.stacks = []        
        
    def get_xyz(self,orbit = 'ASCENDING'):
        if orbit == 'ASCENDING':
            return self.ascending.xyz
        elif orbit == 'DESCENDING':
            return self.descending.xyz

    def get_plh(self,orbit = 'ASCENDING'):
        if orbit == 'ASCENDING':
            return self.ascending.plh,self.ascending.crs
        elif orbit == 'DESCENDING':
            return self.descending.plh,self.ascending.crs
        else:
            raise Exception('Wrong orbit, specify ASCENDING or DESCENDING')

        
    def addStack(self, stack, ovsFactor=1, posFlag = 1, plotFlag=1, **kwargs):
        """Extract measurement timeseries from specified stack
    
        Performed as step-wise procedure:
            (i) radarcoding reflector's geo coordinates (using master frame
                if coreg stack or individual acqisition frames for raw stack)
            (ii) reading SLC crops around radar coordinates from bursts 
                prepared by SNAP
            (iii) oversampling crop by FFT zero-padding using specified 
                'ovsFactor'
            (iv) fitting elliptic paraboloid to oversampled sub-crop of 
                1x1 resolution cells
            (v) witing extracted peak amplitude, radar coordinates, SLC
                into timeseries structure of 'stackType'
        """
        # test if stack geometry viable for reflector and get proper 
        # ellipsoidal coordinates (plh):
        if stack.orbit == 'ASCENDING':
            if self.ascending is not None:
                plh = self.ascending.plh
                xyz = self.ascending.xyz
                crs = self.ascending.crs
            else:
                print('Reflector ' +self.id+ 
                      'not valid for ASCENDING orbit, skipping.')
                return
        elif stack.orbit == 'DESCENDING':
            if self.descending is not None:
                plh = self.descending.plh
                xyz = self.descending.xyz
                crs = self.descending.crs
            else:
                print('Reflector ' +self.id+ 
                      'not valid for DESCENDING orbit, skipping.')
                return
        else:
            print('Not a valid stack geometry, skipping.')
            return
        
        print('Extracting timeseries from stack '+stack.id+' for reflector '
              +self.id)
        # perform data extraction differently for 'raw' or 'coreg stacks:
        if stack.type == 'coreg':
            if posFlag:
                cropSize = coregCropSizePrecise # TODO: remove
                #cropSize = 0.3
            else:
                cropSize = coregCropSizeApprox
            # allocate:
            data = np.zeros(len(stack.acqDates), dtype=stackType)
            # radarcode & crop:
            masterMetadata = stack.masterMetadata
            if 'atmoFile' in kwargs:
                (Azimuth,Range) = radarUtils.radarcode(plh, masterMetadata,
                                                       atmoDir=kwargs['atmoDir'],
                                                       atmoFile=kwargs['atmoFile'],
                                                       crs=crs)
            elif 'atmoDir' in kwargs:
                (Azimuth,Range) = radarUtils.radarcode(plh, masterMetadata,
                                                       atmoDir=kwargs['atmoDir'],
                                                       crs=crs)
            else:
                (Azimuth,Range) = radarUtils.radarcode(plh, masterMetadata,
                                                       crs=crs)   
            # fix inner delay for CAT
            if self.type == 'CAT': # TODO: remove
                Range += 1.5/masterMetadata['rangeSpacing']
            boundingBox = radarUtils.getBoundingBox(Azimuth[0],Range[0],
                                               masterMetadata,cropSize)
            # loop through extracted bursts:
            slcIdx = 0            
            for file,acqDate in zip(stack.files,stack.acqDates):
                #print('Extracting date '+acqDate)
                # get status:
                status = (acqDate > self.startDate and acqDate < self.endDate)
                if ovsFactor == 1:
                    SLC = s1Utils.readSLC(file,boundingBox,method='coreg',
                                        deramp = False)
                    IFG = s1Utils.readIFG(file,boundingBox)
                    metadata = stack.metadata[acqDate]
                    beta0 = np.power(np.abs(SLC),2)/(metadata['beta0']**2)
                    if plotFlag > 1:
                        if not 'outDir' in kwargs:
                            raise Exception('You must specify output directory for plots')
                        plotUtils.plotBeta0(SLC,metadata,
                                            kwargs['outDir'],ovs=False) 
                    
                    # find simple max:
                    idx = np.unravel_index(np.argmax(beta0, axis=None), beta0.shape)
                    # peak measurement:
                    peakSLC = SLC[idx]
                    peakIFG = IFG[idx]
                    RCSdB = radarUtils.beta2RCSdB(beta0[idx],metadata)
                    # coord. differences [m]
                    peakAz = boundingBox[0][0] + idx[0]
                    peakR = boundingBox[1][0] + idx[1]
                    dAz = (Azimuth[0]-peakAz)*metadata['azimuthSpacing']
                    dR = (Range[0]-peakR)*metadata['rangeSpacing']
                    # peak range:
                    # TODO: THis dist is wrong for coreg!
                    dist = (peakR/metadata['RSR'] + metadata['range0time'])*speedOfLight
                    # zero-Doppler acq.time:
                    tAzimuth,_ = radarUtils.radar2time(peakAz,peakR,metadata)
                    acqTime = geoUtils.sod2UTC(tAzimuth)
                    data[slcIdx] = (acqDate,acqTime,status,RCSdB,peakSLC,peakAz,peakR,dAz,dR,dist,peakIFG)
                slcIdx += 1
            
            # get zenith-azimuth-slantRange:
            _,_,satVec = radarUtils.xyz2t(xyz,masterMetadata)
            zas = geoUtils.xyz2zas(xyz,satVec)
            
            # create time series structre
            ts = {'id': stack.id, 'data': data, 'metadata': masterMetadata, 
                  'stack': stack,'type':'coreg','zas':zas}
            self.stacks.append(ts)
            
        elif stack.type == 'raw':
            if posFlag:
                cropSize = rawCropSizePrecise
            else:
                cropSize = rawCropSizeApprox
            #allocate:
            data = np.zeros(len(stack.acqDates), dtype=stackType)
            #% loop through extracted bursts:
            slcIdx = 0
            for file,acqDate in zip(stack.files,stack.acqDates):
                #print('Extracting date '+acqDate)
                # get status:
                status = (acqDate > self.startDate and acqDate < self.endDate)
                # radarcode & crop:
                metadata = stack.metadata[acqDate]
                if 'atmoFile' in kwargs:
                    (Azimuth,Range) = radarUtils.radarcode(plh, metadata,
                                                           atmoDir=kwargs['atmoDir'],
                                                           crs=crs,
                                                           atmoFile=kwargs['atmoFile'])
                    #print(Range)
                elif 'atmoDir' in kwargs:
                    (Azimuth,Range) = radarUtils.radarcode(plh, metadata,
                                                           atmoDir=kwargs['atmoDir'],
                                                           crs=crs)
                else:
                    (Azimuth,Range) = radarUtils.radarcode(plh, metadata,
                                                           crs=crs)
                # fix transponder's inner delay:
                #if self.type == 'CAT':
                #    Range += 1.5/metadata['rangeSpacing']
                if ovsFactor > 1:
                    # extract larger crop to avoid border effects during oversampling:
                    boundingBox = radarUtils.getBoundingBox(Azimuth[0],
                                                   Range[0],metadata,cropSize)
                    # read SLC:
                    SLCderamp = s1Utils.readSLC(file,boundingBox,method='raw',
                                              deramp = True)   
                    SLCovs = radarUtils.oversample(SLCderamp,ovsFactor)
                    # calibrate:
                    beta0 = np.power(np.abs(SLCovs),2)/(metadata['beta0']**2)
                    # TODO: fix NESZ:
                    
                    
                    if plotFlag > 1:
                        if not 'outDir' in kwargs:
                            raise Exception('You must specify output directory for plots')
                        plotUtils.plotBeta0(SLCderamp,metadata,kwargs['outDir'])
                    
                    if status:                    
                        # estimate peak coordinates and beta0:
                        peakAz,peakR,peakBeta0,localAz,localR = radarUtils.estimatePeak(
                            beta0,boundingBox,metadata,ovsFactor,method='fit')                
                        # positional ovs. differences expected vs. observed [m]:
                        dAz = (Azimuth[0]-peakAz)*metadata['azimuthSpacing']
                        dR = (Range[0]-peakR)*metadata['rangeSpacing']                
                        # peak SLC measurement from non-oversampled and not-deramped SLC:
                        #peakSLC = SLCraw[int(round(peakAz)-boundingBox[0][0]),
                        #                 int(round(peakR)-boundingBox[1][0])]
                        #TODO: or perform re-ramping of SLCovs and extract it - to be tested
                        #idx = np.unravel_index(np.argmax(beta0, axis=None), beta0.shape)  
                        dPhase_ovs = s1Utils.getDerampDemodPhase(file,boundingBox,
                                                               ovsFactor)
                        #SLCovs_reramp = utils.reramp(SLCovs,dPhase_ovs)
                        #peakSLC = SLCovs_reramp[idx]
                        # interpolate dPhase and re-ramp:
                        peakSLC = SLCovs[int(localAz),int(localR)]
                        dPhase_interp = radarUtils.interp2d(dPhase_ovs,
                                                            localAz,localR)
                        peakSLC = s1Utils.reramp(peakSLC,dPhase_interp)
                    else:
                        idx = (int(np.round(Azimuth[0])-boundingBox[0][0]),
                               int(np.round(Range[0])-boundingBox[1][0]))
                        peakAz = idx[0]
                        peakR = idx[1]
                        peakBeta0 = beta0[peakAz*ovsFactor,peakR*ovsFactor]
                        # TODO: average:
                        #betaClutt = np.power(np.abs(SLCderamp
                        #                            ),2)/(metadata['beta0']**2)
                        #peakBeta0 = np.mean(betaClutt)
                            
                        peakSLC = SLCovs[peakAz*ovsFactor,peakR*ovsFactor]
                        dAz = float('nan')
                        dR = float('nan')
                    
                    # covert to Radar-Cross-Section (RCS):
                    RCSdB = radarUtils.beta2RCSdB(peakBeta0,metadata)
                
                    # get peak range for later ref.phase calc.:
                    dist = (Range[0]/metadata['RSR'] + 
                            metadata['range0time'])*speedOfLight
                    # TODO: TEMP:
                    #dist = (peakR/metadata['RSR'] + 
                    #        metadata['range0time'])*speedOfLight
                    #dist = (np.round(Range[0])/metadata['RSR'] + 
                    #        metadata['range0time'])*speedOfLight            
                else:
                    boundingBox = radarUtils.getBoundingBox(Azimuth[0],
                                                   Range[0],metadata,
                                                   cropSize)
                    SLCraw = s1Utils.readSLC(file,boundingBox,method='raw',
                                           deramp = False)
                    beta0 = np.power(np.abs(SLCraw),2)/(metadata['beta0']**2)
                    # find simple max:
                    if status:
                        idx = np.unravel_index(np.argmax(beta0, axis=None), 
                                               beta0.shape)
                        if plotFlag > 1:
                            if not 'outDir' in kwargs:
                                raise Exception('You must specify output directory for plots')
                            plotUtils.plotBeta0(SLCraw,metadata,
                                                kwargs['outDir'],ovs=False)                                                
                    else:
                        idx = (int(np.round(Azimuth[0])-boundingBox[0][0]),
                               int(np.round(Range[0])-boundingBox[1][0]))
                    peakSLC = SLCraw[idx]
                    RCSdB = radarUtils.beta2RCSdB(beta0[idx],metadata)    
                    peakAz = boundingBox[0][0] + idx[0]
                    peakR = boundingBox[1][0] + idx[1]
                    dAz = (Azimuth[0]-peakAz)*metadata['azimuthSpacing']
                    dR = (Range[0]-peakR)*metadata['rangeSpacing']
                    dist = (Range[0]/metadata['RSR'] + 
                            metadata['range0time'])*speedOfLight
                    # TODO:
                    #dist = (peakR/metadata['RSR'] + 
                    #        metadata['range0time'])*speedOfLight
                # zero-Doppler acq.time:
                acqTime = geoUtils.sod2UTC(radarUtils.radar2time(peakAz,peakR,
                                                                 metadata)[0])
                # populate data dictionary:
                data[slcIdx] = (acqDate,acqTime,status,RCSdB,peakSLC,peakAz,
                                peakR,dAz,dR,dist,np.nan)                    
                slcIdx += 1
            
            # get zenith-azimuth-slantRange:
            _,_,satVec = radarUtils.xyz2t(xyz,metadata)
            zas = geoUtils.xyz2zas(xyz,satVec)
            # create timeseries structure:
            ts = {'id': stack.id, 'data': data, 'metadata': metadata, 
                  'stack': stack,'type':'raw','zas': zas}
            self.stacks.append(ts)


    def getRCSintegral(self, stackId,ovsFactor=32,**kwargs):
        """TODO
        """
        plh = self.ascending.plh
        crs = self.ascending.crs
        stackIdx = self.getStackIdx(stackId)
        stack = self.stacks[stackIdx]['stack']
        #data = self.stacks[stackIdx]['data']
        
        n = len(stack.files)
        RCSall = np.empty(n)
        SCRall = np.empty(n)
        idx = 0
        for file,acqDate in zip(stack.files,stack.acqDates):
            # get status:
            status = (acqDate > self.startDate and acqDate < self.endDate)
            if status:
                metadata = stack.metadata[acqDate]
                (Azimuth,Range) = radarUtils.radarcode(plh, metadata, crs=crs)                  
                # extract larger crop to avoid border effects during oversampling:
                boundingBox = radarUtils.getBoundingBox(Azimuth[0],
                                               Range[0],metadata,10)
                # read SLC:
                SLCderamp = s1Utils.readSLC(file,boundingBox,method='raw',
                                          deramp = True)
                
                RCSdB,SCR = radarUtils.RCSintegral(SLCderamp,metadata,32)
            else:
                RCSdB = np.nan
                SCR = np.nan
            RCSall[idx] = RCSdB
            SCRall[idx] = SCR            
            idx += 1
            
        return RCSall,SCRall


    def removeStack(self,stackId):
        stackIdx = self.getStackIdx(stackId)
        del self.stacks[stackIdx]


    def updateStack(self, stack, ovsFactor=32, posFlag=1, **kwargs):
        try:
            stackIdx = self.getStackIdx(stack.id)
        except:
            return
        # find out update dates:
        datesProcessed = list(self.stacks[stackIdx]['data']['acqDate'])
        datesUpdate = [q for q in stack.acqDates 
                       if q not in datesProcessed]
        for acqDate in datesUpdate:
            file = stack.files[stack.acqDates.index(acqDate)]
            status = (acqDate > self.startDate and acqDate < self.endDate)
            # radarcode & crop:
            plh,crs = self.get_plh(stack.orbit)
            if stack.type == 'raw':
                if posFlag:
                    cropSize = rawCropSizePrecise
                else:
                    cropSize = rawCropSizeApprox
                metadata = stack.metadata[acqDate]
                if 'atmoDir' in kwargs:
                    (Azimuth,Range) = radarUtils.radarcode(plh, metadata,
                                                           atmoDir=kwargs['atmoDir'],
                                                           crs=crs)
                else:
                    (Azimuth,Range) = radarUtils.radarcode(plh, metadata,
                                                           crs=crs)
                if ovsFactor > 1:
                    boundingBox = radarUtils.getBoundingBox(Azimuth[0],
                                            Range[0],metadata,cropSize)
                    SLCderamp = s1Utils.readSLC(file,boundingBox,method='raw',
                                              deramp = True)   
                    SLCovs = radarUtils.oversample(SLCderamp,ovsFactor)
                    beta0 = np.power(np.abs(SLCovs),2)/(metadata['beta0']**2)
                    if status:
                        (peakAz,peakR,peakBeta0,localAz,localR
                         ) = radarUtils.estimatePeak(
                            beta0,boundingBox,metadata,ovsFactor,method='fit')
                        dAz = (Azimuth[0]-peakAz)*metadata['azimuthSpacing']
                        dR = (Range[0]-peakR)*metadata['rangeSpacing']
                        dPhase_ovs = s1Utils.getDerampDemodPhase(
                            file,boundingBox,ovsFactor)
                        dPhase_interp = radarUtils.interp2d(dPhase_ovs,
                                    localAz,localR)
                        peakSLC = s1Utils.reramp(
                            SLCovs[int(localAz),int(localR)],dPhase_interp)
                    else:
                        idx = (int(np.round(Azimuth[0])-boundingBox[0][0]),
                               int(np.round(Range[0])-boundingBox[1][0]))
                        peakAz = idx[0]
                        peakR = idx[1]
                        peakBeta0 = beta0[peakAz*ovsFactor,peakR*ovsFactor]
                        peakSLC = SLCovs[peakAz*ovsFactor,peakR*ovsFactor]
                        dAz = float('nan')
                        dR = float('nan')
                else:
                     boundingBox = radarUtils.getBoundingBox(Azimuth[0],
                         Range[0],metadata,cropSize)
                     SLCraw = s1Utils.readSLC(file,boundingBox,method='raw',
                                              deramp = False)   
                     beta0 = np.power(np.abs(SLCraw),2)/(metadata['beta0']**2)
                     if status:
                         idx = np.unravel_index(np.argmax(beta0, axis=None), 
                                                beta0.shape)
                     else:
                         idx = (int(np.round(Azimuth[0])-boundingBox[0][0]),
                                int(np.round(Range[0])-boundingBox[1][0]))
                     peakSLC = SLCraw[idx]
                     peakBeta0 = beta0[idx] 
                     peakAz = boundingBox[0][0] + idx[0]
                     peakR = boundingBox[1][0] + idx[1]
                     dAz = (Azimuth[0]-peakAz)*metadata['azimuthSpacing']
                     dR = (Range[0]-peakR)*metadata['rangeSpacing']
                    
                RCSdB = radarUtils.beta2RCSdB(peakBeta0,metadata)
                dist = (Range[0]/metadata['RSR'] + 
                            metadata['range0time'])*speedOfLight
                acqTime = geoUtils.sod2UTC(radarUtils.radar2time(peakAz,peakR,
                                                                 metadata)[0])
                # append to existing data structure:
                self.stacks[stackIdx]['data'] = np.append(
                    self.stacks[stackIdx]['data'],
                    np.array((acqDate,acqTime,status,RCSdB,peakSLC,peakAz,peakR,dAz,dR,
                     dist,np.nan), dtype=stackType))
                self.stacks[stackIdx]['stack'] = stack
            else:
                if posFlag:
                    cropSize = coregCropSizePrecise
                else:
                    cropSize = coregCropSizeApprox
                metadata = stack.metadata[acqDate]
                masterMetadata = stack.masterMetadata
                if 'atmoDir' in kwargs:
                    (Azimuth,Range) = radarUtils.radarcode(plh, masterMetadata,
                                                           atmoDir=kwargs['atmoDir'],
                                                           crs=crs)
                else:
                    (Azimuth,Range) = radarUtils.radarcode(plh, masterMetadata,
                                                           crs=crs)
                boundingBox = radarUtils.getBoundingBox(Azimuth[0],
                    Range[0],masterMetadata,cropSize)
                SLC = s1Utils.readSLC(file,boundingBox,method='coreg',
                                      deramp = False)
                IFG = s1Utils.readIFG(file,boundingBox)
                beta0 = np.power(np.abs(SLC),2)/(metadata['beta0']**2)
                idx = np.unravel_index(np.argmax(beta0, axis=None), 
                                       beta0.shape)
                # peak measurement:
                peakSLC = SLC[idx]
                peakIFG = IFG[idx]
                RCSdB = radarUtils.beta2RCSdB(beta0[idx],metadata)
                # coord. differences [m]
                peakAz = boundingBox[0][0] + idx[0]
                peakR = boundingBox[1][0] + idx[1]
                dAz = (Azimuth[0]-peakAz)*metadata['azimuthSpacing']
                dR = (Range[0]-peakR)*metadata['rangeSpacing']
                # peak range:
                # TODO: THis dist is wrong for coreg!
                dist = (peakR/metadata['RSR'] + metadata['range0time'])*speedOfLight
                # zero-Doppler acq.time:
                tAzimuth,_ = radarUtils.radar2time(peakAz,peakR,metadata)
                acqTime = geoUtils.sod2UTC(tAzimuth)
                # append to existing data structure:
                self.stacks[stackIdx]['data'] = np.append(
                    self.stacks[stackIdx]['data'],
                    np.array((acqDate,acqTime,status,RCSdB,peakSLC,peakAz,peakR,dAz,dR,
                     dist,peakIFG), dtype=stackType))
                self.stacks[stackIdx]['stack'] = stack
        

    def getStackIdx(self,stackId):
        try:
            idx = [stack['id'] for stack in self.stacks].index(stackId)
            return idx
        except:
            print('No such stack')
    
    def detectOutliers(self):
        """Detect outlying samples in all stack timeseries using 
        3*MAD threshold
        """
        if self.stacks: # test if stacks timeseries already extracted
            for stack in self.stacks:
                status = stack['data']['status']
                RCS = stack['data']['RCS']
                RCS_on = stack['data']['RCS'][status]
                lowerThr = np.median(RCS_on)-3*ss.median_absolute_deviation(RCS_on)
                stack['goodIdx'] = np.nonzero((RCS > lowerThr) & status)[0]
        else:
            print('No stack timeseries extracted yet.')
    
    def RCSanalysis(self):
        """Estimate RCS and SCR of reflector as well as clutter using 
        all extracted stack timeseries
        """
        for stack in self.stacks:
            RCS = stack['data']['RCS']
            metadata = stack['metadata']
            # estimate clutter from samples before:
            # TODO: refine this on STATUS!
            clutIdx = np.where(stack['data']['acqDate'] < self.installDate)[0]
            if len(clutIdx) > 4:
                ampl = radarUtils.RCSdB2amp(RCS[clutIdx],metadata)
                # fit Rayleigh distribution:
                mleRayleigh = ss.rayleigh.fit(ampl)
                #stack['clutRCSbefore'] = radarUtils.amp2RCSdB(mleRayleigh[1],metadata)
                stack['clutRCSbefore'] = radarUtils.beta2RCSdB(2*mleRayleigh[1]**2,
                                                               metadata)
            else:
                stack['clutRCSbefore'] = np.nan
            # estimate reflector RCS and SCR:
            ampl = radarUtils.RCSdB2amp(RCS[stack['goodIdx']],metadata)
            # use MLE fitting only for > 20 time series:
            if len(ampl) > 20:
                # fit Rice distribution:
                mleRice = ss.rice.fit(ampl, np.mean(ampl), loc=0,
                                      scale=np.std(ampl)) # add initial guess as constraint
                # get Rice confidence interval
                interval = ss.rice.interval(0.95,mleRice[0],loc=mleRice[1],
                                            scale=mleRice[2])
                RCSinterval = radarUtils.amp2RCSdB(interval,metadata)
                sigRCSmle = np.abs((RCSinterval[1]-RCSinterval[0])/2/2.5)
                # estimate SCR from parameters of Rice distrib.
                #stack['clutRCSafter'] = radarUtils.amp2RCSdB(mleRice[2],metadata) # sigma (scale) parameter
                stack['clutRCSafter'] = radarUtils.beta2RCSdB(2*mleRice[2]**2,
                                                              metadata)
                v = mleRice[0]*mleRice[2]+mleRice[1] # shape parameter of Rice
                stack['reflRCS'] = radarUtils.amp2RCSdB(v,metadata) # v parameter
                #stack['SCR'] = stack['reflRCS']-stack['clutRCSafter']
                stack['SCR'] = 10*np.log10(v**2/(2*(mleRice[2]**2)))
                stack['sigRCS'] = sigRCSmle
            else:
                stack['reflRCS'] = np.nanmean(RCS[stack['goodIdx']])
                stack['sigRCS'] = np.nanstd(RCS[stack['goodIdx']])
                stack['clutRCSafter'] = np.nan
                if np.isnan(stack['clutRCSbefore']):
                    stack['SCR'] = np.nan
                else:
                    stack['SCR'] = stack['reflRCS'] - stack['clutRCSbefore']
            # estimate CRLB predictions:    
            stack['D_A'] = np.std(ampl)/np.mean(ampl)
            stack['sigPhi_SCR'] = (1/np.sqrt(np.power(10,stack['SCR']/10)*2)
                                   /4/np.pi*metadata['wavelength'])
            stack['sigPhi_DA'] = stack['D_A']/4/np.pi*metadata['wavelength']
            # get analytical RCS:
            self.getRCS0(stack['id'])
    
    def getStackStats(self,stackId):
        """Return dictionary containg timeseries statistics for specified 
        stackId
        """    
        stackIdx = self.getStackIdx(stackId)
        stack = self.stacks[stackIdx]
        if not 'SCR' in stack.keys():
            self.detectOutliers()
            self.RCSanalysis()  
        stats = dict((k, stack[k]) for k in ['clutRCS','reflRCS','SCR','sigRCS',
                     'clutRCSbefore','clutRCSafter','sigPhi_DA','sigPhi_SCR'] 
                     if k in stack)
        #calculate Cramer-Rao lower bound for positioning:
        #SCR in linear
        SCR = np.power(10,stats['SCR']/10)
        rangeCRB = np.sqrt(3/2/(np.pi**2)/SCR)*stack['metadata']['rangeResolution']
        azimuthCRB = np.sqrt(3/2/(np.pi**2)/SCR)*stack['metadata']['azimuthResolution']
        rangeSTD = np.std(stack['data']['dR'][stack['goodIdx']])
        azimuthSTD = np.std(stack['data']['dAz'][stack['goodIdx']])
        rangeALE = np.mean(stack['data']['dR'][stack['goodIdx']])
        azimuthALE = np.mean(stack['data']['dAz'][stack['goodIdx']])
        stats['azimuthCRB'] = azimuthCRB
        stats['rangeCRB'] = rangeCRB
        stats['azimuthSTD'] = azimuthSTD
        stats['rangeSTD'] = rangeSTD
        stats['azimuthALE'] = azimuthALE
        stats['rangeALE'] = rangeALE
        stats['RCS0'] = stack['RCS0']
        return stats
    
    def getDphase(self,stackId,masterDate):
        """Return numpy array of differential interferometric phase timeseries
        for specified stackId using specified masterDate as reference
        Also returns timeseries dates in 'YYYYMMDD' format
        """ 
        stackIdx = self.getStackIdx(stackId)
        goodIdx = self.stacks[stackIdx]['goodIdx']
        dates = self.stacks[stackIdx]['data']['acqDate'][goodIdx]
        masterIdx = np.where(dates == masterDate)[0][0]
        SLC = self.stacks[stackIdx]['data']['SLC'][goodIdx]
        # make interferogram w.r.t. master:
        ifg = SLC[masterIdx]*np.conj(SLC)
        # get reference phase:
        dist = self.stacks[stackIdx]['data']['dist'][goodIdx]
        wavelength = self.stacks[stackIdx]['metadata']['wavelength']
        phi0 = np.exp(-1j*(4*np.pi/wavelength * (dist[masterIdx] - dist)))
        ifgCorr = ifg*np.conj(phi0)
        return ifgCorr,dates

    def getDphase_coreg(self,stackId,masterDate):
        """Return numpy array of differential interferometric phase timeseries
        for specified stackId of 'coreg' type (using its master as reference)
        """ 
        stackIdx = self.getStackIdx(stackId)
        stack = self.stacks[stackIdx]
        # TODO: for now use status instead of goodIdx
        goodIdx = self.stacks[stackIdx]['data']['status']
        dates = stack['data']['acqDate'][goodIdx]
        masterIdx = np.where(dates == masterDate)[0][0]
        SLC = stack['data']['SLC'][goodIdx]
        # make interferogram w.r.t. master:
        ifg = SLC[masterIdx]*np.conj(SLC)
        # get reference phase:
        #az = np.mean(stack['data']['azimuth'][stack['goodIdx']])
        #r = np.mean(stack['data']['range'][stack['goodIdx']])
        az = stack['data']['azimuth'][masterIdx]
        r = stack['data']['range'][masterIdx]
        masterMetadata = stack['stack'].metadata[masterDate]
        plh,_ = self.get_plh(stack['metadata']['orbit'])
        xyz = radarUtils.radar2xyz(az,r,plh[2],masterMetadata)
        tAz_master, tR_master = radarUtils.radar2time(az,r,masterMetadata)
        dDist = np.empty(len(dates))
        idx = 0
        for acq in dates:
            metadata = stack['stack'].metadata[acq]
            tA,tR,satVec = radarUtils.xyz2t(xyz,metadata)
            dDist[idx] = (tR_master - tR)*speedOfLight
            idx += 1            
        wavelength = self.stacks[stackIdx]['metadata']['wavelength']
        phi0 = np.exp(-1j*(4*np.pi/wavelength * dDist))
        ifgCorr = ifg*np.conj(phi0)
        return ifgCorr#,dist
    

    def getBtemp(self,stackId,masterDate = None):
        """Return temporal baselines in days for specified stackId 
        and masterDate
        """ 
        stackIdx = self.getStackIdx(stackId)
        #goodIdx = self.stacks[stackIdx]['goodIdx']
        goodIdx = self.stacks[stackIdx]['data']['status']
        acqDates = self.stacks[stackIdx]['data']['acqDate'][goodIdx]
        Btemp = [datetime.datetime.strptime(d,'%Y%m%d').date() for d in acqDates]
        if self.stacks[stackIdx]['type'] == 'coreg':
            masterDate = self.stacks[stackIdx]['metadata']['acqDate'].date()
        elif (self.stacks[stackIdx]['type'] == 'raw') and (masterDate == None):
            print('For raw stacks masterDate must be specified!')
            return
        masterIdx = Btemp.index(masterDate)
        Btemp = np.array([(temp-masterDate).days for temp in Btemp])
        return Btemp,acqDates,masterIdx


    def plotRCS(self,outDir):
        """Plot RCS timeseries for all extracted stack timeries and save 
        as png in specified outDirectory
        """       
        symbolGenerator = plotUtils.symbolGenerator()
        colorGenerator = plotUtils.colorGenerator()
        installDate = dates.datestr2num(self.installDate)
        startDate = dates.datestr2num(self.startDate)
        fig, ax = plt.subplots(figsize=(9,4),dpi=120)        
        anotText = []
        anotColor = []
        for stack in self.stacks:
            color = next(colorGenerator)
            #plotDates = dates.datestr2num(stack['stack'].acqDates)
            plotDates = dates.datestr2num(stack['data']['acqDate'])
            ax.plot(plotDates,stack['data']['RCS'],next(symbolGenerator),
                    color=color,markersize = 2.5,label = stack['id'])        
            # plot estimated reflector RCS:
            ax.hlines(stack['reflRCS'],startDate,np.max(plotDates),
                      color=color,linewidth=1,linestyle = '-.')
            # plot estimated clutter RCS:
            ax.hlines(stack['clutRCSafter'],startDate,np.max(plotDates),
                      color=color,linewidth=1,linestyle = '-.')
            ax.hlines(stack['clutRCSbefore'],np.min(plotDates),startDate,
                      color=color,linewidth=1,linestyle = '-.')            
            ax.annotate('RCS',xy=(np.max(plotDates),stack['reflRCS']-1),
                fontsize = 10,color = color)   
            ax.annotate('clutt',xy=(np.max(plotDates),stack['clutRCSafter']),
                fontsize = 10,color = color)               
            SCR = stack['SCR']
            anotText.append('SCR = ' +f"{SCR:.1f} dB")
            anotColor.append(color)
        
        if stack['RCS0'] is not None:
            RCS0 = np.mean([stack['RCS0'] for stack in self.stacks])
            # plot estimated reflector RCS:
            ax.hlines(RCS0,startDate,np.max(plotDates),
                      color='k',linewidth=0.7,linestyle = '--')
            ax.annotate('RCS0',xy=(np.max(plotDates),RCS0),
                        fontsize = 10,color = 'k')        
        
        # annotate SCR:
        anotPos = stack['reflRCS']
        for txt,col in zip(anotText,anotColor):
            ax.annotate(txt,xy=(np.min(plotDates),anotPos),
                        fontsize = 12,color = col)
            anotPos -= 5
        
        ax.set_ylabel('Apparent RCS [$dBm^2$]')
        start, end = ax.get_xlim()
        stepsize = (end+30-start)/10
        ax.xaxis.set_ticks(np.arange(start, end+30, stepsize))
        ax.xaxis.set_major_formatter(formatter)
        ax.xaxis.set_tick_params(rotation=30, labelsize=10)
        ax.axvline(x = installDate,linestyle=':', color='k')
        ax.axvline(x = startDate,linestyle='--', color='k')
        ax.grid(alpha = 0.3)        
        ax.set_title(self.id,fontsize = 12)
        ax.legend(loc = 'lower right')
        plt.subplots_adjust(hspace=0.5)
        plt.savefig(outDir+os.sep+self.id+'_RCS_TS.png', bbox_inches='tight')
        plt.close()
    
    
    def plotALE(self,stackId,outDir,atmoDir=''):
        """Plot Absolute Location Errors (ALE) for timeseries of specified 
        stackId and save as png in outDirectory
        """ 
        stackIdx = self.getStackIdx(stackId)
        if stackIdx is None:
            return
        stack = self.stacks[stackIdx]     
        if self.ascending != None:
            crs = self.ascending.crs
        else:
            crs = self.descending.crs
        # prepare figure
        fig, axes = plt.subplots(figsize=(8,8))
        axes.axhline(color='k')
        axes.axvline(color='k')        
        startIdx = np.argmax(stack['data']['status'])        
        legendCounter = True
        dAz = np.zeros(len(stack['data']))
        dAz[:] = np.nan
        dR = np.zeros(len(stack['data']))
        dR[:] = np.nan
        # loop through timeseries
        for acqDate in stack['data']['acqDate'][startIdx:]:
            idx = np.where(np.char.find(stack['data']['acqDate'],acqDate)
                        == 0)[0]
            if idx in stack['goodIdx']:
                metadata = stack['stack'].metadata[acqDate]                
                # orbit-wise ECEF apex coordinates:
                xyz = self.get_xyz(stack['metadata']['orbit'])                
                # radar time coordinates:
                (tAzimuth_ETRS,tRange_ETRS,satVec
                 ) = radarUtils.xyz2t(xyz,metadata)                
                # transform from ETRS-89 to ITRS (ITRF2014), t = acqEpoch
                decimalDate = geoUtils.decimalYear(metadata['acqDate'])
                if crs == 'ITRS':
                    xyz = geoUtils.itrf2itrf(xyz,decimalDate)
                else:
                    xyz = geoUtils.etrf2itrf(xyz,decimalDate)
                (tAzimuth_ITRS,tRange_ITRS,satVec
                 ) = radarUtils.xyz2t(xyz,metadata)
                # get SET:
                decimalHour = (metadata['acqDate'].hour
                                + metadata['acqDate'].minute/60 
                                + metadata['acqDate'].second/3600)
                dxyz_set = geoUtils.getSET(xyz,acqDate,decimalHour)
                xyz += dxyz_set
                (tAzimuth_SET,tRange_SET,satVec
                 ) = radarUtils.xyz2t(xyz,metadata)
                # get tropo and iono delay:
                if atmoDir:
                    slantDelay = atmoUtils.getTropoDelay(xyz, satVec, acqDate, 
                                               atmoDir, method='Jolviet')/speedOfLight
                    ionoDelay = atmoUtils.getIonoDelay(xyz, satVec[0], 
                                   metadata['acqDate'], atmoDir=atmoDir)/speedOfLight
                else:                
                    slantDelay = geoUtils.tropoDelay(xyz,satVec)
                    ionoDelay = 0.1/speedOfLight
                #print(slantDelay*speedOfLight)
                slantDelay += 0.15/speedOfLight # TODO
                tRange_TROPO = tRange_SET + slantDelay                    
                tRange_IONO = tRange_TROPO + ionoDelay                
                # S1 residual bistatic correction:
                bistaticAz = s1Utils.bistaticCorrection(tRange_IONO,metadata)
                tAzimuth_BISTATIC = tAzimuth_SET - bistaticAz                    
                # S1 Doppler shift
                doppler = s1Utils.dopplerRgCorrection(tRange_IONO,
                                                    tAzimuth_BISTATIC,metadata)
                tRange_Doppler = tRange_IONO - doppler 
                # S1 FM rate mismash
                FMmismash = s1Utils.FMmismatchCorrection(xyz,tRange_Doppler,
                                                 tAzimuth_BISTATIC,metadata)
                tAzimuth_FM = tAzimuth_BISTATIC + FMmismash
                
                # get IRF peak:
                Azimuth = stack['data']['azimuth'][idx]
                Range = stack['data']['range'][idx]
                tAzimuth, tRange = radarUtils.radar2time(Azimuth,Range,metadata)                
                # convert corrections to metres:
                rangeFactor = metadata['RSR']*metadata['rangeSpacing']
                azimuthFactor = metadata['PRF']*metadata['azimuthSpacing']
                dAz_ETRS = (tAzimuth_ETRS - tAzimuth)*azimuthFactor
                dR_ETRS = (tRange_ETRS - tRange)*rangeFactor
                dAz_ITRS = (tAzimuth_ITRS - tAzimuth)*azimuthFactor
                dR_ITRS = (tRange_ITRS - tRange)*rangeFactor
                dAz_SET = (tAzimuth_SET - tAzimuth)*azimuthFactor
                dR_SET = (tRange_SET - tRange)*rangeFactor
                dAz_TROPO = dAz_SET
                dR_TROPO = (tRange_TROPO - tRange)*rangeFactor
                dAz_IONO = dAz_SET
                dR_IONO = (tRange_IONO - tRange)*rangeFactor
                dAz_BIST = (tAzimuth_BISTATIC - tAzimuth)*azimuthFactor
                dR_BIST = dR_IONO
                dAz_DOPPLER = dAz_BIST
                dR_DOPPLER = (tRange_Doppler - tRange)*rangeFactor
                dAz_FM = (tAzimuth_FM - tAzimuth)*azimuthFactor
                dR_FM = dR_DOPPLER
                # plot all:
                symb = plotUtils.symbolGenerator()
                col = plotUtils.colorGenerator()
                axes.plot(dR_ETRS,dAz_ETRS,next(symb)+next(col),
                          label='Initial' if legendCounter else "")
                axes.plot(dR_ITRS,dAz_ITRS,next(symb)+next(col),
                          label='+ plate m.' if legendCounter else "")
                axes.plot(dR_SET,dAz_SET,next(symb)+next(col),
                          label='+ SET' if legendCounter else "")
                axes.plot(dR_TROPO,dAz_TROPO,next(symb)+next(col),
                          label='+ tropo' if legendCounter else "")
                axes.plot(dR_IONO,dAz_IONO,next(symb)+next(col),
                          label='+ iono' if legendCounter else "")
                axes.plot(dR_BIST,dAz_BIST,next(symb)+next(col),
                          label='+ bistatic' if legendCounter else "")
                axes.plot(dR_DOPPLER,dAz_DOPPLER,next(symb)+next(col),
                          label='+ Doppler' if legendCounter else "")
                #axes.plot(dR_FM,dAz_FM,next(symb)+next(col),
                #          label='+ FM-rate' if legendCounter else "")
                axes.plot(dR_FM,dAz_FM,next(symb)+'k',markersize=10,
                          label='+ FM-rate' if legendCounter else "")
                #axes.plot(0,0,next(symb)+'k',markeredgewidth=5,markersize=20,
                #          label='RP' if legendCounter else "")
                legendCounter = False
                # overwrite dAz/dR:
                dAz[idx] = dAz_FM
                dR[idx] = dR_FM
            
        axes.tick_params(axis='both', which='major', direction='in', 
                         length = 10, labelsize=14)
        axes.set_ylabel(r'$\Delta$ Azimuth [m]',fontsize=14)
        axes.set_xlabel(r'$\Delta$ Range [m]',fontsize=14)
        axes.legend(fontsize = 16)
        axes.set_title(self.id+'_ALE_'+stackId)
        
        plt.subplots_adjust(hspace=0.5)
        plt.savefig(outDir+os.sep+self.id+'_ALE_'+stack['id']+'.png', 
                    bbox_inches='tight')
        plt.close()
        self.stacks[stackIdx]['data']['dAz'] = dAz
        self.stacks[stackIdx]['data']['dR'] = dR
        return dAz, dR


    def plotALE_TS(self,outDir):

        colGenerator = plotUtils.colorGenerator()
        symbGenerator = plotUtils.symbolGenerator()
        
        fig, axes = plt.subplots(2,1,figsize=(10,8),dpi=120)
        
        for stack in self.stacks:
            color = next(colGenerator)
            symb = next(symbGenerator)
            goodIdx = stack['goodIdx']
            dr = stack['data']['dR'][goodIdx]
            da = stack['data']['dAz'][goodIdx]
            acqDates = stack['data']['acqDate'][goodIdx]
            plotDates = dates.datestr2num(acqDates)
            axes[0].plot(plotDates,dr,symb,markersize=8, color=color, 
                     label = stack['id'])
            axes[0].hlines(np.mean(dr),np.min(plotDates),np.max(plotDates),
                  color=color,linewidth=1,linestyle = '-')
            axes[1].plot(plotDates,da,symb,markersize=8, color=color, 
                     label = stack['id'])
            axes[1].hlines(np.mean(da),np.min(plotDates),np.max(plotDates),
                  color=color,linewidth=1,linestyle = '-')
        
        axes[0].set_ylabel('$\Delta$ range [m]',fontsize=15)
        axes[1].set_ylabel('$\Delta$ azimuth [m]',fontsize=15)

        for ax in axes:
            start, end = ax.get_xlim()
            stepsize = (end+50-start)/7
            ax.xaxis.set_ticks(np.arange(start, end+60, stepsize))
            ax.xaxis.set_major_formatter(formatter)
            ax.xaxis.set_tick_params(rotation=30, labelsize=13)
            ax.yaxis.set_tick_params(labelsize=13)
            ax.grid(alpha = 0.3)   
            ax.legend(fontsize=15,loc='right')
        axes[0].set_xticklabels([])
        axes[0].set_title(self.id, fontsize=19, fontweight='bold')
        plt.subplots_adjust(hspace=0.1,wspace=0.5)
        plt.savefig(outDir+os.sep+self.id+'_ALE_TS.png', bbox_inches='tight')
        plt.close()


    def toJSON(self,outDir):
        """Serialize to JSON file in outDirectoy
        """         
        outFile = outDir + os.sep + self.id + '.json'
        with open(outFile, "w") as writeFile:
            json.dump(self, writeFile, cls = crEncoder,indent=2)


    def statsToJSON(self,outDir):
        """All stack timeseries statistics exported to JSON in outDirectory
        """         
        outFile = outDir + os.sep + self.id + '_stats.json'    
        stats = {}
        for stack in self.stacks:
            stats[stack['id']] = self.getStackStats(stack['id'])
        
        with open(outFile, "w") as write_file:
            json.dump(stats, write_file,indent=2)
    
    
    def get_all_timeseries(self):
        """Return string report of all timeseries
        """ 
        outStr = []
        for stack in self.stacks:
            for data in stack['data']:
                acqDate = data['acqDate'] + ' ' + data['acqTime']
                incAngle = stack['zas'][0]*180/np.pi
                RCS = data['RCS']
                dAz = data['dAz']
                dR = data['dR']
                outStr.append(acqDate + ' ' + stack['id'] + ' ' + 
                          'IW'+str(stack['metadata']['swath'])
                          + f" {incAngle:.1f}" + f" {RCS:.1f}" + f" {dAz:.2f}" 
                          + f" {dR:.2f}")
        return outStr


    def print_all_timeseries(self,outDir):
        """Write .txt file containing all timeseries into outDirectory
        """ 
        with open(outDir + os.sep + self.id + '_timeseries.txt', "w") as wf:
            wf.write('Station id: ' +self.id+ '\n')
            wf.write('acqDate acqTime stackId subswath incAngle RCS dAz dR\n')
            #print(*self.get_all_timeseries(),sep = "\n")
            for L in self.get_all_timeseries():
                wf.write(L + '\n')


    def add_insar_ts(self, insarHDF, stackId):
        stack_idx = self.getStackIdx(stackId)
        
        try:
            data = insarUtils.openHDF(insarHDF)
            
            #if 'network3' in data:
            #    network_label = 'network3'
            #    psc_label = 'psc3'
            #else:
            network_label = 'network2'
            psc_label = 'psc2'                

            
            plh_all = data[psc_label + '/plh'][:]
            ts_all = (data[network_label + '/ps_displ'][:] /
                      data[network_label].attrs['m2ph']*1e3)
            var_all = data[network_label+'/ps_var'][:]
            refIdx = data[network_label].attrs['refIdx']
            if ~np.isnan(refIdx):
                ts_all = np.insert(ts_all, refIdx, 0, axis=0)
                var_all = np.insert(var_all, refIdx, 0)
            dates_all = data[network_label].attrs['dates']
            std_res_all = data[network_label+'/stdRes'][:]
            VC_mean = (np.sqrt(data[network_label+'/VC_mean'][:]) /
                       -data[network_label].attrs['m2ph']*1e3)
            
            if 'network3' in data:
                network2 = data['network3']
                psc2 = data['psc3']
                plh_all = np.concatenate((plh_all, psc2['plh'][:]))
                ts_all = np.concatenate((ts_all, 
                                network2['ps_displ'][:]/
                                data[network_label].attrs['m2ph']*1e3))
                var_all = np.append(var_all,
                                    network2['ps_var'][:])
                std_res_all = np.append(std_res_all,
                                        network2['stdRes'][:])
                
            data.close()
        except:
            print('No InSAR HDF data found or specified HDF'+
                  ' does not contain the InSAR network solution.')
            raise
        # extract TS:
        plh = self.get_plh()[0]
        # take index of closest PS:
        idx = np.argmin(cdist(plh_all[:,:2], plh[None,:2]))
        ts = ts_all[idx,:]
        # TODO: add warning if too far
        good_idx = np.intersect1d(
            dates.datestr2num(
                self.stacks[stack_idx]['data']['acqDate'][
                    self.stacks[stack_idx]['goodIdx']
                    ]), 
            dates_all, return_indices=True)[-1]
        # remove outliers:
        plot_dates = dates_all[good_idx]
        plot_ts = ts[good_idx]
        VC = VC_mean[good_idx]
        # stds:
        sig = np.sqrt(var_all[idx])*VC
        std_res = std_res_all[idx]
        insar_ts = {
            'ts': plot_ts,
            'dates' : plot_dates,
            'sigma_ts' : sig,
            'std' : std_res}
        self.stacks[stack_idx]['insar'] = insar_ts
    
    
    def plot_insar_ts(self, out_dir = ''):
        colorGenerator = plotUtils.colorGenerator()
        fig, ax = plt.subplots(figsize=(9,4), dpi=120)
        for stack in self.stacks:
            if 'insar' in stack.keys():
                insar = stack['insar']
                # plot ts with errorbars:
                ax.errorbar(insar['dates'], 
                            insar['ts'], fmt='-^', 
                            yerr = 2.5*insar['sigma_ts'],
                            ecolor = 'gray', label = stack['id'], 
                            color = next(colorGenerator))
        ax.set_ylabel('LOS displacement [mm]', fontsize = 14)
        ax.yaxis.set_tick_params(labelsize=14)
        start, end = ax.get_xlim()
        stepsize = (end+30-start)/10
        ax.xaxis.set_ticks(np.arange(start, end+30, stepsize))
        ax.xaxis.set_major_formatter(formatter)
        ax.xaxis.set_tick_params(rotation=30, labelsize=12)
        #ax.set_ylim(-50,50)
        ax.grid(alpha = 0.35) 
        ax.set_title(self.id, fontsize = 14, fontweight='bold')
        ax.legend(loc = 'best', fontsize=14)
        if out_dir:
            out_fig = (out_dir + os.sep + self.id + '_' + stack['id'] 
                       + '_insar.png')
            plt.savefig(out_fig, bbox_inches='tight')
            plt.close()


    def decomp_insar_ts(self, out_dir = ''):
        """compute vertical/horizontal displacement ts from asc/dsc los displ.
        """
        stacks = self.stacks
        if 'insar' not in stacks[0].keys():
            print('perform InSAR analysis first')
            raise
        # take shortest stack as reference:
        ref_idx = np.argmin([len(s['insar']['dates']) for s in stacks])
        ref_dates = stacks[ref_idx]['insar']['dates']
        # allocate:
        n_stacks = len(stacks)
        n_dates = len(ref_dates)
        los = np.zeros((n_dates, n_stacks))
        sig_los = np.zeros_like(los)
        inc = np.zeros(n_stacks)
        head = np.zeros(n_stacks)
        # prepare decomposition arrays
        for i in np.arange(n_stacks):
            d = stacks[i]['insar']['dates']
            # find closest dates:
            d_idx = np.argmin(np.abs(d[:, np.newaxis] - ref_dates), axis=0)
            los[:,i] = stacks[i]['insar']['ts'][d_idx]
            sig_los[:,i] = stacks[i]['insar']['sigma_ts'][d_idx]
            inc[i] = stacks[i]['zas'][0]
            head[i] = stacks[i]['zas'][1]
        # write decomposition
        vh, sig_vh = insarUtils.decomp(los, sig_los, inc, head)
        self.insar_decomp = {
            'vh' : vh,
            'sig_vh' : sig_vh,
            'dates' : ref_dates}
            
    
    def plot_decomp_ts(self, out_dir = ''):
        try:
            insar_decomp = self.insar_decomp
        except:
            print('Perform InSAR decomposition first.')
            raise
        fig, ax = plt.subplots(figsize=(9,4), dpi=120)
        ax.errorbar(insar_decomp['dates'], 
                    insar_decomp['vh'][0], fmt='-^', 
                    yerr = 2.5*insar_decomp['sig_vh'][0],
                    ecolor = 'gray', label = 'vertical', 
                    color = 'C0')
        ax.errorbar(insar_decomp['dates'], 
                insar_decomp['vh'][1], fmt='-o', 
                yerr = 2.5*insar_decomp['sig_vh'][1],
                ecolor = 'gray', label = 'horizontal (EW)', 
                color = 'C1')
        ax.yaxis.set_tick_params(labelsize=14)
        ax.set_ylabel('Displacement [mm]', fontsize = 14)
        start, end = ax.get_xlim()
        stepsize = (end+30-start)/10
        ax.xaxis.set_ticks(np.arange(start, end+30, stepsize))
        ax.xaxis.set_major_formatter(formatter)
        ax.xaxis.set_tick_params(rotation=30, labelsize=12)
        ax.grid(alpha = 0.35) 
        #ax.set_ylim(-50,50)
        ax.set_title(self.id, fontsize = 14, fontweight='bold')
        ax.legend(loc = 'best', fontsize=14)
        if out_dir:
            out_fig = (out_dir + os.sep + self.id + '_insar_decomp.png')
            plt.savefig(out_fig, bbox_inches='tight')
            plt.close()


class CornerReflector(Reflector):
    
    def __init__(self,reflId,ascending,descending,
                 installDate,startDate,endDate,
                 shape,a0,zenDip,aziDip):
        
        super().__init__(reflId,'CREF',ascending,descending,
                         installDate,startDate,endDate)
        
        self.defGeometry(shape,a0,zenDip,aziDip)

        
    def defGeometry(self,shape,a0,zenDip,aziDip):
        """Specify reflector shape, inner leg-length, zenith and azimuth dip
        """ 
        self.geometry = {
                'shape': shape,
                'a0': a0,
                'zenDip': zenDip,
                'aziDip': aziDip}

    
    def getRCS0(self,stackId):
        """Return analytical RCS of reflector using geometric optics 
        simulation for specified stackId
        """ 
        stackIdx = self.getStackIdx(stackId)
        stack = self.stacks[stackIdx]        
        if not 'RCS0' in stack.keys():
            wavelength = stack['metadata']['wavelength']
            zas = stack['zas']
            if not self.geometry:
                print('First define CR geometry!')
                stack['RCS0'] = np.nan
            else:
                crType = self.geometry['shape']
                a0 = self.geometry['a0']
                zenDip = self.geometry['zenDip']
                aziDip = self.geometry['aziDip']            
                RCS0,dtheta,domega = crUtils.crRCS0(crType,a0,zenDip,aziDip,
                                                    zas,wavelength)
                stack['RCS0'] = RCS0
                stack['dtheta'] = dtheta
                stack['domega'] = domega
        return stack['RCS0']


class Transponder(Reflector):
    
    def __init__(self,reflId,ascending,descending,
                 installDate,startDate,endDate,
                 RCS0):
        
        super().__init__(reflId,'CAT',ascending,descending,
                         installDate,startDate,endDate)
        self.RCS0 = RCS0
    
    def getRCS0(self,stackId):
        """Return expected RCS of transponder
        """
        stackIdx = self.getStackIdx(stackId)
        stack = self.stacks[stackIdx]
        if not 'RCS0' in stack.keys():
            stack['RCS0'] = self.RCS0
        return self.RCS0
    
    def plotRCS_transponder(self,outDir):
                
        symbGenerator = plotUtils.symbolGenerator()
        colGenerator = plotUtils.colorGenerator()
        installDate = dates.datestr2num(self.installDate)
        startDate = dates.datestr2num(self.startDate)
        if self.endDate != '99999999':
            endDate = dates.datestr2num(self.endDate)
        else:
            endDate = np.nan
        # prepare figure:
        fig, ax = plt.subplots(figsize=(9,4),dpi=120)
        
        for stack in self.stacks:
            color = next(colGenerator)
            #plotDates = dates.datestr2num(stack['stack'].acqDates)
            plotDates = dates.datestr2num(stack['data']['acqDate'])
            ax.plot(plotDates,stack['data']['RCS'],next(symbGenerator),
                    color=color,markersize = 2.5,label = stack['id'])
            
            # plot estimated reflector RCS:
            if endDate == np.nan:
                endDate = np.max(plotDates)
            ax.hlines(stack['reflRCS'],startDate,endDate,
                      color=color,linewidth=1,linestyle = '-.')
            ax.hlines(stack['clutRCSbefore'],np.min(plotDates),startDate,
                      color=color,linewidth=1,linestyle = '-.')
                       
        ax.set_ylabel('Apparent RCS [$dBm^2$]')
        start, end = ax.get_xlim()
        #
        #ax.set_xlim(start+90)
        #start, end = ax.get_xlim()
        #
        stepsize = (end+30-start)/10
        ax.xaxis.set_ticks(np.arange(start, end+30, stepsize))
        ax.xaxis.set_major_formatter(formatter)
        ax.xaxis.set_tick_params(rotation=30, labelsize=10)
        ax.axvline(x = installDate,linestyle=':', color='k')
        ax.axvline(x = startDate,linestyle='--', color='k')
        ax.axvline(x = endDate,linestyle='--', color='k')
        ax.grid(alpha = 0.3)        
        ax.set_title(self.id,fontsize = 12)
        ax.legend(loc = 'lower right')
        plt.subplots_adjust(hspace=0.5)
        plt.savefig(outDir+self.id+'_RCS_TS.png', bbox_inches='tight')
        plt.close()

class ReflectorGeometry:
    
    def __init__(self, shape, flipped, a0, 
                 zenDip = 0, aziDip = 0, b0 = 0, aperture = None):
        self.shape = shape
        self.flipped = flipped
        self.a0 = a0
        self.zenDip = zenDip
        self.aziDip = aziDip
        self.b0 = b0
        self.aperture = aperture


class Coordinates:
    
    def __init__(self, latitude, longitude, elevation, 
                 crs = 'ETRS-89', frame = 'ETRF2000', epoch = '2010.0', 
                 epsg = '4258'):
        self.longitude = longitude
        self.latitude = latitude
        self.elevation = elevation
        self.crs = crs
        self.frame = frame
        self.epoch = epoch
        self.epsg = epsg
        self._plh = None
        self._xyz = None
    
    @property
    def plh(self):
        if self._plh is None:
            self._plh = np.array([self.latitude, self.longitude, 
                                  self.elevation])
        return self._plh
    
    @property
    def xyz(self):
        if self._xyz is None:
            # TODO: cope with different CRS:
            self._xyz = geoUtils.plh2xyz(self.plh)
        return self._xyz


class Stack:

    def __init__(self, id, stackDir, sensor, subswath, stackType = 'raw'):
        self.id = id
        self.stackDir = stackDir
        self.sensor = sensor
        self.subswath = subswath
        self.type = stackType
        self._nSLC = None
        self._masterIdx = None
        
    def readData(self):
        if self.sensor == 'Sentinel-1':        
            if self.type == 'raw':
                # list all slaves:
                imgs = glob.glob(self.stackDir + os.sep +"slaves" + os.sep + "*"
                                 +self.subswath+".dim")
                if not imgs:
                    raise Exception('No Sentinel-1 images prepared by SNAP'
                                    +' in specified stack directory')
                print('Reading stack ' +self.id+ ' , ' +
                      str(len(imgs))+ ' SLCs...')
                d = dict()
                for file in imgs:
                    metadata = s1Utils.readMetadata(file)
                    acqStr = metadata['acqDate'].strftime("%Y%m%d")
                    d[acqStr] = metadata.copy()
                self.metadata = OrderedDict(sorted(d.items(), 
                                                   key=lambda t: t[0]))
            elif self.type == 'coreg':
                # determined subswath if none specified:
                if not self.subswath:
                    q = os.listdir(self.stackDir + os.sep +"coreg" + os.sep)
                    if len(q) > 1:
                        #raise Exception('Stack '+self.id+ ' consists of more'
                        #                + ' than single subswaths, specify one.')
                        imgs = glob.glob(self.stackDir + os.sep + "coreg" 
                                         + os.sep + "stack*.dim")
                    else:
                        self.subswath = q[0]   
                else:
                     # list all slaves:
                     imgs = glob.glob(self.stackDir + os.sep 
                                      + "coreg" + os.sep + self.subswath 
                                      + os.sep + "stack*.dim")
                if not imgs:
                    raise Exception('No Sentinel-1 images prepared by SNAP'
                                    +' in specified stack directory')
                self.masterMetadata = s1Utils.readMetadata(imgs[-1],
                                                           MS='master')
                self.masterDate = self.masterMetadata['acqDate'].strftime("%Y%m%d")
                print('Reading stack ' +self.id+ ' , ' +
                      str(len(imgs))+ ' SLCs...')
                d = dict()
                for file in imgs:
                    metadata = s1Utils.readMetadata(file,MS='slave')
                    acqStr = metadata['acqDate'].strftime("%Y%m%d")
                    d[acqStr] = metadata.copy()
                self.metadata = OrderedDict(sorted(d.items(), 
                                                   key=lambda t: t[0]))
            self.acqDates = list(self.metadata.keys())
            self.files = sorted(imgs)
            self.orbit = metadata['orbit']
        else:
            print('Sensor not yet implemented... skipping.')

    def reduce(self, startDate):
        """
        startDate in 'YYYYMMDD' format
        """
        if hasattr(self, 'acqDates'):
            # remove dates before startDate:
            start_idx = next(x for x, val in enumerate(self.acqDates) 
                             if val > startDate) 
            if self.masterDate not in self.acqDates[start_idx:]:
                print('Start date specified after master date. Quiting.')
                raise
            else:
                self.acqDates = self.acqDates[start_idx:]
            # adapt files and metadata dicts:
            self.files = self.files[start_idx:]
            tmp = dict()
            for k,v in self.metadata.items():
               if k > startDate:
                   tmp[k] = v
            self.metadata = tmp
        else:
            print('Perform .readData() first!')

    def remove_idx(self, indices):
        if hasattr(self, 'acqDates'):
            if self.masterIdx in indices:
                print('Cannot delete master date!')
                raise
            else:
                removeDates = np.array(self.acqDates)[np.array(indices)]
                self.acqDates = [j for i, j in enumerate(self.acqDates) 
                                 if i not in indices]
                self.files = [j for i, j in enumerate(self.files) 
                              if i not in indices]
                for k in removeDates:
                    del self.metadata[k]

    @property
    def nSLC(self):
        if self._nSLC is None:
            if hasattr(self, 'files'):
                self._nSLC = len(self.files)
            else:
                print('Perform .readData() first!')
        return self._nSLC
    
    @property
    def masterIdx(self):
        #if self._masterIdx is None:
        if hasattr(self, 'files'):
            self._masterIdx = np.where(np.array(self.acqDates) 
                                       == self.masterDate)[0][0]
        else:
            print('Perform .readData() first!')
        return self._masterIdx


class RawStack(Stack):
    pass        


class CoregStack(Stack):
    pass

        
class TimeSeries:
    pass


#class TimeSeries_insar:
    


class crEncoder(json.JSONEncoder):
     def default(self, obj):
         if isinstance(obj, complex):
             return {"__complex__": True,
                     "real": obj.real,
                     "imag": obj.imag}
         if isinstance(obj, datetime.datetime):
             return {"__datetime__": True,
                     "str": obj.__str__()}
         if isinstance(obj, np.ndarray):
             if obj.dtype == stackType:
                 return {"__stackType__": True,
                         "array": obj.tolist()}
             else:
                 return {"__ndarray__": True,
                         "array": obj.tolist()}
         if isinstance(obj, np.int64):
             return int(obj)
         if isinstance(obj, Coordinates):
             return obj.__dict__
         if isinstance(obj, Stack):
             return obj.__dict__
         if isinstance(obj, Reflector):
             return obj.__dict__
         if isinstance(obj, CornerReflector):
             return obj.__dict__
         if isinstance(obj, Transponder):
             return obj.__dict__
         # Let the base class default method raise the TypeError
         return json.JSONEncoder.default(self, obj)


class crDecoder(json.JSONDecoder):
    def __init__(self, *args, **kwargs):
        json.JSONDecoder.__init__(self, object_hook=self.object_hook, 
                                  *args, **kwargs)
    def object_hook(self, dct):
        if '__complex__' in dct:
            return complex(dct["real"], dct["imag"])
        if "__datetime__" in dct:
            date = datetime.datetime.strptime(dct['str'],
                                              '%Y-%m-%d %H:%M:%S')
            return date
        if "stacks" in dct:
            if dct['type'] == 'CREF':
                cr = CornerReflector(dct['id'],
                                     dct['ascending'],
                                     dct['descending'],
                                     dct['installDate'],
                                     dct['startDate'],
                                     dct['endDate'],
                                     dct['geometry']['shape'],
                                     dct['geometry']['a0'],
                                     dct['geometry']['zenDip'],
                                     dct['geometry']['aziDip'])
            elif dct['type'] == 'CAT':
                cr = Transponder(dct['id'],
                                 dct['ascending'],
                                 dct['descending'],
                                 dct['installDate'],
                                 dct['startDate'],
                                 dct['endDate'],
                                 dct['RCS0'])    
            #cr = Reflector(dct['id'],dct['plh'],dct['installDate'],dct['startDate'])
            #if "type" in dct:
            #    cr.defGeometry(dct['type'],dct['a0'],dct['zenDip'],dct['aziDip'])
            cr.stacks = dct['stacks']
            return cr
        if 'crs' in dct:
            coords = Coordinates(dct['latitude'],
                                 dct['longitude'],
                                 dct['elevation'],
                                 dct['crs'],
                                 dct['frame'],
                                 dct['epoch'],
                                 dct['epsg'])
            return coords
        if 'stackDir' in dct:
            stack = Stack(dct['id'],
                          dct['stackDir'],
                          dct['sensor'],
                          dct['subswath'],
                          dct['type'])
            stack.metadata = dct['metadata']
            stack.acqDates = dct['acqDates']
            stack.files = dct['files']
            stack.orbit = dct['orbit']
            if dct['type'] == 'coreg':
                stack.masterMetadata = dct['masterMetadata']
                stack.masterDate = dct['masterDate']
            return stack
        if "__ndarray__" in dct:
            return np.array(dct["array"])
        if "__stackType__" in dct:
            idx = 0
            data = np.zeros(len(dct["array"]), dtype=stackType)
            for q in dct["array"]:
                #data[idx] = np.array(tuple(q),dtype=stackType)
                data[idx] = tuple(q)
                idx += 1
            return data
        return dct