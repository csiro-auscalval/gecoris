#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""plotUtils

Module containing utility functions for plotting

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

import os
import sys
if sys.stdout.isatty(): # test if interactive shell or not
    import matplotlib
    matplotlib.use('Agg') # change plotting backed accordingly
import matplotlib.pyplot as plt
from matplotlib import dates, gridspec
formatter = dates.DateFormatter('%m-%y') 
import numpy as np
import scipy.stats as ss
# gecoris packages:
from gecoris import radarUtils, geoUtils, s1Utils
# constants:
speedOfLight = 299792458.0
    

def symbolGenerator():
    allSymbols = ['^','o','d','s','X','*','p','v','+']
    for symbol in allSymbols:
        yield symbol

        
def colorGenerator():
    allColors = ['C0','C1','C2','C3','C4','C5','C6','C7']
    for color in allColors:
        yield color


def plotSLC(SLCderamp,SLCovs,acqDate,outDir):
    fig, axes = plt.subplots(1,2,figsize=(12,12))
    im1 = axes[0].imshow(np.abs(SLCderamp))
    im2 = axes[1].imshow(np.abs(SLCovs))
    fig.colorbar(im1,ax=axes[0],fraction=0.046, pad=0.04)
    fig.colorbar(im2,ax=axes[1],fraction=0.046, pad=0.04)
    axes[0].set_title(acqDate)
    plt.savefig(outDir+acqDate+'.png', bbox_inches='tight')
    plt.close()
 

def plotBeta0(SLCderamp,metadata,outDir,ovs=True):
    
    if ovs:
        SLCovs = radarUtils.oversample(SLCderamp,32)
        beta0 = np.power(np.abs(SLCovs),2)/(metadata['beta0']**2)
        
        #idx = np.unravel_index(np.argmax(beta0, axis=None), beta0.shape)
        #beta0peak = beta0[idx]
        #RCSdB = radarUtils.beta2RCSdB(beta0peak,metadata)
        
        acqDate = metadata['acqDate'].isoformat()
        
        fig, axes = plt.subplots(1,2,figsize=(12,12))
        im1 = axes[0].imshow(np.abs(SLCderamp))
        im2 = axes[1].imshow(beta0)
        fig.colorbar(im1,ax=axes[0],fraction=0.046, pad=0.04,label = 'raw [-]')
        fig.colorbar(im2,ax=axes[1],fraction=0.046, pad=0.04,label = 'Beta0 [-]')
        axes[0].set_title(acqDate)
        axes[1].set_title(acqDate)
        #fig.suptitle(acqDate) 
        #axes[1].annotate('RCS = '+f"{RCSdB:.1f} $dBm^2$",xy=(32,32),
        #                 fontsize = 20,color = 'w')
    else:
        beta0 = np.power(np.abs(SLCderamp),2)/(metadata['beta0']**2)
        fig, axes = plt.subplots(figsize=(6,6))
        im = axes.imshow(beta0)
        fig.colorbar(im,ax=axes,fraction=0.046, pad=0.04,label = 'Beta0 [-]')
        acqDate = metadata['acqDate'].isoformat()
        axes.set_title(acqDate)
    
    plt.savefig(outDir+acqDate+'.png')#, bbox_inches='tight')
    plt.close()
    

# def plotALE(cr,stackId,outDir):
#     stackIdx = cr.getStackIdx(stackId)
#     stack = cr.stacks[stackIdx]
    
#     fig, axes = plt.subplots(figsize=(8,8))
#     axes.axhline(color='k')
#     axes.axvline(color='k')
    
#     startIdx = np.argmax(stack['data']['status'])
    
#     legendCounter = True
#     for acqDate in stack['data']['acqDate'][startIdx:]:
#         idx = np.where(np.char.find(stack['data']['acqDate'],acqDate)
#                     == 0)[0]
#         if idx in stack['goodIdx']:
#             metadata = stack['stack'].metadata[acqDate]
            
#             # from ellipsoidal to ECEF:
#             xyz = geoUtils.plh2xyz(cr.plh)
#             (tAzimuth_ETRS,tRange_ETRS,satVec) = radarUtils.xyz2t(xyz,metadata)
            
#             # transform from ETRS-89 to ITRS (ITRF2014), t = acqEpoch
#             decimalDate = geoUtils.decimalYear(metadata['acqDate'])
#             xyz = geoUtils.etrf2itrf(xyz,decimalDate)
#             (tAzimuth_ITRS,tRange_ITRS,satVec) = radarUtils.xyz2t(xyz,metadata)
            
#             # get SET:
#             decimalHour = (metadata['acqDate'].hour
#                             + metadata['acqDate'].minute/60 
#                             + metadata['acqDate'].second/3600)
#             dxyz_set = geoUtils.getSET(xyz,acqDate,decimalHour)
#             xyz += dxyz_set
#             (tAzimuth_SET,tRange_SET,satVec) = radarUtils.xyz2t(xyz,metadata)
            
#             # get tropo delay:
#             slantDelay = geoUtils.tropoDelay(xyz,satVec)
#             tRange_TROPO = tRange_SET + slantDelay    
            
#             # get iono delay
#             ionoDelay = 0.1/speedOfLight # TEMP
#             tRange_IONO = tRange_TROPO + ionoDelay
            
#             # S1 residual bistatic correction:
#             bistaticAz = s1Utils.bistaticCorrection(tRange_IONO,metadata)
#             tAzimuth_BISTATIC = tAzimuth_SET - bistaticAz    
            
#             # S1 Doppler shift
#             doppler = s1Utils.dopplerRgCorrection(tRange_IONO,
#                                                 tAzimuth_BISTATIC,metadata)
#             tRange_Doppler = tRange_IONO + doppler
            
#             # get IRF peak:
#             Azimuth = stack['data']['azimuth'][idx]
#             Range = stack['data']['range'][idx]
#             tAzimuth, tRange = radarUtils.radar2time(Azimuth,Range,metadata)
            
#             # convert to metres:
#             rangeFactor = metadata['RSR']*metadata['rangeSpacing']
#             azimuthFactor = metadata['PRF']*metadata['azimuthSpacing']
#             dAz_ETRS = (tAzimuth_ETRS - tAzimuth)*azimuthFactor
#             dR_ETRS = (tRange_ETRS - tRange)*rangeFactor
#             dAz_ITRS = (tAzimuth_ITRS - tAzimuth)*azimuthFactor
#             dR_ITRS = (tRange_ITRS - tRange)*rangeFactor
#             dAz_SET = (tAzimuth_SET - tAzimuth)*azimuthFactor
#             dR_SET = (tRange_SET - tRange)*rangeFactor
#             dAz_TROPO = dAz_SET
#             dR_TROPO = (tRange_TROPO - tRange)*rangeFactor
#             dAz_IONO = dAz_SET
#             dR_IONO = (tRange_IONO - tRange)*rangeFactor
#             dAz_BIST = (tAzimuth_BISTATIC - tAzimuth)*azimuthFactor
#             dR_BIST = dR_IONO
#             dAz_DOPPLER = dAz_BIST
#             dR_DOPPLER = (tRange_Doppler - tRange)*rangeFactor
            
#             symb = symbolGenerator()
#             col = colorGenerator()
#             axes.plot(dR_ETRS,dAz_ETRS,next(symb)+next(col),
#                       label='ETRS-89' if legendCounter else "")
#             axes.plot(dR_ITRS,dAz_ITRS,next(symb)+next(col),
#                       label='ITRS' if legendCounter else "")
#             axes.plot(dR_SET,dAz_SET,next(symb)+next(col),
#                       label='+ SET' if legendCounter else "")
#             axes.plot(dR_TROPO,dAz_TROPO,next(symb)+next(col),
#                       label='+ tropo' if legendCounter else "")
#             axes.plot(dR_IONO,dAz_IONO,next(symb)+next(col),
#                       label='+ iono' if legendCounter else "")
#             axes.plot(dR_BIST,dAz_BIST,next(symb)+next(col),
#                       label='+ bistatic' if legendCounter else "")
#             axes.plot(dR_DOPPLER,dAz_DOPPLER,next(symb)+next(col),
#                       label='+ Doppler' if legendCounter else "")
#             axes.plot(0,0,next(symb)+'k',markeredgewidth=5,markersize=20,
#                       label='RP' if legendCounter else "")
#             legendCounter = False
        
#     axes.tick_params(axis='both', which='major', direction='in', 
#                      length = 10, labelsize=14)
#     axes.set_ylabel(r'$\Delta$ Azimuth [m]',fontsize=14)
#     axes.set_xlabel(r'$\Delta$ Range [m]',fontsize=14)
#     axes.legend(fontsize = 16)
    
#     plt.subplots_adjust(hspace=0.5)
#     plt.savefig(outDir+os.sep+cr.id+'_ALE_'+stack['id']+'.png', 
#                 bbox_inches='tight')
#     plt.close()


def plotRCS(cr,outDir):
       
    symbGenerator = symbolGenerator()
    colGenerator = colorGenerator()
    installDate = dates.datestr2num(cr.installDate)
    startDate = dates.datestr2num(cr.startDate)
    if cr.endDate != '99999999':
        endDate = dates.datestr2num(cr.endDate)
    else:
        endDate = np.nan
    # prepare figure:
    fig, ax = plt.subplots(figsize=(9,4),dpi=120)
    
    anotText = []
    anotColor = []
    for stack in cr.stacks:
        color = next(colGenerator)
        #plotDates = dates.datestr2num(stack['stack'].acqDates)
        plotDates = dates.datestr2num(stack['data']['acqDate'])
        ax.plot(plotDates,stack['data']['RCS'],next(symbGenerator),
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
    
    RCS0 = np.mean([stack['RCS0'] for stack in cr.stacks])
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
    # vertical lines for installation and end:
    ax.axvline(x = installDate,linestyle=':', color='k')
    ax.axvline(x = startDate,linestyle='--', color='k')
    ax.axvline(x = endDate,linestyle='--', color='k')
    ax.grid(alpha = 0.3)        
    ax.set_title(cr.id,fontsize = 12)
    ax.legend(loc = 'lower right')
    plt.subplots_adjust(hspace=0.5)
    plt.savefig(outDir+cr.id+'.png', bbox_inches='tight')
    plt.close()


def plotRCS_transponder(cr,outDir):

    #formatter = dates.DateFormatter('%m-%y')        
    symbGenerator = symbolGenerator()
    colGenerator = colorGenerator()
    installDate = dates.datestr2num(cr.installDate)
    startDate = dates.datestr2num(cr.startDate)
    if cr.endDate != '99999999':
        endDate = dates.datestr2num(cr.endDate)
    else:
        endDate = np.nan
    # prepare figure:
    fig, ax = plt.subplots(figsize=(9,4),dpi=120)
    
    for stack in cr.stacks:
        color = next(colGenerator)
        #plotDates = dates.datestr2num(stack['stack'].acqDates)
        plotDates = dates.datestr2num(stack['data']['acqDate'])
        ax.plot(plotDates,stack['data']['RCS'],next(symbGenerator),
                color=color,markersize = 2.5,label = stack['id'])
        
        # plot estimated reflector RCS:
        ax.hlines(stack['reflRCS'],startDate,np.max(plotDates),
                  color=color,linewidth=1,linestyle = '-.')      
    
    ax.set_ylabel('Apparent RCS [$dBm^2$]')
    start, end = ax.get_xlim()
    stepsize = (end+30-start)/10
    ax.xaxis.set_ticks(np.arange(start, end+30, stepsize))
    ax.xaxis.set_major_formatter(formatter)
    ax.xaxis.set_tick_params(rotation=30, labelsize=10)
    ax.axvline(x = installDate,linestyle=':', color='k')
    ax.axvline(x = startDate,linestyle='--', color='k')
    ax.axvline(x = endDate,linestyle='--', color='k')
    ax.grid(alpha = 0.3)        
    ax.set_title(cr.id,fontsize = 12)
    ax.legend(loc = 'lower right')
    plt.subplots_adjust(hspace=0.5)
    plt.savefig(outDir+cr.id+'.png', bbox_inches='tight')
    plt.close()


def plotNetworkALE(CRs_plot,outFig):
    
    fig, axes = plt.subplots(2,1,figsize=(8,6),dpi=120)
    idx = 1
    n_stacks = len(CRs_plot[0].stacks)
    for cr in CRs_plot:
        col = colorGenerator()
        for stack in cr.stacks: 
            color = next(col)
            dAz = stack['data']['dAz']
            dR = stack['data']['dR']
            #dAz,dR = getALE(cr,stack['id'])
            dAz = dAz[np.nonzero(dAz)]
            dR = dR[np.nonzero(dR)]
            # boxplot:
            box = axes[0].boxplot(dAz[~np.isnan(dAz)],positions=[idx],
                                      showfliers=False, 
                                      widths = 1.4/n_stacks)
            box2 = axes[1].boxplot(dR[~np.isnan(dAz)],positions=[idx],
                                      showfliers=False, 
                                      widths = 1.4/n_stacks)
            # modify colors
            for _, line_list in box.items():
                for line in line_list:
                    line.set_color(color)
                    line.set_linewidth(1.5)
            for _, line_list in box2.items():
                for line in line_list:
                    line.set_color(color)
                    line.set_linewidth(1.5)
            idx += 2/n_stacks
        idx+=1.5
        
    maxX = len(CRs_plot)*3.5
    axes[0].set_xlim(0,maxX)
    axes[1].set_xlim(0,maxX)
    axes[0].set_xticks(np.arange(1.5,maxX,3.5))
    axes[1].set_xticks(np.arange(1.5,maxX,3.5))
    axes[0].set_xticklabels([])
    axes[1].set_xticklabels([cr.id for cr in CRs_plot])
    axes[0].yaxis.set_tick_params(labelsize=12, direction='in')
    axes[1].yaxis.set_tick_params(labelsize=12, direction='in')
    axes[0].xaxis.set_tick_params(direction='in')
    axes[1].xaxis.set_tick_params(rotation=90, labelsize=12, direction='in')
    axes[0].axhline(y = 0,linestyle=':', color='k')
    axes[1].axhline(y = 0,linestyle=':', color='k')
    axes[0].set_ylabel(r'$\Delta Azimuth \ [m]$',fontsize=14)
    axes[1].set_ylabel(r'$\Delta Range \ [m]$',fontsize=14)
    #axes[0].legend([axes[0].get_lines()[0], axes[0].get_lines()[10]], 
    #               ['asc', 'dsc'], loc='upper right')
    #axes[1].legend([axes[0].get_lines()[0], axes[0].get_lines()[10]], 
    #               ['asc', 'dsc'], loc='upper right')
    axes[1].xaxis.set_tick_params(rotation=25)
    plt.subplots_adjust(hspace=0.5)
    plt.tight_layout()
    plt.savefig(outFig, bbox_inches='tight')
    plt.close()
    

def plotNetworkRCS(CRs,outFig):
    fig, axes = plt.subplots(figsize=(6,4),dpi=120)
    legendCounter1 = True
    legendCounter2 = True
    for cr in CRs:
        col = colorGenerator()
        symb = symbolGenerator()
        for stack in cr.stacks: 
            track = stack['id'][3:6]
            color = next(col)
            symbol = next(symb)
            #RCS = stack['data']['RCS'][stack['goodIdx']]
            #metadata = stack['metadata']
            #ampl = radarUtils.RCSdB2amp(RCS,metadata)
            #mleRice = ss.rice.fit(ampl,loc=0,scale=1)
            #interval = ss.rice.interval(0.95,mleRice[0],loc=mleRice[1],
            #                            scale=mleRice[2])
            #RCSinterval = radarUtils.amp2RCSdB(interval,metadata)
            #sigRCSmle = np.abs((RCSinterval[1]-RCSinterval[0])/2/2.5)
            sigRCSmle = 2.5*stack['sigRCS']
            #if sigRCSmle < 0.05:
            #    sigRCSmle +=0.3
            #print(sigRCSmle)        
            plt.plot(stack['dtheta'],stack['RCS0'],'xk',
                     markeredgewidth=1.8,label=r'$RCS_0$' 
                     if legendCounter1 else "")
            plt.plot(stack['dtheta'],stack['reflRCS'],'o',color=color,
                     marker=symbol,label=r'$\hat{RCS}_{'+track+'}$'
                     if legendCounter2 else "",
                     markersize=6)
            #plt.errorbar(stack['dtheta'],stack['reflRCS'],yerr=sigRCSmle,
            #             color=(0.7, 0.7, 0.7),label=r'$1 \ \sigma$'
            #             if legendCounter1 else "")
            plt.errorbar(stack['dtheta'],stack['reflRCS'],yerr=sigRCSmle,
                         color=color,label=r'$2 \ \sigma$'
                         if legendCounter1 else "",
                         elinewidth=0.5)
            legendCounter1 = False
        legendCounter2 = False
    
    axes.legend(fontsize = 14,loc = 'best')
    axes.set_ylabel(r'RCS [$dBm^2$]',fontsize=14)
    axes.set_xlabel(r'$\theta - \theta_0 [deg]$',fontsize=14)
    #axes.set_xticks(np.arange(-7,-1))
    axes.xaxis.set_tick_params(labelsize=12)
    axes.yaxis.set_tick_params(labelsize=12)
    #axes.axvspan(-7, -6, facecolor='gray', alpha=0.35)
    #axes.axvspan(-3, -2, facecolor='gray', alpha=0.35)
    #axes.set_xlim(-7,-2)
    fig.savefig(outFig)
    plt.close()
    
    
# def plotNetworkSCR(CRs_plot,outFig):
    
#     width = 1
#     idx = 1
#     fig, ax = plt.subplots(figsize=(3,6),dpi=120)    
#     colors = ('#abd9e9','#fdae61','#2c7bb6','#d7191c')
    
#     for cr in CRs_plot:
#         colIdx = 0
#         for stack in cr.stacks: 
#             RCS = stack['data']['RCS'][stack['goodIdx']]
#             metadata = stack['metadata']
#             ampl = radarUtils.RCSdB2amp(RCS,metadata)
#             mleRice = ss.rice.fit(ampl,loc=0,scale=1)
#             interval = ss.rice.interval(0.95,mleRice[0],loc=mleRice[1],
#                                         scale=mleRice[2])
#             RCSinterval = radarUtils.amp2RCSdB(interval,metadata)
#             sigRCSmle = np.abs((RCSinterval[1]-RCSinterval[0])/2)
#             #if sigRCSmle < 0.05:
#             #    sigRCSmle +=0.3
#             sigSCR = sigRCSmle*np.sqrt(2)
            
#             SCRafter = stack['reflRCS']-stack['clutRCSafter']
#             SCRbefore = stack['reflRCS']-stack['clutRCSbefore']
        
#             # plot SCR after with errorbar
#             ax.barh(idx, SCRafter, xerr=sigSCR, align='center', alpha=1, 
#                     ecolor='black', capsize=2,height = width,
#                     color=colors[colIdx])
#             # plot SCR before
#             plt.barh(idx, SCRbefore, height=0.4*width, 
#                      color=colors[colIdx+1], alpha=0.9, label='before')
#             idx += 1
#             colIdx +=2
#         idx +=0.5
    
#     maxY = len(CRs_plot)*2.5
#     #ax.set_xlim(4,ax.get_xlim()[-1])
#     #ax.set_xlim(4,34.454435803759715)
#     ax.set_yticklabels([])
#     ax.set_yticks(np.arange(1.5,maxY,2.5))
#     ax.set_yticklabels([cr.id for cr in CRs_plot])
#     ax.set_ylim(0.2,maxY+0.5)
#     ax.set_xlabel('SCR [dB]',fontsize=14)
#     ax.xaxis.set_tick_params(labelsize=12)
    
#     plt.tight_layout()
#     plt.savefig(outFig, bbox_inches='tight')
#     plt.close()
    
# def plotNetworkSCR_new(CRs_plot,outFig):
#     width = 1
#     idx = 1
#     fig, ax = plt.subplots(figsize=(3,7),dpi=120)    
#     colors = ('#abd9e9','#fdae61','#2c7bb6','#d7191c') 
#     for cr in CRs_plot:
#         colIdx = 0
#         for stack in cr.stacks: 
#             sigSCR = 2*stack['sigRCS']*np.sqrt(2)        
#             #SCRafter = stack['reflRCS']-stack['clutRCSafter']
#             #SCRbefore = stack['reflRCS']-stack['clutRCSbefore']
#             if stack['clutRCSafter'] < 15:
#                 SCRafter = stack['reflRCS']-stack['clutRCSafter']
#             else:
#                 SCRafter = stack['reflRCS']-stack['clutRCSbefore']
#             SCRbefore = stack['RCS0'] - stack['clutRCSbefore']  
#             # plot SCR after with errorbar
#             ax.barh(idx, SCRafter, xerr=sigSCR, align='center', alpha=1, 
#                     ecolor='black', capsize=2,height = width,
#                     color=colors[colIdx])
#             # plot SCR before
#             plt.barh(idx, SCRbefore, height=0.4*width, 
#                      color=colors[colIdx+1], alpha=0.9, label='before')
#             idx += 1
#             colIdx +=2
#         idx +=0.5 
#     maxY = len(CRs_plot)*2.5
#     ax.set_yticklabels([])
#     ax.set_yticks(np.arange(1.5,maxY,2.5))
#     ax.set_yticklabels([cr.id for cr in CRs_plot])
#     ax.set_ylim(0.2,maxY+0.5)
#     ax.set_xlabel('SCR [dB]',fontsize=14)
#     ax.set_xticks(np.arange(0,40,10))
#     ax.xaxis.set_tick_params(labelsize=12)
#     plt.tight_layout()
#     plt.savefig(outFig, bbox_inches='tight')
#     plt.close()   

def plotNetworkSCR_hz(CRs_plot,outFig):
    n_stacks = len(CRs_plot[0].stacks)
    width = 2/n_stacks
    idx = 1
    fig, ax = plt.subplots(figsize=(7,3),dpi=120)    
    #colors = ('#abd9e9','#fdae61','#2c7bb6','#d7191c') 
    colors = ('#d73027',
              '#f46d43',
              '#fdae61',
              '#fee090',
              '#e0f3f8',
              '#abd9e9',
              '#74add1',
              '#4575b4')
    for cr in CRs_plot:
        colIdx = 0
        for stack in cr.stacks: 
            sigSCR = 2*stack['sigRCS']*np.sqrt(2)        
            #SCRafter = stack['reflRCS']-stack['clutRCSafter']
            #SCRbefore = stack['reflRCS']-stack['clutRCSbefore']
            if stack['clutRCSafter'] < 15:
                SCRafter = stack['reflRCS']-stack['clutRCSafter']
            else:
                SCRafter = stack['reflRCS']-stack['clutRCSbefore']
            SCRbefore = stack['RCS0'] - stack['clutRCSbefore']  
            # plot SCR after with errorbar
            ax.bar(idx, SCRafter, yerr=sigSCR, align='center', alpha=1, 
                    ecolor='black', capsize=4,width = width,
                    color=colors[colIdx])
            # plot SCR before
            plt.bar(idx, SCRbefore, width=0.4*width, 
                     color=colors[colIdx+1], alpha=0.9, label='before')
            idx += 2/n_stacks
            colIdx +=2
        idx +=0.5 
    maxY = len(CRs_plot)*2.5
    ax.set_xticklabels([])
    ax.set_xticks(np.arange(1.5,maxY,2.5))
    ax.set_xticklabels([cr.id for cr in CRs_plot])
    ax.set_xlim(0.2,maxY+0.5)
    ax.set_ylabel('SCR [dB]',fontsize=14)
    ax.set_yticks(np.arange(0,40,10))
    ax.yaxis.set_tick_params(labelsize=12)
    ax.xaxis.set_tick_params(rotation=25)
    plt.tight_layout()
    plt.savefig(outFig, bbox_inches='tight')
    plt.close()  


def plotCDF(stations,stackId,outDir):
    for stn in np.arange(len(stations)):
        
        stackIdx = stations[stn].getStackIdx(stackId)
        stack = stations[stn].stacks[stackIdx]
        RCS = stack['data']['RCS']
        metadata = stack['metadata']
        ampl = radarUtils.RCSdB2amp(RCS[stack['goodIdx']],metadata)
        mleRice = ss.rice.fit(ampl,loc=0,scale=1) 
        
        # plot cdf:
        fig, ax = plt.subplots(1, 2,figsize=(12,6))
        vals = ss.rice.ppf([0.001, 0.5, 0.999], mleRice[0],loc=mleRice[1],scale=mleRice[2])
        x = np.linspace(vals[0], vals[2], 500)
        ax[0].plot(x, ss.rice.cdf(x, mleRice[0],loc=mleRice[1],scale=mleRice[2]),
               'r-', lw=9, alpha=0.6, label='Rice cdf')
        # plot data:
        ax[0].plot(np.sort(ampl),np.linspace(0.001,0.999,len(ampl)),'ko',lw=3,
                label='sorted data')
        ax[0].set_xlabel(r'$A$',fontsize=15)
        ax[0].set_ylabel('probability | sample/Nsamples',fontsize=15)
        
        ax[0].plot(x, (ss.rice.pdf(x,mleRice[0],loc=mleRice[1],scale=mleRice[2])-mleRice[1])*
                   mleRice[2],'b--', lw=3, alpha=0.6, 
                   label='Rice pdf')
        
        ax[0].legend(loc='best', frameon=False,fontsize=16)
        ax[0].set_ylim((-0.01,1.01))
        ax[0].tick_params(labelsize=14)
        
        #% plot fitted clutter Rayleigh
        clutIdx = np.where(stack['data']['acqDate'] < stations[stn].installDate)[0]
        ampl = radarUtils.RCSdB2amp(RCS[clutIdx],metadata)
        mleRayleigh = ss.rayleigh.fit(ampl)
        
        vals = ss.rayleigh.ppf([0.001, 0.5, 0.999], loc=mleRayleigh[0],scale=mleRayleigh[1])
        x = np.linspace(vals[0], vals[2], 500)
        ax[1].plot(x, ss.rayleigh.cdf(x,loc=mleRayleigh[0],scale=mleRayleigh[1]),
               'r-', lw=9, alpha=0.6, label='Rayleigh cdf')
        # plot data:
        ax[1].plot(np.sort(ampl),np.linspace(0.001,0.999,len(ampl)),'ko',lw=3,
                label='sorted data')
        ax[1].set_xlabel(r'$A$',fontsize=15)
        #ax[1].set_ylabel('probability | sample/Nsamples',fontsize=13)
        
        ax[1].plot(x, (ss.rayleigh.pdf(x,loc=mleRayleigh[0],
                                   scale=mleRayleigh[1])-mleRayleigh[0])*mleRayleigh[1],
                'b--', lw=3, alpha=0.6, label='Rayleigh pdf')
        
        #clutRCS = radarUtils.amp2RCSdB(mleRayleigh[1],metadata)
        ax[1].legend(loc='best', frameon=False,fontsize=16)
        ax[1].set_ylim((-0.01,1.01))
        ax[1].tick_params(labelsize=14)
        
        plt.subplots_adjust(hspace=0.5)
        plt.savefig(outDir+stations[stn].id+'.png',bbox_inches='tight')


def KS(amp1,amp2,signif):
    """KS 2-sample test
    """
    return ss.ks_2samp(amp1,amp2)[1] > signif


def plotCDFtemporal(stations,stackId,outDir):
    ms = 3
    
    for stn in np.arange(len(stations)):
        
        stackIdx = stations[stn].getStackIdx(stackId)
        
        fig, ax = plt.subplots(1, 2,figsize=(9,4))
        # reflector:
        stack = stations[stn].stacks[stackIdx]
        RCS = stack['data']['RCS']
        metadata = stack['metadata']
        ampl = radarUtils.RCSdB2amp(RCS[stack['goodIdx']],metadata)
        
        half = int(np.floor(len(ampl)/2))
        qrt = int(np.floor(len(ampl)/4))
        ampl1 = ampl[:half]
        ampl2 = ampl[half:]
        ampl3 = ampl[:qrt]
        ampl4 = ampl[qrt:half]
        ampl5 = ampl[half:half+qrt]
        ampl6 = ampl[half+qrt:]
        
        ax[0].plot(np.sort(ampl3),np.linspace(0.001,0.999,len(ampl3)),'.',
                   markersize=ms,label=f'T/4; p = {ss.ks_2samp(ampl3,ampl4)[1]:.2f}')
        ax[0].plot(np.sort(ampl4),np.linspace(0.001,0.999,len(ampl4)),'s',
                   markersize=ms,label=f'T/4; p = {ss.ks_2samp(ampl4,ampl5)[1]:.2f}')
        ax[0].plot(np.sort(ampl5),np.linspace(0.001,0.999,len(ampl5)),'X',
                   markersize=ms,label=f'T/4; p = {ss.ks_2samp(ampl5,ampl6)[1]:.2f}')
        ax[0].plot(np.sort(ampl6),np.linspace(0.001,0.999,len(ampl6)),'d',
                   markersize=ms,label=f'T/4; p = {ss.ks_2samp(ampl6,ampl3)[1]:.2f}')
        ax[0].plot(np.sort(ampl1),np.linspace(0.001,0.999,len(ampl1)),'^',
                   markersize=ms,label=f'T/2; p = {ss.ks_2samp(ampl1,ampl2)[1]:.2f}')
        ax[0].plot(np.sort(ampl2),np.linspace(0.001,0.999,len(ampl2)),'*',
                   markersize=ms,label='T/2')
        #ax[0].plot(np.sort(ampl),np.linspace(0.001,0.999,len(ampl)),'X',
        #           markersize=ms,label = 'T',color='k')    
        ax[0].legend(loc='best', frameon=False,fontsize=11)
        ax[0].set_title('CR')
        ax[0].set_xlabel(r'$A$',fontsize=13)
        ax[0].set_ylabel('sample/Nsamples',fontsize=13)
        # clutter:
        clutIdx = np.where(stack['data']['acqDate'] < stations[stn].installDate)[0]
        ampl = radarUtils.RCSdB2amp(RCS[clutIdx],metadata)
        ampl1 = ampl[:int(np.floor(len(ampl)/2))]
        ampl2 = ampl[int(np.floor(len(ampl)/2)):]
        
        ax[1].plot(np.sort(ampl1),np.linspace(0.001,0.999,len(ampl1)),'^',
                   markersize=ms,label = f'T/2; p = {ss.ks_2samp(ampl1,ampl2)[1]:.2f}')
        ax[1].plot(np.sort(ampl2),np.linspace(0.001,0.999,len(ampl2)),'*',
                   markersize=ms,label = 'T/2')
        #ax[1].plot(np.sort(ampl),np.linspace(0.001,0.999,len(ampl)),'X',
        #           markersize=ms,label = 'T',color='k')
        ax[1].legend(loc='best', frameon=False,fontsize=11)
        ax[1].set_title('clutter')
        ax[1].set_xlabel(r'$A$',fontsize=13)
        #ax[1].set_ylabel('sample/Nsamples',fontsize=13)
        plt.subplots_adjust(hspace=0.5)
        plt.savefig(outDir+stations[stn].id+'.png', bbox_inches='tight')


#%% InSAR plots
    
def ps_ts(dates, ts, model):
    plt.plot(dates, ts)
    plt.plot(dates, model)
    

def plot_psc(psc, outfig = ''):
    plh = psc['plh'][:]*180/np.pi
    D_A = psc['D_A'][:]
    # TODO: log also 1st/2nd order
    print('{:n} PSC.'.format(plh.shape[0]))
    fig, ax = plt.subplots(figsize=(10, 10),dpi=200) 
    q = ax.scatter(plh[:,1], plh[:,0], c=D_A, s=10)
    ax.set_xlabel('Longitude [deg]')
    ax.yaxis.set_tick_params(rotation = 90)
    ax.set_ylabel('Latitude [deg]')
    ax.set_title('D_A of PSC')
    plt.set_cmap('viridis')
    fig.colorbar(q, ax=ax)
    if outfig:
        plt.savefig(outfig, bbox_inches='tight')
        # TODO: add order to title


def plot_network_parms(network, psc, outfig, stdThr = 50):
    # prepare data:
    parms = network['ps_parms'][:]
    res = network['ps_res'][:]
    refIdx = network.attrs['refIdx']
    #plh = geoUtils.xyz2plh(psc['xyz'][:])*180/np.pi
    plh = psc['plh'][:]*180/np.pi
    # remove ref ps
    if not np.isnan(refIdx):
        refLon = plh[refIdx,1]
        refLat = plh[refIdx,0]
        plh = np.delete(plh, refIdx, axis=0)
    vel = parms[:,2]
    resH = parms[:,1]
    std = np.std(res, axis=1)/-network.attrs['m2ph']*1e3
    idx = std < stdThr
    plh = plh[idx,:]
    vel = vel[idx]
    resH = resH[idx]
    std = std[idx]
    # figure:
    fig = plt.figure(figsize=(20, 10)) 
    gs = gridspec.GridSpec(2, 3, height_ratios=[2, 1]) 
    ax0 = plt.subplot(gs[0])
    velPlot = ax0.scatter(plh[:,1], plh[:,0], c=vel*1e3, s = 5)
    fig.colorbar(velPlot,ax=ax0, shrink = 0.6)
    if not np.isnan(refIdx):
        ax0.plot(refLon, refLat, '^k')
    ax0.set_title('Vel [mm/yr]')
    ax1 = plt.subplot(gs[1])
    hPlot = ax1.scatter(plh[:,1], plh[:,0], c=resH, s = 5)
    fig.colorbar(hPlot,ax=ax1, shrink = 0.6)
    if not np.isnan(refIdx):
        ax1.plot(refLon, refLat, '^k')
    ax1.set_title('Res.H [m]')
    ax2 = plt.subplot(gs[2])
    stdPlot = ax2.scatter(plh[:,1], plh[:,0], c=std, s = 5)
    fig.colorbar(stdPlot,ax=ax2, shrink = 0.6)
    if not np.isnan(refIdx):
        ax2.plot(refLon, refLat, '^k')
    ax2.set_title('STD [mm]')
    ax3 = plt.subplot(gs[3])
    ax3.hist(vel*1e3,bins=50)
    ax4 = plt.subplot(gs[4])
    ax4.hist(resH,bins=50)
    ax5 = plt.subplot(gs[5])
    ax5.hist(std,bins=50)
    plt.savefig(outfig, bbox_inches='tight')
    
    
def plot_network_quality(network, outDir = ''):
    if 'VC_mean' in network:
        fig = plt.figure(figsize=(10, 7)) 
        gs = gridspec.GridSpec(2, 2, height_ratios=[1, 0.8]) 
        ax0 = plt.subplot(gs[0,:])
        ax0.bar(network.attrs['dates'], 
                np.sqrt(network['VC_mean'])*180/np.pi,
                width = 3)
        start, end = ax0.get_xlim()
        stepsize = (end+30-start)/10
        ax0.xaxis.set_ticks(np.arange(start, end+30, stepsize))
        ax0.xaxis.set_major_formatter(formatter)
        ax0.xaxis.set_tick_params(rotation=30, labelsize=10)
        ax0.set_title('Mean Variance Components', fontweight='bold')
        ax0.set_ylabel('VC [deg]')
        # axes functions:
        def forward(q):
            return q/180*np.pi / np.abs(network.attrs['m2ph'])*1e3
        def inverse(q):
            return q*np.abs(network.attrs['m2ph'])*180/np.pi /1e3
        secax = ax0.secondary_yaxis('right', functions=(forward, inverse))
        secax.set_ylabel('[mm]')
        
        ax1 = plt.subplot(gs[1,0])
        ax1.hist(network['stdRes'][:],bins = 50)
        ax1.set_title('STD of residuals', fontweight='bold')
        ax1.set_xlabel('STD [mm]')
        ax1.set_ylabel('count')
        
        ax2 = plt.subplot(gs[1,1])
        ax2.hist(network['var'][:],bins = 50)
        ax2.set_title('Variance Factor', fontweight='bold')
        ax2.set_xlabel('$\hat{\sigma}^2$ [-]')
        ax2.set_ylabel('count')
        
        plt.subplots_adjust(hspace=0.1)
        if outDir:
            if 'aps_phase' in network:
                out_fig = outDir + os.sep + 'network_quality_afterAPS.png'
            else:
                out_fig = outDir + os.sep + 'network_quality_beforeAPS.png'
            plt.subplots_adjust(hspace=0.4)
            plt.savefig(out_fig, bbox_inches='tight', dpi=120)
    else:
        print('Perform temporal ambiguity step first.')


def plot_insar_ts(plot_dates, plot_ts, sig, refl_id, stack_id, out_dir = ''):

    fig, ax = plt.subplots(figsize=(9,4),dpi=120)
    ax.errorbar(plot_dates, plot_ts, fmt='-^', yerr=2.5*sig,
                ecolor='gray', label=stack_id, color='C0')
    ax.set_ylabel('LOS displacement [mm]')
    start, end = ax.get_xlim()
    stepsize = (end+15-start)/10
    ax.xaxis.set_ticks(np.arange(start, end+15, stepsize))
    ax.xaxis.set_major_formatter(formatter)
    ax.xaxis.set_tick_params(rotation=30, labelsize=10)
    ax.grid(alpha = 0.35) 
    ax.set_title(refl_id, fontsize = 14, fontweight='bold')
    ax.legend(loc = 'best', fontsize=14)
    if out_dir:
        out_fig = out_dir + os.sep + refl_id + '_' + stack_id + '_insar.png'
        plt.savefig(out_fig, bbox_inches='tight')
        

def plot_network(psc, network, stack):
    
    # parse data:
    azDist = psc['azimuth']*stack.attrs['azimuthSpacing']
    rDist = psc['range']*stack.attrs['rangeSpacing']
    # network types:
    arcs = network['arcs']
    # plot
    f, ax = plt.subplots(figsize=(8,8),dpi=120)
    for (i, j) in arcs:
        ax.plot([rDist[i], rDist[j]],
                [azDist[i], azDist[j]],
                "-r", lw = 0.5)
    ax.set_title('Estimation network')
    ax.set_xlabel('range [m]')
    ax.set_ylabel('azimuth [m]')
    ax.invert_yaxis()
    

def plot_mean_beta0(d, outDir = ''):
    if 'meanBeta' not in d:
        azSize, rgSize = d['master_SLC'].shape
        # allocate:    
        beta0array = np.zeros((azSize, rgSize, len(d['SLC'].keys())))
        idx = 0
        # prepare beta0 time series:
        for k in d['SLC'].keys():
            beta0array[:,:,idx] = (np.power(np.abs(d['SLC'][k][:]),2)/
                                   (d['SLC'][k].attrs['beta0']**2))
            idx += 1
        # temporal average beta0 in dB scale (closer to normal distrib.):    
        meanBeta_dB = np.nanmean(10*np.log10(beta0array+1e-14), axis=2)
        d.create_dataset('meanBeta', data = meanBeta_dB)
    else:
        meanBeta_dB = d['meanBeta']
    
    # plot:
    f, ax = plt.subplots(figsize=(10,10),dpi=300)
    im = ax.imshow(meanBeta_dB, cmap='gray', aspect=5)
    im.set_clim(-20,10)
    c = f.colorbar(im, shrink = 0.3, orientation='horizontal', pad = 0.05)
    c.ax.set_xlabel('Radar brightness [dB]', fontsize=9)
    c.ax.xaxis.set_tick_params(labelsize=8)
    c.ax.xaxis.set_label_position('top')
    ax.xaxis.set_tick_params(labelsize=7)
    ax.xaxis.set_label_position('top')
    ax.xaxis.tick_top()
    ax.yaxis.set_tick_params(labelsize=7, rotation = 90)
    ax.set_ylabel('azimuth [pix]', fontsize=9)
    ax.set_xlabel('range [pix]', fontsize=9)
    if outDir:
        out_fig = outDir + os.sep + 'meanBeta0.png'
        plt.subplots_adjust(hspace=0.2)
        plt.savefig(out_fig, bbox_inches='tight')


def plot_ref_point(psc, network):
    refIdx = network.attrs['refIdx']
    plh = geoUtils.xyz2plh(psc['xyz'][:])*180/np.pi
    refLon = plh[refIdx,1]
    refLat = plh[refIdx,0]
    # plot    
    fig,ax = plt.subplots(figsize=(5, 5))
    ax.scatter(plh[:,1], plh[:,0], s = 5)
    ax.plot(refLon, refLat, '^r', ms  = 10, label='ref.point')
    ax.legend()
    ax.set_xlabel('Longitude [deg]')
    ax.set_ylabel('Latitude [deg]')