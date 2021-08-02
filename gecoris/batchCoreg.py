# -*- coding: utf-8 -*-
"""batchCoreg

Routine to perform geometry-based oregistration and subsequent ESD correction 
of Sentinel-1 slave images prepared by batchPrepare.py in a single-master 
configuration.

Copyright (C) 2021 by R.Czikhardt

Email: czikhardt.richard@gmail.com
Last edit: 25.8.2020

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

print('Geodetic Corner Reflector (In)SAR Toolbox (GECORIS) v.1.0')
print('Copyright (c) 2021 Richard Czikhardt, czikhardt.richard@gmail.com')
print('Dept. of Theoretical Geodesy, Slovak University of Technology')
print('-----------------------------------------------------------------')
print('License: GNU GPL v3+')
print('-----------------------------------------------------------------')


import sys
import glob
import os
import subprocess
import shutil
import numpy as np

def batchCoreg(parmFile):
    # PARAMETERS:
    print('Initializing coregistration...')
    # read input parameters:
    snapPath = parameters['SNAP_parameters']['snapPath']
    repoPath = parameters['SNAP_parameters']['repoPath']
    ram = parameters['SNAP_parameters']['ram']
    cores = parameters['SNAP_parameters']['cores']
    workDir = parameters['Stack_parameters']['workDir']
    swaths = parameters['Stack_parameters']['swaths']
    
    ram = str(ram)
    cores = str(cores*10)
    
    # determine if ESD to be done:
    # - use AOI (latitude) size:
    min_lat = float(parameters['Stack_parameters']['min_lat'])
    max_lat = float(parameters['Stack_parameters']['max_lat'])
    latDif = (max_lat-min_lat)/180*np.pi*6371e3 # approx [m]
    if latDif > 20e3:
        ESD = True
        #if len(swaths) > 1:
        #    graph = 'graphs/CoregGraph_multiswath.xml'
        #    print('Merging multiple sub-swaths.')
        #else:
        graph = os.sep + 'graphs/CoregGraph_ESD.xml'
        graph = os.sep + 'graphs/coregGraph_update.xml'
        print('Including ESD correction in proc.chain.')
    else:
        ESD = False
        graph = os.sep + 'graphs/CoregGraph_smallAOI.xml'
        print('Small AOI, skipping ESD correction in proc.chain.')
    #TODO: multi-swath graph
    
    # if len(swaths) > 1:
    #     outDir = workDir + os.sep + 'coreg' + os.sep
    #     # find out if update requested:
    #     if not os.path.exists(outDir):
    #         os.makedirs(outDir)
    #         update = False
    #         print('Creating single-master stack.')
    #     else:
    #         alreadyCoreg = [slv[6:14] for slv in os.listdir(outDir) 
    #                         if slv.endswith(swaths[0]+'.dim') 
    #                         and slv.startswith('stack')]
    #         update = True
    #         print('Stack already exists, updating...') 
        
    #     # list all slaves dates:
    #     slavesDir = workDir + os.sep +"slaves" + os.sep
    #     imgs = glob.glob(slavesDir + "*"+swaths[0]+".dim")
    #     if update:
    #         imgs = [img for img in imgs if not any(aC in img 
    #                                                for aC in alreadyCoreg)]
    #     imgs.sort()
    #     print('New images in stack:')
    #     print(imgs)
    #     # load master:
    #     master = glob.glob(workDir +os.sep + "master" + os.sep 
    #                        + "*"+swaths[0]+".dim")[0]
    #     # perform coregistration:
    #     for img in imgs:
    #         master1 = master
    #         slave1 = str(img)
    #         master2 = master1.replace(swaths[0],swaths[1])
    #         slave2 = slave1.replace(swaths[0],swaths[1])
    #         master3 = master1.replace(swaths[0],swaths[2])
    #         slave3 = slave1.replace(swaths[0],swaths[2])
    #         cmd = ('gpt -c '+ram+'m -q '+cores+' -x -e ' + '\"' + repoPath 
    #                + graph + '\"' 
    #                + ' -PsourceFile1=' + master1 + ',' + slave1 
    #                + ' -PsourceFile2=' + master2 + ',' + slave2 
    #                + ' -PsourceFile3=' + master3 + ',' + slave3 
    #                + ' -PoutputStack=' + outDir + os.sep + 'stack_' + img[-16:-8] + '.dim'
    #                + ' -PoutputItfg=' + outDir + os.sep + 'itfg_' + img[-16:-8] + '.dim')           
    #         # switch between Linux/Windows:
    #         if os.name == 'nt':
    #             os.system(cmd)
    #         else:
    #             subprocess.call([snapPath+cmd],shell=True)
    #         print('Coregistration of '+img+' finished.')
    #     # imgs.replace(swaths[0],swaths[1])
    #     pass
    # else:
    # loop through subswaths:
    for subswath in swaths:
        
        outDir = workDir + os.sep + 'coreg' + os.sep + subswath
        
        # find out if update requested
        if not os.path.exists(outDir):
            os.makedirs(outDir)
            update = False
            print('Creating single-master stack.')
        else:
            alreadyCoreg = [slv[6:14] for slv in os.listdir(outDir) 
                            if slv.endswith('.dim') and slv.startswith('stack')]
            update = True
            print('Stack already exists, updating...') 
        
        # list all slaves
        slavesDir = workDir + os.sep +"slaves" + os.sep
        imgs = glob.glob(slavesDir + "*"+subswath+".dim")
        
        # if update, reduce by already coregistered slaves
        if update:
            imgs = [img for img in imgs if not any(aC in img for aC in alreadyCoreg)]
        
        imgs.sort()
        print('New images in stack:')
        print(imgs)
        
        # load master
        master = glob.glob(workDir +os.sep + "master" + os.sep + "*"+subswath+".dim")
        
        # perform coregistration:
        for iw in imgs: 
            cmd = ('gpt -c '+ram+'m -q '+cores+' -x -e ' + '\"' + repoPath + graph + '\" ' 
                   + '-PsourceFiles=' + master[0] + ',' + str(iw) 
                        + ' -PoutputStack=' + outDir + os.sep + 'stack_' + iw[-16:-8] + '.dim'
                        + ' -PoutputItfg=' + outDir + os.sep + 'itfg_' + iw[-16:-8] + '.dim')          
            # switch between Linux/Windows:
            if os.name == 'nt':
                os.system(cmd)
            else:
                subprocess.call([snapPath+cmd],shell=True)
            print('Coregistration of '+iw+' finished.')

        # remove corrupted (wrongly coregistered) dates:
        for slv in os.listdir(outDir):
            if slv.endswith('.data'):
                imgDir = outDir+os.sep+slv
                # test if .img file was created:
                if not any(File.endswith(".img") for File in os.listdir(imgDir)):
                    badImg = imgDir[-13:-5]
                    print('Files for date '+ badImg +' corrupted, removing.')
                    shutil.rmtree(imgDir)
                    dimFile = outDir+os.sep+'stack_'+badImg+'.dim'
                    try:
                        os.remove(dimFile)
                    except:
                        dimFile = outDir+os.sep+'itfg_'+badImg+'.dim'
                        os.remove(dimFile)

    # merge sub-swaths:
    if len(swaths) > 1:
        print('Merging sub-swaths...')
        graph = 'graphs/CoregGraph_multiswath.xml'
        outDir = workDir + os.sep + 'coreg' + os.sep
        for stack in os.listdir(outDir+swaths[0]):
            if stack.endswith('.dim') and stack.startswith('stack'):
                iwDirs = [outDir+q+os.sep+stack for q in swaths]
                cmd = ('gpt -c '+ram+'m -q '+cores+' -x -e ' + '\"' 
                       + repoPath + graph + '\" ' 
                       + ' -PsourceFiles=' + ','.join(iwDirs) 
                       + ' -PoutputStack=' + outDir + stack
                       + ' -PoutputItfg=' + outDir + stack.replace('stack','itfg'))
                # switch between Linux/Windows:
                if os.name == 'nt':
                    os.system(cmd)
                else:
                    subprocess.call([snapPath+cmd],shell=True)
                print('Sub-swath merging of '+stack+' finished.')
        # finally remove individual subswaths files:

    print('Done. Thank you for using GECORIS. Do not forget to reference.')

if __name__ == "__main__":  
    # load parameters file:
    parmFile = sys.argv[1]
    with open(parmFile,'r') as inp:
        parameters = eval(inp.read())
    # call function:
    batchCoreg(parameters)
