# -*- coding: utf-8 -*-
"""getMaster

Routine to prepare master for coregistration routine. Either specified by user 
or automatically finds optimal master w.r.t. temporal and perpendicular baslines

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

def getMaster(masterDate,workDir,swaths):
    
    import os
    import glob
    from shutil import copy
    from shutil import copytree
    
###############################################################################
    
    slvDir = workDir + 'slaves/'
    
    # set optimal master if not explicitly specified:
    if not masterDate:
    
        from snappy import jpy
        from snappy import ProductIO
        
        slaves = [slvDir+slv for slv in os.listdir(slvDir) 
                    if slv.endswith(swaths[0]+'.dim')]
        
        #files_list = ['/home/rc/CR_Partizanske/DSC124/slaves/20200116_IW1.dim', '/home/rc/CR_Partizanske/DSC124/slaves/20200122_IW1.dim']
        product_set = []
        for file in slaves: 
            product_set.append(ProductIO.readProduct(file))
        
        InSARStackOverview = jpy.get_type('org.esa.s1tbx.insar.gpf.InSARStackOverview')
        
        optimalMaster = InSARStackOverview.findOptimalMasterProduct(product_set);    
        name = optimalMaster.getName()
        masterDate = name[0:8]    
        print('Optimal master selected: '+masterDate)
        
    # Create master product
    masterDir = workDir+os.sep+'master'
    os.makedirs(masterDir)
    # copy:
    try:
        # list master products for all swaths:
        masterProducts = glob.glob(slvDir+masterDate+"*.dim")
        for masterProduct in masterProducts:
            # copy both .dim file and .data dir
            copy(masterProduct,masterDir)
            productDir = masterProduct[0:-4]+'.data'
            splitPath = os.path.split(productDir)
            copytree(productDir,masterDir+os.sep+splitPath[1])
    except:
        raise Exception("Product corresponding to selected master date does not exist!")    
    #for swath in swaths:
        
        
    print('Directory for master date '+masterDate+' prepared.')

if __name__ == "__main__":
    import sys
    # load parameters file:
    parmFile = sys.argv[1]
    with open(parmFile,'r') as inp:
        parameters = eval(inp.read())
    masterDate = parameters['Stack_parameters']['masterDate']
    workDir = parameters['Stack_parameters']['workDir']
    swaths = parameters['Stack_parameters']['swaths']
    # call function:
    getMaster(masterDate,workDir,swaths)
