********
Examples
********

InSAR processing
----------------



Forcing CR with unstable RCS into PS network
--------------------------------------------

First we perform full InSAR processing of our network.

Then, let's force new PSC layer for a given coordinates:

.. code:: python 

   import numpy as np
   plh = np.array([[
       48.9819620432/180*np.pi,
       21.3368903188/180*np.pi,
       515.6192],[
       48.9806016952/180*np.pi,
       21.3397319418/180*np.pi,
       543.5802
        ]]) 

   HDFfile = '/data/GUDS/CR_VYCHOD/insar_Presov_202109/insar_DSC153.hdf5'

   insarUtils.extract_PSC_manual(HDFfile, plh, psc_name = 'psc_C')

Now we can check on the extracted PSc, e.g. check their ampltide dispersion:

.. code:: python

   data['psc_C/D_A'][:]
   
Now let's connect them to 1st order network via densification routine. We use different
estimation bounds as we know these point undergo significant motion:

.. code:: python

   insarUtils.getPreciseH2PH(data['psc_C'], 
                          parms['outDir'],
                          parms['stackId'])    

   bounds = {'dH' : 5, # [m]
             'vel': 30,  # [mm/year]
             'seasonal': 5
             }

   insarUtils.densify_network(data, 
                               k = 3, 
                               mode_thr = 5, 
                               std_thr = 50,
                               plot_flag = 1,
                               ps_label = 'psc2',
                               network_label = 'network2',
                               new_ps_label = 'psc_C',
                               model = 'polynomial',
                               bounds = bounds)


Optionally check for unwrapping errors:

.. code:: python

   #%% check unwrapping errrors:
   plt.plot(data['network3']['ps_displ'][0])
   plt.plot(data['network3']['ps_displ'][0]-np.pi)
   plt.plot(data['network3']['ps_displ'][0]+np.pi)
