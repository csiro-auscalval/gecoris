***********
Definitions
***********

Terms
-----

stack
   SLC time series per single orbit track

PSc
   Persistent Scatterer (PS) candidates

RCS
   Radar Cross Section


.. _cr_geometry_def:

Corner reflector geometry
-------------------------

.. figure:: _static/CR_geometry_definitions.png
    :align: center
    :figwidth: 600px

    Corner reflector geometry as defined in :download:`reflectors.json <../templates/reflectors.json>` template file.

    
InSAR processing parameters
---------------------------

.. _parms_table:

.. table:: InSAR processing parameters
   
   =====  =====  =======
   A      B      A and B
   =====  =====  =======
   False  False  False
   True   False  False
   False  True   False
   True   True   True
   =====  =====  =======


.. _HDF_datastack:

HDF5 InSAR datastack structure
------------------------------

Full structure of the HDF datastack::

	data
	    - stack
		- {attrs}
		    - masterDate
		    - orbit
		    - wavelength
		    - centerLat
		    - centerLon
		    - centerH
		    - azimuthSpacing
		    - rangeSpacing
		- slcs
		- h2ph
	    - psc
	    	- {attrs}
		- azIdx
		- rgIdx
		- azimuth
		- range
		- xyz
		- plh
		- incAngle
		- satvec
		- SLC
		- IFG
		- D_A
		- (idx)
		- (h2ph)
		- (aps_phase)
	    - psc_B
	    	-//-
	    - psc2
	    	-//-
	    - psc3
	    	-//-
	    - network
		- {attrs}
		    - masterDate
		    - dates
		    - m2ph
		    - n_ps
		    - refIdx
		    - refAz
		    - refR
		- arcs
		- arcPhase
		- B1
		- a_arcs
		- arcPhase_uw
		- arc_parms
		- var
		- VC_mean
		- stdRes
		- (apsPhase)
		- (vario_model)
	    - network2
		-//-
		- a_ps
		- ps_phase_unw
		- ps_parms
		- ps_sig_parms
		- ps_res
		- ps_var
		- ps_displ
		- (RPN)
	    - network_B
		- ps_displ
		- ps_parms
		- ps_res
		- ps_sig_parms
		- ps_var
		
.. _csv_definition:

CSV exchange standard for InSAR
-------------------------------

.. figure:: _static/csv_specification.png
    :align: center
    :figwidth: 600px
