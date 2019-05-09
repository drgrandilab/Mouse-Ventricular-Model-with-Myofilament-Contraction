Matlab code of the Morotti et al. mouse ventricular model (J Physiol. 2014) with integration of the Negroni et
al. contractile model (J Mol Cell Cardiol. 2015).

This model was used for the simulations shown in Surdo et al., Nat Commun. 2017.

________________________________________________________________________________________________________________
Contents:

readme.txt				 	this file

	.m files - model components

morotti_et_al_masterCompute.m		 	loads initial conditions (see below) and runs the simulation
morotti_et_al_masterODEfile.m		 	integrates the following model components:
morotti_et_al_barODEfile.m		 	- beta-adrenergic (PKA) phosphorylation module
morotti_et_al_camkiiODEfile.m		 	- CaMKII phosphorylation module
morotti_et_al_camODEfile.m		 	- CaM module
morotti_et_al_eccODEfile.m		 	- excitation-contraction coupling module (with myofilament)

	.mat files - initial conditions (isotonic contraction at 1-Hz pacing)

yfin_mouse_myofil_isoT_control_1Hz	 	control model (steady-state, Ligtot = 0 uM)
yfin_mouse_myofil_isoT_IBMX_1Hz_120s	 	120 s with Ligtot = 0.1 uM (uMyo, uXBCa & uXBcy = 1)
yfin_mouse_myofil_isoT_ISOall_1Hz_120s	 	120 s with Ligtot = 0.1 uM (uMyo, uXBCa & uXBcy = 0.5)
yfin_mouse_myofil_isoT_ISOxbca_1Hz_120s	 	120 s with Ligtot = 0.1 uM (uXBCa = 0.5; uMyo & uXBcy = 1) 
yfin_mouse_myofil_isoT_ISOxbcy_1Hz_120s	 	120 s Ligtot = 0.1 uM (uXBcy = 0.5; uMyo & uXBCa = 1)  
yfin_mouse_myofil_isoT_ISOtitin_1Hz_120s 	120 s with Ligtot = 0.1 uM (uMyo = 0.5; XBCa & uXBcy = 1)
________________________________________________________________________________________________________________

Note that Ligtot (default value 0 uM) is defined in 'morotti_et_al_masterCompute.m', and uMyo, uXBCa & uXBcy
(default value 1) are defined in 'morotti_et_al_eccODEfile.m'.
________________________________________________________________________________________________________________


References:

J.A. Negroni, S. Morotti, E.C. Lascano, A.V. Gomes, E. Grandi, J.L. Puglisi, D.M. Bers. ß-adrenergic effects on 
cardiac myofilaments and contraction in an integrated rabbit ventricular myocyte model. J Mol Cell Cardiol.
2015 Apr; 81:162-75.

S. Morotti, A.G. Edwards, A.D. McCulloch, D.M. Bers, E. Grandi. A novel computational model of mouse myocyte 
electrophysiology to assess the synergy between Na+ loading and CaMKII. J Physiol. 2014 Mar 15;592(Pt 6):1181-97.

N.C. Surdo, M. Berrera, A. Koschinski, M. Brescia, M.R. Machado, C. Carr, S. Morotti, E. Grandi, P. Wright,
D.M. Bers, J. Gorelik, S. Pantano, M. Zaccolo. FRET biosensor uncovers cAMP nano-domains at ß-adrenergic targets
that dictate precise tuning of cardiac contractility. Nat Commun. 2017.


Please cite the above papers when using this model.