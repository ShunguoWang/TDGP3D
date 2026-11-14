# TDGP3DMT
Trans-Dimensional Gaussian Processes 3D Bayesian Inversion 

Magnetotellurics (MT) images the electrical resistivity of the subsurface. Resistivity is sensitive to fluids, melt, and conductive minerals, so MT is widely used to study crustal structure, volcanoes, geothermal systems, and ore deposits. Most MT data are now interpreted with 3D inversion, which gives a single smooth model. But many different models can fit the same data, and a single model cannot show how well its features are constrained.

Bayesian inversion solves this by producing many models that all fit the data, which shows the uncertainty directly. The method has developed in a few key steps:

1) 1D sampling: Bayesian sampling was first applied to 1D MT with a fixed number of layers (e.g. Grandis et al., 1999; Guo et al., 2011).
2) Trans-dimensional methods: these let the data decide how complex the model should be (e.g. Green, 1995; Malinverno, 2002).
3) Better sampling and models: parallel tempering improved sampling (Sambridge, 2014), and Gaussian processes made flexible, smooth models possible (Ray & Myer, 2019).
4) 2D MT: together, these advances extended Bayesian MT inversion from 1D to 2D (Mandolesi et al., 2018; Blatter et al., 2021).

The step to 3D has been held back by computing cost, because each of the many thousands of forward calculations is expensive.

In this work we present a trans-dimensional Bayesian inversion of 3D MT data. Gaussian processes turn a variable number of nodes into a smooth 3D resistivity model, and every model is tested with a full 3D forward calculation. Parallel computing cuts each calculation to about a few seconds. In a synthetic test with a buried conductive sphere, the method recovers the background resistivity, locates the conductor and its depth, and fits the data to the noise level.

Here, the code is in CODES_TDGP_m0, and the demonstration example is in the folder of 3D_10ohm_V5_kmin2_m0_7day. 
