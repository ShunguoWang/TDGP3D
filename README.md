Introduction

Magnetotellurics (MT) maps the electrical resistivity of the subsurface. Resistivity is sensitive to fluids, melt, and conductive minerals, so MT is widely used to study crustal structure, volcanoes, geothermal systems, and ore deposits. MT data are now usually interpreted with 3D inversion, which produces a single smooth model. However, many different models can fit the same data equally well, and a single model cannot show how well its features are constrained.

Bayesian inversion addresses this problem by producing a large set of models that all fit the data, so the uncertainty can be seen directly. Bayesian MT inversion has developed in three main steps:

1D sampling: Bayesian sampling was first applied to 1D MT models with a fixed number of layers (e.g., Grandis et al., 1999; Guo et al., 2011).
Trans-dimensional methods: these methods let the data decide how complex the model should be (e.g., Green, 1995; Malinverno, 2002), and they were later applied to 1D MT (Mandolesi et al., 2018).
2D MT: building on these advances, Bayesian MT inversion was extended from 1D to 2D (Blatter et al., 2021; Ray, 2021).

Extending this approach to 3D has been limited by computing cost, because sampling requires many thousands of forward calculations, each of which is expensive.

In this work, we present a trans-dimensional Bayesian inversion of 3D MT data. Gaussian processes convert a variable number of nodes into a smooth 3D resistivity model, and each model is evaluated with a full 3D forward calculation. Parallel computing reduces each calculation to a few seconds. In a synthetic test with a buried conductive sphere, the method recovers the background resistivity, locates the conductor and its depth, and fits the data to the noise level.

Code and Data Availability

The inversion code is provided in the folder CODES_TDGP3D_m0, and the synthetic example is provided in the folder 3D_10ohm_V5_kmin2_m0_7day. The 3D forward solver is available on request from the authors of Dai et al. (2022).

References

Blatter, D., Ray, A., & Key, K. (2021). Two-dimensional Bayesian inversion of magnetotelluric data using trans-dimensional Gaussian processes. Geophysical Journal International, 226(1), 548–563.

Dai, S., et al. (2022). [Full reference to be completed.]

Grandis, H., Menvielle, M., & Roussignol, M. (1999). Bayesian inversion with Markov chains—I. The magnetotelluric one-dimensional case. Geophysical Journal International, 138(3), 757–768.

Green, P. J. (1995). Reversible jump Markov chain Monte Carlo computation and Bayesian model determination. Biometrika, 82(4), 711–732.

Guo, R., Dosso, S. E., Liu, J., Dettmer, J., & Tong, X. (2011). Non-linearity in Bayesian 1-D magnetotelluric inversion. Geophysical Journal International, 185(2), 663–675.

Malinverno, A. (2002). Parsimonious Bayesian Markov chain Monte Carlo inversion in a nonlinear geophysical problem. Geophysical Journal International, 151(3), 675–688.

Mandolesi, E., Ogaya, X., Campanyà, J., & Piana Agostinetti, N. (2018). A reversible-jump Markov chain Monte Carlo algorithm for 1D inversion of magnetotelluric data. Computers & Geosciences, 113, 94–105.

Ray, A. (2021). Bayesian inversion using nested trans-dimensional Gaussian processes. Geophysical Journal International, 226(1), 302–326.
