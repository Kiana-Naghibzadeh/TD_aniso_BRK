This code is used to model anisotropic grain growth in an experimentally observed Ni microstructure. 

First, the grain boundary energy distribution for all the interfaces is calculated prior to the main simulation and is stored in an array called "eFull.mat." Throughout the simulation with nucleation of new interfaces, this file will be updated. To evaluate energy values for a given misorientation and normal vector of a boundary, we used the Bulatov-Reed-Kumar (BRK) energy function. The MATLAB function GB5DOF.m used here is provided in the BRK paper:
Vasily V Bulatov, Bryan W Reed, and Mukul Kumar. Grain boundary energy function for fcc metals. Acta Materialia, 65:161–175, 2014.

Next, an anisotropic threshold dynamics simulation for grain growth is performed. The MATLAB function gbm3dParallel.m and some of the functions used therein are modified versions of the isotropic threshold dynamics algorithm provided at http://www.math.lsa.umich.edu/~esedoglu/Research/grains/grains.html (Esedoglu, S.; Otto, F. Threshold dynamics for networks with arbitrary surface tensions. Communications on Pure and Applied Mathematics. 68:5 (2015), pp. 808-864).
