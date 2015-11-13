isQG
==================

For Matlab users
------------------

Matlab users can use invert.py as a blackbox to calculate isQG solution from a set of data:

 - ssd(Ny,Nx): sea surface density (kg/m^3)
 - ssh(Ny,Nx): sea surface height (meter)
 - z(Nz): the vertical coordinate incrementing from bottom to surface, for example z=[-3000 ... -10] (meter). DO NOT place the surface at z=0. N^2 at z=0 does not make physical sense anyway.
 - n2(Nz): the stratification at z points. n2 = -g/rho0 d\rho / dz (/s).
 - lat(Ny): latitude in degree
 - lon(Nx): longitude in degree
 - useanomaly: (optional, default=True), a boolean variable. If true, anomaly fields of ssd and ssh are used in the invertion. The anomaly is the deviation  of the original 2D field from a best-fitted quadratic surface.

The above variables should be prepared and saved in a .mat file. Use the same names all in a lower case. Assume data are saved in datain.mat, use:
  
  python invert.py datain.mat dataout.mat

in command line or 

  run invert.py datain.mat dataout.mat

in iPython

to calculate the inversion and find the output data in dataout.mat (or any file name you specify).

It is written in Python. Try Enthought (https://www.enthought.com).

For Python users
------------------

Read the code. It is straight forward.


