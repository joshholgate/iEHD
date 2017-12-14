# iEHDv2
ideal electrohydrodynamics simulator - VERSION 2.0!!!

1. What is iEHD?

iEHD is a simple implementation of a finite difference, level-set method for solving the "ideal electrohydrodynamic" equations for a conducting liquid in contact with a vacuum, namely the Navier-Stokes equation in the liquid region and the electric Laplace equation in the vacuum region, in axisymmetric (r,z) coordinates. Electrical, surface tension, viscous and gravitational forces are all included. The code is set up to simulate simple test problems for droplets and nearly-flat surfaces. NEW for version 2 - rotational forces now included!



2. Compiling and running iEHDv2.cc

The version of the code given here is now configured for a simulation of a charged droplet at half of the theoretical Rayleigh stability limit which slowly spins up before breaking up. The parameters for each problem are specified at compile time. An uncharged droplet in an electric field can be simulated by replacing lines 19 and 20 with (e.g. at the Taylor limit) #define APPL_E_WEIGHT 0.458258 and #define SELF_E_WEIGHT 0.0. Similarly, spinning is controlled by INIT_L and DLDT - L is the dimensionless angular momentum as defined by Liao and Hill (PRL, 2017).

A flat liquid surface in a perpendicular electric field can be simulated by commenting out #define DROPLET (lines 19+20 have no effect on the surface-in-field simulation).

Before running the code, ensure an "output" folder has been created in the iEHD working directory. Use of optimization flags is strongly recommended for better performance! (I use g++'s miracle-working -O3 flag).

NEW for version 2 - loading initial data from previous simulations! Simply uncomment #define READ_FILE and  follow it with the input data location.



3. Output and data

"macroscopic.dat" provides a timeseries of the volume of the fluid and amplitude of the droplet/surface perturbation.

Fluid and electric variables at each gridpoint and at specific points in time are contained in the files "t#_gridvalues.dat". The gnuplot script "profiles-init.plt" recursively calls "profiles-loop.plt" in order to produce a series of .png plots of the shape of the liquid surface.

NEW for version 2 - 3d rendering of droplets using MayaVi! Couldn't get a single python script which both finds the surface contour of the level-set funstion and renders the drop in 3d using MayaVi... Have written a short bash script to call the two python scripts "find_surface_contour.py" and "3d_render.py"



4. To whom can I complain about iEHD?

Any comments regarding iEHD can be addressed to j.holgate14@imperial.ac.uk.
