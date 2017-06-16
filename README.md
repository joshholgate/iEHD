# iEHD
ideal electrohydrodynamics simulator

1. What is iEHD?

iEHD is a simple implementation of a finite difference, level-set method for solving the "ideal electrohydrodynamic" equations for a conducting liquid in contact with a vacuum, namely the Navier-Stokes equation in the liquid region and the electric Laplace equation in the vacuum region, in axisymmetric (r,z) coordinates. Electrical, surface tension, viscous and gravitational forces are all included. The code is set up to simulate simple test problems for droplets and nearly-flat surfaces.



2. Compiling and running iEHD.cc

The version of the code given here is configured for a simulation of a charged droplet at the theoretical Rayleigh stability limit. The parameters for each problem are specified at compile time. An uncharged droplet in an electric field can be simulated by replacing lines 18 and 19 with (e.g. at the Taylor limit) #define APPL_E_WEIGHT 0.458258 and #define SELF_E_WEIGHT 0.0

If simulating a droplet with charge and an applied field then lines 742-759 should be uncommented; this bit of code ensures that the droplet remains at the centre of the simulation domain rather than being accelerated in the direction of the electric field.

A flat liquid surface in a perpendicular electric field can be simulated by commenting out line 26 (lines 18+19 have no effect on the surface-in-field simulation).

Before running the code, ensure an "output" folder has been created in the iEHD working directory. Use of optimization flags is recommended for better performance (I use g++'s -O3 flag).



3. Output and data

"run_params.txt" is a record of the simulation parameters.

"fluid_params.dat" provides a timeseries of the volume of the fluid and amplitude of the droplet/surface perturbation.

Fluid and electric variables at each gridpoint and at specific points in time are contained in the files "t#_gridvalues.dat". The gnuplot script "profiles-init.plt" recursively calls "profiles-loop.plt" in order to produce a series of .png plots of the shape of the liquid surface.



4. To whom can I complain about iEHD?

Any comments regarding iEHD can be addressed to j.holgate14@imperial.ac.uk.
