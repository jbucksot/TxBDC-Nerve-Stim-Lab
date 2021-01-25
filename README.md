# TxBDC-Nerve-Stim-Lab
Source code which runs the computational model hosted at http://ec2-3-23-162-1.us-east-2.compute.amazonaws.com/
Code was written for use on Ubuntu 18.04.4

Last update: 1/25/2021

Requirements:
  FEniCS Project: https://fenicsproject.org/
  Gmsh 4.0: https://gmsh.info/
  meshio: https://pypi.org/project/meshio/
  pygmsh: https://pypi.org/project/pygmsh/
  
To use:
  1) Edit parameter values in each file to create the desired system.
  2) Run either create_mesh.py or create_IF_mesh.py.
  3) Run solve_fenics.py from the command line. Ensure that the parameter values match the values used in step 2.
  4) Run solve_fibers.py from the command line. Ensure that the parameter values match the values used in steps 2 and 3.
