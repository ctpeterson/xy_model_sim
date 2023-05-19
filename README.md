# XY Model with Python

This is a simple code for performing Monte Carlo simulations of the XY model using Python. The Monte Carlo is performed with Wolff cluster updates. 

To get started, create an XY model object as follows.
 ```
 # Import XY model code
 import xy_model as xym
 
 # "J" coupling (J/k*T)
 J = <value of J>
 
 # Linear dimension of lattice (total number of sites = N*N)
 N = <value of N>
 
 # Create XY model simulation object
 sim_obj = xym.XYSimulation(J, N)
 ```
There are optional arguments that specify the starting lattice, random number seeds, and the locations on disk to save or load spin configurations (if the user chooses to load or save any at all). The options for the "start" are "hot" or "cold", but if "config" is not None, the code will pick up a lattice located in "load_loc" and start from that lattice. To do a single cluster update, simple invoke the "cluster_update" method as follows.
```
# Do a single cluster update
sim_obj.cluster_update()
```
In the "cluster_update" method, you can choose the site to start the cluster from (if "None", then the starting site is random), whether or not you want to save the updated configuration (default is "save = False"), and whether or not you want the method to print some information out once it is done. You can also make simple measurements at any time as follows.
```
# Measure energy
sim_obj.energ()

# Measure magnetization
sim_obj.magn()

# Print out values of energy and magnetization
print(sim_obj.energy, sim_obj.mag)
```
And that's about it! 
