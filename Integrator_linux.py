"""
CMod Astronomical N-body Simulation Project.
Performs a Velocity Verlet time integration of a system, comprising of N particles
(Astronomical bodies) interacting via gravitational forces. It utilises two classes;
Particle3D and ParticleList to achieve this.

Produces a plot of the system's energy as a function of time. Outputs the positions
of each body in the correct format to allow visualisation of the simulation using
VMD (Visual Molecular Dynamics). Prints out Astrophysical Observables, Apo/Peri-apses
and orbital periods of each body around the central body of the system. 
"""


import sys
import numpy as np
from Particle3D import Particle3D
from ParticleList import ParticleList

#Begin main code
def main():

    Read names of input and output files from command line.
    if len(sys.argv)!=4:
        print("Wrong number of arguments.")
        print("Usage: " + sys.argv[0] + " <VMD output file>" + "<system file>" + "<parameters file>")
    else:
        trajectory_file = sys.argv[1] + ".xyz"
        particle_file = sys.argv[2]
        parameters_file = sys.argv[3]
    
    system=ParticleList(trajectory_file,particle_file)
    
    #Opens the file name that contains the parameters of the simulation.
    
    #file_parameters = open(parameters_file, "r")
    
    #line = file_parameters.readline()
    #args = line.split(" ")
    
    #Simulation Parameters and initial conditions.
    #dt=float(args[0]) 
    #numstep=int(args[1]) 
    #n=int(args[2])
    #time= float(args[3])
    
    #file_parameters.close()
    #system.potential_energy()
    #system.tot_kinetic_energy()
    #energy_list = [system.energy()]
    #time_list = [time]
    
   
    
    #Begins the time integration loop by calling appropriate class methods.
    for i in range(numstep):

        system.potential_energy()
        system.tot_kinetic_energy()
        energy = system.energy()
        system.energy_write(i,n)
        system.trajectory(i,n)
        system.apsis()
        system.next_step(dt)
        time+=dt
        time_list.append(time)
        energy_list.append(energy)
        
    
    system.relative_energy_error(energy_list)
    
    #Calls on class method that computes and prints Astrophysical Observables.
    system.observables(numstep,time_list)

    #Converts seconds to hours in time_list.
    time_list = [i/3600 for i in time_list]
    
    #Plots energy against time.
    system.energy_graph(time_list, energy_list)
    
  

main()
