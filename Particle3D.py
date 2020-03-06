"""
Computer Modelling Astronomical N-body Simulation Project. 
Authors: Ioannis Hadjifrangiskou (s1639255) and Michalis Chatzittofi (s1512900) 
Date: 13/03/2019
"""

import numpy as np

class Particle3D(object):

    """
    Class to describe 3D particles.

    Properties:
    position(array) - position along the x,y,z axes
    velocity(array) - velocity along the x,y,z axes
    mass(float) - particle mass

    Methods:
    Formatted output
    Kinetic energy
    First-order velocity update
    Second order position updates
    """

    def __init__(self, file_handle):

        """
        Initialise a Particle3D instance

        :label: label as str
        :position: position as array
        :velocity: velocity as array
        :mass: mass as float
        """
        
        line = file_handle.readline()
        args = line.split(" ")
        label = args[0]
        pos = np.array([args[1],args[2], args[3]], float)
        vel = np.array([args[4],args[5], args[6]], float)
        mass = float(args[7])
        
        self.label = label   
        self.position = pos
        self.velocity = vel
        self.mass = mass

    
    def __str__(self):

        """
        Define output format. Required for VMD. Divides positions by one million for easier recognition by older versions of VMD.
        """

        return str(self.label) + " " +  str(self.position[0]/1.0e6)+ " " +  str(self.position[1]/1.0e6) + " " + str(self.position[2]/1.0e6)
        
    
    def kinetic_energy(self):

        """
        Return kinetic energy as
        1/2*mass*vel^2
        """
        
        return 0.5*self.mass*np.linalg.norm(self.velocity)**2


    def update_velocity(self, dt, acceleration):

        """
        First-order velocity update,
        v(t+dt) = v(t) + dt*a(t)

        :dt: timestep as float
        :acceleration: acceleration of particle as array
        """

        self.velocity = self.velocity + dt*acceleration


    def update_position(self, dt, acceleration):
        
        """
        Second-order position update,
        x(t+dt) = x(t) + dt*v(t) + 1/2*dt^2*a(t)

        :dt: timestep as float
        :acceleration: current acceleration as array
        """

        self.position = self.position + dt*self.velocity + 0.5*dt**2*acceleration


    @staticmethod
    def v_sub(particle1, particle2):
        
        """
        Relative position between two particles.
        Returns the relative position as array. (From particle2 to particle1)
        :particle1: particle1 object.
        :particle2: particle2 object.
        """
        
        dr = particle1.position - particle2.position

        return dr

   
    @staticmethod
    def force_grav(particle1, particle2):

        """
        Method that computes the pairwise gravitational force between two particles.
        G is Newton's Gravitational Constant.
        Returns the force vector of particle1 as an array.
        :particle1: particle1 object.
        :particle2: particle2 object.
        """
        
        G = 6.67408E-11
        r12 = np.linalg.norm(Particle3D.v_sub(particle1, particle2))
        force = -G*particle1.mass*particle2.mass/(r12**3)
        direction_vector = Particle3D.v_sub(particle1, particle2)
        
        return force*direction_vector


    @staticmethod
    def pot_grav(particle1, particle2):

        """
        Method that computes the pairwise gravitational potential between two particles.
        G is Newton's Gravitational Constant
        Returns the potential energy of the interaction as a float
        :particle1: particle1 object.
        :particle2: particle2 object.
        """
        
        G = 6.67408E-11
        r12 = np.linalg.norm(Particle3D.v_sub(particle1, particle2))
        pot_energy = -G*particle1.mass*particle2.mass/r12
        
        return pot_energy
    

    @staticmethod
    def v_sub_mag(particle1, particle2):

        """
        Method that computes and returns the magnitude
        of the distance between two particles
        :particle1: particle1 object.
        :particle2: particle2 object.
        """
        
        return np.linalg.norm(Particle3D.v_sub(particle1, particle2))
        
    

    


        
       
    

