import numpy as np
from Particle3D import Particle3D
import sys
import math as m
import matplotlib.pyplot as pyplot


class ParticleList(object):

    """
    Class that describes lists contaning 3D particle objects.

    Properties:
    Physical Parameters File: Text file that contains the physical parameters of the astronomical bodies.
    Particle List: List that contains the Particle3D properties (Label, Position, Velocity, Mass) of each astronomical body.
    Acceleration List: List that contains the instantaneous accelerations of each particle used to update positions and velocities.
    Distances List: List that contains lists of distances of each body from the Sun (Need Moon separately).
    Trajectory File: Text file which contains the output necessary for VMD use.
    Observables File: Text file which contains the Astrophysical observables.
    Energies File: Text file which contains the kinetic, potential and total energy of the system for different timesteps.
    Labels List: List that contains the labels (names) of each body being simulated.
    y_axis List: List containing lists of y-positions of each body.

    Methods:
    Formatted Outputs.
    Centre of Mass Velocity.
    Position and Velocity Updates.
    Potential, Kinetic, Total Energy, Relative Energy Error.
    Observables (Apo/Peri-apses, Orbital Periods, Keplerian Periods, Semi-Major/Minor axes, Eccentricites).
    """


    
    def __init__(self, vmd_file, particle_file):
        
        """
        Initialises a ParticleList instance.
        Class Variables:
        :file_handle: Text file from which the physical parameters of the astronomical bodies are taken.
        :N: The total number of astronomical bodies simulated.
        :particlelist: The list containing the Particle3D properties of each body.
        :acceleration_list: The list containing the instantaneous accelerations.
        :distances: List containing the lists of distances of each body from the Sun.
        :outfile: Text file in which VMD outputted format will be written.
        :labels: List that contains the labels (names) of each body being simulated.
        :y_axis: List containing lists of y-positions of each body.
        :kin_energy: Kinetic Energy of the system.
        :pot_energy: Potential Energy of the system.
        :tot_energy: Total Energy of the system.
        :energy_error: Maximum deviation of energy from initial value.
        """


    
        #Opens the file name that contains the parameters of the simulation.
    
        file_parameters = open(parameters_file, "r")
    
        line = file_parameters.readline()
        args = line.split(" ")
    
        #Simulation Parameters and initial conditions.
        self.dt=float(args[0]) 
        self.numstep=int(args[1]) 
        self.ne=int(args[2])
        self.time= float(args[3])

        file_parameters.close()


        
         
        file_handle = open(particle_file, "r")
        N = len(file_handle.readlines())
        file_handle.close
        self.file_handle = open(particle_file, "r")
        self.N=N

        
        self.particlelist = []
        self.acceleration_list=[]
        self.labels = []
        self.time_list=[self.time]
        self.energy_error = 0.0 #alo
        self.list_creator()
        

            
        self.distances = []
        self.y_axis = []
        self.velocity_list = []
        
        
        for i in range(self.N-1):
            self.distances.append([])

        for i in range(self.N):
            self.y_axis.append([])
            self.velocity_list.append([])
            
            
        
        #Output files, trajectory, observables, energies
        self.outfile = open(vmd_file, "w")
        
        observables_file = "Observables.txt"
        self.observables_file = open(observables_file, "w")

        energy_file = "Energies.txt"
        self.energy_file = open(energy_file, "w")
        self.energy_file.write("Timestep  Potential Energy (J)      Kinetic Energy (J)       Total Energy (J)" + "\n")
        
    
        
    def __str__(self):
        
        """
        Defines output format.
        Returns particlelist as string, it already has the correct
        format for use by VMD as set in the Particle3D class.
        """
        
        return str(self.particlelist)


  
    def list_creator(self):
        
        """
        Method that fills the particlelist with the Particle3D objects created.
        Fills the labels list with the corresponding names.
        Creates the acceleration_list using the acceleration method.
        Corrects for centre of mass velocity using the com_velocity method.
        """
        
        for i in range(self.N):
            self.particlelist.append(Particle3D(self.file_handle))
            self.labels.append(Particle3D.label(self.particlelist[i]))
          
        self.acceleration_list=self.acceleration()
        self.com_velocity()
        
        

    def com_velocity(self):
        
        """
        Method that calculates the centre of mass velocity of the solar system due to the non-vanishing linear momentum
        of the initial conditions and subtracts it from the velocities of all the bodies to correct for it.
        """
        
        tot_mass=0
        com_momentum=np.array([0,0,0])
        for i in range(self.N):
            tot_mass+=self.particlelist[i].mass
            com_momentum = com_momentum + self.particlelist[i].velocity*self.particlelist[i].mass
        
        for i in range(self.N):
            self.particlelist[i].velocity=self.particlelist[i].velocity-(com_momentum/tot_mass)
        

           
    def next_step(self): #,dt):
        
        """
        Appends y-position of each body to the y_axis list.
        Second order position update using the Particle3D update_position method.
        First order velocity update using the Particle3D update_velocity method.
        The velocity update is done using Verlet. ***(Needs better wording)***
        Argument dt: Timestep set by Parameters text file. (Float)
        """
        for k in range(self.numstep):
            for i in range(self.N):
                self.particlelist[i].update_position(self.dt, self.acceleration_list[i])
                self.y_axis[i].append(self.particlelist[i].position[1]-self.particlelist[0].position[1])
            
            acceleration_list_new=self.acceleration()
        
            for i in range(self.N):
                self.particlelist[i].update_velocity(self.dt, (self.acceleration_list[i]+acceleration_list_new[i])*0.5)
                self.velocity_list[i].append(self.particlelist[i].velocity-self.particlelist[0].velocity)
            self.velocity_list[self.N-1][-1] = self.velocity_list[self.N-1][-1]- self.particlelist[3].velocity
            
            self.acceleration_list=acceleration_list_new
            

        
        
    def trajectory(self,i):
        
        """
        Method that writes the required output for VMD for every n'th (argument) timestep.
        n is taken from the parameters text file.
        Argument i is a dummy variable.
        """
        
        if i % self.ne == 0:
            self.outfile.write(str(self.N) + "\n")
            self.outfile.write("Point = " + str(i/self.ne+1) + "\n")
            for j in range(self.N):
                self.outfile.write(str(Particle3D.__str__(self.particlelist[j])) + "\n")

            
         
    def acceleration(self):
        
        """
        Method that computes instantaneous accelerations of all astronomical bodies
        using Newton's Gravitational Law and Newton's 3rd Law for reduced computation times and appends each to list.
        Returns acceleration_list which contains the instantaneous accelerations for each object.
        """
        
        acceleration_list=np.zeros((self.N,3))
       
        for i in range(self.N):
            j=i+1
            while j<=self.N-1:
                acceleration_list[i]+=Particle3D.force_grav(self.particlelist[i],self.particlelist[j])
                acceleration_list[j]-=acceleration_list[i]
                j+=1
            acceleration_list[i]=acceleration_list[i]/self.particlelist[i].mass
            
        return acceleration_list



    def potential_energy(self):

        """
        Method that computes the potential energy of the system due to gravitational interactions.
        Avoids double-counting of the pairwise interactions.
        Returns pot_energy, the potential energy of the system.
        """

        pot_energy = 0
        for i in range(self.N):
            j = i+1
            while j <= self.N - 1:
                pot_energy += Particle3D.pot_grav(self.particlelist[i], self.particlelist[j])
                j += 1
        
        return pot_energy



    def kinetic_energy(self):

        """
        Method that computes the kinetic energy of the system.
        Returns kin_energy, the kinetic energy of the system.
        """
        
        kin_energy = 0
        for i in range(self.N):
            kin_energy += self.particlelist[i].kinetic_energy()

        return kin_energy

    

    def energy(self):

        """
        Method that computes the total energy of the system.
        Total Energy = Kinetic + Potential.
        Returns Total Energy of the system.
        """
        
        #energy = self.kin_energy + self.pot_energy
        self.energy_list.append(self.kinetic_energy()+self.potential_energy())
        #self.tot_energy = energy
        
        #return energy





    
    def relative_energy_error(self):

        """
        Method that computes the maximum energy fluctuation of the system from the initial value as a percentage.
        """

        relative_energy_error = (self.energy_list[0] - np.min(self.energy_list))*100/self.energy_list[0]

        self.energy_error = relative_energy_error

        
        
    def apsis(self):
        
        """
        Appends the distances of each astronomical body from the Sun
        except the Moon (appends its distance from the Earth) to the corresponding list
        in the distances class variable, using the Particle3D v_sub_mag method.
        """
        
        for i in range(1,self.N-1):
            self.distances[i-1].append(Particle3D.v_sub_mag(self.particlelist[0], self.particlelist[i]))
        self.distances[self.N-2].append(Particle3D.v_sub_mag(self.particlelist[3], self.particlelist[self.N-1]))

        
    
    def moon(self):
        
        """
        Method used exclusively to determine the period of the Moon.
        Computes the relative y positions between the Earth and the Moon.
        """
        
        for i in range(self.numstep):
            self.y_axis[-1][i]=self.y_axis[-1][i]-self.y_axis[3][i]

        

    def observables(self):
        
        """
        Method that computes apo/peri-apsis. It then uses these positions to compute the orbital periods of each body.
        Apses are determined by finding a turning point in the magnitude of distance from the Sun for each body. They are then averaged.
        Eccentricities of the orbits are computed by the formula e = (apoapsis - periapsis)/(apoapsis + periapsis) where apo/peri-apsis correspond to the averaged values respectively.
        Semi-major and Semi-minor axes are computed by the formulae: semi-major = (apoapsis + periapsis)/2 and semi-minor = sqrt(apoapsis*periapsis).
        Orbital periods are determined using a y-axis intercept method and a first order correction. They are then averaged.
        Writes Apo/Peri-apses lists in units of kilometres to an output file.
        Writes Orbital Periods in units of days to the same output file.
        Writes Semi-major and Semi-minor axes in units of kilometeres to the same output file.
        Writes Orbital Eccentricities (dimensionless) to the same file.
        Format:
        Rows - Apoapsis, Periapsis, Period, Semi-Major, Semi-Minor, Eccentricity.
        Columns - Names of bodies. 
        Writes the final Potential, Kinetic and Total Energies and the Relative Energy Error to the same output file.
        Argument numstep is the total number of steps that is read from the parameters text file.
        Argument time_list is the list created to hold the times of each step.
        """
        
        apoapsis = []
        periapsis = []
        y_intercept_time = []
        period_y = np.zeros(self.N-1)
        eccentricity = []
        semi_major = []
        semi_minor = []
        self.moon()
        kepler_periods = []

        #Append empty lists corresponding to each body.
        for i in range(self.N - 1):
            apoapsis.append([])
            periapsis.append([])
            y_intercept_time.append([])
            
            
            
            
        #Compute apo/peri-apses and their times with both methods described. Appends them to corresponding lists. 
        for i in range(self.N-1):
            
            distances = self.distances[i]
            y_axis = self.y_axis[i+1]
            velocity = self.velocity_list[i+1]
            
            
            for j in range(1,numstep-1):
                if distances[j-1] < distances[j] and distances[j+1] < distances[j]:
                    apoapsis[i].append(distances[j])
                    
                    
                if distances[j-1] > distances[j] and distances[j+1] > distances[j]:
                    periapsis[i].append(distances[j])
                    
                    
                  
                if (y_axis[j-1]<0 and y_axis[j]>=0):
                    d_time=y_axis[j]/velocity[j][1]
                    y_intercept_time[i].append(self.time_list[j])

                
               

        
        apoapsis_average = []
        periapsis_average = []
        
        #Averages apoapses and periapses. Converts to km.
        for k in range(len(self.distances)):
            apoapsis_average.append(np.mean(apoapsis[k])/1000)
            
            periapsis_average.append(np.mean(periapsis[k])/1000)
            
        #Computes periods with y-intercept method. Also computes eccentricities, semi-major/minor axes.
        for i in range(self.N-1):
            
            eccentricity.append((apoapsis_average[i] - periapsis_average[i])/(apoapsis_average[i] + periapsis_average[i]))
            semi_major.append((apoapsis_average[i] + periapsis_average[i])/2)
            semi_minor.append(m.sqrt(apoapsis_average[i]*periapsis_average[i]))
            
            for j in range(1,len(y_intercept_time[i])):
                difference = y_intercept_time[i][j]-y_intercept_time[i][j-1]
                period_y[i] = period_y[i]+(difference)
            period_y[i] = period_y[i]/(len(y_intercept_time[i])-1)
                
        #Computes periods using Kepler's Third Law for comparison. The Moon is ommitted since it is a strictly 3-body problem.     
        for i in range(self.N-2):
            kepler_periods.append(m.sqrt(4*m.pi**2*semi_major[i]**3/(6.67408e-20*(self.particlelist[0].mass+self.particlelist[i+1].mass)))/(3600*24))
            
        
        #Converts periods to days.
        period_y = period_y/(24*3600)

        
        
        #Writing observables to the same output file for the y-intercept method of periods.
        for j in range(1,self.N):
            
            print(str(self.labels[j])+ " - "  + " Apoapsis (km): " + str(apoapsis_average[j-1]) + ", Periapsis (km): " + str(periapsis_average[j-1]) + ", Orbital Period (Days): " + str(period_y[j-1]) + ", Eccentricity: " + str(eccentricity[j-1]))
        self.observables_file.write("                  ")
        
        for j in range(1,self.N):
            
            if len(str(self.labels[j])) <= 6:
                self.observables_file.write(str(self.labels[j]) + "           ")
            else:
                self.observables_file.write(str(self.labels[j]) + "         ")
        self.observables_file.write("\n")
        self.observables_file.write("Apoapsis (km):   ")
        
        for j in range(1,self.N):
            self.observables_file.write('%.6e' % apoapsis_average[j-1] + "    ")
            
        self.observables_file.write("\n")
        self.observables_file.write("Periapsis (km):  ")
        
        for j in range(1,self.N):
            self.observables_file.write('%.6e' % periapsis_average[j-1] + "    ")
            
        self.observables_file.write("\n")
        self.observables_file.write("Period (Days):    ")

        
        for j in range(1,self.N):
            
            if len('%.7g' % period_y[j-1]) <= 6:
                self.observables_file.write('%.7g' % period_y[j-1] + "          ")
            else:
                self.observables_file.write('%.7g' % period_y[j-1] + "        ")

        self.observables_file.write("\n")
        self.observables_file.write("Kepler (Days):    ")
        
        for j in range(self.N-2):
            
            if len('%.7g' % kepler_periods[j]) <= 6:
                self.observables_file.write('%.7g' % kepler_periods[j] + "          ")
            else:
                self.observables_file.write('%.7g' % kepler_periods[j] + "        ")
            
        self.observables_file.write("\n")
        self.observables_file.write("Semi-Major (km): ")
        
        for j in range(1,self.N):
            self.observables_file.write('%.6e' % semi_major[j-1] + "    ")

        
            
        self.observables_file.write("\n")
        self.observables_file.write("Semi-Minor (km): ")
        
        for j in range(1, self.N):
            self.observables_file.write('%.6e' % semi_minor[j-1] + "    ")
            
        self.observables_file.write("\n")
        self.observables_file.write("Eccentricity:    ")
        
        for j in range(1,self.N):
            self.observables_file.write('%.3e' % eccentricity[j-1] + "       ")

        for j in range(3):
            self.observables_file.write("\n")
        
        #Writing the final Potential, Kinetic, Total Energies and the Relative Energy error to the output file.
        self.observables_file.write("Final Energies:" + "\n")
        self.observables_file.write("Total Energy = " + str(self.energy_list[-1]) + " J " + "\n")
        self.observables_file.write("Maximum Relative Energy Error = " + str(self.relative_energy_error()) + " % ")


    def energy_graph(self):

        """
        Plots Energy against Time.
        """

        pyplot.title('Total Energy vs Time')
        pyplot.xlabel('Time (s)')
        pyplot.ylabel('Energy (J)')
        pyplot.plot(self.time_list, self.energy_list)
        pyplot.show()
 
