import numpy as np
import matplotlib.pyplot as plt
from pylab import *

user_input = raw_input("Which halo? (no.npy)\n")
print ""
print "Functions:"
print "plot_velocities()"
print "plot_energy_ratio()"
#halos = np.load("halos/mass1473e+12quadrant16shift.npy")
halos = np.load(user_input+'.npy')
h= 0.7
rad = halos[1,:]
mass = halos[2,:]
x = halos[3,:] / h
y = halos[4,:] / h
z = halos[5,:] / h
vx = halos[6,:]
vy = halos[7,:]
vz = halos[8,:]
G = 4.302 * 10**-9  #Mpc M_s^-1 (km/s)^2
H_0 = 70
Potential = []
Kinetic = []
Energy = []
E = 0
system = 1

#KINEMATIC FUNCTIONS
#returns the relative velocity, taking into account hubble expansion
def relative_velocity(xnew, ynew, znew, system):
   vexpx = vx[system] + H_0 * (x[system]-xnew)
   vexpy = vy[system] + H_0 * (y[system]-ynew)
   vexpz = vz[system] + H_0 * (z[system]-znew)
   vrel = ((vexpx - vx[0])**2 + (vexpy - vy[0])**2 + (vexpz - vz[0])**2)**.5
    
   return vrel

#defines center of mass of system and returns components for new origin
def center_of_mass(system):
   totmass = np.sum(mass[0:system])
   compr = rad[0:system] * mass[0:system]
   compx = x[0:system] * mass[0:system]
   compy = y[0:system] * mass[0:system]
   compz = z[0:system] * mass[0:system]
   radnew = sum(compr) / totmass
   xnew = sum(compx) / totmass
   ynew = sum(compy) / totmass
   znew = sum(compz) / totmass
    
   return radnew, xnew, ynew, znew

#Produces energy output of every cluster step until the system is unbound
while (E <= 0):
   system += 1
   rnew, xnew, ynew, znew = center_of_mass(system)
   KE = 0
   PE = 0
    
   for i in range(0, system):
      targmass = mass[i]
      vrel = relative_velocity(xnew, ynew, znew, i)
      KE += .5 * targmass * vrel**2
      #print "targmass:", targmass
      #print "vrel:", vrel

      PE_comp = G * mass[i] * mass[i+1:system] / np.absolute((rad[i] - rad[i+1:system])/h)
      PE -= np.sum(PE_comp)

   E = KE + PE
   print "Potential: ", PE
   print "Kinetic: ", KE
   print "Total Energy: ", E
   print "System: ", system
   print ""

   Kinetic.append(KE)
   Potential.append(PE)
   Energy.append(E)
    
   if (system == halos.shape[1]):
      print "All halos bound"

print "system: ", system

def scale_mass(system):
   maximum = np.log10(np.amax(mass))
   minimum = np.log10(np.amin(mass))
   ranges = maximum - minimum
   boundmass = 50**((np.log10(mass[:system]) - minimum) / ranges + 1)
   unboundmass = 50**((np.log10(mass[system:]) - minimum) / ranges + 1)

   print "These are the bound masses:", boundmass
   return boundmass, unboundmass


def radial_velocity(system):
   vxpol = vx + H_0 * (x - x[0])
   vypol = vy + H_0 * (y - y[0])
   vzpol = vz + H_0 * (z - z[0])
   numerator = vxpol*(x-x[0]) + vypol*(y-y[0]) + vzpol*(z-z[0])
   denominator = ((x-x[0])**2 + (y-y[0])**2 + (z-z[0])**2)**.5

   vradius = np.zeros(numerator.shape[0])
   vradius[1:] = numerator[1:]/denominator[1:]
#   try:
#      np.seterr(all = 'raise')
#      vradius[1:] = numerator[1:]/denominator[1:]
#
#   except FloatingPointError:
#      vradius[0] = 0
#
   print "This is vrad:", vradius
   return vradius

########################BEGINNING OF ANALYSIS##############################

#TO DO: Figure out how to annotate the figures with their halo numbers!
def plot_velocities():
   listno = halos[0,:]
   print "What does system equal?:", system
   vrad = radial_velocity(system)

   fig = plt.figure()
   ax1 = fig.add_subplot(111)
   xax = np.linspace(0, rad[system+200], 100)
   yax = 70 * xax
   ax1.plot(xax, yax, color = 'c', linewidth = 2)

   bmass, ubmass = scale_mass(system)

   target = ax1.scatter(0, vrad[0], s = bmass[0], c = 'black', alpha = 0.7)
   bound = ax1.scatter(rad[1:system], vrad[1:system], s = bmass[1:], c = 'green', alpha = 0.7)
   unbound = ax1.scatter(rad[system:system+200], vrad[system:system+200], s = ubmass[0:200], c = 'red', alpha = 0.7)
   for i in range(0,system+200):
      ax1.annotate('%d'%listno[i], xy = (rad[i], vrad[i]), xytext = (rad[i], vrad[i]), fontweight = 'bold')

   plt.xlim(xmin = 0)
   ax = plt.axes()
   ax.xaxis.grid()
   ax.yaxis.grid()
   ax1.set_xlabel('radius [Mpc]')
   ax1.set_ylabel('velocity [km/s]')
   plt.legend((target, bound, unbound), ('Target Halo', 'Bound Halos', 'Unbound Halos'), scatterpoints = 1, loc = 'upper left', labelspacing = 1)
   massy = np.array_str(mass[0])
   plt.title('Target Halo Mass = ' + massy + " M_sun")
   #plt.savefig(user_input+'velocity.png')
   plt.show()
   plt.close()
   print "Plotting done!"


def plot_energy_ratio():
   T = np.array(Kinetic)
   V = np.array(Potential)
   ratio = T / V
   
   plt.clf()
   ratioo, = plt.plot(rad[0:ratio.shape[0]], np.abs(ratio), 'darkmagenta')
   plt.xlabel('r [Mpc]')
   plt.ylabel('T/V')
   plt.legend([ratioo],["T/V"], loc='best')
   plt.grid()
   massy = np.array_str(mass[0])
   plt.title('Target Halo Mass = ' + massy + " M_sun")
   #plt.savefig(user_input+'ratio.png')
   plt.show()
