import numpy as np

#normalized hubble parameter
h = .7

def find_halo(lbound, ubound, halos, omithalo):
   isSub = halos[9,:]
   mass = halos[1,:]
   xc = halos[3,:]
   yc = halos[4,:]
   zc = halos[5,:]

   #searches for a system of two halos within 1Mpc of each other 
   #who have a mass difference of <10**11
   for target in range(lbound, ubound):
      #print target == any(omithalo)
      if target == any(omithalo):
         print "OMITTED!!"
      else:
         for partner in range(target+1, ubound+1):
            #print "target:", target
            pradius = ((xc[target]-xc[partner])**2 + (yc[target]-yc[partner])**2 + (zc[target]-zc[partner])**2)**.5
            
            mdiff = np.abs(mass[target]-mass[partner])
            
            if pradius < 2.0:
               if target != partner:
                  if isSub[partner] == -1 and isSub[target] == -1:
                     if mdiff < 10**11:
                        print "targetmass: ", mass[target]
                        print "partnermass: ", mass[partner]
                        return target

   #returns -1 if runs off the list
   return -1

def save_halo(radius, halos, halo, quadrant):
   #calculates the relative radius of all halos to the target
   mass = halos[1,:]
   x = halos[3,:]
   y = halos[4,:]
   z = halos[5,:]
   vx = halos[6,:]
   vy = halos[7,:]
   vz = halos[8,:]
   print vz
   
   rad = ((x-x[halo])**2 + (y-y[halo])**2 + (z-z[halo])**2)**.5

   #finds where the relative radius is less than inputted Mpc
   boundhalos = np.where(rad < radius)
   print "boundhalos:", boundhalos

   #constructs desired components
   hlist = np.arange(0, rad[boundhalos].shape[0])
   rf = rad[rad<radius]
   mf = mass[rad<radius]
   xf = x[rad<radius]
   yf = y[rad<radius]
   zf = z[rad<radius]
   vxf = vx[rad<radius]
   vyf = vy[rad<radius]
   vzf = vz[rad<radius]
   listy = [hlist, rf, mf, xf, yf, zf, vxf, vyf, vzf]
   #creates a single numpy array and saves it into "halos" directory 
   output = np.array(listy)
   print output.shape
   properorder = np.argsort(output[1])
   newhalos = output[:,properorder]
   newhalos[0] = np.arange(0, newhalos.shape[1])
   targg = str(newhalos[2,0])
   np.save("halos/mass"+targg[0]+targg[2:]+"quadrant"+str(quadrant)+".npy", newhalos)


def construct(radius, lmass, umass, halos, quadrant):
   h = 0.7
   lbound = ubound = 0
   halolist = []

   #sets up lower and upper bounds of the list
   for i in range(0,halos.shape[1]):
      if (halos[1,i] < lmass):
         lbound = i
      elif (halos[1,i] < umass):
         ubound = i
      else:
         break

   #sweeps through all the indices  in order to find a halo which
   #matches necessary conditions and then saves them
   omithalo = []
   while (lbound < ubound):
      #print "lbound:", lbound
      #print "ubound:", ubound
      #print "to go", ubound-lbound
   
      halo = find_halo(lbound, ubound, halos, omithalo)

      lbound += 1   

      #signifies the whole list was searched to no avail
      if (halo == -1):
         continue
      else:
         omithalo.append(halo)

      save_halo(radius, halos, halo, quadrant)

      #continues to loop after conditional b/c there may be two higher
      #mass halos which share a fractional difference with each other      


####################MAIN#######################
for i in range(0, 20):
   print "BEGINNING OF TRIAL NUMBER", i+1, "###################################################"
   if (i < 10):
      halos = np.load("quadrants/halos_0214_0"+str(i)+".npy")
   else:
      halos = np.load("quadrants/halos_0214_"+str(i)+".npy")     

   construct(10, 10**12, 10**13, halos, i)

