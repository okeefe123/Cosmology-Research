
# coding: utf-8

# In[1]:

import numpy as np


# In[34]:
input = raw_input("choose halo (no.npy)\n")
halos = np.load(input+'.npy')


# In[48]:

halos.shape[1]
#bigarray = np.arange(0,halos.shape[1]) + 1
#halos[0] = bigarray
print halos[0]
partner = np.argmax(halos[2,1:200])
# In[49]:

#mass
#print halos[2,0:10]
#print halos[:,0]
#print halos[:,4]


# In[50]:

halos[:,[0,partner]] = halos[:,[partner,0]]


# In[51]:

#print halos[:,0]
#print halos[:,4]


# In[52]:

H_0 = 70
x = halos[3,:]
y = halos[4,:]
z = halos[5,:]

# In[55]:

#Reconfigure the radii such that they reflect the new mass as target halo
halos[1] = ((x-x[0])**2 + (y-y[0])**2 + (z-z[0])**2)**.5

#print halos[0,0:10]
#print halos[1,0:10]


# In[56]:

properorder = np.argsort(halos[1])
haloshift = halos[:,properorder]
np.save(input+'shift.npy', haloshift)


# In[ ]:



