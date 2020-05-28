import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return -2* (pow(x[0], 2)) +20*x[0]- (pow(x[1], 2))+16*x[1]-3* (pow(x[2], 4)) -48* (pow(x[2], 3)) -288* (pow(x[2], 2)) -768*x[2]-878
#       -2* (pow(x[0], 2)) +20*x[0]- (pow(x[1], 2))+16*x[1]-3* (pow(x[2], 4)) -48* (pow(x[2], 3)) -288* (pow(x[2], 2)) -768*x[2]-878
x,y=np.mgrid[-20:20:0.1, -20:20:0.1]
fig,ax0=plt.subplots()
fig,ax1=plt.subplots()
fig,ax2=plt.subplots()
lev_region=[-10000,-9000,-6000,-4000,-2000,-1500,-1000,-750,-500,-250,-150,-100,-20,-10,-5]

cs0=ax0.contour(x, y, f([x,y,0]),levels=lev_region,colors='k')
cs1=ax1.contour(x, y, f([x,0,y]),levels=lev_region)
cs2=ax2.contour(x, y, f([0,x,y]),levels=lev_region)

ax0.clabel(cs0)
ax1.clabel(cs1)
ax2.clabel(cs2)

inf=[]
with open("ans1.dat") as tr:
    for line in tr:
        inf.append([float(line.split(',')[0]),float(line.split(',')[1]),float(line.split(',')[2])])
inf=np.array(inf)
ax0.plot(inf[:,0],inf[:,1],"go:")
ax1.plot(inf[:,0],inf[:,2],"go:")
ax2.plot(inf[:,1],inf[:,2],"go:")

ax0.set_title('Linii urovneniya & \n traektoriay poiska \n x1,x2')
ax1.set_title('Linii urovneniya & \n traektoriay poiska \n x1,x3')
ax2.set_title('Linii urovneniya & \n traektoriay poiska \n x2,x3')

plt.show()
