import scipy as sy
import pylab as plb 
import matplotlib 
from pylab import *

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np

from scipy.interpolate import interp2d
import matplotlib.mlab as ml

from mpl_toolkits.axes_grid1 import make_axes_locatable

from matplotlib.mlab import griddata
from matplotlib.ticker import AutoMinorLocator

from matplotlib.colors import LogNorm
#-----------------------------------------------------------------------

matplotlib.rc('font', family='serif')

matplotlib.rc('text.latex', 
              preamble=[r'\usepackage[T1]{fontenc}',
                        r'\usepackage{amsmath}',
                        r'\usepackage{txfonts}',
                        r'\usepackage{textcomp}'])

matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.unicode'] = True   
#-----------------------------------------------------------------------
#---------------------------------------------------------------
#-----------------------------------------------------------
data = plb.loadtxt('band1_interpolated.txt')

kxn1= data[:,0]
kyn1= data[:,1]
dispn1= data[:,2]
#-----------------------------------------------------------Sep2test5.txt
dispn1=dispn1*1000.0
#-----------------------------------------------------------
#-----------------------------------------------------------
data = plb.loadtxt('Sep2test5.txt')

kxedge= data[:,0]
kyedge= data[:,1]
#-----------------------------------------------------------
# Definition of Matrices:
#-----------------------------------------------------------
dim1=500
dim2=500

k=0

A = np.zeros((dim2,dim1))
xmat = np.zeros((dim2,dim1))
ymat = np.zeros((dim2,dim1))


for j in range(0,dim2):
    for i in range(0,dim1):

        A[j, i] = dispn1[k]
        xmat[j, i] = kxn1[k]
        ymat[j, i] = kyn1[k]
        k = k + 1
            
#-----------------------------------------------------------


x = kxn1 #np.arange(1, 10)
y = kyn1 #x.reshape(-1, 1)
z = dispn1
h = A #x * y

#x = np.arange(1, 10)
#y = x.reshape(-1, 1)
#h = x * y

#cs = plt.contourf(h, levels=[-0.8, -0.6, -0.4])
#    colors=['#808080', '#A0A0A0', '#C0C0C0'], extend='both')
#cs.cmap.set_over('red')
#cs.cmap.set_under('blue')
#cs.changed()

cm = plt.cm.get_cmap('tab20b')

CS = plt.contour(ymat,xmat,h,20, cmap = cm, linewidth=5)

#CS = plt.contour(ymat,xmat,h,[2.0,2.3,2.6], cmap = cm, linewidth=5)


#plt.clabel(CS, fontsize=9, inline=1)

#cbar = plt.colorbar()

plt.plot(kyedge, kxedge,'k-', linewidth=2)


#gca().set_xlim(0.5,0.75)
#gca().set_ylim(0.15,0.18)

savefig('contour.png')
plt.show()
