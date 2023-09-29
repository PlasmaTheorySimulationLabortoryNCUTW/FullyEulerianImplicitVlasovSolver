import h5py
import numpy as np
import io
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.style.use("presentation.mplstyle")
from matplotlib import cm, colors,patches
import mpl_toolkits.mplot3d.axes3d
from PIL import Image
from scipy.signal import hilbert,savgol_filter
from PIL import Image
def envelope_fig(data,t,frames,x_left,x_right):
    bio = io.BytesIO()
    plt.close()
    x        = np.linspace(0,4*np.pi,257)
    fig = plt.figure()
    plt.imshow(data.T,origin = 'lower',aspect = 'auto',norm=colors.SymLogNorm(linthresh = 1e-5,vmin=1e-3,base=10),extent=[0,4*np.pi,-10 ,10], cmap = 'gray')
    plt.colorbar()
    plt.xlabel(r'$x/\lambda_{De}$')
    plt.ylabel(r'$v/v_{th}$')
    plt.title(r'$t\omega_{pe}=$'+str(round(t,1)))
    plt.savefig(bio,facecolor = 'white')
    frames.append(Image.open(bio))
listf    = open('file_list_f.txt', 'r')
namelist = listf.read().split()
frames   = []
dump     = 4
x        = np.linspace(0,4*np.pi,257)
dx       = x[2]-x[1]
v        = np.linspace(-10,10,256)
fv_t     = np.zeros((len(namelist)//dump+1,256))
t        = np.zeros(len(namelist)//dump+1)
for i,name in enumerate(namelist):
    filename = "./data/"+name
    with h5py.File(filename, "r") as f:
        fvx             = np.array(f['fuex']).reshape((257,256))
        axis            = np.array(f['axis'])
        t[i//dump]      = axis[11]
        envelope_fig(fvx,t[i//dump],frames,x[ind]-200,x[ind]+600)
frames[0].save('dis_L.gif',save_all = True,append_images = frames[1:],duration=500,loops = 0)
