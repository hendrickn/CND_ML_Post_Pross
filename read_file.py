from scipy.io import loadmat
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import svd
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import copy
import math
from matrix import*
from scipy.signal import hilbert
import time
import cProfile
from numba import vectorize, cuda

#tf.debugging.set_log_device_placement(True)

y=[i for i in range(520)]
x=[i for i in range(37)]
dt=9.62916600773235e-12
c=3*10**8/math.sqrt(5.63)
def match_filter(Res) :
    #cp = cProfile.Profile()
    #cp.enable()
    dx = 0.02
    nip = 37
    t_sample=np.array([i*dt for i in range(520)])
    rang = t_sample*c/2
    antx =[i*dx for i in range(nip)]
    pix_x = antx
    pix_y = rang
    pixels = np.zeros((len(pix_y), len(pix_x)))
    delay = np.zeros(pixels.shape)
    f = 2e9
    a = math.sqrt(2) * math.pi * f
    dd=-2 * math.sqrt(math.exp(1) / (2* a ** 2))*( a ** 2 )
    for ii in range(nip) :
        t0 = time.time()
        print('Pulse '+str(ii+1)+ ' of '+str(nip)+'\n')
        for jj in range(520) :
            for k in range(nip) :
                t = math.sqrt((antx[k] - pix_x[ii])** 2 + pix_y[jj] ** 2) * 2 / c
                delay[:, k] = t_sample-t
            b = delay-1/f
            #d = -2 * math.sqrt(math.exp(1) / (2* a ** 2))*( a ** 2 )* b
            #A=np.multiply(Res,np.multiply(d,np.exp(-(a**2)*(np.multiply(b,b)))))
            A=dot_mul(Res,dot_mul(dd* b, dot_exp(-(a**2)*(dot_mul(b,b)))))
            pixels[jj, ii] = np.sum(A)
        print(time.time()-t0)
    #cp.disable()
    #cp.print_stats()
    return pixels


def time_crop(Bscan) :
    BBscan=copy.deepcopy(Bscan)
    Y=[entropy(BBscan,k) for k in range(len(BBscan))]
    Yp = [(Y[i + 1] - Y[i]) / dt for i in range(len(Y) - 1)]
    ypmax=np.max(Yp)
    i0=Yp.index(ypmax)
    for i in range(i0+1) :
        for k in range(len(Bscan[0])) :
            BBscan[i,k]=0

    j=i0
    while abs(Yp[j])<10**(8) :
        for k in range(len(Bscan[0])):
            BBscan[j, k] = 0
        j+=1
    return BBscan

def read_B_scan(i) :
    Matlab = loadmat('3DGR\\raw_data')
    return Matlab['raw_data'][:, :, i]


def plot_SVD(B_scan) :
    U,s,V=svd(B_scan)
    X=[i for i in range(len(s))]
    plt.bar(X,s)
    return U,s,V


def PCA_data(Bscan,num) :
    U,s,V=svd(Bscan)
    I=np.zeros(Bscan.shape)
    for  i in range(num) :
        I[i][i]=s[i]
    G=np.subtract(Bscan,np.dot(np.dot(U,I),V))
    return G


def mean_data(Bscan) :
    BBscan=copy.deepcopy(Bscan).conj()
    print(Bscan.max())
    Mean=np.array([[np.mean(i)]*(len(Bscan[0])) for i in Bscan])
    print(Mean.max())
    BBscan=np.subtract(BBscan,Mean)
    return BBscan


def plot_color_map(Bscan1,Bscan2) :
    BBscan1=np.array([[Bscan1[i,j] for j in range(len(Bscan1[0]))] for i in range(len(Bscan1)-1,-1,-1)])
    BBscan2=np.array([[Bscan2[i,j] for j in range(len(Bscan2[0]))] for i in range(len(Bscan2)-1,-1,-1)])
    levels = MaxNLocator(nbins=1000).tick_values(BBscan1.min(), BBscan1.max())
    fig, (ax0,ax1) = plt.subplots(ncols=2)
    cmap = plt.get_cmap('seismic')
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    im = ax0.pcolormesh(x, y, BBscan1, cmap=cmap, norm=norm)
    fig.colorbar(im, ax=ax0)
    ax0.set_title('Bscan1')
    levels2=MaxNLocator(nbins=1000).tick_values(BBscan2.min(), BBscan2.max())
    norm2=BoundaryNorm(levels2, ncolors=cmap.N, clip=True)
    imm=ax1.pcolormesh(x,y,BBscan2,cmap=cmap,norm=norm2)
    fig.colorbar(imm,ax=ax1)
    ax1.set_title('Bscan2')


def epsilon(Bscan) :
    BBBscan=np.zeros(Bscan.shape)
    for i in range(len(Bscan)) :
        Norm=np.sum(dot_mul(Bscan,Bscan)[i,:])
        if Norm<10**(-11) :
            for j in range(len(Bscan[0])) :
                BBBscan[i,j]=0
        else :
            for j in range(len(Bscan[0])) :
                BBBscan[i,j]=(abs(Bscan[i,j])**2)/Norm
    return BBBscan

def entropy(Bscan,k) :
    BG=epsilon(Bscan)
    ent=np.sum([-BG[k,i]*math.log(BG[k,i]) if BG[k,i]!=0 else 0 for i in range(len(Bscan[0]))])
    return ent

def gate(Bscan,alpha) :
    #t0=ent.index(np.min(ent))
    BCscan=np.zeros(Bscan.shape)
    threshold=alpha*math.log(len(Bscan[0]))
    for t in range(len(Bscan)) :
        if entropy(Bscan,t)<threshold :
            BCscan[t,:]=copy.deepcopy(Bscan[t,:])
    return BCscan

def epsw(Bscan,alpha) :
    BCscan=np.zeros(Bscan.shape)
    for i in range(len(BCscan)) :
        for j in range(len(BCscan[0])) :
            BCscan[i,j]=epsilon(Bscan)[i,j]*gate(i,Bscan,alpha)
    return BCscan


#plot_color_map(read_B_scan(19),epsilon(read_B_scan(19)))
Bscan=read_B_scan(9)
X=[i for i in range(len(Bscan))]
BBscan=np.array([[Bscan[i][j] if i>125 else 0 for j in range(37)]for i in range(520)])
#plot_color_map(read_B_scan(19),epsilon(BBscan))
# Y=[entropy(Bscan,k) for k in range(len(Bscan))]
# Yp=[(Y[i+1]-Y[i])/dt for i in range(len(Y)-1)]
# Xp=[i for i in range(len(Y)-1)]
# plt.subplot(1,2,1)
# plt.plot(X[150:],Y[150:])
# plt.subplot(1,2,2)
# plt.plot(Xp,Yp)
#plot_color_map(read_B_scan(1),abs(hilbert(match_filter(PCA_data(read_B_scan(1),2)))))
aa=[0.6, 0.7, 0.8]
# for a in aa :
#     plot_color_map(read_B_scan(15),abs(hilbert(match_filter(gate(read_B_scan(15),a)))))
#         #
#         # plot_color_map(read_B_scan(9),time_crop((read_B_scan(9))))
# plt.show()






# plot_color_map(read_B_scan(17),PCA_data(read_B_scan(17),2))
#
# plot_color_map(read_B_scan(17),epsilon(read_B_scan(17)))
# plt.show()

#
# plot_color_map(read_B_scan(17),abs(hilbert(match_filter(PCA_data(read_B_scan(17),2)))))
#
# plot_color_map(read_B_scan(17),abs(hilbert(match_filter(epsilon(read_B_scan(17))))))
