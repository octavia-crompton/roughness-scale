import numpy as np

def reshape_topo(topo, p):
    xx, yy, zz = topo[:,0],topo[:,1],topo[:,2],
    xx = xx.reshape(p.Nxcell, p.Nycell)
    yy = yy.reshape(p.Nxcell, p.Nycell)
    zz = zz.reshape(p.Nxcell, p.Nycell)
    
    return xx, yy, zz


def xgrad(zed):
    return zed[1:, :] - zed[:-1, :]

def ygrad(zed):
    return zed[:, 1:] - zed[:, :-1]

def comp_y_slope(zed, p):
    return np.mean(zed[:, 1:] - zed[:, :-1])/p.dy

def equivalent_y(zed, yy, p):
    
    slope = comp_y_slope(zed, p)
    
    planar = yy*slope
    
    return planar

def detrend_y(zed, yy, p):
    
    planar = equivalent_y(zed, yy, p)

    yp = zed - planar
    
    return yp

def plot_y_transct(zed):
    plt.plot(zed[1])
    plt.plot(zed[-1])    
    plt.plot(zed[int(zed.shape[0]/2)], '--') 
    
