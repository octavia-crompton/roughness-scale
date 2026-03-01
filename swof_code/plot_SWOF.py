import matplotlib.pylab as plt
import pandas as pd
import numpy as np
from matplotlib import cm
import matplotlib.colors as colors
from mpl_toolkits.mplot3d import Axes3D
import os
from shapely import geometry
import matplotlib as mpl
mpl.rcParams['figure.dpi']= 100
mpl.rcParams['figure.figsize']= (5,3)
import cv2
from matplotlib.ticker import FormatStrFormatter
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import seaborn as sns

BrBG_big = cm.get_cmap('BrBG_r', 512)
browns = ListedColormap(BrBG_big(np.linspace(0.5, 0.95, 256)), 'browns')

mako_big = cm.get_cmap('mako_r', 512)
mako_r = ListedColormap(mako_big(np.linspace(0.1, 0.9, 256)), 'mako_r')

mako_big = cm.get_cmap('mako', 512)
MAKO = ListedColormap(mako_big(np.linspace(0.15, 0.9, 256)), 'mako')

BrBG_big = cm.get_cmap('BrBG_r', 512)
browns = ListedColormap(BrBG_big(np.linspace(0.5, 0.95, 256)), 'browns')

blues_big = cm.get_cmap('Blues', 512)
blues = ListedColormap(blues_big(np.linspace(0.1, 0.9, 256)))

terrain = cm.get_cmap('BrBG_r', 512)
terrain = ListedColormap(terrain(np.linspace(0.46,0.9, 256)), 'terrain')

GnBu = cm.get_cmap('GnBu', 512)
GnBu = ListedColormap(GnBu(np.linspace(0.1,0.9, 256)), 'GnBu')

cividis = cm.get_cmap('cividis', 512)
cividis = ListedColormap(cividis(np.linspace(0.,0.9, 256)), 'cividis')

sns.set_context('notebook', font_scale=1,
                rc={'lines.linewidth': 1})

long_names = {'So': 'Slope angle, So', 
        'p': 'rain intensity, p (cm/hr)',
        'tr': 'rain duration, tr (min)',
        'alpha_b': 'bare soil n, alpha_b',         
        'alpha_v': 'vegetation n, alpha_v',                  
        'Ks_b': 'bare soil Ksat, Ks_b (cm/hr)',
        'Ks_v': 'vegetation Ksat, Ks_v (cm/hr)',
        'fV'   : 'vegetation fraction, fV',
        'sigma'   : 'vegetation lengthscale, sigma',
        'l'   : 'hillslope length, l (m)',
        'dx'   : 'grid size, dx (m)',
        'seed'   : 'random seed',
        'nchecks'   : 'number of checks',
        'patch_scale'   : 'patch scale',
        }


def get_colors():
    
    fig= plt.figure(figsize = (0, 0))
    ax = plt.gca()#visible = 0);
    
    global blue
    global orange
    global green

    l = ax.plot(np.arange(10), visible = 0)
    blue = l[0].get_color()

    l = ax.plot(np.arange(10), visible = 0)
    orange = l[0].get_color()

    l = ax.plot(np.arange(10), visible = 0)
    green = l[0].get_color()

get_colors()
def mako_r(n):
    return sns.color_palette('mako_r', n)        


def mako(n):
    return sns.color_palette('mako', n)        


def print_input_params(sim_list):
    '''

    Parameters
    ----------
    directory_name
    '''
    s = pd.DataFrame(sim_list)
    for col in s:
        vals= str(s[col].unique()).strip(']').strip('[')
        vals = ', '.join(vals.split(' ')).replace(', ,', ', ').replace('  ', ' ')
        if col in long_names and col != 'seed':
            print( '{0}: {1}'.format(long_names[col], vals) )
        elif col != 'seed':
            print( '{0}: {1}'.format(col, vals) )


def print_input_params_sample(sim_list):
    '''

    Parameters
    ----------
    directory_name
    '''


    s = pd.DataFrame(sim_list)
    for col in s:
        vals= [ str(s[col].min()), str(s[col].max())]
        vals = ', '.join(vals)
        
        # print( '{0}: {1}'.format(names[col], vals) )
        if col in names:
            print( '{0}: {1}'.format(names[col], vals) )
        else:
            print( '{0}: {1}'.format(col, vals) )
        

def plot_IhU(sim, ucut = 0, dcut = 0, suptitle = True):
    '''
    '''
    # visualize roughness and end of rain water height for a single storm instance
    # downslope is to the right
    
    fig, axes = plt.subplots(1, 3, figsize=(20, 4))
    xc, yc, zc, infl_2d = trim(sim, sim.infl_2d, ucut, dcut)
    xc, yc, zc, h_max = trim(sim, sim.h_max, ucut, dcut)
    xc, yc, zc, U_max = trim(sim, sim.vc[sim.i_tr], ucut, dcut)
    
    p = axes[0].pcolormesh(infl_2d*100, cmap = 'Greens')
    plt.colorbar(p, ax = axes[0])
    axes[0].set_title('Cumulative infiltration $I$')

    p = axes[1].pcolormesh(h_max*100, cmap = 'Blues')
    plt.colorbar(p, ax = axes[1])
    axes[1].set_title('maximum $h$ (cm)')

    p = axes[2].pcolormesh(U_max*100, cmap = 'mako')
    axes[2].set_title('maximum $U$ (cm/s)')
    plt.colorbar(p, ax = axes[2])

    if suptitle:
        t = plt.suptitle(sim.name)
        t.set_y(0.0)
        t.set_x(0.2)

    return fig


def plot_zhU(sim, ucut = 0, dcut = 0, suptitle = True):
    '''
    '''
    # visualize roughness and end of rain water height for a single storm instance
    # downslope is to the right
    
    fig, axes = plt.subplots(1, 3, figsize=(20, 4))
    xc, yc, zc, zc = trim(sim, sim.zc, ucut, dcut)
    xc, yc, zc, h_max = trim(sim, sim.h_max, ucut, dcut)
    xc, yc, zc, U_max = trim(sim, sim.vc[sim.i_tr], ucut, dcut)
    
    p = axes[0].pcolormesh(zc)
    plt.colorbar(p, ax = axes[0])
    axes[0].set_title('Topography $z$')

    p = axes[1].pcolormesh(h_max*100, cmap = 'Blues')
    plt.colorbar(p, ax = axes[1])
    axes[1].set_title('maximum $h$ (cm)')

    p = axes[2].pcolormesh(U_max*100, cmap = 'mako')
    axes[2].set_title('maximum $U$ (cm/s)')
    plt.colorbar(p, ax = axes[2])

    if suptitle:
        t = plt.suptitle(sim.name)
        t.set_y(0.0)
        t.set_x(0.2)

    return fig

def plot_vIU(sim, suptitle = True):
    '''
    '''
    # visualize vegetation pattern and end of rain water height for a single storm instance
    # downslope is to the right
    
    fig, axes = plt.subplots(1, 3, figsize=(20, 4))

    p = axes[0].pcolormesh(sim.veg*1., cmap = 'Greens')
    plt.colorbar(p, ax = axes[0])
    axes[0].set_title('Vegetation cover')


    p = axes[1].pcolormesh(sim.infl_2d*100, cmap = 'Blues')
    plt.colorbar(p, ax = axes[1])
    axes[1].set_title('Cumulative infiltration $I$')

    p = axes[2].pcolormesh(sim.h_max*100, cmap = 'Blues')
    plt.colorbar(p, ax = axes[2])
    axes[2].set_title('maximum $h$ (cm)')

    if suptitle:
        t = plt.suptitle(sim.name)
        t.set_y(0.0)
        t.set_x(0.2)

    return fig


def plot_Ihq(sim, suptitle = True):
    '''
    '''
    # visualize roughness and end of rain water height for a single storm instance
    # downslope is to the right
    
    fig, axes = plt.subplots(1, 3, figsize=(20, 4))
    p = axes[0].pcolormesh(sim.infl_2d*100, cmap = 'Greens')
    plt.colorbar(p, ax = axes[0])
    axes[0].set_title('Cumulative infiltration $I$')

    p = axes[1].pcolormesh(sim.h_max*100, cmap = 'Blues')
    plt.colorbar(p, ax = axes[1])
    axes[1].set_title('maximum $h$ (cm)')

    p = axes[2].pcolormesh(sim.q_max*1e4, cmap = 'mako')
    axes[2].set_title('maximum $q$ (cm^2/s)')
    plt.colorbar(p, ax = axes[2])

    if suptitle:
        t = plt.suptitle(sim.name)
        t.set_y(0.0)
        t.set_x(0.2)

    return fig

##### Hysteresis functions #########


def plot_obs_hysteresis(sim, suptitle = True):
    '''
    Q(L) vs dQ(L)/dt hysteresis
    '''
    f = int(10/np.unique(np.diff(sim.time))[-1].round(2))
    dq_1d = np.abs(np.diff(np.diff(sim.Vol_bound_tot[::f])))
    dq_1d = np.append( [0], dq_1d)
    
    q_1d = np.diff(sim.Vol_bound_tot[::f])

    fig = plt.figure(figsize = (5, 3))

    plt.scatter(q_1d/np.max(q_1d), dq_1d/np.max(dq_1d), c = sim.time[::f][1:]/60, cmap = 'mako')
    plt.colorbar(label = 'time')

    if suptitle:
        t = plt.suptitle(sim.name, fontsize = 12)
        t.set_y(-0.12)
        t.set_x(0.2)

    poly = geometry.Polygon(zip(q_1d/q_1d.max(), dq_1d/dq_1d.max()))        
    plt.title(r'hysteresis $\eta_{1}$= {0:.2}'.format( poly.area, '{obs}'), fontsize = 12)
    plt.xlabel(r'$ q(L) $')
    plt.ylabel(r'$ dq(L)/dt $')
    
    return fig


def plot_hydro_hysteresis(sim, suptitle = True):
    '''
    <h> vs <q> hysteresis
    '''
    
    f = int(20/np.unique(np.diff(sim.time))[-1].round(2))
    h_1d = sim.Vol_of_tot[::f]

    q_1d = np.diff(sim.Vol_bound_tot[::f])
    q_1d = np.append( [0],q_1d)
    fig = plt.figure(figsize = (5, 3))

    plt.scatter(h_1d/np.max(h_1d), q_1d/np.max(q_1d), c = sim.time[::f]/60, cmap = 'mako')
    plt.colorbar(label = 'time')

    if suptitle:
        t = plt.suptitle(sim.name, fontsize = 12)
        t.set_y(-0.12)
        t.set_x(0.2)

    try:
        plt.title('hysteresis = {0:.2}'.format(sim.hydro_hyst_orient), fontsize = 12)
    except:
        poly = geometry.Polygon(zip(h_1d/h_1d.max(), q_1d/q_1d.max()))        
        plt.title('hysteresis = {0:.2}'.format( poly.area), fontsize = 12)
    
    plt.xlabel(r'$\langle h \rangle$')
    plt.ylabel(r'$q(L)$')
    
    return fig

def plot_dh_hysteresis(sim, suptitle = True):
    '''
    d<h>/dt vs q hysteresis
    '''
    
    f = int(10/np.unique(np.diff(sim.time))[-1].round(2))
    h_1d = np.diff(sim.Vol_of_tot[::f])
    q_1d = np.diff(sim.Vol_bound_tot[::f])

    fig = plt.figure(figsize = (5, 3))

    plt.scatter( q_1d/np.max(q_1d), h_1d/np.max(h_1d), c = sim.time[::f][1:]/60, 
                s =  sim.time[::f][1:]/10,
                cmap = 'mako')
    plt.colorbar(label = 'time')

    if suptitle:
        t = plt.suptitle(sim.name, fontsize = 12)
        t.set_y(-0.12)
        t.set_x(0.2)

    plt.ylabel(r'$d(\langle h \rangle)/dt$')
    plt.xlabel(r'$q(L)$')
    
    return fig



def plot_q_hysteresis(sim, suptitle = True):
    '''
    h vs q hysteresis
    '''
    h_1d = sim.hc.mean(1).mean(1)
    u_1d = sim.uc.mean(1).mean(1)
    v_1d = sim.vc.mean(1).mean(1)
    q_1d = sim.qc.mean(1).mean(1)    
    fig = plt.figure(figsize = (5, 3))
    plt.scatter(h_1d/h_1d.max(), q_1d/q_1d.max(), c = sim.t/60, cmap = 'mako')
    plt.colorbar(label = 'time')
    
    if suptitle:
        t = plt.suptitle(sim.name, fontsize = 12)
        t.set_y(-0.12)
        t.set_x(0.2)

    poly = geometry.Polygon(zip(h_1d/h_1d.max(), q_1d/q_1d.max()))        
    plt.title('hysteresis = {0:.2}'.format( poly.area), fontsize = 12)
    plt.xlabel(r'$\langle h \rangle$')
    plt.ylabel(r'$\langle q \rangle$')
    return fig


def plot_v_hysteresis(sim, suptitle = True):
    '''
    h vs q hysteresis
    '''
    h_1d = sim.hc.mean(1).mean(1)
    u_1d = sim.uc.mean(1).mean(1)
    v_1d = sim.vc.mean(1).mean(1)
    q_1d = sim.qc.mean(1).mean(1)    
    fig = plt.figure(figsize = (5, 3))
    plt.scatter(h_1d/h_1d.max(), v_1d/v_1d.max(), c = sim.t/60, cmap = 'mako')
    plt.colorbar(label = 'time')
    
    if suptitle:
        t = plt.suptitle(sim.name, fontsize = 12)
        t.set_y(-0.12)
        t.set_x(0.2)

    poly = geometry.Polygon(zip(h_1d/h_1d.max(), v_1d/v_1d.max()))        
    plt.title('hysteresis = {0:.2}'.format( poly.area), fontsize = 12)
    plt.xlabel(r'$\langle h \rangle$')
    plt.ylabel(r'$\langle v \rangle$')
    return fig


def hydrographs(summary):
    fig = plt.figure(figsize = (5,3))
    for key in summary.index:    

        s = summary.loc[key]

        plt.plot( np.diff(s.Vol_bound_tot[::60]) , label = np.round(s.fV, 2)       )
    return fig

######################## Wrap plots ############################################

def loop_summary(function, summary, name, out_dir, **kwargs):
    '''
    apply the same plotting functionality to all simulations in summary
    '''  
    for key in summary.index:
        
        sim = summary.loc[key]
        fig = function(sim, **kwargs);        

        path = save_fig(fig, key, name, out_dir)    

    return os.path.dirname(path)

def save_fig(fig, key, name, out_dir):
    '''
    '''
    figdir = os.path.join(out_dir, 'figures')
    if os.path.isdir(figdir) == 0:
        os.mkdir(figdir)

    key = key.replace('/',',')
    simdir = os.path.join(figdir, key)
    if os.path.isdir(simdir) == 0:
        os.mkdir(simdir)
        
    fig.savefig(get_path(key, name, out_dir), bbox_inches = 'tight')

    path = get_path(key, name, out_dir)
    return (os.path.dirname(path))



def save_summ_fig(fig, name, base_dir):
    '''
    '''
    figdir = os.path.join(base_dir, 'figures')
    if os.path.isdir(figdir) == 0:
        os.mkdir(figdir)

    path = os.path.join(base_dir, 'figures',   name + '.png')
        
    fig.savefig(path, bbox_inches = 'tight')


def get_path(key, name, base_dir):

    return os.path.join(base_dir, 'figures', key,   name + '.png')


############################# General plots ###################################


def plot_hydrographs(summary,  ax=None, label = '', time_res = 30, c = 'b', alpha = 0.3):
    '''
    hydrographs (flux3) with negative values removed

    Parameters
    ----------
    ax
    nonzero
    summary : pd.DataFrame
        SVE simulations

    '''
    if not ax:
        fig, ax = plt.subplots(1, figsize = (7, 3))
    else:
        fig = plt.gcf()

    for key in summary.index:
        sim = summary.loc[key]

        f = int(max(time_res/np.unique(np.diff(sim.time))) )
        dt =  np.mean(np.diff(sim.time[::f]))
        if label == '':
            ax.plot(sim.time[::f][1:]/60, np.diff(sim.Vol_bound_tot[::f])/dt*3.6e5/sim.l/sim.L, c = c, alpha = alpha)
        else:
            ax.plot(sim.time[::f][1:]/60, np.diff(sim.Vol_bound_tot[::f])/dt*3.6e5/sim.l/sim.L,  
                    label = sim[label], c = c, alpha = alpha)

    ax.set_xlabel('minutes')
    ax.set_ylabel('cm/hr')
    if label:
        plt.legend()
        
    return fig, ax
######################## 3D plots ############################################


def fix_3D_axes(ax):
    '''
    Set up 3D axes (and get rid of colored axis planes)
    '''
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.xaxis.pane.set_edgecolor('w')
    ax.yaxis.pane.set_edgecolor('w')
    ax.zaxis.pane.set_edgecolor('w')
    ax.set_xticks([])
    ax.set_zticks([])
    ax.set_yticks([])
    ax.grid(False)

    return ax

def format_plot(xfld, yfld, cfld, subset, ax = None,
                inds = None, vmax = 80, vmin = 0, 
                legend = 0, alpha = 1, legend_title = None,
                cmap = cm.ocean):
    
    if ax is None:
        ax = plt.gca()
        
    sns.scatterplot(data = subset,
                    x = xfld, y = yfld,
                    hue = cfld, palette = cmap, ax = ax, s = 70,
                    alpha = alpha, vmin = vmin, vmax = vmax)

    ax.set_xlabel(format_name(xfld))
    ax.set_ylabel(format_name(yfld))
    
    if legend:
        ax.legend(title = legend_title)
    else: 
        ax.legend([])
        ax.get_legend().remove()
        

def format_colorbar(vmax, cmap, cbaxes, fmt, label = 'cm'  , alpha = 1, vmin = 0):
    
    from matplotlib.ticker import FuncFormatter
    
    a = np.array([[np.round(vmin*100,2), np.round(vmax*100, 2)]])
    plt.figure(figsize=(0, 0))
    img = plt.imshow(a, cmap=cmap)
    plt.gca().set_visible(False)
    cbar = plt.colorbar(img, cax = cbaxes, format=FuncFormatter(fmt),  alpha = alpha, label = label)
    return cbar
    
def plot_surface(sim, fld, title = '', color = cm.Greens,
                plot_veg = False, erode = False, vmin = None, vmax = None,
                 ax = '', dcut = 0, ucut = 0, stack = 0, alpha = 0.5):
    '''
    '''
    if ax == '':
        fig = plt.figure(figsize=(5,5))
        ax = fig.add_subplot(111, projection='3d')

    ax = fix_3D_axes(ax)
    norm = plt.Normalize()
    xc, yc, zc, veg = trim(sim, sim.veg, dcut, ucut, stack )
    xc, yc, zc, Ic = trim(sim, sim[fld], dcut, ucut, stack )

    if vmin == None:
        vmin = Ic.min()
    if vmax == None:
        vmax = Ic.max()        

    I_norm = colors.Normalize( vmin = vmin,
                             vmax =   vmax
                            )
    I_colors = color(I_norm(Ic))
    
    if dcut == 0:
        veg = np.array(sim['veg']).reshape( sim.Nxcell, sim.Nycell)[:, ucut:]
    else:
        veg = np.array(sim['veg']).reshape( sim.Nxcell, sim.Nycell)[:, ucut:-dcut]
        
    veg_norm = colors.Normalize(vmin =  0. ,
                                vmax =  1.5
                              )
    if erode == True:
        veg = veg - cv2.erode(veg, np.ones((3,3)), 1)
    veg_colors = cm.Greens(veg_norm(veg ))

    ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))

    if plot_veg:
        
        ax.plot_surface(xc, yc , zc ,
                        facecolors = veg_colors, rstride=1, cstride=1,
                        linewidth=0, antialiased=True, shade=False, alpha=1)

    ax.plot_surface(xc, yc , zc ,
                    facecolors=I_colors, rstride=1, cstride=1,
                    linewidth=0, antialiased=True, shade=False, alpha=alpha)   
  
    ax.view_init(25, 20)

    t = ax.set_title(title, fontsize  = 14)
    t.set_y(0.9)
    
    return  ax

def trim(sim, array, dcut = 0, ucut = 0, stack= 0, offset = 0):
    
    if dcut == 0:
        xc = np.array(sim.xc)[:, ucut:]
        yc = np.array(sim.yc)[:, ucut:]
        zc = np.array(sim.zc)[:, ucut:]
        array = array[:, ucut:] + offset
        
    else:      
        xc = np.array(sim.xc)[:, ucut:-dcut]
        yc = np.array(sim.yc)[:, ucut:-dcut]
        zc = np.array(sim.zc)[:, ucut:-dcut]
        array = array.reshape( sim.Nxcell, sim.Nycell)[:, ucut:-dcut]  + offset
        
    if stack == 1:
        xc = np.vstack((xc - sim.dx/2, xc + sim.dx/2))
        yc = np.vstack((yc,yc))
        zc = np.vstack((zc,zc))
        array = np.vstack((array,array))
    return xc,yc,zc, array


def rescale (zed):
    ''' 
    rescales an input array from 0 to 1
    ''' 
    return (zed- zed.min())/(zed.max() - zed.min())


def plot_h_surface(sim, hscale = 1, dcut = 0, ucut = 0, title = None, project_veg = 0):
    
    fig = plt.figure(figsize=(5, 5))

    ax = fig.add_subplot(111, projection='3d')
    ax = fix_3D_axes(ax)

    xc, yc, zc, veg = trim(sim, sim.veg, dcut, ucut )
    xc, yc, zc, hc = trim(sim, sim.hc[sim.i_tr], dcut, ucut )
    remove_micro = 1
    if remove_micro == 1:
        zc = rescale(1-yc)*sim.So*sim.l    
    zc = zc - zc.min()
    

    if project_veg:
        
        veg_norm = colors.Normalize(0., 1.1)
        veg_colors = cm.Greens(veg_norm(veg ))    

        ax.plot_surface(xc, yc, zc*0, facecolors = veg_colors , 
                            rstride = 1, cstride = 1, alpha = 0.4,
                            linewidth=0, antialiased=True, shade=False)        

    else:
        veg_norm = colors.Normalize(0., 1.1)
        veg_contour = veg - cv2.erode(veg, np.ones((3,3)), 1)
        veg_colors = cm.Greens(veg_norm(veg_contour ))    

        ax.plot_surface(xc, yc, zc, facecolors = veg_colors , 
                                rstride = 1, cstride = 1, alpha = 1,
                                linewidth=0, antialiased=True, shade=False)    
    
    
    h_norm = colors.Normalize(vmin= hc.min(), vmax= hc.max()*1.2)
    h_colors = cm.Blues(h_norm(hc)) 

    cbaxes = fig.add_axes([0.85, 0.25, 0.03, 0.45]) 
    fmt = lambda x, pos: '{:.1f}'.format(x)
    format_colorbar( hc.max(), 'Blues', cbaxes, fmt, label = '$h$ (cm)' )

    plot =  ax.plot_surface(xc, yc, zc + hc, facecolors = h_colors, 
                        rstride=1, cstride=1,cmap = cm.Blues,alpha = 0.8,
                        linewidth=0,antialiased=True, shade=False),


    ax.view_init(25, 20)
    ax.set_title(title)

    return ax

def plot_U_surface(sim, hscale = 1, dcut = 0, ucut = 0, title = None ):
    
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(111, projection='3d')
    ax = fix_3D_axes(ax)

    xc, yc, zc, veg = trim(sim, sim.veg, dcut , ucut )
    xc, yc, zc, U = trim(sim, sim.vc[sim.i_tr], dcut, ucut )
    xc, yc, zc, hc = trim(sim, sim.hc[sim.i_tr], dcut, ucut )    


    veg_norm = colors.Normalize(0., 1.1)
    veg = veg - cv2.erode(veg, np.ones((3,3)), 1)
    veg_colors = cm.Greys(veg_norm(veg ))    

    ax.plot_surface(xc, yc, zc, facecolors = veg_colors , 
                            rstride = 1, cstride = 1, alpha = 1,
                            linewidth=0, antialiased=True, shade=False)        

    x, y = np.where(U > 1)
    xx = x[x< sim.Nxcell]
    U[xx , y] =   (U[ xx - 1, y] + U[ xx - 1, y])/2

    xx = x[x ==  sim.Nxcell]
    if len(xx > 0):
        U[xx , y] =  U[ xx - 1, y] 

    
    Umax = U.max()*1.05
    Umin = U.min()*0.95
    U_norm = colors.Normalize(vmin= U.min(),  vmax=Umax)  
    U_colors = cm.cividis(U_norm(U))

    cbaxes = fig.add_axes([0.85, 0.25, 0.03, 0.45]) 
    fmt = lambda x, pos: '{:.1f}'.format(x)
    format_colorbar(Umax, cm.cividis, cbaxes, fmt, label = '$|U|$ (cm/s)'  , vmin = Umin)

    plot =  ax.plot_surface(xc, yc, zc + hc, facecolors = U_colors, 
                        rstride=1, cstride=1,cmap = cm.cividis,alpha = 0.8,
                        linewidth=0,antialiased=True, shade=False),

    ax.view_init(25, 20)
    ax.set_title(title)

    return ax    

def plot_I_surface(sim, hscale = 1, dcut = 0, ucut = 0, project_veg = 0, title = None, vmax = None , ax = ''):

    if ax == '':
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(111, projection='3d')


    ax = fix_3D_axes(ax)

    xc, yc, zc, veg = trim(sim, sim.veg, dcut , ucut )
    xc, yc, zc, infl = trim(sim, sim.infl_2d, dcut, ucut )
    xc, yc, zc, hc = trim(sim, sim.hc[sim.i_tr], dcut, ucut )    
    
    remove_micro = 1
    if remove_micro == 1:
        zc = rescale(1-yc)*sim.So*sim.l    
    zc = zc - zc.min()
    

    if project_veg:
        
        veg_norm = colors.Normalize(0., 1.1)
        veg_colors = cm.Greens(veg_norm(veg ))    

        ax.plot_surface(xc, yc, zc*0, facecolors = veg_colors , 
                            rstride = 1, cstride = 1, alpha = 0.4,
                            linewidth=0, antialiased=True, shade=False)        

    elif project_veg == 2:
        veg_norm = colors.Normalize(0., 1.1)
        veg_contour = veg - cv2.erode(veg, np.ones((3,3)), 1)
        veg_colors = cm.Greens(veg_norm(veg_contour ))    

        ax.plot_surface(xc, yc, zc, facecolors = veg_colors , 
                                rstride = 1, cstride = 1, alpha = 1,
                                linewidth=0, antialiased=True, shade=False)  
        
    if vmax:
        I_norm = colors.Normalize(vmin= infl.min(),  vmax=vmax)
    else:      
        I_norm = colors.Normalize(vmin= infl.min(),  vmax=infl.max())  
        
    I_colors = GnBu(I_norm(infl))

    cbaxes = fig.add_axes([0.85, 0.25, 0.03, 0.45]) 
    fmt = lambda x, pos: '{:.1f}'.format(x)

    if vmax:
        format_colorbar( vmax, GnBu, cbaxes, fmt, label = '$I$ (cm)' )
    else:
        format_colorbar( infl.max(), GnBu, cbaxes, fmt, label = '$I$ (cm)' )
        
    plot =  ax.plot_surface(xc, yc, zc + hc, facecolors = I_colors, 
                        rstride=1, cstride=1,cmap = GnBu,alpha = 0.8,
                        linewidth=0,antialiased=True, shade=False),


    ax.view_init(25, 20)
    ax.set_title(title)

    return ax    
## Plots for debug mode outputs


name_format = {
         'Vol_of_tot' : 'overland flow', 
         'Vol_inf_tot' : 'infiltrated',
         'Vol_rain_tot' : 'rain',    
         'Vol_bound_tot' : 'boundary'}


def plot_check_vol(df):
    '''
    Plots volume track plot
    
    Boundary = negative equals inflow     
    '''
    fig = plt.figure(figsize = (6, 3.5));
    for fld in ['Vol_of_tot', 'Vol_inf_tot', 'Vol_rain_tot', 'Vol_bound_tot']:
        plt.plot(df.time, df[fld], label = name_format[fld])

    plt.legend()
    plt.xlabel('time (s)')
    plt.ylabel('volume (m$^3$)')
    
    plt.tight_layout()
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    return fig

def plot_fluxes(df):
    '''
    Plots volume track plot
    '''
    fig = plt.figure(figsize = (6, 3.5));
    for fld in ['Vol_of_tot', 'Vol_inf_tot', 'Vol_rain_tot', 'Vol_bound_tot']:
        plt.plot(df.time[:-1], np.diff(df[fld]), label = name_format[fld])
    plt.legend()
    plt.xlabel('time (s)')
    plt.ylabel('volume (m$^3$)')
    
    plt.tight_layout()
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    return fig


def plot_outflow(b, p, check_vol = None, results = None):
    '''
    '''
    fig = plt.figure(figsize = (7, 4));

    time = b['time']

    plt.plot(time,  - b['left'], label = 'left')
    plt.plot(time,  b['top'], label = 'top')    

    plt.plot(time,  b['right'], label = 'right')
    plt.plot(time,  - b['bottom'], label ='bottom')    

    outflow = - b.left - b.bottom + b.right + b.top
    plt.plot(time, outflow,  ':',label = 'outflow',)

    if check_vol is not None:
        plt.plot(check_vol.time, check_vol.Vol_bound_tot, '--')    
    
    if results is not None:
        plt.title('Boundary outflow volume: {0:.2f} m$^2$'.format(results.boundary_vol))
    
    plt.xlabel('time (s)')
    plt.ylabel('m$^3$')    
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))    
    plt.tight_layout()
    
    return fig


def plot_flux_LR_BT(LR, BT):
    '''
    Plots volume track plot
    '''
    fig = plt.figure(figsize = (5, 3));

    plt.plot(LR.time, LR['left flux (m^2/s)'], label= 'left')
    plt.plot(LR.time, LR['right flux (m^2/s)'], label= 'right')
    plt.plot(BT.time, BT['bottom flux (m^2/s)'], label = 'bottom')
    plt.plot(BT.time, BT['top flux (m^2/s)'], label = 'top')

    plt.legend()
    plt.xlabel('time')
    
    plt.tight_layout()
    
    return fig    



def colormap(array,ax = None):
    '''

    '''
    if not ax:
        plt.figure(figsize = (5,3))
        ax = plt.gca()
  
    g = ax.pcolormesh(np.fliplr(array), cmap = 'Blues')
    plt.colorbar(g, ax = ax)
    ax.set_xticks([])
    ax.set_yticks([])



def plot_Re(summary):
    '''
    Plot 
    '''
    b =np.linspace(0, np.max(summary.Re), 20)
    plt.hist(summary.Re, b, label = 'mean $Re(t_r)$', alpha = 0.5)
    plt.hist(summary.Re_all, b, label = 'mean $Re$', alpha = 0.5)
    plt.legend()
    plt.xlabel('$Re$')


############################  Tracer visuals ###########


def infiltrated_locations(sim, positions, 
            title = 'Locations of infiltrated particles', 
            dcut = 0, ucut = 0, 
            alpha = 0.1, stack = 0, scale = 3):
    '''
    Illustrating where the tracers go
    '''
    l = sim.l - dcut - ucut

    plt.figure(figsize = (scale*l/sim.L, scale))

    if stack == 0:
        plt.pcolormesh(sim.yc, sim.xc, sim.hc[sim.i_tr], cmap = 'Blues')

    elif stack == 1:
        xc = np.vstack((sim.xc - sim.dx/2, sim.xc + sim.dx/2))
        yc = np.vstack((sim.yc,sim.yc))
        hc = np.vstack((sim.hc[sim.i_tr], sim.hc[sim.i_tr]))
        plt.pcolormesh(yc, xc, hc, cmap = 'Blues')
        
    inds = []
    
    for i in range(len(positions)):
        p = positions.iloc[i]
        xo = np.array(p.xo)
        yo = np.array(p.yo)    
        
        #     if  p.escape_in_domain == 1:
        #         plt.plot( p.y_f,p.x_f,  'y.', alpha = 1)   

        if  p.t_f < positions.t_f.max() and p.y_f < sim.l - 1:
            plt.plot( p.y_f,p.x_f,  'y.', alpha = alpha)       

        if  p.escape==0:
            plt.plot( p.y_f,p.x_f,  'y.', alpha = alpha)               

    plt.xticks([])
    plt.yticks([])
    plt.xlim(ucut, sim.l-dcut)
    plt.ylim(0, sim.L)
    plt.title(title)


def plot_3D_trajectories(sim, positions, escape = '', ax = '', trapped = '', title = '', dcut = 1,  ucut = 70, min_curve = 0., point = 1, stack = 0, alpha = 0.1):
    '''
    '''
    if dcut == 0:
        dcut = 1
    if ax == '':
        fig = plt.figure(figsize= (6,5))
        ax =    fig.add_subplot(111, projection='3d');

    if title == 1:
        title =  r'$p$={0}, $\sigma_y$={1},  $F_v = {2}$'.format(
        sim.p, sim.sigma, sim.fV)

    ax = plot_surface(sim, 'veg', title, ucut = ucut, dcut = dcut, alpha = 0.5,vmax = 1.5, color = cm.Greens, plot_veg = 0, ax = ax, stack = stack);

    for ind in range(len(positions)):

        p = positions.iloc[ind]
        xo = np.array(p.xo)
        yo = np.array(p.yo)    
        xo = xo[yo <= sim.l-dcut*sim.dx] 
        yo = yo[yo <= sim.l-dcut*sim.dx] 
        
        if len(xo) > 0 and p.escape == 1:
            xo = np.concatenate((xo , [xo[-1]]))
            yo = np.concatenate((yo , [sim.l-dcut*sim.dx]))
        
        xo = xo[yo >= ucut*sim.dx]
        yo = yo[yo >= ucut*sim.dx]
        zo = np.array([sim.zc[int(xo[i]/sim.dx), min(int(yo[i]/sim.dx), sim.Nycell-1)] for i in range(len(xo))])
        
        if escape != '':

            if p.escape == 1 and np.std(xo) > min_curve:
                ax.plot(xo,  yo, zo, escape, alpha = alpha)        

                # if point == 1:
                #     ax.plot(xo[-1],  yo[-1], zo[-1], escape +  '.', ms= 5, alpha = 0.4)                          

        if trapped != '':

            if p.escape == 0 and np.std(xo) > min_curve:
                ax.plot(xo,  yo, zo, trapped, alpha = alpha)  ;    
                
                if point == 1:
                    ax.plot(xo[-1],  yo[-1], zo[-1], trapped +  '.', ms = 5, alpha = 0.4) ;                         

    # ax.set_xlim(ucut*sim.dx, sim.l-dcut*sim.dx)
    ax.view_init(25, 20);

    return ax



def plot_3D_trajectories_I(sim, positions, escape = '', 
                         trapped = '', title = '', dcut = 30, ucut = 70,
                         min_curve = 0., point = 1, stack = 0, alpha = 0.1):
    '''
    '''
    fig = plt.figure(figsize= (10, 4));
    ax =   fig.add_subplot(111, projection='3d');

    t = ax.set_title(title, fontsize = 10)
    t.set_y(0.9)
    _ = plot_surface(sim, 'infl_2d',  title, ucut = ucut, dcut = dcut, alpha = 0.8,
                 color = cm.Blues, plot_veg = 0, ax = ax, stack = stack);
    
    for ind in range(len(positions)):

        p = positions.iloc[ind]
        xo = np.array(p.xo)
        yo = np.array(p.yo)    
        xo = xo[yo <= sim.l-dcut*sim.dx]
        yo = yo[yo <= sim.l-dcut*sim.dx]
        
        xo = xo[yo >= ucut*sim.dx]
        yo = yo[yo >= ucut*sim.dx]
        zo = np.array([sim.zc[int(xo[i]/sim.dx), int(yo[i]/sim.dx)] for i in range(len(xo))])
        
        if escape != '':

            if p.escape == 1 and np.std(xo) > min_curve:
                ax.plot(xo,  yo, zo, escape, alpha = alpha)        

                if point == 1:
                    ax.plot(xo[-1],  yo[-1], zo[-1], escape +  '.', ms= 5, alpha = 0.4)                          

        if trapped != '':

            if p.escape == 0 and np.std(xo) > min_curve:
                ax.plot(xo,  yo, zo, trapped, alpha = alpha)      
                
                if point == 1:
                    ax.plot(xo[-1],  yo[-1], zo[-1], trapped +  '.', ms = 5, alpha = 0.4)                          
        
    cbaxes = fig.add_axes([0.65, 0.25, 0.015, 0.45]) 
    fmt = lambda x, pos: '{:.1f}'.format(x)
    format_colorbar( sim.infl_2d.max(), cm.Blues, cbaxes, fmt )
    
def plot_3D_points(sim, positions, escape = '', ax = '',
                         trapped = '', title = '', dcut = 0, ucut = 0,
                         min_curve = 0., point = 1, stack = 0, alpha = 0.05):
    '''
    '''
    fig = plt.figure(figsize= (10, 4));
    
    if ax == '':
        ax =    fig.add_subplot(111, projection='3d');


    if title == 1:
        title =  r'$p$={0}, $\sigma_y$={1},  $F_v = {2}$'.format(
        sim.p, sim.sigma, sim.fV)

    _ = plot_surface(sim, 'veg',  title, ucut = ucut, dcut = dcut,
                 color = cm.Greens, plot_veg = 0, ax = ax, stack = stack);

    for ind in range(len(positions)):

        p = positions.iloc[ind]
        xo = np.array(p.xo)
        yo = np.array(p.yo)    
        xo = xo[yo <= sim.l-dcut*sim.dx]
        yo = yo[yo <= sim.l-dcut*sim.dx]
        
        xo = xo[yo >= ucut*sim.dx]
        yo = yo[yo >= ucut*sim.dx]
        zo = np.array([sim.zc[int(xo[i]/sim.dx), int(yo[i]/sim.dx)] for i in range(len(xo))])
        
        if escape != '':

            if ind < 500 and p.escape == 1 and np.std(xo) > min_curve:
                ax.plot(xo,  yo, zo, escape, alpha = alpha)        

                if point == 1:
                    ax.plot(xo[0],  yo[0], zo[0], escape +  '.', ms= 5, alpha = 1)                          

        if trapped != '':

            if ind < 500 and p.escape == 0 and np.std(xo) > min_curve:
                ax.plot(xo,  yo, zo, trapped, alpha = alpha)      
                
                if point == 1:
                    ax.plot(xo[0],  yo[0], zo[0], trapped +  '.', ms = 5, alpha = 1)                          

    ax.view_init(25, 30);


def plot_trajectories(sim, positions, title= '', min_curve = 0, ucut = 0, dcut = 0):
    '''
    '''
    dx = sim.dx
    l = (sim.l - dcut*dx - ucut*dx)/sim.L
    plt.figure(figsize = (l*2., 2.))

    inds = []
    for i in range(len(positions)):
        
        p = positions.iloc[i]
        xo = np.array(p.xo)
        yo = np.array(p.yo)    
 

        if p.escape == 1 and np.std(xo) > min_curve:
            plt.plot(yo,  xo, 'b', alpha = 0.1)        

        if p.escape == 0 and np.std(xo) > min_curve:
            plt.plot(yo,  xo, 'y', alpha = 0.1)        

    plt.pcolormesh(sim.yc, sim.xc, sim.veg, cmap = 'Greens', alpha = 0.5)
    plt.title(title)
    plt.xticks([])
    plt.yticks([])

    plt.gca().set_xlim(ucut*sim.dx, sim.l - dcut*sim.dx - sim.dx)


def streamplot(sim, title = '', ucut = 0, dcut = 0):
    '''
    Varying color along a streamline    
    '''
    Re = sim.vc*sim.hc/1e-6
    l = sim.l/sim.dx - dcut - ucut
    plt.figure(figsize = (l/sim.L*sim.dx*2., 2.))
    ax = plt.gca()
    ax.pcolormesh(sim.yc[0], sim.xc[:, 0], sim.veg, cmap = 'Greens', vmin = 0, vmax = 1, alpha = 0.2)

    strm = ax.streamplot(sim.yc[0], sim.xc[:, 0], sim.vc[sim.i_tr], sim.uc[sim.i_tr], 
                         color = np.sqrt(1 + Re[sim.i_tr]), 
                         linewidth = 1, #+ Re[sim.i_tr]/Re.max(),
                         cmap = 'Blues', arrowstyle='-',
                         density = 1
                        )
    
    ax.set_title(title)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlim(ucut, sim.l-dcut*sim.dx)

def plot_Re(summary):
    '''
    '''
    b =np.linspace(0, np.max(summary.Re), 20)
    plt.hist(summary.Re, b, label = 'mean $Re(t_r)$', alpha = 0.5)
    plt.hist(summary.Re_all, b, label = 'mean $Re$', alpha = 0.5)
    plt.legend()
    plt.xlabel('$Re$')


#### Trouble shooting / assessment measures for tracers
from source_functions_1p3 import *

def plot_tracer_error(summary):
    '''
    Histogram of tracer errors
        plot_tracer_error(summary)
    '''
    bins = np.arange(-0.05, 0.055, 0.006)
    plt.hist( summary.infl_frac - summary.tracer_IF , 20, alpha = 0.5, label= 'rain')
    plt.xlabel('Tracer error: $C_L - C$')
    # plt.title('Tracer runoff coefficient error')
    plt.ylabel('Count')

    
def compare_hydrographs(sim):
    '''
    Compare real and tracer hydrographs 
    '''
    positions, recap = integrate_positions(sim, N = 1000, t_min = 1)
    a, b = plt.histogram(positions.query('escape == 1').t_esc, np.arange(0, (sim.tr+60)*60,120))
    plt.plot(sim['t'][:len(sim.hydro)] - 120, sim.hydro[:len(sim.t)]/sim.hydro.max())
    plt.plot(b[1:], a/a.max())


############# Visual helpers - formatting ################
names = { 'fV' : 'F_v',
            'dC_par' : 'C_{||} - C_{\perp}',
             'phi_veg' : '\\phi_{V}', 
             'sigma' : r'\\sigma_x',
             'excess' : r'(p - K_s)',
             'excess_tr' : r'(p - K_s) t_R',
             'aniso' : r'\\sigma_y/\\sigma_x',
             'Re_all': 'Re', 
             'offset' : '\\Delta x',
             # 'divert' : r'\\xi' ,
             'Ks_v'  : 'K_v',
             'K_v'  : 'K_b',
             'tr'  : 't_R',

             'alpha_b'  : r'\\alpha_b',
             'alpha_v'  : r'\\alpha_v',
             'Slope': 'S_o',
             'Precipitation': 'p',
             'p-Ks_v'  : 'p-K_v',
             'hydro' : 'Q_L',
             'final_qL' : 'Q_L(end)',
             'flashy' : 't_{rise}',

             'Delta_flashy' : '\\Delta t_{rise}',
             'curve' : '\\eta',
             'divert' : 'L_{esc} - L_{infl}',      
             'path_avg' : 'L',
             ## Stresses   
             'tau_v' : '\\tau_v',
             'tau_b' : '\\tau_b',
             'tau_b_std' : 'std(\\tau_b)',         
             'tau_edge' : '\\tau_{edge}',
             'tau_edge_a' : '\\tau_{approach}',
             ## Velocities
             'h_edge' : 'h_{edge}',
             'V_edge' : 'V_{edge}',
             'U_edge' : 'U_{edge}',
             't_HD'  : 't_{HD}',
             'U_te'  : 'U_{te}',
             't_te'  : 't_{te}',
             't_fall'  : 't_{fall}',
             't_t'  : 't_{t}',
             'U_avg'  : 'U_{avg}',
             
             'U_tr_mean' : 'mean(U)',
             'U_max' : 'U',
             'U_std' : 'std(U)',
             'V_std' : 'std(V)', 
             'V_tr_mean' : 'mean( V)',
             'std(V/U)' : 'std(V/U)',          
             'V_tr_std' : 'std(V_{t_R})',
             'U_tr_std' : 'std(U_{t_R})',

             'Fr_mean' : 'Fr',
             'Fr_max' : 'max(Fr)',

             ## scenarios
             'Delta_C' : '\\Delta C',
             'Delta_IF' : '\\Delta IF',
             'Delta_I_v' : '\\Delta I_v',
             'Delta_I_v_tr' : '\\Delta I_v(rise)',
             'Delta_I_v_rec' : '\\Delta I_v(rec)',
             'Delta_U' : '\\Delta U_{max}',
             'Delta_U_edge' : '\\Delta U_{edge}',
             
             'Delta_V' : '\\Delta V_{max}  ',
             'Delta_V_edge' : '\\Delta V_{edge} ',
             'Delta_h' : '\\Delta h_{max}',
             'Delta_h_edge' : '\\Delta h_{edge}',
             'Delta_curve' : '\\Delta \\eta',
             'Delta_V_std' : '\\Delta std(V)', 
             'Delta_tau_b' : '\\Delta \\tau_b',
             'Delta_tau_edge' : '\\Delta \\tau_{edge}',

             'path_IF' : 'L_{infl}',
             'path_source' : 'L_{esc}',
             'infl_frac' : 'IF',
             # 'I_v' : 'run–on \\ fraction, I_v',  
             'IF_v' : 'IF_v',
             'I_r' : 'I_{cc}/I_{v}',
             'I_v_tr' : 'I_v(rise)',
             'I_v_rec' : 'I_v(rec)',
             'positive' : 'I_v/I_b',

             ## tracers
             'tracer_err' : 'IF - IF_L',
             'abs_err' : '|IF - IF_L|',
             'std_LS' : 'std(L_{esc})',
             'std_Linfl' : 'std(L_{infl})',
             
             ## flume
             'p_equiv' : 'p_{equiv}',

             ## run around
             'std_flowlines' : 'std(FL)',
             'diverted_mean' : 'diverted',
             'Delta_diverted_mean' : '\\Delta diverted',
             'captured_mean' : 'captured',
             'Delta_captured_mean' : '\\Delta captured',

             'dthetacoef' : '\\Delta \\theta',

             'Litter_total' : 'Litter',
             'Velocity_mean' : 'U',

            'K_invRe' : 'K',
            'K_invf' : 'K',
            'Foliar_total' : 'F_{foliar}',
            'Litter_total' : 'F_{litter}',
            'Rock_total' : 'F_{rock}',
            'Basal_total' : 'F_{basal}',
            'Ground_cover_total' : 'F_{ground-cover}',
            'Soil_total' : 'F_{soil}',
            'Canopy_gap_total' : 'F_{canopy-gap}',
            'Canopy_gap_average' : 'F_{canopy-gap-avg}',

            'Basal_gap_total' : 'F_{basal-gap}',
            'Basal_gap_average' : 'F_{basal-gap-avg}',

            'Foliar_forb' : 'F_{foliar-forb}',
            'Foliar_grass' : 'F_{foliar-grass}',
            'Foliar_shrub' : 'F_{foliar-shrub}',
            'Foliar_total' : 'F_{foliar}',
            'Litter_protected' : 'F_{litter-protected}',
            'Litter_unprotected' : 'F_{litter-unprotected}',

            'Rock_protected' : 'F_{rock-protected}',
            'Rock_unprotected' :'F_{rock_unprotected}' , 
            'Rock_umprotected' :'F_{rock_unprotected}' , 

            'Basal_protected' : 'F_{basal-protected}',
            'Basal_unprotected' : 'F_{basal-unprotected}',            

            'Soil_protected' : 'F_{soil-protected}',
            'Soil_unprotected' : 'F_{soil-unprotected}',

            'Ground_cover_protected' : 'F_{ground-cover-protected}',
            'Ground_cover_unprotected' : 'F_{ground-cover-unprotected}',

             'm1' : 'T',
             'f_dw' : 'f_{dw}',
             'dw' : 'f_{E}',         
             'nrmse_dw' : 'NRMSE_{dw}',
             'nrmse_invf' : 'NRMSE_{K}',
             'nrmse_invRe' : 'NRMSE_{K}',
             'r2_dw' : 'r^2_{dw}',
             'r2_invf' : 'r^2_{K}',
             'r2_invRe' : 'r^2_{K}',         
}


def format_name(fld, updates = {}):
    names.update(updates)
    fld = fld.replace("_reg", "_{reg}")
    fld = fld.replace("_avg", "_{avg}").replace("_cal", "_c")
    if fld in names:
        return  '${0}$'.format(names[fld]).replace("\\\\", "\\").replace("D_", "\Delta \ ")
    else: 
        return '${0}$'.format(fld).replace("\\\\", "\\").replace("D_", "\Delta ")




def reformat_to_numeric(summary, fld):
    '''
    For colorbars 
    '''
    count = 1
    for val in summary[fld].unique():
        
        inds = summary[summary[fld] == val].index
        summary.loc[inds, 'case'] = summary.loc[inds, fld].replace(val, str(count))

        count += 1
    return summary

def make_title(sim, cols = ['fV', 'p', 'Ks_v'], prec = 2 ):
    '''
    '''
    try:
        sim[cols] = sim[cols].astype(float)
    except:
        pass
    title = ', '.join(['{0}={1}'.format(names[k], np.round(sim[k],prec)) if k in 
        names else '{0}={1}'.format(k, np.round(sim[k], prec)) for k in cols])
    title = '${0}$'.format(title)  
    
    return title.replace('\\\\', '\\')


def plot_directly_connected(sim):

    L_dc = np.zeros_like(sim['patch_LB'])
    for i in range(sim.Nxcell):
        l = int(sim.patch_LB[i, -1])
        L_dc[i, -l:] = l

    sim['L_dc_array']  = L_dc

    fig = plt.figure(figsize= (10, 4))

    ax =    fig.add_subplot(121, projection='3d');
    _ = plot_surface(sim,'L_dc_array', 'Directly connected $L_{dc}$ ',
                 color = cm.BrBG_r, plot_veg = 0, ax = ax,  ucut = 0, dcut = 0)

    ax =    fig.add_subplot(122, projection='3d');
    _ = plot_surface(sim,'patch_LB', 'Flowlength $L_{b}$',
                 color =  cm.BrBG_r, plot_veg = 0, ax = ax,  ucut = 0, dcut = 0)


def format_plot(summary, xfld, yfld, cfld, cmap = cm.ocean, 
                alpha = 0.5, err_tol = 0.1, vmax = None, vmin = None,
                s = None, edgecolor = None):
    '''
    Will create a color-coded scatter plot for tracer fields
    continuous colorbar!
    '''
    if err_tol:
        inds = summary.query('abs_err < {0}'.format(err_tol)).index
    else: 
        inds = summary.index
    c = plt.scatter(summary.loc[inds][xfld].astype(float), 
                summary.loc[inds][yfld].astype(float), 
                c = summary.loc[inds][cfld].astype(float),  cmap = cmap, visible = 0, vmin = vmin, vmax = vmax)

    plt.scatter(summary.loc[inds][xfld].astype(float),  
                summary.loc[inds][yfld].astype(float),  
                c = summary.loc[inds][cfld].astype(float),  cmap = cmap, alpha =alpha,
                vmin = vmin, edgecolor = edgecolor,
                 vmax = vmax, s = s)

    plt.xlabel(format_name(xfld))
    plt.ylabel(format_name(yfld))
    plt.colorbar(c ,label = format_name(cfld))

def format_plot_sns(subset, xfld, yfld, cfld, ax = None,  vmax = 80, vmin = 0, 
                legend = 0, alpha = 1, legend_title = None, cmap = cm.ocean):
    """
    updated: includes discrete legend
    """
    if ax is None:
        ax = plt.gca()
        
    sns.scatterplot(data = subset.loc[inds],
                    x = xfld, y = yfld,
                    hue = cfld, palette = cmap, ax = ax, s = 40,
                    alpha = alpha, vmin = vmin, vmax = vmax)

    ax.set_xlabel(format_name(xfld))
    ax.set_ylabel(format_name(yfld))
    
    if legend:
        ax.legend(title = legend_title)
    else: 
        ax.legend([])
        ax.get_legend().remove()
        

###################################
# Equations to fit a power law....

import statsmodels.api as sm
import statsmodels.formula.api as smf


def fit_power_law(subset, target = 'curve', 
                cols = ['fV', 'phi_veg', 'p', ], train_test = 0):
    '''
    Fit equation:
       log (y) = a log(x_1) + b log(x_2)   
    
    '''
    from sklearn.model_selection import train_test_split

    data = np.log(subset[[target]].astype(float))
    
    for col in cols:
        data[col] = np.log(subset[col].astype(float))

    if train_test == 1:
        
        X_train, X_test, y_train, y_test = train_test_split(
            data, data[target], test_size=0.3, random_state=1)

        train = pd.concat([X_train, y_train])
    
        res = smf.ols(formula= '{0} ~ {1}'.format(target, ' + '.join(cols)), data=X_train).fit()
    
    else:
        res = smf.ols(formula= '{0} ~ {1}'.format(target, '+'.join(cols)), data=data).fit()

    data['y_hat'] = res.predict(data[cols]) 

    return res

######## Power law fit code... #######

def make_equation(res):
    '''
    returns a power law equation    
    '''

    inds = abs(res.params.round(2)) > 5e-2
    res.params = res.params[inds]
    
    result = list(zip( res.params.index, res.params))    
    #equation = '1/{0:.0f} \ '.format(1/np.exp(result[0][1])) + \
    if 'Intercept' in res.params:
        equation =  '*'.join(['{0}^{1:.2f}'.format(r[0],r[1]) for r in result[1:]])
        equation = '{0:.2f} \ '.format(np.exp(res.params.Intercept)) + equation
    else:
        equation =  '*'.join(['{0}^{1:.2f}'.format(r[0],r[1]) for r in result])
    equation = equation.replace("_reg", "_{reg}")
    equation = equation.replace("_avg", "_{avg}")


    local_names = { 
              'Re_all': 'Re',
              'Slope': 'So',
              'Precipitation': 'p',
              'path_avg' : 'L',
              'dw_guess' : 'f_{guess}',
             'dw' : 'f',
              'p_Ks': '(p-K_v)'}

    for s in res.params.index:
        
        if s in local_names:

            equation = equation.replace(s, local_names[s])
            
        elif s in names:

            equation = equation.replace(s, names[s])

    equation = equation.replace('^', '^{').replace('*', '} ') + '}'    
    return r'${0}$'.format(equation).replace('\\\\', '\\')



def make_prediction(res, cols, subset):
    ''' 
    returns prediction for a power law fit
    '''
    x = np.log(subset[cols].astype(float))
    logy = res.predict(x)
    y = np.exp(logy)
    return y

def wrap_fit(summary, cols = ['D',  'excess', 'fV'], target = 'path_avg'):
    '''
    '''
    subset = summary.copy()
    
    
    res = fit_power_law(subset, target = target, cols =  cols)
    resid = make_prediction(res, cols, subset ) - subset[target]

    equation = make_equation(res)
    subset[target + '_hat'] = make_prediction(res, cols, subset)
    # rint (res.rsquared_adj, np.mean(np.abs(resid**2)))
    
    return res, subset, equation


############ Linear fit code ###########

###################################
# Equations to fit a linear law....

import statsmodels.api as sm
import statsmodels.formula.api as smf


def fit_linear(subset, target = 'curve', 
                cols = ['fV', 'phi_veg', 'p', ], train_test = 0):
    '''
    Fit equation:
       log (y) = a log(x_1) + b log(x_2)   
    
    '''
    from sklearn.model_selection import train_test_split
 
    
    data = (subset[[target]+ cols].astype(float))
    

    if train_test == 1:
        
        X_train, X_test, y_train, y_test = train_test_split(
            data, data[target], test_size=0.3, random_state=1)

        train = pd.concat([X_train, y_train])
    
        res = smf.ols(formula= '{0} ~ {1}'.format(target, ' + '.join(cols)), data=X_train).fit()
    
    else:
        res = smf.ols(formula= '{0} ~ 1 + {1}'.format(target, '+'.join(cols)), data=data).fit()

    data['y_hat'] = res.predict(data[cols]) 

    return res


def make_linear_equation(res):
    '''
    returns a linear law equation    
    '''

    inds = abs(res.params.round(2)) > 1e-2
    res.params = res.params[inds]
    
    result = list(zip( res.params.index, res.params))    
    #equation = '1/{0:.0f} \ '.format(1/np.exp(result[0][1])) + \
    equation =  '+'.join(['{1:.1f}\ {0}'.format(r[0],r[1]) for r in result[1:]])

    local_names = { 
              'Re_all': 'Re',
              'Slope': 'So',
              'Precipitation': 'p',
              'path_avg' : 'L',
              'amplitude' : 'A',
              'p_Ks': '(p-K_v)'}

    for s in res.params.index:
        
        if s in local_names:

            equation = equation.replace(s, local_names[s])
            
        elif s in names:

            equation = equation.replace(s, names[s])
            
            equation = equation.replace('D_', '\Delta ')            

    #equation = equation.replace('^', '^{').replace('*', '} ') + '}'    
    equation = equation.replace('+-', '-')
    return r'${0}$'.format(equation).replace('\\\\', '\\')


def make_linear_prediction(res, cols, subset):
    ''' 
    returns prediction for a linear law fit
    '''
    x = (subset[cols].astype(float))
    y = res.predict(x)

    return y

def wrap_linear_fit(summary, cols = ['D',  'excess', 'fV'], target = 'path_avg'):
    '''
    '''
    subset = summary.copy()
    
    
    res =  fit_linear(subset, target = target, cols =  cols)
    resid = make_linear_prediction(res, cols, subset ) - subset[target]

    equation = make_linear_equation(res)
    
    return res, subset, equation


#### GRL fits ########## ##### ##### 

from scipy import stats
from scipy.optimize import curve_fit

def fit_exponent_grl(subset):
    
    x = subset.path_source.astype(float)
    y = 1 - subset.infl_frac.astype(float)

    pars, cov = curve_fit(f=exponent, xdata=x, ydata=y, p0=[0, 0], bounds=(-np.inf, np.inf))
    
    # Get the standard deviations of the parameters (square roots of the # diagonal of the covariance)
    stdevs = np.sqrt(np.diag(cov))
    # Calculate the residuals
    res = y - exponent(x, *pars)

    a = pars[0]
    b = pars[1]
    rmse = np.sqrt(np.mean(res**2))

    print('y = b + e^(a x)')
    coefs = ['a', 'b']
    for i, p, var in zip(range(n), pars, np.diag(cov)):
        sigma = var**0.5
        print ('{0}: {1:.2f} [{2:.2f} ]'.format(coefs[i], p,
                                      sigma*tval))


    stdevs = np.sqrt(np.diag(cov))
    return x, y, a, b, rmse, stdevs,res

def exponent(x, a, b):
    return b*np.exp(x*a)



def linear(x, a, b):
    
    return b + a*x

def fit_linear_grl(subset):
    
    x =  np.log(subset.path_source.astype(float))
    y = np.log((1 - subset.infl_frac.astype(float)))

    pars, cov = curve_fit(f=linear, xdata=x, ydata=y, p0=[0, 0], bounds=(-np.inf, np.inf))
    # Get the standard deviations of the parameters (square roots of the # diagonal of the covariance)
    stdevs = np.sqrt(np.diag(cov))
    # Calculate the residuals
    res = np.exp(y) - np.exp(linear(x, *pars))

    kappa = pars[0]
    c = pars[1]
    rmse = np.sqrt(np.mean(res**2))

    rmse = np.sqrt(np.mean(res**2))

    alpha = 0.05 # 95% confidence interval = 100*(1-alpha)
    n = len(y)    # number of data points
    p = len(pars) # number of parameters
    
    dof = max(0, n - p) # number of degrees of freedom
    tval = stats.t.ppf(1.0-alpha/2., dof) 

    print('log(y) = k log(x) + c')
    coefs = ['k', 'c']
    for i, p, var in zip(range(n), pars, np.diag(cov)):
        sigma = var**0.5
        print ('{0}: {1:.2f} [{2:.2f} ]'.format(coefs[i], p,
                                      sigma*tval))

    stdevs = np.sqrt(np.diag(cov))
    return x, y, kappa, c, rmse, stdevs, res    



# https://kitchingroup.cheme.cmu.edu/blog/2013/02/12/Nonlinear-curve-fitting-with-parameter-confidence-intervals/
def fit_power_grl(subset):
    
    df_grouped_by = subset.groupby(['path_source','C'])
    df_balanced = df_grouped_by.apply(lambda x: x.sample(df_grouped_by.size().min()).reset_index(drop=True))
    df_balanced = df_balanced.droplevel(['path_source', 'C'])
    subset  =  df_balanced
    x =    subset.path_source.astype(float)
    y = 1 - subset.infl_frac.astype(float)

    pars, cov = curve_fit(f=power_law, xdata=x, ydata=y, p0=[0, 0], bounds=(-np.inf, np.inf))
    # Get the standard deviations of the parameters (square roots of the # diagonal of the covariance)
    stdevs = np.sqrt(np.diag(cov))
    # Calculate the residuals
    res = y - power_law(x, *pars)

    kappa = pars[0]
    c = pars[1]
    rmse = np.sqrt(np.mean(res**2))

    alpha = 0.05 # 95% confidence interval = 100*(1-alpha)
    n = len(y)    # number of data points
    p = len(pars) # number of parameters
    
    dof = max(0, n - p) # number of degrees of freedom
    tval = stats.t.ppf(1.0-alpha/2., dof) 


    print('y = c x^k')
    coefs = ['k', 'c']
    for i, p, var in zip(range(n), pars, np.diag(cov)):
        sigma = var**0.5
        print ('{0}: {1:.2f} [{2:.2f}]'.format(coefs[i], p,
                                      sigma*tval))
                
    stdevs = np.sqrt(np.diag(cov))        
    for i, p,var in zip(range(n), pars, np.diag(cov)):
     # computed as 1.96 * standard deviation error
        CI = stdevs[i]*1.96
        #print ('p{0}: {1:.2f} [{2:.2f} ]'.format(i, p, CI))

    return x, y, kappa, c, rmse, stdevs,res

def power_law(x, a, b):
    return b*np.power(x, a)    