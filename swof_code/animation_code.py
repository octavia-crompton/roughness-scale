import numpy as np  # Numerical operations
import matplotlib.pyplot as plt  # Plotting
import cv2  # OpenCV for image processing
plt.rcParams['savefig.bbox'] = 'tight'  # Savefig config
import matplotlib.animation as animation  # Animation support
import matplotlib.colors as colors  # Color normalization
from matplotlib import cm  # Colormap handling
import matplotlib.gridspec as gridspec  # Grid layout for subplots
import cmocean  # Oceanographic colormaps
from plot_SWOF import *  # Project-specific plotting utilities

def fix_axes(ax):
    """
    Configure 3D axes for cleaner visuals (remove colored panes, ticks, grid).
    """
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


def make_filename(sim, cols = [ 'p', 'Ks_v']):
    """
    Generate a filename string from simulation parameters.
    """
    title = ','.join(['{0}={1}'.format(k, sim[k]) for k in cols])
    title = '{0}'.format(title)  
    if "scenario" in  sim:
        return  title  +  "," + sim.scenario  
    else: 
        return title

def animation_flds(sim,  dt_anim = 10, trim = 100, truncate = True, max_U = 0.2):
    """
    Extract and preprocess simulation fields for animation.
    Returns: processed arrays and animation parameters.
    """
    fps = int(60/sim['dt']) # frames per second
    freq = max(1, int(dt_anim/sim['dt']))  # animation frame frequency
    hydro_freq =  max(1, int(dt_anim))     # hydrograph frequency
    U = np.sqrt(sim.uc**2 + sim.vc**2)[::freq, :, :trim]  # velocity magnitude
    V = sim.uc[::freq, :, :trim]  # u-component
    I = sim.I[::freq, :, :trim]   # infiltration
    # Patch edge values for U
    t, x, y = np.where(U > max_U)
    xx = x[x< sim.Nxcell]
    U[t, xx , y] =   U[t, xx - 1, y] 
    xx = x[x ==  sim.Nxcell]
    if len(xx > 0):
        U[t, xx , y] =  U[t, xx - 1, y] 

    hc =  sim.hc[::freq, :, :trim]
    t = sim.t[::freq]
    hydro = sim.hydro[::freq][:len(t)]
    
    if len(hydro) < len(t):
        hydro = np.insert(hydro, -1, 0)

    yc = sim.yc[:, :trim]
    xc = sim.xc[:, :trim] 
    zc = sim.zc[:, :trim]
    
    # remove microtopography
    
    zc = zc - zc.min() 
    
    if truncate:
        frn = (hc.mean(1).mean(1)*100 > 0.0001).sum()
    else:
        frn = len(hc) - 1 # (hc.mean(1).mean(1)*100 > 0).sum()

    return fps, freq, hydro_freq, U,V, I, hc, t, hydro, xc, yc, zc,frn


def format_colorbar(fld, cmap, cbaxes, fmt, label = "cm"  , vmin = 0):
    """
    Draw a colorbar for a given field and colormap.
    """
    from matplotlib.ticker import FuncFormatter
    a = np.array([[np.round(vmin*100, 2), np.round(fld*100, 2)]])
    plt.figure(figsize=(0, 0))
    img = plt.imshow(a, cmap=cmap)
    plt.gca().set_visible(False)
    cbar = plt.colorbar(img, cax = cbaxes, format=FuncFormatter(fmt), label = label)
    


class animation_class:
    """
    Container for FullSWOF simulation results and derived fields for animation.
    """
    def __init__(self, sim, dt_anim = 10, trim = 100):
        # Animation parameters
        fps = int(sim['dt']/3) # frame per sec
        freq = max(1, int(dt_anim/sim['dt']))
        hydro_freq =  max(1, int(dt_anim))
        # Extract fields
        U = np.sqrt(sim.uc**2 + sim.vc**2)[::freq, :, :trim]
        V = sim.uc[::freq, :, :trim]
        I = sim.I[::freq, :, :trim]    
        hc =  sim.hc[::freq, :, :trim]
        t = sim.t[::freq]
        hydro = sim.hydro[::freq][:len(U)]
        yc = sim.yc[:, :trim]
        xc = sim.xc[:, :trim] 
        zc = sim.zc[:, :trim]
        zc = zc - zc.min() 
        # Compute frame number
        frn = (hc.mean(1).mean(1)*100 > 0.01).sum()
        # Store fields
        self.fps = fps
        self.freq = freq
        self.hydro_freq = hydro_freq
        self.U = U
        self.V = V        
        self.I = I
        self.t = t
        self.hc = hc        
        self.hydro = hydro
        self.xc = xc        
        self.yc = yc
        self.zc = zc
        self.frn = frn 
        # Vegetation mask and colors
        veg_norm = colors.Normalize(0., 1.1)
        self.veg = sim.veg - cv2.erode(sim.veg, np.ones((3,3)), 1)
        self.veg_colors = cm.Greens(veg_norm(self.veg ))   


############ Static visuals... ################

def rescale (zed):
    """ 
    rescales an input array from 0 to 1
    """ 
    return (zed- zed.min())/(zed.max() - zed.min())


def static_hUq(sim, infl = 0, h_scale = 1, 
               cols =  ['p', 'Ks_v','fV', 'tr' ],  dt_anim = 60, 
               trim = 100, truncate = 0):
    """
    """
    def hydro_axis(ax, t, hydroh):
        """
        Set up axis for hydrograph
        """
        t = t/60.

        q = hydro
        ax.set_xlim(0, t[frn-1])
        ax.set_ylim(-0.01, np.max(q)*1.1)
        ax.set_xlabel("time (min)")
        ax.set_ylabel(" cm/hr")
        ax.set_title("Hydrograph", y = 1.38)

        return ax

    
    fps, freq, hydro_freq, U, V, I, hc, t, hydro, xc, yc, zc, frn = animation_flds(
        sim, dt_anim, trim, truncate)

    
    # remove microtopography    
    zc = rescale(1-yc)*sim.So*sim.l    
    Umax = np.percentile(U[np.isnan(U) == 0], 99.)*1.2
    hmax = hc.ravel().max()*1.2

    ### Set up figure
    fig = plt.figure(figsize = (17, 5))
    plt.subplots_adjust(wspace = -0.1)

    gs = gridspec.GridSpec(4, 9)
    ax = fig.add_subplot(gs[:, 0:3], projection='3d')
    ax = fix_axes(ax)

    ax2 = fig.add_subplot(gs[:, 3:6], projection='3d')
    ax2 = fix_axes(ax2)

    ax3 = fig.add_subplot(gs[1:3, 7:])
    ax3 = hydro_axis(ax3, t, hydro)

    ax.set_title("Depth", y =1)
    ax2.set_title("Velocity", y = 1)    

    ### Add vegetation
    veg_norm = colors.Normalize(0., 1.1)

    veg_contour = sim.veg - cv2.erode(sim.veg, np.ones((3,3)), 1)
    veg_colors = cm.Greens(veg_norm(veg_contour ))    
    
    contour = 0

    if contour == 1:
        ax.plot_surface(xc, yc, zc, facecolors = veg_colors , 
                            rstride = 1, cstride = 1, alpha = 0.8,
                            linewidth=0, antialiased=True, shade=False)      

        ax2.plot_surface(xc, yc, zc, facecolors = veg_colors , 
                            rstride = 1, cstride = 1, alpha = 0.8,
                            linewidth=0, antialiased=True, shade=False)      
        
    if infl == 0:
        veg = sim.veg
        veg_colors = cm.Greens(veg_norm(veg))    

        ax.plot_surface(xc, yc, zc*0, facecolors = veg_colors , 
                                rstride = 1, cstride = 1, alpha = 0.5,
                                linewidth=0, antialiased=True, shade=False)      

        ax2.plot_surface(xc, yc, zc*0, facecolors = veg_colors , 
                                rstride = 1, cstride = 1, alpha = 0.5,
                                linewidth=0, antialiased=True, shade=False)      
        

    i_tr = np.where(t < sim.tr*60)[-1][-1]
    h_norm = colors.Normalize(vmin= .0,  vmax = hmax)
    h_colors = cm.Blues(h_norm(hc[i_tr]))

    U_norm = colors.Normalize(vmin= -0.01,  vmax = Umax)  
    U_colors = cm.cividis(U_norm(U[i_tr])) 

    I_norm = colors.Normalize(vmin= -0.0,  vmax = I.max())  
    I_colors = GnBu(I_norm(I[-1])) 
    
    if infl:
        ax.plot_surface(xc, yc, zc*0 , facecolors = I_colors, 
                            rstride = 1, cstride = 1, alpha = 1,
                            linewidth=0, antialiased=True, shade=False)      

        ax2.plot_surface(xc, yc, zc*0, facecolors = I_colors, 
                            rstride = 1, cstride = 1, alpha = 1,
                            linewidth=0, antialiased=True, shade=False)          


    cbaxes = fig.add_axes([.355, 0.25, 0.01, 0.45]) 
    fmt = lambda x, pos: '{:.1f}'.format(x)
    format_colorbar( hc.max(), "Blues", cbaxes, fmt , label = "cm" )

    cbaxes = fig.add_axes([0.61, 0.25, 0.01, 0.45])  
    fmt = lambda x, pos: '{:.0f}'.format(x)
    format_colorbar( Umax, cm.cividis, cbaxes, fmt, label = "cm/s" )

    if infl: 
        cbaxes = fig.add_axes([0.10, 0.25, 0.01, 0.45])  
        fmt = lambda x, pos: '{:.0f}'.format(x)
        format_colorbar( I.max(), GnBu, cbaxes, fmt, label = "$I$ (cm)" )

    plot =  [ax.plot_surface(xc, yc, zc + hc[i_tr]*h_scale, facecolors = h_colors, 
                        rstride=1, cstride=1,cmap = GnBu, alpha = 0.8,
                        linewidth=0,antialiased=True, shade=False),
             ax2.plot_surface(xc, yc,  zc + hc[i_tr]*h_scale, facecolors = U_colors, 
                        rstride=1, cstride=1, alpha = 0.8,
                        linewidth=0,antialiased=True, shade=False),
             ax3.plot(t/60., hydro,
                           '#2b8cbe', animated = True)   
            ]

    ax2.view_init(25, 20)
    ax.view_init(25, 20)
    ax.set_zlim(0, (zc + hc*h_scale*2).max())
    ax2.set_zlim(0, (zc + hc*h_scale*2).max())

    
    plt.suptitle(make_title(sim,cols),fontsize = 12)


def static_hIq(sim, infl = 0, h_scale = 1, 
               cols =  ['p', 'Ks_v','fV', 'tr' ],  
               dt_anim = 60, trim = 100, truncate = 0):
    """
    """
    def hydro_axis(ax, t, hydroh):
        """
        Set up axis for hydrograph
        """
        t = t/60.

        q = hydro
        ax.set_xlim(0, t[frn])
        ax.set_ylim(-0.01, np.max(q)*1.1)
        ax.set_xlabel("time (min)")
        ax.set_ylabel(" cm/hr")
        ax.set_title("Hydrograph", y = 1.38)

        return ax

    fps, freq, hydro_freq, U, V, I, hc, t, hydro, xc, yc, zc, frn = animation_flds(
        sim, dt_anim, trim, truncate)
    veg = sim.veg[:, :trim]
    
    # remove microtopography 
    planar  = rescale(1-yc)*sim.So*sim.l    
    micro = zc - planar    
    zc = planar
    Umax = np.percentile(U[np.isnan(U) == 0], 99.)*1.2
    hmax = hc.ravel().max()*1.2

    ### Set up figure
    fig = plt.figure(figsize = (17, 5))
    plt.subplots_adjust(wspace = -0.1)

    gs = gridspec.GridSpec(4, 9)
    ax = fig.add_subplot(gs[:, 0:3], projection='3d')
    ax = fix_axes(ax)

    ax2 = fig.add_subplot(gs[:, 3:6], projection='3d')
    ax2 = fix_axes(ax2)

    ax3 = fig.add_subplot(gs[1:3, 7:])
    ax3 = hydro_axis(ax3, t, hydro)

    ax.set_title("Depth", y =1)
    ax2.set_title("Infiltration", y = 1)    

    ### Add vegetation
    veg_norm = colors.Normalize(0., 1.1)

    veg_contour = veg - cv2.erode(veg, np.ones((3,3)), 1)
    veg_colors = cm.Greens(veg_norm(veg_contour ))    
    
    contour = 0
    
    if contour == 1:
        ax.plot_surface(xc, yc, zc, facecolors = veg_colors , 
                            rstride = 1, cstride = 1, alpha = 0.8,
                            linewidth=0, antialiased=True, shade=False)      

        ax2.plot_surface(xc, yc, zc, facecolors = veg_colors , 
                            rstride = 1, cstride = 1, alpha = 0.8,
                            linewidth=0, antialiased=True, shade=False)      
        
        
    project_veg = 0
    if project_veg == 1:
        # for projection to horizontal
        veg_norm = colors.Normalize(0., 1.1)
        veg_colors = cm.Greens(veg_norm(veg ))    

        ax.plot_surface(xc, yc, zc*0, facecolors = veg_colors , 
                            rstride = 1, cstride = 1, alpha = 0.4,
                            linewidth=0, antialiased=True, shade=False)        

    elif project_veg == 2:
        # for incline
        ax.plot_surface(xc, yc, zc, facecolors = veg_colors , 
                                rstride = 1, cstride = 1, alpha = 1,
                                linewidth=0, antialiased=True, shade=False)  
    project_micro = 1
    
    if project_micro == 1:    
       
        micro_norm = colors.Normalize(- micro.max(), micro.max())
        micro_colors = cm.terrain(micro_norm(micro ))    

        ax.plot_surface(xc, yc, micro/5, facecolors = micro_colors , 
                            rstride = 1, cstride = 1, alpha = 0.4,
                            linewidth=0, antialiased=True, shade=False) 
        
        ax2.plot_surface(xc, yc, micro/5, facecolors = micro_colors , 
                            rstride = 1, cstride = 1, alpha = 0.4,
                            linewidth=0, antialiased=True, shade=False)         

        cbaxes = fig.add_axes([0.10, 0.25, 0.01, 0.45])  
        fmt = lambda x, pos: '{:.0f}'.format(x)
        format_colorbar( micro.max(),"terrain", cbaxes, fmt, label = "terrain (cm)" , vmin = 0)

        
    i_tr = np.where(t < sim.tr*60)[-1][-1]
    h_norm = colors.Normalize(vmin= .0,  vmax = hmax)
    h_colors = cm.Blues(h_norm(hc[i_tr]))

    I_norm = colors.Normalize(vmin= -0.0,  vmax = I.max())  
    I_colors = GnBu(I_norm(I[-1])) 


    cbaxes = fig.add_axes([.355, 0.25, 0.01, 0.45]) 
    fmt = lambda x, pos: '{:.1f}'.format(x)
    format_colorbar( hc.max(), "Blues", cbaxes, fmt , label = "$h$ (cm)" )

    cbaxes = fig.add_axes([0.61, 0.25, 0.01, 0.45])  
    fmt = lambda x, pos: '{:.0f}'.format(x)
    format_colorbar( I.max()*100, GnBu, cbaxes, fmt, label = "$I$ (cm)" )


    plot =  [ax.plot_surface(xc, yc, zc + hc[i_tr]*h_scale, facecolors = h_colors, 
                        rstride=1, cstride=1,cmap = GnBu, alpha = 0.8,
                        linewidth=0,antialiased=True, shade=False),
             ax2.plot_surface(xc, yc,  zc + hc[i_tr]*h_scale, facecolors = I_colors, 
                        rstride=1, cstride=1, alpha = 0.8,
                        linewidth=0,antialiased=True, shade=False),
             ax3.plot(t/60., hydro,
                           '#2b8cbe', animated = True)   
            ]

    ax2.view_init(25, 20)
    ax.view_init(25, 20)
    ax.set_zlim(0, (zc + hc*h_scale*2).max()*1.3)
    ax2.set_zlim(0, (zc + hc*h_scale*2).max()*1.3)

    
    plt.suptitle(make_title(sim,cols),fontsize = 12)
    
    

##### Code for comparing scenarios


class animation_class:
    """
    Contains the FullSWOF results
    """
    def __init__(self, sim, dt_anim = 10, trim = 151, scale_micro = True ):
        
        fps = int(sim['dt']/3) # frame per sec

        freq = max(1, int(dt_anim/sim['dt']))

        hydro_freq =  max(1, int(dt_anim))
        U = np.sqrt(sim.uc**2 + sim.vc**2)[::freq, :, :trim]
        V = sim.uc[::freq, :, :trim]
        I = sim.I[::freq, :, :trim]    

        hc =  sim.hc[::freq, :, :trim]
        t = sim.t[::freq]
        hydro = sim.hydro[::freq][:len(U)]

        yc = sim.yc[:, :trim]
        xc = sim.xc[:, :trim] 
        zc = sim.zc[:, :trim]
        zc = zc - zc.min() 

        if scale_micro == 1:
            
            planar  = rescale(1-yc)*sim.So*sim.l    
            micro = zc - planar
            zc = planar + micro/10

            self.micro = micro
            self.planar = planar
    

        frn = (hc.mean(1).mean(1)*100 > 0.01).sum()
        
        self.fps = fps
        self.freq = freq
        self.hydro_freq = hydro_freq
        self.U = U
        self.V = V        
        self.I = I
        self.t = t
        self.hc = hc        
        self.hydro = hydro
        self.xc = xc        
        self.yc = yc
        self.zc = zc
        self.frn = frn 
        

        veg_norm = colors.Normalize(0., 1.1)
        self.veg = sim.veg - cv2.erode(sim.veg, np.ones((3,3)), 1)
        self.veg_colors = cm.Greens(veg_norm(self.veg ))   
    

    
def static_visual(sim1, sim2, sim3, plot_infl = 1,
                  titles = ["With mound", "Without mound", "Difference"]):
    """
    """
    h_scale = 3

    a1 = animation_class(sim1, dt_anim = 30)
    a2 = animation_class(sim2, dt_anim = 30)
    a3 = animation_class(sim3, dt_anim = 30)

    
    frn = np.max([a2.frn, a1.frn, a3.frn])
    fps = a3.fps

    fig = plt.figure(figsize = (16, 5))
    plt.subplots_adjust(wspace = -0.1)
    
    if True:
        
        import matplotlib.gridspec as gridspec
        gs = gridspec.GridSpec(4, 9)
        ax1 = fig.add_subplot(gs[:, 0:3], projection='3d')
        ax1 = fix_axes(ax1)

        ax2 = fig.add_subplot(gs[:, 3:6], projection='3d')
        ax2 = fix_axes(ax2)

        ax3 = fig.add_subplot(gs[:, 6:9], projection='3d')
        ax3 = fix_axes(ax3)

        ax1.set_title(titles[0], y = 1)
        ax2.set_title(titles[1], y = 1)
        ax3.set_title(titles[2], y = 1)

    if False:
        ax1.plot_surface(a1.xc, a1.yc , a1.zc*0, facecolors = a1.veg_colors , 
                                rstride = 1, cstride = 1, alpha = 1,
                                linewidth=0, antialiased=True, shade=False)      

        ax2.plot_surface(a2.xc, a2.yc , a2.zc*0, facecolors = a2.veg_colors , 
                                rstride = 1, cstride = 1, alpha = 1,
                                linewidth=0, antialiased=True, shade=False)      

        ax3.plot_surface(a3.xc, a3.yc , a3.zc*0, facecolors = a3.veg_colors , 
                                rstride = 1, cstride = 1, alpha = 1,
                                    linewidth=0, antialiased=True, shade=False)      

    i_tr = np.where(sim1.t < sim1.tr*60)[-1][-1]

       
    project_micro = 1
    
    if project_micro == 1:    
        
        
        micro_max = max([a1.micro.max(), a2.micro.max(), a3.micro.max()])
        micro_min = max([a1.micro.min(), a2.micro.min(), a3.micro.min()])
        micro_norm = colors.Normalize(micro_min, micro_max)

        micro_c1 = cm.terrain(micro_norm(a1.micro ))    
        micro_c2 = cm.terrain(micro_norm(a2.micro ))            

        ax1.plot_surface(a1.xc, a2.yc, a1.micro/5, facecolors = micro_c1 , 
                            rstride = 1, cstride = 1, alpha = 0.4,
                            linewidth=0, antialiased=True, shade=False) 
        
        ax2.plot_surface(a1.xc, a2.yc, a2.micro/5, facecolors = micro_c2 , 
                            rstride = 1, cstride = 1, alpha = 0.4,
                            linewidth=0, antialiased=True, shade=False)         

        cbaxes = fig.add_axes([0.09, 0.25, 0.01, 0.45])  
        fmt = lambda x, pos: '{:.0f}'.format(x)
        format_colorbar( micro_max,"terrain", cbaxes, fmt, label = "terrain (cm)" , vmin = micro_min)
        
 
    if plot_infl == 1:
        
        
        Imax = np.max([a1.I.ravel().max(), a2.I.ravel().max(),a3.I.ravel().max()])
        I_norm = colors.Normalize(vmin = 0.0,  vmax=Imax)  
        DI_norm = colors.Normalize(vmin = a3.I.ravel().min(), vmax =  a3.I.ravel().max())

        I_colors1 = GnBu(I_norm(a1.I[-1])) 
        I_colors2 = GnBu(I_norm(a2.I[-1]))
        I_colors3 = GnBu(DI_norm(a3.I[-1]))

        cbaxes =  fig.add_axes([0.61, 0.25, 0.01, 0.45])  
        fmt = lambda x, pos: '{:.1f}'.format(x)
        format_colorbar( Imax, GnBu, cbaxes, fmt,
                        label = "$I$ (cm)" , vmin = 0)

        cbaxes = fig.add_axes([0.89, 0.25, 0.01, 0.45])  
        fmt = lambda x, pos: '{:.2f}'.format(x)
        format_colorbar( a3.I.ravel().max(), GnBu, 
                        cbaxes, fmt, label = "$\Delta I$ (cm)" ,
                        vmin = a3.I.ravel().min()*100)



        plot =  [ax1.plot_surface(a1.xc, a1.yc,  a1.zc, facecolors = I_colors1, 
                            rstride=1, cstride=1,cmap = GnBu, alpha = 1.,
                            linewidth=0,antialiased=True, shade=False),
                 ax2.plot_surface(a2.xc, a2.yc,  a2.zc, facecolors = I_colors2, 
                            rstride=1, cstride=1, alpha = 1.,
                            linewidth=0,antialiased=True, shade=False),
                 ax3.plot_surface(a3.xc, a3.yc,  a3.zc, facecolors = I_colors3, 
                            rstride=1, cstride=1, alpha = 1.,
                            linewidth=0,antialiased=True, shade=False), 
                ]

    else:
        
        hmax = np.max([a1.hc.ravel().max(), a1.hc.ravel().max(),a3.hc.ravel().max()])
        hmin = np.min([a1.hc.ravel().min(), a1.hc.ravel().min(),a3.hc.ravel().min()])    

        h_norm = colors.Normalize(vmin= hmin,   vmax= hmax)
        h_colors1 = GnBu(h_norm(a1.hc[i_tr]))
        h_colors2 = GnBu(h_norm(a2.hc[i_tr]))
        h_colors3 = GnBu(h_norm(a3.hc[i_tr]))

        ax1.plot_surface(a1.xc, a1.yc , a1.zc , facecolors = h_colors1, 
                            rstride = 1, cstride = 1, alpha = 1,
                            linewidth=0, antialiased=True, shade=False)      

        ax2.plot_surface(a1.xc, a1.yc , a1.zc, facecolors = h_colors2, 
                            rstride = 1, cstride = 1, alpha = 1,
                            linewidth=0, antialiased=True, shade=False)   

        ax3.plot_surface(a1.xc, a1.yc , a1.zc , facecolors = h_colors3, 
                            rstride = 1, cstride = 1, alpha = 1,
                            linewidth=0, antialiased=True, shade=False)   


        cbaxes = fig.add_axes([0.89, 0.25, 0.01, 0.45])  
        fmt = lambda x, pos: '{:.1f}'.format(x)
        format_colorbar( hmax, GnBu, cbaxes, fmt, label = "$h$ (cm)" )

 
    ax1.set_zlim(0, (a1.zc*1.2 ).max())
    ax2.set_zlim(0, (a2.zc*1.2 ).max())
    ax3.set_zlim(0, (a3.zc*1.2).max())

    ax1.view_init(25, 20)
    ax2.view_init(25, 20)
    ax3.view_init(25, 20)

    # plt.suptitle(sim.name, fontsize = 14)
    return fig

