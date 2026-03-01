from plot_SWOF import *
from scipy import interpolate

site_names = {'EM': 'Audubon research ranch',
 'PC': 'Audubon research ranch',
 'ER1': 'Empire Ranch',
 'ER2': 'Empire Ranch',
 'ER3': 'Empire Ranch',
 'ER4G': 'Empire Ranch',
 'ER4S': 'Empire Ranch',
 'ER5': 'Empire Ranch',
 'PCE': 'Porter Canyon',
 'PCW': 'Porter Canyon',
 'AB': 'San Rafael Valley',
 'SA': 'San Rafael Valley',
 'TA': 'San Rafael Valley',
 'WI': 'San Rafael Valley',
 'K2': 'WGEW',
 'K3': 'WGEW',
 'CR': 'WGEW',
 'LH1': 'WGEW',
 'LH2': 'WGEW',
 'LH3': 'WGEW',
 'YG1': 'Young, AZ',
 'YG2': 'Young, AZ',
 'YG3': 'Young, AZ'}

# colors
blue = (0.2980392156862745, 0.4470588235294118, 0.6901960784313725)
orange = (0.8666666666666667, 0.5176470588235295, 0.3215686274509804)
green = (0.3333333333333333, 0.6588235294117647, 0.40784313725490196)


def nrmse(x, y):
    """
    """
    return np.sqrt(np.mean((x-y)**2))/np.mean(x)


### Correction factors

def get_c_q(m):  
    """
    c_q = <h><U>/q_L
    """
    return (m+1)**2 *(2 + m)**(-1) * (2*m + 1)**(-1)

def get_c_q_a(m, a):
    """
    c_q = <h>_a<U>_a/q_L
    """
    return (m+1)**2*(2 + m)**(-1)*(2*m + 1)**(-1)*(
        1-a**((2*m+1)/(m+1)))*(1-a**((m+2)/(m+1)))/(1-a)**2    

def get_c_r(m):
    """
    c_r = r_o / <r>
    """
    return (1+m)**(1 - m)*(2*m + 1)**(-1)*(2 + m)**(m) 

def get_c_r_a(m,a):
    """
    c_r = r_o / <r>_a
    """
    return (1+m)**(1 - m)*(2*m + 1)**(-1)*(2 + m)**(m)*\
        (1 - a)**(m-1)*( 1 - a**((2*m + 1)/(1+m)))*(1 - a**((2 + m)/(m+1)))**(-m)

def get_c_U(m, a):
    """
    c_U = <U>_a / <U> 
    """
    return (1 - a)**(-1)*( 1 - a**((2*m + 1)/(1+m)) )

def get_correction(m, a):
    """
    combine correction factors
    """
    return get_c_U(m, a)**(1+m)*get_c_r(m)*get_c_q(m)**m #get_c_U(m, a)*get_c_r(m)*get_c_q(m)**m

### Naming conventions

def r_to_Knf(summary):
    """
    convert simulation scheme to abbreviation
    """
    name = summary.scheme.iloc[0].lower()
    if name == 'laminar_dw':
        return "K_o"
    if name == 'manning':
        return "n_o"
    if name == 'dw':
        return "f_o"    
    else:
        return "r_o"


def m_to_name(m):
    """
    exponent to name
    """
    if  m == 0:
        return 'Distributed'
    if m == 0.5:
        return 'Darcy Weisbach'
    if m == 0.66:
        return 'Mannings'
    if m == 2/3:
        return 'Mannings'
    if m == 2:
        return 'Laminar'

def names_to_m(name):
    """
    name to exponent m
    """
    if name ==  'Distributed':
        return 0
    if name == 'Darcy Weisbach':
        return 0.5
    if name == 'Mannings':
        return 0.66
    if name == 'Laminar':
        return 2

### Interpolate gridded output

def interp_6(sim, fld, add = 1):
    """
    Estimated fld at x = 6
    """
    ind6 = int(6//sim.dx) - 1

    if len(sim[fld].shape) == 2:    

        if sim.l == 6:  

            return sim[fld][:, ind6]

        else:   
            return (sim[fld][:, ind6] + sim[fld][:, ind6+add])[0]/2
    else:
        
        if sim.l == 6:  

            return sim[fld][:, 0, ind6] 

        else:   
            
            return (sim[fld][:, 0, ind6] + sim[fld][:, 0, ind6+add])/2        
        
def interp_27_to_6(sim, fld, add = 1):
    """
    Computes avg(fld) between x = 2.7 and x = 6
    """
    ind6 = int(6//sim.dx) -1
    if sim.dx == 1:
        ind27 = round(2.7/sim.dx)-1
    else:
        ind27 = int(2.7//sim.dx)-1
        
    if len(sim[fld].shape) == 2:

        if sim.l == 6:  
            return sim[fld][:, ind27:ind6+1] 

        else:   
            return (sim[fld][:, ind27:ind6+1] + sim[fld][:, ind27+add:ind6+1+add])/2
        
    else:
       
        if sim.l == 6:  

            return sim[fld][:, 0, ind27:ind6+1].mean(1)

        else:  

            return (sim[fld][:, 0, ind27:ind6+1] + sim[fld][:, 0, ind27+add:ind6+1+add]).mean(1)/2
        
def interp_to_6(sim, fld, add = 1, start = 0):
    """
    Computes avg(fld) between start and x = 6
    """    
    
    ind6 = int(6//sim.dx) - 1
    
    if len(sim[fld].shape) == 2:    
        if sim.l == 6:  
            return sim[fld][:, start:ind6+1] 

        else:   
            return (sim[fld][:, start:ind6+1] + sim[fld][:, start+add:ind6+add+ 1])/2
    else:        
        if sim.l == 6:  
            return sim[fld][:, 0,  start:ind6+1].mean(1) 

        else:             
            return (sim[fld][:,0, start:ind6+1] + sim[fld][:, 0, start+add:ind6+1+add]).mean(1)/2

### Plot utilities

from PIL import Image

def load_PIL_img(label, multiple = 1):
    
    rootdir = '/Users/octaviacrompton/Google_Drive_quatratavia/data/polyakov2017rainfall/Site photos/labeled'
    
    expdir = "{0}/{1:.0f}".format(rootdir, label)
    names = os.listdir(expdir)
    names = [n for n in names if 'JPG' in n]    
    img_path = os.path.join(expdir, names[0])
    
    img = Image.open(img_path)

    if img.size[1] > img.size[0]:
        img = img.rotate(-90, expand=True)

    basewidth = 500

    wpercent = (basewidth/float(img.size[0]))
    hsize = int((float(img.size[1])*float(wpercent)))
    img = img.resize((basewidth,hsize))
    img = np.array(img)


    if len(names) > 1 and multiple == 1:
        img_path = os.path.join(expdir, names[1])
        img2 = Image.open(img_path)
        
        if img2.size[1] > img2.size[0]:
            img2 = img2.rotate(-90, expand=True)
        wpercent = (basewidth/float(img2.size[0]))
        
        #hsize = int((float(img2.size[1])*float(wpercent)))
        img2 = img2.resize((basewidth,hsize))
        img2 = np.array(img2)
        
        img = np.concatenate((img, img[:, :200]*0 + 255, img2, img[:, :100]*0 + 255), axis = 1)
    
    return img

         
def trim_df(df):
    """
    remove indices where rain intensity increases
    """
    inds =   np.where((np.diff(df.Precipitation, prepend = 0) != 0))[0] 
    # drop first 2 indices after precip increase
    inds = np.hstack((inds +3, inds+1, inds+2 ))
    df = df.drop(df.index[inds])
    # require time > 5, rainfall > 0
    return df[(df.Time > 5) & (df.Precipitation >0) &   (np.diff(df.Precipitation, prepend = 0) == 0)]
        
    
def visualize_case_row(folder, allsum, meta, criteria = 'w_eff', r = 'K', r_fld = 'K_avg'):
    
    
    label = float(folder.split("-")[-1])
    ID = folder.split("-")[0]
    case = meta.query("label == {0}".format(label)).iloc[0]

    rootdir = '/Users/octaviacrompton/Google_Drive_quatratavia/data/polyakov2017rainfall/rain'
    df = pd.read_csv('{0}/{1}.csv'.format(rootdir, folder.strip("_n")))
    df['Infiltration'] = df['Precipitation']  - df['Runoff'] 
    dft = trim_df(df)
    row = meta.query("folder == '{0}'".format(folder)).iloc[0]

    summary = allsum[folder]
    points = df[['Runoff', 'Velocity_mean']].dropna()
    
    
    best = summary.sort_values(criteria).tail(10)
    sim = summary.sort_values(criteria).tail().iloc[-1]

    ## Prepare U - Q comparison
    vc6 = interp_6(sim, 'vc')
    hc6 = interp_6(sim, 'hc')

    hydro = hc6*vc6*3.6e5/6       

    v = interp_27_to_6(sim, 'vc')
    if len(hydro) > len(v):
        v = np.insert(v, -1, 0)
    f = interpolate.interp1d(hydro, v) 

    inds = np.where(( points.Runoff/10 < hydro.max()) & ( points.Runoff/10 > hydro.min()))[0]

    Velocity = f(points.Runoff.iloc[inds]/10)  
    mpl.rcParams['figure.dpi'] = 300
    
    fig, axes = plt.subplots(1,3, figsize = (18,4), 
                             sharey = False, sharex = False,  constrained_layout = True)    
    axes = axes.ravel()
    fig.subplots_adjust(hspace = 0.4)
    fig.subplots_adjust(wspace = 0.2)
    
    axes[2].scatter( interp_27_to_6(sim, 'vc')*100, hydro[:], color= 'grey', alpha = 0.4)
    axes[2].scatter(Velocity*100,  points.Runoff.iloc[inds]/10, color= 'grey',
                   s = 80, alpha = 1, label = "SVE model")
    axes[2].scatter(points.Velocity_mean*100, points.Runoff/10, 
                   color = orange, marker = 's', s = 80, label = "observed")  
    axes[2].set_xlim(0, np.percentile(sim.vc, 99)*100)
    axes[2].legend(loc = 'upper left')
    axes[2].set_ylabel("$q_L$ (cm/hr)") 
    axes[2].set_xlabel("$U_E$ (cm/s)")

    vc6 = interp_6(sim, 'vc')
    hc6 = interp_6(sim, 'hc')
    hydro = hc6*vc6*3.6e5/6              

    axes[1].plot(df.Time, df.Precipitation/10, 'o', ms = 5, label = "rainfall $p$")
    axes[1].plot(df.Time, df.Precipitation/10, '-', color = 'k', alpha= 0.1)

    axes[1].plot(df.Time, df.Runoff/10, 'o-', ms = 5, label = "runoff $q_L$")    
    axes[1].plot(sim.t/60, hydro, '-', color= 'grey', alpha = 1, label = " model $q_L$") 

    #axes[1].plot(dft.Time, dft.Infiltration/10, 'o', ms = 5, alpha = 0.5, label = "infiltration $I$")    
    #axes[1].plot(dft.Time, dft.Infiltration/10, ':', color = 'grey', alpha = 0.5)
    
    
    axes[1].set_xlabel("$t$ (min)")
    axes[1].legend(loc = 'upper left')
    #axes[1].legend(loc='center left', bbox_to_anchor=(0.95, 0.5))
    max_y = max(axes[2].get_ylim()[1],axes[1].get_ylim()[1])*1.25
    
    axes[1].set_ylim(axes[2].get_ylim()[0], max_y)
    axes[2].set_ylim(axes[2].get_ylim()[0], max_y)
    axes[1].set_ylabel("cm/hr") 
   
    import calendar
    title = "{0} ({1}),  {4:s} {3:.0f},  plot = {2:.0f},".format(
        site_names[ID], ID, case.Plot, case.Year,  calendar.month_abbr[int(case.Month)])
    
    axes[0].text(-.02, 1.15, title, transform=axes[0].transAxes, ha="left", va="top",
                fontsize = 18, style = 'normal');
    
    c_c = get_correction(float(sim.m), 2.7/6)

    stats_str = folder + "; " + row.Vegetation + ", " + row.Soil_texture + ", condition = " + row.Condition
    
    stats_str += "\n" + "$K_s$ = {0:.2f} cm/hr".format(sim.Ks_v) +\
                "; $dI/dt={0:.2f}$ cm/hr [{1:.2f}, {2:.2f}] ".format(row.I_trend, row.I_trend_low, row.I_trend_high)
    
    stats_str += "\n\n" + "estimate {3} = {0:.2f}, calibrate {3} = {1:.2f} +/- {2:.2f} ".format(
        case[r_fld]*c_c, sim.alpha_v, best.alpha_v.std() , r)
    
    stats_str += "; " + "tested range = [{0:.2f},{1:.2f}]".format( summary[r].min(), summary[r].max())
    
    stats_str += "\n" +"$w_e$= {0:.2f}, $U_e$ =  {1:.2f}, NSE = {2:.2f}".format(
        sim.w_eff, sim.U_eff, sim.NSE)
    
    img = load_PIL_img(label, 0)

    gs = axes[0].get_gridspec()

    axes[0].axis('off')

    axbig = fig.add_subplot(gs[:1])
    axbig.imshow(img)
    axbig.axis('off') 
    
    axes[0].text(-.05 , -0.3, stats_str, transform=axes[0].transAxes, ha="left", va="top",
               style = 'normal');
    
    return stats_str

def visualize_case(folder, allsum, meta, r = 'K', criteria = 'w_eff', r_fld = 'K_avg'):
    """
    """
    df = pd.read_csv(
        '/Users/octaviacrompton/Google_Drive_quatratavia/data/polyakov2017rainfall/rain/{0}.csv'.format(
            folder.strip("_n")))

    summary = allsum[folder]
    
    label = folder.split("-")[2]

    ID = folder.split("-")[0]
    print (folder, label, ID)
    case = meta.query("label == {0}".format(label)).iloc[0]

    best = summary.sort_values(criteria).tail(10)    
    points = df[['Runoff', 'Velocity_mean']].dropna()

    ## Prepare U - Q comparison
    sim = summary.sort_values(criteria).tail().iloc[-1]
    
    vc6 = interp_6(sim, 'vc')
    hc6 = interp_6(sim, 'hc')
        
    hydro = hc6*vc6*3.6e5/6       

    f = interpolate.interp1d(hydro[:], sim.vc[:, 0, 5:].mean(1))

    v = interp_27_to_6(sim, 'vc')
    
    if len(hydro) > len(v):
        v = np.insert(v, -1, 0)
        
    f = interpolate.interp1d(hydro, v) 

    inds = np.where(( points.Runoff/10 < hydro.max()) & ( points.Runoff/10 > hydro.min()))[0]

    Velocity = f(points.Runoff.iloc[inds]/10) 
     
    # Make figure
    fig, axes = plt.subplots(1,2, figsize = (10, 3), sharey = True)

    axes[0].scatter( interp_27_to_6(sim, 'vc')*100, hydro[:], color= 'grey',
                    alpha = 0.4)
    axes[0].scatter(Velocity*100,  points.Runoff.iloc[inds]/10,  
                    color= 'grey',  s = 80, alpha = 1, label = "SVE model")
    axes[0].scatter(points.Velocity_mean*100,  points.Runoff/10, 
                   color = orange, marker = 's',   s = 80, label = "observed")  
    
    axes[0].legend()
    axes[0].set_ylabel("$q_L$ (cm/hr)", fontsize = 14)
    axes[0].set_xlabel("$U_E$ (cm/s)", fontsize = 14)
    axes[0].set_title("{0} ({1}), plot = {2:.0f}, year = {3:.0f}".format(
        site_names[ID], ID, case.Plot, case.Year), fontsize = 14)
    
    axes[1].plot(df.Time, df.Precipitation/10, 'o', alpha = 0.7, label = "rainfall $p$")
    axes[1].plot(df.Time, df.Runoff/10, '.-', alpha = 0.6, label = "runoff $q_L$")    
    axes[1].plot(sim.t/60, hydro, ':', color = 'grey', alpha = 1, label = "SVE model $q_L$")   
    axes[1].plot(df.Time, df.Precipitation/10, 'k', alpha = 0.1)
    
    axes[1].set_xlabel("$t$ (min)")
    axes[1].legend(loc='center left', bbox_to_anchor=(0.95, 0.5), fontsize= 15)
    
    c_c = get_correction(float(sim.m), 2.7/6)

    stats_str = "estimate {3} = {0:.2f}, calibrate {3} = {1:.2f} +/- {2:.2f} ".format(
        case[r_fld]*c_c, sim.alpha_v, best.alpha_v.std() , r)
    stats_str += "\n" + "tested range = [{0:.2f},{1:.2f}]".format(summary[r].min(), summary[r].max())
    stats_str += "\n" + "analytic c_c = {0:.2f}, calibration c_c = {1:.2f}".format(c_c,
                                                                     best.alpha_v.mean()/case[r_fld] )
    stats_str += "\n" +"w_eff = {0:.2f}, U_eff = {1:.2f}, NSE = {2:.2f}".format(
        sim.w_eff, sim.U_eff, sim.NSE)

    print( stats_str)
    
# def visualize_case_im(folder, allsum, meta, criteria = 'w_eff', r = 'K', r_fld = 'K_avg'):

#     label = float(folder.split("-")[-1])
#     ID = folder.split("-")[0]
#     case = meta.query("label == {0}".format(label)).iloc[0]

#     rootdir = '/Users/octaviacrompton/Google_Drive_quatratavia/data/polyakov2017rainfall/rain'
#     df = pd.read_csv('{0}/{1}.csv'.format(rootdir, folder.strip("_n")))
    
#     summary = allsum[folder]
#     points = df[['Runoff', 'Velocity_mean']].dropna()

#     best = summary.sort_values(criteria).tail(10)
#     sim = summary.sort_values(criteria).tail().iloc[-1]

#     ## Prepare U - Q comparison
#     vc6 = interp_6(sim, 'vc')
#     hc6 = interp_6(sim, 'hc')

#     hydro = hc6*vc6*3.6e5/6       

#     v = interp_27_to_6(sim, 'vc')
#     if len(hydro) > len(v):
#         v = np.insert(v, -1, 0)
#     f = interpolate.interp1d(hydro, v) 

#     inds = np.where(( points.Runoff/10 < hydro.max()) & ( points.Runoff/10 > hydro.min()))[0]

#     Velocity = f(points.Runoff.iloc[inds]/10)  

#     mpl.rcParams['figure.dpi'] = 200

#     fig, axes = plt.subplots(2,2, figsize = (10,7), sharey = True, sharex = False,  constrained_layout = True)    
#     axes = axes.ravel()
#     fig.subplots_adjust(hspace = 0.4)
#     fig.subplots_adjust(wspace = 0.1)
    
#     axes[0].scatter( interp_27_to_6(sim, 'vc')*100, hydro[:], color= blue,
#                     alpha = 0.4)
#     axes[0].scatter(Velocity*100,  points.Runoff.iloc[inds]/10,  
#                     color= blue,
#                    s = 80, alpha = 1, label = "SVE model")
#     axes[0].scatter(points.Velocity_mean*100, points.Runoff/10, 
#                    color = orange, marker = 's', s = 80, label = "observed")  
#     axes[0].set_xlim(0, np.percentile(sim.vc, 99)*100)
#     axes[0].legend()
#     axes[0].set_ylabel("$q_L$ (cm/hr)", fontsize = 14)
#     axes[0].set_xlabel("$U_E$ (cm/s)", fontsize = 14)

#     vc6 = interp_6(sim, 'vc')
#     hc6 = interp_6(sim, 'hc')
#     hydro = hc6*vc6*3.6e5/6              

#     axes[1].plot(sim.t/60, hydro, '-', alpha = 1, label = "SVE model")   
#     axes[1].plot(df.Time, df.Runoff/10, '.-', label = "observed $q_L$")    
#     axes[1].plot(df.Time, df.Precipitation/10, '.', label = "rainfall")
#     axes[1].plot(df.Time, df.Precipitation/10, '-', alpha = 0.1, c = 'k')    
    
#     axes[1].set_xlabel("$t$ (min)")
#     axes[1].legend(loc='center left', bbox_to_anchor=(0.95, 0.5), fontsize= 15)

#     t = plt.suptitle("{0} ({1}), plot = {2:.0f}, year = {3:.0f}".format(
#         site_names[ID], ID, case.Plot, case.Year), fontsize = 14)
#     t.set_y(0.93)

#     c_c = get_correction(float(sim.m), 2.7/6)

#     stats_str = "estimate {3} = {0:.2f}, calibrate {3} = {1:.2f} +/- {2:.2f} ".format(
#         case[r_fld]*c_c, sim.alpha_v, best.alpha_v.std() , r)
#     stats_str += "\n" + "tested range = [{0:.2f},{1:.2f}]".format( summary[r].min(), summary[r].max())
#     stats_str += "\n" + "analytic c_c = {0:.2f}, calibration c_c = {1:.2f}".format(c_c,
#                                                                      best.alpha_v.mean()/case[r_fld] )
#     stats_str += "\n" +"w_eff = {0:.2f}, U_eff = {1:.2f}, NSE = {2:.2f}".format(
#         sim.w_eff, sim.U_eff, sim.NSE)
    
#     img = load_PIL_img(label)

#     gs = axes[2].get_gridspec()
        
#     for ax in axes[2:]:
#         ax.remove()

#     axes[2].axis('off')

#     axbig = fig.add_subplot(gs[2:])
#     axbig.imshow(img)
#     axbig.axis('off') 
    
#     return stats_str
def visualize_case_im(folder, allsum, meta, criteria = 'w_eff', r = 'K', r_fld = 'K_avg'):

    label = float(folder.split("-")[-1])
    ID = folder.split("-")[0]
    case = meta.query("label == {0}".format(label)).iloc[0]

    rootdir = '/Users/octaviacrompton/Google_Drive_quatratavia/data/polyakov2017rainfall/rain'
    df = pd.read_csv('{0}/{1}.csv'.format(rootdir, folder.strip("_n")))
    
    summary = allsum[folder]
    points = df[['Runoff', 'Velocity_mean']].dropna()

    best = summary.sort_values(criteria).tail(10)
    sim = summary.sort_values(criteria).tail().iloc[-1]

    ## Prepare U - Q comparison
    vc6 = interp_6(sim, 'vc')
    hc6 = interp_6(sim, 'hc')

    hydro = hc6*vc6*3.6e5/6       

    v = interp_27_to_6(sim, 'vc')
    if len(hydro) > len(v):
        v = np.insert(v, -1, 0)
    f = interpolate.interp1d(hydro, v) 

    inds = np.where(( points.Runoff/10 < hydro.max()) & ( points.Runoff/10 > hydro.min()))[0]

    Velocity = f(points.Runoff.iloc[inds]/10)  

    mpl.rcParams['figure.dpi'] = 200

    fig, axes = plt.subplots(2,2, figsize = (10,7), sharey = True, sharex = False,  constrained_layout = True)    
    axes = axes.ravel()
    fig.subplots_adjust(hspace = 0.4)
    fig.subplots_adjust(wspace = 0.1)
    
    axes[0].scatter( interp_27_to_6(sim, 'vc')*100, hydro[:], color= 'grey',
                    alpha = 0.4)
    axes[0].scatter(Velocity*100,  points.Runoff.iloc[inds]/10,  
                    color= 'grey',
                   s = 80, alpha = 1, label = "SVE model")
    axes[0].scatter(points.Velocity_mean*100, points.Runoff/10, 
                   color = orange, marker = 's', s = 80, label = "observed")  
    axes[0].set_xlim(0, np.percentile(sim.vc, 99)*100)
    axes[0].legend()
    axes[0].set_ylabel("$q_L$ (cm/hr)", fontsize = 14)
    axes[0].set_xlabel("$U_E$ (cm/s)", fontsize = 14)

    vc6 = interp_6(sim, 'vc')
    hc6 = interp_6(sim, 'hc')
    hydro = hc6*vc6*3.6e5/6              

    axes[1].plot(df.Time, df.Precipitation/10, '.', label = "rainfall $p$")
    axes[1].plot(df.Time, df.Runoff/10, '.-', label = "runoff $q_L$")    
    axes[1].plot(sim.t/60, hydro, '-', color= 'grey', alpha = 1, label = "SVE model $q_L$")   
    axes[1].plot(df.Time, df.Precipitation/10, ':', alpha = 0.5, c = 'grey')    
    
    axes[1].set_xlabel("$t$ (min)")
    axes[1].legend(loc='center left', bbox_to_anchor=(0.95, 0.5), fontsize= 15)

    import calendar
    title = "{0} ({1}),  {4:s} {3:.0f},  plot = {2:.0f},".format(
        site_names[ID], ID, case.Plot, case.Year,  calendar.month_abbr[int(case.Month)])
    

    t = plt.suptitle(title, fontsize = 14)
    t.set_y(0.93)

    c_c = get_correction(float(sim.m), 2.7/6)
    print (label)


    img = load_PIL_img(label)

    gs = axes[2].get_gridspec()
        
    for ax in axes[2:]:
        ax.remove()

    axes[2].axis('off')

    axbig = fig.add_subplot(gs[2:])
    axbig.imshow(img)
    axbig.axis('off') 

    stats_str = folder + "; " + row.Vegetation + ", " + row.Soil_texture + ", condition = " + row.Condition
 

    stats_str += "\n" + "$K_s$ = {0:.2f} cm/hr".format(sim.Ks_v) +\
                "; $dI/dt={0:.2f}$ cm/hr [{1:.2f}, {2:.2f}] ".format(row.I_trend, row.I_trend_low, row.I_trend_high)

    stats_str += "\n\n" + "estimate {3} = {0:.2f}, calibrate {3} = {1:.2f} +/- {2:.2f} ".format(
        case[r_fld]*c_c, sim.alpha_v, best.alpha_v.std() , r)
    
    stats_str += "; " + "tested range = [{0:.2f},{1:.2f}]".format( summary[r].min(), summary[r].max())
    

    # stats_str += "\n" + "analytic c_c = {0:.2f}, calibration c_c = {1:.2f}".format(c_c,
    #                                                                 best.alpha_v.mean()/case[r_fld] )
    stats_str += "\n" +"$w_e$= {0:.2f}, $U_e$ =  {1:.2f}, NSE = {2:.2f}".format(
        sim.w_eff, sim.U_eff, sim.NSE)
    

    return stats_str


 
def get_steady( df):
    """
    """
    
    first = np.where(df.Runoff > 0)[0][0]
    df = df.iloc[first:]

    # Remove gaps in rain
    if False:
        after_break = np.where((df.Precipitation == 0) & (df.Time < 60))[0][0]
        df = df.iloc[after_break:]
    
        
    # inds = True for 1st value of each new/increased precip value
    inds = np.diff(df.Precipitation, prepend = df.Precipitation.iloc[0]) != 0

    # inds = integers indicating final row before precip increase
    inds = np.where(inds)[0] - 1
    
    # concatenate with the rows before index, and the rows before that
    inds = sorted(np.concatenate((inds, inds - 1)) )

    # limit the data to these points
    steady = df.iloc[inds].reset_index()
    
    rain_end = np.where(np.diff(steady.Precipitation, prepend = 0) < 0)[0]
    steady = steady.drop(rain_end, axis = 0)
    
    steady = steady.query("Time > 15 and Precipitation > 0")
    
    steady['Infiltration'] = steady['Precipitation']  - steady['Runoff'] 
    steady['Time_hr'] = steady['Time']/60
    
    return steady

def group_by_rain(df, lag = 1):
    """
    Functions to estimate Ksat, infiltration
    """
    if False:
        after_break = np.where((df.Precipitation == 0) & (df.Time < 60))[0][0]
        df = df.iloc[after_break:]
        
    df = df.query("Precipitation > 0 and Time > 15 and Runoff > 0").reset_index()
    
    # Drop datapoint at the time of the precip increase and the following lag measurements
    inds = np.where(np.diff(df.Precipitation, prepend = df.Precipitation.iloc[0]) > 0 )[0] 

    if lag == 2:
        inds = np.concatenate((inds, inds+ 1))
        
    elif lag == 3:
        inds = np.concatenate((inds, inds+ 1, inds+ 2))

    elif lag == 4:
        inds = np.concatenate((inds, inds+ 1, inds+ 2, inds+ 3))

    df = df.drop(index = sorted(inds)).reset_index()
    
    df = df[np.diff(df.Runoff, append = 0) < 10]

    # trim the front and end of the timeseries
    first = np.where(np.diff(df.Runoff)  < 10)[0][0] 
    last = np.where(df.Precipitation  > 0)[0][-1]  + 1
    grouped = df.iloc[first:last]
    
    #grouped = grouped.groupby("Precipitation").mean().reset_index()[['Time', 'Precipitation', 'Runoff']]
    grouped['Infiltration'] = grouped['Precipitation']  - grouped['Runoff']
    grouped['Time_hr'] = grouped['Time']/60
    return grouped

# def group_by_rain(df, lag = 0):
#     """
#     Functions to estimate Ksat, infiltration
#     """

#     # Drop datapoint at the time of the precip increase and the following lag measurements
#     inds = np.where(np.diff(df.Precipitation, prepend = df.Precipitation.iloc[0]) > 0 )[0] 

#     if lag == 1:
#         inds = np.concatenate((inds, inds+ 1))
        
#     elif lag == 2:
#         inds = np.concatenate((inds, inds+ 1, inds+ 2))

#     elif lag == 3:
#         inds = np.concatenate((inds, inds+ 1, inds+ 2,  inds+ 3))

    
#     df = df.reset_index().drop(index = sorted(inds)  )
    
#     # trim the front and end of the timeseries
#     first = np.where(df.Precipitation  > 0)[0][0]  + lag +2
#     last = np.where(df.Precipitation  > 0)[0][-1]  + 1
#     grouped = df.iloc[first:last]
    
#     #grouped = grouped.groupby("Precipitation").mean().reset_index()[['Time', 'Precipitation', 'Runoff']]
#     grouped['Infiltration'] = grouped['Precipitation']  - grouped['Runoff']
#     grouped['Time_hr'] = grouped['Time']/60
#     return grouped

# def get_steady(df):
#     """
#     """
#     tr = df.Time[df.Precipitation > 0].iloc[-1]
#     # Remove gaps in rain
#     first = np.where(df.Runoff > 0)[0][0]
#     df = df.iloc[first:]
#     df = df.query("(Precipitation > 0 and Time > 2)  or Time > {0} ".format(tr))    
#     # inds = True for 1st value of each new/increased precip value
#     inds = np.diff(df.Precipitation, prepend = df.Precipitation.iloc[0]) != 0

#     # inds = integers indicating final row before precip increase
#     inds = np.where(inds)[0] - 1
    
#     # concatenate with the rows before index, and the rows before that
#     inds = sorted(np.concatenate((inds, inds - 1)) )

#     # limit the data to these points
#     steady = df.reset_index().iloc[inds]
#     steady['Infiltration'] = steady['Precipitation']  - steady['Runoff'] 
#     steady['Time_hr'] = steady['Time']/60
    
#     return steady