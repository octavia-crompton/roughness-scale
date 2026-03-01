import numpy as np
import pandas as pd


def insert_fld(summary, fld):
    """
    Add object field to pandas dataframe
    """
    summary[fld] = 0
    summary[fld] = summary[fld].astype(object)
    
    return summary 

def increment(indirect, i, j, sim):
    """ 
    march up the hillslope starting at the outlet
    """
    patch_LB = sim.patch_LB.astype(int)
    patch_LV = sim.patch_LV.astype(int)
    veg = sim.veg.astype(int)
    
    vegL = patch_LV[i, j]

    # length of the next upslope bare soil patch is...
    upslope = np.maximum(- 1 + j - vegL, - sim.Nycell)

    bareL = patch_LB[i, upslope]
        
    survives =  bareL*(sim.p-sim.Ks_b) + vegL*(sim.p - sim.Ks_v)
    
    indirect[i, j] = survives
    j = j - (bareL + vegL)  - 1

    return indirect, j


def get_source(sim):
    """
    Estimate the potential source area 
    """
    indirect = np.zeros_like(sim.veg)
    direct = np.zeros_like(sim.veg)
    source = np.ones_like(sim.veg)

    patch_LB = sim.patch_LB.astype(int)
    patch_LV = sim.patch_LV.astype(int)
    veg = sim.veg.astype(int)

    for i in range(sim.Nxcell):

        if veg[i, -1] == 1:
            j = -1
            for iterate in range(1000):    
                indirect, j = increment(indirect, i, j, sim)
                
                if j < - sim.Nycell:
                    break

        elif veg[i, -1] == 0:

            bareL =  patch_LB[i,-1]
            direct[i, -1] = bareL*(sim.p-sim.Ks_b) # m2 cm/hr
            
            j = np.maximum(- 2 - bareL, - sim.Nycell)

            for iterate in range(1000):    
                indirect, j = increment(indirect, i, j, sim)

                if j < - sim.Nycell:
                    break

    for i in range(sim.Nxcell):
        strip = indirect[i, :]
        ext = np.where(strip < 0)[0]
        if len(ext) > 0:
            source[i, :ext[-1]+2] = 0
            indirect[i, :ext[-1]+1] = 0    
        
    return indirect, direct, source

def add_source(summary):
    """
    Adds connectivity metrics to the dataframes:

    """
    for fld in ["upslope_array", "direct_array", "indirect_array"]:
        
        summary = insert_fld(summary, fld)
    
    for key in summary.index:
        sim = summary.loc[key]
        
        indirect, direct, source = get_source(sim)
        summary.at[key, "source_UI"] = source.mean()
        summary.at[key, "upslope_array"] = source
        summary.at[key, "direct"] = direct.sum()/(sim.p*sim.Nxcell*sim.Nycell)        
        # summary.at[key, "direct_array"] = direct
        summary.at[key, "indirect"] = indirect.sum()/(sim.p*sim.Nxcell*sim.Nycell)               
        # summary.at[key, "indirect_array"] = indirect                
        summary.at[key, "runoff_potential"] =  (indirect.sum() + direct.sum())/(sim.p*sim.Nxcell*sim.Nycell)

    return summary


##### Particle functions 

def intx(inpt, Nxcell):
    i = int(round(inpt))
    if i > Nxcell - 1:
        i = i- 1
    return i

def inty(inpt, Nycell):
    i = int(round(inpt))
    if i > Nycell - 1:
        i = i - 1
    return i

def tracer_init(sim, N = 1000, initialize = False, t_min = 1, y_min = 1, seed = None):
    """
    initialize the particle tracers: starting times between 1 and the end of storm 

    to : time index
    y_init, x_init : y, x coordinates
    """
    if seed:

        np.random.seed(seed)
    else: 
        np.random.seed(0)
    to = np.random.randint( t_min, len(sim.t[sim.t <= sim.tr*60]),  N)

    if sim.L > 2:
        x_init = np.random.randint( 1e-4, sim.L, N) - 0.001 + np.random.rand(N)
    else: 
        x_init =  np.random.rand(N)*sim.L 
    if sim.l > 2:
        y_init = np.random.randint( sim.dx*y_min, sim.l, N) - 0.001 + np.random.rand(N)           
    else: 
        y_init = (np.random.rand(N)*(sim.l-sim.dx) + sim.dx ) 
        
    
    if initialize == 'top':

        y_init = np.random.randint(1, 2,  N)
        x_init = np.linspace( 1e-4, sim.L, N) - 0.001
        # (np.random.rand(N)*sim.dx*4)

    elif initialize == 'top_t0':

        to = np.random.randint(1, 2,  N)
        y_init = ((np.random.rand(N))*sim.dx)        

    elif initialize == "t0":
        
        to = np.random.randint(1, 5,  N)
        y_init = ((np.random.rand(N))*sim.dx)             

    elif initialize == "flume":

        to = np.random.randint(1, 5,  N)
        y_init = ((np.random.rand(N))*sim.dx)             


    elif initialize == "along":
        # at the center of the flume... 
        x_init = np.random.randint( sim.L/2, sim.L/2+1,  N) - 0.001

    x_init[x_init >=  sim.L - sim.dx] = sim.L - sim.dx
    y_init[y_init >= sim.l - sim.dx] = sim.l - sim.dx

    dx = sim.dx
    xo = x_init
    yo = y_init


    Nx = sim.Nxcell
    Ny = sim.Nycell

    # update depth and velocity

    ho = [sim.hc[to[i], int(xo[i]/dx), min(int(yo[i]/dx), Ny - 1)] for i in range(len(x_init))]
    vo = [sim.vc[to[i], int(xo[i]/dx), min(int(yo[i]/dx), Ny - 1)] for i in range(len(x_init))]
    uo = [sim.uc[to[i], int(xo[i]/dx), min(int(yo[i]/dx), Ny - 1)] for i in range(len(x_init))]
    veg = [sim.veg[int(xo[i]/dx),   min(int(yo[i]/dx), Ny - 1)] for i in range(len(x_init))]

    loop = np.ones_like(x_init)*0

    # label particles
    label = np.arange(len(ho)).astype(int)

    positions = {}
    for l in label:
        positions[l] = {"xo" : [xo[l]], 
                        "yo" : [yo[l]],
                        "ho" : [ho[l]],
                        "uo" : [uo[l]],
                        "vo" : [vo[l]],
                        "t" :  [sim.t[to[l]]],
                        "loop" : [loop[l]],
                        "veg" : [veg[l]]
                       }
    return label, ho, vo, uo, xo, yo, loop, positions, to


def step(label, ho, vo, uo, xo, yo, loop, positions, to, sim, h_min, t_scale, periodic = 0):
    """
    """
    dx = sim.dx
    Nx = sim.Nxcell
    Ny = sim.Nycell    
    
    dt = t_scale*sim['dt']
    Ks = sim.Ks

    Ks = np.diff(sim.I, axis = 0, append=0)/(sim['dt'])
    Ks[-1] = sim.Ks

    vo = np.array([sim.vc[to[i], int(xo[i]/dx), min(int(yo[i]/dx), Ny - 1)] for i in range(len(xo))])
    uo = np.array([sim.uc[to[i], int(xo[i]/dx), min(int(yo[i]/dx), Ny - 1)] for i in range(len(xo))])    
    ho = np.array([sim.hc[to[i], int(xo[i]/dx), min(int(yo[i]/dx), Ny - 1)] for i in range(len(xo))])
    
    ## Add option for impermeable soils here...
    Ko = np.array([Ks[to[i], int(xo[i]/dx), min(int(yo[i]/dx), Ny- 1)] for i in range(len(xo))])

    veg = np.array([sim.veg[int(xo[i]/dx), min(int(yo[i]/dx), Ny- 1)] for i in range(len(xo))])

    rain_input = sim.p*dt/3.6e5*(np.array(to) <= sim.tr*60/dt)
    p_infl = np.random.rand(len(xo)) > Ko*dt/( h_min + ho + rain_input)
    # p_infl[ho < h_min] = 0
    # p_infl = np.random.rand(len(xo)) > Ko*dt/(h_min + ho)

    xnew = xo + dt*np.array(uo)
    ynew = yo + dt*np.array(vo) 
    xo = xnew
    yo = ynew

    if periodic == 0:
        yo[yo >= sim.l] = sim.l  - dx*0.01
    else:
        loop[yo >= sim.l] = loop[yo >= sim.l] + sim.l         
        yo[yo >= sim.l] = yo[yo >= sim.l] - sim.l
    
    #  needs fixing if slope curvature! 
    xo[xo>=sim.L] = sim.L - dx
    #xo[xo>=sim.L] =  xo[xo>=sim.L] - 2*dt*np.array(uo[xo>=sim.L])

    xo[xo < 0] =  dx
    #xo[xo < 0] =  xo[xo < 0] - 2*dt*np.array(uo[xo < 0])
    
    label  = label[p_infl]  
    ho  = ho[p_infl]  
    yo  = yo[p_infl]  
    xo  = xo[p_infl]  
    vo  = vo[p_infl]  
    uo  = uo[p_infl]  
    to  = to[p_infl] 
    veg = veg[p_infl] 
    loop  = loop[p_infl]      

    for i in range(len(xo)):
        l = label[i]
        positions[l]["xo"].append(xo[i])
        positions[l]["yo"].append(yo[i])
        positions[l]["ho"].append(ho[i])
        positions[l]["uo"].append(uo[i])
        positions[l]["vo"].append(vo[i])    
        positions[l]["loop"].append(loop[i])            
        positions[l]["t"].append(sim.t[to[i]])
        positions[l]["veg"].append(veg[i])

    to = to + t_scale
    to[to >=len(sim.t)] = len(sim.t) - 1

    return label, ho, vo, uo, xo, yo, loop, positions, to

def describe_positions(positions, sim):
    """

    """
    recap = pd.Series(name = sim.name)

    positions['y_f'] = [y[-1] for y in positions['yo']]
    positions['veg_f'] = [veg[-1] for veg in positions['veg']]
    positions['veg_bin'] = [int(np.sum(positions.iloc[i].veg) > 0) for i in range(len(positions))]

    positions['y_i'] = [y[0] for y in positions['yo']]
    positions['x_f'] = [y[-1] for y in positions['xo']]
    positions['x_i'] = [y[0] for y in positions['xo']]    

    positions['t_f'] = [t[-1] for t in positions['t']]
    positions['t_i'] = [t[0] for t in positions['t']]
    positions['std_x'] = [np.std(xo) for xo in positions['xo']]    
    positions['range_x'] = [np.max(xo) - np.min(xo) for xo in positions['xo']]    
    positions['loop'] = [l[-1] for l in positions['loop']]
    positions['t_esc'] = np.nan

    if sim.veg_type == 'dot':
        positions['up_out'] = 0
        positions['captured'] = 0
        
        x, y = np.where(sim.veg == 1)
        
        query = "((x_i > {0} and x_i < {0} + 1) or (x_i < {1} and x_i > {1} - 1) ) and y_i < {2}".format(
            x.max()*sim.dx, x.min()*sim.dx, y.min()*sim.dx - sim.r)
        
        inds= positions.query(query).index
        
        positions.loc[inds, 'up_out'] = 1

        inds = positions.query("(up_out == 1) and (veg_bin== 1)").index
        positions.loc[inds, 'captured'] = 1

        positions['up_in'] = 0
        x, y = np.where(sim.veg == 1)
        query = "x_i <= {0} and x_i >= {1} and y_i < {2}".format(
            x.max()*sim.dx, x.min()*sim.dx, y.min()*sim.dx - sim.r )
        inds= positions.query(query).index

        positions.loc[inds, 'up_in'] = 1

        positions[ 'diverted'] = 0
        inds = positions.query("(up_in == 1) and (veg_bin == 0) and y_f > {0}".format(y.max()*sim.dx)).index
        positions.loc[inds, 'diverted'] = 1

    for t in positions.index:

        pos = positions.loc[t]
        dx = sim.dx
        xo = pos.xo
        yo = pos.yo

        escape_time = np.where(np.array(yo) >= sim.l - dx*0.1)[0]
        if len(escape_time) > 0:
            positions.at[t, 't_esc'] = (pos.t[escape_time[0]])

        # residence in vegetated sites
        residence = np.array([sim.veg[int(xo[i]/sim.dx), min(int(yo[i]/dx), sim.Nycell- 1)]
                 for i in range(len(xo))]).mean()
        positions.at[t,"residence"] = residence
        
        # get arc length
        apts = np.array([xo, yo]).T 
        lengths = np.sqrt(np.sum(np.diff(apts, axis=0)**2, axis=1)) # Length between corners
        total_length = np.sum(lengths)
        positions.at[t, 'length'] = total_length/sim.l
        
    escape =  positions['y_f'] >= sim.l - sim.dx*0.1
    if sim.Tbound == 4:
        escape = escape*0
    positions['escape'] = escape

    escape_in_domain = positions.query("(t_f == {0}) and (y_f < {1})".format(sim.t.max(), sim.l-0.1*sim.dx))
    positions["escape_in_domain"]  = 0
    positions.loc[escape_in_domain.index, "escape_in_domain"]  = 1
    if sim.Tbound != 4:
        positions['time_to_escape'] = (positions['t_esc'][escape] - positions['t_i'][escape])/60
    else:
        positions['time_to_escape'] = np.nan
    
    positions['time_to_infiltrate'] = (positions['t_f'][escape == 0] - positions['t_i'][escape==0])/60

    distance = (positions['y_f'] - positions['y_i'] + positions.loop)/sim.l
    positions["distance"] = distance
 
    recap['tracer_C'] =  np.nanmean(positions.escape)
    recap['tracer_C_inc'] = (len(escape_in_domain) + len(escape))/len(positions)
    
    recap['tracer_IF'] = 1 - np.nanmean(positions.escape)
    recap['tracer_IF_inc'] = 1 - (len(escape_in_domain) + len(escape))/len(positions)     
    
    tracer_err = sim['infl_frac'] - (1 - np.nanmean(positions.escape))
    tracer_err_inc = sim['infl_frac'] - (1 - (len(escape_in_domain) + len(escape))/len(positions))
        
    tracer_source = distance[positions.escape == 1]
    
    recap["path_source"]  = np.mean(tracer_source)
    recap["path_source_med"]  = np.median(tracer_source)

    recap["path_avg"] = distance.mean()

    recap["path_IF"] = distance[positions.escape == 0].mean()
    recap["path_IF_med"] = np.median(distance[positions.escape == 0])

    recap["curve"] = np.mean(positions['std_x'])
    recap["eta"] = np.mean(positions['range_x'])        

    recap['length_source'] = positions.length[positions.escape == 1].mean()
    recap['length_IF'] = positions.length[positions.escape == 0].mean()
    recap['max_length'] = np.max(positions.length)

    recap['tracer_err'] = sim['infl_frac'] - recap['tracer_IF']
    recap['tracer_err_inc'] = sim['infl_frac'] - recap['tracer_IF_inc']
    recap['abs_err'] = np.abs(recap['tracer_err'])
    recap['divert'] = recap.path_source - recap.path_IF

    recap['vcount'] = np.sum([np.sum(positions.iloc[i].veg) > 0 for i in range(len(positions))])
    
    if sim.veg_type == 'dot':    
        l = []
        for i in range(20):               
            s = positions.sample(len(positions)//2, random_state= i)
            l.append( s['captured'].mean())
        
        recap['captured_mean'] = positions['captured'].mean()
        recap['captured_std'] = np.std(l)

        l = []
        for i in range(20):               
            s = positions.sample(len(positions)//2, random_state= i)
            l.append( s['diverted'].mean())
        
        recap['diverted_mean'] = positions['diverted'].mean()
        recap['diverted_std'] = np.std(l)

    return  positions, recap

def integrate_positions(sim, N = 1000, initialize = False, t_min = 0, y_min = 0, 
            h_min = 0, t_scale = 1, periodic = 0, seed = None):

    label, ho, vo, uo, xo, yo, loop, positions, to = tracer_init(sim = sim, N = N, initialize = initialize, 
        t_min = t_min, y_min = y_min, seed = seed)

    # start = int(60/dt) + 1
    for t in range(t_min + 1, len(sim.t)+1):        
        label, ho, vo, uo, xo, yo, loop, positions, to = step(label, ho, vo, uo, xo, yo, loop, 
                                                                positions, to, sim, h_min, t_scale, periodic)
    
    positions = pd.DataFrame(positions).T     

    positions, recap = describe_positions(positions, sim)
    
    return positions, recap

def run_tracers(summary, N = 1000, initialize = 0, t_min = 1,  y_min = 1, h_min = 5e-4, t_scale = 1):
    """
    adding tracers
    """
    cols = [ 'y_f', 'y_i', 'x_f', 'x_i', 't_f', 't_i', 'veg_bin',
            'std_x',  't_esc', 'residence', 'length','escape', 'escape_in_domain', 
            'time_to_escape', 'time_to_infiltrate', 'distance']
    
    tracers = pd.DataFrame()
    
    for key in summary.index: 

        sim = summary.loc[key]  

        if sim.Tbound == 4:
            periodic = 1
        else: 
            periodic = 0

        positions, recap = integrate_positions(sim, N = N, initialize = initialize, 
            t_min = t_min,  y_min = y_min, h_min = h_min, 
            t_scale = t_scale, periodic = periodic)

        recap['N'] = N
        recap['t_min'] = t_min
        recap['h_min'] = h_min
        recap['y_min'] = y_min
        recap['periodic'] = periodic
        recap['t_scale'] = t_scale

        for fld in cols: 
            recap.loc[fld] = np.array(positions[fld])
        

        recap['key'] = key
        tracers = pd.concat([tracers, recap], axis = 1)
        print (key, np.round(recap['tracer_err'] , 3))

    return tracers.T

def add_tracers(summary, tracers):
    
    summary = pd.concat([summary, tracers], axis = 1)
    return summary


def delete_tracers(summary, tracers):

    tracer_cols = tracers.columns
    
    summary = summary.drop(tracer_cols, axis =  1)

    return summary



def error_dist(sim, count, t_min = 1, y_min = 1, h_min = 5e-4, N = 1000, t_scale = 1,
               initialize = "", periodic = 0):
    """
    See if the worst case was a fluke
    """
    dist = pd.DataFrame()
    
    for i in range(count):
        np.random.seed(i)
        positions, recap = integrate_positions(sim, N = N, t_min = t_min,  h_min = h_min,
            y_min  = y_min, t_scale = t_scale, periodic = periodic, initialize = initialize, seed = i)
        # dist = dist.append(recap)
        dist = pd.concat([dist, recap], axis = 0)

    return dist


def add_tracer_lengthscales(summary):
    """
    U_t : mean tracer velocity
    U_te : mean velocity of tracers that escape
    x_t : mean displacement of tracers
    x_te : mean displacement of tracers that escape
    t_t : mean transport time of tracers
    t_te : mean transport time of tracers that escape

    U_b : mean Eulerian velocity in bare soil areas
    U_v : mean Eulerian velocity in vegetated areas
    U_avg : mean Eulerian velocity 

    t_HD : hydrograph duration
    """
    for key in summary.index:

        sim = summary.loc[key]
        #sim.t_esc[sim.t_esc.isna()] = sim.t_f

        inds = [np.where(np.isnan(sim['t_esc']))[0]]
        np.array(sim['t_esc'])[inds] = np.array(sim['t_f'])[inds]

        U_t = (sim.length*sim.l/(np.minimum(sim.t_f, sim.t_esc) - sim.t_i))  # m    
        # U_t = U_t[U_t<100]

        summary.at[key, 'U_t'] = np.mean(U_t[U_t<100])
        summary.at[key, 'U_te'] = np.mean(U_t[sim.escape == 1])

        summary.at[key, 'x_t'] = np.nanmean(sim.distance)
        summary.at[key, 'x_te'] = np.nanmean(sim.distance[sim.escape == 1])

        summary.at[key, 't_t'] = (np.minimum(sim.t_f, sim.t_esc) - sim.t_i).mean()  # m   
        summary.at[key, 't_te'] = (np.minimum(sim.t_f, sim.t_esc) - sim.t_i)[sim.escape == 1].mean()  # m       

        summary.at[key, 'U_b'] = sim.vc[1:sim.i_tr, sim.veg == 0].mean()
        summary.at[key, 'U_v'] = sim.vc[1:sim.i_tr, sim.veg == 1].mean()
        summary.at[key, 'U_avg'] = sim.vc[1:].mean()    

        summary.at[key, 'U_avg_l'] = sim.vc[1:, :, -1].mean()    

        
        if len(sim.t) > len(sim.hydro):
            t = sim.t[:-1]
        else:
            t = sim.t
        try:
            summary.at[key, 't_HD'] = t[sim.hydro > 2e-2][-1]      
        except:
            summary.at[key, 't_HD'] = np.nan
        summary['t_fall'] = summary['t_HD'] - summary['tr']
        
    return summary    
