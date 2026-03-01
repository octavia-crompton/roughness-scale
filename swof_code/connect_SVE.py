"""

All units in m, m/s, s, etc.
"""
import numpy as np
from plot_SWOF import *
from gw_functions import *
from scipy.optimize import minimize
import time
from scipy import interpolate

from sklearn.model_selection import train_test_split
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import LinearRegression   

def get_GW_param(sim,  a = 2/3, threshold = 2e-7, ntstep = 100, use_IF = 1):
    """
    Extract GW input parameters from an SVE simulation `sim`

    Parameters:
    ----------
    sim : dict
        simulation dictionary from SVE model

    Returns:
    -------
    param : dict
        parameter dictionary needed to run GW model, with components:

        ksatV : Ks (m/s)
        rain : rain intensity (m/s)
        So : slope (m/m)
        Kr : generalized roughness parameter
             So^eta/alpha
        a, m : flow regime exponents  (m = a+1)
            convention used in GW
        t_rain : duration (s)
        L : hillslope length (m)
        Ao : sorptivity (m/s**.5)
        t_pond : time of ponding (s)

    """
    if use_IF == 1:

        ksatV = sim.infl_frac*sim.p/3.6e5
    
    elif sim.Ks.mean() > 0:
        ksatV = sim.Ks.mean()
    else:
        ksatV = 0

    rain = sim.p/3.6e5
    So = sim.So

    #   alpha_eff = get_alpha_avg(sim)
    #   Kr = sim.So ** (1/2.)/get_alpha_avg(sim)
    
    # fit a 
    x = sim.hc.mean(1).mean(1)
    y = sim.qy.mean(1).mean(1)
    inds = y > 1e-7

    x = np.log(x[inds])
    y = np.log(y[inds])
    Kr = np.exp(y - (a+1)*x).mean()
    alpha_eff = sim.So**(1/2.)/Kr
    
    t_rain = sim.t_rain
    L = sim.l

    dt = 1.0

    t_pond =  SVE_t_pond(sim, threshold)

    # Ao = np.sqrt(2 * t_pond * rain * (rain - ksatV) ** 2 /
    #                  (rain - ksatV / 2.))

    # delta_theta = sim.theta_s - sim.theta_i

    param = {'ksatV': ksatV,
             'rain': rain,
             'So': So,
             'Kr': Kr,
             'a': a,
             'm': a + 1,
             't_rain': t_rain,
             'L': L,
             'dt': dt,
             'ntstep': ntstep,
             'alpha_eff' : alpha_eff,
             'Ao': 0,
             't_pond': t_pond,
             }

    return param

def SVE_t_pond(sim, threshold = 1e-3):
    """
    Extract time of ponding from an SVE simulation
    """
    Vol_of_tot = sim.Vol_of_tot
    ponded_inds = np.where(Vol_of_tot > threshold)[0]
    
    if len(ponded_inds) > 0:
        t_pond = sim.time[np.min(ponded_inds[0])]
    else:
        t_pond = np.inf

    if t_pond == sim["dt"]:
        t_pond = 0

    return t_pond

#### Calibration plots

def plot_GW_SVE( sim, res):
    """
    
    """
    fig, ax = plt.subplots(1, 1, figsize = (7, 3))

    ax.plot(res.t/60, res.q_GW/res.L*3.6e5, 'k--', label = "GW")

    plt.xlabel("min")
    ax.plot(res.t/60, res.q_SVE/res.L*3.6e5, label = "SVE")
    
    ax.legend()
    ax.set_ylabel("runoff (cm/hr)")
    ax.set_xlabel("time (min)")
    
    t = plt.suptitle( sim.name, fontsize = 12)
    t.set_y(-0.08)
    t.set_x(0.3)

    return fig


def plot_GW_SVE(sim, res):
    """
    
    """
    fig, ax = plt.subplots(1, 1, figsize = (7, 3))

    ax.plot(res.t/60, res.q_GW/res.L*3.6e5, 'k--', label = "GW")

    plt.xlabel("min")
    ax.plot(res.t/60, res.q_SVE/res.L*3.6e5, label = "SVE")
    
    ax.legend()
    ax.set_ylabel("runoff (cm/hr)")
    ax.set_xlabel("time (min)")
    
    t = plt.suptitle( sim.name, fontsize = 12)
    t.set_y(-0.08)
    t.set_x(0.3)

    return fig

def plot_uncalibrated(sim):
    """
    Plot GW results without calibration
    """
    param = get_GW_param(sim)
    try:
        res = compute_GW(sim, param)
    except:
        param["ksatV"] = sim.Ks.mean()/2 
        res = compute_GW(sim, param)    
    q_nrmse, Q_nrmse, IF_percent_diff =  get_GW_error(sim, res)
    print ("Uncalibrated error is {0:.2f}".format(Q_nrmse))

    res = compute_GW(sim, param)
    q_nrmse, Q_nrmse, IF_percent_diff =  get_GW_error(sim, res)

    plt.figure(figsize = (5,3))
    plt.plot(res.t/60, res.q_GW*3.6e5/sim.l, label = "GW")
    plt.plot(res.t/60, res.q_SVE*3.6e5/sim.l, label = "SVE")
    plt.legend()
    plt.xlabel("time (min)")
    plt.ylabel("q (cm/hr)")    
    plt.xlim(0, sim.tr*3)

def compare_opti_GW_SVE(  res, res2, label1, label2):
    """
    Compare two optimizations
    """
    fig, ax = plt.subplots(1, 1, figsize = (7, 3))
    
    ax.plot(res.t_avg/60, res.q_SVE/res.l*3.6e5, label = "SVE")
    ax.plot(res.t_o/60, res.q_o/res.l*3.6e5, 'g--', label = label1 + " {0:.2f}".format(res.Q_err_o))
    ax.plot(res2.t_o/60, res2.q_o/res.l*3.6e5, '--', label = label2 + " {0:.2f}".format(res2.Q_err_o))

    ax.legend()
    ax.set_ylabel("runoff (cm/hr)")
    ax.set_xlabel("time (min)")
    
    t = plt.suptitle( res.name, fontsize = 12)
    t.set_y(-0.08)
    t.set_x(0.2)
    ax.set_title("Q(NRMSE) = {0:.2f}".format(res.Q_err_o))
    
    return fig

def plot_opti_GW_SVE( res):
    """
    Optimized GW check
    """
    fig, ax = plt.subplots(1, 1, figsize = (7, 3))
    
    ax.plot(res.t_avg/60, res.q_SVE/res.l*3.6e5, label = "SVE")
    ax.plot(res.t_avg/60, res.q_avg/res.l*3.6e5, 'k--', label = "GW")
    ax.plot(res.t_o/60, res.q_o/res.l*3.6e5, 'g--', label = "optimized GW")

    ax.legend()
    ax.set_ylabel("runoff (cm/hr)")
    ax.set_xlabel("time (min)")
    
    t = plt.suptitle( res.name, fontsize = 12)
    t.set_y(-0.08)
    t.set_x(0.2)
    ax.set_title("Q(NRMSE) = {0:.2f}".format(res.Q_err_o))
    
    qmax = np.max([res.q_SVE.max(), res.q_avg.max(), res.q_o.max()])
    plt.ylim(0, qmax*3.6e5/res.l*1.01)

############################# Calibration functions ####################################

from scipy.stats.mstats import gmean
import numpy as np

# Functions for estimating 'effective' hillslope parameters

roughness_names = {"alpha_ss" :  r"$\langle n \rangle_{SS}$",
                   "alpha_avg" :  r"$\langle n \rangle$",
                   "alpha_gmean" : r"$\langle n \rangle_{g}$",                   
                  }

def get_alpha_ss(sim):
    """
    estimate effective roughness using a sum of squares approach 
    """
    return np.sqrt(sim.alpha_v**2*sim.fV**2+ sim.alpha_b**2*(1-sim.fV)**2)


def get_alpha_avg(sim):
    """
    estimate effective roughness using a weighted average
    """
    return  np.mean(get_alpha(sim).ravel())

def get_alpha_gmean(sim):
    """
    estimate effective roughness using a geometric average
    """
    return  gmean(get_alpha(sim).ravel())

def get_alpha(sim):
    """
    get the simulation alpha field
    """
    alpha = sim.veg.copy().astype(float)
    alpha[sim.veg < 1] = sim.alpha_b
    alpha[sim.veg == 1] = sim.alpha_v
    
    return alpha

def get_Ks(sim):
    """
    get the simulation Ks field 
    """
    Ks = sim.veg.copy().astype(float)
    Ks[sim.veg < 1] = sim.Ks_b/3.6e5
    Ks[sim.veg == 1] = sim.Ks_v/3.6e5
    
    return Ks

def update_alpha(param, sim, avg_func):
    """
    Update GW parameters - alpha and Ksat - with the weighted average
    
    resistors in series
    """
    param["alpha_eff"] = avg_func(sim)
    param["Kr"] = sim.So**0.5/param["alpha_eff"]
    
    return param    


def update_Ks(param, sim, avg_func = np.mean):
    """
    Update GW parameters - alpha and Ksat - with the weighted average
    
    resistors in series
    """
    param["ksatV"] = avg_func(get_Ks(sim).ravel())
    param["Kr"] = sim.So**0.5/param["alpha_eff"]
    
    return param    


def update_gmean(param, sim):
    """
    Update GW parameters - alpha - with the geometric mean
    """
    param["alpha_eff"] = gmean(get_alpha(sim).ravel())
    param["Kr"] = sim.So**0.5/param["alpha_eff"]
    return param


def update_ss_mean(param, sim):
    """
    Update GW parameter - alpha  - with the sum of squares mean  
    """
    param["alpha_eff"] = np.sqrt(sim.alpha_v**2*sim.fV + sim.alpha_b**2*(1-sim.fV))
    param["Kr"] = sim.So**0.5/param["alpha_eff"]
    
    return param    


# Calibration / inverse model functions

def compute_GW(sim, param):
    """
    """
    res = Comparison_function2(param)
    res = pd.Series(res)

    from scipy import interpolate
    if sim.l > 10:
        dt = 15
    else:
        dt = 1
        
    d = int(dt/np.diff(sim.time)[0])  
    t_sve = sim.time[::d]  
    Q_sve = sim["Vol_bound_tot"][::d]

    tmax = max(res.t[-1], t_sve[-1])
    tmin = max(res.t[0], t_sve[0])

    res["tmax"] = tmax
    res["tmin"] = tmin

    if res.t[-1] < tmax: 

        times = np.arange(np.round(res.t[-1])+1, tmax+2)

        res["t"] = np.hstack((res["t"], times))
        res["q"] = np.hstack((res["q"], 0*times))    

    if t_sve[-1] < tmax: 

        times = np.arange(np.round(t_sve[-1])+1, tmax+2)

        t_sve = np.hstack((t_sve, times))
        Q_sve = np.hstack((Q_sve, Q_sve[-1]*np.ones_like(times)))    
    
    f = interpolate.interp1d(res.t, res.q)
    t_GW = np.arange(tmin, tmax)
    q_GW = f(t_GW)      #  m2/s

    res["t"] = t_GW[1:]
    res["q_GW"] = q_GW[1:]


    f = interpolate.interp1d(t_sve, Q_sve/sim.L)
    Q_sve = f(t_GW)      
    res["q_SVE"] = np.diff(Q_sve)
    res["Q_GW"] = np.cumsum(q_GW)  
    res["Q_SVE"] = Q_sve

    return res


def compare_GW_SVE(sim, res):
    """
    Interpolates and computes error metrics
    """
    d = int(10/np.diff(sim.time)[0])  
    t_sve = sim.time[::d]  
    Q_sve = sim["Vol_bound_tot"][::d]

    tmax = max(res.t[-1], t_sve[-1])
    tmin = max(res.t[0], t_sve[0])

    res["tmax"] = tmax
    res["tmin"] = tmin

    if res.t[-1] < tmax: 
        times = np.arange(np.round(res.t[-1])+1, tmax+2)

        res["t"] = np.hstack((res["t"], times))
        res["q"] = np.hstack((res["q"], 0*times))    

    if t_sve[-1] < tmax: 
        
        times = np.arange(np.round(res.t[-1])+1, tmax+2)
        t_sve = np.hstack((t_sve, times))
        Q_sve = np.hstack((Q_sve, Q_sve[-1]*np.ones_like(times)))    

    f = interpolate.interp1d(res.t, res.q)
    t_GW = np.arange(tmin, tmax)
    q_GW = f(t_GW) #  m2/s
    
    f = interpolate.interp1d(t_sve, Q_sve/sim.L)
    Q_sve = f(t_GW)      
    q_sve = np.diff(Q_sve, prepend = 0)  # m2/s

    Q_GW = np.cumsum(q_GW)  

    inds =  np.where((q_GW > 0.2/3.6e5) &  (q_sve >  0.2/3.6e5))[0]

    q_nrmse = np.zeros_like(q_sve)
    q_nrmse[inds] = (q_sve[inds] - q_GW[inds])**2
    q_nrmse = np.sqrt(q_nrmse[inds].mean())/(Q_sve.max())*100

    Q_nrmse = np.zeros_like(q_sve).astype(float)    
    Q_nrmse[inds] = (Q_sve[inds] - Q_GW[inds])**2
    Q_nrmse = np.sqrt(Q_nrmse[inds].mean())/(Q_sve.max())*100    

    if sim.infl_frac == 0:
        IF_diff = 0
    else:
        IF_diff = (sim.infl_frac - res.infl_frac)
    
    res['q_nrmse'] = q_nrmse
    res['Q_nrmse'] = Q_nrmse
    res['IF_diff'] = IF_diff
    res['t_GW'] = t_GW
    res['Q_GW'] = Q_GW
    res['Q_sve'] = Q_sve
    res['Q_GW'] = Q_GW
    res['Q_sve'] = Q_sve    
    
    return res
    

def get_GW_error(sim, res ):
    """
    Computes error metrics only
    """
    q_SVE = res.q_SVE
    q_GW = res.q_GW    
    
    Q_SVE = res.Q_SVE
    Q_GW = res.Q_GW    
    
    inds =  np.where((q_GW > 0.2/3.6e5) &  (q_SVE >  0.2/3.6e5))[0]
    
    q_nrmse = np.zeros_like(q_SVE)
    q_nrmse[inds] = (q_SVE[inds] - q_GW[inds])**2
    q_nrmse = np.sqrt(q_nrmse[inds].mean())/(Q_SVE.max())*100

    Q_nrmse = np.zeros_like(q_SVE).astype(float)    
    Q_nrmse[inds] = (Q_SVE[inds] - Q_GW[inds])**2
    Q_nrmse = np.sqrt(Q_nrmse[inds].mean())/(Q_SVE.max())*100    

    if sim.infl_frac == 0:
        IF_diff = 0
    else:
        IF_diff = (sim.infl_frac - res.infl_frac)

    return q_nrmse, Q_nrmse, IF_diff


def objective(p, param, sim):
    """
    Penalizes solutions with no domain 2 output
    """
    param["ksatV"] = p[0]
    param["alpha_eff"] = p[1]
    param["Kr"] = sim.So**0.5/param["alpha_eff"]

    res = compute_GW(sim, param)
    if len(res['q_d2']) == 0:
        return 100
    q_nrmse, Q_nrmse, IF_diff = get_GW_error(sim, res)
    return abs(Q_nrmse)    


def objective_imp(p, param, sim):
    """
    Penalizes solutions with no domain 2 output
    """
    param["alpha_eff"] = p[0]
    param["Kr"] = sim.So**0.5/param["alpha_eff"]

    res = compute_GW(sim, param)
    if len(res['q_d2']) == 0:
        return 100    
    q_nrmse, Q_nrmse, IF_diff = get_GW_error(sim, res)
    return abs(Q_nrmse)    



def objective_a(p, param, sim):
    """
    Penalizes solutions with no domain 2 output
    also adjusts the exponent a
    """
    param["ksatV"] = p[0]
    param["alpha_eff"] = p[1]
    param["Kr"] = sim.So**0.5/param["alpha_eff"]
    
    param["a"] = p[2]
    param["m"] = p[2] +1   

    res = compute_GW(sim, param)

    if len(res['q_d2']) ==0:
        return 100
    
    q_nrmse, Q_nrmse, IF_diff = get_GW_error(sim, res)
    return abs(Q_nrmse)    


def optimized_GW(result, param, sim):
    """
    Takes result of the nelder mead optimization and returns
    """
    opti_NM = result["final_simplex"][0][1]    
    opti_param = param.copy()
    
    opti_param["alpha_eff"] = opti_NM[1]
    opti_param["Kr"] = sim.So**0.5/opti_NM[1]
    opti_param["ksatV"] =  opti_NM[0]
    
    if len(opti_NM) > 2:
        opti_param["a"] =  opti_NM[2]
        opti_param["m"] =  opti_NM[2]+1

    opti_res = Comparison_function2(opti_param)
    opti_res = pd.Series(opti_res)
    opti_res["error"] = result.fun
    opti_res["success"] = 1 if result.message == "Optimization terminated successfully." else 0

    d = int(10/np.diff(sim.time)[0])  
    t_sve = sim.time[::d]  
    Q_sve = sim["Vol_bound_tot"][::d]

    tmax = max(opti_res.t[-1], t_sve[-1])
    tmin = max(opti_res.t[0], t_sve[0])

    opti_res["tmax"] = tmax
    opti_res["tmin"] = tmin

    if opti_res.t[-1] < tmax: 
        times = np.arange(np.round(opti_res.t[-1])+1, tmax+2)

        opti_res["t"] = np.hstack((opti_res["t"], times))
        opti_res["q"] = np.hstack((opti_res["q"], 0*times))    
    
    if t_sve[-1] < tmax: 
        times = np.arange(np.round(opti_res.t[-1])+1, tmax+2)

        t_sve = np.hstack((t_sve, times))
        Q_sve = np.hstack((Q_sve, Q_sve[-1]*np.ones_like(times)))    

    f = interpolate.interp1d(opti_res.t, opti_res.q)
    t_GW = np.arange(tmin, tmax)
    q_GW = f(t_GW)  #  m2/s

    opti_res["t"] = t_GW[1:]
    opti_res["q_GW"] = q_GW[1:]
    
    f = interpolate.interp1d(t_sve, Q_sve/sim.L)
    Q_sve = f(t_GW)      
    opti_res["q_SVE"] = np.diff(Q_sve)  # m2/s

    opti_res["Q_GW"] = np.cumsum(q_GW)  
    opti_res["Q_SVE"] = Q_sve    


    return opti_res

def optimized_GW_imp(result, param, sim):
    """
    Takes result of the nelder mead optimization and returns
    """
    opti_NM = result["final_simplex"][0][0]

    opti_param = param.copy()
    opti_param["alpha_eff"] = opti_NM[0]
    
    if len(opti_NM) > 1:
        opti_param["a"] =  opti_NM[1]
        opti_param["m"] =  opti_NM[1]+1
        
    opti_param["Kr"] = sim.So**0.5/opti_NM[0]

    opti_res = Comparison_function2(opti_param)
    opti_res = pd.Series(opti_res)
    opti_res["error"] = result.fun
    opti_res["success"] = 1 if result.message == "Optimization terminated successfully." else 0
    
    if sim.l > 10:
        time_res = 15
    else:
        time_res = 1
    d = int(time_res/np.diff(sim.time)[0])  
    t_sve = sim.time[::d]  
    Q_sve = sim["Vol_bound_tot"][::d]

    tmax = max(opti_res.t[-1], t_sve[-1])
    tmin = max(opti_res.t[0], t_sve[0])

    opti_res["tmax"] = tmax
    opti_res["tmin"] = tmin

    if opti_res.t[-1] < tmax: 
        times = np.arange(np.round(opti_res.t[-1])+1, tmax+2)

        opti_res["t"] = np.hstack((opti_res["t"], times))
        opti_res["q"] = np.hstack((opti_res["q"], 0*times))    
    
    if t_sve[-1] < tmax: 
        times = np.arange(np.round(opti_res.t[-1])+1, tmax+2)

        t_sve = np.hstack((t_sve, times))
        Q_sve = np.hstack((Q_sve, 0*times))    
    

    f = interpolate.interp1d(opti_res.t, opti_res.q)
    t_GW = np.arange(tmin, tmax)
    q_GW = f(t_GW)      #  m2/s

    opti_res["t"] = t_GW[1:]
    opti_res["q_GW"] = q_GW[1:]
    
    f = interpolate.interp1d(t_sve, Q_sve/sim.L)
    Q_sve = f(t_GW)      
    opti_res["q_SVE"] = np.diff(Q_sve)

    opti_res["Q_GW"] = np.cumsum(q_GW)  
    opti_res["Q_SVE"] = Q_sve    



    return opti_res



def calibrate_to_GW_a(summary, maxiter = 60, a = 2/3., use_a = 1, use_IF = 0):
    """
    Calibrate SVE simulations to GW solution, using nelder mead
    exponent a is prescribed 

    """
    calibrated = {}
    initial_start = time.time()
    for ind, key in enumerate(summary.index):

        print (ind, key,)
        start = time.time()
        sim = summary.loc[key]

        if np.max(sim.Vol_bound_tot/sim.Vol_rain_tot) < 0.02:
            continue
        else:
            calibrated[key] = {}

        for var in ['p', 'So', 'Ks_v', 'Ks_b', 'fV']:
            calibrated[key][var] = sim[var]
        calibrated[key]['infl_frac_sve'] = sim['infl_frac']           

        if use_a == 1:
            a = get_exp(sim)

        param = get_GW_param(sim, a, threshold = 2e-7, ntstep = 100, use_IF = use_IF)


        try:
            res = compute_GW(sim, param)


        except IndexError:
            try:
                param["ksatV"] = param["ksatV"]/2
                res = compute_GW(sim, param)
        
            except:
                print ("no initial solution")
                pass

        for var in ['Ao', 't_pond', 'q_SVE'] :   
            calibrated[key][var] = res[var]
        
        calibrated[key]['l'] = res["L"]
        
        for var in  ['t', 'infl_frac'] :  
            calibrated[key][var + '_avg'] = res[var]   
        
        calibrated[key]['q_avg'] = res['q_GW']           
        calibrated[key]['q_SVE'] = res['q_SVE']     
        calibrated[key]['Q_avg'] = res['Q_GW']           
        calibrated[key]['Q_SVE'] = res['Q_SVE']         


        calibrated[key]['alpha_avg'] = get_alpha_avg(sim) 
        calibrated[key]['Ks_avg'] = np.mean(get_Ks(sim).ravel())*3.6e5    

        calibrated[key]['alpha_gm'] = gmean(get_alpha(sim).ravel())
        calibrated[key]['Ks_gm'] = gmean(get_Ks(sim).ravel())*3.6e5

        q_nrmse, Q_nrmse, IF_percent_diff =  get_GW_error(sim, res)

        calibrated[key]['q_err_avg'] = np.round(q_nrmse, 2)
        calibrated[key]['Q_err_avg'] = np.round(Q_nrmse, 2)
        calibrated[key]['IF_err_avg']  = np.round(IF_percent_diff, 2)

        fatol = 0.01
        
        try:
            param["ksatV"] =  sim.Ks.mean() 
            result = minimize(objective, [param["ksatV"], param["alpha_eff"], param['a']], args = (param, sim),
                          method='nelder-mead',
                          options= {"maxiter": maxiter, "return_all" : False, "fatol" : fatol, "xatol" : 1, 
                                    "adaptive" : True}) 

        except:
            try:
                param["ksatV"] = sim.Ks.mean()/2 # gmean(get_Ks(sim).ravel())
                result = minimize(objective, [param["ksatV"], param["alpha_eff"], param['a']], args = (param, sim),
                          method='nelder-mead',
                          options= {"maxiter": maxiter, "return_all" : False, "fatol" : fatol, "xatol" : 1, 
                                    "adaptive" : True}) 

                res = compute_GW(sim, param)
        
            except:
                print ("not optimized")
                continue


        calibrated[key]['fatol'] = fatol
        calibrated[key]['maxiter'] = maxiter 

        calibrated[key]['nit'] =  result.nit    
        calibrated[key]["opti_err"] = result.final_simplex[-1][0]
        opti_res = optimized_GW(result, param, sim)    
        q_nrmse, Q_nrmse, IF_percent_diff =  get_GW_error(sim, opti_res)

        calibrated[key]['alpha_o'] = opti_res.alpha_eff
        calibrated[key]['Ks_o'] =     opti_res.ksatV*3.6e5
        calibrated[key]['a'] = opti_res.a

        calibrated[key]['q_o'] = opti_res['q_GW']
        calibrated[key]['Q_o'] = opti_res['Q_GW']
        calibrated[key]['t_o'] = opti_res['t']
        calibrated[key]['infl_frac_o'] = opti_res['infl_frac']   
        calibrated[key]['success'] = opti_res['success']       

        calibrated[key]['q_err_o'] = np.round(q_nrmse, 2)
        calibrated[key]['Q_err_o'] = np.round(Q_nrmse, 2)    
        calibrated[key]['IF_err_o']  = np.round(IF_percent_diff, 2)
        calibrated[key]['opti_time']  = time.time() - start
        
    calibrated = pd.DataFrame(calibrated).T        
  
    calibrated["IF_err_avg"] = calibrated.infl_frac_sve - calibrated.infl_frac_avg
    calibrated["IF_err_o"] = calibrated.infl_frac_sve - calibrated.infl_frac_o
    calibrated["IF_err_o_abs"] = np.abs(calibrated["IF_err_o"])

    return calibrated


def calibrate_to_GW(summary, maxiter = 60, a = 0., use_a = 1, use_IF = 0):
    """
    Calibrate SVE simulations to GW solution, using nelder mead
    """

    calibrated = {}
    initial_start = time.time()
    for ind, key in enumerate(summary.index):

        print (ind, key,)
        start = time.time()
        sim = summary.loc[key]

        if np.max(sim.Vol_bound_tot/sim.Vol_rain_tot) < 0.02:
            print ("skipping", key)
            continue
        else:
            calibrated[key] = {}

        for var in ['p', 'So', 'Ks_v', 'Ks_b', 'fV']:
            calibrated[key][var] = sim[var]
        
        calibrated[key]['infl_frac_sve'] = sim['infl_frac']           
        
        if use_a == 1:
            a = get_exp(sim)

        param = get_GW_param(sim, a, threshold = 2e-7, ntstep = 100, use_IF = use_IF)
        
        try:
            res = compute_GW(sim, param)

        except IndexError:
            try:
                param["ksatV"] = gmean(get_Ks(sim).ravel())
                res = compute_GW(sim, param)
        
            except:
                print ("no initial solution")
                pass

        for var in  ['Ao', 't_pond', 'q_SVE'] :   
            calibrated[key][var] = res[var]
        
        calibrated[key]['l'] = res["L"]
        
        for var in  ['t', 'infl_frac'] :  
            calibrated[key][var + '_avg'] = res[var]   

        calibrated[key]['a'] = a
        
        calibrated[key]['q_avg'] = res['q_GW']           
        calibrated[key]['q_SVE'] = res['q_SVE']     
        calibrated[key]['Q_avg'] = res['Q_GW']           
        calibrated[key]['Q_SVE'] = res['Q_SVE']         


        calibrated[key]['alpha_avg'] = get_alpha_avg(sim) 
        calibrated[key]['Ks_avg'] = np.mean(get_Ks(sim).ravel())*3.6e5    

        calibrated[key]['alpha_gm'] = gmean(get_alpha(sim).ravel())
        calibrated[key]['Ks_gm'] = gmean(get_Ks(sim).ravel())*3.6e5

        q_nrmse, Q_nrmse, IF_percent_diff =  get_GW_error(sim, res)

        calibrated[key]['q_err_avg'] = np.round(q_nrmse, 2)
        calibrated[key]['Q_err_avg'] = np.round(Q_nrmse, 2)
        calibrated[key]['IF_err_avg']  = np.round(IF_percent_diff, 2)

        fatol = 0.01
        
        try:
            param = get_GW_param(sim, a, threshold = 2e-7, ntstep = 100, use_IF = use_IF)
            result = minimize(objective, [param["ksatV"], param["alpha_eff"]], args = (param, sim),
                          method='nelder-mead',
                          options= {"maxiter": maxiter, "return_all" : False, "fatol" : fatol, "xatol" : 1, 
                                    "adaptive" : True}) 

        except:
            try:
                param["ksatV"] = sim.Ks.mean()/2 
                result = minimize(objective, [param["ksatV"], param["alpha_eff"]], args = (param, sim),
                          method='nelder-mead',
                          options= {"maxiter": maxiter, "return_all" : False, "fatol" : fatol, "xatol" : 1, 
                                    "adaptive" : True}) 

                res = compute_GW(sim, param)
        
            except:
                print ("not optimized")
                continue


        calibrated[key]['fatol'] = fatol
        calibrated[key]['maxiter'] = maxiter 

        calibrated[key]['nit'] =  result.nit    
        calibrated[key]["opti_err"] = result.final_simplex[-1][0]
        opti_res = optimized_GW(result, param, sim)    
        q_nrmse, Q_nrmse, IF_percent_diff =  get_GW_error(sim, opti_res)

        calibrated[key]['alpha_o'] = opti_res.alpha_eff
        calibrated[key]['Ks_o'] =     opti_res.ksatV*3.6e5

        calibrated[key]['q_o'] = opti_res['q_GW']
        calibrated[key]['Q_o'] = opti_res['Q_GW']
        calibrated[key]['t_o'] = opti_res['t']
        calibrated[key]['infl_frac_o'] = opti_res['infl_frac']   
        calibrated[key]['success'] = opti_res['success']       

        calibrated[key]['q_err_o'] = np.round(q_nrmse, 2)
        calibrated[key]['Q_err_o'] = np.round(Q_nrmse, 2)    
        calibrated[key]['IF_err_o']  = np.round(IF_percent_diff, 2)
        calibrated[key]['opti_time']  = time.time() - start
        
    calibrated = pd.DataFrame(calibrated).T        
  
    calibrated["IF_err_avg"] = calibrated.infl_frac_sve - calibrated.infl_frac_avg
    calibrated["IF_err_o"] = calibrated.infl_frac_sve - calibrated.infl_frac_o
    calibrated["IF_err_o_abs"] = np.abs(calibrated["IF_err_o"])

    return calibrated


def calibrate_to_GW_imp(summary, maxiter = 60, a = 2/3., use_a = 1):
    """
    """
    calibrated = {}
    initial_start = time.time()
    for ind, key in enumerate(summary.index):

        start = time.time()
        sim = summary.loc[key]

        if np.max(sim.Vol_bound_tot/sim.Vol_rain_tot) < 0.02:
            continue
        else:
            calibrated[key] = {}

        for var in ['p', 'So', 'Ks_v', 'Ks_b', 'fV']:
            calibrated[key][var] = sim[var]

        if use_a == 1:
            a = get_exp(sim)

        param = get_GW_param(sim, a, threshold = 2e-7, ntstep = 100, use_IF = 0)

        try:
            res = compute_GW(sim, param)

        except:
            print ("no initial solution")
            pass

        for var in  ['L',  'Ao', 't_pond', 'q_SVE'] :   
            calibrated[key][var] = res[var]

        for var in  ['t'] :  
            calibrated[key][var + '_avg'] = res[var]   

        calibrated[key]['q_avg'] = res['q_GW']           
        calibrated[key]['q_SVE'] = res['q_SVE']     
        calibrated[key]['Q_avg'] = res['Q_GW']           
        calibrated[key]['Q_SVE'] = res['Q_SVE']         

        ### add alpha eff
        calibrated[key]['alpha_avg'] = get_alpha_avg(sim) 
        calibrated[key]['alpha_gm'] = gmean(get_alpha(sim).ravel())

        q_nrmse, Q_nrmse, IF_percent_diff =  get_GW_error(sim, res)

        calibrated[key]['q_err_avg'] = np.round(q_nrmse, 2)
        calibrated[key]['Q_err_avg'] = np.round(Q_nrmse, 2)
        fatol = 0.01
        try:
            result = minimize(objective_imp, [param["alpha_eff"]], args = (param, sim),
                          method='nelder-mead',
                          options= {"maxiter": maxiter, "return_all" : False, "fatol" : fatol, "xatol" : 1, 
                                    "adaptive" : True}) 
    
        except:
                print ("not optimized")
                continue

        calibrated[key]['fatol'] = fatol
        calibrated[key]['maxiter'] = maxiter 

        calibrated[key]['nit'] =  result.nit    
        calibrated[key]["opti_err"] = result.final_simplex[-1][0]
        opti_res = optimized_GW_imp(result, param, sim)    
        q_nrmse, Q_nrmse, IF_percent_diff =  get_GW_error(sim, opti_res)

        calibrated[key]['alpha_o'] = opti_res.alpha_eff
        calibrated[key]['Ks_o'] =     opti_res.ksatV*3.6e5

        calibrated[key]['q_o'] = opti_res['q_GW']
        calibrated[key]['Q_o'] = opti_res['Q_GW']
        calibrated[key]['t_o'] = opti_res['t']
        calibrated[key]['success'] = opti_res['success']       

        calibrated[key]['q_err_o'] = np.round(q_nrmse, 2)
        calibrated[key]['Q_err_o'] = np.round(Q_nrmse, 2)    
        calibrated[key]['opti_time']  = time.time() - start
        
    calibrated = pd.DataFrame(calibrated).T        

    return calibrated


def update_summary(summary, compare):
    """
    Add results from compare / calibrated to summary

        summary = update_summary(summary, compare)    
    """
    summary["Ks_o"] = compare.Ks_o
    summary["alpha_o"] = compare.alpha_o    
    summary["Ks_vo"] = ((compare.Ks_o - summary.Ks_b*(1-summary.fV))/summary.fV)
    summary["alpha_vo"] = ((compare.alpha_o**2 - summary.alpha_b**2*(1-summary.fV))/summary.fV)**0.5

    summary["Q_err_o"] = compare["Q_err_o"]
    summary["IF_err_o"] = compare["IF_err_o"]
    summary["IF_err_o_abs"] = compare["IF_err_o_abs"]

    summary["Delta_s"] = summary["Ks_v"] - summary["Ks_vo"]
    summary["Delta_h"] = summary["alpha_v"] - summary["alpha_vo"]

    summary["psi_s"] =  summary["Ks_v"]/summary["Ks_vo"] 
    summary["psi_h"] =  summary["alpha_v"]/summary["alpha_vo"] 
    return summary



