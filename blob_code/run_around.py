import cv2
import numpy as np
import pandas as pd

def total_stress(sim):
    """
    rho*g*h*S_o
    """
    tau_t = 1000*9.8*sim.hc*sim.So

    if sim.fric == 6:
        tau_t[:, sim.veg == 1] = tau_t[:, sim.veg == 1]*(1-sim.phi_veg)
    
    return tau_t 

def kirst_f(sim):
    """
    """
    Re = sim.vc*sim.hc/1e-6
    f = Re/24
    f[Re>48] = 0.5
    f[Re<1] = 24
    return f

def bed_stress(sim):
    """
    """
    f = kirst_f(sim)
    tau_b = f/8*1000*sim.vc**2
    return tau_b


def stress_measures(summary):
    """
    tau_edge : vegetation patch edge
    tau_v : vegetation patch     
    tau_b : bare soil areas
 
    """
    for key in summary.index:
        sim = summary.loc[key]

        import cv2
        kernel = np.ones((3,3),np.uint8)

        dilation = cv2.dilate(sim.veg,kernel,iterations = 2)

        boundary = dilation-sim.veg
        boundary[sim.xc > sim.L/2] = 0
        x, y = np.where(boundary == 1)

        tau_bed = bed_stress(sim)
        tau_t = total_stress(sim)

        summary.at[key, 'tau_edge'] = tau_t[sim.i_tr, x, y].mean()            

        summary.at[key, 'tau_edge'] = tau_t[sim.i_tr, x, y].mean()    
        summary.at[key, 'tau_v'] = tau_t[sim.i_tr, sim.veg == 1].mean()
        veg = sim.veg.copy()

        summary.at[key, 'tau_b'] = tau_t[sim.i_tr, veg == 0].mean()    
        summary.at[key, 'tau_b_std'] = tau_t[sim.i_tr, veg == 0].std()
        summary.at[key, 'tau_bed_b'] = tau_bed[sim.i_tr, veg == 0].mean()        

        summary.at[key, 'V_edge'] = sim.uc[sim.i_tr, x, y].mean()
        summary.at[key, 'h_edge'] = sim.hc[sim.i_tr, x, y].mean()
        summary.at[key, 'U_edge'] = sim.vc[sim.i_tr, x, y].mean()

        if 'offset' in sim:

            boundary[sim.yc > sim.l/sim.dx/2 + sim.offset  ] = 0
            x_approach, y_approach = np.where(boundary == 1)    

            summary.at[key, 'tau_edge_a'] = tau_t[sim.i_tr, x_approach, y_approach].mean()    
    
    return summary


def compare_cases(summary, sim_list, scenario, use_tracers = 1):
    """
    Compares the 'both' scenario to 'hydraulics only' or 'hydrology only'
    """
    compare = pd.DataFrame(index = summary.query("scenario == 'both'").index)    

    for key in compare.index:
        
        sim = summary.loc[key]
      
        new_key = key.replace("scenario-both", 'scenario-' + scenario)
        try:
            sim_h = summary.loc[new_key]
        except:
            sim_h = summary.loc[new_key.replace("Ks_v-{0}".format(sim.Ks_v), "Ks_v-{0}".format(summary.query("scenario == 'both'").Ks_v.min()))]
        
        compare.at[key, 'pair_key'] = new_key

        if scenario == 'hydrology':
            assert  sim_h.Ks_v == sim.Ks_v    
            assert  sim_h.Ks.sum() == sim.Ks.sum()
        cols = ['tr', 'fV', 'p', 'phi_veg', 'alpha_v', 'alpha_b',
               'Ks_v', 'r', "C", "IF",
               'offset', 'I_v', 'I_v_tr', 'I_v_rec', 'IF_v', 
               # 'tau_b', 'V_edge', 'U_edge', 'tau_edge', 
               'flashy',  'excess', 'K_v'] 

        if type(sim_list) == list:
            sim_list = pd.DataFrame(sim_list)

        cols = cols + list(sim_list.columns)
        cols = list(set(cols))
        
        if use_tracers == 1:        
            cols = cols + ['divert', 'curve', 'vcount', 'diverted_mean', 'captured_mean']

        for col in cols:
            if col in sim:
                compare.at[key, col] = sim[col]

        compare.at[key, 'p-Ks_v'] = sim.p - sim.Ks_v      
        
        compare.at[key, 'Delta_C'] = (sim.C - sim_h.C)*100
        
        compare.at[key, 'Delta_U'] = (np.mean(sim.vc[sim.i_tr]) - 
                                      np.mean(sim_h.vc[sim_h.i_tr]))/np.mean(sim.vc[sim.i_tr])*100
        
        compare.at[key, 'Delta_h'] = (sim.h_max.mean() - sim_h.h_max.mean())/sim.h_max.mean()*100
        
        cols = [# 'tau_b', 'tau_v', 'tau_edge',
                 'I_v', 'I_v_tr', 'I_v_rec', 'IF_v', 'flashy', "C", "IF"]
        
        if use_tracers == 1:        
            cols = cols + ['divert', 'curve', 'vcount', 'diverted_mean', 'captured_mean']

        for col in cols:  
            
            compare.at[key, 'Delta_' + col ] = (sim[col] -  sim_h[col])

        # for col in ['V_edge', 'h_edge', 'U_edge', 'V_edge',  'V_std']:  
            
        #    compare.at[key, 'Delta_' + col ] = (sim[col] -  sim_h[col])
            
        veg = sim.veg.copy()
        veg[sim.xc > sim.L/2] = 0
        veg = veg == 1
        sim.vc[sim.vc > 1] = np.nan
    
        # compare.at[key, 'tau_n'] = sim.tau_t*(1/(1-sim.phi_veg)-1)**0.5
    
    return compare



### Unused snippets ####################################

def upper_approach(sim):

    upper = pd.DataFrame()

    edges = (np.diff(sim.veg, 1, axis = 0, prepend= 0) )
    midpoint = sim.l/sim.dx/2 + sim.offset 

    temp = np.where(edges == -1)
    temp = pd.DataFrame(temp, index = ['x', 'y']).T
    temp = temp.query("y <= {0}".format(midpoint))


    temp['offset'] = 0
    temp['V'] = sim.uc[sim.i_tr, temp['x'], temp['y']]*100
    temp['h'] = sim.hc[sim.i_tr, temp['x'], temp['y']]*100    
    temp['U'] = sim.vc[sim.i_tr, temp['x'], temp['y']]*100        
    upper = upper.append(temp)
    
    for offset in range(3):
        temp['x'] = temp['x']+1
        temp['offset'] += 1
        temp['V'] = sim.uc[sim.i_tr, temp['x'], temp['y']]*100
        temp['h'] = sim.hc[sim.i_tr, temp['x'], temp['y']]*100    
        temp['U'] = sim.vc[sim.i_tr, temp['x'], temp['y']]*100        
        
        upper = upper.append(temp)
    
    return upper

def lower_approach(sim):

    lower = pd.DataFrame()

    edges = (np.diff(sim.veg, 1, axis = 0, prepend= 0) )
    midpoint = sim.l/sim.dx/2 + sim.offset 

    temp = np.where(edges == 1)
    temp = pd.DataFrame(temp, index = ['x', 'y']).T
    temp = temp.query("y <= {0}".format(midpoint ))
    temp['x'] = temp['x']-3
    
    temp['offset'] = -2
    temp['V'] = sim.uc[sim.i_tr, temp['x'], temp['y']]*100
    temp['h'] = sim.hc[sim.i_tr, temp['x'], temp['y']]*100    
    temp['U'] = sim.vc[sim.i_tr, temp['x'], temp['y']]*100        
    lower = lower.append(temp)
    
    for offset in range(2):
        temp['x'] = temp['x']+1
        temp['offset'] += 1
        temp['V'] = sim.uc[sim.i_tr, temp['x'], temp['y']]*100
        temp['h'] = sim.hc[sim.i_tr, temp['x'], temp['y']]*100    
        temp['U'] = sim.vc[sim.i_tr, temp['x'], temp['y']]*100        
        
        lower = lower.append(temp)

    return lower

