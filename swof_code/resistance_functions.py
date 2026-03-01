"""
Glossary:

alpha : coefficient in:

    U = h^m So^w alpha^-1
"""
import numpy as np
from scipy.optimize import fsolve


def get_exp(sim):
    if sim.scheme == "mann":
        a = 2/3.
    elif sim.scheme == "cyl":
        a = 0       
    elif sim.scheme == "cyl-dyn":
        a = 0               
    elif sim.scheme == "james":
        a = 0    
    elif sim.scheme == 'DW':
        a = 1/2.
        
    return a
    
################  Function from the kinematic resistance paper ###############

def get_alpha_2(So, qo,  s1, s2):
    """
    result from the roughness paper
    needs reformating to be more useful.
    """
    exponent = s2.eta*(s1.a+1)- s1.eta*(s2.a+1)
    
    return (s1.alpha**(s2.a + 1)*So**exponent*qo**(s2.a - s1.a))**(1/(s1.a+1))
   
def get_n_from_cyl(So, qo, alpha):
    """
    From the resistance schemes formulation
    ** Also need .... an upslope matching scheme
    A Giraldez and Woolhiser scheme
    **
    """
    return (alpha**(5/3)*So**(-1/3)*qo**(2/3))

def get_cyl_from_n(So, qo,  n):

    return (n**(3./5)*So**(1./5)*qo**(-2/5))
              

################  Cylinder functions ###############

    
def get_Rv(phi, D):
    """
    Vegetation hydraulic radius 
    
    Note that D and phi are not independent parameters -- only need to vary one. 
   
    phi: stem density
    D : stem diameter
    
    Increase with increasing step density, decreasing phi
    """    
    R_v = np.pi/4*(1-phi)/phi*D

    return R_v


def get_alpha_o(D, phi):
    """
    Non dynamic part of the cylinder array scheme, 
    i.e. alpha = alpha_o sqrt(Cd) 

    Delta_S : spacing between stems
    D : stem diameter
    phi : veg density
    """
    Delta_S2 = np.pi*D**2/(4*phi)

    return np.sqrt(1/Delta_S2)*np.sqrt(D/((1-phi)*2*9.8))


def get_RE_v(phi, D, U):
    """
    get vegetation reynolds number
    
    phi: stem density
    D : stem diameter
    """    
    RE_v = np.pi/4*(1-phi)/phi*D

    return U*RE_v/1e-6


def get_Cd_h(phi, D, q_o, h):    
    """
    Cd for cylinder array scheme (Wang 2018), 
    with q and h as input arguments
    """
    R_v = np.pi/4*(1-phi)*D/phi
    U = q_o/h
    Re_v = U*R_v/1e-6
    Cd = 50*(Re_v+ 1)**-0.43 + 0.7*(1-np.exp(-Re_v/15000))
    
    return Cd



def get_Cd(phi, D, U):
    """
    Cd for cylinder array scheme (Wang 2018)
    """
    RE_v = get_RE_v(phi, D, U)
    
    return 50*(RE_v+1)**-0.43 + 0.7*(1-np.exp(-RE_v/15000.))


def get_Cd_A(phi, D, U):
    """
    Cd for sheltering part of cylinder array scheme
    """
    re = get_RE_v(phi, D, U)
    
    return 50*re**-0.43 


def get_Cd_B(phi, D, U):
    """
    Cd for blocking part of cylinder array scheme
    """
    re = get_RE_v(phi, D, U)
    
    return 0.7*(1-np.exp(-re/15000))



def get_alpha_dynC(phi, D, U):
    """
    Estimate resistance given vegetation characteristics and flow velocity
    Wang 2018 with dynamic Cd
    """
    Cd =  get_Cd(phi, D, U)
    
    Delta_S2 = np.pi*D**2/4/phi

    alpha_o = np.sqrt((D)/(2*9.8*(1-phi)*Delta_S2)) 
    alpha = np.sqrt(Cd)*alpha_o
    
    return alpha

def balance_cyl_fixC(h, So, q, Cd, phi, D):
    """
    For the fied Cd cylinder array resistance formulation, 
    implicitly solve for the depth for which kinematic wave approximation holds

    """
    U = q/h
    
    alpha = get_alpha_fixC(Cd, phi, D)
    
    return U - So**0.5/alpha

def balance_cyl_dynC(h, So, q, phi, D):
    """
    For the dynamic cylinder array resistance formulation, 
    implicitly solve for the depth for which kinematic wave approximation holds
    """
    U = q/h
    
    alpha = get_alpha_dynC(phi, D, U)
    
    return U - So**0.5/alpha

def get_alpha_fixC(Cd, phi, D):
    """
    Estimate resistance given vegetation characteristics and drag coefficient
    for cylinder array scheme
    """
    Delta_S2 = np.pi*D**2/(4*phi)

    alpha = np.sqrt((Cd*D)/(2*9.8*(1-phi)*Delta_S2)) 
    
    return alpha


def equiv_Cd( U, So, phi, D):
    """
    For a given flow velocity and vegetation parameters, determinine an equivalent drag coefficient
    """
    alpha_o = get_alpha_o(D, phi)

    Cd = So/(alpha_o**2*U**2)

    return Cd


##### James now in the mix  ####

# def get_f(f, k, h, U):
#     """ 
#     Colebrook White
#     """
#     Re = h*U/1.14e-6
#     return 1/np.sqrt(f) + 2*np.log(k/(12*h) + 2.51/(Re*np.sqrt(f)))


# def fsolve_f(k, h, U):
#     """
#     implicitly solve for f
#     """
#     try:
#         f = fsolve(get_f, 0.001, args = (k, h, U))[0]
#     except:
#         f = [fsolve(get_f, 0.001, args = (k, h[i], U[i]))[0] for i in range(len(h))]
#         f = np.array(f) 
#     return f

# def get_alpha_James(h, q, k, phi, D):
#     """
#     Estimate alpha given roughness k and dynamic Cd, and known vegetation properties    
#     """
#     U = q/h

#     Cd = get_Cd(phi, D, U) 
#     denom = (1-phi)*9.8*h

#     Delta_S2 = np.pi*D**2/4/phi
#     try:
#         f = fsolve(get_f, 0.001, args = (k, h, U))[0]
#     except:
#         f = [fsolve(get_f, 0.001, args = (k, h[i], U[i]))[0] for i in range(len(h))]
#         f = np.array(f)
#     num = f/8 + 2*h*Cd*phi/(np.pi*D)

#     A = np.sqrt(num/denom)

#     return A


def get_alpha_James_fix(h, f, Cd,  phi, D):
    """
    Estimate alpha given fixed DW f and fixed Cd, and known vegetation properties
    """
    
    denom = (1-phi)*9.8*h
    
    num = f/8 + 2*h*Cd*phi/(np.pi*D)
    
    A = np.sqrt(num/denom)
    
    return A

def get_alpha_James_dynC(h, q, f, phi, D):
    """
    Estimate alpha given fixed DW f and dynamic Cd, and known vegetation properties    
    """
    U = q/h
    
    Cd = get_Cd(phi, D, U) 

    denom = (1-phi)*9.8*h

    Delta_S2 = np.pi*D**2/4/phi

    num = f/8 + 2*h*Cd*phi/(np.pi*D)
    
    A = np.sqrt(num/denom)
    
    return A


# def get_alpha_James_dynf(h, q, k, phi, D, Cd):
#     """
#     Estimate alpha given roughness k and dynamic Cd, and known vegetation properties    
#     """
#     U = q/h

#     # Cd = get_Cd(phi, D, U) 
#     denom = (1-phi)*9.8*h

#     Delta_S2 = np.pi*D**2/4/phi
#     try:
#         f = fsolve(get_f, 0.001, args = (k, h, U))[0]
#     except:
#         f = [fsolve(get_f, 0.001, args = (k, h[i], U[i]))[0] for i in range(len(h))]
#         f = np.array(f)
#     num = f/8 + 2*h*Cd*phi/(np.pi*D)

#     A = np.sqrt(num/denom)

#     return A

### Get kinematic depth for the James schemes, with and without dynamic resistance formulations

def balance_James_fix(h, q, f, Cd, phi, D, So): 
    """
    For fixed Cd and fixed f, 
    implicitly solve for the depth for which kinematic wave approximation holds

    """
    U = q/h
    
    alpha = get_alpha_James_fix(h, f, Cd, phi, D)
    
    return U - So**0.5/alpha

def balance_James_dynC(h, q, f, phi, D, So): 
    """
    
    For dynamic Cd and fixed f, 

    implicitly solve for the depth for which kinematic wave approximation holds

    """
    U = q/h

    alpha = get_alpha_James_dynC(h, q, f, phi, D)
  
    return U - So**0.5/alpha



# def balance_James_dynf(h, q, k, Cd, phi, D, So): 
#     """
    
#     For fixed Cd and dynamic f, 
#     implicitly solve for the depth for which kinematic wave approximation holds

#     """
#     U = q/h

#     # what f do we predict for this velocity
#     f = fsolve_f(k, h, U)

#     alpha = get_alpha_James_dynf(h, q, k, phi, D, Cd)
  
#     return U - So**0.5/alpha


def balance_James_dyn(h, q, phi, D, So): 
    """
    For the dynamic cylinder array resistance formulation, 
    implicitly solve for the depth for which kinematic wave approximation holds
    """
    U = q/h

    # could add laminar for lower Re    
    f = 0.5    
    alpha = get_alpha_James_dynC(h, q, f, phi, D)
  
    return U - So**0.5/alpha


def kirst(q):
    
    if type(q) in [float, np.float64, int]:
        q = np.array([q])
    
    Re = q/1e-6
    
    f = 24/Re
    f[Re>48] = 0.5
    f[Re<1] = 24
    
    return f


def balance_James_kirst(h, q,  phi, D, So): 
    """
    estimate f given k
    
    For the dynamic cylinder array resistance formulation, 
    implicitly solve for the depth for which kinematic wave approximation holds

    """
    U = q/h

    f =  kirst(q)
    alpha = get_alpha_James_dynC(h, q, f, phi, D)
  
    return U - So**0.5/alpha


######### Equivalent resistance ###################


def equiv_n(U, So, q):
    """
    Determine n assuming So=Sf,
    for a given velocity and depth

    U = So^1/2 h^2/3 /n
    q = So^1/2 h^5/3 / n
    """
    h = q/U
    
    n = h**(2/3.)*So**0.5/U
    
    return n 


def equiv_DW_ff(U, So, q):
    """
    Determine f assuming So=Sf,
    for a given velocity and depth

    f = 8gh So / U^2
    """
    h = q/U
    
    f = 8*g*h/U**2
    
    return f 


######### Compare hydrographs between schemes ###################

# def get_QNRMSE(sim, sim2):
#     n = min(len(sim.Vol_bound_tot), len(sim2.Vol_bound_tot))
#     Q1 = sim.Vol_bound_tot[:n]
#     Q2 = sim2.Vol_bound_tot[:n]

#     Q_nrmse = np.zeros_like(Q1).astype(float)    
#     Q_nrmse = (Q1-Q2)**2
#     Q_nrmse = np.sqrt(Q_nrmse.mean())/((sim.p-sim.Ks_v)*sim.L)*100   
#     return Q_nrmse

