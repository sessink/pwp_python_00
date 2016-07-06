"""
File to contain ice model functions to go with PWP.py

"""

import numpy as np
import seawater as sw
import matplotlib.pyplot as plt
from IPython.core.debugger import Tracer
import timeit
import os
import PWP_helper as phf

debug_here = Tracer()

#define constants
L_ice = 333e3 #Latent heat of fusion (J kg-1) - from Hyatt 2006
rho_ice = 920. #density of ice (kg/m3) - from Hyatt 2006
k_ice = 2. #thermal conductivity of sea ice (W m-1 K-1) - from Thorndike 1992
c_ice = 2e6/rho_ice #heat capacity of sea ice (J kg-1 K-1) - from Thorndike 1992/Hyatt 2006
c_sw = 4183. #heat capacity of seawater (J kg-1 K-1) - from Hyatt 2006
sal_ice = 5. #salinity of sea ice (PSU) - from Hyatt 2006

melt_lyrs = 5 #TODO: play with this. PWP has issues with super thin, freshwater lenses.
thin_ice = 0.005 #m    

# def melt_all_ice(h_ice, temp_ice_surf, sal_sw, temp_sw):
#
#     print "Ice has completely melted..."
#
#     #compute freshening due to melt
#     dsal_sw_melt = h_ice*(sal_sw[:melt_lyrs].mean()-sal_ice)/(dz*melt_lyrs)
#     sal_sw[:melt_lyrs] = sal_sw[:melt_lyrs]+dsal_sw_melt
#
#     #compute temp loss due to melt
#     F_i = h_ice*rho_ice*L_ice #energy used to melt ice
#     dT_surf = F_i/(dz*rho_sw0*c_sw)
#     temp_sw[0] = temp_sw[0]-dT_surf
#
#     #reset ice
#     h_ice = 0.
#     temp_ice_surf = np.nan
#
#     #debug_here()
#     return h_ice, temp_ice_surf, temp_sw, sal_sw

def iceGrowthModel_ode(t, y, F_atm, F_ocean, temp_swfz, k_ice):
    
    """   
    # dT/dt =  (F_atm + F_ocean)/(c*rho*h)   (1)
    # dh/dt = -(k*(T-T_fz)/h - F_ocean)/L    (2)
    # where T is surface temp
        
    # in the model, let y1 = T and y2 = h   
    """
    
    dy0 = 2*(F_atm+F_ocean)/(c_ice*y[1]*rho_ice)
    dy1 = (-k_ice*(y[0]-temp_swfz)/y[1] - F_ocean)/L_ice
    
    return [dy0, dy1]
    

def get_ocean_ice_heat_flux(temp_sw, sal_sw, rho_sw, params):    
    
    rho_sw0 = rho_sw[0] #density of seawater in model layer 1 (kg/m3)
    temp_swfz = sw.fp(sal_sw[0], p=1) #freezing point of seawater at given salinity
    
    dT_surf = temp_sw[0] - temp_swfz
    F_sw = dT_surf*params['dz']*rho_sw0*c_sw #J per m2 per time step
    F_sw_dt = F_sw/params['dt'] #in W/m2

    return F_sw_dt
    
def melt_ice(h_ice_i, temp_ice_surf, temp_sw, sal_sw, rho_sw, q_net, dz, dt, source):
    
    """
    Function to melt ice
    
    source: 'ocean' or 'atm'
    
    q_net - heat flux (W/m2)
    
    only way ice temp changes here is if ALL the ice melts. Otherwise, the ice will be in 
    dis-equilibrium until next time step.
    
    """
    
    temp_swfz = sw.fp(sal_sw[0], p=1) #freezing point of seawater at given salinity
    
    if source == 'atm':
        q_net = q_net*dt
    
    
    #print "Applying ocean feedback..."  
    ice_melt_potential = q_net/(rho_ice*L_ice) #amount of ice that can be melted with available heat flux

    #check if heat flux can completely melt current ice 
    if ice_melt_potential >= h_ice_i:
        print "%s heat flux has melted all ice..." %source
        #reset ice thickness and temperature
        temp_ice_surf = np.nan
        h_ice_f = 0.
        h_ice_melt = h_ice_i
        
    else:
        #ice survives only melt a portion of the ice with total ice melt potential
        h_ice_melt = ice_melt_potential
        h_ice_f = h_ice_i - h_ice_melt
    
    #compute freshening due to melt:
    dsal_sw_melt = h_ice_melt*(sal_sw[:melt_lyrs].mean()-sal_ice)/(dz*melt_lyrs)
    sal_sw[:melt_lyrs] = sal_sw[:melt_lyrs]+dsal_sw_melt
    
    #energy used to melt ice
    q_ice = h_ice_melt*rho_ice*L_ice 
    
    if source == 'ocean':
        #compute ocean surface temp loss due to melt:
        rho_sw0 = rho_sw[0] #density of seawater in model layer 1 (kg/m3)
        dT_surf = q_ice/(dz*rho_sw0*c_sw)
        temp_sw[0] = temp_sw[0]-dT_surf
    
        if ice_melt_potential < h_ice_i:
            assert np.abs(temp_sw[0]-temp_swfz)<0.001, "Something went wrong. Ocean surface not at freezing temp after incomplete ice melt."
            
        print "F_sw: %.2f, h_i: %.2f, dh: -%.4f" %(q_net/dt, h_ice_i, h_ice_melt)
        
    elif source == 'atm' and q_net>q_ice:
        
        #use excess atmospheric heat to warm ocean
        qnet_rem = q_net - q_ice
        dT_surf = qnet_rem/(dz*rho_sw0*c_sw)
        temp_sw[0] = temp_sw[0]-dT_surf #shouldn't this be positive temp_sw[0]+dT_surf ???
        
        print "F_atm: %.2f, h_i: %.2f, dh: -%.4f" %(q_net/dt, h_ice_i, h_ice_melt)
    
    
    return h_ice_f, temp_ice_surf, temp_sw, sal_sw
        
   
def create_initial_ice(h_ice_i, temp_ice_surf, temp_sw, sal_sw, rho_sw, dz):
    
    temp_swfz = sw.fp(sal_sw[0], p=1) #freezing point of seawater at given salinity
    rho_sw0 = rho_sw[0] #density of seawater in model layer 1 (kg/m3)
    
    print "initiating ice growth..."
    ##create ice according to Hyatt 2006
    #first, create a thin layer of sea ice (eqn. 5.11, Hyatt 2006)
    h_ice_f = rho_sw0*c_sw*dz*(temp_swfz-temp_sw[0])/(rho_ice*L_ice)
    
    #compute salinity change in top layer due to brine rejection (eqn. 5.12, Hyatt 2006)
    dsal_sw = h_ice_f*(sal_sw[0]-sal_ice)/dz
    sal_sw[0] = sal_sw[0]+dsal_sw
    
    #set ice to freezing temp of water
    temp_ice_surf_f = temp_swfz
    
    #set ocean surface to freezing temp
    temp_sw[0] = temp_swfz
    
    
    return h_ice_f, temp_ice_surf_f, temp_sw, sal_sw
    
    
def modify_thin_ice(h_ice_i, temp_ice_surf_i, temp_sw, sal_sw, rho_sw, F_atm, F_ocean, params):
    
    """
    Use this method to:
    
    1. modify thin ice up to some threshold. Let's say 1cm (0.01m).
    2. melt thick or thin ice if available heat flux is sufficient. ODE solver cannot handle this.
    
    """
    
    temp_swfz = sw.fp(sal_sw[0], p=1) #freezing point of seawater at given salinity
    rho_sw0 = rho_sw[0] #density of seawater in model layer 1 (kg/m3)
    
    #available heat: heat available to warm and/or melt ice:
    total_available_heat = (F_atm+F_ocean)*params['dt']  
    #warming_sink: heat gain required to bring ice to freezing point:
    warming_heat_sink = (temp_ice_surf_i-temp_swfz)*c_ice*h_ice_i 
    #melting_heat_sink: heat gain required to completely melt ice of thickness, after it has warmed to the freezing point:
    melting_heat_sink = h_ice_i*rho_ice*L_ice
    #total ice heat sink: warming_heat_sink + melting_heat_sink
    total_ice_heat_sink = warming_heat_sink + melting_heat_sink
    
    if total_available_heat>0:
        
        print "melting thin ice..."
        
        if total_available_heat<= warming_heat_sink:
            #use heat warm the ice
            temp_ice_surf_f = temp_ice_surf_i + 2*available_heat/(c_ice*h_ice_i*rho_ice)
            #no change in ice thickness, since ice was not warmed to freezing temp
            dh_ice = 0
            #reset ocean surface to freezing temp:
            temp_sw[0] = temp_swfz
            
        else:
            
            """
            In this scenario total_available_heat exceeds heat needs to warm ice to freezing temp 
            Hence, use heat left over from warming to melt the ice and possibly warm the ocean.
            """
            
            available_heat_for_melt = total_available_heat-warming_heat_sink
            
            if available_heat_for_melt <= melting_heat_sink:
                
                """
                In this scenario, there is enough heat to warm the ice to its freezing temp, 
                but not enough to completely melt it.
                """
                
                #compute change in ice thickness
                dh_ice = -available_heat_for_melt/(rho_ice*L_ice)
                #reset ocean surface temp:
                temp_sw[0] = temp_swfz
                #set ice surface temp to freezing
                temp_ice_surf_f = temp_swfz
                
            else:
                
                """
                In this scenario, there is enough heat to warm the ice to its freezing temp, 
                AND completely melt it. That is, total_available_heat > total_ice_heat_sink
                """
                
                #first find available ATMOSPHERIC heat.
                available_ATM_heat_for_ocean_warming = F_atm*params['dt'] - (total_ice_heat_sink-F_ocean*params['dt'])
                
                #completely melt ice and use leftover ATMOSPHERIC heat to warm the ocean
                dh_ice = -h_ice_i
                dT_surf = available_ATM_heat_for_ocean_warming/(params['dz']*rho_sw0*c_sw)
                temp_sw[0] = temp_sw[0]+dT_surf
                
                #set ice surface temp to freezing
                temp_ice_surf_f = np.nan #since there is no ice
                
                print "Ice has completely melted."
                
        
    else:
        
        """
        if available heat is negative, grow ice according to Hyatt 2006 (eqn. 5.11). 
        
        With this approach, we assume that ice is essentially transparent to radiative
        heat fluxes. The intent is to mitigate the rapid ice growth that occurs when 
        we apply large negative heat fluxes to very thin ice.
        
        With this model, ice continues to grow/accumulate until it exceeds the thin ice
        threshold. 
        
        """
        
        print "growing thin ice..."
    
        #compute ocean surface temp change due to heat loss
        dT_surf = total_available_heat/(params['dz']*rho_sw0*c_sw)
        temp_sw[0] = temp_sw[0]+dT_surf
        
        #keep ice at the freezing point (might be redundant)
        temp_ice_surf_f = temp_swfz 
        
        #if there is enough heat to cool the ocean surface past its freezing point, use the remaining heat flux to create ice
        if temp_sw[0]+dT_surf < temp_swfz:
            
            temp_sw[0] = temp_swfz #set ocean to freezing point.
            dh_ice = rho_sw0*c_sw*params['dz']*(temp_swfz-temp_sw[0])/(rho_ice*L_ice)
            
        else:
            #this is the case where the heat loss is not enough to cool the ocean past the freezing point.
            dh_ice = 0
            
    #get final ice thickness
    h_ice_f = h_ice_i + dh_ice
  
    #compute salinity change in top layer due to brine rejection or ice melt (eqn. 5.12, Hyatt 2006)
    dsal_sw = dh_ice*(sal_sw[0]-sal_ice)/params['dz']
    sal_sw[0] = sal_sw[0]+dsal_sw
    
    print "F_atm: %.2f, F_ocean: %.2f, h_i=%.4f, dh: %.4f" %(F_atm, F_ocean, h_ice_i, dh_ice)
    
    #debug_here()
    return h_ice_f, temp_ice_surf_f, temp_sw, sal_sw
    
    
            

def modify_existing_ice(temp_ice_surf_i, h_ice_i, temp_sw, sal_sw, rho_sw, F_atm, F_ocean, params):
    """
    input:
    
    temp_sw - ocean temp profile
    sal_sw - ocean salinity profile
    rho_sw - ocean density profile
    z - vertical coordinates
    dt - time step in seconds
    F_atm - net surface heating (positive into the ocean)
    temp_ice_surf_i - initial surface ice temperature
    h_ice_i - initial ice thickness
    
    """
    
    from scipy.integrate import odeint
    import scipy.integrate as spi
    
    #TODO: incorporate ice fraction
    
    switch_algorithm = False
    
    #derive/set additional parameters
    rho_sw0 = rho_sw[0] #density of seawater in model layer 1 (kg/m3)
    temp_swfz = sw.fp(sal_sw[0], p=1) #freezing point of seawater at given salinity
    temp_ice_base = temp_swfz
    #F_atm_rem = 0. #heat left over after melting ice
    # h_ice_melt = 0. #ice melt due to ocean feedback
    # F_sw = 0. #ocean heat flux
    dz = params['dz']
    dt = params['dt']
    
    assert h_ice_i >=0., "Error! negative ice thickness. Something went terribly wrong!"
    
    #available heat: heat required to warm and/or melt ice:
    available_heat = (F_atm+F_ocean)*params['dt']  
    #warming_sink: heat gain required to bring ice to freezing point:
    warming_heat_sink = (temp_ice_surf_i-temp_swfz)*c_ice*h_ice_i 
    #melting_heat_sink: heat gain required to completely* melt ice of thickness, after it has warmed to the freezing point:
    melting_heat_sink = (h_ice_i-thin_ice)*rho_ice*L_ice #*or melt ice to the point where it becomes thin
    #total ice heat sink: warming_heat_sink + melting_heat_sink
    total_ice_heat_sink = warming_heat_sink + melting_heat_sink


    if h_ice_i <= thin_ice:
        #modeling the conductive heat flux through thin ice is numerically unstable
        #use simplified approach that treats ice as transparent, except when positive heating is applied.
        h_ice_f, temp_ice_surf_f, temp_sw, sal_sw = modify_thin_ice(h_ice_i, temp_ice_surf_i, temp_sw, sal_sw, rho_sw, F_atm, F_ocean, params)
        
        return h_ice_f, temp_ice_surf_f, temp_sw, sal_sw
        
        
    else:
        
        if available_heat > total_ice_heat_sink:
        
            # if there is enough avaiable heat to create thin ice or completely melt the ice,
            # use thin ice algorithm. This is to circumvent the numerical issues that arise with very thin ice.
            # (this if-block is separated from the previous for the sake of clarity.)
        
            h_ice_f, temp_ice_surf_f, temp_sw, sal_sw = modify_thin_ice(h_ice_i, temp_ice_surf_i, temp_sw, sal_sw, rho_sw, F_atm, F_ocean, params)
        
            return h_ice_f, temp_ice_surf_f, temp_sw, sal_sw
        
        print "modifying thick ice..."
        
        #if the ice is "thick", we can treat is using the Thorndike algorithm

        #define initial conditions and time step for ice growtg model
        y_init = np.array([temp_ice_surf_i, h_ice_i])
        t_end = dt  
        t_ice_model = np.linspace(0, t_end, 5e4)
        dt_ice_model = t_ice_model[1]-t_ice_model[0]

        #load and intialize ice growth model
        ode =  spi.ode(iceGrowthModel_ode)
        #ode.set_integrator('vode', nsteps=500, method='adams') # BDF method suited to stiff systems of ODEs
        ode.set_integrator('lsoda')
        ode.set_initial_value(y_init, 0)
        ode.set_f_params(F_atm*dt_ice_model, F_ocean*dt_ice_model, temp_swfz, k_ice*dt_ice_model)
        ts = []
        ys = []
    
        #debug_here()

        #run ice growth model
        while ode.successful() and ode.t < t_end:
            ode.integrate(ode.t + dt_ice_model)
            ts.append(ode.t)
            ys.append(ode.y)
            temp_ice = ode.y[0] 
            h_ice = ode.y[1]
            
            if temp_ice >=temp_swfz or h_ice<=thin_ice:
                print "ice has become too thin or warm. May become numerically unstable. Switching to thin ice algorithm..."
                switch_algorithm=True

                #abort and switch to thin ice algorithm.??
                break


        

        if switch_algorithm:
            
            #use simplified approach that treats ice as transparent, except when positive heating is applied.
            h_ice_f, temp_ice_surf_f, temp_sw, sal_sw = modify_thin_ice(h_ice_i, temp_ice_surf_i, temp_sw, sal_sw, rho_sw, F_atm, F_ocean, params)
        
            return h_ice_f, temp_ice_surf_f, temp_sw, sal_sw
            
            
        else:
            t = np.vstack(ts)
            y = np.vstack(ys)
            temp_ice_surf_f = y[-1,0] 
            h_ice_f = y[-1,1]

        #compute ice temp and thickness change
        dh_ice = h_ice_f - h_ice_i
        dtemp_ice_surf = temp_ice_surf_f - temp_ice_surf_i

        print "F_atm: %.2f, F_ocean: %.2f, h_i=%.2f, dh: %.4f" %(F_atm, F_ocean, h_ice_i, dh_ice)

        #check if all ice has melted
        assert h_ice_f >= 0, "Something weird has happened. Negative ice thickness..."
        # if h_ice_f <= h_ice_i:
        #
        #     print "Warning: ice model produced negative ice thickness. Not good..."
        #     print "Reseting to zero ice..."
        #
        #     h_ice, temp_ice_surf, temp_sw, sal_sw = melt_all_ice(h_ice, temp_ice_surf, sal_sw, temp_sw)
        #
        #     return h_ice, temp_ice_surf, temp_sw, sal_sw, qnet_rem
    

        #compute salinity change in top layer due to brine rejection (or ice melt?) (eqn. 5.12, Hyatt 2006)
        num_lyrs = 5 #TODO: play with this. PWP has issues with super thin, freshwater lenses.
        dsal_sw = dh_ice*(sal_sw[:num_lyrs].mean()-sal_ice)/(dz*num_lyrs)
        sal_sw[:num_lyrs] = sal_sw[:num_lyrs]+dsal_sw

        
        return h_ice_f, temp_ice_surf_f, temp_sw, sal_sw
    
    
    # if h_ice_i == 0.:
    #     print "initiating ice growth..."
    #     ##create ice according to Hyatt 2006
    #     #first, create a thin layer of sea ice (eqn. 5.11, Hyatt 2006)
    #     h_ice_f = rho_sw0*c_sw*dz*(temp_swfz-temp_sw[0])/(rho_ice*L_ice)
    #
    #     #compute salinity change in top layer due to brine rejection (eqn. 5.12, Hyatt 2006)
    #     dsal_sw = h_ice_f*(sal_sw[0]-sal_ice)/dz
    #     sal_sw[0] = sal_sw[0]+dsal_sw
    #
    #     #set ice to freezing temp of water
    #     temp_ice_surf_f = temp_swfz
    #
    #     #set ocean surface to freezing temp
    #     temp_sw[0] = temp_swfz

    #elif h_ice_i > 0:
        
        
           
        # h_ice_i = h_ice #initial ice thickness
        # temp_ice_surf_i = temp_ice_surf #initial surface ice temperature
        
        
        # if ocean_fb_ON:
        #     # Find ocean sensible heat flux.
        #     # That is, the heat flux required to bring surface temp to freezing.
        #     dT_surf = temp_sw[0] - temp_swfz
        #     F_sw = dT_surf*dz*rho_sw0*c_sw
        #     F_sw_dt = F_sw/dt
        #

        #check if ocean heat heat flux can completely melt current ice
        # h_ice_melt = F_sw/(rho_ice*L_ice) #amount of ice that can be melted
        # if h_ice_melt >= h_ice:
        #
        #     h_ice, temp_ice_surf, temp_sw, sal_sw =  melt_all_ice(h_ice, temp_ice_surf, sal_sw, temp_sw)
        #     return h_ice, temp_ice_surf, temp_sw, sal_sw, qnet_rem
            
        # else:
        #
        #     temp_sw[0] = temp_swfz

                # print "Ocean mixing has melted ice..."
                #
                # #compute freshening due to melt
                # # dsal_sw_melt = h_ice*(sal_sw[0]-sal_ice)/dz
                # # sal_sw[0] = sal_sw[0]-dsal_sw_melt
                # #
                #
                # dsal_sw_melt = h_ice*(sal_sw[:melt_lyrs].mean()-sal_ice)/(dz*melt_lyrs)
                # sal_sw[:melt_lyrs] = sal_sw[:melt_lyrs]+dsal_sw_melt
                #
                # #compute temp loss due to melt
                # F_i = h_ice*rho_ice*L_ice #energy used to melt ice
                # dT_surf = F_i/(dz*rho_sw0*c_sw)
                # temp_sw[0] = temp_sw[0]-dT_surf
                #
                # #reset ice
                # h_ice = 0.
                # temp_ice_surf = np.nan

                #debug_here()
            #     return h_ice, temp_ice_surf, temp_sw, sal_sw, qnet_rem
            #
            # else:
            #
            #     temp_sw[0] = temp_swfz
                

        #debug_here()
        #h_ice_i = h_ice #initial ice thickness after ocean melt
        

        


    

        

    
        
        
        
        # plt.figure()
        # plt.subplot(211)
        # plt.plot(y[:,0]); plt.ylabel('Ice temperature (C)')
        # plt.subplot(212)
        # plt.plot(y[:,1]); plt.ylabel('Ice thickness (m)')
        
         
    
    

    