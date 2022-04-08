#!/usr/bin/env python3

kb = 8.617e-5   # Boltzmann in eV/K

import argparse
import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.integrate import odeint
from labellines import labelLine, labelLines


# set LaTeX font
font = {'family' : 'Deja Vu Sans',
        'serif' : ['Computer Mordern Roman'],
        'size' : 16,}
plt.rc('font',**font)
plt.rc('text',usetex=True)

def TST(nu0,E,T):
    return nu0*np.exp(-E/kb/T)

def RE_inter(y, t, nui1,nu1i,nui2,nu2i,l1,l2,l_0,c_1):
    c_i, c_2 = y
    dydt = [+l1*c_1*(1-c_i)*nu1i/l_0
            +l2*c_2*(1-c_i)*nu2i/l_0
            -nui1*c_i
            -nui2*c_i,                   # dci/dt
            1.0/l2*(nui2*c_i*l_0
                -l2*c_2*(1-c_i)*nu2i)   # dc2/dt
            ]
    return dydt

def Jac_inter(y, t, nui1,nu1i,nui2,nu2i,l1,l2,l_0,c_1):
    c_i, c_2 = y
    jac  = [[-l1*c_1*nu1i/l_0-l2*c_2*nu2i/l_0-nui1-nui2,l2*(1-c_i)*nu2i/l_0],
            [1.0/l2*(nui2*l_0+l2*c_2*nu2i),1.0/l2*(-l2*(1-c_i)*nu2i)]]
    return jac


if __name__=="__main__":
    
    default = argparse.ArgumentDefaultsHelpFormatter
    parser  = argparse.ArgumentParser(formatter_class=default)
    parser.add_argument('-c1',metavar='c1',type=float,default=1e20,help='concentration in material 1 m-3')
    parser.add_argument('-E',metavar='E',type=str,default='0.3,0.8,1.0,0.4',help='the 4 energy barriers of the model as follow: E1->i,Ei->1,Ei->2,E2->i (eV)')
    parser.add_argument('-l',metavar='l',type=str,default='1e-10,1e-10',help='I-I in mat 1 and 2 as follow:l1,l2 (m)')
    parser.add_argument('-nu0',metavar='nu0',type=float,default=1e13,help='pre-exp factor of the different process (s-1)')
    parser.add_argument('-T',metavar='T',type=float,default=500,help='Temperature at which to solve the system (K)')
    parser.add_argument('-t',metavar='t',type=float,default=10,help='Total physical time of the simulation (s)')
    parser.add_argument('-dt',metavar='dt',type=float,default=0.001,help='time step for the resolution of the equation (s)')
    parser.add_argument('-n_i',metavar='n_i',type=float,default=1e20,help='Total concentration of interface site at the interface (m-2)')
    parser.add_argument('-atol',metavar='atol',type=float,default=1e-10,help='absolute tolerance for resolution')
    parser.add_argument('-rtol',metavar='rtol',type=float,default=1e-15,help='relative tolerance for resolution')
    parser.add_argument('--plot',action='store_true',default=False,help='plot output')

    c1 = 1e18       # cm in W in steady-state
    E_list = [[0.39,1.00,0.53,0.39],
              [0.39,1.05,0.58,0.39],
              [0.39,1.10,0.63,0.39],
              [0.39,1.15,0.68,0.39],
              [0.39,1.20,0.73,0.39],
              [0.39,1.25,0.78,0.39],
              [0.39,1.30,0.83,0.39],
              [0.39,1.35,0.88,0.39],
              [0.39,1.40,0.93,0.39],
              [0.39,1.45,0.98,0.39],
              [0.39,1.50,1.03,0.39],
              [0.39,1.55,1.08,0.39],
              [0.39,1.60,1.13,0.39]]
    l1 = 110e-12
    l2 = 65e-12
    n_i = 1e19      # interface site m-2
    nu0 = 1e13
    T   = 475
    t_fin = [1e-1,5e-1,1e0,5e0 ,1e1,5e1,1e2,5e2,1e3,5e3,1e4,5e4,1e5]
    
    dt = [1e-4,1e-4,1e-3,1e-3,1e-2,1e-2,1e-1,1e-1,1e0,1e0,1e1,1e1,1e1]
    
    a_tol = 1e-10
    r_tol = 1e-15



    c_interface  = []
    c_mat_2      = []
    c_int_eq     = []
    tau_interface = []
    tau_mat_2     = []
    t_simu        = []

    for r in range(len(E_list)):
        N_t =int( t_fin[r]/dt[r])
        t  = np.linspace(0,t_fin[r],num=N_t)
        t_simu.append(t)
        E_1i = E_list[r][0]
        E_i1 = E_list[r][1]
        E_i2 = E_list[r][2]
        E_2i = E_list[r][0]

        nu_1i = TST(nu0,E_1i,T)
        nu_i1 = TST(nu0,E_i1,T)
        nu_i2 = TST(nu0,E_i2,T)
        nu_2i = TST(nu0,E_2i,T)

        # Definition of the normalisation concentration (from the ceoncentration of interface site)
        l_0 = 1.0/np.sqrt(n_i)  # Distance in m
        n_0 = 1.0/l_0/l_0/l_0   # concentration in m-3

        # initialisation
        y0 = [0,0]   # initialisation y = [c1,ci,c2] normalized to n_0 or n_i
        sol = odeint(RE_inter, y0, t,args=(nu_i1,nu_1i,nu_i2,nu_2i,l1,l2,l_0,c1/n_0),Dfun=Jac_inter,mxstep=5000,rtol=r_tol,atol=a_tol)

        c_i = sol[:,0]
        c_2 = sol[:,1]
        c_interface.append(c_i)
        c_mat_2.append(c_2)
        #   STEADY-STATE Value
        c_i_eq = 1.0/(1.0+n_i*nu_i1/(l1*c1*nu_1i))
        print('At steady-state, ci_eq = ','%6.4e'%c_i_eq+' => ci = ','%3.2e'%float(c_i_eq*n_i)+' m-2')
        c_2_eq = l1/l2*nu_i2/nu_2i*nu_1i/nu_i1
        print('At steady-state, c2/c1 = ','%6.4e'%c_2_eq+' => c2 = ','%3.2e'%float(c_2_eq*c1)+' m-3')
        c_int_eq.append(c_i_eq)
        #   Time constant (to reach 95 % of the max)
        if max(c_i)>0.95*c_i_eq:
            t_i = t[abs(c_i-0.95*c_i_eq).tolist().index(min(abs(c_i-0.95*c_i_eq)))]
            print('t_95_i = '+'%3.2e'%t_i+' s')
        else:
            t_i = max(t)
            print('Simulation ime not enough to reach steady state i, increase -t')
        if max(c_2)>0.95*c_2_eq*c1/n_0:
            t_2 = t[abs(c_2*n_0/c1-0.95*c_2_eq).tolist().index(min(abs(c_2*n_0/c1-0.95*c_2_eq)))]
            print('t_95_2 = '+'%3.2e'%t_2+' s')
        else:
            t_2 = max(t)
            print('Simulation time not enough to reach steady state 2, increase -t')
        tau_interface.append(t_i)
        tau_mat_2.append(t_2)



    plt.figure(figsize=(6,15),edgecolor='None',facecolor='None')

    plt.subplot(311)
    E_short = [0.53,0.98,1.03,1.08,1.13]
    E_i2 = [E[2]/kb/T for E in E_list]
    color = "black"
    labels = [r'$< 21.5$', r'$23.9$', r'$25.2$', r'$26.4$', r'$27.6$']
    locs = np.flip(np.logspace(-4 + np.log10(8), 0 + np.log10(1), num=7))
    lines = []
    for i, E in enumerate(E_short[::-1]):
        r = E_i2.index(E/kb/T)
        label = labels[::-1][i]
        lines.append(plt.loglog(t_simu[r]/tau_interface[r],c_mat_2[r]*n_0/c1,linewidth=2,color=color, label=label, zorder=i+2))
        if E!=0.93/kb/T:
            labelLines(lines[-1], xvals=[locs[i], locs[i+1]], zorder=i+2, fontsize=9)
    
    plt.annotate('(a)',xy=(1.5e-3,0.8e5))

    plt.ylim(1e2,2e5)
    plt.xlim(1e-3,100)

    ax = plt.gca()
    ax.set_position([0.15,0.71,0.80,0.25])
    
    plt.xlabel(r'$t/\tau_\mathrm{i}$')
    plt.ylabel(r'$c_2/c_1$')


#    plt.xlim(1e-2,1e0)

    plt.subplot(312)
    for r in range(len(E_list)):
        plt.loglog(t_simu[r]/tau_interface[r],c_interface[r]/c_int_eq[r],linewidth=2,color='k')

    
    plt.ylabel(r'$\theta_\mathrm{i}/\theta_\mathrm{i}^\mathrm{eq}$')
    plt.xlabel(r'$t/\tau_\mathrm{i}$')

    plt.annotate('(b)',xy=(1.5e-3,6e-1))
    plt.xlim(1e-3,100)
    ax = plt.gca()
    ax.set_position([0.15,0.39,0.80,0.25])

#    plt.xlim(1e-2,1e0)

    plt.subplot(325)
    x = [E[2]/kb/T for E in E_list]
    plt.semilogy(x,tau_interface,'o',label=r'$\tau_\mathrm{i}$',markersize=6,markeredgecolor=(0,0,0.),markerfacecolor='None',markeredgewidth=2)
    plt.semilogy(x,tau_mat_2,'d',label=r'$\tau_2$',markersize=6,markeredgecolor="tab:orange",markerfacecolor='None',markeredgewidth=2)
    plt.legend(loc='lower right')
    plt.annotate('(c)',xy=(15,3e02))

    plt.xlabel(r'$E_\mathrm{i\rightarrow2}/k_\mathrm{B}T$')
    plt.ylabel(r'$\tau^{95\%}$ (s)')
    plt.yticks([1e-2,2e-2,3e-2,4e-2,6e-2,7e-2,8e-2,9e-2,
        1e-1,2e-1,1e-1,4e-1,5e-1,6e-1,7e-1,8e-1,9e-1,
        1e-0,2e-0,3e-0,4e-0,5e-0,6e-0,7e-0,8e-0,9e-0,
        1e1,2e1,3e1,4e1,5e1,6e1,7e1,8e1,9e1,
        1e2,2e2,3e2,4e2,5e2,6e2,7e2,8e2,9e2,
        1e3])
    plt.ylim(1e-2,1e3)
    ax = plt.gca()
    ax.set_position([0.15,0.06,0.32,0.25])
    # plt.xticks([0.5,0.7,0.9,1.1])

    plt.subplot(326)
    x = [E[2]/kb/T for E in E_list]
    plt.semilogy(x,c_int_eq,'s',markersize=6,markeredgecolor=(0,0,0),markerfacecolor='None',markeredgewidth=2)
    plt.ylim(1e-5,1e0)
    plt.annotate('(d)',xy=(15,3e-1))

    plt.ylabel(r'$\theta_\mathrm{i}^\mathrm{eq}$')
    plt.xlabel(r'$E_\mathrm{i\rightarrow 2}/k_\mathrm{B}T$')

    plt.yticks([1e-5,2e-5,3e-5,4e-5,6e-5,7e-5,8e-5,9e-5,
        1e-4,2e-4,3e-4,4e-4,5e-4,6e-4,7e-4,8e-4,9e-4,
        1e-3,2e-3,3e-3,4e-3,5e-3,6e-3,7e-3,8e-3,9e-3,
        1e-2,2e-2,3e-2,4e-2,5e-2,6e-2,7e-2,8e-2,9e-2,
        1e-1,2e-1,3e-1,4e-1,5e-1,6e-1,7e-1,8e-1,9e-1,
        1e-0])
    # plt.xticks([0.5,0.7,0.9,1.1])

    ax = plt.gca()
    ax.set_position([0.62,0.06,0.32,0.25])

#   plt.figure(figsize=(5,5))
#   plt.plot(t,c_1*n_0)
##    if plot:
#        plt.figure(figsize=(10,5))
#    #   Plot of C_2
#        plt.subplot(121)
#        plt.semilogy(t,c_2*n_0/c1,color=(0.8,0,0),linewidth=2)
#        plt.plot(t,np.ones(len(t))*c_2_eq,'--',color='k')
#        plt.annotate(r'$\frac{c_2^\mathrm{eq}}{c_1}$='+'%2.2e'%c_2_eq,xy=(max(t)*0.5,0.2*c_2_eq))
#        plt.annotate(r'$\tau_\mathrm{2}$='+'%2.2e'%t_2+' s',xy=(max(t)*0.5,0.1*c_2_eq))
#        plt.yticks([0,c_2_eq],['0',r'$\frac{c_2^\mathrm{eq}}{c_1}$'])
#        plt.xlim(0,max(t))
#        plt.ylim(0,1.1*c_2_eq)
#        plt.xlabel('time (s)')
#        plt.ylabel(r'$\frac{c_2}{c_1}$')
#        ax = plt.gca()
#        ax.set_position([0.08,0.15,0.40,0.8])
#    
#    #   Plot of C_interface
#        plt.subplot(122)
#        plt.plot(t,c_i,color=(0,0.5,0),linewidth=2)
#        plt.plot(t,np.ones(len(t))*c_i_eq,'--',color='k')
#        plt.annotate(r'$\theta_i^\mathrm{eq}$='+'%2.2e'%c_i_eq,xy=(max(t)*0.5,0.2*c_i_eq))
#        plt.annotate(r'$\tau_\mathrm{i}$='+'%2.2e'%t_i+' s',xy=(max(t)*0.5,0.1*c_i_eq))
#        plt.yticks([0,c_i_eq],['0',r'$\theta_i^\mathrm{eq}$'])
#        plt.xlim(0,max(t))
#        plt.ylim(0,1.1*c_i_eq)
#        plt.xlabel('time (s)')
#        plt.ylabel(r'$\theta_i$')
#        ax = plt.gca()
#        ax.set_position([0.57,0.15,0.40,0.8])
    plt.savefig("kinetic_interface.pdf")
    # plt.show()



