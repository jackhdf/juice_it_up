# Numerical method : Runge-Kutta Cash-Karp or Bulirsch-Stoer (both adaptive step size)
# Simulates the motion of particles under the influence of various forces
# Command line : python Adaptive_RK_CK_particle_moons_jupiter_traj.py model_parameters_jupiter text_files/Initial_Conditions_Callisto.txt 1000 0 0
# simulates 1000 particles starting from Callisto with job_ID = 0 and script_ID = 0


from math import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import interp1d
from time import *
import spiceypy as spice
import cProfile
import sys
import importlib



plt.close("all")

## returns position of a planet in a given reference frame (defined by a global variable) at ephemeris time t with respect to a certain origin
def spice_GetPosition(id,origin,et): # the last three arguments are optional
    state, ltime = spice.spkezr(id,et,frame,'NONE',origin) # Return the position of a target body relative to an observing body
    return state[0:3]

## returns velocity of a planet in a given reference frame (defined by a global variable) at ephemeris time t with respect to a certain origin
def spice_GetVelocity(id,origin,et): # the last three arguments are optional
    state, tlime = spice.spkezr(id,et,frame,'NONE',origin) # Return the state (position and velocity) of a target body relative to an observing body
    return state[3:]


## uses SPICE to transform position and velocity contained in y from the the Jupiter Inertial Frame (JIF) to the IAU_JUPITER body-fixed frame
def spice_rot2inert(y,t):
    frame_from = 'IAU_JUPITER'
    frame_to = 'JUICE_JUPITER_IF_J2000' # 'JUICE_JSM'
    et = t

    rot_mat = spice.pxform(frame_from, frame_to, et)  # rotation matrix from IAU_JUPITER frame to JIF
    w = np.array([0.,0.,w_jup])
    transfo_mat = spice.rav2xf(rot_mat,-w)
    yout = np.dot(transfo_mat,y)
    return yout


## uses SPICE to transform position and velocity contained in y from the IAU_JUPITER body-fixed frame to the Jupiter Inertial Frame (JIF) t
def spice_inert2rot(y,t):
    frame_to = 'IAU_JUPITER'
    frame_from = 'JUICE_JUPITER_IF_J2000' # 'JUICE_JSM'
    et = t

    rot_mat = spice.pxform(frame_from, frame_to, et)  # rotation matrix from IAU_JUPITER frame to JIF
    w = np.array([0.,0.,w_jup])
    transfo_mat = spice.rav2xf(rot_mat,w)
    yout = np.dot(transfo_mat,y)

    return yout


## function with returns the value of k.T (k: Boltzmann constant, T: temperature) which is a function of the distance of the particle (used in the model of plasma drag)
def fun_kT(dist):
    f_kT = interp1d(abs_logkT,np.exp(ord_logkT),kind='linear')
    return(f_kT(dist))
    ## plot values
    #xnew = np.linspace(3.8,170,500)
    #plt.plot(abs_logkT,np.exp(ord_logkT),'o',xnew,f_kT(xnew),'-')
    #plt.show()

## function with returns the value of N (number of ion particle / volume) which is a function of the distance of the particle (used in the model of plasma drag)
def fun_N(dist):
    f_N = interp1d(abs_logN0,np.exp(ord_logN0),kind='linear')
    return(f_N(dist))
    ## plot values
    #xnew = linspace(3.8,170,500)
    #plt.plot(abs_logN0,np.exp(ord_logN0),'o',xnew,f_N(xnew),'-')
    #plt.show()

## returns the reference distance ([Rj]) for the plasma drag as a function of distance
def fun_H(dist):
    if dist > 1.:
        if dist <= 3.8:
            H = 1
        elif dist <= 5.5:
            H = 0.2
        elif dist <= 7.9:
            H = 0.2
        elif dist <= 20:
            H = (1.82 - 0.041*(dist/Rj))
        elif dist <= 170:
            H = 1
    else:
        H = 0
    return H

## returns 1 if particle is in the sun, and 0 if it is in the shadow of Jupiter
# y is position and velocity in IAU_JUPITER frame
# t is ephemeris time [s]
def srp(y,t):
    origin = 'SUN'
    vect_sun_jupiter = spice_GetPosition('JUPITER', origin, t)

    # Given the position of a particle wrt to Jupiter in IAU_JUPITER frame, compute particle-sun vector and test whether we have an occultation by Jupiter
    vect_jupiter_dust = y[0:3]  # position of a particle wrt to Jupiter in IAU_JUPITER frame (body fixed frame)
    vect_sun_dust = vect_sun_jupiter + vect_jupiter_dust  # Compute particle position vector dust wrt to Sun

    ## Compute closest distance between vect_sun_dust line and Jupiter ellipsoid
    radii = spice.bodvrd('JUPITER', "RADII", 3)[1]  # Get ellipsoid radii in body fixed reference frame (IAU_JUPITER) : Fetch from the kernel pool the double precision values of an item associated with a body.
    direction = - vect_sun_dust / np.sqrt(vect_sun_dust[0] ** 2 + vect_sun_dust[1] ** 2 + vect_sun_dust[2] ** 2)  # Get nearest point on ellipsoid to the line defined by the particle-sun line

    # vect_jupiter_dust defines the point from which we look at the sun along the direction 'dust-sun'
    # then we check whether the dust-sun line intersects Jupiter = check whether 'Jupiter is between the dust grain and the sun'
    p_near, dist = spice.npedln(radii[0], radii[1], radii[2], vect_jupiter_dust, direction)  # Find nearest point on a triaxial ellipsoid to a specified line and the distance from the ellipsoid to the line.
    # there is an error in the above line when vec_jupiter_dust is nan. it is nan because y becomes a nan vector (the values of position and velocity are too big because particle escapes the system)

    if dist == 0: return 1 # in the sun
    else: return 0 # in the shadow of Jupiter

## associated Legendre polynomials (normalized) which are used for the VIPAL model of the magnetic field
# x = cos(theta) and n and m are order coefficients for the VIPAL model
def P(x,n,m):
    if n >= m:
        coeff = ( (2*factorial(n-m)) / factorial(n+m) )**(1/2) # normalization coefficient

        if n == 0:
            if m == 0:
                return coeff * 1

        elif n == 1:
            if m == 0:
                return coeff * x
            elif m == 1:
                return coeff * (-1) * (1 - x**2)**(1/2)

        elif n == 2:
            if m == 0:
                return coeff * (1/2) * (3*(x**2) - 1)
            elif m == 1:
                return coeff * (-3*x) * (1 - x**2)**(1/2)
            elif m == 2:
                return coeff * 3 * (1 - x ** 2)

        elif n == 3:
            if m == 0:
                return coeff * (1 / 2) * (5 * (x ** 3) - 3*x)
            elif m == 1:
                return coeff * (-3/2) * (5*(x ** 2) - 1) * (1-x**2)**(1/2)
            elif m == 2:
                return coeff * 15*x * (1 - x ** 2)
            elif m == 3:
                return coeff * (-15) * (1 - x ** 2) ** (3 / 2)

        elif n == 4:
            if m == 0:
                return coeff * (1/8) * (35*(x**4) - 30*(x**2) + 3)
            elif m == 1:
                return coeff * (-5/2) * (7*(x**3) - 3*x) * (1 - x**2)**(1/2)
            elif m == 2:
                return coeff * (15/2) * (7*(x**2) - 1) * (1 - x**2)
            elif m == 3:
                return coeff * (-105*x) * (1 - x ** 2) ** (3 / 2)
            elif m == 4:
                return coeff * 105 * (1 - x ** 2)**2

    elif n < m: return 0 # n < m (this happens when we call P(x,n-1,m) (see definition of P(x,n,m) on wikipedia to verify that it equals 0 when n<m)

## VIPAL model of the magnetic field around Jupiter, expressed in body-fixed coordinates (IAU_JUPITER)
# see article for complete description of VIPAL model : "dynamics and distribution of Jovian dust ejected from the Galilean satellites" by Xiaodong Liu
def VIPAL_B(r, th, phi):
    x = cos(th) # argument of Legendre function
    Bvec = np.array([0.,0.,0.]) # initialization

    for n in range(1, degree + 1):  # increment over n from 1 to 5 (innermost summation sign)
        for m in range(n+1): # increment over m from 0 to n (innermost summation sign)
            comp_r = (n+1) * P(x,n,m) * (g[m,n]*cos(m*phi) + h[m,n]*sin(m*phi)) # radial component of magnetic field
            comp_theta = (1/sin(th)) * (np.sqrt((n+m)*(n-m)) * P(x,(n-1),m) - n*cos(th)*P(x,n,m)) * (g[m,n]*cos(m*phi) + h[m,n]*sin(m*phi)) # colatitude component of magnetic field
            comp_phi = (m/sin(th)) * P(x,n,m) * (g[m,n]*sin(m*phi) - h[m,n]*cos(m*phi)) # longitude component of magnetic field
            Bvec = Bvec + ((Rj/r)**(n+2)) * np.array([comp_r, comp_theta, comp_phi])

    return Bvec

## returns the charge [Coulomb] of a particle (which depends on its size, mass and distance from Jupiter, NOT if it is in the sun or in the shadow of Jupiter)
def fun_Q(y):
    # data points from graph of potential as a function of distance from Jupiter
    absc = np.array([1,2,3,4,5,6,7,8,9,10])
    ord = np.array([1.1, 4.4, 6.5, 4.14, -8.77, 1.4, 1.4, 1.5, 1.7, 1.87])

    dist = sqrt(y[0]**2 + y[1]**2 + y[2]**2) / Rj # [Rj]
    index = np.abs(absc-dist).argmin() # index of the closest value in abs to dist
    potential = ord[index]

    num_electrons = 700*potential*(r_dust*1e6) # number of electrons
    charge = (1.60217662e-19)*num_electrons
    #chargebis = 4*pi*epsilon0*r_dust*potential # more precise
    return charge


## fictitious centrifugal force for the rotating frame
def acc_centrifugal(y):
    w = np.array([0., 0., w_jup])
    #y_inert = spice_rot2inert(y,t)
    temp = np.cross(w, y[0:3])
    acc = - np.cross(w, temp) # here the position and velocity are expressed in the rotating frame
    return acc


## fictitious Coriolis force for the rotating frame
def acc_coriolis(y):
    w = np.array([0., 0., w_jup]) # angular velocity vector
    coriolis = - 2 * np.cross(w, y[3:]) # here the position and velocity are expressed in the rotating frame
    return coriolis


## plasma drag model
# nh is a matrix containing the number density for different heavy ions
# each row (5) contains values for a given distance from jupiter: |1.0|3.8|5.5|7.9|20|170| [Rj]
# mh is a vector containing the mass for different heavy ions
# r_dust is in meters, density is in kg/m^3
# dist is in km and vel_rel in km/s
def acc_PDdrag(y):
    dist = np.sqrt(y[0]**2 + y[1]**2 + y[2]**2) / Rj # distance of particle from Jupiter [Rj]
    vel_rel = y[3:] #velocity of the particle relative to the plasma (we suppose that the plasma has zero velociy)
    norm_vel = np.sqrt(vel_rel[0]**2 + vel_rel[1]**2 + vel_rel[2]**2) # norm of the velocity of the particle relative to the plasma

    l = atan(y[1]/y[0]) # longitude of the particle (with respect to the x-axis of the body-fixed frame: arbitrary but no more precisions in the article)
    if abs((y[2] / (dist * Rj))) <= float(1): # asin is valid
        Lambda = asin(y[2]/(dist*Rj)) # latitude [radians]
    elif (y[2]/(dist*Rj)) > float(1): # to treat the cases when the z-component is so big that y[2] is bigger than dist*Rj
        Lambda = pi/2
    elif (y[2]/(dist*Rj)) < -float(1):
        Lambda = - pi / 2

    # calculates the number density of different species of heavy ions [O+ | O++ | S+ | S++ | S+++ | Na+]
    if dist > 1.:
        if dist <= 3.8:
            H0 = fun_H(dist) # scalar [Rj]
            lambda_c = tan(alpha[0]) * cos(l - l0)
            my_nh = gk[0,:] * (N0*exp(r0[0]/dist - ((dist/H0-1)*(Lambda-lambda_c))**2)) # vector giving the value of N_k for different types of ions

        elif dist <= 5.5:
            N = fun_N(dist) # scalar [Rj]
            H0 = fun_H(dist) # scalar [Rj]
            kT = fun_kT(dist) # scalar [Joules]
            H = H0*((kT/E0[1])**(1./2)) # scalar [Rj]
            z0 = dist*tan(alpha[1])*cos(l-l0) # [Rj]
            my_nh = gk[1,:] * N * exp(-((dist*Lambda - z0)/H)**2)

        elif dist <= 7.9:
            N = fun_N(dist) # scalar [Rj]
            H0 = fun_H(dist) # scalar [Rj]
            kT = fun_kT(dist) # scalar [Joules]
            H = H0 * ((kT / E0[2]) ** (1. / 2)) # scalar [Rj]
            z0 = dist*tan(alpha[0])*cos(l-l0) # [Rj]
            my_nh = gk[2, :] * N * exp(-((dist * Lambda - z0) / H) ** 2)

        elif dist <= 20:
            N = fun_N(dist)  # scalar [Rj]
            H0 = fun_H(dist)  # scalar [Rj]
            z0 = (((7*dist-16)/30)*Rj) * cos(l-l0) # [km]
            H = (1.82-0.041*dist)*Rj  # scalar [km]
            #kT = E0[3] - E1 * exp(-((dist*Rj * Lambda - z0) / H) ** 2)  # scalar [Joules]
            my_nh = gk[3, :] * N * exp(-((dist*Rj * Lambda - z0) / H) ** 2)

        elif dist <= 170:
            N = fun_N(dist)  # scalar [Rj]
            H = fun_H(dist)  # scalar [Rj]
            z0 = r0[1] * tan(alpha[2]) * cos(l - l0*ratio*(dist-r0[1]))
            my_nh = gk[4, :] * N * exp(-((dist * Lambda - z0) / H) ** 2) # scalar [Rj]

        elif dist > 170:
            my_nh = np.zeros((6)) # beyond 170*Rj, the plasma drag is zero

    else : my_nh = np.zeros((6)) # no plasma drag inside of Jupiter

    acc = (- ((3*norm_vel)/(4*r_dust*density*1e-6)) * np.dot(my_nh, mh)) * vel_rel # using speed in rotating frame (formula from Liu is in inertial frame, but doesn't change anything)

    return acc


## Poynting-Robertson drag (only acts when in the sun)
def acc_PRdrag(y,t,flag_srp):
    if flag_srp == 1:  # in the sun
        origin = 'SUN'
        vect_sun_jupiter = spice_GetPosition('JUPITER', origin, t)
        vect_jupiter_dust = y[0:3]  # position of a particle wrt to Jupiter in IAU_JUPITER frame (body fixed frame)
        vect_sun_dust = vect_sun_jupiter + vect_jupiter_dust  # Compute particle position vector dust wrt to Sun
        direction = - y[3:] / np.sqrt(y[3]**2 + y[4]**2 + y[5]**2) # direction of PR drag vector : it is pointing in the opposite direction of the velocity vector
        R = np.sqrt(vect_sun_dust[0]**2 + vect_sun_dust[1]**2 + vect_sun_dust[2]**2) # distance of particle from the sun

        acc = (1/m_dust) * (((r_dust**2)*(1e-6)*Ls)/(4*(c**2))) * np.sqrt(mu_sun/(R**5)) * direction # [km/s^2] (be careful that r_dust is initially expressed in meters) formula is from Wikipedia
        return acc

    elif flag_srp == 0:  # in the shadow of Jupiter
        return np.array([0.,0.,0.])


## Lorentz force
def acc_lorentz(y):
    ## Constant magnetic field pointing towards the positive z-axis of IAU_JUPITER frame
    #W = np.array([0.,0.,w_jup]) # angular velocity vector
    #B = np.array([0.,0.,B_magn])# suppose a constant magnetic field parallele to z-axis in body-fixed frame

    ## VIPAL model of magnetic field
    radius, colat, lon = spice.recsph(y[0:3])  # r is distance from origin, colat is angle from the positive z-axis, lon is longitude in radians (between -pi and pi) (ALL in rotating frame)
    B = (4.2e-4)*VIPAL_B(radius, colat, lon)  # this function returns the vector of the magnetic field at position y [Tesla]

    # calculation of Lorentz acceleration (equation in Liu has the vector r in inertial frame)
    Q = fun_Q(y)
    lo = (Q/m_dust) * (np.cross(y[3:], B))

    return lo


def acc_grav_sun(y,t,flag_srp):
    origin = 'SUN'
    vect_sun_jupiter = spice_GetPosition('JUPITER', origin, t)

    ## Given the position of a particle wrt to Jupiter in IAU_JUPITER frame, compute particle-sun vector and test whether we have an occultation by Jupiter
    vect_jupiter_dust = y[0:3]  # position of a particle wrt to Jupiter in IAU_JUPITER frame (body fixed frame)
    vect_sun_dust = vect_sun_jupiter + vect_jupiter_dust  # Compute particle position vector dust wrt to Sun
    norm_dist = np.sqrt(vect_sun_dust[0]** 2 + vect_sun_dust[1] ** 2 + vect_sun_dust[2] ** 2) # norm of the distance between particle and sun
    acc_grav_sun = (-mu_sun / norm_dist ** 3) * vect_sun_dust[0:3] # here the position y is expressed in the rotating frame.
    # in reality we should be using y expressed in the inertial frame for the "real" forces (like those due to gravity), but here we only need the distance from the origin which is the same in both frames

    if flag_srp == 1:  # in the sun
        return (1-beta)*acc_grav_sun # we incorporate solar radiation pressure here
    elif flag_srp == 0:  # in the shadow of Jupiter
        return acc_grav_sun


## acceleration of the particle due to Jupiter
def acc_grav_primary(y):
    norm_r2 = np.sqrt(y[0] ** 2 + y[1] ** 2 + y[2] ** 2)
    direction = - y[0:3] / norm_r2
    acc = (mu / (norm_r2**2)) * direction

    return acc

## acceleration of the particle due to Jupiter taking into account oblateness of the planet
def acc_grav_primary_oblate(y):
    r = np.sqrt(y[0] ** 2 + y[1] ** 2 + y[2] ** 2) # norm of the position
    acc = ((-mu / r ** 3) * y[0:3]) * (1 + ((3 * J2 * (a ** 2)) / (2 * (r ** 2))) * (1 - (5 * (y[2] ** 2)) / (r ** 2)))  # verifier l'indice r(2) du code de Nico
    acc[2] = acc[2] - (mu * 3 * J2 * (a ** 2) * y[2]) / (r ** 5)

    return acc


## acceleration generated by the four moons on the particle
def acc_moons(y):
    dist = np.sqrt((y[0,0]-y[1:,0]) ** 2 + (y[0,1]-y[1:,1]) ** 2 + (y[0,2]-y[1:,2]) ** 2) # distances between the particle and the four moons (4x1)
    mat1 = np.vstack((y[0,0:3], y[0,0:3], y[0,0:3], y[0,0:3])) # position of particle (4x3)
    mat2 = mat1 - y[1:,0:3] # position of particle - position of each of the four moons (4x3)
    temp = (-G*mass_moons / (dist ** 3)) # scalar value for each moon (4x1)

    accs = np.vstack((temp[0]*mat2[0,:], temp[1]*mat2[1,:], temp[2]*mat2[2,:], temp[3]*mat2[3,:])) # the four vectors of acceleration pulling on the particle (4x3)

    net_acc = np.array([sum(accs[:,0]), sum(accs[:,1]), sum(accs[:,2])]) # we sum the accelerations of the four moons along the three x,y and z
    return net_acc


def f(t, y, moon_flag):
    # initialization
    ag = ags = ace = aco = alo = apr = apd = amoons = np.zeros((3))

    ## calculation of different kinds of accelerations (y is a vector)
    ## we are calculating trajectory of moon, so we only consider the gravity of jupiter
    if moon_flag == 1 :
        f1 = y[3:]  # derivative of position = velocity for the particle
        ag = acc_grav_primary(y)
        ace = acc_centrifugal(y)
        aco = acc_coriolis(y)

        f2 = ag + aco + ace # never more forces acting on the moons

    ## we are calculating trajectory of dust particle : y is a matrix
    ## to use the exact positions of the moons at each call of the function f ( because if we only updated the positions of the moons in the driver function, then there would errors in the positions in the moons when calling f
    elif moon_flag == 0:
        f1 = y[0,3:]  # derivative of position = velocity for the particle
        flag_srp = srp(y[0, :], t)  # determines if particle is in occultation with Jupiter
        ag = acc_grav_primary(y[0,:]) #acc_grav_primary_oblate(y[0,:]) # acc_grav_primary(y[0,:])
        ace = acc_centrifugal(y[0,:])
        aco = acc_coriolis(y[0,:])
        ags = acc_grav_sun(y[0,:], t, flag_srp)
        alo = acc_lorentz(y[0,:])
        apr = acc_PRdrag(y[0,:],t,flag_srp)
        apd = acc_PDdrag(y[0,:])
        amoons = acc_moons(y) # acceleration exerted by the moons on the dust particle

        f2 = ag + ags + ace + aco + alo + apr + apd + amoons

    acc = np.concatenate((f1, f2), 0) # RHS of ODE system

    all_accs = np.vstack((ag, ags, ace, aco, alo, apr, apd, amoons)) # ag, ags, ace, aco, alo, apr, apd, amoons (in that order always)

    return acc, all_accs


##calculate position and velocity of four galilean moons in body-fixed Jupiter frame at time t
def moons(y, t, h, moon_flag):
    yout = np.zeros((4,6)) # position and velocity of four moons
    origin = 'JUPITER'
    yout[0,:] = np.array(np.append(spice_GetPosition('IO', origin, t), spice_GetVelocity('IO', origin, t))) # initial state of Io
    yout[1,:] = np.array(np.append(spice_GetPosition('EUROPA', origin, t), spice_GetVelocity('EUROPA', origin, t))) # initial state of Europa
    yout[2,:] = np.array(np.append(spice_GetPosition('GANYMEDE', origin, t), spice_GetVelocity('GANYMEDE', origin, t))) # initial state of Ganymede
    yout[3,:] = np.array(np.append(spice_GetPosition('CALLISTO', origin, t), spice_GetVelocity('CALLISTO', origin, t))) # initial state of Callisto

    ## using rk_ck which works great, ie same as SPICE (rk4 works too)
    # yout[0,:], yerr, acc = rk_ck(t, y[1,:], h, moon_flag)
    # yout[1,:], yerr, acc = rk_ck(t, y[2,:], h, moon_flag)
    # yout[2,:], yerr, acc = rk_ck(t, y[3,:], h, moon_flag)
    # yout[3,:], yerr, acc = rk_ck(t, y[4,:], h, moon_flag)

    return yout

## RK4
def rk4(t, y, h,moon_flag):
    k1 = h * f(t, y, moon_flag)[0]
    k2 = h * f(t+h/2, y+0.5*k1, moon_flag)[0]
    k3 = h * f(t+h/2, y+0.5*k2, moon_flag)[0]
    k4 = h * f(t+h, y+k3, moon_flag)[0]

    y = y + (k1 + 2*k2 + 2*k3 + k4) / 6
    return y
    #return y[0,:]


## for Bulirsch-Stoer : polynomial interpolation to evalute n functions at t = 0 by fitting a polynomial to a
# sequence of estimates with progressively smaller values t = t_est, and corresponding function
# vectors yest
# iest is the number in the sequence of calls
# tp is a vector that stores the time intervals that were accomplished
# yz is the extrapolated function values
# dy is their estimated error
# WARNING : if the numerical method is slow, try using rational function extrapolation (see page 731 of Numerical Recipes)
def interp(iest, t_est, yest):
    c = np.zeros((6))
    t = np.zeros((kmaxx))
    d = np.zeros((6, kmaxx))
    t[iest-1] = t_est # save current independent variable
    dy = yest
    yz = yest

    if iest == 1: # store first estimate in first column
        d[:,0] = yest # why change d if it is not being returned by interp, and d is not used in stepper, so what's the point of d ?
    else:
        c = yest
        #for k1 in range(0,iest): # we cannot use (iest-1) as an upper bound because then if iest = 2, the code will never go into the for loop
        for k1 in range(1, iest): # k1 = 1..iest-1
            delta = 1. / (t[iest-k1-1] - t_est)
            f1 = t_est * delta
            f2 = t[iest-k1-1] * delta
            for j in range(6): # propagate tableau one diagonal more
                q = d[j,k1-1]
                d[j,k1-1] = dy[j]
                delta = c[j] - q
                dy[j] = f1*delta
                c[j] = f2*delta
                yz[j] = yz[j] + dy[j]
        d[:,iest-1] = dy

    return yz, dy


# ## modified midpoint step : integrates a system of ODE from time t to t+htot using nstep substeps
# # y is a matrix here, with the first line containing the state of the particle
# # ts is the current time
# # htot is the total step to be made (H in Numerical Recipes)
# # nstep is the number of substeps to be used.
# def mmid_bis(y, ts, htot, nstep):
#     h = htot / nstep # stepsize this trip
#     ym = y # z0 here, but in general it is z_m-1 (matrix)
#     yn = y + h*np.vstack((f(ts, y, 0)[0], y[1:,:])) # first step z1, but in general it is z_m+1 (matrix)
#
#     t = ts + h
#     h2 = 2.*h
#     yout = np.vstack((f(t, yn, 0)[0], yn[1:,:])) # temporary storage of RHS of ODE system (matrix)
#     for n in range(1,nstep): # general step
#         swap = ym + h2*yout # matrix
#         ym = yn # matrix
#         yn = swap # matrix
#         t = t + h
#         yout = np.vstack((f(t, yn, 0)[0], yn[1:,:])) # matrix
#
#     yout = 0.5 * (ym[0,:] + yn[0,:] + h*yout[0,:]) # we take the first line of the matrix
#
#     return yout


## modified midpoint step : integrates a system of ODE from time t to t+htot using nstep substeps
# y is a matrix here, with the first line containing the state of the particle
# ts is the current time
# htot is the total step to be made (H in Numerical Recipes)
# nstep is the number of substeps to be used.
# NOTE: in the matrices z_prev, z_current and z_future, the only valide values are contained in the first line (values for the particle
def mmid(y, ts, htot, nstep):
    h = htot / nstep # stepsize this trip
    z_prev = y # z_0
    z_current = z_prev + h*np.vstack((f(ts, z_prev, 0)[0], y[1:,:])) # z_1
    z_future = y # initialization
    # z_prev = z_m-1 ; z_current = z_m ; z_future = z_m+1

    for m in range(1,nstep): # general step
        z_future_ = z_prev + 2*h*np.vstack((f(ts+m*h, z_current, 0)[0], y[1:,:]))

    yout = 0.5 * (z_future[0,:] + z_current[0,:] + h*(f(ts+htot, z_current, 0)[0])) # we take the first line of the matrix
    return yout


# ALGORITHM ROUTINE
# one iteration of Cash-Karp
# yerr is estimate of local truncation error
# yout is incremented variable y
def rk_ck(t, y, h, moon_flag):
    k1 = h * f(t, y, moon_flag)[0]
    k2 = h * f(t+h*B[1,0], y+B[1,1]*k1, moon_flag)[0]
    k3 = h * f(t+h*B[2,0], y+B[2,1]*k1+B[2,2]*k2, moon_flag)[0]
    k4 = h * f(t+h*B[3,0], y+B[3,1]*k1+B[3,2]*k2+B[3,3]*k3, moon_flag)[0]
    k5 = h * f(t+h*B[4,0], y+B[4,1]*k1+B[4,2]*k2+B[4,3]*k3+B[4,4]*k4, moon_flag)[0]
    k6 = h * f(t+h*B[5,0], y+B[5,1]*k1+B[5,2]*k2+B[5,3]*k3+B[5,4]*k4+B[5,5]*k5, moon_flag)[0]

    accelerations = np.zeros((8,3)) # used to store values of different types of acceleration

    if moon_flag == 1: # y is a vector
        yout = y + B[6,1]*k1 + B[6,2]*k2 + B[6,3]*k3 + B[6,4]*k4 + B[6,5]*k5 + B[6,6]*k6 # incremented variable yout, using fifth-order RK
        y4 = y + B[7,1]*k1 + B[7,2]*k2 + B[7,3]*k3 + B[7,4]*k4 + B[7,5]*k5 + B[7,6]*k6 # embedded fourth-order RK
    elif moon_flag == 0: # y is a matrix
        yout = y[0,:] + B[6, 1] * k1 + B[6, 2] * k2 + B[6, 3] * k3 + B[6, 4] * k4 + B[6, 5] * k5 + B[6, 6] * k6  # incremented variable yout, using fifth-order RK
        y4 = y[0,:] + B[7, 1] * k1 + B[7, 2] * k2 + B[7, 3] * k3 + B[7, 4] * k4 + B[7, 5] * k5 + B[7, 6] * k6  # embedded fourth-order RK
        accelerations = f(t + h, np.vstack((yout, y[1:])), 0)[1]

    yerr = yout - y4 # truncation error : difference between fourth and fifth order method

    return yout, yerr, accelerations # vectors



# # Bulirsch-Stoer step with monitoring of local truncation error to ensure accuracy and adjust stepsize
# # htry is the stepsize to be attempted
# # yscal is the vector against which the error is scaled
# # hnext is normally the same value as htry
# def stepper_bs(y, tt, htry, yscal, t1, moonflag, hnext):
#     first = 1
#     epsold = -1.
#     a = np.zeros((imaxx))
#     alf = np.zeros((kmaxx+1,kmaxx+1))
#     nseq = [0, 2, 4, 6, 8, 10, 12, 14, 16, 18]
#     exitflag = 0
#     #d = np.zeros((6,kmaxx))
#     err = np.zeros((kmaxx)) # err[k] is the estimate of a new stepsize H_k to obtain convergence in a fixed column k
#     #t = np.zeros((kmaxx))
#     ysav = yseq = yerr = np.zeros((6))
#     kmax_local = tnew = kopt = 0 # not sure that this has to be initialized
#     km = 2
#     # we must also initialize kopt and kmax_local and ysav
#
#     if eps != epsold: # a new tolerance, so reinitialize
#         hnext = tnew = -1e29 # "impossible" values
#         eps1 = S1*eps
#         a[0] = nseq[0] + 1 # check for indices (when C says a[1] then it probably means the first element)
#         for k in range(kmaxx):
#             a[k+1] = a[k] + nseq[k+1] # compute work coefficients A_k
#         for iq in range(1,kmaxx+1): # 1..kmaxx
#             for k in range(iq): # 0..(iq-1)
#                 alf[k,iq] = eps1 ** ((a[k] - a[iq]) / ((a[iq] - a[0] + 1.) * (2*k - 1))) # compute alpha(k,q)
#         epsold = eps
#         temp = 0
#         for kopt in range(1,kmaxx): # determine optimal row number for convergence
#             if a[kopt+1] > a[kopt] * alf[kopt-1,kopt]:
#                 temp = kopt
#                 break # gets out of this for loop
#         kmax_local = temp
#
#     h = htry
#     ysav = y # save the starting values
#     if tt != tnew or h != hnext: # a new stepsize or a new integration
#         first = 1 # restablish the order window
#         kopt = kmax_local
#     reduct = 0
#     bool = True # for the infinite loop below
#     errmax = 0. # to expand the scope of errmax
#
#     while bool: # infinite loop
#         for k in range(1,kmax_local+1): # evaluate the sequence of modified midpoint integrations
#             tnew = tt + h
#             if tnew == tt:
#                 print('Step size underflow in stepper routine')
#             yseq = mmid(ysav, tt, h, nseq[k])
#             test = np.sqrt(h/nseq[k]) # squared, since error series is even
#             y, yerr = interp(k, test, yseq) # perform extrapolation: we fit a polynomial to data points (test,yseq)
#             if k != 1: # compute normalized error estimate epsilon(k)
#                 errmax = tiny
#                 for i in range(6):
#                     errmax = max(errmax, abs(yerr[i]/yscal[i]))
#                 errmax = errmax / eps
#                 km = k - 1
#                 err[km] = (errmax/S1) ** (1./(2*km+1))
#             if k != 1 and (k >= kopt - 1 or first): # in order window
#                 if errmax < 1.: # converged
#                     exitflag = 1
#                     break # get out of the for loop
#                 if k == kmax_local or k == kopt + 1: # check for possible stepsize reduction
#                     red = S2 / err[km]
#                     break
#                 elif k == kopt and (alf[kopt-1,kopt] < err[km]):
#                     red = 1. / err[km]
#                     break
#                 elif kopt == kmax_local and alf[km, kmax_local-1] < err[km]:
#                     red = alf[km,kmax_local-1] * (S2/err[km])
#                     break
#                 elif alf[km,kopt] < err[km]:
#                     red = alf[km,kopt-1] / err[km]
#                     break
#         if exitflag:
#             break # get out of for loop # bool = False # to break out of the infinite loop
#         red = min(red, redmin) # reduce stepsize by at least REDMIN and at most redmax
#         red = max(red, redmax)
#         h = h * red
#         reduct = 1
#         # try again (and start again at the top of the infinite loop
#
#     tt = tnew # successful step taken
#     hdid = h
#     first = 0
#     wrkmin = 1e35
#
#     for kk in range(0,km): # compute optimal row for convergence and corresponding stepsize
#         fact = max(err[kk], scalmx)
#         work = fact*a[kk+1]
#         if work < wrkmin:
#             scale = fact
#             wrkmin = work
#             kopt = kk + 1
#
#     hnext = h / scale
#     if kopt >= k and kopt != kmax_local and not reduct: # check for possible order increase, but not if stepsize was just reduced
#         fact = max(scale / alf[kopt-1,kopt], scalmx)
#         if a[kopt+1]*fact <= wrkmin:
#             hnext = h / fact
#             kopt = kopt + 1
#
#     acc = f(tt, np.vstack((y,ysav[1:,:])), 0)[1] # to obtain accelerations
#     yout = y
#
#     return yout, tt, hnext, hdid, acc

## WORKING VERSION OF BULIRSCH STOER
def stepper_bs_bis(y, tt, htry, yscal, t1, moonflag, hnext):
    first = 1
    epsold = -1.
    a = np.zeros((imaxx)) # A_k
    alf = np.zeros((kmaxx,kmaxx)) # alpha(k,q) k = 1..8 ; q = 1..8 and k < q
    nseq = [0, 2, 4, 6, 8, 10, 12, 14, 16, 18]
    exitflag = 0
    #d = np.zeros((6,kmaxx))
    err = np.zeros((kmaxx)) # err[k] is the estimate of a new stepsize H_k to obtain convergence in a fixed column k
    #t = np.zeros((kmaxx))
    ysav = yseq = yerr = np.zeros((6))
    #kmax_local = tnew = kopt = 0 # only needs to be initialized if eps == epsold
    #km = 0 # in the case that k = 1:
    red = 0 # if kmaxx_local = 1, then the infinite loop is truly infinite

    if eps != epsold: # a new tolerance, so reinitialize
        hnext = tnew = -1e29 # "impossible" values
        eps1 = S1*eps
        a[0] = nseq[1] + 1
        for k in range(kmaxx):
            a[k+1] = a[k] + nseq[k+2] # compute work coefficients A_k : # equation (16.4.6) of Numerical Recipes
        for iq in range(1,kmaxx+1): # 1..kmaxx (because starting iq at zero does nothing interesting)
            for k in range(0,iq): # 0..(iq-1)
                alf[k,iq-1] = eps1 ** ((a[k+1] - a[iq]) / ((a[iq] - a[0] + 1.) * (2*k + 3))) # compute alpha(k,q) : # equation (16.4.10) of Numerical Recipes
        epsold = eps
        temp = 2 # temporary variable which represents the largest allowed column, kmax_local (we set temp to 1 so that we can at least use the first column of the tableau
        for kopt in range(1,kmaxx):  # determine optimal row number for convergence (kopt is q in Numerical Recipes)
            #print a[kopt + 1], (a[kopt] * alf[kopt - 1, kopt])
            if a[kopt+1] > (a[kopt] * alf[kopt-1, kopt]):  # condition (16.4.13) of Numerical Recipes
                temp = kopt
                break  # gets out of this for loop
        kmax_local = temp # note : indexing of

    h = htry
    ysav = y # save the starting values
    if tt != tnew or h != hnext: # a new stepsize or a new integration
        first = 1 # restablish the order window
        kopt = kmax_local
    reduct = 0
    bool = True # for the infinite loop below
    errmax = 0. # to expand the scope of errmax

    while bool: # infinite loop
        for k in range(1,kmax_local+1): # evaluate the sequence of modified midpoint integrations
            tnew = tt + h
            if tnew == tt:
                print('Step size underflow in stepper routine')
            yseq = mmid(ysav, tt, h, nseq[k])
            test = np.sqrt(h/nseq[k]) # squared, since error series is even
            y, yerr = interp(k, test, yseq) # perform extrapolation: we fit a polynomial to data points (test,yseq)
            if k != 1: # compute normalized error estimate epsilon(k)
                errmax = tiny
                for i in range(6):
                    errmax = max(errmax, abs(yerr[i]/yscal[i]))
                errmax = errmax / eps
                km = k - 1
                err[km-1] = (errmax/S1) ** (1./(2*km+1))
            if k != 1 and (k >= (kopt - 1) or first): # in order window
                if errmax < 1.: # converged
                    exitflag = 1
                    break # get out of the for loop
                if k == kmax_local or k == kopt + 1: # check for possible stepsize reduction
                    red = S2 / err[km-1]
                    break
                elif k == kopt and (alf[kopt-2,kopt-1] < err[km-1]):
                    red = 1. / err[km-1]
                    break
                elif kopt == kmax_local and alf[km-1, kmax_local-2] < err[km-1]:
                    red = alf[km-1,kmax_local-2] * (S2/err[km-1])
                    break
                elif alf[km-1,kopt-1] < err[km-1]:
                    red = alf[km-1,kopt-2] / err[km-1]
                    break
        if exitflag:
            break # get out of for loop # bool = False # to break out of the infinite loop
        red = min(red, redmin) # reduce stepsize by at least REDMIN and at most redmax
        red = max(red, redmax)
        h = h * red
        reduct = 1
        # try again (and start again at the top of the infinite loop

    tt = tnew # successful step taken
    hdid = h
    first = 0
    wrkmin = 1e35

    for kk in range(0,km): # compute optimal row for convergence and corresponding stepsize
        fact = max(err[kk], scalmx)
        work = fact*a[kk+1]
        if work < wrkmin:
            scale = fact
            wrkmin = work
            kopt = kk + 1

    hnext = h / scale
    if kopt >= k and kopt != kmax_local and not reduct: # check for possible order increase, but not if stepsize was just reduced
        fact = max(scale / alf[kopt-1,kopt], scalmx)
        if a[kopt+1]*fact <= wrkmin:
            hnext = h / fact
            kopt = kopt + 1

    acc = f(tt, np.vstack((y,ysav[1:,:])), 0)[1] # to obtain accelerations
    yout = y

    return yout, tt, hnext, hdid, acc




# STEPPER ROUTINE : Fifth-order Runge-Kutta step with monitoring of local truncation error to ensure accuracy and adjust stepsize
# CAREFUL : the exponents exp_shrink and exp_grow are negative so the fraction in the book is flipped
# we do not control the accumulation of local errors at each time step
# yscal is a vector against which the error is scaled
# htry is the stepsize to be attempted
# eps is required accuracy
# hdid is the stepsize that is actually accomplished
# hnext is the estimated next stepsize
# Delta_0 is the desired accuracy (Delta_0 is the error produced after taking a step h0)
# Delta_1 is the error produced after taking a step h1
# if Delta_1 > Delta_0, we need to shrink the time step
# if Delta_1 < Delta_0, we need to grow the time step

def stepper_rkck(y, t, htry, yscal, t1, moon_flag):
    yerr = np.zeros((6))
    h = htry # set stepsize to the initial trial value
    errmax = 0. # maximum error for all components of the ODE system (errmax = Delta_1 / yscal)

    while True : # infinite loop
        yout, yerr, acc = rk_ck(t, y, h,0) # take one step (yerr = Delta_1)

        errmax = 0. # evaluate accuracy
        for i in range(6):
            errmax = max(errmax, abs(yerr[i] / yscal[i])) # find the biggest error for each component of the ODE system
        errmax = errmax / eps # scale relative to required tolerance
        if errmax <= 1. : # step succeeded. compute next iteration by going back in the driver routine. if errmax <= 1, then Delta_0 >= Delta_1. so we must grow the size of the step, which will be done outside of the infinite loop
            break
        htemp = S*h*(errmax**exp_shrink) # if errmax > 1, then Delta_0 < Delta_1. the current stepsize is too big, so reduce stepsize. htemp is the size of the next step to be attempted.
        if h >= 0.: # truncation error too large, reduce stepsize by no more than a factor of ten (while keeping the stepsize above htemp)
            h = max(htemp,0.1*h)
        else: # we have a negative time step, thus we are doing backwards integration. Taking the min, means taking the biggest step in absolute value
            h = min(htemp, 0.1*h)
        tnew = t + h
        if tnew == t: # means that the stepsize is zero
            print('stepsize underflow in stepper routine')

    if errmax > errcon:
        hnext = S*h*(errmax**exp_grow) # we grow the next stepsize since we now have that  Delta_0 >= Delta_1
    else:
        hnext = 5. * h # no more than a factor of 5 increase

    hdid = h
    t = t + h

    return yout, t, hnext, hdid, acc



## a particle has hit a moon if the distance between the position vector of the particle and the position vector of the moon is smaller than the radius of the moon
def hitMoons(y):
    for i in range(np.size(mass_moons)):
        distance = np.sqrt((y[0,0] - y[i+1,0])**2 + (y[0,1] - y[i+1,1])**2 + (y[0,2] - y[i+1,2])**2) # distance between particle and center of the moon
        if distance <= radius_moons[i]:
            return True # particle hit the moon

    return False

## returns True if the particle has a trajectory that traverses a moon's surface
## we draw a straight line between the position y of the particle at time t  and the position of the particle for the previous iteration
## we discretize this line in smaller segments. For each endpoint of the segments, we call hitMoons to check if the particle has traversed the moons
## the discretization of the line is done such that each segment is no longer than the diameter of the smallest moon
## all calculations are done in the rotating frame
## we suppose that the velocity of the particle is constant between yprev and ycurrent
## t_prev is the time at yprev
def traverseMoons(yprev, ycurrent, t_prev, t_current):
    # only do the traverse check if we are reasonably close to any moon at the previous and current iteration
    distance_from_any_moon = 1000 # distance under which we check if the particle might traverse the moon
    distances_prev = distances_current = np.zeros((np.size(mass_moons)))
    pos_prev = yprev[0, 0:3]  # position of particle at the previous iteration
    pos_current = ycurrent[0, 0:3]  # position of particle at the current iteration
    for i in range(np.size(mass_moons)):
        distances_prev[i] = np.sqrt((yprev[0, 0] - yprev[i + 1, 0]) ** 2 + (yprev[0, 1] - yprev[i + 1, 1]) ** 2 + (yprev[0, 2] - yprev[i + 1, 2]) ** 2)  # distance between particle and center of the moon
        distances_current[i] = np.sqrt((ycurrent[0, 0] - ycurrent[i + 1, 0]) ** 2 + (ycurrent[0, 1] - ycurrent[i + 1, 1]) ** 2 + (ycurrent[0, 2] - ycurrent[i + 1, 2]) ** 2)  # distance between particle and center of the moon

    if (distances_prev < distance_from_any_moon).any() or (distances_current < distance_from_any_moon).any() : # the particle is too close to at least one moon at either the previous iteration or the current iteration
        print('')
        print('Approaching a moon...')
        delta_x = pos_current[0] - pos_prev[0] # positive if particle is advancing towards positive x ; negative otherwise
        delta_y = pos_current[1] - pos_prev[1]
        delta_z = pos_current[2] - pos_prev[2]
        min_radius = min(radius_moons)
        seg_num_x = abs(ceil(delta_x / min_radius)) # number of segments in the x-dimension
        seg_num_y = abs(ceil(delta_y / min_radius)) # number of segments in the x-dimension
        seg_num_z = abs(ceil(delta_z / min_radius)) # number of segments in the x-dimension
        max_seg_num = max(seg_num_x, seg_num_y, seg_num_z) # find the maximum number of segments of the 3 dimensions. this number is then the number of segments in all three dimensions
        seg_length_x = delta_x / max_seg_num # x-length of the segment (can be negative)
        seg_length_y = delta_y / max_seg_num # y-length of the segment
        seg_length_z = delta_z / max_seg_num # z-length of the segment
        points_x = pos_prev[0] + seg_length_x * np.linspace(0, max_seg_num, max_seg_num+1) # x-coordinates of the 3D points that define the segments
        points_y = pos_prev[1] + seg_length_y * np.linspace(0, max_seg_num, max_seg_num+1) # y-coordinates of the 3D points that define the segments
        points_z = pos_prev[2] + seg_length_z * np.linspace(0, max_seg_num, max_seg_num+1) # z-coordinates of the 3D points that define the segments
        all_points = np.vstack((points_x, points_y, points_z)) # position of the particle at different times between pos_prev and pos_current (dimension : 3 x (max_seg_num + 1))
        t = t_prev + (abs(t_current-t_prev) / max_seg_num) * np.linspace(0, max_seg_num, max_seg_num+1) # vector containing times at each endpoint of segments (assuming that the velocity of the particle is constant)

        for i in range(int(max_seg_num)+1):
            ytemp = np.zeros((5,6)) # stores the position of particle at i iterations after pos_prev, as well as the positions of the moons at i iterations after pos_prev
            ytemp[0,:] = np.append(all_points[:,i], np.zeros((1,3))) # we do not need velocity of particle here
            ytemp[1,:] = np.array(np.append(spice_GetPosition('IO', origin, t[i]), spice_GetVelocity('IO', origin, t[i])))
            ytemp[2,:] = np.array(np.append(spice_GetPosition('EUROPA', origin, t[i]), spice_GetVelocity('EUROPA', origin, t[i])))
            ytemp[3,:] = np.array(np.append(spice_GetPosition('GANYMEDE', origin, t[i]), spice_GetVelocity('GANYMEDE', origin, t[i])))
            ytemp[4,:] = np.array(np.append(spice_GetPosition('CALLISTO', origin, t[i]), spice_GetVelocity('CALLISTO', origin, t[i])))

            bool = hitMoons(ytemp)
            if bool == True:
                return True # one of the endpoints of the segments is inside the radius of one of the moons
        return False

    else :
        return False # the particle is too far away from any moon


# DRIVER ROUTINE : Runge-Kutta driver with adaptive stepsize control. Integrate starting values ystart[1..nvar] from x1 to x2 with accuracy eps
# ystart is the initial value of y
# eps is the desired accuracy
# n is number of variables/equations
# t1 and t2 are starting and end times of integration
# h1 is to be set as a guessed first stepsize
# hmin is the minimum allowed stepsize (can be zero)
# nok and nbad are respectively the number of good and bad (but retried and fixed) steps taken
# kount is the current number of time steps
# kmax is the maximum number of steps that can be stored
# tp and yp stores results (t,y) at approximate time intervals dxsav
# results are only stored at intervals greater than dxsav
# t is not incremented in for loop because its value is determined at each output of the stepper routine

def driver(ystart, t1, t2):
    t = t1
    tp = np.zeros((1,kmax))[0]
    yp = np.zeros((5,6,kmax))
    accp = np.zeros((8,3,kmax))
    orbitp = np.zeros((3,kmax))
    yscal = np.zeros((5,6))[0]
    hdid = hnext = 0 # initialized here for scope

    if (t2-t1) >= 0.:
        h = abs(h1) # forward integration
    else:
        h = -abs(h1) # backwards integration

    nok = nbad = kount = 0
    y = ystart
    acc = np.zeros((8,3))

    if kmax > 0: # if the maximum number of steps that can be stored is valid
        tsav = t - dtsav*2. # assures storage of first step

    for nstep in range(itermax): # take at most itermax steps
        yscal = abs(y[0,:]) + abs(h*f(t,y,0)[0]) + tiny # a matrix:  cfr p. 719 of Numerical Recipes : scaling used to monitor accuracy. This general-purpose choice can be modified if need be
        if kmax > 0 and kount < kmax-1 and abs(t-tsav) > abs(dtsav):
            tp[kount] = t # save first and intermediate results
            yp[:,:,kount] = y
            accp[:,:,kount] = acc

            state = spice_rot2inert(y[0,:],t) # transforms to JIF
            elts = spice.oscelt(state,t,mu) # orbital elements
            apoapsis = elts[0] / (1 - elts[1])  # apoapsis distance = perifocal distance / (1-eccentricity)
            semimaj = (elts[0] + apoapsis) / 2
            # v = sqrt(y[0, 3] ** 2 + y[0, 4] ** 2 + y[0, 5] ** 2)  # instantaneous speed
            # dist = sqrt(y[0, 0] ** 2 + y[0, 1] ** 2 + y[0, 2] ** 2)  # dist from origin
            # E = (v**2)/2 - mu/dist # specific orbital energy
            # semimaj = - mu / (2*E) # semi-major axis
            orbitp[:, kount] = np.array([semimaj, elts[2], elts[0]])   # stores semi-major axis, inclination and perifocal distance

            tsav = t
            kount = kount + 1
        if (t+h-t2)*(t+h-t1) > 0.: # term 1 and term 2 are of the same sign if LHS is positive
            h = t2-t # if stepsize can overshoot, decrease stepsize


        if integrator == 'rkck': # RK Cash-Karp
            y[0,:], t, hnext, hdid, acc = stepper_rkck(y, t, h, yscal, t1,0) # we tried with step h, and then actually accomplished step hdid
            y[1:, :] = moons(y, t, hdid,1)  # update position and speed of moons at time t (after executing this line, the matrix y is entirely updated
        elif integrator == 'bs': # Burlisch-Stoer
            y[0, :], t, hnext, hdid, acc = stepper_bs_bis(y, t, h, yscal, t1,0, hnext)  # we tried with step h, and then actually accomplished step hdid
            y[1:, :] = moons(y, t, hdid,1)  # update position and speed of moons at time t (after executing this line, the matri


        if hdid == h:
            nok = nok + 1 # good step was taken
        else:
            nbad = nbad + 1 # bad step was taken

        hitJupiter = sqrt(y[0,0]**2+y[0,1]**2+y[0,2]**2) <= Rj # True if the particle has hit Jupiter
        hitMoon = hitMoons(y)
        traverseMoon = traverseMoons(yp[:,:,kount-1], y, tp[kount-1], t)
        outsideHill = sqrt(y[0,0]**2+y[0,1]**2+y[0,2]**2) > hill_radius

        if (t-t2)*(t2-t1) >= 0. or hitJupiter or hitMoon or traverseMoon or outsideHill: # are we done? ie LHS is positive if t > t2 , ie we are done (second condition is to see if we hit Jupiter's surface
            #ystart = y # no real point to this line
            if kmax: # enters the loop if kmax != 0
                tp[kount] = t # save final step
                yp[:,:,kount] = y
                accp[:, :, kount] = acc

                state = spice_rot2inert(y[0, :], t)  # transforms to JIF
                elts = spice.oscelt(state, t, mu)  # orbital elements
                apoapsis = elts[0] / (1-elts[1]) # apoapsis distance = perifocal distance / (1-eccentricity)
                semimaj = (elts[0] + apoapsis)/2
                #v = sqrt(y[0, 3] ** 2 + y[0, 4] ** 2 + y[0, 5] ** 2)  # instantaneous speed
                #dist = sqrt(y[0, 0] ** 2 + y[0, 1] ** 2 + y[0, 2] ** 2)  # dist from origin
                #E = (v**2)/2 - mu/dist # specific orbital energy
                #semimaj = - mu / (2*E) # semi-major axis
                orbitp[:, kount] = np.array([semimaj, elts[2], elts[0]]) #np.array([elts[0], elts[2]])  # stores perifocal distance and inclination

                kount = kount + 1
            if hitJupiter:
                print('Simulation stopped because the particle hit Jupiter')
                print('')
            elif hitMoon or traverseMoon:
                print('Simulation stopped because the particle hit a moon')
                print('')
            elif outsideHill:
                print('Simulation stopped because the particle escaped the Hill sphere')
                print('')
            return tp, yp, accp, orbitp, kount, nok, nbad # normal exit
        #if abs(hnext) <= hmin :
        #    print('Step size too small in driver routine')
        h = hnext

    print('Maximum number of iterations reached in driver routine : returning results')
    return tp, yp, accp, orbitp, kount, nok, nbad  # Exit after maximum number of iterations reached


def load_initial_conditions(ic_file):

    ## read file of initial conditions for particles
    file = open(ic_file,'r')  # file contains initial positions and velocities of particles in IAU_JUPITER frame at time t1
    lines = file.readlines()  # array of lines (each line represents a particle)
    file.close()
    column = np.size(lines[0].split())
    y0_particles = np.zeros((num_particles, 6))  # matrix of initial conditions for particles

    ## INITIAL POSITION OF PARTICLE AND FOUR GALILEAN MOONS USING SPICE
    ystart = np.zeros((5, 6))  # 5 bodies to simulate : 1 particle + 4 moons

    ## initial conditions in ROTATING frame
    # ystart[0,:] will be initialized later
    ystart[1, :] = np.array(np.append(spice_GetPosition('IO', origin, t1), spice_GetVelocity('IO', origin, t1)))  # initial state of Io
    ystart[2, :] = np.array(np.append(spice_GetPosition('EUROPA', origin, t1), spice_GetVelocity('EUROPA', origin, t1)))  # initial state of Europa
    ystart[3, :] = np.array(np.append(spice_GetPosition('GANYMEDE', origin, t1), spice_GetVelocity('GANYMEDE', origin, t1)))  # initial state of Ganymede
    ystart[4, :] = np.array(np.append(spice_GetPosition('CALLISTO', origin, t1), spice_GetVelocity('CALLISTO', origin, t1)))  # initial state of Callisto

    ## initial conditions in INERTIAL frame (to integrate in inertial frame (without centrifugal and coriolis force))
    # ystart[0,:] = np.concatenate((pos_inert, v_orb_inert)) # initial state vector for particle expressed in IAU_JUPITER frame (body-fixed)
    # ystart[1,:] = np.array(np.append(spice_GetPosition('IO', origin, t1), spice_GetVelocity('IO', origin, t1))) # initial state of Io
    # ystart[2,:] = np.array(np.append(spice_GetPosition('EUROPA', origin, t1), spice_GetVelocity('EUROPA', origin, t1))) # initial state of Europa
    # ystart[3,:] = np.array(np.append(spice_GetPosition('GANYMEDE', origin, t1), spice_GetVelocity('GANYMEDE', origin, t1))) # initial state of Ganymede
    # ystart[3,:] = np.array(np.append(spice_GetPosition('GANYMEDE', origin, t1), spice_GetVelocity('GANYMEDE', origin, t1))) # initial state of Ganymede
    # ystart[4,:] = np.array(np.append(spice_GetPosition('CALLISTO', origin, t1), spice_GetVelocity('CALLISTO', origin, t1))) # initial state of Callisto

    ## determine initial conditions of particles
    for i in range(num_particles):
        data = np.zeros((column))  # initialize
        str_list = lines[i].split()  # splits the 7 numbers of a line into a list of strings
        for j in range(column):
            data[j] = float(str_list[j])  # makes an array with the 7 values
        pos_inert = data[1:4]  # we suppose that all particles are emitted at the center of a given moon
        v_orb_inert = data[4:]
        W = np.array([0., 0., w_jup])  # angular velocity vector
        v_orb_rot = v_orb_inert - np.cross(W, pos_inert)
        pos_rot = spice_inert2rot(np.append(pos_inert, v_orb_inert), t1)[0:3]
        y0_particles[i, :] = np.append(pos_rot, v_orb_rot)

    return y0_particles, ystart


def simulate(y0_particles, ystart):
    tps = np.zeros((num_particles, itermax))
    yps = np.zeros((5, 6, itermax, num_particles))
    accps = np.zeros((8, 3, itermax, num_particles))
    orbitps = np.zeros((3, itermax, num_particles))
    kounts = np.zeros((num_particles))

    # loop over all particles
    t_start = time()
    for i in range(num_particles):
        ystart[0, :] = y0_particles[i, :]  # initial condition for one particle
        tp, yp, accp, orbitp, kount, nok, nbad = driver(ystart, t1, t2)  # trajectory of one particle
        print('')
        print('-- Trajectory of particle %d calculated' % (i + 1))
        print('Actual duration of simulation : %f days' % ((tp[kount - 1] - t1) / (24 * 3600)))
        tps[i, :] = tp
        yps[:, :, :, i] = yp
        accps[:, :, :, i] = accp
        orbitps[:, :, i] = orbitp
        kounts[i] = kount
    t_end = time()
    kounts = kounts.astype(int)

    print('')
    print('Real time of whole simulation : %f minutes' % ((t_end - t_start) / 60))
    return yps, tps, accps, orbitps, kounts


def save_outputs(yps, tps, accps, orbitps, kounts, job_ID, script_ID):
    np.save('output_state_'+str(script_ID)+'_'+str(job_ID), yps)
    np.save('output_time_'+str(script_ID)+'_'+str(job_ID), tps)
    np.save('output_acc_'+str(script_ID)+'_'+str(job_ID), accps)
    np.save('output_orbitelem_'+str(script_ID)+'_'+str(job_ID), orbitps)
    np.save('output_kounts_'+str(script_ID)+'_'+str(job_ID), kounts)









#############################################################################################
# MAIN SCRIPT
#############################################################################################
# Three-body problem: Jupiter, Sun and particle
# state vector is of dimension 6
# all units are expressed in SI units, except distances which are expressed in kilometers


if __name__ == "__main__":
    print('')
    print('Initializing simulation...')

    ## GLOBAL CONSTANTS
    global job_ID, num_particles, metakernel, integrator, frame, traj_flag, origin, tiny, hill_radius, G, m1, m_sun, mu, mu_sun, Rj, a, AU, d, w_sun, w_jup, J2, B_magn, epsilon0, Sdot, Ls, Sdot, c, mass_moons, m_dust, r_dust, beta, density, itermax, kmax, h1, hmin, dtsav, eps, exp_shrink, exp_grow, S, errcon, B, g, h, degree, abs_logN0, ord_logN0, abs_logkT, ord_logkT, N0, r0, alpha, l0, gk, E0, E1, ratio, mh, kmaxx, imax, S1, S2, redmax, redmin, scalmx, radius_moons

    ## ARGUMENTS FROM COMMAND LINE
    args = sys.argv # arguments from the command line
    parameters_filename = args[1]
    ic_file = args[2]
    num_particles = int(args[3]) # number of particles to be simulated
    job_ID = int(args[4])
    script_ID = int(args[5])

    from model_parameters_jupiter import * ## importlib.import_module(filename)
    spice.furnsh(metakernel)

    y0_particles, ystart = load_initial_conditions(ic_file)

    print('Model parameters and initial conditions of particles loaded')

    ## SIMULATION
    #pr = cProfile.Profile()  # profiler
    #pr.enable()
    yps, tps, accps, orbitps, kounts = simulate(y0_particles, ystart)
    #pr.disable()

    ## save results of simulation
    save_outputs(yps, tps, accps, orbitps, kounts, job_ID, script_ID)
    print('')
    print('Results of simulation saved')

    ## profiler
    #pr.dump_stats('profile_'+str(script_ID)+'_'+str(job_ID)+'.cprof')
    #print('Profiling saved')
    print('DONE')

# plot_figs(yps, tps, accps, orbitps, kounts, t1, row, particleID)  # plot trajectories of the particles


