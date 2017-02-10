######################################################################
###################  MODEL PARAMETERS  ###############################
######################################################################
import spiceypy as spice
from math import *
import numpy as np


## USER PARAMETERS
duration = 20*24*3600 # duration of the simulation in real time [s]
traj_flag = True  # if this flag is True then the trajectories of the particles are plotted in blue
itermax = 100000 # max number of iterations (increase if increase of the simulation time)
particleID = 0  # ID of the particle for which to show energy and acceleration plots
shuffle_particles = 1 # if set to 1, the program will shuffle the particles that will be simulated from the initial conditions file. If set to 0, the program will simulate the particles in the order that they appear in the initial conditions file
show_plots = 1 # 1 if we want to see the plots at the end of the simulation

## INTEGRATOR PARAMETERS
integrator = 'rkck' # 'bs' for Burlisch-Stoer
frame = 'IAU_JUPITER' # when integrating in body-fixed frame
# frame = 'JUICE_JUPITER_IF_J2000' # when integrating in inertial frame (but then plots are not valid)
origin = 'JUPITER' # origin of the system
metakernel = 'text_files/Traj_juice_141a_generic_cfg_linux.txt' # metakernel file name
spice.furnsh(metakernel)

## FORCES TO BE APPLIED: 1 is applied, 0 is not applied
# if all the following parameters are set to 0, then only primary gravity of Jupiter will be considered
central_body_oblateness = 1
sun_gravity = 1
moon_gravity = ['IO', 'EUROPA', 'GANYMEDE', 'CALLISTO'] ## List with the names of the moons IN ORDER OF INCREASING DISTANCE from Jupiter, to be taken into account (should be SPICE compatible names)
sun_radiation_pressure = 1
sun_pr_drag = 1
lorentz = 1 # Here we will have to see how to trigger the use of the B field model and particle potential for Jupiter or Saturnbut for now we should be able to switch on or off this acceleration from the config file - the charge of the particle will be derived from the potential and density values
plasma_drag = 1 # same comment than for Lorentz

## TIME PARAMETERS
#t_init = '2016-01-01T00:00:01.816092' # time in Julian calendar
#t1 = spice.str2et(t_init) # ephemeris time : number of seconds past the J2000 epoch: initial time [s]
t1 = 5.0487847e+08
t2 = t1 + duration # final ephemeris time [s]

## MODEL PARAMETERS
tiny = 1e-30
G = 6.67384e-20  # gravitational constant [km^3 * kg^-1 * s^-2]
m_jup = 1.8986e27  # mass of Jupiter [kg]
m_sun = 1.9891e30  # mass of sun [kg]
mu = G * m_jup
mu_sun = G * m_sun
Rj = 71492.  # radius of Jupiter [km]
a = 71492  # equatorial radius of Jupiter [km]
AU = 1.496e8  # astronomical unit [km]
semimaj = 5.20260 * AU # semi-major axis of orbit of Jupiter [kg]
e = 0.048498 # excentricity of Jupiter's orbit
d = 5 * AU  # distance of the sun from Jupiter [km]
w_sun = (2 * pi) / (4335.3545 * 24 * 3600)  # angular velocity of Jupiter around Sun
w_jup = (2 * pi) / (9 * 3600 + 55 * 60 + 27)  # angular velocity of the rotation of Jupiter on itself
J2 = 0.014736  # Jupiter's J2 value [m^5 . s^-2] but not sure about the units (although it works with this value)
B_magn = 1e-12  # for constante magnetic field [Tesla]
epsilon0 = 8.85418782e-12  # vacuum permittivity [Farad/m]
Sdot = 1361  # average total solar irradiance at 1 AU [Watt/m^2]
Ls = 4 * pi * (d ** 2) * Sdot  # solar luminosity [microWatt = kg*km^2/s^3] at a distance of 5 AU from the Sun (initially solar luminosity is in W/m^2, so we multiply this value by the surface of a sphere defined by a radius of 5 AU)
c = 299792.458  # speed of light [km/s]
mass_moons = np.array([8.9319e22, 4.8e22, 1.4819e23, 1.0759e23])  # mass of four galilean moons in the same order as below [kg]
radius_moons = np.array([1830., 1560.8, 2631.2, 2410.3])  # radius of the moons of Jupiter [km]
num_zones = np.size(mass_moons) + 1  # number of zones
r = 5*Rj # initial distance from Jupiter
hill_radius = semimaj*(1-e)*(m_jup/(3*m_sun))**(1/3) # Hill radius of Jupiter [km]

## PARTICLE PARAMETERS
r_dust = 1e-6  # radius of the particle [m]
density = 1500  # kg/m^3
m_dust = (4. / 3) * (r_dust ** 3) * density  # [kg]
beta = 2.278887e-02  # solar radiation pressure coefficient
particle_potential_fixed = 5 ## in V fixed value if not chosen as model dependent
particle_potential_model_dependent = 0 ## if set to 1, we have to use a model for the potential (a table with the potential as function of the particle position in the system  TBC)

## PLASMA DRAG MODEL PARAMETERS
# each row (5) contains values for a given distance from jupiter: |1.0|3.8|5.5|7.9|20|170| [Rj]. Coefficients (see article, table 6) for O+, O++, S+, S++, S+++, Na+ (rows are for different distances of the particle from Jupiter)
eV = 1.60218e-19  # [1 eV = 1.60218e-19 Joules]
abs_logN0 = np.array([3.8, 4.9, 5.1, 5.3, 5.5, 5.65, 5.8, 5.9, 6.4, 7.4, 7.9, 10, 20, 60, 100, 170])  # Jovicentric distance [Rj]
ord_logN0 = np.array([1.55, 2.75, 2.91, 3.27, 2.88, 3.57, 3.31, 3.35, 3.18, 2.78, 2.25, 1.48, 0.2, -2, -2,-3])  # electron density [cm^-3]
abs_logkT = abs_logN0
ord_logkT = eV * np.array([1.67, -0.31, -0.18, 0.37, 0.92, 1.15, 1.33, 1.54, 1.63, 1.67, 1.75, 2., 2., 2., 2., 2.])  # [Joules]
N0 = 4.65  # N0 value for inner plasmasphere around Jupiter [cm^-3] from table 6 of the article
r0 = np.array([7.68, 20.])  # value for inner plasmasphere region and outer disc [Rj]
alpha = (np.array([7., 7., 10.77]) / 360) * 2 * pi  # inclination of the magnetic axis [radians]
l0 = (21 / 360) * 2 * pi  # reference longitude  [radians] (valid for all five regions around Jupiter)
gk = np.array([[0.2, 0.02, 0.7, 0.03, 0., 0.],  # inner plasmasphere
               [0.2, 0.02, 0.7, 0.03, 0., 0.],  # cool torus
               [0.06, 0.08, 0.24, 0.25, 0.01, 0.01],  # warm torus
               [0.06, 0.08, 0.24, 0.25, 0.01, 0.01],  # inner disc
               [0.07, 0.06, 0.06, 0.26, 0.06, 0.05]])  # outer disc
E0 = np.array([0., 1., 1., 100., 100.]) * eV  # values for the 5 regions [Joules]
E1 = 85 * eV  # valid for the last two outer regions [Joules]
ratio = 0.9  # [degrees/Rj]
uma = 1.660538921e-27  # atomic mass unit [kg]
m_elec = 9.10938356e-31  # mass of one electron [kg]
m_ion1 = 15.9994 * uma - 1 * m_elec  # O+
m_ion2 = 15.9994 * uma - 2 * m_elec  # O++
m_ion3 = 32.065 * uma - 1 * m_elec  # S+
m_ion4 = 32.065 * uma - 2 * m_elec  # S++
m_ion5 = 32.065 * uma - 3 * m_elec  # S+++
m_ion6 = 22.98976928 * uma - 1 * m_elec  # Na+
mh = np.array([m_ion1, m_ion2, m_ion3, m_ion4, m_ion5, m_ion6])  # = ion mass [kg/ion] (formula from wikipedia)

# VIPAL MODEL PARAMETERS (MAGNETIC FIELD)
# matrix of Schmidt coefficients (m=0,...,5 rows and n=1,...,5 columns)
degree = 4  # normally 5, but we are missing associated normalized Legendre polynomials for order 5
g = np.array([[4.2, 0.6441, -0.1058, -0.7466, -0.066],
              [-0.6975, -0.8672, -0.59, 0.3282, 0.0737],
              [0., 0.9598, 0.6322, -0.338, -0.1711],
              [0., 0., 0.4671, 0.1826, -0.1793],
              [0., 0., 0., -0.1429, -0.0077],
              [0., 0., 0., 0., -0.0740]])
# matrix of Schmidt coefficients (m=0,...,5 rows and n=1,...,5 columns)
h = np.array([[0., 0., 0., 0., 0.], # the line with the exponent m=0 does not exist because in the formula of B, we see that sin(m.phi) = 0 if m = 0
              [0.1973, -0.4041, -0.231, 0.3283, 0.2065],
              [0., 0.603, 0.516, -0.2131, -0.1167],
              [0., 0., -0.1131, -0.0606, -0.0288],
              [0., 0., 0., -0.0486, -0.0050],
              [0., 0., 0., 0., -0.2279]])

## NUMERICAL METHOD PARAMETERS
kmax = itermax # max number of steps that can be stored (if kmax = 0 there is no intermediate storage)
h1 = 100 # guessed first stepsize [s]
hmin = 1 # minimum allowed stepsize (can be zero)
dtsav = 1 # results y are only stored at intervals greater than  (not very useful in practice)
kmaxx = 8 # maximum row number used in the extrapolation (used in Burlisch-Stoer) : if kmaxx == 8 then the value of n = 16
imaxx = kmaxx + 1 # used in Burlisch-Stoer
S1 = 0.25 # safety factor (used in Burlisch-Stoer)
S2 = 0.7 # safety factor (used in Burlisch-Stoer)
redmax = 1e-5 # Maximum factor for stepsize reduction (used in Burlisch-Stoer)
redmin = 0.7 # Minimum factor for stepsize reduction (used in Burlisch-Stoer)
scalmx = 0.1 # 1/scalmx is the maximum factor by which a stepsize can be increased (used in Burlisch-Stoer)
eps = 1e-8 # eps is a tolerance level
exp_shrink = - 0.25 # use this exponent when reducing the time step
exp_grow = - 0.2 # use this exponent when increasing the time step
S = 0.9 # safety factor
errcon = 1.89e-4 # equals (5/SAFETY) raised to the power (1/PGROW), see use below.equals (5/SAFETY) raised to the power (1/PGROW), see use below.
B = np.array([[0., 0., 0., 0., 0., 0., 0.],
              [1. / 5, 1. / 5, 0., 0., 0., 0., 0.],
              [3. / 10, 3. / 40, 9. / 40, 0., 0., 0., 0.],
              [3. / 5, 3. / 10, -9. / 10, 6. / 5, 0., 0., 0.],
              [1., -11. / 54, 5. / 2, -70. / 27, 35. / 27, 0., 0.],
              [7. / 8, 1631. / 55296, 175. / 512, 575. / 13824, 44275. / 110592, 253. / 4096, 0.],
              [0., 37. / 378, 0., 250. / 621, 125. / 594, 0., 512. / 1771],
              [0., 2825. / 27648, 0., 18575. / 48384, 13525. / 55296, 277. / 14336, 1. / 4]])  # Butcher tableau for Cash-Karp method (better than RKF)




## EXAMPLES OF INITIAL CONDITIONS FOR ONE PARTICLE
## INITIAL VELOCITY OF PARTICLE WRT IAU_JUPITER FRAME
## circular orbit in equatorial plane with position along x-axis and initial velocity along y-axis
#pos_inert = np.array([r,0.,0.]) # this position is actually in body fixed frame, but since we consider that the inertial frame is the body-fixed frame at t1, then it's the same  ##rot2inert(np.array([r,0.,0.,0.,v_orb_inert,0.]), t1, w_jup)[0:3]
#v_orb_inert = np.array([0., np.sqrt(mu/r), 0.]) # orbital velocity [km/s] in inertial frame (in equatorial plane

## circular orbit in equatorial plane with position along x-axis and initial velocity along y-axis (particle is trying to hit a moon from the back)
## be carful when using initial condition t1-10 for particle then it seems in the inertial frame like the particle is hitting the moon from the front, but actually in the rotating frame it is hitting the moon from the back
#pos_inert = spice_GetPosition('IO', 'JUPITER', t1-50) # this position is actually in body fixed frame, but since we consider that the inertial frame is the body-fixed frame at t1, then it's the same  ##rot2inert(np.array([r,0.,0.,0.,v_orb_inert,0.]), t1, w_jup)[0:3]
#v_orb_rot = spice_GetVelocity('IO', 'JUPITER', t1-50)

## polar orbit with position along z-axis and initial velocity along y-axis
# pos_inert = np.array([0., 0., -r]) # this position is actually in body fixed frame, but since we consider that the inertial frame is the body-fixed frame at t1, then it's the same  ##rot2inert(np.array([r,0.,0.,0.,v_orb_inert,0.]), t1, w_jup)[0:3]
# v_orb_inert = np.array([0., np.sqrt(mu/r), 0.]) # orbital velocity [km/s] in inertial frame (in equatorial plane

## slight inclination
# pos_inert = np.array([r,0.,0.]) # this position is actually in body fixed frame, but since we consider that the inertial frame is the body-fixed frame at t1, then it's the same  ##rot2inert(np.array([r,0.,0.,0.,v_orb_inert,0.]), t1, w_jup)[0:3]
# v_orb_inert = np.sqrt(mu/r) * np.array([0., 0.5, 0.5]) # orbital velocity [km/s] in inertial frame (in equatorial plane


## elliptical orbit in equatorial plane with position along x-axis and initial velocity along y-axis : with J2 gravity term, there should be a precession of orbit
# pos_inert = np.array([r,0.,0.]) # this position is actually in body fixed frame, but since we consider that the inertial frame is the body-fixed frame at t1, then it's the same  ##rot2inert(np.array([r,0.,0.,0.,v_orb_inert,0.]), t1, w_jup)[0:3]
# dist = np.sqrt(pos_inert[0]**2 + pos_inert[1]**2 + pos_inert[2]**2) # initial distance from Jupiter
# v_orb_inert = np.array([0., np.sqrt(mu/dist)/1.7, 50.]) # orbital velocity [km/s] in inertial frame (in equatorial plane

## hyperbolic orbit in equatorial plane
# pos_inert = np.array([r,0.,0.]) # this position is actually in body fixed frame, but since we consider that the inertial frame is the body-fixed frame at t1, then it's the same  ##rot2inert(np.array([r,0.,0.,0.,v_orb_inert,0.]), t1, w_jup)[0:3]
# dist = np.sqrt(pos_inert[0]**2 + pos_inert[1]**2 + pos_inert[2]**2) # initial distance from Jupiter
# escape_vel = np.sqrt((2*mu)/dist) # magnitude of escape velocity
# direction = np.array([-r,3*Rj,0]) / np.sqrt(r**2 + 9*(Rj**2)) # direction of initial velocity (pointing at 2*Rj from center of Jupiter)
# v_orb_inert = escape_vel * direction

## hyperbolic flight in polar orbit about Jupiter
# pos_inert = np.array([0.,0.,-r]) # this position is actually in body fixed frame, but since we consider that the inertial frame is the body-fixed frame at t1, then it's the same  ##rot2inert(np.array([r,0.,0.,0.,v_orb_inert,0.]), t1, w_jup)[0:3]
# dist = np.sqrt(pos_inert[0]**2 + pos_inert[1]**2 + pos_inert[2]**2) # initial distance from Jupiter
# escape_vel = np.sqrt((2*mu)/dist) # magnitude of escape velocity
# direction = np.array([3*Rj,0,r]) / np.sqrt(r**2 + 9*(Rj**2)) # direction of initial velocity (pointing at 2*Rj from center of Jupiter)
# v_orb_inert = escape_vel * direction

## test case : 1 micron diameter of particle, 1e-15 kg for the mass of particle, with initial position at one moon radius of Io on the leading side (longitude = 90 degrees)
##  (ie ejected from the surface of Io) with magnitude of initial velocity equal to circular velocity of Io + escape velocity at the surface of Io, in the direction of the circular velocity of Io
# origin = 'JUPITER'
# moon = 'IO'
# state_moon = np.array(np.append(spice_GetPosition(moon, origin, t1), spice_GetVelocity(moon, origin, t1))) # initial state of Io
# vel_moon = np.sqrt(state_moon[3]**2 + state_moon[4]**2 + state_moon[5]**2)
# radius_moon = spice.bodvrd(moon, "RADII", 3)[1][1]
# mass_moon = mass_moons[0] #spice.bodvrd(moon, "MU",1)
# escape_vel_moon = np.sqrt((2*G*mass_moon)/radius_moon) # magnitude of escape velocity at surface of Io
# mag_vel = vel_moon + escape_vel_moon # magnitude of initial velocity of particle
# direction_vel = state_moon[3:] / vel_moon # direction of initial velocity of particle
# pos_inert = state_moon[0:3] + radius_moon*direction_vel
# v_orb_rot = mag_vel * direction_vel

## transformation to body-fixed frame
#W = np.array([0.,0.,w_jup]) # angular velocity vector
#v_orb_rot = v_orb_inert - np.cross(W,pos_inert)
#initial_state = spice_inert2rot(np.concatenate((pos_inert,v_orb_inert)),t1) # the initial velocity makes it a retrograde orbit which causes the particle to exit the system
