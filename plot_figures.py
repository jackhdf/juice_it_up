## Plots figures
## command line : python3 plot_figures.py 1000 10 0 text_files/Initial_Conditions_Callisto.txt
## this command plots 1000 particules that were simulated on 10 jobs using script_ID = 0
# using the argument "all" instead of 1000 will plot all particles in the initial conditions file


from math import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import spiceypy as spice
import sys
import os

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

# plot figures in Jupiter body-fixed frame
def plot_figs(yps, tps, accps, orbitps, kounts, t1, row, particleID, script_ID):
    # print('Plotting...')

    # os.makedirs('Figures/')

    ## Converting position and velocity of particle from IAU_JUPITER to JIF
    yps_inert = np.zeros((5, 6, itermax, row))  # initialization
    for k in range(np.size(yps, axis=3)):  # over the particles
        for i in range(kounts[k]):  # over the iterations
            for j in range(np.size(yps, axis=0)):  # over the particle and moons
                yps_inert[j, :, i, k] = spice_rot2inert(yps[j, :, i, k], tps[k, i])

    # translates distances into Jovian radii
    #yps = yps / Rj
    #yps_inert = yps_inert / Rj

    ## sphere and disc for each moon of Jupiter
    # theta = np.linspace(0, 2 * pi, 100)
    # disc_0_x = Rj*np.cos(theta)  # Jupiter
    # disc_0_y = Rj*np.sin(theta)
    # disc_1_x = radius_moons[0] * np.cos(theta)  # Io
    # disc_1_y = radius_moons[0] * np.sin(theta)
    # disc_2_x = radius_moons[1] * np.cos(theta)  # Europa
    # disc_2_y = radius_moons[1] * np.sin(theta)
    # disc_3_x = radius_moons[2] * np.cos(theta)  # Ganymede
    # disc_3_y = radius_moons[2] * np.sin(theta)
    # disc_4_x = radius_moons[3] * np.cos(theta)  # Callisto
    # disc_4_y = radius_moons[3] * np.sin(theta)
    #
    # phi, theta = np.mgrid[0:pi:25j, 0:2 * pi:25j]
    # sphere_0_x = Rj*np.sin(phi) * np.cos(theta)  # Jupiter
    # sphere_0_y = Rj*np.sin(phi) * np.sin(theta)
    # sphere_0_z = Rj*np.cos(phi)
    # sphere_1_x = radius_moons[0] * np.sin(phi) * np.cos(theta)
    # sphere_1_y = radius_moons[0] * np.sin(phi) * np.sin(theta)
    # sphere_1_z = radius_moons[0] * np.cos(phi)
    # sphere_2_x = radius_moons[1] * np.sin(phi) * np.cos(theta)
    # sphere_2_y = radius_moons[1] * np.sin(phi) * np.sin(theta)
    # sphere_2_z = radius_moons[1] * np.cos(phi)
    # sphere_3_x = radius_moons[2] * np.sin(phi) * np.cos(theta)
    # sphere_3_y = radius_moons[2] * np.sin(phi) * np.sin(theta)
    # sphere_3_z = radius_moons[2] * np.cos(phi)
    # sphere_4_x = radius_moons[3] * np.sin(phi) * np.cos(theta)
    # sphere_4_y = radius_moons[3] * np.sin(phi) * np.sin(theta)
    # sphere_4_z = radius_moons[3] * np.cos(phi)
    theta = np.linspace(0, 2 * pi, 100)
    disc_0_x = np.cos(theta)  # Jupiter
    disc_0_y = np.sin(theta)
    disc_1_x = (radius_moons[0] / Rj) * np.cos(theta)  # Io
    disc_1_y = (radius_moons[0] / Rj) * np.sin(theta)
    disc_2_x = (radius_moons[1] / Rj) * np.cos(theta)  # Europa
    disc_2_y = (radius_moons[1] / Rj) * np.sin(theta)
    disc_3_x = (radius_moons[2] / Rj) * np.cos(theta)  # Ganymede
    disc_3_y = (radius_moons[2] / Rj) * np.sin(theta)
    disc_4_x = (radius_moons[3] / Rj) * np.cos(theta)  # Callisto
    disc_4_y = (radius_moons[3] / Rj) * np.sin(theta)

    phi, theta = np.mgrid[0:pi:25j, 0:2 * pi:25j]
    sphere_0_x = np.sin(phi) * np.cos(theta)  # Jupiter
    sphere_0_y = np.sin(phi) * np.sin(theta)
    sphere_0_z = np.cos(phi)
    sphere_1_x = (radius_moons[0] / Rj) * np.sin(phi) * np.cos(theta)
    sphere_1_y = (radius_moons[0] / Rj) * np.sin(phi) * np.sin(theta)
    sphere_1_z = (radius_moons[0] / Rj) * np.cos(phi)
    sphere_2_x = (radius_moons[1] / Rj) * np.sin(phi) * np.cos(theta)
    sphere_2_y = (radius_moons[1] / Rj) * np.sin(phi) * np.sin(theta)
    sphere_2_z = (radius_moons[1] / Rj) * np.cos(phi)
    sphere_3_x = (radius_moons[2] / Rj) * np.sin(phi) * np.cos(theta)
    sphere_3_y = (radius_moons[2] / Rj) * np.sin(phi) * np.sin(theta)
    sphere_3_z = (radius_moons[2] / Rj) * np.cos(phi)
    sphere_4_x = (radius_moons[3] / Rj) * np.sin(phi) * np.cos(theta)
    sphere_4_y = (radius_moons[3] / Rj) * np.sin(phi) * np.sin(theta)
    sphere_4_z = (radius_moons[3] / Rj) * np.cos(phi)

    ## 3D plot of trajectory (INERTIAL)
    fig1 = plt.figure(1)
    axes = Axes3D(fig1)
    ## plot and scatters of initial points for particles
    for i in range(row):
        if traj_flag:
            axes.plot(yps_inert[0, 0, 0:kounts[i], i], yps_inert[0, 1, 0:kounts[i], i], yps_inert[0, 2, 0:kounts[i], i],
                      c='b')
        axes.scatter(yps_inert[0, 0, 0, i], yps_inert[0, 1, 0, i], yps_inert[0, 2, 0, i], c='b')
        axes.scatter(yps_inert[0, 0, kounts[i] - 1, i], yps_inert[0, 1, kounts[i] - 1, i],
                     yps_inert[0, 2, kounts[i] - 1, i], c='r')

    ## plot of the trajectories of moons (we take the number of iterations kounts[0] of the first particle for the plot of the trajectories of the moons)
    axes.plot(yps_inert[1, 0, 0:kounts[0], 0], yps_inert[1, 1, 0:kounts[0], 0], yps_inert[1, 2, 0:kounts[0], 0], c='g')
    axes.plot(yps_inert[2, 0, 0:kounts[0], 0], yps_inert[2, 1, 0:kounts[0], 0], yps_inert[2, 2, 0:kounts[0], 0], c='r')
    axes.plot(yps_inert[3, 0, 0:kounts[0], 0], yps_inert[3, 1, 0:kounts[0], 0], yps_inert[3, 2, 0:kounts[0], 0], c='c')
    axes.plot(yps_inert[4, 0, 0:kounts[0], 0], yps_inert[4, 1, 0:kounts[0], 0], yps_inert[4, 2, 0:kounts[0], 0], c='m')

    ## final points of the trajectories of moons
    axes.scatter(yps_inert[1, 0, kounts[0] - 1, 0], yps_inert[1, 1, kounts[0] - 1, 0],
                 yps_inert[1, 2, kounts[0] - 1, 0], c='k')
    axes.scatter(yps_inert[2, 0, kounts[0] - 1, 0], yps_inert[2, 1, kounts[0] - 1, 0],
                 yps_inert[2, 2, kounts[0] - 1, 0], c='k')
    axes.scatter(yps_inert[3, 0, kounts[0] - 1, 0], yps_inert[3, 1, kounts[0] - 1, 0],
                 yps_inert[3, 2, kounts[0] - 1, 0], c='k')
    axes.scatter(yps_inert[4, 0, kounts[0] - 1, 0], yps_inert[4, 1, kounts[0] - 1, 0],
                 yps_inert[4, 2, kounts[0] - 1, 0], c='k')
    axes.scatter(0., 0., 0., c='k')
    axes.scatter(yps_inert[1, 0, 0, 0], yps_inert[1, 1, 0, 0], yps_inert[1, 2, 0, 0], c='b')
    axes.scatter(yps_inert[2, 0, 0, 0], yps_inert[2, 1, 0, 0], yps_inert[2, 2, 0, 0], c='b')
    axes.scatter(yps_inert[3, 0, 0, 0], yps_inert[3, 1, 0, 0], yps_inert[3, 2, 0, 0], c='b')
    axes.scatter(yps_inert[4, 0, 0, 0], yps_inert[4, 1, 0, 0], yps_inert[4, 2, 0, 0], c='b')

    axes.plot_wireframe(sphere_0_x, sphere_0_y, sphere_0_z, color='k')
    axes.plot_wireframe(yps_inert[1, 0, 0, 0] + sphere_1_x, yps_inert[1, 1, 0, 0] + sphere_1_y,
                        yps_inert[1, 2, 0, 0] + sphere_1_z, color='k')  # kounts[0]-1
    axes.plot_wireframe(yps_inert[2, 0, 0, 0] + sphere_2_x, yps_inert[2, 1, 0, 0] + sphere_2_y,
                        yps_inert[2, 2, 0, 0] + sphere_2_z, color='k')
    axes.plot_wireframe(yps_inert[3, 0, 0, 0] + sphere_3_x, yps_inert[3, 1, 0, 0] + sphere_3_y,
                        yps_inert[3, 2, 0, 0] + sphere_3_z, color='k')
    axes.plot_wireframe(yps_inert[4, 0, 0, 0] + sphere_4_x, yps_inert[4, 1, 0, 0] + sphere_4_y,
                        yps_inert[4, 2, 0, 0] + sphere_4_z, color='k')

    axes.set_xlabel('X-axis [Rj]')
    axes.set_ylabel('Y-axis [Rj]')
    axes.set_zlabel('Z-axis [Rj]')
    axes.set_label('3D plot of trajectory (INERTIAL FRAME)')
    axes.set_aspect('equal', 'datalim')
    fig1.savefig('Figures/traj_3d_inertial_script_' + str(script_ID) + '.png')

    ## 3D plot of trajectory (ROTATING)
    # fig2 = plt.figure(2)
    # axes = Axes3D(fig2)
    # axes.plot(yp[0,0, 0:kount], yp[0,1, 0:kount], yp[0,2, 0:kount], c='b')
    # axes.plot(yp[1,0, 0:kount], yp[1,1, 0:kount], yp[1,2, 0:kount], c='g')
    # axes.plot(yp[2,0, 0:kount], yp[2,1, 0:kount], yp[2,2, 0:kount], c='r')
    # axes.plot(yp[3,0, 0:kount], yp[3,1, 0:kount], yp[3,2, 0:kount], c='c')
    # axes.plot(yp[4,0, 0:kount], yp[4,1, 0:kount], yp[4,2, 0:kount], c='m')
    # axes.scatter(yp[0,0, kount-1], yp[0,1, 0], yp[0,2, kount-1], c='k')
    # axes.scatter(yp[1,0, kount-1], yp[1,1, 0], yp[1,2, kount-1], c='k')
    # axes.scatter(yp[2,0, kount-1], yp[2,1, 0], yp[2,2, kount-1], c='k')
    # axes.scatter(yp[3,0, kount-1], yp[3,1, 0], yp[3,2, kount-1], c='k')
    # axes.scatter(yp[4,0, kount-1], yp[4,1, 0], yp[4,2, kount-1], c='k')
    # axes.scatter(0., 0., 0., c='k')
    # axes.set_xlabel('X-axis [Rj]')
    # axes.set_ylabel('Y-axis [Rj]')
    # axes.set_zlabel('Z-axis [Rj]')
    # axes.set_label('3D plot of trajectory (ROTATING FRAME)')
    # axes.set_aspect('equal', 'datalim')


    ## Projection of trajectory on X-Y plane (ROTATING)
    # fig3 = plt.figure(3)
    # ## plot and scatters of initial points for particles
    # if traj_flag:
    #     for i in range(row):
    #         plt.plot(yps[0, 0, 0:kounts[i], i], yps[0, 1, 0:kounts[i], i], c='b')
    #         plt.scatter(yps[0, 0, 0, i], yps[0, 1, 0, i], c='b')
    #         #plt.scatter(yps[0, 0, kounts[i] - 1, i], yps[0, 1, kounts[i] - 1, i], c='k')
    #
    # ## plot of the trajectories of moons (we take the number of iterations kounts[0] of the first particle for the plot of the trajectories of the moons)
    # plt.plot(yps[1, 0, 0:kounts[0], 0], yps[1, 1, 0:kounts[0], 0], c='g')
    # plt.plot(yps[2, 0, 0:kounts[0], 0], yps[2, 1, 0:kounts[0], 0], c='r')
    # plt.plot(yps[3, 0, 0:kounts[0], 0], yps[3, 1, 0:kounts[0], 0], c='c')
    # plt.plot(yps[4, 0, 0:kounts[0], 0], yps[4, 1, 0:kounts[0], 0], c='m')
    #
    # ## final points of the trajectories of moons
    # plt.scatter(yps[1, 0, kounts[i] - 1, 0], yps[1, 1, kounts[i] - 1, 0], c='k')
    # plt.scatter(yps[2, 0, kounts[i] - 1, 0], yps[2, 1, kounts[i] - 1, 0], c='k')
    # plt.scatter(yps[3, 0, kounts[i] - 1, 0], yps[3, 1, kounts[i] - 1, 0], c='k')
    # plt.scatter(yps[4, 0, kounts[i] - 1, 0], yps[4, 1, kounts[i] - 1, 0], c='k')
    # plt.scatter(0., 0., c='k')
    # plt.scatter(yps[1, 0, 0, 0], yps[1, 1, 0, 0], c='b')
    # plt.scatter(yps[2, 0, 0, 0], yps[2, 1, 0, 0], c='b')
    # plt.scatter(yps[3, 0, 0, 0], yps[3, 1, 0, 0], c='b')
    # plt.scatter(yps[4, 0, 0, 0], yps[4, 1, 0, 0], c='b')
    #
    # plt.plot(disc_0_x, disc_0_y, c='k')  # Jupiter
    # plt.plot(yps[1, 0, kounts[i] - 1, 0] + disc_1_x, yps[1, 1, kounts[i] - 1, 0] + disc_1_y, c='k')
    # plt.plot(yps[2, 0, kounts[i] - 1, 0] + disc_2_x, yps[2, 1, kounts[i] - 1, 0] + disc_2_y, c='k')
    # plt.plot(yps[3, 0, kounts[i] - 1, 0] + disc_3_x, yps[3, 1, kounts[i] - 1, 0] + disc_3_y, c='k')
    # plt.plot(yps[4, 0, kounts[i] - 1, 0] + disc_4_x, yps[4, 1, kounts[i] - 1, 0] + disc_4_y, c='k')
    # plt.xlabel('X-axis [Rj]')
    # plt.ylabel('Y-axis [Rj]')
    # plt.title('Projection of trajectory on X-Y plane (ROTATING FRAME)')
    # plt.axes().set_aspect('equal', 'datalim')


    ## Projection of trajectory on X-Y plane (INERTIAL)
    fig4 = plt.figure(4)
    ## plot and scatters of initial points for particles
    for i in range(row):
        if traj_flag:
            plt.plot(yps_inert[0, 0, 0:kounts[i], i], yps_inert[0, 1, 0:kounts[i], i], c='b')
        plt.scatter(yps_inert[0, 0, 0, i], yps_inert[0, 1, 0, i], c='b')
        plt.scatter(yps_inert[0, 0, kounts[i] - 1, i], yps_inert[0, 1, kounts[i] - 1, i], c='r')

    ## plot of the trajectories of moons (we take the number of iterations kounts[0] of the first particle for the plot of the trajectories of the moons)
    plt.plot(yps_inert[1, 0, 0:kounts[0], 0], yps_inert[1, 1, 0:kounts[0], 0], c='g')
    plt.plot(yps_inert[2, 0, 0:kounts[0], 0], yps_inert[2, 1, 0:kounts[0], 0], c='r')
    plt.plot(yps_inert[3, 0, 0:kounts[0], 0], yps_inert[3, 1, 0:kounts[0], 0], c='c')
    plt.plot(yps_inert[4, 0, 0:kounts[0], 0], yps_inert[4, 1, 0:kounts[0], 0], c='m')

    ## final points of the trajectories of moons
    plt.scatter(yps_inert[1, 0, kounts[0] - 1, 0], yps_inert[1, 1, kounts[0] - 1, 0], c='k')
    plt.scatter(yps_inert[2, 0, kounts[0] - 1, 0], yps_inert[2, 1, kounts[0] - 1, 0], c='k')
    plt.scatter(yps_inert[3, 0, kounts[0] - 1, 0], yps_inert[3, 1, kounts[0] - 1, 0], c='k')
    plt.scatter(yps_inert[4, 0, kounts[0] - 1, 0], yps_inert[4, 1, kounts[0] - 1, 0], c='k')
    plt.scatter(0., 0., c='k')
    plt.scatter(yps_inert[1, 0, 0, 0], yps_inert[1, 1, 0, 0], c='b')
    plt.scatter(yps_inert[2, 0, 0, 0], yps_inert[2, 1, 0, 0], c='b')
    plt.scatter(yps_inert[3, 0, 0, 0], yps_inert[3, 1, 0, 0], c='b')
    plt.scatter(yps_inert[4, 0, 0, 0], yps_inert[4, 1, 0, 0], c='b')

    plt.plot(disc_0_x, disc_0_y, c='k')  # Jupiter
    plt.plot(yps_inert[1, 0, 0, 0] + disc_1_x, yps_inert[1, 1, 0, 0] + disc_1_y, c='k')
    plt.plot(yps_inert[2, 0, 0, 0] + disc_2_x, yps_inert[2, 1, 0, 0] + disc_2_y, c='k')
    plt.plot(yps_inert[3, 0, 0, 0] + disc_3_x, yps_inert[3, 1, 0, 0] + disc_3_y, c='k')
    plt.plot(yps_inert[4, 0, 0, 0] + disc_4_x, yps_inert[4, 1, 0, 0] + disc_4_y, c='k')
    plt.title('Projection of trajectory on X-Y plane (INERTIAL FRAME)')
    plt.xlabel('X-axis [Rj]')
    plt.ylabel('Y-axis [Rj]')
    plt.axes().set_aspect('equal', 'datalim')
    fig4.savefig('Figures/traj_xy_inertial_script_' + str(script_ID) + '.png')

    ## Projection of trajectory on X-Z plane (INERTIAL)
    fig7 = plt.figure(7)
    for i in range(row):
        if traj_flag:
            plt.plot(yps_inert[0, 0, 0:kounts[i], i], yps_inert[0, 2, 0:kounts[i], i], c='b')
        plt.scatter(yps_inert[0, 0, 0, i], yps_inert[0, 2, 0, i], c='b')
        plt.scatter(yps_inert[0, 0, kounts[i] - 1, i], yps_inert[0, 2, kounts[i] - 1, i], c='r')

    ## plot of the trajectories of moons (we take the number of iterations kounts[0] of the first particle for the plot of the trajectories of the moons)
    plt.plot(yps_inert[1, 0, 0:kounts[0], 0], yps_inert[1, 2, 0:kounts[0], 0], c='g')
    plt.plot(yps_inert[2, 0, 0:kounts[0], 0], yps_inert[2, 2, 0:kounts[0], 0], c='r')
    plt.plot(yps_inert[3, 0, 0:kounts[0], 0], yps_inert[3, 2, 0:kounts[0], 0], c='c')
    plt.plot(yps_inert[4, 0, 0:kounts[0], 0], yps_inert[4, 2, 0:kounts[0], 0], c='m')

    ## final points of the trajectories of moons
    plt.scatter(yps_inert[1, 0, kounts[0] - 1, 0], yps_inert[1, 2, kounts[0] - 1, 0], c='k')
    plt.scatter(yps_inert[2, 0, kounts[0] - 1, 0], yps_inert[2, 2, kounts[0] - 1, 0], c='k')
    plt.scatter(yps_inert[3, 0, kounts[0] - 1, 0], yps_inert[3, 2, kounts[0] - 1, 0], c='k')
    plt.scatter(yps_inert[4, 0, kounts[0] - 1, 0], yps_inert[4, 2, kounts[0] - 1, 0], c='k')
    plt.scatter(0., 0., c='k')
    plt.scatter(yps_inert[1, 0, 0, 0], yps_inert[1, 2, 0, 0], c='b')
    plt.scatter(yps_inert[2, 0, 0, 0], yps_inert[2, 2, 0, 0], c='b')
    plt.scatter(yps_inert[3, 0, 0, 0], yps_inert[3, 2, 0, 0], c='b')
    plt.scatter(yps_inert[4, 0, 0, 0], yps_inert[4, 2, 0, 0], c='b')

    plt.plot(disc_0_x, disc_0_y, c='k')  # Jupiter
    plt.plot(yps_inert[1, 0, 0, 0] + disc_1_x, yps_inert[1, 2, 0, 0] + disc_1_y, c='k')
    plt.plot(yps_inert[2, 0, 0, 0] + disc_2_x, yps_inert[2, 2, 0, 0] + disc_2_y, c='k')
    plt.plot(yps_inert[3, 0, 0, 0] + disc_3_x, yps_inert[3, 2, 0, 0] + disc_3_y, c='k')
    plt.plot(yps_inert[4, 0, 0, 0] + disc_4_x, yps_inert[4, 2, 0, 0] + disc_4_y, c='k')

    plt.xlabel('X-axis [Rj]')
    plt.ylabel('Z-axis [Rj]')
    plt.title('Projection of trajectory on X-Z plane (INERTIAL FRAME)')
    plt.axes().set_aspect('equal', 'datalim')
    fig7.savefig('Figures/traj_xz_inertial_script_' + str(script_ID) + '.png')

    ## Projection of trajectory on Y-Z plane (INERTIAL)
    fig8 = plt.figure(8)
    for i in range(row):
        if traj_flag:
            plt.plot(yps_inert[0, 1, 0:kounts[i], i], yps_inert[0, 2, 0:kounts[i], i], c='b')
        plt.scatter(yps_inert[0, 1, 0, i], yps_inert[0, 2, 0, i], c='b')
        plt.scatter(yps_inert[0, 1, kounts[i] - 1, i], yps_inert[0, 2, kounts[i] - 1, i], c='r')

    ## plot of the trajectories of moons (we take the number of iterations kounts[0] of the first particle for the plot of the trajectories of the moons)
    plt.plot(yps_inert[1, 1, 0:kounts[0], 0], yps_inert[1, 2, 0:kounts[0], 0], c='g')
    plt.plot(yps_inert[2, 1, 0:kounts[0], 0], yps_inert[2, 2, 0:kounts[0], 0], c='r')
    plt.plot(yps_inert[3, 1, 0:kounts[0], 0], yps_inert[3, 2, 0:kounts[0], 0], c='c')
    plt.plot(yps_inert[4, 1, 0:kounts[0], 0], yps_inert[4, 2, 0:kounts[0], 0], c='m')

    ## final points of the trajectories of moons
    plt.scatter(yps_inert[1, 1, kounts[0] - 1, 0], yps_inert[1, 2, kounts[0] - 1, 0], c='k')
    plt.scatter(yps_inert[2, 1, kounts[0] - 1, 0], yps_inert[2, 2, kounts[0] - 1, 0], c='k')
    plt.scatter(yps_inert[3, 1, kounts[0] - 1, 0], yps_inert[3, 2, kounts[0] - 1, 0], c='k')
    plt.scatter(yps_inert[4, 1, kounts[0] - 1, 0], yps_inert[4, 2, kounts[0] - 1, 0], c='k')
    plt.scatter(0., 0., c='k')
    plt.scatter(yps_inert[1, 1, 0, 0], yps_inert[1, 2, 0, 0], c='b')
    plt.scatter(yps_inert[2, 1, 0, 0], yps_inert[2, 2, 0, 0], c='b')
    plt.scatter(yps_inert[3, 1, 0, 0], yps_inert[3, 2, 0, 0], c='b')
    plt.scatter(yps_inert[4, 1, 0, 0], yps_inert[4, 2, 0, 0], c='b')

    plt.plot(disc_0_x, disc_0_y, c='k')  # Jupiter
    plt.plot(yps_inert[1, 1, 0, 0] + disc_1_x, yps_inert[1, 2, 0, 0] + disc_1_y, c='k')
    plt.plot(yps_inert[2, 1, 0, 0] + disc_2_x, yps_inert[2, 2, 0, 0] + disc_2_y, c='k')
    plt.plot(yps_inert[3, 1, 0, 0] + disc_3_x, yps_inert[3, 2, 0, 0] + disc_3_y, c='k')
    plt.plot(yps_inert[4, 1, 0, 0] + disc_4_x, yps_inert[4, 2, 0, 0] + disc_4_y, c='k')

    plt.xlabel('Y-axis [Rj]')
    plt.ylabel('Z-axis [Rj]')
    plt.title('Projection of trajectory on Y-Z plane (INERTIAL FRAME)')
    plt.axes().set_aspect('equal', 'datalim')
    fig8.savefig('Figures/traj_yz_inertial_script_' + str(script_ID) + '.png')

    # ## Projection of trajectory on X-Z plane (ROTATING)
    # fig9 = plt.figure(9)
    # plt.plot(yp[0, 0, 0:kount], yp[0, 2, 0:kount], c='b')
    # plt.plot(yp[1, 0, 0:kount], yp[1, 2, 0:kount], c='g')
    # plt.plot(yp[2, 0, 0:kount], yp[2, 2, 0:kount], c='r')
    # plt.plot(yp[3, 0, 0:kount], yp[3, 2, 0:kount], c='c')
    # plt.plot(yp[4, 0, 0:kount], yp[4, 2, 0:kount], c='m')
    # plt.scatter(yp[0, 0, 0], yp[0, 2, 0], c='k')
    # plt.scatter(yp[1, 0, 0], yp[1, 2, 0], c='k')
    # plt.scatter(yp[2, 0, 0], yp[2, 2, 0], c='k')
    # plt.scatter(yp[3, 0, 0], yp[3, 2, 0], c='k')
    # plt.scatter(yp[4, 0, 0], yp[4, 2, 0], c='k')
    # plt.scatter(0., 0., 0., c='c')
    # plt.xlabel('X-axis [Rj]')
    # plt.ylabel('Y-axis [Rj]')
    # plt.title('Projection of trajectory on X-Z plane (ROTATING FRAME)')
    # plt.axes().set_aspect('equal', 'datalim')
    #
    # ## Projection of trajectory on Y-Z plane (ROTATING)
    # fig10 = plt.figure(10)
    # plt.plot(yp[0, 1, 0:kount], yp[0, 2, 0:kount], c='b')
    # plt.plot(yp[1, 1, 0:kount], yp[1, 2, 0:kount], c='g')
    # plt.plot(yp[2, 1, 0:kount], yp[2, 2, 0:kount], c='r')
    # plt.plot(yp[3, 1, 0:kount], yp[3, 2, 0:kount], c='c')
    # plt.plot(yp[4, 1, 0:kount], yp[4, 2, 0:kount], c='m')
    # plt.scatter(yp[0, 1, kount-1], yp[0, 2, kount-1], c='k')
    # plt.scatter(yp[1, 1, kount-1], yp[1, 2, kount-1], c='k')
    # plt.scatter(yp[2, 1, kount-1], yp[2, 2, kount-1], c='k')
    # plt.scatter(yp[3, 1, kount-1], yp[3, 2, kount-1], c='k')
    # plt.scatter(yp[4, 1, kount-1], yp[4, 2, kount-1], c='k')
    # plt.scatter(0., 0., 0., c='k')
    # plt.xlabel('X-axis [Rj]')
    # plt.ylabel('Y-axis [Rj]')
    # plt.title('Projection of trajectory on Y-Z plane (ROTATING FRAME)')
    # plt.axes().set_aspect('equal', 'datalim')



    # Energy of particle
    #yps_inert = yps_inert * Rj * 1000
    distance = np.sqrt(yps_inert[0, 0, 0:kounts[particleID], particleID] ** 2 + yps_inert[0, 1, 0:kounts[particleID],
                                                                                particleID] ** 2 + yps_inert[0, 2,
                                                                                                   0:kounts[particleID],
                                                                                                   particleID] ** 2)  # distance of particle from center of the system
    v2_norm = np.sqrt(yps_inert[0, 3, 0:kounts[particleID], particleID] ** 2 + yps_inert[0, 4, 0:kounts[particleID],
                                                                               particleID] ** 2 + yps_inert[0, 5,
                                                                                                  0:kounts[particleID],
                                                                                                  particleID] ** 2)  # speed
    kinetic = (v2_norm ** 2) / 2
    potential = - mu / distance
    energy_tot = kinetic + potential
    # theory
    # kinetic_th = mu/(2*r) * np.ones((np.size(energy_tot)))
    # potential_th = -mu/(r) * np.ones((np.size(energy_tot)))
    # energy_tot_th = -mu/(2*r) * np.ones((np.size(energy_tot)))

    fig5 = plt.figure(5)
    handle1, = plt.plot((tps[particleID, 0:kounts[particleID]] - t1) / 3600, energy_tot, label='Total energy', c='r')
    handle2, = plt.plot((tps[particleID, 0:kounts[particleID]] - t1) / 3600, kinetic, label='Kinetic energy', c='b')
    handle3, = plt.plot((tps[particleID, 0:kounts[particleID]] - t1) / 3600, potential, label='Potential energy', c='g')
    # handle4, = plt.plot((tp[1:kount + 1] - t1) / 3600, energy_tot_th, label='Total energy (theory)', c='r', linestyle = 'dashed')
    # handle5, = plt.plot((tp[1:kount + 1] - t1) / 3600, kinetic_th, label='Kinetic energy (theory)', c='b', linestyle = 'dashed')
    # handle6, = plt.plot((tp[1:kount + 1] - t1) / 3600, potential_th, label='Potential energy (theory)', c='g', linestyle = 'dashed')
    plt.legend(handles=[handle1, handle2, handle3])
    plt.xlabel('Time [hours]')
    plt.ylabel('Energy [J]')
    plt.title('Energy of particle')
    fig5.savefig('Figures/energy_script_' + str(script_ID) + '.png')

    ## Plot of accelerations due to different forces as a function of time
    acc_mag = np.zeros((8, kounts[particleID]))
    for i in range(kounts[particleID]):
        for j in range(8):
            acc_mag[j, i] = sqrt(accps[j, 0, i, particleID] ** 2 + accps[j, 1, i, particleID] ** 2 + accps[
                j, 2, i, particleID] ** 2)  # ag, ags, ace, aco, alo, apr, apd, amoons

    fig6 = plt.figure(6)
    handle1, = plt.plot((tps[particleID, 0:kounts[particleID]] - t1) / 3600, acc_mag[0, :], label='Gravity of Jupiter',
                        c='r')
    handle2, = plt.plot((tps[particleID, 0:kounts[particleID]] - t1) / 3600, acc_mag[1, :],
                        label='Gravity of Sun (including SRP)', c='b')
    # handle3, = plt.plot((tps[particleID, 0:kounts[particleID]] - t1) / 3600, acc_mag[2,:], label='Centrifugal force (fictitious)', c='k')
    # handle4, = plt.plot((tps[particleID, 0:kounts[particleID]] - t1) / 3600, acc_mag[3,:], label='Coriolis force (fictitious)', c='k')
    handle5, = plt.plot((tps[particleID, 0:kounts[particleID]] - t1) / 3600, acc_mag[4, :], label='Lorentz force',
                        c='m')
    handle6, = plt.plot((tps[particleID, 0:kounts[particleID]] - t1) / 3600, acc_mag[5, :],
                        label='Poynting-Robertson drag', c='y')
    handle7, = plt.plot((tps[particleID, 0:kounts[particleID]] - t1) / 3600, acc_mag[6, :], label='Plasma drag', c='g')
    handle8, = plt.plot((tps[particleID, 0:kounts[particleID]] - t1) / 3600, acc_mag[7, :], label='Gravity of moons',
                        c='c')
    plt.legend(handles=[handle1, handle2, handle5, handle6, handle7, handle8], loc='best')
    plt.xlabel('Time [hours]')
    plt.ylabel('Acceleration [km/s^2]')
    plt.title('Different accelerations on the particle')
    fig6.savefig('Figures/accelerations_script_' + str(script_ID) + '.png')

    # ## plot of semi-major axis
    # fig11 = plt.figure(11)
    # plt.plot((tps[particleID, 0:kounts[particleID]] - t1) / 3600, orbitps[0, 0:kounts[particleID], particleID])
    # plt.xlabel('Time [hours]')
    # plt.ylabel('semi-major axis of particle [km]')
    # plt.title('Semi-major axis of particle')
    # fig11.savefig('Figures/semi_major_axis_vs_time_script_' + str(script_ID) + '.png')
    #
    # ## plot of inclination of trajectory
    # fig12 = plt.figure(12)
    # plt.plot((tps[particleID, 0:kounts[particleID]] - t1) / 3600, orbitps[1, 0:kounts[particleID], particleID])
    # plt.xlabel('Time [hours]')
    # plt.ylabel('Inclination of trajectory of particle [radians]')
    # plt.title('Inclination of trajectory of particle ')
    # fig12.savefig('Figures/inclination_vs_time_script_' + str(script_ID) + '.png')

    ## plot of the three components for a given acceleration vector as a function of distanc
    # distance_sorted = sorted(distance)
    # sort_index = np.argsort(distance) # returns the indices of the sorted distance vector
    # accp_sorted = np.zeros((8,3,kmax))
    # for i in range(kounts[particleID]):
    #     accp_sorted[:,:,i] = accps[:,:,sort_index[i],particleID]
    #
    # index = 0 # indicates which acceleration we are considering : ag, ags, ace, aco, alo, apr, apd, amoons
    # fig13 = plt.figure(13)
    # handle1, = plt.plot(distance_sorted, accp_sorted[index, 0, 0:kounts[particleID]], label='x-component of acceleration vector', c='r')
    # handle2, = plt.plot(distance_sorted, accp_sorted[index, 1, 0:kounts[particleID]], label='y-component of acceleration vector', c='b')
    # handle3, = plt.plot(distance_sorted, accp_sorted[index, 2, 0:kounts[particleID]], label='z-component of acceleration vector', c='g')
    # handle4, = plt.plot(distance_sorted, np.sqrt(accp_sorted[index, 0, 0:kounts[particleID]]**2 + accp_sorted[index, 1, 0:kounts[particleID]]**2 +accp_sorted[index, 2, 0:kounts[particleID]]**2), label='magnitude of acceleration vector', c='k')
    # plt.legend(handles=[handle1, handle2, handle3, handle4], loc = 'best')
    # plt.xlabel('Distance [km]')
    # plt.ylabel('Acceleration [km/s^2]')
    # plt.title('Components for a given acceleration vector')
    #
    # index = 4  # indicates which acceleration we are considering : ag, ags, ace, aco, alo, apr, apd, amoons
    # fig14 = plt.figure(14)
    # handle1, = plt.plot(distance_sorted, accp_sorted[index, 0, 0:kounts[particleID]], label='x-component of acceleration vector', c='r')
    # handle2, = plt.plot(distance_sorted, accp_sorted[index, 1, 0:kounts[particleID]], label='y-component of acceleration vector', c='b')
    # handle3, = plt.plot(distance_sorted, accp_sorted[index, 2, 0:kounts[particleID]], label='z-component of acceleration vector', c='g')
    # handle4, = plt.plot(distance_sorted, np.sqrt(accp_sorted[index, 0, 0:kounts[particleID]]**2 + accp_sorted[index, 1, 0:kounts[particleID]]**2 +accp_sorted[index, 2, 0:kounts[particleID]]**2), label='magnitude of acceleration vector', c='k')
    # plt.legend(handles=[handle1, handle2, handle3, handle4], loc='best')
    # plt.xlabel('Distance [km]')
    # plt.ylabel('Acceleration [km/s^2]')
    # plt.title('Components for a given acceleration vector')

    # ## plot of the apocenter as a function of the distance
    # apocenter = 2*orbitps[0,:,particleID] - orbitps[2,:,particleID]
    # #apocenter_sorted = np.zeros((1,kmax))
    # #for i in range(kount):
    # #    apocenter_sorted[i] = apocenter[sort_index[i]]
    # fig15 = plt.figure(15)
    # plt.plot((tps[particleID, 0:kounts[particleID]] - t1) / 3600, apocenter[0:kounts[particleID]])
    # plt.xlabel('Time [hours]')
    # plt.ylabel('Apocenter distance [km]')
    # plt.title('Apocenter distance')

    if show_plots == 1:
        plt.ion()
        plt.show()
        plt.pause(0.001)
        input("Press [enter] to close plots and exit.")


if __name__ == "__main__":
    print("Plotting...")
    from model_parameters_jupiter import *  ## importlib.import_module(filename)
    spice.furnsh(metakernel)

    args = sys.argv  # arguments from the command line
    ic_file = args[4]

    try:
        total_num_particles = int(args[1])
    except ValueError:
        file = open(ic_file,'r')  # file contains initial positions and velocities of particles in IAU_JUPITER frame at time t1
        lines = file.readlines()  # array of lines (each line represents a particle)
        file.close()
        total_num_particles = np.size(lines)

    num_jobs = int(args[2])  # number of jobs used in the simulation
    script_ID = int(args[3])  # if we are calling several times super_script.py (for example to simulate initial conditions from different moons), this identifier is useful

    output_acc = np.load('output_acc_'+str(script_ID)+'_'+str(0)+'.npy')
    output_kounts = np.load('output_kounts_' + str(script_ID) + '_' + str(0) + '.npy')
    output_orbitelem = np.load('output_orbitelem_' + str(script_ID) + '_' + str(0) + '.npy')
    output_state = np.load('output_state_' + str(script_ID) + '_' + str(0) + '.npy')
    output_time = np.load('output_time_' + str(script_ID) + '_' + str(0) + '.npy')

    ## Collate all the output numpy arrays from all the jobs into one numpy array, so that plot_figures.py can plot results for a given script_ID (so for a given set of initial conditions
    for job_ID in range(num_jobs-1):
        temp = np.load('output_acc_'+str(script_ID)+'_'+str(job_ID+1)+'.npy')
        output_acc = np.concatenate((output_acc, temp),axis=3)

        temp = np.load('output_kounts_' + str(script_ID) + '_' + str(job_ID + 1) + '.npy')
        output_kounts = np.concatenate((output_kounts, temp), axis=0)

        temp = np.load('output_orbitelem_' + str(script_ID) + '_' + str(job_ID + 1) + '.npy')
        output_orbitelem = np.concatenate((output_orbitelem, temp), axis=2)

        temp = np.load('output_state_' + str(script_ID) + '_' + str(job_ID + 1) + '.npy')
        output_state = np.concatenate((output_state, temp), axis=3)

        temp = np.load('output_time_' + str(script_ID) + '_' + str(job_ID + 1) + '.npy')
        output_time = np.concatenate((output_time, temp), axis=0)

    #print('Output files regrouped')

    plot_figs(output_state, output_time, output_acc, output_orbitelem, output_kounts, t1, total_num_particles, particleID, script_ID)
    print('Figures from simulation ' + str(script_ID) + ' saved')
