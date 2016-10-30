## creats a movie of the trajectories of the particles
## command line : python movie_maker.py 20 1 0
## makes a movie from 20 particles simulated with 1 job using script_ID = 0

from math import *
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
from mpl_toolkits.mplot3d import Axes3D
import spiceypy as spice
import sys


def spice_rot2inert(y,t):
    frame_from = 'IAU_JUPITER'
    frame_to = 'JUICE_JUPITER_IF_J2000' # 'JUICE_JSM'
    et = t

    rot_mat = spice.pxform(frame_from, frame_to, et)  # rotation matrix from IAU_JUPITER frame to JIF
    w = np.array([0.,0.,w_jup])
    transfo_mat = spice.rav2xf(rot_mat,-w)
    yout = np.dot(transfo_mat,y)
    return yout



def collate(script_ID, num_jobs):
    output_acc = np.load('output_acc_' + str(script_ID) + '_' + str(0) + '.npy')
    output_kounts = np.load('output_kounts_' + str(script_ID) + '_' + str(0) + '.npy')
    output_orbitelem = np.load('output_orbitelem_' + str(script_ID) + '_' + str(0) + '.npy')
    output_state = np.load('output_state_' + str(script_ID) + '_' + str(0) + '.npy')
    output_time = np.load('output_time_' + str(script_ID) + '_' + str(0) + '.npy')

    ## Collate all the output numpy arrays from all the jobs into one numpy array, so that plot_figures.py can plot results for a given script_ID (so for a given set of initial conditions
    for job_ID in range(num_jobs - 1):
        temp = np.load('output_acc_' + str(script_ID) + '_' + str(job_ID + 1) + '.npy')
        output_acc = np.concatenate((output_acc, temp), axis=3)

        temp = np.load('output_kounts_' + str(script_ID) + '_' + str(job_ID + 1) + '.npy')
        output_kounts = np.concatenate((output_kounts, temp), axis=0)

        temp = np.load('output_orbitelem_' + str(script_ID) + '_' + str(job_ID + 1) + '.npy')
        output_orbitelem = np.concatenate((output_orbitelem, temp), axis=2)

        temp = np.load('output_state_' + str(script_ID) + '_' + str(job_ID + 1) + '.npy')
        output_state = np.concatenate((output_state, temp), axis=3)

        temp = np.load('output_time_' + str(script_ID) + '_' + str(job_ID + 1) + '.npy')
        output_time = np.concatenate((output_time, temp), axis=0)

    print('Output files regrouped')
    return output_acc, output_kounts, output_orbitelem, output_state, output_time


if __name__ == "__main__":
    from model_parameters_jupiter import *  ## importlib.import_module(filename)
    spice.furnsh(metakernel)

    args = sys.argv  # arguments from the command line
    total_num_particles = int(args[1])
    num_jobs = int(args[2])  # number of jobs used in the simulation
    script_ID = int(args[3])  # if we are calling several times super_script.py (for example to simulate initial conditions from different moons), this identifier is useful

    # collate results from all the jobs together
    output_acc, output_kounts, output_orbitelem, output_state, output_time = collate(script_ID, num_jobs)

    # UTILISER PYTABLES POUR CREER LES VIDEOS PLUS RAPIDEMENT

    ## create movie
    FFMpegWriter = manimation.writers['ffmpeg']
    metadata = dict(title='Movie Test', artist='Matplotlib', comment='Movie support!')
    writer = FFMpegWriter(fps=15, metadata=metadata)
    number_of_frames = max(output_kounts) # equals the maximum number of iterations for all the particles
    print(str(number_of_frames) + ' frames to render')

    ## Projection of trajectory on X-Y plane (INERTIAL)
    fig = plt.figure()
    ## Converting position and velocity of particle from IAU_JUPITER to JIF
    yps_inert = np.zeros((5, 6, itermax, total_num_particles))  # initialization
    for k in range(np.size(output_state, axis=3)):  # over the particles
        for i in range(output_kounts[k]):  # over the iterations
            for j in range(np.size(output_state, axis=0)):  # over the particle and moons
                yps_inert[j, :, i, k] = spice_rot2inert(output_state[j, :, i, k], output_time[k, i])

    ## sphere and disc for each moon of Jupiter
    theta = np.linspace(0, 2 * pi, 100)
    disc_0_x = Rj * np.cos(theta)  # Jupiter
    disc_0_y = Rj * np.sin(theta)
    disc_1_x = radius_moons[0] * np.cos(theta)  # Io
    disc_1_y = radius_moons[0] * np.sin(theta)
    disc_2_x = radius_moons[1] * np.cos(theta)  # Europa
    disc_2_y = radius_moons[1] * np.sin(theta)
    disc_3_x = radius_moons[2] * np.cos(theta)  # Ganymede
    disc_3_y = radius_moons[2] * np.sin(theta)
    disc_4_x = radius_moons[3] * np.cos(theta)  # Callisto
    disc_4_y = radius_moons[3] * np.sin(theta)

    phi, theta = np.mgrid[0:pi:25j, 0:2 * pi:25j]
    sphere_0_x = Rj * np.sin(phi) * np.cos(theta)  # Jupiter
    sphere_0_y = Rj * np.sin(phi) * np.sin(theta)
    sphere_0_z = Rj * np.cos(phi)
    sphere_1_x = radius_moons[0] * np.sin(phi) * np.cos(theta)  # Jupiter
    sphere_1_y = radius_moons[0] * np.sin(phi) * np.sin(theta)
    sphere_1_z = radius_moons[0] * np.cos(phi)
    sphere_2_x = radius_moons[1] * np.sin(phi) * np.cos(theta)  # Jupiter
    sphere_2_y = radius_moons[1] * np.sin(phi) * np.sin(theta)
    sphere_2_z = radius_moons[1] * np.cos(phi)
    sphere_3_x = radius_moons[2] * np.sin(phi) * np.cos(theta)  # Jupiter
    sphere_3_y = radius_moons[2] * np.sin(phi) * np.sin(theta)
    sphere_3_z = radius_moons[2] * np.cos(phi)
    sphere_4_x = radius_moons[3] * np.sin(phi) * np.cos(theta)  # Jupiter
    sphere_4_y = radius_moons[3] * np.sin(phi) * np.sin(theta)
    sphere_4_z = radius_moons[3] * np.cos(phi)

    ## plot of the trajectories of moons (we take the number of iterations kounts[0] of the first particle for the plot of the trajectories of the moons)
    plt.plot(yps_inert[1, 0, 0:output_kounts[0], 0], yps_inert[1, 1, 0:output_kounts[0], 0], c='g')
    plt.plot(yps_inert[2, 0, 0:output_kounts[0], 0], yps_inert[2, 1, 0:output_kounts[0], 0], c='r')
    plt.plot(yps_inert[3, 0, 0:output_kounts[0], 0], yps_inert[3, 1, 0:output_kounts[0], 0], c='c')
    plt.plot(yps_inert[4, 0, 0:output_kounts[0], 0], yps_inert[4, 1, 0:output_kounts[0], 0], c='m')

    ## final points of the trajectories of moons
    plt.scatter(yps_inert[1, 0, output_kounts[0] - 1, 0], yps_inert[1, 1, output_kounts[0] - 1, 0], c='k')
    plt.scatter(yps_inert[2, 0, output_kounts[0] - 1, 0], yps_inert[2, 1, output_kounts[0] - 1, 0], c='k')
    plt.scatter(yps_inert[3, 0, output_kounts[0] - 1, 0], yps_inert[3, 1, output_kounts[0] - 1, 0], c='k')
    plt.scatter(yps_inert[4, 0, output_kounts[0] - 1, 0], yps_inert[4, 1, output_kounts[0] - 1, 0], c='k')
    plt.scatter(0., 0., c='k')
    plt.scatter(yps_inert[1, 0, 0, 0], yps_inert[1, 1, 0, 0], c='b')
    plt.scatter(yps_inert[2, 0, 0, 0], yps_inert[2, 1, 0, 0], c='b')
    plt.scatter(yps_inert[3, 0, 0, 0], yps_inert[3, 1, 0, 0], c='b')
    plt.scatter(yps_inert[4, 0, 0, 0], yps_inert[4, 1, 0, 0], c='b')

    plt.plot(disc_0_x, disc_0_y, c='k')  # Jupiter
    plt.plot(yps_inert[1, 0, output_kounts[0] - 1, 0] + disc_1_x, yps_inert[1, 1, output_kounts[0] - 1, 0] + disc_1_y, c='k')
    plt.plot(yps_inert[2, 0, output_kounts[0] - 1, 0] + disc_2_x, yps_inert[2, 1, output_kounts[0] - 1, 0] + disc_2_y, c='k')
    plt.plot(yps_inert[3, 0, output_kounts[0] - 1, 0] + disc_3_x, yps_inert[3, 1, output_kounts[0] - 1, 0] + disc_3_y, c='k')
    plt.plot(yps_inert[4, 0, output_kounts[0] - 1, 0] + disc_4_x, yps_inert[4, 1, output_kounts[0] - 1, 0] + disc_4_y, c='k')

    plt.title('Projection of trajectory on X-Y plane (INERTIAL FRAME)')
    plt.xlabel('X-axis [km]')
    plt.ylabel('Y-axis [km]')
    plt.axes().set_aspect('equal', 'datalim')
    plt.xlim(-hill_radius, hill_radius)
    plt.ylim(-hill_radius, hill_radius)

    ## x-y plot
    print('')
    print('Rendering x-y plot...')
    with writer.saving(fig, "movie_xy.mp4", number_of_frames):
        for i in range(number_of_frames):
            for j in range(total_num_particles):
                plt.scatter(yps_inert[0, 0, i, j], yps_inert[0, 1, i, j], c='b')
            writer.grab_frame()
            if (i % 10) == 0:
                print(str(i) + ' frames rendered...')

    print('')
    print('Video rendered')

    ## x-z plot
    print('')
    print('Rendering x-z plot...')
    with writer.saving(fig, "movie_xz.mp4", number_of_frames):
        for i in range(number_of_frames):
            for j in range(total_num_particles):
                plt.scatter(yps_inert[0, 0, i, j], yps_inert[0, 1, i, j], c='b')
            writer.grab_frame()
            if (i % 10) == 0:
                print(str(i) + ' frames rendered...')

    print('')
    print('Video rendered')

    # y-z plot
    print('')
    print('Rendering y-z plot...')
    with writer.saving(fig, "movie_yz.mp4", number_of_frames):
        for i in range(number_of_frames):
            for j in range(total_num_particles):
                plt.scatter(yps_inert[0, 0, i, j], yps_inert[0, 1, i, j], c='b')
            writer.grab_frame()
            if (i % 10) == 0:
                print(str(i) + ' frames rendered...')

    print('')
    print('All videos rendered')