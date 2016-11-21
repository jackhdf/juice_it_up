## Run this script to start several jobs
## If one wants to do a quick simulation, one can run Adaptive_RK_CK_particle_moons_jupiter_traj.py from the command line
## command line : python super_script.py text_files/Initial_Conditions_Callisto.txt 1000 10 0
## runs super_script.py with 1000 particles from Callisto using 10 jobs (with each 100 particles) script_ID = 0
# using the argument "all" instead of 1000 will simulate all particles in the initial conditions file


import numpy as np
from math import *
import sys
import subprocess as sp
import os

if __name__ == "__main__":

    ## Starts several jobs on command line for a given file of initial conditions
    args = sys.argv  # arguments from the command line
    ic_file = args[1]

    try:
        num_particles = int(args[2])
    except ValueError:
        file = open(ic_file,'r')  # file contains initial positions and velocities of particles in IAU_JUPITER frame at time t1
        lines = file.readlines()  # array of lines (each line represents a particle)
        file.close()
        num_particles = np.size(lines)

    num_jobs = int(args[3]) # number of jobs
    script_ID = int(args[4]) # if we are calling several times super_script.py (for example to simulate initial conditions from different moons), this identifier is useful

    file = open(ic_file,'r')  # file contains initial positions and velocities of particles in IAU_JUPITER frame at time t1
    lines = file.readlines()  # array of lines (each line represents a particle)
    file.close()
    total_num_particles = np.size(lines) # total number of particles in a file of initial conditions
    num_particles_per_job = ceil(num_particles/num_jobs)
    diff = abs(num_particles_per_job*num_jobs - num_particles) # the last job will have diff less particles to simulate than the other jobs

    if num_particles > total_num_particles or num_jobs > num_particles:
        print('Arguments are not valid')
        sys.exit(0) # stops executing code

    print('Script ' + str(script_ID) + ' launched...')

    proc = [None] * num_jobs # empty list of subprocesses
    for job_ID in range(num_jobs):
        if job_ID < (num_jobs-1):
            print('Job ' + str(job_ID) + ' initialized')
            command = ['python', 'Adaptive_RK_CK_particle_moons_jupiter_traj.py', ic_file, str(num_particles_per_job), str(job_ID), str(script_ID), 'multi']
            proc[job_ID] = sp.Popen(command) #, shell = True, stdout = log, stderr = sp.PIPE)
        elif job_ID == (num_jobs-1):
            print('Job ' + str(job_ID) + ' initialized')
            command = ['python', 'Adaptive_RK_CK_particle_moons_jupiter_traj.py', ic_file, str(num_particles_per_job-diff), str(job_ID), str(script_ID), 'multi']
            proc[job_ID] = sp.Popen(command) #, shell = True, stdout = log, stderr = sp.PIPE)

    for i in range(num_jobs):
        proc[i].wait()
    print('All jobs are done')
    os.system('python plot_figures.py ' + str(num_particles) + ' ' + str(num_jobs) + ' ' + str(script_ID) + ' ' + ic_file)
    #os.system(['python', 'plot_figures.py', str(num_particles), str(num_jobs), str(script_ID), ic_file])

    print('')
    print('----------------------------')
    print('ALL DONE')
    print('')
    print('Now, take a cup of coffee, sit back, relax and enjoy the results !')