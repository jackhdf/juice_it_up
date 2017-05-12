The following code simulates the motion of a dust particle in the Jovian system under the influence of the following forces:

1) Jupiter gravity with J2 term
2) Sun gravity
3) Solar radiation pressure
4) Lorentz force
5) Poynting-Robertson drag
5) Plasma drag
6) Four galilean moons' gravity

The integrator is a Runge-Kutta Cash-Karp scheme with adaptive step. The integrator can go forward and backward in time. 
We integrate in the IAU_JUPITER frame (body-fixed frame). The kernel pool is described in the file "Traj_juice_141a_generic_cfg_linux.txt"

The file "model_parameters_jupiter.py" is used as a configuration file for the simulation. All physical parameters of the model are present as well
as the forces that the user chooses to apply. The user can also change the initial time, the duration of the simulation as well as which particle in 
the initial conditions file to simulate.

To run a simple simulation of one particle with the parameters specified in the above file, one can type the following command:


python3 Adaptive_RK_CK_particle_moons_jupiter_traj.py ./text_files/Initial_Conditions_test.txt 1 0 0 solo

All results are plotted by default and data from the simulation is also saved for later use.
All the forces, the integrator and the initialization of the simulation are done in "Adaptive_RK_CK_particle_moons_jupiter_traj.py"
If one wants to simulate 1000 particles on 10 simulataneous jobs specified for example in "Initial_Conditions_Callisto.txt", one can use the following command:


python3 super_script.py text_files/Initial_Conditions_Callisto.txt 1000 10 0
