rm -f sim_ode.mod ../bin/sim_ode
gfortran sim_ode.f90 dlsode.f -o ../bin/sim_ode
