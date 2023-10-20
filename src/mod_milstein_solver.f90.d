mod_milstein_solver.mod = $(src/mod_milstein_solver.f90)
$(src/mod_milstein_solver.f90) += $(iso_fortran_env.mod)
$(src/mod_milstein_solver.f90) += $(mod_random_number_generator.mod)
