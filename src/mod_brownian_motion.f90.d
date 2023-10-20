mod_brownian_motion.mod = $(src/mod_brownian_motion.f90)
$(src/mod_brownian_motion.f90) += $(mod_random_number_generator.mod)
$(src/mod_brownian_motion.f90) += $(iso_fortran_env.mod)
