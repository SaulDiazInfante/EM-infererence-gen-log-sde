mod_path_sampler.mod = $(src/mod_path_sampler.f90)
$(src/mod_path_sampler.f90) += $(mod_stochastic_logistic_model.mod)
$(src/mod_path_sampler.f90) += $(mod_brownian_motion.mod)
$(src/mod_path_sampler.f90) += $(iso_fortran_env.mod)
