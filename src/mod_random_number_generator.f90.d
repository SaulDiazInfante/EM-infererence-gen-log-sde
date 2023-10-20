mod_random_number_generator.mod = $(src/mod_random_number_generator.f90)
$(src/mod_random_number_generator.f90) += $(mkl_vsl_type.mod)
$(src/mod_random_number_generator.f90) += $(mkl_vsl.mod)
$(src/mod_random_number_generator.f90) += $(iso_fortran_env.mod)
$(src/mod_random_number_generator.f90) += mkl_vsl.f90
