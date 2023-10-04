!> @ingroup modules
!> @author E. Lince-Gomez, F. Baltazar-Larios, S. Diaz-Infante
!> @brief This module contains functions to handle various aspects of 
!> one-dimensional Brownian motion. 
!> @see Higham & Kloeden 2021
module mod_brownian_motion_sampler
  use iso_fortran_env, only : int32, real32
  use stdlib_random, only: random_seed
  use stdlib_stats_distribution_normal, only: norm => rvs_normal
  implicit none
contains
!> @brief Returns a realization of a standard Brownian motion given a l
!> homogenoun stenci with N operational points 
!! standar gfortran random number generator, the returned value is a real 32.
!> @param [in] ix  Dummy paramer for the seed initialization for the random
!> generator 
    subroutine get_brownian_path(n_resolution, r_operative_factor,&
        time_horizon, brownian_path)
        implicit none
        integer(int32), intent(in) :: n_resolution, r_operative_factor
        real(real32), intent(in) :: time_horizon
        real(real32), intent(out) :: &
            brownian_path(r_operative_factor * n_resolution)
        integer(int32) :: n_operative, i, seed_get, seed_put
        real(real32) h_op, h_res, xi, browinian_operative_inc, &
            browinian_resolution_inc
        ! random number initialization
        seed_put = 1234567
        call random_seed(seed_put, seed_get)
        brownian_path(:) = 0.0
        xi = 0.0
        if (r_operative_factor > 1 ) then 
            h_op = r_operative_factor * h_res
            n_operative = n_resolution / r_operative_factor
        else if (r_operative_factor == 1) then 
            h_op = h_res
            n_operative = n_resolution
        end if
        
        do i=2, n_resolution
            xi = norm()
            browinian_resolution_inc = sqrt(h_res) * xi
            brownian_path(i) =  brownian_path(i -1 ) + browinian_resolution_inc
        end do
    end subroutine get_brownian_path
end module mod_random_number_generator
