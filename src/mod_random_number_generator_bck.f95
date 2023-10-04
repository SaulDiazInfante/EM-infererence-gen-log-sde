module mod_random_number_generator
  ! use iso_fortran_env, only : int32, real32
  implicit none
contains
!> @brief Returns a random variables with uniform distribution using the
!! standar gfortran random number generator, the returned value is a real 32.
!> @param [in] ix  Dummy paramer for the seed initialization for the random
!> generator 
real function unif_(ix)
!  integer(int32), intent(in):: ix
!  real(real32) unif_, x
!    call random_seed()
!    call random_number(x)
!    unif_ = x
!    return
    print*, ix
  end function unif_
! 
!> @brief Implementation of the Box-muller algorithm to genrate Gaussian random
!>  variables
!> @param [in ] seed  initial value for the uniform random generator
!> @todo: Try other random-number generators, such like mersene an others
!  real(real32) function boxmuller()
!    integer(int32) iset
!    real fac, gset, rsq, v1, v2, seed
!    save iset, gset
!    data iset/0/
!    1 if (iset.eq.0) then
!        call random_number(seed)
!        v1 = 2. * seed - 1.0
!        call random_number(seed)
!        v2 = 2. * seed - 1.0
!        rsq = v1 ** 2 + v2 ** 2
!        if (rsq.ge.1..or.rsq.eq.0)goto 1
!            fac = sqrt(-2.0 * log(rsq) / rsq)
            !gset = v1 * fac
!            boxmuller = v2 * fac
!            iset = 1
!        else
!            boxmuller=gset
!            iset=0
!        endif
!      return
!  end function boxmuller
!> @brief Returns a number with Gaussian distribution using the Box-Muller
!> algortihm.
!> @details This function can be ommited in the case of a step forward 
!>  implelemtation
!> @param [in] seed a int32 with the initial value for the randim uniform 
!> generator
!> @todo: Implement this furnction such that returns a realization path of the
!> standard Brownian motion
!  real(real32) function normalvar()  result(r_x)
!    r_x = boxmuller()
!  end function normalvar
end module mod_random_number_generator
