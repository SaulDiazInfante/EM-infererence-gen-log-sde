program mod_montecarlo_path_sampler_test
    use mod_brownian_motion
    use mod_path_sampler
    use ifport
    use iso_fortran_env, only: &
            stdin => input_unit, &
            stdout => output_unit, &
            stderr => error_unit
    
    implicit none
    integer :: n_montecarlo, n_resolution, r_operative_factor, user_seed
    integer ::  date_time(8)
    character*10 b(3)
    character (len = 100) :: file_name
    character(len=8) :: fmt
    logical :: debug
    real(kind=8), allocatable :: brownian_path(:), sample_path(:), mc_sample(:, :)
    integer(kind=4) :: n_operative, i
    real(kind=8) h_op, h_res, time_horizon
    real(kind=8) :: x_0, alpha, m, sigma
!
    time_horizon = 32.0
    n_resolution = 4096
    n_montecarlo = 1000
    r_operative_factor = 64
    n_operative = 0
    h_op = 0.0
    user_seed=123456
    debug = .false.
    x_0 = 0.1
    alpha = 0.8
    m = 1.5
    sigma = 0.1
    debug = .false.
!   termalization phase
    call get_brownian_path(&
            &n_resolution, &
            &r_operative_factor, &
            &time_horizon, &
            &brownian_path, &
            &h_op, &
            &n_operative, &
            &user_seed, &
            &debug &
            &)
    allocate(mc_sample(n_montecarlo, n_operative + 1))
!    mc_sample(:, :) = 0.0
!    print *, "head(X_mc)", mc_sample(1:5, 1:10)
    call date_and_time(b(1), b(2), b(3), date_time)
    write (file_name, &
            &"('../data/output_mc_', &
            &I4.4, '-', &
            &I2.2, '-', &
            &I2.2, '_', &
            &I2.2, ':', &
            &I2.2, ':'  &
            &I2.2, ':', &
            &I3.3, &
            &'.csv') &
            ") date_time(1), date_time(2), date_time(3),&
            &date_time(5), date_time(6), date_time(7), date_time(8)
    open (unit=6, carriagecontrol='fortran')
    do i=1, n_montecarlo
        user_seed = i
        debug = .false.
        call get_brownian_path(&
            &n_resolution, &
            &r_operative_factor, &
            &time_horizon, &
            &brownian_path, &
            &h_op, &
            &n_operative, &
            &user_seed, &
            &debug &
        &)
! To tag the current path
!        print *, '(===)', i, ' ', trim(file_name)
       ! sample_path = [(0.0, i=0, n_operative)]
        call get_sample_path(&
            &brownian_path, &
            &n_operative, &
            &h_op, &
            &x_0, &
            &alpha, &
            &m, &
            &sigma, &
            &sample_path, &
            &debug&!, &
           ! &file_name&
        &)
        mc_sample(i, :) = sample_path
!        print *, "(===):", i
!        print *, "head(X_mc[i, :] )", mc_sample(i, 1:10)
        call progress(i)
    end  do
    call save_mc_csv(mc_sample, h_op, n_montecarlo, n_operative, file_name)
! TODO:write a function to save a csv
! write(1, fmt=csv_fmt) 0, sep, 0.0, sep, x_0
end program mod_montecarlo_path_sampler_test

subroutine progress(j)
    implicit none
    integer, intent(in) :: j
    integer::i, k
    character(len=60) :: bar
    bar="???%[....................................................]"
    i = 0
    do k=0, j
        if  (MOD(k, 20)  == 0) then
            i = i+1
            write(unit=bar(1:3),fmt="(i3)") (2 * i - 2)
            bar(5 + i: 5+i+1)= "#>"
        end if
    enddo
    ! print the progress bar.
    write(unit=6, fmt="(a1, a1, a60)") '+',char(60), bar
    return
end subroutine progress

subroutine save_mc_csv(mc_sample, h_op, n_montecarlo, n_operative, filename)
    integer(kind=4), intent(in) :: n_montecarlo, n_operative
    real(kind=8), intent(in) ::  mc_sample(n_montecarlo, n_operative), h_op
    character (len = 100), intent(in), optional :: filename
    integer(kind=4) :: i, j, fileunit
    character (len = 100) :: default_filename, header
    character(len=:), allocatable :: csv_fmt, header_fmt
    character(len=1)  :: sep = ','
    default_filename = '../data/mc_sample_.csv'
    header_fmt ='(a4, a, a4, a, a12, a, a12)'
    csv_fmt = '(i4, a, i4, a, f12.8, a, f12.8)'
    header = 'idx' // sep // 'j' // sep  // 't_j' // sep // 'X(t_{j})'
    if(present(filename)) then
        open(newunit=fileunit, file=trim(filename))
        print *, filename
    else
        open(newunit=fileunit, file=trim(default_filename))
    end if
    write(fileunit, header_fmt) header
    do i=1, n_montecarlo
        do j=1, n_operative + 1
            write(fileunit, fmt=csv_fmt) i, sep, (j - 1), sep, (j - 1) * h_op, sep, mc_sample(i, j)
        end do
    end do
end subroutine  save_mc_csv