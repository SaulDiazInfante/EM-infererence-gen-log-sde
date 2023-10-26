program mod_montecarlo_path_sampler_test
    use mod_brownian_motion
    use mod_path_sampler
    implicit none
    integer :: n_montecarlo, n_resolution, r_operative_factor, user_seed
    integer ::  date_time(8)
    character*10 b(3)
    character (len = 100) :: file_name
    character(len=8) :: fmt
    logical :: debug
    real(kind=8), allocatable :: brownian_path(:), sample_path(:)
    integer(kind=4) :: n_operative, i
    real(kind=8) h_op, h_res, time_horizon
    real(kind=8) :: x_0, alpha, m, sigma
!
    time_horizon = 32.0
    n_resolution = 4096
    n_montecarlo = 10
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
    !># TODO: make a sampler of 10000 ptahs
    do i=1, n_montecarlo
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
        print *, '(===)', i, ' ', trim(file_name)
        sample_path = [(0.0, i=0, n_operative)]
        call get_sample_path(&
            &brownian_path, &
            &n_operative, &
            &h_op, &
            &x_0, &
            &alpha, &
            &m, &
            &sigma, &
            &sample_path, &
            &debug, &
            &file_name&
        &)
    end  do
end program mod_montecarlo_path_sampler_test
