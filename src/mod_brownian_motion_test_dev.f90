program mod_brownian_motion_test
    use mod_brownian_motion, only: get_brownian_path
    implicit none
    integer(kind=4) :: n_resolution, r_operative_factor
    real(kind=8), allocatable :: brownian_path(:)
    integer(kind=4) :: n_operative, i, seed_get, seed_put
    real(kind=8) h_op, h_res, xi, browinian_operative_inc, &
        browinian_resolution_inc, time_horizon
    time_horizon = 16.0
    n_resolution = 2048
    r_operative_factor = 1
    n_operative = 2048
    get_brownian_path(&
        n_resolution, &
        r_operative_factor, &
        time_horizon, &
        123456, &
        brownian_path &
    )
    print*,"Hello Brownian path"
end program mod_brownian_motion_test
