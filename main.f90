!> @author
!> Sam Hatfield, AOPP, University of Oxford
!> @brief
!> An observing system simulation experiment (OSSE) using the Lorenz '63 model
!> and 4DVar.
!> Based on code by Amos Lawless, University of Reading.
program lorenz96_4dvar
    use params
    use lorenz96, only: run_model
    use utils, only: time_seed, randn
    use io, only: output
    use assim, only: calc_cost, calc_cost_grad

    implicit none

    real(dp), dimension(tstep,n_x) :: truth = 0.0_dp
    real(dp), dimension(tstep,n_x) :: best_guess
    real(dp) :: obs(tstep/freq,n_x)
    real(dp) :: diagn(max_iterations,1)
    real(dp) :: l(n_x), f, norm, initial(n_x)
    real(dp) :: time(tstep)
    integer :: i, j = 1

    ! Check whether observation frequency divides into total number of timsteps
    if (mod(tstep, freq) /= 0) then
        stop 'Number of timesteps must be a multiple of the observation frequency'
    end if

    ! Seed random number generator
    call time_seed

    ! Generate time axis
    time = (/ (real(i)*h, i = 0, tstep-1) /)

    ! Run truth
    truth = run_model(tstep, (/ (randn(1.0_dp, 0.1_dp), i = 1, n_x) /))

    ! Calculate observations
    do i = 1, tstep, freq
        last = i
        do j = 1, n_x
            obs(1+i/freq,j) = randn(truth(i,j), sqrt(obs_var))
        end do
    end do

    ! Output truth and observations
    call output(time, truth, "truth.txt")
    call output(time, obs, "obs.txt", freq)

    ! Set initial best guess
    initial = (/ (randn(1.0_dp, 5.0_dp), i = 1, n_x) /)

    ! Perform minimisation
    do j = 1, max_iterations
        ! Compute cost of current best guess
        best_guess = run_model(tstep, initial)
        diagn(j,1) = calc_cost(tstep, best_guess, obs)

        ! Output first guess
        if (j == 1) then
            call output(time, best_guess, "first_guess.txt")
        end if

        ! Compute gradient of cost function
        l = calc_cost_grad(tstep, best_guess, obs)

        ! Compute norm of cost function gradient
        norm = sqrt(sum(l**2))

        ! Normalise gradient vector
        l = l/norm

        ! Update initial state estimate at beginning of window
        initial = initial - 0.5_dp*l
    end do

    ! Output final best guess
    call output(time, best_guess, "final_guess.txt")

    ! Output diagnostics
    call output((/ (real(i,dp), i = 1, max_iterations) /), diagn, "diagnostics.txt")
end program lorenz96_4dvar
