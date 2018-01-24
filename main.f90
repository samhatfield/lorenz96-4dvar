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

    real(dp), dimension(n_x,tstep) :: truth = 0.0_dp
    real(dp), dimension(n_x,tstep) :: best_guess
    real(dp) :: obs(n_x,tstep/freq)
    real(dp) :: cost, diagn(1,max_iterations)
    real(dp) :: grad(n_x), f, norm, initial(n_x)
    real(dp) :: time(tstep)
    integer :: i, j, iters

    ! Dummy variables for gradient descent algorithm
    real(dp) :: d(n_x), grad_old(n_x), w(n_x)

    ! Flags to control gradient descent algorithm behaviour
    integer :: printflags(2) = (/ -1, 2 /), flag = 0, rest = 0, method = 3

    ! Variables required by gradient descent subroutine call, but that aren't actually used
    real(dp) :: eps = 1.0d-5
    logical :: finish = .false.

    ! Check whether observation frequency divides into total number of timsteps
    if (mod(tstep, freq) /= 0) then
        stop 'Number of timesteps must be a multiple of the observation frequency'
    end if

    ! Seed random number generator
    call time_seed

    ! Generate time axis
    time = (/ (real(i)*h, i = 0, tstep-1) /)

    ! Spin up truth
    truth(:,1) = (/ (randn(0.0_dp, 5.0_dp), i = 1, n_x) /)
    do i = 1, spin_up
        truth(:,:1) = run_model(1, truth(:,1))
    end do

    ! Run truth
    truth = run_model(tstep, truth(:,1))

    ! Calculate observations
    do i = 1, tstep, freq
        last = i
        do j = 1, n_x
            obs(j, 1+i/freq) = randn(truth(j,i), sqrt(obs_var))
        end do
    end do

    ! Output truth and observations
    call output(time, truth, "truth.txt")
    call output(time, obs, "obs.txt", freq)

    ! Set initial best guess
    initial = (/ (truth(i,1) + randn(0.0_dp, 1.0_dp), i = 1, n_x ) /)

    ! Perform minimisation
    iters = 1
    do
        ! Compute cost of current best guess
        if (flag == 0 .or. flag == 1) then
            best_guess = run_model(tstep, initial)
            cost = calc_cost(tstep, best_guess, obs)
            diagn(1,iters) = cost
    
            ! Output first guess
            if (iters == 0) then
                call output(time, best_guess, "first_guess.txt")
            end if
    
            ! Compute gradient of cost function
            grad = calc_cost_grad(tstep, best_guess, obs)
        end if

        ! Use gradient descent algorithm to step towards minimum of cost function
        call cgfam(n_x, initial, cost, grad, d, grad_old, printflags, eps, w, flag, rest, method, finish)
        if (flag <= 0 .or. iters >= max_iterations) exit
        if (flag == 1) iters = iters + 1
        if (flag == 2) then
            if (cost < 0.2_dp * n_obs * obs_var * n_x) then
                finish = .true.
            end if
        end if
    end do

    write (*,'(A5,I5,A11)') 'Took ', iters, ' iterations'
    write (*,'(A11,F9.2)') 'Final cost ', cost

    ! Output final best guess
    call output(time, best_guess, "final_guess.txt")

    ! Output diagnostics
    call output((/ (real(i,dp), i = 1, max_iterations) /), diagn, "diagnostics.txt")
end program lorenz96_4dvar
