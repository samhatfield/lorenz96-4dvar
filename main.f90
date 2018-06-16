!> @author
!> Sam Hatfield, AOPP, University of Oxford
!> @brief
!> An observing system simulation experiment (OSSE) using the Lorenz '96 model
!> and 4DVar.
!> Based on code by Amos Lawless, University of Reading.
program lorenz96_4dvar
    use params
    use lorenz96, only: run_model
    use utils, only: time_seed, randn
    use minimisation, only: minimise
    use io, only: output, input
    use rp_emulator

    implicit none

    real(dp), dimension(n_x,tstep) :: truth
    real(dp), dimension(n_x,tstep) :: guess_traj
    real(dp), dimension(n_x,tstep/freq) :: obs, innov
    real(dp), dimension(n_x) :: initial, initial_del, final_del
    real(dp) :: time(tstep)
    real(dp) :: diagn(1,out_iter*max_iterations) = -1.0_dp, diagn_inner(1,max_iterations)
    integer :: i, j, num_iters = 1, num_iters_inner
    logical :: files_exist(3)

    !==========================================================================
    ! Setup
    !==========================================================================

    RPE_DEFAULT_SBITS = sbits

    ! Check whether observation frequency divides into total number of timsteps
    if (mod(tstep, freq) /= 0) then
        stop 'Number of timesteps must be a multiple of the observation frequency'
    end if

    ! Seed random number generator
    call time_seed

    ! Generate time axis
    time = (/ (real(i)*h, i = 0, tstep-1) /)

    !==========================================================================
    ! Generate truth and observations and initial guess
    !==========================================================================
    
    ! If all of the files already exist, use them for the truth, observations and first guess
    inquire(file="truth.txt", exist=files_exist(1))
    inquire(file="obs.txt", exist=files_exist(2))
    inquire(file="first_guess.txt", exist=files_exist(3))
    if (all(files_exist)) then
        call input(truth, "truth.txt")
        call input(obs, "obs.txt")
        call input(guess_traj, "first_guess.txt")
    ! Else create them all from scratch
    else
        ! Spin up truth
        truth(:,1) = (/ (randn(0.0_dp, 5.0_dp), i = 1, n_x) /)
        do i = 1, spin_up
            truth(:,:1) = run_model(1, truth(:,1))
        end do

        ! Run truth
        truth = run_model(tstep, truth(:,1))

        ! Calculate observations
        do i = 1, tstep, freq
            do j = 1, n_x
                obs(j, 1+i/freq) = randn(truth(j,i), sqrt(obs_var))
            end do
        end do

        ! Output truth and observations
        call output(time, truth, "truth.txt")
        call output(time, obs, "obs.txt", freq)

        ! Set initial best guess
        initial = (/ (truth(i,1) + randn(0.0_dp, init_err), i = 1, n_x ) /)

        ! Output first guess
        guess_traj = run_model(tstep, initial)
        call output(time, guess_traj, "first_guess.txt")
    end if

    ! Store the time index of the last observation
    do i = 1, tstep, freq
        last = i
    end do

    !==========================================================================
    ! Perform minimisation
    !==========================================================================

    ! Perform outer loop
    do j = 1, out_iter
        ! Generate nonlinear forecast
        guess_traj = run_model(tstep, initial)

        ! Calculate innovations
        innov = obs - guess_traj(:,1:tstep:freq)

        ! Initial increment for inner loop
        initial_del = 0.0_dp

        ! Perform inner loop
        call minimise(initial_del, innov, guess_traj, final_del, diagn_inner)

        ! If loop ended early
        if (any(diagn_inner(1,:) < 0.0_dp)) then
            num_iters_inner = minloc(diagn_inner(1,:), dim=1)-1
        else
            num_iters_inner = max_iterations
        end if

        ! Store diagnostics
        diagn(1,num_iters:num_iters+num_iters_inner-1) = diagn_inner(1,:)
        num_iters = num_iters + num_iters_inner

        ! Update best guess
        initial = initial + final_del
    end do

    !==========================================================================
    ! Finish
    !==========================================================================

    ! Output final best guess
    guess_traj = run_model(tstep, initial)
    call output(time, guess_traj, "final_guess.txt")

    ! Output diagnostics
    call output((/ (real(i,dp), i = 1, size(diagn(1,:)) ) /), diagn, "diagnostics.txt")

    ! write (*,'(A5,I5,A11)') 'Took ', num_iters-1, ' iterations'
    ! write (*,'(A11,F9.2)') 'Final cost ', diagn(1,num_iters-1)
    print *, diagn(1,num_iters-1), num_iters
end program lorenz96_4dvar
