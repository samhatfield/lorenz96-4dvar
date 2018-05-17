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
    use io, only: output

    implicit none

    real(ap), dimension(n_x,tstep) :: truth = 0.0_ap
    real(ap), dimension(n_x,tstep) :: guess_traj
    real(ap) :: obs(n_x,tstep/freq)
    real(ap), dimension(n_x) :: initial, best_guess
    real(ap) :: time(tstep)
    real(ap) :: diagn(1,max_iterations)
    integer :: i, j

    !==========================================================================
    ! Setup
    !==========================================================================

    ! Check whether observation frequency divides into total number of timsteps
    if (mod(tstep, freq) /= 0) then
        stop 'Number of timesteps must be a multiple of the observation frequency'
    end if

    ! Seed random number generator
    call time_seed

    ! Generate time axis
    time = (/ (real(i)*h, i = 0, tstep-1) /)

    !==========================================================================
    ! Generate truth and observations
    !==========================================================================
    
    ! Spin up truth
    truth(:,1) = (/ (randn(0.0_ap, 5.0_ap), i = 1, n_x) /)
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

    !==========================================================================
    ! Set initial guess
    !==========================================================================

    ! Set initial best guess
    initial = (/ (truth(i,1) + randn(0.0_ap, init_err), i = 1, n_x ) /)

    ! Output first guess
    guess_traj = run_model(tstep, initial)
    call output(time, guess_traj, "first_guess.txt")
    
    !==========================================================================
    ! Perform minimisation
    !==========================================================================

    call minimise(initial, obs, best_guess, diagn)

    !==========================================================================
    ! Finish
    !==========================================================================

    ! Output final best guess
    guess_traj = run_model(tstep, best_guess)
    call output(time, guess_traj, "final_guess.txt")

    ! Output diagnostics
    call output((/ (real(i,ap), i = 1, max_iterations) /), diagn, "diagnostics.txt")

    write (*,'(A5,I5,A11)') 'Took ', size(pack(diagn(1,:), diagn(1,:) >= 0.0_ap)), ' iterations'
    write (*,'(A11,F9.2)') 'Final cost ', diagn(1,minloc(diagn(1,:))-1)
end program lorenz96_4dvar
