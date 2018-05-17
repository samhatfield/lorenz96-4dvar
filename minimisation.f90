!> @author
!> Sam Hatfield, AOPP, University of Oxford
!> @brief
!> Contains minimisation subroutine.
module minimisation
    use params
    use assim, only: calc_cost_grad, calc_cost
    use lorenz96, only: run_model
    use io, only: output

    implicit none

    private
    public minimise

contains
    !> @brief
    !> Find minimum of cost function given an initial guess using the conjugate
    !> gradient algorithm.
    !> @param[in] initial_del the initial guess
    !> @param[in] innov the innovations
    !> @param[in] guess_traj the trajectory around which the minimisation is
    !> performed
    !> @param[inout] del the computed minimum
    !> @param[inout] diagn array containing diagnostics of minimisation
    !> performance
    subroutine minimise(initial_del, innov, guess_traj, del, diagn)
        real(ap), intent(in) :: initial_del(n_x), innov(n_x,tstep/freq), guess_traj(n_x,tstep)
        real(ap), intent(inout) :: del(n_x), diagn(1,max_iterations)
        real(ap) :: cost, grad(n_x)
        integer :: iters, i

        ! Dummy variables for conjugate gradient algorithm
        real(ap) :: d(n_x), grad_old(n_x), w(n_x)
    
        ! Flags to control conjugate gradient algorithm behaviour
        integer :: printflags(2) = (/ -1, 2 /), flag, rest, method = 3
    
        ! Variables required by conjugate gradient subroutine call, but that aren't actually used
        real(ap) :: eps = 1.0d-5
        logical :: finish

        ! Initialise minimisation parameters
        flag = 0
        rest = 0
        iters = 1
        finish = .false.
        del = initial_del
        diagn = -1.0_ap

        do
            ! Compute cost of current best guess
            if (flag == 0 .or. flag == 1) then
                cost = calc_cost(tstep, guess_traj, del, innov)
                diagn(1,iters) = cost
        
                ! Compute gradient of cost function
                grad = calc_cost_grad(tstep, guess_traj, del, innov)
            end if
    
            ! Use gradient descent algorithm to step towards minimum of cost function
            call cgfam(n_x, del, cost, grad, d, grad_old, printflags, eps, w, flag, rest, method, finish)
            if (flag <= 0 .or. iters >= max_iterations) exit
            if (flag == 1) iters = iters + 1
            if (flag == 2) then
                if (maxval(abs(grad)) <= 0.2_ap) then
                    finish = .true.
                end if
            end if
        end do
    end subroutine minimise
end module minimisation
