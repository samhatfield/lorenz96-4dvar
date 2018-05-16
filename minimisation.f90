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
    !> @param[in] initial the initial guess
    !> @param[in] obs the observations
    !> @param[inout] best_guess the computed minimum
    !> @param[inout] diagn array containing diagnostics of minimisation
    !> performance
    subroutine minimise(initial, obs, best_guess, diagn)
        real(ap), intent(in) :: initial(n_x), obs(n_x,tstep/freq)
        real(ap), intent(inout) :: best_guess(n_x,tstep), diagn(1,max_iterations)
        real(ap) :: cost, grad(n_x)
        integer :: iters = 1, i

        ! Dummy variables for conjugate gradient algorithm
        real(ap) :: d(n_x), grad_old(n_x), w(n_x)
    
        ! Flags to control conjugate gradient algorithm behaviour
        integer :: printflags(2) = (/ -1, 2 /), flag = 0, rest = 0, method = 3
    
        ! Variables required by conjugate gradient subroutine call, but that aren't actually used
        real(ap) :: eps = 1.0d-5
        logical :: finish = .false.

        diagn = -1.0_dp

        do
            ! Compute cost of current best guess
            if (flag == 0 .or. flag == 1) then
                best_guess = run_model(tstep, initial)
                cost = calc_cost(tstep, best_guess, obs)
                diagn(1,iters) = cost
        
                ! Compute gradient of cost function
                grad = calc_cost_grad(tstep, best_guess, obs)
            end if
    
            ! Use gradient descent algorithm to step towards minimum of cost function
            call cgfam(n_x, initial, cost, grad, d, grad_old, printflags, eps, w, flag, rest, method, finish)
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
