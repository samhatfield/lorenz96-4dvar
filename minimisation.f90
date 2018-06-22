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
        real(dp), intent(in) :: initial_del(n_x), innov(n_x,tstep/freq), guess_traj(n_x,tstep)
        real(dp), intent(inout) :: del(n_x), diagn(1,max_iterations)
        real(dp) :: cost, grad(n_x)
        integer :: i, iter
        common /runinf/iter

        ! Dummy variables for conjugate gradient algorithm
        real(dp) :: d(n_x), grad_old(n_x), w(n_x)
    
        ! Flags to control conjugate gradient algorithm behaviour
        integer :: printflags(2) = (/ -1, 2 /), flag, rest, method = 3
    
        ! Variables required by conjugate gradient subroutine call, but that aren't actually used
        real(dp) :: eps = 1.0d-5
        logical :: finish

        ! Initialise minimisation parameters
        flag = 0
        rest = 0
        finish = .false.
        del = initial_del
        diagn = -1.0_dp

        do
            ! Compute cost of current best guess
            if (flag == 0 .or. flag == 1) then
                cost = calc_cost(tstep, guess_traj, del, innov)

                if (flag == 0) diagn(1,1) = cost

                ! Compute gradient of cost function
                grad = calc_cost_grad(tstep, guess_traj, del, innov)
            end if
    
            ! Use gradient descent algorithm to step towards minimum of cost function
            call cgfam(n_x, del, cost, grad, d, grad_old, printflags, eps, w, flag, rest, method, finish)
            if (flag <= 0 .or. iter >= max_iterations) exit
            if (flag == 2) then
                ! Save cost in diagnostics array
                diagn(1,iter+1) = cost
                if (maxval(abs(grad)) <= 1.0_dp) then
                    finish = .true.
                end if
            end if
        end do
    end subroutine minimise
end module minimisation
