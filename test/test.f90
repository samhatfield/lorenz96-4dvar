program test
    use params
    use lorenz96, only: run_model, run_tangent_linear, run_adjoint
    use assim, only: calc_cost, calc_cost_grad
    use utils, only: randn
    use rp_emulator

    implicit none
    
    RPE_DEFAULT_SBITS = sbits

    call time_seed

    call test_tl
    call test_adj
    call test_grad
contains
    subroutine test_tl
        integer, parameter :: nsteps = 10
        real(dp), dimension(n_x) :: pert_orig, orig_state, pert, pert_state, diff
        real(dp), dimension(n_x,nsteps) :: orig_traj, pert_traj, pert_traj_tl
        real(dp) :: rel_err
        real(dp) :: gamma = 100_dp
        integer :: i

        print *, '============================================================'
        print *, 'Tangent linear test'
        print *, '============================================================'
        print *, 'The relative error should tend towards zero as the magnitude'
        print *, 'of the perturbation decreases.'
        print *, '============================================================'

        pert_orig = (/ ( randn(0.0_dp, 1.0_dp), i = 1, n_x ) /)

        write (*, '(A10, A16)') 'Magnitude', 'Relative error'
        do i = 1, 9
            gamma = gamma * 0.1_dp

            orig_state = (/ (1.0_dp, i = 1, n_x) /)
            pert = gamma * pert_orig
            pert_state = orig_state + pert

            orig_traj = run_model(nsteps, orig_state)
            pert_traj = run_model(nsteps, pert_state)
            pert_traj_tl = run_tangent_linear(nsteps, orig_traj, pert)

            diff = pert_traj(:,nsteps) - orig_traj(:,nsteps)
            rel_err = 100.0_dp * norm2(diff - pert_traj_tl(:,nsteps)) &
                & / norm2(pert_traj_tl(:,nsteps))

            write (*,'(E10.1, F16.9)') gamma, rel_err
        end do
        print *, '============================================================'
    end subroutine test_tl

    subroutine test_adj
        integer, parameter :: nsteps = 10
        real(dp), dimension(n_x) :: pert_1, pert_2, adj_traj
        real(dp), dimension(n_x,nsteps) :: traj, fwd_traj
        real(dp) :: forward_product, adjoint_product
        integer :: i

        print *, ''
        print *, '============================================================'
        print *, 'Adjoint model test'
        print *, '============================================================'
        print *, 'The inner product of the integration of the first'
        print *, 'perturbation with the second perturbation should equal the'
        print *, 'inner product of the first perturbation with the second'
        print *, 'perturbation put through the adjoint model.'
        print *, '<Mdx, dy> = <dx, M^Tdyz>'
        print *, '============================================================'

        ! Generate random unit vector
        call random_number(pert_1)
        call random_number(pert_2)

        ! Get nonlinear trajectory
        traj = run_model(nsteps, (/ (1.0_dp, i = 1, n_x) /))
        
        ! Compute tangent linear evolution of first perturbation
        fwd_traj = run_tangent_linear(nsteps, traj, pert_1)
        
        ! Compute adjoint evolution of second perturbation
        adj_traj = run_adjoint(nsteps, traj, pert_2)

        ! Compute LHS and RHS of above identity
        forward_product = dot_product(fwd_traj(:,nsteps), pert_2)
        adjoint_product = dot_product(pert_1, adj_traj)

        write (*, '(3A20)') 'Forward product', 'Adjoint product', 'Difference'
        write (*,'(3F20.16)') forward_product, adjoint_product, forward_product - adjoint_product
        print *, '============================================================'
    end subroutine test_adj

    subroutine test_grad
        integer, parameter :: n_samples = 14
        real(dp), dimension(n_x) :: initial, pert, grad1, grad2, rand_unit
        real(dp) :: grad1_norm(n_x), alpha = 1.0_dp, obs(n_x,tstep/freq), innov(n_x,tstep/freq)
        real(dp) :: cost1, cost2
        real(dp) :: traj(n_x,tstep), del(n_x)
        integer :: i, j

        print *, ''
        print *, '============================================================'
        print *, 'Gradient test'
        print *, '============================================================'
        print *, 'The relative error should tend towards one as the magnitude'
        print *, 'of the perturbation decreases.'
        print *, '============================================================'

        initial = (/ ( randn(0.0_dp, 1.0_dp), i = 1, n_x ) /)
        pert = (/ ( randn(0.0_dp, 1.0_dp), i = 1, n_x ) /)

        ! Unperturbed input to cost function
        del = 0.1_dp

        ! Compute trajectory to linearise about
        traj = run_model(tstep, initial)

        ! Calculate observations
        do i = 1, tstep, freq
            last = i
            do j = 1, n_x
                obs(j, 1+i/freq) = randn(traj(j,i), sqrt(obs_var))
            end do
        end do

        ! Compute innovations
        innov = obs - traj(:,1:tstep:freq)

        ! Compute cost function and gradient at unperturbed input value
        cost1 = calc_cost(tstep, traj, del, innov)
        grad1 = calc_cost_grad(tstep, traj, del, innov)

        do i = 1, n_samples
            ! Perturb input to cost function
            initial = del + alpha*pert

            ! Compute cost at perturbed input
            cost2 = calc_cost(tstep, traj, initial, innov)

            write (*,'(E10.1, F18.11)') alpha, (cost2 - cost1)/(dot_product(initial, grad1))

            alpha = alpha * 0.1_dp
        end do
        print *, '============================================================'
    end subroutine test_grad

    subroutine output(time_axis, output_array, filename, stride_in)
        real(dp), intent(in) :: time_axis(:), output_array(:,:)
        character(len=*), intent(in) :: filename
        integer, optional :: stride_in
        integer :: i, stride

        if (present(stride_in)) then
            stride = stride_in
        else
            stride = 1
        end if

        open(1, file=filename)
        do i = 1, size(output_array, 1)
            write (1,*) time_axis(i*stride), output_array(i,:)
        end do
        close(1)
    end subroutine output

    subroutine time_seed()
      integer :: i, n, clock
      integer, allocatable :: seed(:)

      call random_seed(size = n)
      allocate(seed(n))

      call system_clock(count=clock)

      seed = clock + 37 * (/ (i - 1, i = 1, n) /)
      call random_seed(put = seed)

      deallocate(seed)
    end subroutine
end program test
