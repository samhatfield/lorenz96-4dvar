!> @author
!> Sam Hatfield, AOPP, University of Oxford
!> @brief
!> Contains functions for integrating the Lorenz '96 system, along with its
!> tangent linear and adjoint models.
module lorenz96
    use params
    use rp_emulator

    implicit none

    private
    public run_model, run_tangent_linear, run_adjoint

    ! Model parameters
    real(dp) :: F = 8.0_dp

contains
    !> @brief
    !> The full nonlinear ODE.
    !> @param[in] x the state vector at which to evaluate the ODE
    !> @return dRdT the evaluated ODE
    function dRdT(x)
        type(rpe_var), intent(in) :: x(n_x)
        type(rpe_var) :: dRdT(n_x)

        dRdT = (cshift(x, 1) - cshift(x, -2))*cshift(x, -1) - x + rpe_literal(F)
    end function dRdT

    !> @brief
    !> The Jacobian of the ODE multiplied by the given perturbation.
    !> @param[in] x the state vector at which to evaluate the Jacobian
    !> @param[in] dx the perturbation vector
    !> @return jacob the evaluated Jacobian multiplied by the given
    !> perturbation
    function jacob(x, dx)
        type(rpe_var), intent(in) :: x(n_x), dx(n_x)
        type(rpe_var) :: jacob(n_x)

        jacob = (cshift(x, 1) - cshift(x, -2))*cshift(dx, -1) &
            & + (cshift(dx, 1) - cshift(dx, -2))*cshift(x, -1) - dx
    end function jacob

    !> @brief
    !> The adjoint of the Jacobian of the ODE multiplied by the given
    !> perturbation.
    !> @param[in] x the state vector at which to evaluate the adjoint of
    !> the Jacobian
    !> @param[in] dx_a the perturbation vector
    !> @return jacob the evaluated adjoint of the Jacobian multiplied by the
    !> given perturbation
    function jacob_a(x, dx_a)
        type(rpe_var), intent(in) :: x(n_x), dx_a(n_x)
        type(rpe_var) :: jacob_a(n_x)

        jacob_a = (cshift(x, 2) - cshift(x, -1))*cshift(dx_a, 1) &
            & - cshift(x, 1)*cshift(dx_a, 2) + cshift(x, -2)*cshift(dx_a, -1) - dx_a
    end function jacob_a

    !> @brief
    !> Run full nonlinear model using modified Euler scheme.
    !> @param[in] tstep the length of the integration
    !> @param[in] in the initial condition
    !> @return out the computed trajectory
    function run_model(tstep, in_array_dp) result(out_array_dp)
        integer, intent(in) :: tstep
        real(dp), intent(in) :: in_array_dp(n_x)
        real(dp), dimension(n_x,tstep) :: out_array_dp
        type(rpe_var), dimension(n_x,tstep) :: out_array
        type(rpe_var), dimension(n_x) :: k1, k2, k3, k4
        type(rpe_var) :: half, two, six, h_
        integer :: i

        half = 0.5_dp; two = 2.0_dp; six = 6.0_dp
        h_ = h

        ! First step
        out_array(:,1) = in_array_dp

        do i = 1, tstep-1
            k1 = h_*dRdT(out_array(:,i))
            k2 = h_*dRdT(out_array(:,i) + half*k1)
            k3 = h_*dRdT(out_array(:,i) + half*k2)
            k4 = h_*dRdT(out_array(:,i) + k3)

            out_array(:,i+1) = out_array(:,i) + (k1 + two*k2 + two*k3 + k4)/six
        end do

        out_array_dp = out_array%val
    end function run_model

    !> @brief
    !> Run tangent linear model.
    !> @param[in] tstep the length of the integration
    !> @param[in] in the trajectory around which to perform the integration
    !> @param[in] din the initial condition
    !> @return dout the computed trajectory
    function run_tangent_linear(tstep, in_array_dp, din_dp) result(dout_dp)
        integer, intent(in) :: tstep
        real(dp), intent(in) :: in_array_dp(n_x,tstep), din_dp(n_x)
        real(dp) :: dout_dp(n_x,tstep)
        type(rpe_var) :: in_array(n_x,tstep), dout(n_x,tstep)
        type(rpe_var), dimension(n_x) :: k1, k2, k3, dk1, dk2, dk3, dk4
        type(rpe_var) :: half, two, six, h_
        integer :: i

        half = 0.5_dp; two = 2.0_dp; six = 6.0_dp
        h_ = h

        in_array = in_array_dp

        ! First step
        dout(:,1) = din_dp

        do i = 1, tstep-1
            k1 = h_*dRdT(in_array(:,i))
            k2 = h_*dRdT(in_array(:,i) + half*k1)
            k3 = h_*dRdT(in_array(:,i) + half*k2)
            dk1 = h_*jacob(in_array(:,i), dout(:,i))
            dk2 = h_*jacob(in_array(:,i) + half*k1, dout(:,i) + half*dk1)
            dk3 = h_*jacob(in_array(:,i) + half*k2, dout(:,i) + half*dk2)
            dk4 = h_*jacob(in_array(:,i) + k3, dout(:,i) + dk3)

            dout(:,i+1) = dout(:,i) + (dk1 + two*dk2 + two*dk3 + dk4)/six
        end do

        dout_dp = dout%val
    end function run_tangent_linear

    !> @brief
    !> Run adjoint model.
    !> @param[in] tstep the length of the integration
    !> @param[in] in the trajectory around which to perform the integration
    !> @param[in] din_a the initial condition
    !> @return dout_a the computed trajectory
    function run_adjoint(tstep, in_array_dp, din_a_dp) result(dout_a_dp)
        integer, intent(in) :: tstep
        real(dp), intent(in) :: in_array_dp(:,:), din_a_dp(n_x)
        real(dp) :: dout_a_dp(n_x)
        type(rpe_var) :: in_array(size(in_array_dp,1),size(in_array_dp,2)), dout_a(n_x)
        type(rpe_var), dimension(n_x) :: k1, k2, k3, dk1_a, dk2_a, dk3_a, dk4_a
        type(rpe_var) :: half, two, six, h_
        integer :: i

        half = 0.5_dp; two = 2.0_dp; six = 6.0_dp
        h_ = h

        in_array = in_array_dp

        ! First step
        dout_a = din_a_dp

        do i = tstep-1, 1, -1
            k1 = h_*dRdT(in_array(:,i))
            k2 = h_*dRdT(in_array(:,i) + half*k1)
            k3 = h_*dRdT(in_array(:,i) + half*k2)
            dk1_a = h_*jacob_a(in_array(:,i) + k3, dout_a)
            dk2_a = h_*jacob_a(in_array(:,i) + half*k2, dout_a + half*dk1_a)
            dk3_a = h_*jacob_a(in_array(:,i) + half*k1, dout_a + half*dk2_a)
            dk4_a = h_*jacob_a(in_array(:,i), dout_a + dk3_a)

            dout_a = dout_a + (dk1_a + two*dk2_a + two*dk3_a + dk4_a)/six
        end do

        dout_a_dp = dout_a
    end function run_adjoint
end module lorenz96
