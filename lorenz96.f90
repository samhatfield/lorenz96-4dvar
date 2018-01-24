!> @author
!> Sam Hatfield, AOPP, University of Oxford
!> @brief
!> Contains functions for integrating the Lorenz '96 system, along with its
!> tangent linear and adjoint models.
module lorenz96
    use params

    implicit none

    private
    public run_model, run_tangent_linear, run_adjoint

    ! Model parameters
    real(ap) :: F = 20.0_ap

contains
    !> @brief
    !> The full nonlinear ODE.
    !> @param[in] x the state vector at which to evaluate the ODE
    !> @return dRdT the evaluated ODE
    function dRdT(x)
        real(ap), intent(in) :: x(n_x)
        real(ap) :: dRdT(n_x)

        dRdT = (cshift(x, 1) - cshift(x, -2))*cshift(x, -1) - x + F
    end function dRdT

    !> @brief
    !> The Jacobian of the ODE multiplied by the given perturbation.
    !> @param[in] x the state vector at which to evaluate the Jacobian
    !> @param[in] dx the perturbation vector
    !> @return jacob the evaluated Jacobian multiplied by the given
    !> perturbation
    function jacob(x, dx)
        real(ap), intent(in) :: x(n_x), dx(n_x)
        real(ap) :: jacob(n_x)

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
        real(ap), intent(in) :: x(n_x), dx_a(n_x)
        real(ap) :: jacob_a(n_x)

        jacob_a = (cshift(x, 2) - cshift(x, -1))*cshift(dx_a, 1) &
            & - cshift(x, 1)*cshift(dx_a, 2) + cshift(x, -2)*cshift(dx_a, -1) - dx_a
    end function jacob_a

    !> @brief
    !> Run full nonlinear model using modified Euler scheme.
    !> @param[in] tstep the length of the integration
    !> @param[in] in the initial condition
    !> @return out the computed trajectory
    function run_model(tstep, in_array) result(out_array)
        integer, intent(in) :: tstep
        real(ap), intent(in) :: in_array(n_x)
        real(ap), dimension(n_x, tstep) :: out_array
        real(ap), dimension(n_x) :: k1, k2, k3, k4
        integer :: i

        ! First step
        out_array(:,1) = in_array

        do i = 1, tstep-1
            k1 = h*dRdT(out_array(:,i))
            k2 = h*dRdT(out_array(:,i) + 0.5_ap*k1)
            k3 = h*dRdT(out_array(:,i) + 0.5_ap*k2)
            k4 = h*dRdT(out_array(:,i) + k3)

            out_array(:,i+1) = out_array(:,i) + (k1 + 2.0_ap*k2 + 2.0_ap*k3 + k4)/6.0_ap
        end do
    end function run_model

    !> @brief
    !> Run tangent linear model.
    !> @param[in] tstep the length of the integration
    !> @param[in] in the trajectory around which to perform the integration
    !> @param[in] din the initial condition
    !> @return dout the computed trajectory
    function run_tangent_linear(tstep, in_array, din) result(dout)
        integer, intent(in) :: tstep
        real(ap), intent(in) :: in_array(n_x, tstep), din(n_x)
        real(ap) :: dout(n_x,tstep)
        real(ap), dimension(n_x) :: k1, k2, k3, dk1, dk2, dk3, dk4
        integer :: i

        ! First step
        dout(1,:) = din

        do i = 1, tstep-1
            k1 = h*dRdT(in_array(:,i))
            k2 = h*dRdT(in_array(:,i) + 0.5_ap*k1)
            k3 = h*dRdT(in_array(:,i) + 0.5_ap*k2)
            dk1 = h*jacob(in_array(:,i), dout(:,i))
            dk2 = h*jacob(in_array(:,i) + 0.5_ap*k1, dout(:,i) + 0.5_ap*dk1)
            dk3 = h*jacob(in_array(:,i) + 0.5_ap*k2, dout(:,i) + 0.5_ap*dk2)
            dk4 = h*jacob(in_array(:,i) + k3, dout(:,i) + dk3)

            dout(:,i+1) = dout(:,i) + (dk1 + 2.0_ap*dk2 + 2.0_ap*dk3 + dk4)/6.0_ap
        end do
    end function run_tangent_linear

    !> @brief
    !> Run adjoint model.
    !> @param[in] tstep the length of the integration
    !> @param[in] in the trajectory around which to perform the integration
    !> @param[in] din_a the initial condition
    !> @return dout_a the computed trajectory
    function run_adjoint(tstep, in_array, din_a) result(dout_a)
        integer, intent(in) :: tstep
        real(ap), intent(in) :: in_array(:,:), din_a(n_x)
        real(ap) :: dout_a(n_x)
        real(ap), dimension(n_x) :: k1, k2, k3, dk1_a, dk2_a, dk3_a, dk4_a
        integer :: i

        ! First step
        dout_a = din_a

        do i = tstep-1, 1, -1
            k1 = h*dRdT(in_array(:,i))
            k2 = h*dRdT(in_array(:,i) + 0.5_ap*k1)
            k3 = h*dRdT(in_array(:,i) + 0.5_ap*k2)
            dk1_a = h*jacob_a(in_array(:,i) + k3, dout_a)
            dk2_a = h*jacob_a(in_array(:,i) + 0.5_ap*k2, dout_a + 0.5_ap*dk1_a)
            dk3_a = h*jacob_a(in_array(:,i) + 0.5_ap*k1, dout_a + 0.5_ap*dk2_a)
            dk4_a = h*jacob_a(in_array(:,i), dout_a + dk3_a)

            dout_a = dout_a + (dk1_a + 2.0_ap*dk2_a + 2.0_ap*dk3_a + dk4_a)/6.0_ap
        end do
    end function run_adjoint
end module lorenz96
