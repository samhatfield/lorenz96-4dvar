!> @author
!> Sam Hatfield, AOPP, University of Oxford
!> @brief
!> Contains experiment parameters.
module params
    use, intrinsic :: iso_fortran_env

    implicit none

    integer, parameter :: dp = real64
    integer, parameter :: sp = real32

    ! Write values of X variables in output file?
    logical, parameter :: output_x = .true.

    ! State vector properties
    integer, parameter :: n_x = 40

    real(dp), parameter :: h = 0.01_dp
    real(dp), parameter :: wind_len = 1.2_dp
    integer, parameter :: tstep = wind_len/h
    integer, parameter :: freq = 10
    integer, parameter :: n_obs = (tstep - 1) / freq + 1
    real(dp), parameter :: obs_var = 2.0_dp
    integer, parameter :: max_iterations = 3000
    integer :: last
end module params
