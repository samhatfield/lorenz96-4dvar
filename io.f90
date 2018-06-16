!> @author
!> Sam Hatfield, AOPP, University of Oxford
!> @brief
!> Contains functions for I/O.
module io
    use params

    implicit none

    private
    public output, input

contains
    !> @brief
    !> Outputs the given 2D array with the given time axis as the first column.
    !> @param[in] time_axis an array of time values corresponding to each row
    !> of the output array
    !> @param[in] output_array the array containing the data to output
    !> @param[in] filename the output file name
    !> @param[in] stride_in an optional argument allowing only every nth value
    !> to be output
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
        do i = 1, size(output_array, 2)
            write (1,*) time_axis(1+(i-1)*stride), output_array(:,i)
        end do
        close(1)
    end subroutine output
    
    !> @brief
    !> Reads in the data from the given file name into the given 2D array.
    !> @param[in] input_array the array containing the data to output
    !> @param[in] filename the input file name
    subroutine input(input_array, filename)
        real(dp), intent(inout) :: input_array(:,:)
        character(len=*), intent(in) :: filename
        ! For storing the first column, which is not used (it's the timestep)
        real(dp) :: dummy
        integer :: i
        
        open(1, file=filename)
        do i = 1, size(input_array, 2)
            read(1,*) dummy, input_array(:,i)
        end do
        close(1)
    end subroutine
end module io
