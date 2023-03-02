module testing
    ! Module with subroutines for testing quickly that everything works.
    ! This is for developers use and shuld not be imported in the final version of the main program.
contains

    subroutine testMatrix(var)
        ! Prints the values of a 3xN matrix in a readable way.
        ! Useful for quickly see r/v/F matrices.
        implicit none
        double precision, intent(in), dimension(:,:) :: var
        integer :: i
    
        do i=1,size(var,dim=2)
            print*,var(:,i)
        end do
    end subroutine testMatrix

    subroutine lowercase(s1)
        ! Transforms the input string to lowercase.
        character(*), intent(inout) :: s1
        character(len(s1)):: s2
        character          :: ch
        integer,parameter  :: DUC = ichar('A') - ichar('a')
        integer            :: i
        
        do i = 1,len(s1)
           ch = s1(i:i)
           print*,ch
           if (ch >= 'A'.AND.ch <= 'Z') ch = char(ichar(ch)-DUC)
           s2(i:i) = ch
        end do
        s1 = s2
    end subroutine lowercase

end module testing