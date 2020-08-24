module numericalMethods
    implicit none

contains

    subroutine matrixPLU(a, p, l, u)
        real, dimension(:,:) :: a, p, l, u
        real :: max
        integer :: i = 0, j = 0, k = 0, n = 0, index = 0
        n = size(a, 1)
        
        ! Initializing matrices
        u(:,:) = a(:,:)
        do i = 1, n
            do j = 1, n
                if (i == j) then
                    l(i,j) = 1
                    p(i,j) = 1
                else
                    l(i,j) = 0
                    p(i,j) = 0
                end if
                ! u(i,j) = a(i,j)
            end do
        end do

        do k = 1, n-1
            do i = k, n-1
                ! Find maximum element in row
                max = u(i,k)
                if (u(i+1,k) > max) then
                    index = i+1 
                    max = u(i+1,k)
                end if
            end do
            ! Swap i and k rows in U
            call swapper(u, i, k, k, n, 0)
            ! Swap i and k rows till k-1 column in L
            call swapper(l, i, k, 1, k-1, 0)
            ! Swap i and k rows in P
            call swapper(p, i, k, 1, n, 0)
            do j = k+1, n 
                l(j,k) = u(j,k)/u(k,k)
                u(j,k:n) = u(j,k:n) - l(j,k)*u(k,k:n)
            end do
        end do
    end subroutine

    subroutine swapper(matrix, in1, in2, start, end, axis)
        real, dimension(:,:) :: matrix
        integer :: in1, in2, start, end, axis

        if (axis == 0) then
            matrix(in1,start:end) = matrix(in1,start:end) + matrix(in2,start:end)
            matrix(in2,start:end) = matrix(in1,start:end) - matrix(in2,start:end)
            matrix(in1,start:end) = matrix(in1,start:end) - matrix(in2,start:end)
        else
            matrix(start:end,in1) = matrix(start:end,in1) + matrix(start:end,in2)
            matrix(start:end,in2) = matrix(start:end,in1) - matrix(start:end,in2)
            matrix(start:end,in1) = matrix(start:end,in1) - matrix(start:end,in2)
        end if

    end subroutine

end module numericalMethods