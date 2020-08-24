module numericalMethods
    implicit none

contains

    subroutine matrixPLU(a, p, l, u)
        real, dimension(:,:) :: a, p, l, u
        real :: max
        integer :: i = 0, j = 0, k = 0, n = 0, index = 1
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
            end do
        end do
        do k = 1, n-1
            max = u(k,k)
            do i = k, n-1
                ! Find maximum element in column below current diagonal entry
                if (abs(u(i+1,k)) >= max) then
                    ! Store index
                    index = i+1 
                    ! Update maximum variable
                    max = abs(u(i+1,k))
                end if
            end do
            ! Swap i and k rows in U
            call swapper(u, index, k, k, n, 0)
            ! Swap i and k rows till k-1 column in L
            call swapper(l, index, k, 1, k-1, 0)
            ! Swap i and k rows in P
            call swapper(p, index, k, 1, n, 0)
            do j = k+1, n 
                ! Divide (j,k) entry by (k,k) entry
                l(j,k) = u(j,k)/u(k,k)
                ! Subtract normalised element multiplied by row after kth column from original row
                u(j,k:n) = u(j,k:n) - l(j,k)*u(k,k:n)
            end do
        end do
    end subroutine matrixPLU

    subroutine swapper(matrix, i, j, start, end, axis)
        real, dimension(:,:) :: matrix
        integer :: i, j, start, end, axis

        if (axis == 0) then
            matrix(i,start:end) = matrix(i,start:end) + matrix(j,start:end)
            matrix(j,start:end) = matrix(i,start:end) - matrix(j,start:end)
            matrix(i,start:end) = matrix(i,start:end) - matrix(j,start:end)
        else
            matrix(start:end,i) = matrix(start:end,i) + matrix(start:end,j)
            matrix(start:end,j) = matrix(start:end,i) - matrix(start:end,j)
            matrix(start:end,i) = matrix(start:end,i) - matrix(start:end,j)
        end if

    end subroutine swapper

    subroutine bisection(f, a, b, tol)
        real, external :: f
        real :: a, b, c, tol
        do while (b - a > tol)
            c = (a + b)/2
            if (f(c) == 0) then
                b = c
                a = b
            else if (f(c)*f(a) < 0) then
                b = c
            else
                a = c
            end if
        end do 
    end subroutine bisection

    subroutine gramSchmidtOrthoNorm(u)
        real, dimension(:,:) :: u
        real, dimension(size(u, 1), size(u, 2)) :: v
        integer :: i, j, n
        n = size(u)
        do i = 1, n
            v(:,i) = u(:,i)
            do j = 1, n
                v(:,i) = v(:,i) - dot_product(v(:,j), v(:,i))*v(:,j)
            end do
            v(:,i) = v(:,i)/norm2(v(:,i))   
            u(:,i) = v(:,i)
        end do
    end subroutine gramSchmidtOrthoNorm

end module numericalMethods

pure function func(x)
    real, intent(in) :: x ! input
    func = x**3 + 1
end function func

program testCases
    use numericalMethods

    real, dimension(4,4) :: a, p, l, u, vec_list = 0
    integer :: i, j
    real :: x1 = -3, x2 = 3

    interface
        function func(i) result(j)
            real, intent(in) :: i
            real :: j
        end function func
    end interface

    ! Matrix PLU factorisation
    print *, "Matrix PLU factorisation:"
    print *, "A = "
    do i = 1,4
        do j = 1,4
            a(i,j) = i + j - 4
        end do
        a(4,1) = 0
        print *, a(i,:)
    end do
    call matrixPLU(a, p, l, u)
    print *, "P ="
    do i = 1,4
        print *, p(i,:)
    end do
    print *, "L ="
    do i = 1,4
        print *, l(i,:)
    end do
    print *, "U ="
    do i = 1,4
        print *, u(i,:)
    end do

    print *, "Bisection method:"
    call bisection(func, x1, x2, 1e-5)
    print *, "x1 =", x1, "x2 =", x2

    print *, "Gram-Schmidt Orthonormalisation:"
    do j = 1, 4
        call random_number(vec_list(:,j))
    end do
    print *, "Vector list:"
    do i = 1, 4
        print *, vec_list(i,:)
    end do
    call gramSchmidtOrthoNorm(vec_list)
    print *, "Orthonormalised:"
    do i = 1, 4
        print *, vec_list(i,:), "WRONG?"
    end do

end program testCases
