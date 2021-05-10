! program no. 6: linear strength source
! -------------------------------------
! this program finds the pressure distribution on an arbitrary airfoil
! by representing the surface as a finite number of source panels with
! linear strength (alpha=0, neumann b.c., program by steven yon, 1989).
module numericalMethods
    implicit none

contains

    subroutine matrx(a,n,g)
        real, dimension(n,n) :: a, temp
        real, dimension(:) :: g
        integer :: n
        integer :: i = 0, j = 0, k = 0, l = 0, p = 0, p2 = 0

        ! Initialisation
        do i = 1, n-1
            g(i) = 0
        end do

        ! Convert coefficient matrix to upper triangular form
        do i = 1, n-1
            5 if (abs(a(i,i)) .lt. 1e-7) goto 9

            p = a(i,i)
            do j = i, n
                a(i,j) = a(i,j)/p
            end do

            do k = i + 1, n - 1
                p2 = a(k,i)
                do l = i, n
                    a(k,l) = a(k,l) - p2*a(i,l)
                end do
            end do
        end do

        ! Back substitute triangularised matrix to get values of solution vector
        do i = n - 1, 1, -1
            g(i) = a(i,n)
            do j = 1, n-1
                a(i,i) = 0
                g(i) = g(i) - a(i,j)*g(j)
            end do
        end do

        return

        ! order matrix so that diagonal coefficients are non-zero and stop if matrix is singular
        9 if (i .ne. n - 1) then
            do j = 1,n
                temp(i,j) = a(i,j)
                a(i,j) = a(i+1,j)
                a(i+1,j) = temp(i,j)
            end do
            goto 5
        else
            goto 10
        endif

        10 write(6,*) "No solution"
        stop
    end

    subroutine linearsource(airfoil, a)
        character(32), intent(in) :: airfoil
        real ep(400,2),ept(400,2),pt1(400,2),pt2(400,2)
        real co(400,2),b(400,400),g(400),v(400)
        real a(400, 400)
        real th(400)
        integer :: i = 0, j = 0, m = 0, n = 0
        real :: xt, zt, x2t, z2t, x, z, x2, z2, r1, r2, th1, th2, al, xx
        real :: u1l, u2l, w1l, w2l, u1, u2, w1, w2, dz, dx, holda, holdb, vel, va, cp, ans1

        ! open(8,file='cpd2.dat',status='new')
        open(9,file=airfoil,status='old')
        write(6,*) 'enter number of panels'
        read(5,*) m

        n=m+1
        al=0

        ! read in the panel end points
        do i=1,m+1
            read(9,*) ept(i,1), ept(i,2)
        end do

        ! convert paneling to clockwise
        do i=1,n
            ep(i,1)=ept(n-i+1,1)
            ep(i,2)=ept(n-i+1,2)
        end do

        ! establish coordinates of panel end points
        do i=1,m
            pt1(i,1)=ep(i,1)
            pt2(i,1)=ep(i+1,1)
            pt1(i,2)=ep(i,2)
            pt2(i,2)=ep(i+1,2)
        end do

        ! find panel angles th(j)
        do i=1,m
            dz    =pt2(i,2)-pt1(i,2)
            dx    =pt2(i,1)-pt1(i,1)
            th(i) =atan2(dz,dx)
        end do

        th(m+1)=0

        ! establish collocation points
        do i=1,m
            co(i,1)=(pt2(i,1)-pt1(i,1))/2+pt1(i,1)
            co(i,2)=(pt2(i,2)-pt1(i,2))/2+pt1(i,2)
        end do

        write(6,*) 'enter x coord. of wake point'

        read(5,*) xx

        co(m+1,1)=xx
        co(m+1,2)=0

        ! establish influence coefficients
        do i=1,m+1
            do j=1,m
            ! convert collocation point to local panel coords.
            xt  = co(i,1)  - pt1(j,1)
            zt  = co(i,2)  - pt1(j,2)
            x2t = pt2(j,1) - pt1(j,1)
            z2t = pt2(j,2) - pt1(j,2)
            x   =  xt * cos(th(j)) + zt * sin(th(j))
            z   = -xt * sin(th(j)) + zt * cos(th(j))
            x2  =  x2t* cos(th(j)) + z2t* sin(th(j))
            z2  = 0
            ! find r1, r2, th1, th2
            r1  = sqrt(x**2+z**2)
            r2  = sqrt((x-x2)**2+z**2)
            th1 = atan2(z,x)
            th2 = atan2(z,x-x2)
            ! compute velocity components as functions of
            ! sigma1 and sigma2. these velocities are in
            ! the jth reference frame.
            if(i.eq.j) then
                u1l = 0.15916
                w1l = -0.5*(x-x2)/x2

                u2l = -0.15916
                w2l =  0.5*(x)/x2
            else
                u1l=  ((x2 - z * (th2-th1)) + (x - x2) * log(r2/r1))/(6.28319*x2)
                w1l=  -(z * log(r2/r1)      + (x - x2) * (th2-th1))/(6.28319*x2)

                u2l= -((x2 - z * (th2-th1)) + x * log(r2/r1))/(6.28319*x2)
                w2l=   (z * log(r2/r1)      + x * (th2-th1))/(6.28319*x2)
            end if
            ! transform the local velocities into the global
            ! reference frame.
            u1 =  u1l * cos(-th(j)) + w1l * sin(-th(j))
            u2 =  u2l * cos(-th(j)) + w2l * sin(-th(j))
            w1 = -u1l * sin(-th(j)) + w1l * cos(-th(j))
            w2 = -u2l * sin(-th(j)) + w2l * cos(-th(j))
            ! compute the coefficients of sigma in the
            ! influence matrix
            if(j.eq.1) then
                a(i,1) = -u1 * sin(th(i)) + w1 * cos(th(i))
                holda  = -u2 * sin(th(i)) + w2 * cos(th(i))
                b(i,1) =  u1 * cos(th(i)) + w1 * sin(th(i))
                holdb  =  u2 * cos(th(i)) + w2 * sin(th(i))
            else if(j.eq.m) then
                a(i,m) = -u1 * sin(th(i)) + w1 * cos(th(i)) + holda
                a(i,n) = -u2 * sin(th(i)) + w2 * cos(th(i))
                b(i,m) =  u1 * cos(th(i)) + w1 * sin(th(i)) + holdb
                b(i,n) =  u2 * cos(th(i)) + w2 * sin(th(i))
            else
                a(i,j) = -u1 * sin(th(i)) + w1 * cos(th(i)) + holda
                holda  = -u2 * sin(th(i)) + w2 * cos(th(i))
                b(i,j) =  u1 * cos(th(i)) + w1 * sin(th(i)) + holdb
                holdb  =  u2 * cos(th(i)) + w2 * sin(th(i))
            end if
            end do
            a(i,n+1) = sin(th(i))
        end do
        n=m+2

        ! if(m.eq.10) then
        !     do i=1,11
        !         write(6,10) a(i,1),a(i,2),a(i,3),a(i,4),a(i,5),a(i,6),a(i,7),a(i,8),a(i,9),a(i,10),a(i,11)
        !     end do
        ! end if
        ! solve for the solution vector of source strengths
        call matrx(a,n,g)
        ! convert source strengths into tangential
        ! velocities along the airfoil surface and cp's
        ! on each of the panels.

        continue

        n=m+1
        do i=1,m
            vel=0
            do j=1,n
                vel=vel+b(i,j)*g(j)
            end do
            v(i)=vel+cos(al)*cos(th(i))+sin(al)*sin(th(i))
        end do

        write(6,*) ' '
        write(6,*) 'smooth the velocity distribution?'
        write(6,*) '1 = yes'
        write(6,*) '2 = no'
        read(5,*) ans1

        do i=2,m
            if(ans1.eq.1) then
                va=(v(i)+v(i-1))/2
                cp=1-va**2
                write(8,*) pt1(i,1), ' ,',cp
            else
                cp=1-v(i)**2
                write(8,*) co(i,1),' ,',cp
            end if
        end do

        write(6,*) ' '
        write(6,*) 'lift coefficient = 0'
        stop
    end 

end module numericalMethods

program airfoil
    use numericalMethods
    real a(400, 400)
    call linearsource("../coordinates/cosinepanels.dat", a)
    
    do i = 1, m+1
        print *, a(i,:), "WRONG?"
    end do
end program airfoil
! format(,f6.2,1x,f5.2,1x,f5.2,1x,f5.2,1x,f5.2,1x,f5.2,1x,f5.2,1x,f5.2,1x,f5.2,1x,f5.2,1x,f5.2)