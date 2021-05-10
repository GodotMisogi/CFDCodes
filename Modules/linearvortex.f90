! PROGRAM No. 7: LINEAR STRENGTH VORTEX
! -------------------------------------
! THIS PROGRAM FINDS THE PRESSURE DISTRIBUTION ON AN ARBITRARY AIRFOIL
! BY REPRESENTING THE SURFACE AS A FINITE NUMBER OF VORTEX PANELS WITH
! LINEAR STRENGTH (NEUMANN B.C., PROGRAM BY STEVEN YON, 1989).

real ep(400,2),ept(400,2),pt1(400,2),pt2(400,2)
real co(400,2),a(400,400),b(400,400),g(400)
real th(400),dl(400)

open(8,file='cplv.dat',status='new')
open(9,file='afoil2.dat',status='old')
write(6,*) 'enter number of panels'
read(5,*) m

n=m+1
write(6,*) 'enter angle of attack in degrees'
read(5,*) alpha

al=alpha/57.2958
! read in the panel end points
do i=1,m+1
    read(9,*) ept(i,1), ept(i,2)
end do

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
    dz=pt2(i,2)-pt1(i,2)
    dx=pt2(i,1)-pt1(i,1)
    th(i)=atan2(dz,dx)
end do
! establish collocation points
do i=1,m
    co(i,1)=(pt2(i,1)-pt1(i,1))/2+pt1(i,1)
    co(i,2)=(pt2(i,2)-pt1(i,2))/2+pt1(i,2)
end do
! establish influence coefficients
do i=1,m
    do j=1,m
        ! convert collocation point to local panel coords.
        xt  = co(i,1)  - pt1(j,1)
        zt  = co(i,2)  - pt1(j,2)
        x2t = pt2(j,1) - pt1(j,1)
        z2t = pt2(j,2) - pt1(j,2)
        x   =  xt *cos(th(j)) + zt *sin(th(j))
        z   = -xt *sin(th(j)) + zt *cos(th(j))
        x2  =  x2t*cos(th(j)) + z2t*sin(th(j))
        z2  = 0
        ! save panel lengths for lift coeff. calc.
        if(i.eq.1) then
            dl(j)=x2
        end if
        ! find r1, r2, th1, th2
        r1  = sqrt(x**2+z**2)
        r2  = sqrt((x-x2)**2+z**2)
        th1 = atan2(z,x)
        th2 = atan2(z,x-x2)
        ! compute velocity componants as functions of
        ! gamma1 and gamma2. these velocities are in
        ! the jth reference frame.
        if(i.eq.j) then
            u1l = -0.5*(x-x2)/(x2)
            w1l = -0.15916

            u2l = 0.5*(x)/(x2)
            w2l = 0.15916
        else
            u1l = -(z * log(r2/r1) + (x - x2) * (th2-th1))/(6.28319*x2)
            w1l = -((x2 - z * (th2-th1)) + (x - x2) * log(r2/r1))/(6.28319*x2)

            u2l =  (z * log(r2/r1) + x * (th2-th1))/(6.28319*x2)
            w2l =  ((x2 - z * (th2-th1)) + x * log(r2/r1))/(6.28319*x2)
        end if
    
        ! transform the local velocities into the
        ! global reference frame.
        u1 =  u1l * cos(-th(j)) + w1l * sin(-th(j))
        u2 =  u2l * cos(-th(j)) + w2l * sin(-th(j))
        w1 = -u1l * sin(-th(j)) + w1l * cos(-th(j))
        w2 = -u2l * sin(-th(j)) + w2l * cos(-th(j))
        ! compute the coefficients of gamma in the
        ! influence matrix.
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
    a(i,n+1)=cos(al)*sin(th(i))-sin(al)*cos(th(i))
end do

! add the kutta condition
a(n,1)=1
a(n,n)=1

if(m.eq.10) then
    do i=1,11
        write(6,10) a(i,1),a(i,2),a(i,3),a(i,4),a(i,5),a(i,6),a(i,7), a(i,8),a(i,9),a(i,10),a(i,11)
    end do
end if

n=n+1
! solve for the solution vector of vortex strengths
call matrx(a,n,g)
! convert vortex strengths into tangential
! velocities along the airfoil surface and cp's
! on each of the panels.
200 continue
n=m+1
cl=0
do i=1,m
    vel=0
    do j=1,n
        vel=vel+b(i,j)*g(j)
    end do
    v=vel+cos(al)*cos(th(i))+sin(al)*sin(th(i))
    cl=cl+v*dl(i)
    cp=1-v**2
    write(8,*) co(i,1),' ,',cp
end do

write(6,*) ' '
write(6,*) 'lift coefficient=',cl

stop
10 format(/,f6.2,1x,f5.2,1x,f5.2,1x,f5.2,1x,f5.2,1x,f5.2,1x,f5.2,1x,f5.2,1x,f5.2,1x,f5.2,1x,f5.2)
end

subroutine matrx(a,n,g)
    real, dimension(n,n) :: a, temp
    real, dimension(:) :: g

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
