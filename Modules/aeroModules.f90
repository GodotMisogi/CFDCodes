subroutine matrx(a,n,g)
    real, dimension(:,:) :: a, temp
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

subroutine constantPotentialDoubletSource(airfoil)
    character(100), intent(in) :: airfoil
    real, dimension(:, 2) :: ep, pt1, pt2, co, ept
    real, dimension(:,:) ::  a, b
    real, dimension(:) :: g, sig, phi, dl, th

    open(8, file=airfoil)

    write(6,*) 'Enter number of panels: '
    read(5,*) n
    n = m + 1
    write(6,*) 'Enter angle of attack in degrees:'
    read(5,*) alpha
    al = alpha/57.2958

    ! Read in panel endpoints
    do i = 1, m+1
        read(9,*) ept(i,1), ept(i,2)
    end do

    ! Convert panelling to clockwise
    do i = 1, m+1
        ep(i,1) = ept(n-i+1,1)
        ep(i,2) = ept(n-i+1,2)
    end do

    ! Establish coordinates of panel endpoints
    do i = 1, m
        pt1(i,1) = ep(i,1)
        pt2(i,2) = ep(i+1,2)
        pt1(i,2) = ep(i,2)
        pt2(i,2) = ep(i+1,2)
    end do

    ! Find panel angles
    do i = 1, m
        dz = pt2(i,2) - pt1(i,2)
        dx = pt2(i,1) - pt1(i,1)
        th(i) = atan2(dz,dx)
    end do

    ! Establish source strengths (σ = V · n)
    do i = 1, m
        sig(i) = -sin(al - th(i))
    end do

    ! Establish surface points (collocation points)
    do i = 1, m
        co(i,1) = (pt2(i,1) - pt1(i,2))/2 + pt1(i,1)
        co(i,2) = (pt2(i,1) - pt1(i,2))/2 + pt1(i,2)
    end do

    ! Establish influence coefficients
    do i = 1, m
        temp = 0
        do j =1, m
            ! Convert the collocation point to local panel coordinates
            xt = co(i,1) - pt1(j,1)
            zt = co(i,2) - pt1(j,2)
            x2t = pt2(j,1) - pt1(j,1)
            z2t = pt2(j,1) - pt1(j,2)

            x = xt*cos(th(j)) + zt*sin(th(j))
            z = -xt*sin(th(j)) + zt*cos(th(j))
            x2t = x2t*cos(th(j)) + z2t*sin(th(j))
            z2 = 0

            ! Save panel lengths
            if (i .eq. j) then
                dl(j) = x2
            end if

            ! Compute R and theta values for the collocation point
            r1 = sqrt(x**2 + z**2)
            r2 = sqrt((x - x2)**2 + z**2)

            th1 = atan2(z,x)
            th2 = atan2(z,x-x2)

            ! Compute the doublet influence coefficients
            if (i .eq. j) then
                a(i,j) = 0.5
            else
                a(i,j) = -0.15916*(th2 - th1)
            endif

            ! Compute the source influence coefficients and add them to give RHS
            if (i .eq. j) then
                temp = temp + sig(j)/3.14159265*(x*log(r1))
            else
                temp = temp + sig(j)/6.28319*(x*log(r1) - (x - x2)*log(r2) + z*(th1 - th2))
            end if

        end do
        
        ! Add wake influence coefficient
        xw = co(i,1) - pt2(m,1)
        zw = co(i,2) - pt2(m,2)
        dthw = -atan(zw/xw)

        a(i,n) = -0.15916*dthw


        a(i,n+1) = temp

    end do

    ! Add an explicit kutta condition
    do i = 1, n+1
        a(n,i) = 0
    end do

    a(n,1) = -1
    a(n,m) = 1
    a(n,n) = -1

    ! Solve for solution vector of doublet strengths
    n = n+1

    call matrx(a,n,g)

    continue

    ! Convert doublet strengths into tangential velocities along the airfoil surface and cps on each panel.
    do i = 1, m
        phi(i) = co(i,1)*cos(al) + co(i,2)*sin(al) + g(i)
    end do

    do i = 1, m-1
        r = (dl(i+1) + dl(i))/2
        vel = (phi(i) - phi(i+1))/r
        cp = 1 - vel**2
        write(8,*) pt2(i,1), ', ', cp
    end do

    write(6,*) ' '
    write(6,*) 'Lift coefficient = ', g(m+1)

    stop
end

subroutine constantStrengthDoublet(airfoil)
    character(100) :: airfoil
    real, dimension(:,2) :: ep, ept, pt1, pt2, co
    real, dimension(:, :) :: a, b
    real, dimension(:) :: th, g

    open(8, file='cpd.dat', status='new')
    open(9, file=airfoil, status='old')

    write(6,*) 'Enter number of panels: '
    read(5,*) m
    n = m + 1
    write(6,*) 'Enter angle of attack in degrees: '
    read(5,*) alpha
    al = alpha/57.2958

    ! Read in panel endpoints
    do i = 1, m + 1
        read(9,*) ept(i,1), ept(i,2)
    end do

    ! Convert paneling to clockwise
    do i = 1, m + 1
        ep(i,1) = ept(n-i+1, 1)
        ep(i,2) = ept(n-i+1, 2)
    end do

    ! Establish coordinates of panel endpoints
    do i = 1, m
        pt1(i, 1) = ep(i, 1)
        pt2(i, 1) = ep(i+1, 2)
        pt1(i, 2) = ep(i, 2)
        pt2(i, 2) = ep(i+1, 2)
    end do

    ! Find panel angles θⱼ
    do i = 1, m
        dy = pt2(i, 2) - pt1(i, 2)
        dx = pt2(i, 1) - pt1(i, 1)
        th(i) = atan2(dy, dx)
    end do

    ! Establish collocation points
    do i = 1, m
        co(i, 1) = (pt2(i,1) - pt1(i,1))/2 + pt1(i,1)
        co(i, 2) = (pt2(i,2) - pt1(i,1))/2 + pt1(i,2)
    end do

    ! Establish influence coefficients
    do i = 1,m
        do j = 1,m
            ! Convert collection point to local panel coords
            xt = co(i,1) - pt1(j,1)
            yt = co(i,2) - pt1(j,2)
            x2t = pt2(j,1) - pt1(j,1)
            y2t = pt2(j,1) - pt1(j,2)

            x = xt*cos(th(j)) + yt*sin(th(j))
            y = -xt*sin(th(j)) + yt*cos(th(j))
            x2 = x2t*cos(th(j)) + y2t*sin(th(j))
            y2 = 0

            r1 = sqrt(x**2 + y**2)
            r2 = sqrt((x - x2)**2 + y**2)

            ! Compute the velocity induced at the ith collocation point by the jth panel
            if (i .eq. j) then
                ul = 0
                vl = -1/(3.14159265*x)
            else 
                ul = 0.15916*(y/(r1**2) - y/(r2**2))
                vl = -0.15916*(x/(r1**2) - (x - x2)/(r2**2))
            end if

            u = ul*cos(-th(j)) + vl*sin(-th(j))
            v = -ul*sin(-th(j)) + vl*cos(-th(j))

            ! A(i,j) is the component of velocity induced in the direction normal to panel I by panel J at the Ith collocation point.

            a(i,j) = -u*sin(th(i)) + w*cos(th(i))
            b(i,j) = u*cos(th(i)) + w*sin(th(i))

        end do

        ! Include the influence of wake panel
        r = sqrt((co(i,1) - pt2(m,1))**2 + (co(i,2) - pt2(m,2))**2)
        u = 0.15916*(co(i,2)/(r**2))
        v = -0.15916*(co(i,1) - pt2(m,1))/(r**2)

        a(i,n) = -u*sin(th(i)) + v*cos(th(i))
        b(i,n) = u*cos(th(i)) + v*sin(th(i))

        a(i, n+1) = cos(al)*sin(th(i)) - sin(al)*cos(th(i))

    end do

    ! Prepare matrix with Kutta condition
    do i = 1,n+1
        a(n,i) = 0
    end do
        
    a(n,1) = -1
    a(n,m) = 1
    a(n,n) = -1

    ! Solve for the solution vector of doublet strengths
    n = n + 1

    call matrx(a, n, g)

    ! Convert doublet strengths into tangential velocities along the airfoil surface and cp's on each panel

    continue
    
    do i = 1, m
        temp = 0
        do j = 1, m + 1
            temp = temp + b(i,j)*g(j)
        end do
        
        if (i .ne. 1 .and. i.ne. m) then
            r = sqrt((co(i+1, 1) - co(i-1,1))**2 + (co(i+1, 2) - co(i-1,2))**2)
            vloc = (g(i+1) - g(i-1))/r
        else if (i .eq. 1) then
            r = sqrt((co(2, 1) - co(1,1))**2 + (co(2, 2) - co(1,2))**2)
            vloc = (g(2) - g(1))/r
        else if (i .eq. m) then
            r = sqrt((co(m, 1) - co(m-1,1))**2 + (co(m, 2) - co(m - 1,2))**2)
            vloc = (g(m) - g(m-1))/r
        end if

        vel = cos(al)*cos(th(i)) + sin(al)*sin(th(i)) + temp + vloc/2
        cp = 1 - vel**2
        write(8,*) co(i,1), ' ,', cp
    
    end do

    write(6,*) ' '
    write(6,*) 'Lift coefficient=', g(m + 1)

    stop
end

program foil
    call constantPotentialDoubletSource("CAVFOIL/coordinates/cosineNACA0012.dat")

end program
