program main

    implicit none

    ! set m = n = e^(-1), t = 0.
    real (kind = 8) :: m, n, dt, gamma
    m = 1.0
    n = 1.0
    gamma = 1.0
    dt = 0.02

    call l3(m,n,dt,gamma)
    print *, "m(dt/2) = ", m
    print *, "exp(-0.01) = ", exp(-0.01) 
end


subroutine l3 (m,n,dt,gamma)

! L3 numerically solves the ODE system m_t = -(m^2+n^2)^(\gamma-1)*m & 
! n_t = -(m^2+n^2)^(\gamma-1)*n over a small time interval dt/2 with the given initial value of m and n.
!
! Parameters
!   INPUT
!       -m,n, initial values of m and n.
!       -dt, 2 * length of the time interval.
!       -gamma, a constant parameter.
!   OUTPUT
!       -m,n, values of m, n after the time interval.

    implicit none

    real (kind = 8) :: m, n, dt, gamma
    real (kind = 8) :: a, b, at, bt

    a = m**2
    b = n**2
    at = a * exp(-0.5 * dt * ((a+b)**(gamma-1) + ((a+b)**(1-gamma) - dt * (1-gamma))**(-1)))
    bt = b * exp(-0.5 * dt * ((a+b)**(gamma-1) + ((a+b)**(1-gamma) - dt * (1-gamma))**(-1)))

    m = sqrt(at)
    n = sqrt(bt)

    return
end subroutine