program main
    implicit none
    integer p, i, j, l, flag
    real h, abs_err, rel_err
    real, allocatable :: a(:,:), x(:), phi(:), f(:)
    integer, allocatable :: v(:)
    
    do p = 1, 14
        h = 1.0/real(2**p)
        l = 2**p - 1
        allocate(a(l,l), x(l), f(l), v(l), phi(l))

        do i = 1, l
            x(i) = i * h
            f(i) = (1+ 4*x(i) + 2*x(i)**2 - x(i)**4) * exp(x(i))
            a(i,i) = (2 + (h*x(i))**2)
            v(i) = i
            phi(i) = (1-x(i)**2) * exp(x(i))
        end do

        if (l > 1) then
            do j = 1, l-1
                a(j,j+1) = -1
                a(j+1,j) = -1
            end do
        end if

        f(1) = f(1) + (l+1)**2
        a = a / h**2

        call sgesv(l,1,a,l,v,f,l,flag)
        
        abs_err = maxval(abs(phi-f))
        rel_err = abs_err / h**2
        ! write(*,'(es20.6)') h
        write ( *, '(a,es14.4,a,es14.4,a,es14.4)' ) '  h = ', h, ',  abs_err = ', abs_err, ',  rel_err = ', rel_err
        deallocate(a,x,v,f,phi)
    end do


end