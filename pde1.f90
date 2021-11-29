program main

    implicit none

    integer, parameter :: rk = kind ( 1.0D+00 )
    integer nx, ny, nz_num, mr
    integer, parameter :: itr_max = 1
    real (kind = 8) dx, dy
    real (kind = 8), allocatable, dimension(:, :) :: a,b,c,s,p
    integer h, x_bd, y_bd

    do h = 1, 5
        nx = 10 * 2**h - 1
        ny = 10 * 2**h - 1
        mr = min(nx*ny, 1000)
        nz_num = (3 * ny - 2) * (3 * nx - 2)
        dx = real(2,8) / (10 * 2**h)
        dy = real(2,8) / (10 * 2**h)
        x_bd = int(0.15 * nx)
        y_bd = int(0.0125/2 * ny)

        ! write(*, '(/a,2g14.6)') 'dx & dy: ', dx, dy

        allocate(a(nx,ny), b(nx,ny), c(nx,ny), s(nx,ny), p(nx,ny))

        a(1:nx,1:ny) = 0.1D+00
        b(1:nx,1:ny) = 0.0D+00
        c(1:nx,1:ny) = 0.1D+00
        a(1:x_bd,floor(real(ny)/2 - y_bd):ceiling(real(ny)/2 + y_bd)) = 1.1D+00
        s(1:nx,1:ny) = 1.0D+00

        call pde_solver (a,b,c,s,nx,ny,nz_num,dx,dy,itr_max,mr,p,x_bd,y_bd)

        deallocate(a,b,c,s,p)
    end do

end

! PDE_SOLVER solves the pde in the first task numerically using GMRES method.
!
! Parameters
!   INPUT
!       -a(nx,ny), b(nx,ny), c(nx,ny), s(nx,ny), matrices of coefficients in the pde.
!       -nz_num = (3 * ny - 2) * (3 * nx - 2), the max number of nonzero elements of the lhs matrix.
!       -nx, ny, numbers of meshes in x and y coordinates.
!       -dx, dy, mesh size in x and y coordinate.
!       -mr, the max number of inner iterations. 0 < mr <= nx*ny
!   OUTPUT
!       -p(nx,ny), numerical solution of pde in the matrix form.
subroutine pde_solver (a, b, c, s, nx, ny, nz_num, dx, dy, itr_max, mr, p, x_bd, y_bd)

    implicit none

    integer, intent(in) :: nx, ny, nz_num, itr_max, mr, x_bd, y_bd
    real (kind = 8), intent(in) :: dx, dy
    real (kind = 8), intent(in), dimension(nx, ny) :: a, b, c, s
    real (kind = 8), intent(out), dimension(nx, ny) :: p
    integer i, j, k, m
    integer exact_num, num
    integer, allocatable :: i_lhs(:), j_lhs(:), i_raw(:), j_raw(:) ! indices of nonzero elements of the lhs matrix 
    real (kind = 8), allocatable :: lhs(:), lhs_raw(:) ! nonzero elements of the lhs matrix
    real (kind = 8), allocatable, dimension(:,:) :: t, v, u ! submatrices in the k th row
    real (kind = 8), allocatable, dimension(:) :: rhs, p_estimate, p_exact
    integer test
    real ( kind = 8 ) tol_abs
    real ( kind = 8 ) tol_rel
    real ( kind = 8 ) p_error, rel_error

!
!  Set the matrix.
!
    allocate(t(nx,nx), v(nx,nx), u(nx,nx))
    allocate(lhs_raw(nz_num), i_raw(nz_num), j_raw(nz_num))
    num = 1 ! NUM counts the number of elements on the main diagonal, superdiagonal and subdiagonal.
    exact_num = 0 ! EXACT_NUM counts the number of nonzero elements.

    do k = 1, ny
        do j = 1, nx

            t(j,j) = -1 * (a(j,k) / dx**2 + c(j,k) / dy**2)
            u(j,j) = 0.5 * c(j,k) / dy**2
            v(j,j) = 0.5 * c(j,k) / dy**2
            
            if (j /= 1) then 
                t(j,j) = t(j,j) - a(j-1,k) / (2 * dx**2)
                u(j,j) = u(j,j) - b(j-1,k) / (8*dx*dy)
                v(j,j) = v(j,j) + b(j-1,k) / (8*dx*dy)
            end if
            if (k /= 1) then
                t(j,j) = t(j,j) - c(j,k-1)/ (2 * dy**2)
                v(j,j) = v(j,j) + c(j,k-1) / (2 * dy**2)
            end if
            if (j /= nx) then
                t(j,j) = t(j,j) - a(j+1,k) / (2 * dx**2)
                u(j,j) = u(j,j) + b(j+1,k) / (8*dx*dy)
                v(j,j) = v(j,j) - b(j+1,k) / (8*dx*dy)
            end if
            if (k /= ny) then
                t(j,j) = t(j,j) - c(j,k+1) / (2 * dy**2)
                u(j,j) = u(j,j) + c(j,k+1) / (2 * dy**2)
            end if
            !     t(j,j) = -1 * ((a(j,k) + (a(j+1,k)+a(j-1,k))/2) / dx**2 + (c(j,k) + (c(j,k+1)+c(j,k-1))/2) / dy**2)
            !     u(j,j) = 0.5 * ((c(j,k) + c(j,k+1)) / dy**2 + (b(j+1,k) - b(j-1,k)) / (4*dx*dy))
            !     v(j,j) = 0.5 * ((c(j,k) + c(j,k-1)) / dy**2 - (b(j+1,k) - b(j-1,k)) / (4*dx*dy))
        end do

        do i = 1, nx-1

            t(i,i+1) = (a(i+1,k) + a(i,k)) / (2*(dx**2))
            t(i+1,i) = (a(i,k) + a(i+1,k)) / (2*(dx**2))
            u(i,i+1) = (b(i,k) + b(i+1,k) / 2) / (4*dx*dy)
            u(i+1,i) = -1 * (b(i+1,k) + b(i,k) / 2) / (4*dx*dy)
            v(i,i+1) = -1 * (b(i,k) + b(i+1,k) / 2) / (4*dx*dy)
            v(i+1,i) = (b(i+1,k) + b(i,k) / 2) / (4*dx*dy)

            if (k /= 1) then
                t(i,i+1) = t(i,i+1) - b(i,k-1) / (8*dx*dy)
                t(i+1,i) = t(i+1,i) + b(i+1,k-1) / (8*dx*dy)
                v(i,i+1) = v(i,i+1) - b(i,k-1) / (8*dx*dy)
                v(i+1,i) = v(i+1,i) + b(i+1,k-1) / (8*dx*dy)
            end if
            if (k /= ny) then
                t(i,i+1) = t(i,i+1) + b(i,k+1) / (8*dx*dy)
                t(i+1,i) = t(i+1,i) - b(i+1,k+1) / (8*dx*dy)
                u(i,i+1) = u(i,i+1) + b(i,k+1) / (8*dx*dy)
                u(i+1,i) = u(i+1,i) - b(i+1,k+1) / (8*dx*dy)
            end if
            ! t(i,i+1) = (a(i+1,k) + a(i,k)) / (2*(dx**2)) + (b(i,k+1) - b(i,k-1)) / (8*dx*dy)
            ! t(i+1,i) = (a(i,k) + a(i+1,k)) / (2*(dx**2)) - (b(i+1,k+1) - b(i+1,k-1)) / (8*dx*dy)
            ! u(i,i+1) = (b(i,k) + (b(i,k+1) + b(i+1,k)) / 2) / (4*dx*dy)
            ! u(i+1,i) = -1 * (b(i+1,k) + (b(i+1,k+1) + b(i,k)) / 2) / (4*dx*dy)
            ! v(i,i+1) = -1 * (b(i,k) + (b(i+1,k) + b(i,k-1)) / 2) / (4*dx*dy)
            ! v(i+1,i) = (b(i+1,k) + (b(i+1,k-1) + b(i,k)) / 2) / (4*dx*dy)
        end do

        do m = 1, nx
            i = (k-1) * nx + m
            i_raw(num) = i

            if (k /= 1) then
                if (m /= 1) then
                    j_raw(num) = (k-2) * nx + m - 1
                    lhs_raw(num) = v(m,m-1)
                    if (abs(lhs_raw(num)) > 1.0D-08) then
                        exact_num = exact_num + 1
                    end if
                    num = num + 1
                    i_raw(num) = i
                end if
                j_raw(num) = (k-2) * nx + m
                lhs_raw(num) = v(m,m)
                if (abs(lhs_raw(num)) > 1.0D-08) then
                    exact_num = exact_num + 1
                end if                
                num = num + 1
                i_raw(num) = i
                if (m /= nx) then
                    j_raw(num) = (k-2) * nx + m + 1
                    lhs_raw(num) = v(m,m+1)
                    if (abs(lhs_raw(num)) > 1.0D-08) then
                        exact_num = exact_num + 1
                    end if                    
                    num = num + 1
                    i_raw(num) = i
                end if
            end if

            if (m /= 1) then
                j_raw(num) = (k-1) * nx + m - 1
                lhs_raw(num) = t(m,m-1)
                if (abs(lhs_raw(num)) > 1.0D-08) then
                    exact_num = exact_num + 1
                end if                
                num = num + 1
                i_raw(num) = i
            end if
            j_raw(num) = (k-1) * nx + m
            lhs_raw(num) = t(m,m)
            if (abs(lhs_raw(num)) > 1.0D-08) then
                exact_num = exact_num + 1
            end if               
            num = num + 1
            i_raw(num) = i
            if (m /= nx) then
                j_raw(num) = (k-1) * nx + m + 1
                lhs_raw(num) = t(m,m+1)
                if (abs(lhs_raw(num)) > 1.0D-08) then
                    exact_num = exact_num + 1
                end if                   
                num = num + 1
                i_raw(num) = i
            end if 
            
            if (k /= ny) then
                if (m /= 1) then
                    j_raw(num) = k * nx + m - 1
                    lhs_raw(num) = u(m,m-1)
                    if (abs(lhs_raw(num)) > 1.0D-08) then
                        exact_num = exact_num + 1
                    end if                       
                    num = num + 1
                    i_raw(num) = i
                end if
                j_raw(num) = k * nx + m
                lhs_raw(num) = u(m,m)
                if (abs(lhs_raw(num)) > 1.0D-08) then
                    exact_num = exact_num + 1
                end if                   
                num = num + 1
                i_raw(num) = i
                if (m /= nx) then
                    j_raw(num) = k * nx + m + 1
                    lhs_raw(num) = u(m,m+1)
                    if (abs(lhs_raw(num)) > 1.0D-08) then
                        exact_num = exact_num + 1
                    end if                       
                    num = num + 1
                    i_raw(num) = i
                end if                
            end if

        end do
    end do

    deallocate(u,v,t)

!
!  Use compressed row storage for the LHS matrix.
!

    allocate(lhs(exact_num),i_lhs(exact_num), j_lhs(exact_num))
    k = 1
    do m = 1, nz_num
        if(abs(lhs_raw(m)) > 1.0D-08) then
            lhs(k) = lhs_raw(m)
            i_lhs(k) = i_raw(m)
            j_lhs(k) = j_raw(m)
            k = k + 1
        end if
    end do

    deallocate(lhs_raw, i_raw, j_raw)
!
!  Set the right hand side. (This part has no use in the test stage)
!
    allocate(rhs(nx*ny))

    do i = 1, ny*nx
        j = i / nx + 1
        k = mod(i, nx)
        rhs(i) = -1 * s(k,j)
    end do
!
!  Set the exact solution and corresponding rhs vector.
!   p(x,y) = x^2 + y^2    
!
    allocate(p_exact(nx*ny))

    do i = 1, nx*ny
        j = i / nx + 1
        k = mod(i, nx) 
        p_exact(i) = (k * dx)**2 + (j * dy - 1)**2
    end do

    rhs(1:nx*ny) = 0.4
    do i = 1, x_bd
        do j = ny/2 - y_bd, ny/2 + y_bd
            rhs(j*nx+i) = rhs(j*nx+i) + 2
        end do
    end do

!
!  Set the initial solution estimate.
!
    test = 1
    allocate(p_estimate(nx*ny))
    p_estimate(:) = 0.0D+00
    p_error = maxval ( abs ( ( p_exact(:) - p_estimate(:) )**2 ) )

    tol_abs = 1.0D-08
    tol_rel = 1.0D-08

    write ( *, '(a)' ) ' '
    ! write ( *, '(a,i8)' ) '  Test ', test
    ! write ( *, '(a,i8)' ) '  Matrix order N = ', nx*ny
    ! write ( *, '(a,i8)' ) '  Inner iteration limit = ', mr
    ! write ( *, '(a,i8)' ) '  Outer iteration limit = ', itr_max
    write ( *, '(a,g14.6)' ) '  Initial P_ERROR = ', p_error

    call mgmres_st ( nx*ny, exact_num, i_lhs, j_lhs, lhs, p_estimate, rhs, itr_max, mr, &
        tol_abs, tol_rel )

    deallocate(i_lhs, j_lhs, lhs, rhs)

    write(*, '(a,g14.6)') '  h^2 = ', dx**2
    p_error = maxval ( abs ( ( p_exact(:) - p_estimate(:) )**2 ) )
    write ( *, '(a,g14.6)' ) '  Final P_ERROR = ', p_error

    rel_error = p_error / sqrt(sum(p_exact(:)**2))
    ! write ( *, '(a,g14.6)' ) '  Relative P_ERROR = ', rel_error

    do i = 1, nx
        do j = 1, ny
            p(i,j) = p_estimate((j-1)*nx + i)
        end do
    end do

    deallocate(p_exact, p_estimate)

    return
end