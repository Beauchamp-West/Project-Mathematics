program main

    implicit none

    integer nx, ny, nz_num, mr, h, i, j, k
    integer, parameter :: itr_max = 1
    real (kind = 8) dx, dy, tol_abs, tol_rel, x, y, p_error
    real (kind = 8), allocatable, dimension(:, :) :: a,b,c,s
    real (kind = 8), allocatable, dimension(:) :: p, p_exact, rhs

    tol_abs = 1.0D-08
    tol_rel = 1.0D-08

    ! write ( *, '(a)' ) 'test'
    open(1, file='task1.txt')

    do h = 1, 3
        nx = 10 * 2**h - 1
        ny = 10 * 2**h - 1
        mr = min(nx*ny, 1000)
        nz_num = 9 * ny * nx
        dx = real(2,8) / (10 * 2**h)
        dy = real(2,8) / (10 * 2**h)

        allocate(a(nx,ny), b(nx,ny), c(nx,ny), s(nx,ny), p(nx*ny), p_exact(nx*ny), rhs(nx*ny))

        a(1:nx,1:ny) = 0.1D+00
        b(1:nx,1:ny) = 0.0D+00
        c(1:nx,1:ny) = 0.1D+00
        s(1:nx,1:ny) = 1.0D+00

        !
        ! Set the exact solution.
        !
        do i = 0, nx*ny-1
            k = i / nx + 1          ! y coordinate
            j = mod(i, nx) + 1      ! x coordinate
            x = j * dx
            y = k * dy - 1
            p_exact(i+1) = x * (x-2) * (y-1) * (y+1)
            rhs(i+1) = (x**2 - 2 * x + y**2 - 1) / 5
        end do
        s = -1 * reshape(rhs,[nx,ny]) ! Only used for test
        p_error = maxval(abs(p_exact))

        ! write ( *, '(a)' ) ' '
        ! write ( *, '(a,g14.6)' ) '  Initial P_ERROR = ', p_error

        call pde_solver (a,b,c,s,nx,ny,nz_num,dx,dy,itr_max,mr,p,tol_abs,tol_rel)

        write(1, '(a,g14.2)') '  h = ', dx
        p_error = maxval(abs(p_exact - p))
        write (1, '(a,g14.6/)' ) '  Final P_ERROR = ', p_error

        deallocate(a,b,c,s,p,p_exact,rhs)
    end do

    close(1)
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
!       -itr_max, the max number of outer iterations.
!       -tol_abs, tol_rel, absolute and relative tolerance.
!   OUTPUT
!       -p(nx*ny), numerical solution of pde in the matrix form.
subroutine pde_solver (a, b, c, s, nx, ny, nz_num, dx, dy, itr_max, mr, p, tol_abs, tol_rel)

    implicit none

    integer, intent(in) :: nx, ny, nz_num, itr_max, mr
    real (kind = 8), intent(in) :: dx, dy, tol_abs, tol_rel
    real (kind = 8), intent(in), dimension(nx, ny) :: a, b, c, s
    real (kind = 8), intent(out), dimension(nx*ny) :: p
    integer i, j, k, i_kj, row, col
    integer, allocatable :: i_lhs(:), j_lhs(:) ! indices of nonzero elements of the lhs matrix 
    real (kind = 8), allocatable :: lhs(:) ! nonzero elements of the lhs matrix
    real (kind = 8), allocatable, dimension(:,:) :: a_pad, b_pad, c_pad ! a,b,c padding with zeros
    real (kind = 8), allocatable, dimension(:) :: rhs, p_estimate

    !
    !  Set the matrix.
    !
    allocate(a_pad(nx+2,ny+2), b_pad(nx+2,ny+2), c_pad(nx+2,ny+2))
    allocate(lhs(nz_num), i_lhs(nz_num), j_lhs(nz_num))

    call add_zeros(a,nx,ny,a_pad)
    call add_zeros(b,nx,ny,b_pad)
    call add_zeros(c,nx,ny,c_pad)

    do i = 0, nz_num-1
        row = i / 9 + 1
        col = mod(i,9) + 1
        k = (row-1) / nx + 1
        j = mod(row-1,nx) + 1

        !
        !  Use compressed row storage for the LHS matrix.
        !
        i_lhs(i+1) = k * (nx+2) + j + 1
        if (col == 1) then
            j_lhs(i+1) = (k-1) * (nx+2) + j
            lhs(i+1) = (b_pad(j+1,k+1) + (b_pad(j,k+1) + b_pad(j+1,k)) / 2) / (4*dx*dy)
        elseif (col == 2) then
            j_lhs(i+1) = (k-1) * (nx+2) + j + 1
            lhs(i+1) = (c_pad(j+1,k+1) + c_pad(j+1,k)) / (2*dy**2) - (b_pad(j+2,k+1) - b_pad(j,k+1)) / (8*dx*dy)
        elseif (col == 3) then
            j_lhs(i+1) = (k-1) * (nx+2) + j + 2
            lhs(i+1) = -1 * (b_pad(j+1,k+1) + (b_pad(j+2,k+1) + b_pad(j+1,k)) / 2) / (4*dx*dy)
        elseif (col == 4) then
            j_lhs(i+1) = k * (nx+2) + j
            lhs(i+1) = (a_pad(j+1,k+1) + a_pad(j,k+1)) / (2*dx**2) - (b_pad(j+1,k+2) - b_pad(j+1,k)) / (8*dx*dy)
        elseif (col == 5) then
            j_lhs(i+1) = k * (nx+2) + j + 1 
            lhs(i+1) = -1 * (a_pad(j+1,k+1) + (a_pad(j+2,k+1) + a_pad(j,k+1)) / 2) / (dx**2) - (c_pad(j+1,k+1) + &
            (c_pad(j+1,k+2) + c_pad(j+1,k)) / 2) / (dy**2)
        elseif (col == 6) then
            j_lhs(i+1) = k * (nx+2) + j + 2
            lhs(i+1) = (a_pad(j+1,k+1) + a_pad(j+2,k+1)) / (2*dx**2) + (b_pad(j+1,k+2) - b_pad(j+1,k)) / (8*dx*dy)
        elseif (col == 7) then
            j_lhs(i+1) = (k+1) * (nx+2) + j
            lhs(i+1) = -1 * (b_pad(j+1,k+1) + (b_pad(j,k+1) + b_pad(j+1,k+2)) / 2) / (4*dx*dy)
        elseif (col == 8) then
            j_lhs(i+1) = (k+1) * (nx+2) + j + 1  
            lhs(i+1) = (c_pad(j+1,k+2) + c_pad(j+1,k+1)) / (2*dy**2) + (b_pad(j+2,k+1) - b_pad(j,k+1)) / (8*dx*dy)
        elseif (col == 9) then
            j_lhs(i+1) = (k-1) * (nx+2) + j + 2
            lhs(i+1) = (b_pad(j+1,k+1) + (b_pad(j+2,k+1) + b_pad(j+1,k+2)) / 2) / (4*dx*dy)
        end if

        ! d(j,j) = -1 * ((a(j,k) + (a(j+1,k)+a(j-1,k))/2) / dx**2 + (c(j,k) + (c(j,k+1)+c(j,k-1))/2) / dy**2)
        ! u(j,j) = 0.5 * ((c(j,k) + c(j,k+1)) / dy**2 + (b(j+1,k) - b(j-1,k)) / (4*dx*dy))
        ! l(j,j) = 0.5 * ((c(j,k) + c(j,k-1)) / dy**2 - (b(j+1,k) - b(j-1,k)) / (4*dx*dy))
        ! d(i,i+1) = (a(i+1,k) + a(i,k)) / (2*(dx**2)) + (b(i,k+1) - b(i,k-1)) / (8*dx*dy)
        ! d(i+1,i) = (a(i,k) + a(i+1,k)) / (2*(dx**2)) - (b(i+1,k+1) - b(i+1,k-1)) / (8*dx*dy)
        ! u(i,i+1) = (b(i,k) + (b(i,k+1) + b(i+1,k)) / 2) / (4*dx*dy)
        ! u(i+1,i) = -1 * (b(i+1,k) + (b(i+1,k+1) + b(i,k)) / 2) / (4*dx*dy)
        ! l(i,i+1) = -1 * (b(i,k) + (b(i+1,k) + b(i,k-1)) / 2) / (4*dx*dy)
        ! l(i+1,i) = (b(i+1,k) + (b(i+1,k-1) + b(i,k)) / 2) / (4*dx*dy)
    end do

    deallocate(a_pad,b_pad,c_pad)

    !
    !  Set the right hand side.
    !
    allocate(rhs((nx+2)*(ny+2)))
    rhs = 0.0D+00

    do i = 0, ny*nx-1
        k = i / nx + 1
        j = mod(i, nx) + 1
        i_kj = k * (nx+2) + j + 1
        rhs(i_kj) = -1 * s(j,k)
    end do

    !
    !  Set the initial solution estimate.
    !
    allocate(p_estimate((nx+2) * (ny+2)))
    p_estimate(:) = 0.0D+00
    
    call mgmres_st ( (nx+2)*(ny+2), nz_num, i_lhs, j_lhs, lhs, p_estimate, rhs, itr_max, mr, &
        tol_abs, tol_rel )

    deallocate(i_lhs, j_lhs, lhs, rhs)

    do k = 1, ny
        do j = 1, nx
            p(j+(k-1)*nx) = p_estimate(k*(nx+2) + j + 1)
        end do
    end do

    deallocate(p_estimate)

    return
end subroutine

! ADD_ZEROS add zeros around the matrix.
!
! Parameters
!
!     -INPUT
!         x,y, size of the matrix.
!         a(x,y), the original matrix.

!     -OUTPUT
!         b(x+1,y+1), the matrix after adding edges.
subroutine add_zeros(a,x,y,b)

    implicit none
    integer :: x, y, i
    real (kind = 8), dimension(x, y) :: a
    real (kind = 8), dimension(x+2, y+2) :: b

    b(1,:) = 0.0D+00
    b(x+2,:) = 0.0D+00
    do i = 1, x
        b(i+1,:) = [0.0D+00, a(i,:), 0.0D+00]
    end do

    return
end subroutine