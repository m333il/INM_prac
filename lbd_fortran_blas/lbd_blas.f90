program lbd_algorithm
    implicit none
    double precision, external :: dnrm2, ddot
    external :: dcopy, dgemv
    !gfortran lbd.f90 -Og -g -Wall -Wextra -pedantic -fimplicit-none -fcheck=all -fbacktrace
    character(len=32) :: n_str
    integer :: n
    integer :: i, j
    double precision, dimension(:, :), allocatable :: A, AT, B, tU, tV, UT, VT, U, V, UUT, VVT, eye
    double precision, dimension(:), allocatable :: alpha, betta, v1
    double precision, dimension(:, :), allocatable :: AV, UB
    double precision :: time, pnorm

    call get_command_argument(1, n_str)
    read(n_str, *) n

    allocate(A(n, n), AT(n, n))
    allocate(v1(n))

    call random_number(A)
    A = A * 100

    call random_number(v1)
    pnorm = dnrm2(n, v1, 1)
    v1 = v1 / pnorm

    allocate(alpha(n))
    allocate(betta(n))
    alpha = 0.0
    betta = 0.0

    allocate(tU(n, n))
    allocate(tV(n, n + 1))
    tU = 0.0
    tV = 0.0
    tV(:, 1) = v1

    time = lbd(A, alpha, betta, tU, tV, AT, n)
    write(*, *) "Time (msec) = ", time

    allocate(U(n, n), V(n, n))
    do j = 1, n
        do i = 1, n
            U(i, j) = tU(i, j)
            V(i, j) = tV(i, j)
        end do
    end do

    allocate(UT(n, n), VT(n, n))
    call transpose(U, UT, n)
    call transpose(V, VT, n)
    
    allocate(B(n, n))
    B = 0.0
    do i = 1, n - 1
        B(i, i) = alpha(i)
        B(i, (i + 1)) = betta(i)
    end do
    B(n, n) = alpha(n)

    allocate(AV(n, n), UB(n, n))
    AV = 0.0
    UB = 0.0
    call dgemm('N', 'N', n, n, n, 1.0d0, A, n, V, n, 0.0d0, AV, n)
    call dgemm('N', 'N', n, n, n, 1.0d0, U, n, B, n, 0.0d0, UB, n)
    call daxpy(n * n, -1.0d0, AV, 1, UB, 1)
    write(*, *) "norm(AV - UB) = ", dnrm2(n * n, UB, 1)

    allocate(UUT(n, n), VVT(n, n))
    UUT = 0.0
    VVT = 0.0
    call dgemm('N', 'T', n, n, n, 1.0d0, U, n, U, n, 0.0d0, UUT, n)
    call dgemm('N', 'T', n, n, n, 1.0d0, V, n, V, n, 0.0d0, VVT, n)
    do i = 1, n
        UUT(i, i) = UUT(i, i) - 1.0
        VVT(i, i) = VVT(i, i) - 1.0
    end do
    allocate(eye(n, n))
    eye = 0.0
    do i = 1, n
        eye(i, i) = 1.0
    end do
    write(*, *) "norm(I - UU^T) = ", dnrm2(n * n, UUT, 1)
    write(*, *) "norm(I - VV^T) = ", dnrm2(n * n, VVT, 1)

    deallocate(A, AT)
    deallocate(B, alpha, betta)
    deallocate(U, V)
    deallocate(UT, VT)
    deallocate(tU, tV)
    deallocate(AV, UB)
    deallocate(v1)
    deallocate(UUT, VVT)
    deallocate(eye)
    
contains

subroutine transpose(A, AT, n)
    implicit none
    double precision, allocatable, dimension(:, :), intent(in) :: A
    double precision, allocatable, dimension(:, :), intent(inout) :: AT
    integer, intent(in) :: n
    integer :: i, j

    do i = 1, n
        do j = 1, n
            AT(i, j) = A(j, i)
        end do
    end do

end subroutine transpose

subroutine reorthog(M, it, n)
    double precision, allocatable, dimension(:, :), intent(inout) :: M
    integer, intent(in) :: it, n
    integer :: i
    double precision :: a
    do i = 1, it - 1
        a = ddot(n, M(:, i), 1, M(:, it), 1)
        call daxpy(n, -a, M(:, i), 1, M(:, it), 1)
    end do

end subroutine reorthog   

function lbd(A, alpha, betta, U, V, AT, n) result(time)
    double precision, allocatable, dimension(:, :), intent(inout) :: A, AT, U, V
    double precision, allocatable, dimension(:), intent(inout) :: alpha, betta
    integer, intent(in) :: n
    integer :: i, k
    double precision :: t_start, t_finish, time
    double precision, allocatable, dimension(:) :: r, p
    
    allocate(r(n), p(n))
    r = 0.0
    p = 0.0
    call cpu_time(t_start)
 
    do k = 1, n
        if (k == 1) then
            call dgemv('N', n, n, 1.0d0, A, n, V(:, k), 1, 0.0d0, r, 1)
        else
            call dcopy(n, U(:, k - 1), 1, r, 1)
            call dgemv('N', n, n, 1.0d0, A, n, V(:, k), 1, -betta(k - 1), r, 1)
        end if
        U(:, k) = r
        call reorthog(U, k, n)
        alpha(k) = dnrm2(n, U(:, k), 1)
        U(:, k) = U(:, k) / alpha(k)

        p = V(:, k)
        call dcopy(n, V(:, k), 1, p, 1)
        call dgemv('T', n, n, 1.0d0, A, n, U(:, k), 1, -alpha(k), p, 1)
        V(:, k + 1) = p
        call reorthog(V, k + 1, n)
        betta(k) = dnrm2(n, V(:, k + 1), 1)
        V(:, k + 1) = V(:, k + 1) / betta(k)
    end do
    call cpu_time(t_finish)
    
    time = t_finish - t_start
    time = time * 1000
    deallocate(r, p)
end function lbd

end program lbd_algorithm
