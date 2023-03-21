program lbd_algorithm
    implicit none
    !gfortran lbd.f90 -Og -g -Wall -Wextra -pedantic -fimplicit-none -fcheck=all -fbacktrace
    character(len=32) :: n_str
    integer :: n
    integer :: i, j
    double precision, dimension(:, :), allocatable :: A, AT, B, tU, tV, UT, VT, U, V, UUT, VVT, eye
    double precision, dimension(:), allocatable :: alpha, betta, v1
    double precision, dimension(:, :), allocatable :: AV, UB
    double precision :: time

    call get_command_argument(1, n_str)
    read(n_str, *) n

    allocate(A(n, n), AT(n, n))
    allocate(v1(n))

    call random_number(A)
    A = A * 100

    call random_number(v1)
    v1 = v1 / norm(v1, n)

    allocate(alpha(n))
    allocate(betta(n))
    alpha = 0.0
    betta = 0.0

    allocate(tU(n, n))
    allocate(tV(n, n + 1))
    tU = 0.0
    tV = 0.0
    do i = 1, n
        tV(i, 1) = v1(i)
    end do
    
    time = lbd(A, alpha, betta, tU, tV, AT, n)
    write(*, *) "Time (msec) = ", time

    allocate(U(n, n), V(n, n))
    do i = 1, n
        do j = 1, n
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

    call matmul(A, V, AV, n)
    call matmul(U, B, UB, n)

    write(*, *) "norm(AV - UB) = ", mat_norm(AV, UB, n)

    allocate(UUT(n, n), VVT(n, n))
    UUT = 0.0
    VVT = 0.0
    call matmul(U, UT, UUT, n)
    call matmul(V, VT, VVT, n)

    allocate(eye(n, n))
    eye = 0.0
    do i = 1, n
        eye(i, i) = 1.0
    end do
    write(*, *) "norm(I - UU^T) = ", mat_norm(eye, UUT, n)
    write(*, *) "norm(I - VV^T) = ", mat_norm(eye, VVT, n)

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

function mat_norm(A, B, n) result(res)
    implicit none
    double precision, allocatable, dimension(:, :), intent(in) :: A, B
    integer, intent(in) :: n
    double precision :: res
    integer :: i, j
    res = 0.0
    do i = 1, n
        do j = 1, n
            res = res + (A(i, j) - B(i, j)) * (A(i, j) - B(i, j))
        end do
    end do
    res = sqrt(res)
end function mat_norm

subroutine matmul(A, B, C, n)
    implicit none
    double precision, allocatable, dimension(:, :), intent(in) :: A, B
    double precision, allocatable, dimension(:, :), intent(inout) :: C
    integer, intent(in) :: n
    integer :: k, j, i

    do j = 1, n
        do k = 1, n
            do i = 1, n
                C(i, j) = C(i, j) + A(i, k) * B(k, j)
            end do
        end do
    end do

end subroutine matmul

function norm(vec, n) result (res)
    implicit none
    double precision :: res
    double precision, dimension(:), intent(in) :: vec
    integer, intent(in) :: n
    integer :: i
    res = 0.0

    do i = 1, n
        res = res + vec(i) * vec(i)
    end do
    res = sqrt(res)
end function norm

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

function dot(a, b, n) result (res)
    implicit none
    double precision :: res
    double precision, dimension(:), intent(in) :: a, b
    integer, intent(in) :: n
    integer :: i
    res = 0.0

    do i = 1, n
        res = res + a(i) * b(i)
    end do
end function dot

subroutine reorthog(M, it, n)
    double precision, allocatable, dimension(:, :), intent(inout) :: M
    integer, intent(in) :: it, n
    integer :: i
    double precision :: a
    do i = 1, it - 1
        a = dot(M(:, i), M(:, it), n)
        do j = 1, n
            M(j, it) = M(j, it) - a * M(j, i)
        end do
    end do

end subroutine reorthog   

subroutine multiply_substract(A, V, U, betta, M, it, n)
    double precision, allocatable, dimension(:, :), intent(in) :: A
    double precision, allocatable, dimension(:, :), intent(inout) :: M
    double precision, dimension(:), intent(in) :: V, U
    double precision, intent(in) :: betta
    integer, intent(in) :: n, it
    integer :: i, j

    do j = 1, n
        do i = 1, n
            M(i, it) = M(i, it) + A(i, j) * V(j)
        end do
    end do

    do i = 1, n
        M(i, it) = M(i, it) - betta * U(i)
    end do

end subroutine multiply_substract

function lbd(A, alpha, betta, U, V, AT, n) result(time)
    double precision, allocatable, dimension(:, :), intent(inout) :: A, AT, U, V
    double precision, allocatable, dimension(:), intent(inout) :: alpha, betta
    integer, intent(in) :: n
    integer :: i, k
    double precision :: t_start, t_finish, time
    double precision, allocatable, dimension(:) :: u0

    call transpose(A, AT, n)
    call cpu_time(t_start)

    do k = 1, n
        if (k == 1) then
            allocate(u0(n))
            u0 = 0.0
            call multiply_substract(A, V(:, k), u0, betta(k), U, k, n)
            deallocate(u0)
        else
            call multiply_substract(A, V(:, k), U(:, k - 1), betta(k - 1), U, k, n)
            !write(*, *) 1
        end if
        call reorthog(U, k, n)
        alpha(k) = norm(U(:, k), n)
        do i = 1, n
            U(i, k) = U(i, k) / alpha(k)
        end do

        call multiply_substract(AT, U(:, k), V(:, k), alpha(k), V, k + 1, n)
        call reorthog(V, k + 1, n)
        betta(k) = norm(V(:, k + 1), n)

        do i = 1, n
            V(i, k + 1) = V(i, k + 1) / betta(k)
        end do
    end do
    call cpu_time(t_finish)
    
    time = t_finish - t_start
    time = time * 1000

end function lbd

end program lbd_algorithm
