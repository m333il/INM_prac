program ldb_algorithm
    implicit none
    !gfortran ldb.f90 -Og -g -Wall -Wextra -pedantic -fimplicit-none -fcheck=all -fbacktrace
    character(len=32) :: n_str
    integer :: n
    integer :: i, j
    double precision, dimension(:), allocatable :: A, AT, B, alpha, betta, tU, tV, UT, VT, U, V, v1
    double precision, dimension(:), allocatable :: AV, UB
    double precision :: time

    call get_command_argument(1, n_str)
    read(n_str, *) n

    allocate(A(0 : n * n - 1), AT(0 : n * n - 1))
    allocate(v1(0 : n - 1))

    call random_number(A)
    A = A * 100

    call random_number(v1)
    v1 = v1 / norm(v1, 0, n)

    allocate(alpha(0 : n))
    allocate(betta(0 : n))
    alpha = 0.0
    betta = 0.0

    allocate(tU(0 : n * (n + 1) - 1))
    allocate(tV(0 : n * (n + 2) - 1))
    tU = 0.0
    tV = 0.0
    do i = 0, n - 1
        tV(n + i) = v1(i)
    end do
    
    time = lbd(A, alpha, betta, tU, tV, AT, n)
    write(*, *) time

    allocate(U(0 : n * n - 1), V(0 : n * n - 1))
    do i = 0, n - 1
        do j = 0, n - 1
            U(i + j * n) = tU(i + (j + 1) * n)
            V(i + j * n) = tV(i + (j + 1) * n)
        end do
    end do

    allocate(UT(0 : n * n - 1), VT(0 : n * n - 1))
    call transpose(U, UT, n)
    call transpose(V, VT, n)
    
    allocate(B(0 : n * n - 1))
    B = 0.0
    do i = 0, n - 2
        B(i * n + i) = alpha(i + 1)
        B(i + (i + 1) * n) = betta(i + 1)
    end do
    B(n - 1 + (n - 1) * n) = alpha(n)

    allocate(AV(0 : n * n - 1), UB(0 : n * n - 1))
    AV = 0.0
    UB = 0.0

    call matmul(A, V, AV, n)
    call matmul(U, B, UB, n)

    write(*, *) mat_norm(AV, UB, n)

    deallocate(A, AT)
    deallocate(B, alpha, betta)
    deallocate(U, V)
    deallocate(UT, VT)
    deallocate(tU, tV)
    deallocate(AV, UB)
    deallocate(v1)
    
contains

function mat_norm(A, B, n) result(res)
    implicit none
    double precision, allocatable, dimension(:), intent(in) :: A, B
    integer, intent(in) :: n
    double precision :: res
    integer :: i, j
    res = 0.0
    do i = 0, n - 1
        do j = 0, n - 1
            res = res + (A(i + j * n) - B(i + j * n)) * (A(i + j * n) - B(i + j * n))
        end do
    end do
    res = sqrt(res)
end function mat_norm

subroutine matmul(A, B, C, n)
    implicit none
    double precision, allocatable, dimension(:), intent(in) :: A, B
    double precision, allocatable, dimension(:), intent(inout) :: C
    integer, intent(in) :: n
    integer :: k, j, i

    do j = 0, n - 1
        do k = 0, n - 1
            do i = 0, n - 1
                C(i + j * n) = C(i + j * n) + A(i + k * n) * B(k + j * n)
            end do
        end do
    end do

end subroutine matmul

function norm(vec, lvec, n) result (res)
    implicit none
    double precision :: res
    double precision, allocatable, dimension(:), intent(in) :: vec
    integer, intent(in) :: n, lvec
    integer :: i
    res = 0.0

    do i = 0, n - 1
        res = res + vec(lvec + i) * vec(lvec + i)
    end do
    res = sqrt(res)
end function norm

subroutine transpose(A, AT, n)
    implicit none
    double precision, allocatable, dimension(:), intent(in) :: A
    double precision, allocatable, dimension(:), intent(inout) :: AT
    integer, intent(in) :: n
    integer :: i, j

    do i = 0, n - 1
        do j = 0, n - 1
            AT(i * n + j) = A(j * n + i)
        end do
    end do

end subroutine transpose

subroutine multiply_substract(A, V, lv, U, lu, betta, vec, lvec, n)
    double precision, allocatable, dimension(:), intent(in) :: A, V, U
    double precision, intent(in) :: betta
    integer, intent(in) :: n, lv, lu, lvec
    double precision, allocatable, dimension(:), intent(inout) :: vec
    integer :: i, j

    do j = 0, n - 1
        do i = 0, n - 1
            vec(lvec + i) = vec(lvec + i) + A(i + j * n) * V(lv + j)
        end do
    end do

    do i = 0, n - 1
        vec(lvec + i) = vec(lvec + i) - betta * U(lu + i)
    end do

end subroutine multiply_substract

function lbd(A, alpha, betta, U, V, AT, n) result(time)
    double precision, allocatable, dimension(:), intent(inout) :: A, AT, alpha, betta, U, V
    integer, intent(in) :: n
    integer :: i, k
    double precision :: t_start, t_finish, time

    call transpose(A, AT, n)
    call cpu_time(t_start)

    do k = 1, n
        call multiply_substract(A, V, k * n, U, (k - 1) * n, betta(k - 1), U, k * n, n)
        alpha(k) = norm(U, k * n, n)
        !write(*, *) alpha(k)
        do i = 0, n - 1
            U(i + k * n) = U(i + k * n) / alpha(k)
        end do

        call multiply_substract(AT, U, k * n, V, k * n, alpha(k), V, (k + 1) * n, n)
        betta(k) = norm(V, (k + 1) * n, n)

        do i = 0, n - 1
            V(i + (k + 1) * n) = V(i + (k + 1) * n) / betta(k)
        end do
    end do
    call cpu_time(t_finish)
    
    time = t_finish - t_start

end function lbd

end program ldb_algorithm
