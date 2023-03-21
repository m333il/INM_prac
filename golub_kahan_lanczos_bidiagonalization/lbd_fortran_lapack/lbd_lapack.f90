program lbd_algorithm
    implicit none
    !gfortran lbd.f90 -Og -g -Wall -Wextra -pedantic -fimplicit-none -fcheck=all -fbacktrace
    character(len=32) :: n_str
    integer :: n
    double precision, allocatable :: A(:, :), D(:), E(:), tauq(:), taup(:), work(:), lwork_arr(:)
    double precision :: t_start, t_end, time
    integer :: lwork, info
    call get_command_argument(1, n_str)
    read(n_str, *) n

    allocate(A(n, n))
    allocate(D(n))
    allocate(E(n - 1))
    allocate(tauq(n))
    allocate(taup(n))
    lwork = -1
    allocate(lwork_arr(1))
    call random_number(A)

    call dgebrd(n, n, A, n, D, E, tauq, taup, lwork_arr, lwork, info)
    lwork = lwork_arr(1)
    allocate(work(lwork))
    
    call cpu_time(t_start)
    call dgebrd(n, n, A, n, D, E, tauq, taup, work, lwork, info)
    call cpu_time(t_end)
    
    time = (t_end - t_start) * 1000
    write(*, *) "Time (msec) = ", time

    deallocate(A)
    deallocate(D, E)
    deallocate(tauq, taup)
    deallocate(work)
    deallocate(lwork_arr)
end program lbd_algorithm
