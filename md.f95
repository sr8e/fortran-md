program md
    use MD_STRUCT
    use POTENTIALS, only: set_morse_params
    use IO

    implicit none

    !function
    type(Particle), dimension(:), allocatable :: particles
    !params
    real, parameter :: dt = 0.001 !
    integer, parameter :: stepmax = 100
    integer, parameter :: seed = 1891109
    !potential
    real, parameter :: eps = 0.2703, alpha = 1.1646, r0 = 3.253
    real, parameter :: cutoff = 12.0
    !mass
    real, parameter :: mass = 26.98

    integer step, index
    !total energy
    real, dimension(3) :: toteng  ! pe, ke, temp

    ! i/o
    integer :: itraj = 1, ieng = 2

    integer :: iload = 10, idump = 11
    open(iload, FILE="fcc_al.atom")
    call read_atom(iload, particles, box_size)
    close(iload)

    open(itraj, FILE="trajectory.txt")
    open(ieng, FILE="energy.txt")

    !module initialize
    !random
    call setup_random(seed)
    !simulation box
    is_periodic = [.true., .true., .true.]
    !potentials
    call set_morse_params(eps, alpha, r0)

    ! set init vel
    call init_velocity_temp(particles, 300.0)
    ! set init force
    call iterate_pair_f(particles, cutoff)
    call bulk_copy_f(particles)

    do step = 0, stepmax
        call velo_verlet_x(particles, dt)
        call iterate_pair_f(particles, cutoff)
        call velo_verlet_v(particles, dt)
        call bulk_copy_f(particles)

        toteng = total_energy(particles)
        write(ieng, '(i5, ",", *(1pe0.8, :, ","))') step, toteng, sum(toteng(1:2))
    end do

    close(itraj)
    close(ieng)

    open(idump, FILE="al_after.cfg")
    call write_cfg(idump, particles, box_size)
    close(idump)

end program md

