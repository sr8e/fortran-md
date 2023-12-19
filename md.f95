program md
    use MD_STRUCT
    use POTENTIALS, only: set_morse_params

    implicit none

    !function
    type(Particle), dimension(2) :: particles
    !params
    real, parameter :: dt = 0.001 !
    integer, parameter :: stepmax = 100
    !potential
    real, parameter :: eps = 0.2703, alpha = 1.1646, r0 = 3.253
    !mass
    real, parameter :: mass = 26.98

    integer step, index
    !total energy
    real, dimension(2) :: toteng  ! pe, ke

    ! i/o
    integer :: ifile = 1
    open(ifile, FILE="trajectory.txt")

    call set_morse_params(eps, alpha, r0)
    call bulk_set_mass(particles, mass)

    ! set init pos
    particles(1)%x(1) = -1.0
    particles(2)%x(1) = 1.0
    ! set init force
    call iterate_pair_f(particles)
    call bulk_copy_f(particles)

    do step = 0, stepmax
        call velo_verlet_x(particles, dt)
        call iterate_pair_f(particles)
        call velo_verlet_v(particles, dt)
        call bulk_copy_f(particles)

        write(ifile, '(i5, ",", *(1pe0.8, :, ","))') step, particles(1)%x, particles(2)%x
        toteng = total_energy(particles)
    end do

    close(ifile)
end program md

