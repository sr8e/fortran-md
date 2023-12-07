module MD_STRUCT
    implicit none
    type Particle
        real :: mass = 26.98 ! g/mol
        real, dimension(3) :: x = 0, v = 0, f = 0, f_next = 0
    end type Particle

    contains

    subroutine bulk_copy_f(p_set)

        type(Particle), dimension(:), intent(inout) :: p_set
        integer index
        do index = 1, size(p_set)
            p_set(index)%f = p_set(index)%f_next
            p_set(index)%f_next = 0
        end do
    end subroutine bulk_copy_f

    subroutine velo_verlet_x(p_set, dt)
        ! update coordinate

        type(Particle), dimension(:), intent(inout) :: p_set
        type(Particle) tmp
        integer index, dimen
        real diff
        real, intent(in) :: dt 

        do index=1, size(p_set)
            tmp = p_set(index)
            do dimen=1, 3
                diff = tmp%v(dimen) * dt + acc(tmp%mass, tmp%f(dimen)) * dt ** 2 / 2
                p_set(index)%x(dimen) = tmp%x(dimen) + diff
            end do
        end do
    end subroutine velo_verlet_x
                
    subroutine velo_verlet_v(p_set, dt)
        ! update velocity

        type(Particle), dimension(:), intent(inout) :: p_set
        type(Particle) tmp
        integer index, dimen
        real diff
        real, intent(in) :: dt

        do index = 1, size(p_set)
            tmp = p_set(index)
            do dimen = 1, 3
                p_set(index)%v(dimen) = tmp%v(dimen) + &
                    (tmp%f(dimen) + tmp%f_next(dimen)) * dt / tmp%mass / 2
            end do
        end do
    end subroutine velo_verlet_v

    subroutine iterate_pair_f(p_set)

        type(Particle), dimension(:), intent(inout) :: p_set
        integer i, j
        
        do i = 1, size(p_set) - 1
            do j = i + 1, size(p_set)
                call add_force(p_set(i), p_set(j))
                call add_force(p_set(j), p_set(i))
            end do
        end do
    end subroutine iterate_pair_f

    subroutine add_force(pi, pj)

        type(Particle), intent(inout) :: pi
        type(Particle), intent(in) :: pj
        real, dimension(3) :: distance
        real r
        integer index
        distance = pj%x - pi%x
        r = d_euc(distance)
        do index = 1, 3
            pi%f_next(index) = pi%f_next(index) + f_morse(r) * distance(index) / r
        end do
    end subroutine add_force

    function d_euc(d_vec)
        real, dimension(3) :: d_vec
        real d_euc 
        d_euc = sqrt(d_vec(1) ** 2 + d_vec(2) ** 2 + d_vec(3) ** 2)
    end function d_euc

    function f_morse(r)
        ! unit: eV/A
        real, parameter :: eps = 0.2703, alpha = 1.1646, r0 = 3.253
        real r, f_morse
        f_morse = -2.0 * alpha * eps * ( &
        exp(-2.0 * alpha * (r - r0)) - exp(-alpha * (r - r0)) &
        )
    end function f_morse

    function acc(mass, force)
        ! unit: mass -> g/mol, force -> eV/A, acc -> A/ps^-2
        ! 1 eV = e [J] = e [kg (m/s)^2] = 10^-1 e [g (A/ps)^2]

        real, parameter :: NA = 6.02214076e+23, e = 1.6021892e-19
        real mass, force, acc
        acc = force / mass * e * NA / 10.0
    end function acc

end module MD_STRUCT

program md
    use MD_STRUCT

    implicit none

    !function
    real f, a
    type(Particle), dimension(2) :: particles
    !params
    real, parameter :: dt = 0.001 !
    integer, parameter :: stepmax = 100
    integer step, index

    ! i/o
    !open(ifile, FILE="")
    !read(ifile, *) line
    !close(ifile)

    ! set init pos
    particles(1)%x(1) = -1.0
    particles(2)%x(1) = 1.0
    ! set init force
    call iterate_pair_f(particles)
    call bulk_copy_f(particles)

    do step = 0, stepmax
        call velo_verlet_x(particles, dt)
        call iterate_pair_f(particles)
        call bulk_copy_f(particles)
        call velo_verlet_v(particles, dt)
        write(6, *) step, particles(1)%x, particles(2)%x
    end do

end program md

