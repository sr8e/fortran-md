module MD_STRUCT
    use UTILS
    use POTENTIALS, only: f, p

    implicit none

    type Particle
        real mass ! g/mol
        character(2) elem
        real, dimension(3) :: x = 0, v = 0, f = 0, f_next = 0
        real :: ke = 0, pe = 0
    end type Particle

    real, dimension(3) :: box_size
    logical, dimension(3) :: is_periodic ! boundary condition

    contains

    subroutine bulk_set_mass(p_set, mass)
        type(Particle), dimension(:), intent(inout) :: p_set
        real, intent(in) :: mass
        integer index

        do index = 1, size(p_set)
            p_set(index)%mass = mass
        end do
    end subroutine bulk_set_mass

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
                diff = tmp%v(dimen) * dt + acc(tmp%mass, tmp%f(dimen)) * dt ** 2 / 2.0
                p_set(index)%x(dimen) = tmp%x(dimen) + diff
                if (is_periodic(dimen)) then
                    p_set(index)%x(dimen) = periodic(p_set(index)%x(dimen), box_size(dimen))
                end if
            end do
        end do
    end subroutine velo_verlet_x

    subroutine velo_verlet_v(p_set, dt)
        ! update velocity
        type(Particle), dimension(:), intent(inout) :: p_set
        type(Particle) tmp
        integer index, dimen
        real, intent(in) :: dt

        do index = 1, size(p_set)
            tmp = p_set(index)
            do dimen = 1, 3
                p_set(index)%v(dimen) = tmp%v(dimen) + &
                    acc(tmp%mass, (tmp%f(dimen) + tmp%f_next(dimen))/2.0) * dt
            end do
            p_set(index)%ke = convert_unit(tmp%mass * d_euc(tmp%v) ** 2 / 2.0, .true.)
        end do
    end subroutine velo_verlet_v

    subroutine iterate_pair_f(p_set, cutoff)
        type(Particle), dimension(:), intent(inout) :: p_set
        real, intent(in) :: cutoff
        integer i, j
        integer s, t, u
        ! reset potential energy
        do i = 1, size(p_set)
            p_set(i)% pe = 0
        end do

        ! rc must be less than length of shortest edge
        if (cutoff .ge. minval(box_size)) then
            write(0, *) "simulation box is too small (must be larger than r_c)"
            stop
        end if

        do i = 1, size(p_set)
            do j = 1, size(p_set)
                if (i .eq. j) then
                    cycle
                end if
                ! iterate image cell
                do s = -1, 1
                    if (.not. is_periodic(1) .and. s .ne. 0 ) then
                        cycle
                    end if
                    do t = -1, 1
                        if (.not. is_periodic(2) .and. t .ne. 0) then
                            cycle
                        end if
                        do u = -1, 1
                            if (.not. is_periodic(3) .and. u .ne. 0) then
                                cycle
                            end if

                            call add_force(p_set(i), p_set(j), [s, t, u], cutoff)
                        end do
                    end do
                end do
            end do
        end do
    end subroutine iterate_pair_f

    subroutine add_force(pi, pj, cell_ofs, cutoff)
        ! calculate force of particle j to particle i
        type(Particle), intent(inout) :: pi
        type(Particle), intent(in) :: pj
        integer, dimension(3), intent(in) :: cell_ofs
        real, intent(in) :: cutoff
        real, dimension(3) :: distance
        real r
        integer index
        distance = pj%x + box_size * cell_ofs - pi%x
        r = d_euc(distance)
        if (r .gt. cutoff) then
            return
        end if

        pi%pe = pi%pe + p(r) / 2.0
        do index = 1, 3
            pi%f_next(index) = pi%f_next(index) + f(r) * distance(index) / r
        end do
    end subroutine add_force

    function total_energy(p_set)
        real, dimension(2) :: total_energy
        type(Particle), dimension(:) :: p_set
        integer i
        total_energy = 0
        do i = 1, size(p_set)
            total_energy(1) = total_energy(1) + p_set(i)%pe
            total_energy(2) = total_energy(2) + p_set(i)%ke
        end do
    end function total_energy

    !potentials

end module MD_STRUCT
