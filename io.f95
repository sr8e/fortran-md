module IO
    use MD_STRUCT
    implicit none

    contains

    subroutine read_atom(ifile, p_set, box_size)
        integer, intent(in) :: ifile
        type(Particle), dimension(:), allocatable, intent(out) :: p_set
        real, dimension(3), intent(out) :: box_size
        integer i, n

        read(ifile, *) n

        ! allocate
        allocate(p_set(n))

        read(ifile, *) box_size

        do i = 1, n
            read(ifile, *) p_set(i)%mass, p_set(i)%elem, p_set(i)%x, p_set(i)%v
        end do
    end subroutine read_atom

    subroutine write_cfg(ifile, p_set, box_size)
        integer, intent(in) :: ifile
        type(Particle), dimension(:), intent(in) :: p_set
        real, dimension(3), intent(in) :: box_size
        integer i, j
        real len
        type(Particle) p
        real, dimension(3) :: reg_x

        write(ifile, '("Number of particles = ", i0)') size(p_set)
        write(ifile, '("A = 1 Angstrom")')
        do i = 1, 3
            do j = 1, 3
                if (i == j) then
                    len = box_size(i)
                else
                    len = 0
                end if
                write(ifile, '("H0(", i0, ",", i0, ") = ", e0.8, " A")') i, j, len
            end do
        end do
        write(ifile, '("entry_count = 6")')

        do i = 1, size(p_set)
            p = p_set(i)
            ! regularize position
            do j = 1, 3
                reg_x(j) = p%x(j) / box_size(j)
            end do
            write(ifile, '(e0.8)') p%mass
            write(ifile, '(a)') p%elem
            write(ifile, '(*(e0.8, :, " "))') reg_x, p%v
        end do
    end subroutine write_cfg

end module IO
