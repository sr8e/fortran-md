
program lattice
    open(1, FILE="fcc_al.atom")
    call create_lattice(1, "fcc", 26.98, "Al", 4.05, [3, 3, 3])
    close(1)
end program lattice

subroutine create_lattice(ifile, lat, mass, elem, a, duplicate)
    character(3), intent(in) :: lat
    character(2), intent(in) :: elem
    real, intent(in) :: mass, a
    integer, dimension(3), intent(in) :: duplicate
    integer n, i, j, k, l
    real, dimension(3) :: box_size
    real, dimension(:, :), allocatable :: unit

    select case(lat)
        case("fcc")
            n = 4 * product(duplicate)
            box_size = a * duplicate
            write(ifile, '(i0)') n
            write(ifile, '(*(1pe0.8, :, " "))') box_size
            ! unit cell
            allocate(unit(4, 3))
            unit(1, :) = [0.0, 0.0, 0.0]
            unit(2, :) = [0.5, 0.5, 0.0]
            unit(3, :) = [0.0, 0.5, 0.5]
            unit(4, :) = [0.5, 0.0, 0.5]

            do i = 0, duplicate(1) - 1
                do j = 0, duplicate(2) - 1
                    do k = 0, duplicate(3) - 1
                        do l = 1, 4
                            write(ifile, '(1pe0.8, " ", a, " ", *(1pe0.8, :, " "))') &
                                mass, elem, a*([i, j, k] + unit(l, :)), 0.0, 0.0, 0.0
                        end do
                    end do
                end do
            end do
    end select
end subroutine create_lattice
