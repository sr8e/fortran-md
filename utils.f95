module UTILS
    implicit none

    real, parameter :: NA = 6.02214076e+23, e = 1.6021892e-19, pi = 3.1415926536

    contains

    function convert_unit(val, to_eV)
        ! convert eV <-> (g/mol)(A/ps)^2
        ! 1 eV = e [J] = e [kg (m/s)^2] = 10^-1 e NA [(g/mol)(A/ps)^2])
        logical to_eV
        real convert_unit, val
        if (to_eV) then
            convert_unit = val * 10.0 / (e * NA)
        else
            convert_unit = val * (e * NA) / 10.0
        end if
    end function convert_unit

    function dot(v1, v2)
        ! calculate inner product
        real, dimension(3) :: v1, v2
        real dot
        dot = sum(v1 * v2)
    end function dot

    function d_euc(d_vec)
        ! calculate euclid norm of vector
        real, dimension(3) :: d_vec
        real d_euc
        d_euc = sqrt(dot(d_vec, d_vec))
    end function d_euc

    function acc(mass, force)
        ! calculate acceleration of mass by applied force
        ! unit: mass -> g/mol, force -> eV/A, acc -> A/ps^-2

        real mass, force, acc
        acc = convert_unit(force / mass, .false.)
    end function acc

    function periodic(p, q)
        real periodic, p, q
        periodic = mod(p, q)
        if (p .lt. 0) then
            periodic = periodic + q
        end if
    end function periodic

    ! random
    subroutine setup_random(seed)
        integer, intent(in) :: seed
        integer n
        integer, dimension(:), allocatable :: seed_array

        ! get seed array size (compiler dependent)
        call random_seed(size=n)

        ! allocate seed array and put specified value
        allocate(seed_array(n))
        seed_array = seed
        call random_seed(put=seed_array)
    end subroutine setup_random

    function normdist(mu, var)
        implicit none
        ! generate random number which follows N(mu, var)
        ! by Box-Muller's method
        real mu, var
        real normdist, u1, u2
        call random_number(u1)
        call random_number(u2)
        normdist = sqrt(var) * sqrt(-2.0*log(1-u1))*cos(2*pi*u2) + mu
    end function

end module UTILS
