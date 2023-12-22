module UTILS
    implicit none

    real, parameter :: NA = 6.02214076e+23, e = 1.6021892e-19

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

end module UTILS
