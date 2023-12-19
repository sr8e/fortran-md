module POTENTIALS
    implicit none

    integer pot_id
    ! morse
    real eps, alpha, r0

    contains

    subroutine set_morse_params(e, a, r)
        real, intent(in) :: e, a, r
        ! set potential type
        pot_id = 1

        eps = e
        alpha = a
        r0 = r
    end subroutine set_morse_params

    function f(r)
        real f, r
        select case(pot_id)
            case(1)
                f = f_morse(r)
            case default
                f = 0
        end select
    end function f

    function p(r)
        real p, r
        select case (pot_id)
            case(1)
                p = p_morse(r)
            case default
                p = 0
        end select
    end function p

    function f_morse(r)
        ! unit: eV/A
        real r, f_morse
        f_morse = -2.0 * alpha * eps * ( &
            exp(-2.0 * alpha * (r - r0)) - exp(-alpha * (r - r0)) &
        )
    end function f_morse

    function p_morse(r)
        ! unit: eV
        real r, p_morse
        p_morse = eps * (exp(-2.0 * alpha * (r - r0)) - 2.0 * exp(-alpha * (r - r0)))
    end function p_morse
end module POTENTIALS
