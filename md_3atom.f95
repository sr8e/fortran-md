program md
    use MD_STRUCT
    use POTENTIALS, only: set_morse_params
    use IO

    implicit none

    !params
    integer step
    integer, parameter :: stepmax = 1000
    !potential
    real, parameter :: eps = 0.2703, alpha = 1.1646, r0 = 3.253
    !mass
    real, parameter :: mass = 26.98

    !total energy
    real, dimension(3) :: toteng  ! pe, ke, temp

    ! i/o
    integer :: itraj = 1, ieng = 2, iload = 10

    dt = 0.001  ! ps
    cutoff = 12.0  ! A

    open(iload, file="models/3atom.atom")
    call read_atom(iload)
    close(iload)
    init_box = box_size

    open(itraj, FILE="trajectory.txt")
    open(ieng, FILE="energy.txt")

    !module initialize
    !simulation box
    is_periodic = [.false., .false., .false.]
    !potentials
    call set_morse_params(eps, alpha, r0)

    ! set init force
    call velo_verlet_step(.true.)

    do step = 0, stepmax
        call velo_verlet_step(.false.)
        write(itraj, '(i5, ",", *(1pe0.8, :, ","))') step, p_set(1)%x, p_set(2)%x, p_set(3)%x
        toteng = total_energy()
        write(ieng, '(i5, ",", *(1pe0.8, :, ","))') step, toteng, sum(toteng(1:2))
    end do

    close(itraj)
    close(ieng)

end program md

