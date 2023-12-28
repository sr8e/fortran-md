program md
    use MD_STRUCT
    use POTENTIALS, only: set_morse_params
    use IO

    implicit none

    !params
    integer step
    integer, parameter :: stepmax = 1000
    integer, parameter :: seed = 1891109
    !potential
    real, parameter :: eps = 0.2703, alpha = 1.1646, r0 = 3.253
    !mass
    real, parameter :: mass = 26.98

    !total energy
    real, dimension(3) :: toteng  ! pe, ke, temp

    ! i/o
    integer :: ieng = 2, iload = 10, idump = 11

    dt = 0.001
    cutoff = 12.0

    open(iload, file="models/fcc_al.atom")
    call read_atom(iload)
    close(iload)

    open(ieng, file="energy.txt")

    !module initialize
    !random
    call setup_random(seed)
    !simulation box
    is_periodic = [.true., .true., .true.]
    !potentials
    call set_morse_params(eps, alpha, r0)

    ! set init vel
    call init_velocity_temp(300.0)
    ! set init force
    call velo_verlet_step(.true.)

    do step = 0, stepmax
        call velo_verlet_step(.false.)
        
        toteng = total_energy()
        write(ieng, '(i5, ",", *(1pe0.8, :, ","))') step, toteng, sum(toteng(1:2))
    end do

    close(ieng)

    open(idump, file="al_after.cfg")
    call write_cfg(idump)
    close(idump)

end program md

