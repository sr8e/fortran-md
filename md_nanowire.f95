program md
    use MD_STRUCT
    use POTENTIALS, only: set_morse_params
    use IO

    implicit none

    !params
    integer step
    integer, parameter :: stepmax = 2000
    !potential
    real, parameter :: eps = 0.2703, alpha = 1.1646, r0 = 3.253
    !mass
    real, parameter :: mass = 26.98

    !total energy
    real, dimension(3) :: toteng  ! pe, ke, temp
    real :: prev_eng = 0

    !relax
    real :: erate = 0.01
    integer :: dstep = 0, dstepmax = 100  ! final nominal strain = erate * dstepmax

    ! i/o
    integer :: ieng = 2, istrain = 3, iload = 10, idump = 11
    character(4) fn_step  ! for output file

    dt = 0.001
    cutoff = 12.0

    open(iload, file="models/fcc_al.atom")
    call read_atom(iload)
    close(iload)
    init_box = box_size

    open(ieng, FILE="energy.txt")
    open(istrain, file="strain.txt")

    !module initialize
    !simulation box
    is_periodic = [.false., .false., .true.]
    !potentials
    call set_morse_params(eps, alpha, r0)

    ! set init force
    call velo_verlet_step(.true.)
    do step = 0, stepmax
        call velo_verlet_step(.false.)
        
        ! relax
        call relax_quickmin()
        toteng = total_energy()
        write(ieng, '(i5, ",", *(1pe0.8, :, ","))') step, toteng, sum(toteng(1:2))

        if(abs(sum(toteng(1:2))-prev_eng) .lt. 1.0e-4) then
            write(istrain, '(i0, ",", *(1pe0.8, :, ","))') step, dstep * erate, sum(toteng(1:2))
            call deform(erate, dstep)
            dstep = dstep + 1
            if (dstep .eq. dstepmax) then
                exit
            end if
        end if
        prev_eng = sum(toteng(1:2))
        if(mod(step, 100) .eq. 0) then
            print *, step, dstep
            write(fn_step, '(i0.4)') dstep
            open(idump, file="al_deform_"//fn_step//".cfg")
            call write_cfg(idump)
            close(idump)
        end if
    end do

    if (dstep .ne. dstepmax) then
        print *, "deformation has not completed!"
    end if

    close(ieng)
    close(istrain)



end program md

