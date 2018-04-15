!
! File:   sim_ode.f
! Authors: Xubin Li & Anil Korkut @ MDACC
! Created on March 27, 2017, 5:05 AM
!

!
! computer simulation of network dynamics (steady state) from perturbation
! Hopfield network models and perturbation
! model X[i]'=eps[i]*tanh(Sigma(W[i,j]*X[j])+U[i])-apha[i]*X[i]
!

module sim_ode
    implicit none
    double precision, allocatable, save :: u(:), w(:), alpha(:), eps(:)

contains

    subroutine conv_int_str(num, str)
        implicit none
        integer num
        character*(*) :: str

        integer ndigits, tmp
        character(len = 10) :: format_string

        if (num < 0) then
            write(*, *) "ERROR: negative integer can not be converted."
            stop
        end if
        ndigits = 1
        tmp = num / 10
        do while (tmp /= 0)
            tmp = tmp / 10
            ndigits = ndigits + 1
        end do
        if (ndigits > 10) then
            write(*, *) "ERROR: the integer is too large ( >= 10^10 ) ."
            stop
        end if

        format_string(:) = ' '
        format_string = '(I' // char(48 + ndigits) // ')'

        str(:) = ' '
        write (str, format_string) num
        str = trim(adjustl(str))

        return
    end subroutine

    subroutine str_split(str, word)
        implicit none
        character*(*) :: str
        character*(*) :: word(:)
        integer :: pos1, pos2, pos3, n
 
        pos1 = 1
        pos2 = 0
        pos3 = 0
        n = 0
        word(:) = ""
        
        do
            pos2 = index(str(pos1:), " ")
            pos3 = index(str(pos1:), char(9))
            if (pos3 /= 0) then
                if (pos2 == 0) then
                    pos2 = pos3
                else 
                    if (pos2 > pos3) pos2 = pos3
                end if
            end if
          
            if (pos2 == 0) then
                n = n + 1
                word(n) = str(pos1:)
                exit
            else if (pos2 /= 1) then
                n = n + 1
                word(n) = str(pos1:(pos1 + pos2 - 2))
            end if
            pos1 = pos2 + pos1
        end do
    end subroutine

    ! ode solver interface
    subroutine sim(neq, x, dlsode_istate)
        implicit none
        external fex, jex
        integer neq
        double precision x(*)
        integer dlsode_istate

        integer iopt, itask, itol, mf, liw, lrw
        double precision atol, rtol, t, tout
        integer, allocatable :: iwork(:)
        double precision, allocatable :: rwork(:)

        itol = 1
        rtol = 1.0D-20
        atol = 1.0D-10

        t = 0.0D0
        tout = 1.0D-10

        itask = 1
        dlsode_istate = 1
        iopt = 0

        mf = 10
        liw = 20
        lrw = 20 + 16 * neq
        allocate(iwork(liw))
        allocate(rwork(lrw))

        do while (.true.)
            call dlsode(fex, neq, x, t, tout, itol, rtol, atol, itask, &
            dlsode_istate, iopt, rwork, lrw, iwork, liw, jex, mf)
            tout = tout * 2.0
            if (dlsode_istate > 0) then
                cycle
            else
                exit
            end if
        end do

        if (dlsode_istate /= -1) then
            write(*, *) "Warning: dlsode may not converge due to dlsode_istate = ", dlsode_istate
        end if

        deallocate(iwork)
        deallocate(rwork)
        return
    end subroutine

end module

! ode of netwrok dynamics
subroutine fex(neq, t, x, xdot)
    use sim_ode
    implicit none
    integer neq
    double precision t, x(*), xdot(*)
    
    integer i

    do i = 1, neq
        xdot(i) = eps(i) &
        *tanh(sum(w(((i - 1) * neq + 1):((i - 1) * neq + neq)) * x(1:neq)) + u(i)) &
        -alpha(i) * x(i)
    end do

    return
end subroutine

! dummy subroutine; copied form the example in dlsode.f
subroutine jex(neq, t, y, ml, mu, pd, nrpd)
    double precision pd, t, y
    dimension y(3), pd(nrpd, 3)
    pd(1, 1) = -.04d0
    pd(1, 2) = 1.d4 * y(3)
    pd(1, 3) = 1.d4 * y(2)
    pd(2, 1) = .04d0
    pd(2, 3) = -pd(1, 3)
    pd(3, 2) = 6.d7 * y(2)
    pd(2, 2) = -pd(1, 2) - pd(3, 2)
    return
end

program main
    use sim_ode
    implicit none

    integer n_nodes
    double precision, allocatable :: x(:)
    integer dlsode_istate

    character(len = 10000) :: line
    character(len = 100), dimension(100) :: buf_str
    integer io_status

    character(len = 100) :: file_input, file_pert, file_out, model_file_prefix, model_file_suffix
    integer model_id_beg, model_id_end
    integer, allocatable :: model_ids(:)
    character(len = 100) :: model_id_str
    integer unit_id, unit_pert, unit_out, unit_model
    character(len = 100000) :: my_fmt

    integer n_models, n_nodes_pert, i, j, k, m, n
    integer, allocatable :: pert_index(:)
    character(len = 100) :: pert_name
    
    integer :: num_args
    character(len = 50), dimension(:), allocatable :: args
    
    ! program discription
    write(*, *) "Predicting steady states of network nodes under in silico perturbations"
    write(*, *) "Authors: Xubin Li & Anil Korkut @ MDACC"
    write(*, *) "Version: 20180327"
    write(*, *) ""

    ! parse command line
    num_args = command_argument_count()
    allocate(args(num_args))

    do i = 1, num_args
        call get_command_argument(i, args(i))
        args(i) = trim(adjustl(args(i)))
    end do

    if (num_args == 0) then
        write(*, *) "DESCRIPTION"
        write(*, *) "-i input", char(9), "Run your prediction"
        write(*, *) "-h", char(9), char(9), "Help"
        write(*, *) ""
        stop
    end if
    n = 1
    do while (.true.)
        if (n > num_args) exit
        if (args(n) == "-i") then
            file_input = args(n + 1)
            n = n + 1
        else if (args(n) == "-h") then
            write(*, *) "An example of input:"
            write(*, *) "n_nodes", char(9), "99"
            write(*, *) "file_pert", char(9), "pert_insilico.txt"
            write(*, *) "model_file_prefix", char(9), "model_"
            write(*, *) "model_file_suffix", char(9), ".txt"
            write(*, *) "model_id_beg", char(9), "1"
            write(*, *) "model_id_end", char(9), "80"
            write(*, *) ""
            stop
        end if
        n = n + 1
    end do
    deallocate(args)

    ! parse input file
    unit_id = 51
    open(unit_id, file = file_input)
    buf_str(:) = ""
    do while (.true.)
        read(unit_id, '(A)', iostat = io_status) line
        if (io_status > 0) then
            write(*, *) "Error: check ", file_input, " something was wrong"
            stop
        else if (io_status < 0) then
            exit
        end if

        line = trim(adjustl(line))
        if (line(1:1) == "#" .or. line(1:1) == "") then
            cycle
        else
            call str_split(line, buf_str)
            if (buf_str(1) == "n_nodes") then
                read(buf_str(2), *) n_nodes
            else if (buf_str(1) == "file_pert") then
                read(buf_str(2), *) file_pert
            else if (buf_str(1) == "model_file_prefix") then
                read(buf_str(2), *) model_file_prefix
            else if (buf_str(1) == "model_file_suffix") then
                read(buf_str(2), *) model_file_suffix
            else if (buf_str(1) == "model_id_beg") then
                read(buf_str(2), *) model_id_beg
            else if (buf_str(1) == "model_id_end") then
                read(buf_str(2), *) model_id_end
            end if
        end if
    end do
    close(unit_id)
    ! print
    write(*, *) "Parameters and file names from input:"
    write(*, *) "n_nodes", char(9), n_nodes
    write(*, *) "file_pert", char(9), file_pert
    write(*, *) "model_file_prefix", char(9), model_file_prefix
    write(*, *) "model_file_suffix", char(9), model_file_suffix
    write(*, *) "model_id_beg", char(9), model_id_beg
    write(*, *) "model_id_end", char(9), model_id_end
    write(*, *) ""
    
    ! note of results
    write(*, *) "Note: the prediction results will be written in files with prefix 'predict_'"
    write(*, *) "The 1st column is model ID"
    write(*, *) "The 2nd column is istate from calling DLSODE (where istate = -1 indicates convergence)"
    write(*, *) "The remaining columns are final values of x's"
    write(*, *) ""

    ! model ids
    n_models = model_id_end - model_id_beg + 1
    allocate(model_ids(n_models))
    do i = 1, n_models
        model_ids(i) = model_id_beg + i - 1
    end do

    ! allocate w, x, u
    allocate(w(n_nodes * n_nodes))
    w(:) = 0.0
    allocate(x(n_nodes))
    allocate(u(n_nodes))
    x(:) = 0.0
    u(:) = 0.0

    ! alpha and eps
    allocate(alpha(n_nodes))
    allocate(eps(n_nodes))
    alpha(:) = 1.0
    eps(:) = 1.0

    ! process each pert condictions and predict response
    unit_pert = 51
    open(unit_pert, file = file_pert)
    do while (.true.)
        read(unit_pert, '(A)', iostat = io_status) line
        if (io_status > 0) then
            write(*, *) "Error: check ", file_pert, " something was wrong"
            stop
        else if (io_status < 0) then
            exit
        end if

        line = trim(adjustl(line))
        if (line(1:1) == "#" .or. line(1:1) == "") cycle

        ! parse line
        call str_split(line, buf_str)
        
        ! assign
        read(buf_str(1), *) n_nodes_pert
        allocate(pert_index(n_nodes_pert))
        u(:) = 0.0
        
        read(buf_str(2:(2*n_nodes_pert + 2)), *) &
        (pert_index(k), k = 1, n_nodes_pert), &
        (u(pert_index(k)), k = 1, n_nodes_pert), &
        pert_name
        
        ! convert inhibited percentage to a log2 unit
        do i = 1, n_nodes
            u(i) = log10(1.0 - u(i)) * 3.322
        end do
        pert_name = trim(adjustl(pert_name))

        ! open file for out
        unit_out = 91
        open(unit = unit_out, file = "predict_" // trim(pert_name) // ".txt")
    
        ! loop over each model
        do m = 1, n_models
            write(*, *) "Predict ", trim(pert_name), char(9), "Model ", model_ids(m)

            ! read
            unit_model = 52
            call conv_int_str(model_ids(m), model_id_str)
            open(unit_model, file = trim(model_file_prefix) // trim(model_id_str) // trim(model_file_suffix))
            do i = 1, n_nodes
                read(unit_model, *) (w((i - 1) * n_nodes + j), j = 1, n_nodes)
            end do
            do i = 1, n_nodes
                read(unit_model, *) alpha(i)
            end do
            do i = 1, n_nodes
                read(unit_model, *) eps(i)
            end do
            close(unit_model)

            ! calculation
            x(:) = 0.0
            call sim(n_nodes, x, dlsode_istate)
            write(my_fmt, '( "(I10, 1X, I10, 1X,", I10, "(F9.3, 1X))" )') n_nodes
            write(unit_out, my_fmt) model_ids(m), dlsode_istate, (x(i), i = 1, n_nodes)

            call flush()
        end do

        ! deallocate
        deallocate(pert_index)
        close(unit_out)
        
        write(*, *) ""
    end do

    ! deallocate
    deallocate(w)
    deallocate(x)
    deallocate(u)
    deallocate(alpha)
    deallocate(eps)
    deallocate(model_ids)
    
    write(*, *) "Done!"
end program
