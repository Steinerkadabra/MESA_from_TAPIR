program hello
    real(8) :: MESA_Teff, MESA_end_time
    real(kind=8), allocatable, dimension(:) :: time, rate
    integer :: nlines, i, ierr



    nlines = -1
    open(unit = 10, file = "accr_rate_noburst.dat")
    do 
        read(10,*,iostat=io)
        if (io/=0) exit
            nlines = nlines + 1 
    end do 
    close(10)


    open(unit = 10, file = "accr_rate_noburst.dat")

    allocate(time(nlines))
    allocate(rate(nlines))

    read(10,*)
    read(10,*)
    do i = 1, nlines -2
         read(10,*) time(i),  rate(i)
    end do
    close(10)


    do i = 1, nlines
         time(i) = time(i)/3.1558149984d7
         rate(i) = rate(i)/1.9892d33*3.1558149984d7
    end do

    do i = 1, nlines-2
	ierr = 0
	
	if (time(i) > 205195.87628815800) then
	 
        write(*,*) 'start_time', time(i), 'end_time', time(i+1), 'start_rate', rate(i), 'end_rate', rate(i+1)
        call start_MESA(rate(i), rate(i+1) , time(i), time(i+1), MESA_Teff, MESA_end_time, ierr)
        write(*,*) 'resulting MESA Teff:', MESA_Teff
        end if
	if (ierr == 1) STOP "MESA has not reached end time"
    end do
end program


subroutine start_MESA(start_rate, end_rate, start_time, end_time, MESA_Teff, MESA_end_time, ierr)
    real(8), intent(in) :: start_rate, end_rate, start_time, end_time
    real(8), intent(out) :: MESA_Teff, MESA_end_time
    integer, intent(out) :: ierr
    character(len = 300) :: run_string
    character(len=26) :: start_rate_str, end_rate_str, start_time_str, end_time_str
    character(len=10) :: rn_tapir
    character(len=12) :: output_file

    ierr = 0
    rn_tapir = './rn_tapir'
    output_file = '> output.txt'

    write(start_rate_str , *) start_rate
    write(end_rate_str , *) end_rate
    write(start_time_str , *) start_time
    write(end_time_str , *) end_time

    run_string = rn_tapir//start_rate_str//end_rate_str//start_time_str//end_time_str//output_file

    call execute_command_line(run_string, wait = .true.)


    open (2, file = 'MESA_for_TAPIR.txt', status = 'old')
    read(2,*) MESA_Teff, MESA_end_time
    close(2)

    if (abs(MESA_end_time - end_time) > 1d-10) then
        ierr = 1
    end if

 end subroutine start_MESA
