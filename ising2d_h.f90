! define J = 1
module extension
    implicit none
    integer,private :: fdigit = 0
contains
    ! called after all outprocessing(any) is completed
    ! value / max : % (value out of max)
    ! isflushed : it is used when you plan to flush. I mainly want to output result at the same time as the progress bar
    subroutine pbout(value, max, isflushed)
        integer,intent(in) :: value, max
        double precision,save :: rate = 0d0, time = 0d0
        double precision dt, dr, estimate
        integer,parameter :: digit = 50
        integer remain(2), values(8), i
        character(len=20) :: FMT
        logical,intent(in) :: isflushed
        if (isflushed) then
            write (*, *)
        end if

        ! check the number of digits of the variable max and create a format specifier
        write (FMT, '(a, i0, a)') "(2(i", int(log10(real(max))) + 1, ",a))"

        ! Get date and time
        call date_and_time(values=values)
        dr = rate   ! store the previously executed rate in dr for the time being
        dt = time   ! store the previously executed time in dt for the time being

        ! Time update : how much time has passed since 9 o clock today, calculations over a day can not be dispalyed.
        time = ((values(5)*60d0+values(6))*60d0+values(7))*1d3+values(8)
        dt = time - dt

        ! percentage update
        rate = dble(value) / dble(max)
        dr = rate - dr

        ! Remaining time calculation (milliseconds)
        estimate = (1d0 - rate) * (dt / dr)

        ! min/sec表記に変換
        remain(1) = int(estimate/6d4) ! minutes
        remain(2) = int((estimate-remain(1)*6d4)*1d-3) ! seconds

        write (*, FMT, advance='no') value, " / ", max, " ["
        do i = 1, int(digit*rate)
            write (*, '(a)', advance='no') "="
        end do
        write (*, '(a)', advance='no') ">"
        do i = int(digit*rate) + 1, digit
            write (*, '(a)', advance='no') "-"
        end do

        if (isflushed) then
            fdigit = 2 * int(log10(real(max))) + 75
        end if
        write (*, '(a, f7.2, a, i3, a, i2, a, a)', advance='no') "]", 100d0*rate, "% ", remain(1), "m", remain(2), "s", char(13)

        if (value == max) then
            write (*, *)
        end if
    end subroutine

    ! ony needed if isfushed us set to true in pbout()
    !call before the first outout processing(Any) starts. Not required when false
    subroutine pbflush()
        character(len=9) FMT
        if (fdigit == 0) then
            return
        end if

        write (FMT, '(a, i0, a)') "(", fdigit, "x, 3a)"
        write (*, FMT, advance='no') char(13), char(8), char(13)
    end subroutine
end module

module subprog
    implicit none
    double precision :: h = 0.8d0
contains
    ! inits : initialize the spin 
    ! init  : 0 or 1
    ! s     : an array that stores the spin direction of each lattice
    ! n     : number of all grid points

    ! init = 0 -> All spins are set to random number
    ! init != 0 -> All spins are set to 0
    subroutine inits(init, s, n)
        integer, intent(in) :: init, n
        double precision, intent(out) :: s(n)
        integer i
        double precision rnd(n)
        
        if(init == 0) then
            call random_number(rnd)
            do i = 1, n
                if(rnd(i) > 0.5d0) then
                    s(i) = 1
                else
                    s(i) = -1
                endif
            enddo
        else
            do i = 1, n
                s(i) = 1
            enddo
        endif
    end subroutine inits
    
    ! ising : update the spin configurstion。
    ! beta  : basically inverse of temperature
    ! s     : array of spin
    ! n     : number of all grid points
    ! hit   : number of reserved spins
    ! ne    : array that stores adjacent coordinate adress of each spin
    subroutine ising(beta, s, n, hit, ne)
        double precision, intent(in) :: beta
        integer, intent(in) :: n, ne(n,4)
        double precision, intent(inout) :: s(n)
        double precision, intent(out) :: hit
        integer i, m, nhit
        double precision rnd1(n), rnd2(n), ds
        
        nhit = 0
        ! set random numbers in array rand1 and rand 2
        call random_number(rnd1)
        call random_number(rnd2)

        do i = 1, n
            ! selet random grid points ad convert to straight numbers
            m = int(dble(n)*rnd1(i)) + 1
            ! in rare cases the adress may exceed the range so in that case set it to m
            if(m > n) then
                m = n
            endif
            
            !calculate the energy if a random grid point is flipped
            ds = 2*s(m)*(beta*(s(ne(m,1))+s(ne(m,2))+s(ne(m,3))+s(ne(m,4)))+h) 
            
            ! if the generated random number value is smaller than the enrgy calculated above, perform spin flip
            if(rnd2(i) < exp(-dble(ds))) then
                ! count the number of flipped spins
                nhit = nhit + 1
                ! invert the spin
                s(m) = -1.0d0*s(m)
            endif
        enddo
        !find the ratio of number of flipped spin to the total
        hit = dble(nhit) / dble(n)
    end subroutine ising

    ! magnet    : calculate magnetix=zation
    ! mgn       : array that store spin
    ! s         : direction of lattice spin
    ! n         : total numbe of spins
    subroutine magnet(mgn, s, n)
        integer, intent(in) :: n
        double precision, intent(in) :: s(n)
        double precision, intent(out) :: mgn
        integer i, m
        
        ! examine the global magnetization by suming the global spin configuration
        m=0
        do i=1,n
            m=m+int(s(i))
        enddo
        ! examie the ratio of total magnetization .i.e. magnetization per spin
        mgn=dble(m)/dble(n)
    end subroutine magnet

    subroutine hamiltonian(energy, s, n, ne,beta)
        integer, intent(in) :: n, ne(n,4)
        double precision, intent(in) :: beta
        double precision, intent(in) :: s(n)
        double precision, intent(out) :: energy
        double precision delta
        integer i
        energy = 0.0d0
        delta = 0.0d0
        do i = 1, n
            delta = - 0.5*beta*s(i)*(s(ne(i,1))+s(ne(i,2))+s(ne(i,3))+s(ne(i,4)))-h*s(i)
            energy = energy + delta
        end do
        energy = energy / dble(n)
    end subroutine hamiltonian
      
    ! List vector： Compute adresses of adjacent grid points
    ! ne(m,direction) 
    ! m : street number (any where between 1 and L**2)
    ! direction : 1->right, 2-> Up, 3->Left, 4-> down
    ! n1 : Length of One turn in x direction
    ! n2 : Length of turn in the y direction
    subroutine list_vector(ne, n1, n2)
        integer, intent(in) :: n1, n2
        integer, intent(out) :: ne((n1*n2),4)
        integer m1, m2, i1, i2

        ! i2=y:1 ~ n2
        do i2 = 1, n2
            ! m2:(y-1)*n1 = (y-1)L
            ! m2 is pivot of y
            m2 = (i2-1)*n1
            ! i1=x:1 ~ n1
            do i1 = 1, n1
               
                m1 = m2 + i1
               
                ne(m1,1) = ( m2 + mod(i1,n1) ) + 1

                
                ne(m1,2) = mod(i2,n2)*n1 + i1

              
                ne(m1,3) = m2 + mod(i1-2+n1,n1) + 1

                
                ne(m1,4) = mod(i2-2+n2,n2)*n1 + i1
            enddo
        enddo
    end subroutine list_vector
end module subprog

program ising2d
    use extension
    use subprog
    implicit none
    double precision, allocatable :: s(:)
    integer, allocatable :: ne(:,:)
    double precision :: beta, hit, mgn, thit, tmgn, tmgn2, energy, t_energy1, t_energy2,tmgn4,bc
    double precision :: specific_heat, suscepibility
    integer :: n1=96, n2=96, nbeta=100
    integer :: init=0, therm=1000, inter=10, niter=50000
    integer :: n, i, j, nb
    open(10, file='mc_96_0.8.txt')
    n=n1*n2
    allocate (s(n), ne(n,4))
    
    call random_seed
    
    call list_vector(ne, n1, n2)
    
    
    write(*,200) n1, n2, init, therm, inter, niter, h

200 format('size:',2I3,' hot(1)/cold(0):',I2, &
    ' Therm:',I5,' Int:',I5,' NConf:',I5,' h:', F6.2)
    
    
    do nb=0,nbeta
        ! beta = 1 / T
        beta=0.01d0+nb*0.0079d0
        
        
        call inits(init, s, n)
        ! thit = Total hits
        ! tmgn = Total magnetization
        thit=0.0
        tmgn=0.0
        tmgn2 = 0.0d0
        t_energy1 = 0.0d0
        t_energy2 = 0.0d0
        tmgn4 = 0.0d0
        bc = 0
   
        do i=1,therm
            call ising(beta, s, n, hit, ne)
        enddo
        
        do i=1,niter
            
            do j=1,inter
                call ising(beta, s, n, hit, ne)
              
                thit=thit+hit
            enddo
            
            call magnet(mgn, s, n)
          
            tmgn = tmgn + abs(mgn)
            tmgn2 = tmgn2 + mgn*mgn
            tmgn4 = tmgn4 + tmgn2*tmgn2
            call hamiltonian(energy, s, n, ne,beta)
            t_energy1 = t_energy1 + energy
            t_energy2 = t_energy2 + energy * energy
        enddo
        !
        tmgn=tmgn/dble(niter)
        tmgn2 = tmgn2 / dble(niter)
        tmgn4 = tmgn4/dble(niter)
        thit=thit/dble(inter*niter)
        t_energy1 = t_energy1 / dble(niter)
        t_energy2 = t_energy2 / dble(niter)
        specific_heat = n * (t_energy2 - t_energy1*t_energy1)*beta
        suscepibility = n * (tmgn2 - tmgn*tmgn)
        bc = 1- (tmgn4/3*(tmgn2)**2)

        call pbflush()
        write(*,300) beta, tmgn, thit, t_energy1, specific_heat, suscepibility
        call pbout(nb, nbeta, .true.)

        write(10, *) beta, tmgn, specific_heat, suscepibility
    enddo
    close(10)
300 format('beta=',f8.5,' mgn=',f9.6,' hit=',f9.5,' energy=',f0.5,' specific heat=', f0.5,' suscepibility=', f0.5)
end program ising2d
