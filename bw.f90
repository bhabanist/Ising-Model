program ising_model

  implicit none
  
  integer :: i, ic, nm
  real(kind=8) :: h, bmin, xa, xb, b, xx, fx, df, zro,chi,deno
  integer, parameter :: q = 4
  !character(len=20), parameter :: filename = "bw_0.1.txt"
  
  h = 0.1d0
  bmin = 0.01d0
  nm = 100
  
  ! cordination number q
  zro = 1.0d-10
  
  ! Open output file
  open(unit=1, file="mbw_0.1.txt")

  do i = 0, nm
    xa = 0.0d0
    xb = 1.0d0
    b = bmin + i * 0.0079d0
    ic = 0
    
123 continue
    ic = ic + 1
    xx = (xa + xb) / 2.0d0
    fx = tanh(q * b * xx + h)
    df = fx - xx
    if (abs(df) <= zro) then
      goto 234
    else if (df < 0.0d0) then
      xb = xx
      goto 123
    else
      xa = xx
      goto 123
    endif
    
234 continue
deno = (COSH(h + q*b*xx))**2 - 4*b
chi = 1/deno
    write(1, *) b, xx,chi
  end do
  
  ! Close output file
  close(unit=1)

end program ising_model

