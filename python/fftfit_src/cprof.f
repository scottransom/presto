C  @(#)cprof.f    1.1 9/7/90
      subroutine cprof(y,nmax,nh,c,amp,pha)

C  Compute FFT of profile in array y(nmax), and return amplitude and phase
C  in arrays amp(nh) and pha(nh).  Note that nh=nmax/2, and that the DC term
C  is returned in c(0), fundamental in c(1), ..., Nyquist freq in c(nh).

      parameter(MAXSAM=8192)
      real*4 y(nmax),amp(nh),pha(nh)
      complex c(0:nh),temp(MAXSAM)

      do 10 i=1,nh
 10      temp(i)=cmplx(y(2*i-1),y(2*i))
      call ffft(temp,nmax,1,1)
      c(0)=temp(1)
      do 20 i=1,nh
         c(i)=temp(i+1)
         amp(i)=cabs(c(i))
         pha(i)=0.
 20      if(amp(i).gt.0.) pha(i)=aimag(clog(c(i)))
      return
      end
