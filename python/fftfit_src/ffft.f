C @(#)ffft.f    1.1 9/7/90
      subroutine ffft(d,npts,isign,ireal)

C  Fourier transform of length npts=2**k, performed in place.
C  Input data in array d, treated as complex if ireal=0, and as real if ireal=1.
C  In either case the transform values are returned in array d, treated as
C  complex. The DC term is d(1), and d(npts/2+1) is the term at the Nyquist
C  frequency.  The basic algorithm is the same as Norm Brenner's FOUR1, and
C  uses radix-2 transforms.

C  J. H. Taylor, Princeton University.

      parameter(MAXSAM=8192)
      complex d(MAXSAM),t,w,wstep,tt,uu
      data pi/3.14159265/

C  Shuffle the data to bit-reversed order.

      imax=npts/(ireal+1)
      irev=1
      do 5 i=1,imax
         if(i.ge.irev) go to 2
         t=d(i)
         d(i)=d(irev)
         d(irev)=t
 2       mmax=imax/2
 3       if(irev.le.mmax) go to 5
         irev=irev-mmax
         mmax=mmax/2
         if(mmax.ge.1) go to 3
 5       irev=irev+mmax

C  The radix-2 transform begins here.

         api=isign*pi/2.
         mmax=1
 6       istep=2*mmax
         wstep=cmplx(-2.*sin(api/mmax)**2,sin(2.*api/mmax))
         w=1.
         do 9 m=1,mmax

C  This in the inner-most loop -- optimization here is important!
            do 8 i=m,imax,istep
               t=w*d(i+mmax)
               d(i+mmax)=d(i)-t
 8             d(i)=d(i)+t

 9             w=w*(1.+wstep)
         mmax=istep
         if(mmax.lt.imax) go to 6
         if(ireal.eq.0) return

C  Now complete the last stage of a doubled-up real transform.

         jmax=imax/2 + 1
         wstep=cmplx(-2.*sin(isign*pi/npts)**2,sin(isign*pi/imax))
         w=1.0
         d(imax+1)=d(1)

         do 10 j=1,jmax
            uu=cmplx(real(d(j))+real(d(2+imax-j)),aimag(d(j)) -
     +           aimag(d(2+imax-j)))
            tt=w*cmplx(aimag(d(j))+aimag(d(2+imax-j)),-real(d(j)) +
     +           real(d(2+imax-j)))
            d(j)=uu+tt
            d(2+imax-j)=conjg(uu-tt)
 10         w=w*(1.+wstep)

      return
      end
