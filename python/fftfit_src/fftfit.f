      subroutine fftfit(prof,s,phi,nmax,shift,eshift,snr,esnr,b,errb
     $     ,ngood)

C  Fourier transform domain routine for determining pulse TOAs.
C  Input data:
C    prof(nmax)    profile
C    s(nh)        standard profile amplitude
C    phi(nh)         standard profile phase
C    nmax        length of prof

C  Outputs:
C    shift        shift required to align std with prof, in bins
C    eshift            uncertainty in shift
C    snr        signal to noise ratio of prof
C    esnr        uncertainty in snr

C  Method:
C  It is assumed that prof(j)=a + b*std(j-tau) + noise(j), where the
C  (j-tau) subscript is supposed to mean that the time-shift between the
C  observed and standard profiles need not be an integer number of bins.
C  The algorithm is a straightforward Chi-squared minimization with respect
C  to a, b, and tau.  (The result for a, being of no physical or instrumental
C  interest, is not actually evaluated -- though it easily could be.)
C  First and second partial derivatives of Chisqr with respect to b and tau
C  are computed, and used to evaluate s, snr, and their "standard errors,"
C  or one-sigma uncertainties.  The only special trick is that the expression
C  for the best-fit value of tau defines it implicitly, and cannot be solved
C  analytically.  It is solved numerically instead, finding the minimum of
C  Chisqr near a best guess from a CCF at 32 lags done in the Fourier domain.

C  Also note that it may
C  be desirable to return the scale factor b relating prof to std, instead
C  of snr.  In that case you could also return the noise estimate rms.

      parameter (twopi=6.2831853,MAXSAM=8192)
      real*4 prof(MAXSAM),p(MAXSAM/2),theta(MAXSAM/2)
      real*4 s(MAXSAM/2),phi(MAXSAM/2),r(MAXSAM/2),tmp(MAXSAM/2)
      complex cp(0:MAXSAM/2)
      logical low,high

C        if (phi(1).ne.0) then
C          print *,' Phase of fundamental not zero, check .hm file'
C          stop
C        end if

      nh=nmax/2
      sum=0.
      do 1 i=nh/2+1,nh
         sum=sum+s(i)
 1    continue
      ave=2*sum/nh

      do 2 i=1,nh
         if(s(i).lt.ave) go to 3
 2    continue
 3    ngood=i-1

      call cprof(prof,nmax,nh,cp,p,theta)
      do 10 k=1,nh
         tmp(k)=p(k)*s(k)
C     write(13,1010) k,s(k),phi(k),p(k),theta(k)
 1010    format(i5,4f10.5)
 10      r(k)=theta(k)-phi(k)

      fac=nmax/twopi

      call fccf(tmp,r,shift)

C  The "DO 60" loop solves the transcendental equation yielding the best-fit
C  value of tau.

      tau=shift
      nsum0=min(16,ngood/4)
      do 60 nsum=nsum0,ngood
         dtau=0.2/nsum
         edtau=0.01/nsum
         if (nsum.gt.(nh/2.+.5)) edtau=1.e-4

         ntries=0
         low=.false.
         high=.false.
 50      ftau=dchisqr(tau,tmp,r,nsum)
         ntries=ntries+1
         if(ftau.lt.0.0) then
            a=tau
            fa=ftau
            tau=tau+dtau
            low=.true.
         else
            b=tau
            fb=ftau
            tau=tau-dtau
            high=.true.
         end if
         if (ntries.gt.100) then
            shift=0.
            eshift=999.
            snr=0.
            esnr=0.
            return
         end if
         if (low.neqv.high) go to 50
         tau=zbrent(a,b,fa,fb,edtau,tmp,r,nsum)
 60   continue

      s1=0.
      s2=0.
      s3=0.
      do 80 k=1,ngood
         cosfac=cos(-r(k)+k*tau)
         s1=s1 + tmp(k)*cosfac
         s2=s2 + s(k)**2
 80      s3=s3 + k**2 *tmp(k)*cosfac
      b=s1/s2
      s1=0.
      do 90 k=1,ngood
         sq=p(k)**2-2.*b*p(k)*s(k)*cos(r(k)-k*tau)+(b*s(k))**2
 90      s1=s1+sq
      rms=sqrt(s1/ngood)
      errb=rms/sqrt(2.0*s2)
      errtau=0.0                ! ###
      if(s3.gt.0.0) errtau=rms/sqrt(2.0*b*s3)
      snr=2.0*sqrt(2.0*nh)*b/rms
      shift=fac*tau
      eshift=fac*errtau
      esnr=snr*errb/b

      return
      end
