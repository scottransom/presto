	subroutine rdascat(lu,npsr)
	include '../include/rdascat.h'

	open(lu,file='../lib/ascat1.dat',status='old')

	do 100 i=1,NP
	read(lu,1001,end=999) jname(i),bname(i),ra2000(i),ra1950(i),rae(i)
1001	format(a12,a8,2f16.12,d10.2)
	read(lu,1002) dec2000(i),dec1950(i),dece(i),n5,n4,n3,n2,n1,n0,
     +    dmin(i),dmax(i),dist(i),ndflag(i),lcode(i),ucode(i)
1002	format(2f16.12,d10.2,6i1,3x,3f8.2,i2,1x,2a1)
	nscode(i)=n5*32768+n4*4096+n3*512+n2*64+n1*8+n0
	read(lu,1003) ldeg(i),bdeg(i),pmra(i),pmrae(i),pmdec(i),
     +    pmdece(i),posepoch(i)
1003	format(2f10.4,2(f10.4,f8.4),f10.2)
	read(lu,1004) p(i),pe(i),pdot(i),pdote(i),epoch(i)
1004	format(f22.19,d10.2,d16.8,d10.2,f14.6)
	read(lu,10042) f2(i),f2e(i),f3(i),f3e(i),tau(i),ntauflag(i),
     +    t408(i),n2,n1,n0
10042	format(d14.6,d10.2,d14.6,d10.2,f6.2,i2,f8.1,3x,3i1)
	ntype(i)=n2*64 + n1*8 + n0
	read(lu,1005) dm(i),dme(i),rm(i),rme(i),we(i),w50(i),w10(i)
1005	format(f14.5,f12.5,2f12.3,3f10.3)
	read(lu,1006) s400(i),s600(i),s1400(i),modcode(i),
     +    limcode(i),distmod(i),lum(i),
     +    bsurf(i),age(i),edot(i),ibin(i)
1006	format(3f8.2,2i2,f8.4,4f10.3,3x,i1)
	if(ibin(i).ge.1) then
	  read(lu,1007) pb(i),pbe(i),a1(i),a1e(i),om(i),ome(i)
1007	  format(f18.12,d10.2,f15.8,d10.2,f12.6,f10.6)
	  read(lu,1008) omdot(i),omdote(i),e(i),ee(i),t0(i),t0e(i)
1008	  format(f12.6,f10.6,f12.9,d10.2,f18.10,d10.2)
          read(lu,1009) gamma(i),gammae(i),pbdot(i),pbdote(i),
     +      si(i),sie(i),r(i),re(i)
1009	  format(2(d12.4,d10.2),f8.4,f6.2,d12.4,d10.2)
	endif
	if(ibin(i).ge.2) then
	  read(lu,1007) pb2(i),pb2e(i),a12(i),a12e(i),om2(i),om2e(i)
	  read(lu,1008) omdot2(i),omdot2e(i),e2(i),e2e(i),t02(i),t02e(i)
          read(lu,1009) gamma2(i),gamma2e(i),pbdot2(i),pbdot2e(i),
     +      si2(i),si2e(i),r2(i),r2e(i)
	endif
100	continue

999	npsr=i-1
	return
	end
