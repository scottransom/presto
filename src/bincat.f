	program bincat

C  Reads pulsar catalog data from ascat1.dat (on lu 23), and writes it
C  to an unformatted binary file ascat.bin (on lu 24)

	include '../include/rdascat.h'

	open(24,file='../lib/bincat.dat',status='unknown',
     +    form='unformatted')

	call rdascat(23,npsr)

	write(24) jname,bname,ra2000,ra1950,
     +    rae,dec2000,dec1950,dece,nscode,
     +    dmin,dmax,dist,ndflag,lcode,ucode,
     +    ldeg,bdeg,pmra,pmrae,pmdec,
     +    pmdece,posepoch,p,pe,pdot,pdote,
     +    f2,f2e,f3,f3e,epoch,dm,dme,
     +    rm,
     +    rme,we,w50,w10,s400,s600,s1400,
     +    tau,ntauflag,t408,ntype,
     +    modcode,limcode,distmod,lum,bsurf,age,edot,ibin,
     +    pb,pbe,a1,a1e,om,ome,
     +    omdot,omdote,e,ee,t0,t0e,
     +    gamma,gammae,pbdot,pbdote,si,sie,
     +    r,re,
     +    pb2,pb2e,a12,a12e,om2,om2e,
     +    omdot2,omdot2e,e2,e2e,t02,t02e,
     +    gamma2,gamma2e,pbdot2,pbdot2e,si2,
     +    si2e,r2,r2e,npsr

	end
