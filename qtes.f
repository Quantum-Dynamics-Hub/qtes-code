	 program main
!---------------------------------------------------------------
!cccc MonteCarlo sampling, Gaussian deviates, matrix QPOT cccc
!--------- NOTE on the POTENTIAL -------------------------------
!--------- to make system calls to analytical potentials -------
!----  change 'call derivs' to 'call derivs0,1,2,3 -------------
!--  derivs() calls dftb and has some dftb keywords and Tel ---- 
!--  derivx() is a template for external call of the potential -
!---------------------------------------------------------------
	implicit real*8 (a-h,o-z)
        include 'bohm.inc'
        include 'mpif.h'
        complex*16 im,ov
!	parameter(ndx=250,npx=100000)
        character*2 asymbol(ndx)
        character*60 fstat
        dimension x(ndx,npx),p(ndx,npx),s(npx),u(npx),du(nqx,npx),
     &		w(npx),dx(ndx),xmn(ndx),q(ndx,ndx),v(npx),dv(ndx,npx)
        integer iat(ndx)
        dimension tim(npx)
        parameter(BohrtoA=0.52918d0,Hartoev=27.211d0)
        common/wf/a0(nqx,nqx),q0(ndx),p0(ndx),isr
        common/prop/dt,dt2,tau(ndx),tau2(ndx),am(ndx),
     &          ipot,ipf,nqd,kmax,kout,krst
        common/molecule/asymbol,nat
        common/chk/xsv(ndx*npx),psv(ndx,npx),ssv(npx),wsv(npx),ist
        character dum*80, name1*1,name2*2, name3*3,name4*4
        !parameter (MDV=10000000)
        !real*8  VV(MDV)             !scratch array / buffer
        character infile*80, outfile*80, blabel*70
        logical ex
        real*8 ts,tf
        
        integer myid, ierr, numprocs, root

        call mpi_init(ierr)
        call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
        call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
        
        include 'dftb/mpiio.inc'
        infile='IN_dftb'
        outfile='out.dftb'


11     format(I1)
12     format(I2)
13     format(I3)
14     format(I4)

        if (myid.eq.0) then 
          write(*,*)'############################################'
          write(*,*)
          write(*,*)' Parallel version of bohmd. Ver(17/08/2016)'
          write(*,*)
          write(*,*)'############################################'
          endif 

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!------------- read the input file -------------------------------------
        call read_input_mpi(nd,np,iat,pow,myid,numprocs,ist,x,p,s,w,t
     $       ,isr)

        call mpi_bcast(ist,1,MPI_INTEGER,0,MPI_COMM_WORLD,
     $        ierr)
        call mpi_bcast(isr,1,MPI_INTEGER,0,MPI_COMM_WORLD,
     $        ierr)
        call mpi_bcast(x,ndx*np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,
     $        ierr)
        call mpi_bcast(p,ndx*np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,
     $        ierr)
        call mpi_bcast(s,np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,
     $        ierr)
        call mpi_bcast(w,np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,
     $        ierr)

!------ setup dftb, read  SK files     -------------------------------     

        if (myid.eq.0)  then 
          inquire(file=infile,exist=ex)
          if(ex)then
              open(istdin,file=infile,status='old')
            else
              write(*,*)myid,' cannot find file', infile
              write(*,*)' stop!'
              call mpi_abort()
            endif
          open(istdout,file=outfile,status='unknown')    !<-----  to be deleted
          endif
        call dftb_init_mpi(ntype)                        !<--  just  returns  ntype

!-------------------- open files --------------------------------------

        if (myid.eq.0) then 
          if(ist.eq.0) then
           write(*,*)  'INITIALIZE QTs',ist,kmax,np,numprocs
           fstat='unknown'
          else
           write(*,*)  'READ QTs from file',ist,kmax,np,numprocs
           fstat='old'
          endif
          open(21,file='trajt')
          open(22,file='dftbt')
        endif 

!----------------------------------------------------------------------
!--------------------- initialize trajectories ------------------------

        if(ist.eq.0) then
        if (myid.eq.0) then
          open( 10,file=      'traj',status=fstat)
          open( 11,file=      'aver',status=fstat)
          open( 12,file=      'qt.dat',status=fstat)
          open( 14,file=      'pfit',status=fstat)
          open( 15,file=       'qtx',status=fstat)
          open( 16,file=       'qty',status=fstat)
          open( 17,file=       'qtz',status=fstat)
          open( 19,file=      'ener',status=fstat)
          open( 25,file=       'clx',status=fstat)
          open( 26,file=       'cly',status=fstat)
          open( 27,file=       'clz',status=fstat)
          open( 29,file=       'wft',status=fstat)
!          open(30,file='cm',status=fstat)
!          open(31,file='bnd',status=fstat)
          open( 13,file=       'cor',status=fstat)
        endif
        call mysecond(ts)
         call init_traj(nd,iat,np,pow,x,p,s,w,u,du,v,
     &                  dv,myid,numprocs)
        call mysecond(tf)
        write(1000+myid,*) kmax,np,(x(1,k),k=1,np)
        if (myid.eq.0) write(21,*) ts,tf,tf-ts
        call flush(21)

!------------------ compute average positions of atoms --------------
        call aver(np,nd,nqd,t,w,s,x,p,myid,numprocs)
        call over(nd,np,x,t,p,s,w,v,dv,u,0,myid,numprocs)
!----------------- print traj #2  -----------------------------------
        if (myid.eq.0.and.np.ge.2) then
         write(10,1000) t,(x(k,2),k=1,4*nqd),s(2)
         flush(10)
        endif

!------------------ read trajectories ------------------------
        else
        if (myid.eq.0) then
          open( 10,file=      'traj', ACCESS = 'APPEND',status=fstat)
          open( 11,file=      'aver', ACCESS = 'APPEND',status=fstat)
          open( 12,file=      'qt.dat', ACCESS = 'APPEND',status=fstat)
          open( 14,file=      'pfit', ACCESS = 'APPEND',status=fstat)
          open( 15,file=       'qtx', ACCESS = 'APPEND',status=fstat)
          open( 16,file=       'qty', ACCESS = 'APPEND',status=fstat)
          open( 17,file=       'qtz', ACCESS = 'APPEND',status=fstat)
          open( 19,file=      'ener', ACCESS = 'APPEND',status=fstat)
          open( 25,file=       'clx', ACCESS = 'APPEND',status=fstat)
          open( 26,file=       'cly', ACCESS = 'APPEND',status=fstat)
          open( 27,file=       'clz', ACCESS = 'APPEND',status=fstat)
          open( 29,file=       'wft', ACCESS = 'APPEND',status=fstat)
          open( 13,file=       'cor', ACCESS = 'APPEND',status=fstat)
        endif
        call mysecond(ts)
        if(myid.eq.0) write(*,*) 'READ QTs:Ntr, kmax, procs',np,kmax,numprocs
         call set_traj(nd,iat,np,pow,x,p,s,w,u,du,v,
     &                  dv,myid,numprocs)
        call mysecond(tf)
        if (myid.eq.0) write(21,*) ts,tf,tf-ts
        call flush(21)
        if(myid.eq.0)  write(*,*) 'NUMBER OF TIMESTEPS ON RESTART', kmax
        endif

!------------------ compute average positions of atoms --------------
!        call aver(np,nd,nqd,t,w,s,x,p,myid,numprocs)
!        call over(nd,np,x,t,p,s,w,v,u,myid,numprocs)
!----------------- print traj #2  -----------------------------------
!        if (myid.eq.0) then
!         write(10,1000) t,(x(k,2),k=1,4*nqd),s(2)
!         flush(10)
!        endif
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!---------------- dynamics loop -------------------------------------
!--------------------------------------------------------------------
	do 300 i=1,kmax
         if(myid.eq.0) write(*,*) 't=',t,i
	 do 301 ii=1,kout
          call trot(np,nd,iat,t,w,x,p,s,v,dv,u,du,myid,numprocs)
          call over(nd,np,x,t,p,s,w,v,dv,u,i,myid,numprocs)
301      continue
         call aver(np,nd,nqd,t,w,s,x,p,myid,numprocs)
         if (myid.eq.0) then
          write(10,1000) t,(x(k,2),k=1,4*nqd),s(2)
          flush(10)
          open(32,file='traj.sav')
          write(32,*) t
          do n=1,np
           write(32,1008) (x(k,n),k=1,nd),w(n)
           write(32,1008) (p(k,n),k=1,nd),s(n)
           write(12,1008) (x(k,n),k=1,nd),w(n)
           write(12,1008) (p(k,n),k=1,nd),s(n)
          enddo
          close(32)
         endif
300	continue
!-----------  print final positions ---------------------------------
        nwr=10000
        if(np.lt.nwr) nwr=np
        npr=nd
        if(npr.gt.20) npr=20

        if (myid.eq.0) then
          write(*,*) 'PRINTED IN ener: time,tot_en,q_pot,
     &          cl_pot,kin_qm,kin_cl'
          write(*,*) 'PRINTING POSITIONS UP TO ',npr,' DIMENSIONS'
          do  i=1,nwr
           write(29,1000) (x(j,i),j=1,npr)
          enddo 
          
        write(*,*) 'SG: CLOSING FILES'
          close(10)
          close(11)
          close(12)
          close(14)
          close(15)
          close(16)
          close(17)
          close(19)
          close(25)
          close(26)
          close(27)
          close(29)
!          close(30)
!          close(31)
          close(13)
          close(istdin)
          close(istdout)      !<----  to be deleted
        write(*,*) 'SG: CLOSED FILES'
          endif 
c         close(95+myid)
!--------------------------------------------------------------------
1000    format(200(e12.5,1x))
1008    format(700(e23.16,1x))
        !call mpi_barrier(MPI_COMM_WORLD,ierr)
        !write(*,*) 'MPI_BARRIER ERR',ierr
        call mpi_finalize(ierr)
        write(*,*) 'MPI_FINALIZE ERR',ierr
        stop
        end
!-------------------------------------------------------------
!-------------------------------------------------------------



!-------------------------------------------------------------
!-------------------------------------------------------------
        subroutine init_traj(nd,iat,np,pow,x,p,s,w,u,du,v,dv,
     $             myid,numprocs)
!------- read input file: ndx, ndq, npx are from bohm.inc ----
        implicit real*8 (a-h,o-z)
        include 'bohm.inc'
        include 'mpif.h'
        real*8 vlocal(npx), dvlocal(ndx,npx) 
        character*2 asymbol(ndx)
        integer iat(ndx),idum(ndx)
        real*8 x(ndx,npx),p(ndx,npx),s(npx),w(npx),dx(ndx)
        real*8 v(npx),dv(ndx,npx),u(npx),du(nqx,npx),q(ndx),pq(ndx)
!  dver(0)       : ehrenfest averaged potential over the wavefunction
!  dver(1:ndx)   : ehnrefest averaged forces    over the wavefunction
        real*8 dver(0:ndx)
        real*8 ta,tb
        common/wf/a0(nqx,nqx),q0(ndx),p0(ndx),isr
        common/prop/dt,dt2,tau(ndx),tau2(ndx),am(ndx),
     &          ipot,ipf,nqd,kmax,kout,krst
        common/molecule/asymbol,nat
        common/pot_xyz/ qc(3,ndx),ql,qls,dw,sc,esh,dch,zh,ch,bnd
!----------- normalize and sample WF  ------------------------------------
        tc=0d0
        pi=4d0*atan(1d0)
        vol=1.d0
        an=1d0
        do 91 i=1,nd
         q(i)=q0(i)
         pq(i)=0d0
         dx(i)=1.d0/dsqrt(4.d0*a0(i,i))
         if(i.le.nqd) then
          vol=vol*dx(i)*2d0
          an=an*dsqrt(2d0*a0(i,i)/pi)
         endif
91      continue
!------------ print the  wavepacket parameters ----------------------------
        if (myid.eq.0) then 
          write(*,*) 'INITITAL WAVEFUNCTION: sampling cutoff = ',pow
          write(*,*) 'width:'
	  write(*,*) (a0(i,i),i=1,nqd)
	  write(*,*) 'position: '
          write(*,*) (q0(i),i=1,nqd)
  	  write(*,*) 'momentum: '
          write(*,*) (p0(i),i=1,nqd)
          write(*,*) 'GASDEV RANDOM SAMPLING USING ',np,' POINTS'
          write(*,*) 'dx',(dx(i),i=1,nqd)
          write(*,*) 'VOL=',vol*an,an
          endif
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!----------------------- set up trajectories-----------------------------
        anrm=0.d0
        en0=0.d0
        epot=0.d0
        ekin=0d0
        do n=1,nqd
         idum(n)=-n
        enddo
!---------- random sampling; Gaussian deviates <<<<<<<<<<<<<<<<<<
        ic=0
        do 101 i=1,np*50
         ww=0d0
         do 102 k=1,nqd
          dev=gasdev(idum(k))
          if(isr.eq.1) then
          q(k)=q0(k)+dx(k)*dev
!-----------sampling for momentum(a mimic to quantum potential)-------
          if(k.eq.(nqd+1))then
          pq(k)=p0(k)+2d0*a0(k,k)*dx(k)*dev
          else
          pq(k)=p0(k)
          endif
!--------------------------------------------------------------------
          elseif(isr.eq.0) then
          q(k)=q0(k)
          pq(k)=p0(k)
          endif
          ww=ww+2d0*a0(k,k)*(q(k)-q0(k))**2
102      continue
         ww=exp(-ww)
         if(ww.gt.10**(-pow)) then
          ic=ic+1
         endif       
         if(ic.gt.np) goto 105
         s(ic)=0d0
         w(ic)=1d0/(np*1d0)
         do 103 k=1,nd
          p(k,ic)=pq(k)
          x(k,ic)=q(k)
103      continue
c--------------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
101     continue
!---------------------------------------------------------------
105     if (myid.eq.0) write(*,*) 'FINISHED SAMPLING',ic,np,'OUT OF',i
        if(ic.lt.np) np=ic
c distribute work
        nwork = (np + numprocs-1) /numprocs
        ip0 = myid*nwork  + 1 
        ip1 = myid*nwork  + nwork
        if(ip1.gt.np)  ip1 = np     ! the last one has the least work to do  
        
!---------------- compute classical potential ------------------
        nat1 = nd/3
         call mysecond(ta)
        do 112 ip =ip0,ip1
         !call derivs0(nd,iat,x(1,ip),vlocal(ip),dvlocal(1,ip))       !<--- old(fit)
c         call derivs_ho(nd,nqd,x(1,ip),vlocal(ip),dvlocal(1,ip))       !<--- old(fit)
         call derivs(nat1,iat,x(1,ip),vlocal(ip),dvlocal(1,ip))       !<--- new(dftb)
         vlocal(ip) = vlocal(ip) -esh                                 !<-- this probably should  be moved  to myid=0
         
112     continue
         call mysecond(tb)
        if (myid.eq.0) write(22,*) ta,tb,tb-ta
        call flush(22)
c
c  now collect all potentials and from  all Mpi tasks         
        nlocal = ip1-ip0+1
        call mpi_gather(vlocal(ip0),nlocal,MPI_DOUBLE_PRECISION,
     $        v,nlocal,MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD,ierr)   

        call mpi_gather(dvlocal(1,ip0),ndx*nlocal,MPI_DOUBLE_PRECISION,
     $        dv,ndx*nlocal,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

        call mpi_barrier(MPI_COMM_WORLD,ierr)

        call mpi_Bcast(v,np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        
        call mpi_Bcast(dv,ndx*np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,
     $        ierr)
!--------------------------------------------------------
!--------------------- Ehrenfest averaging <<<<<<<<<<<<<<
!--------------------- added 09/24/2012 SG --------------
        if(ipf.eq.0) then
         do j=0,nd
          dver(j)=0d0
         enddo
         do i=1,np
          dver(0)=dver(0)+w(i)
          do j=nqd+1,nd
           dver(j)=dver(j)+dv(j,i)*w(i)
          enddo
         enddo
         do i=1,np
          do j=nqd+1,nd
           dv(j,i)=dver(j)/dver(0)
          enddo
         enddo
        endif
!------------------------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!-------- s(0) is absorbed into the sampling --------------------
        do 100 i=1,np
         ek=0d0
         do 108 n=1,nd
108      ek=ek+p(n,i)**2/2d0/am(n)
         epot=epot+v(i)*w(i)
         ekin=ekin+ek*w(i)
         anrm=anrm+w(i)
100     continue
         
        t=0d0
c        call qpot(np,nd,t,w,x,s,u,du,qp,myid,numprocs)
        call cqf(np,nd,t,w,x,s,u,du,qp,myid,numprocs)
        if (myid.eq.0) then 
          write(*,*) 'NORM=',anrm
          write(*,*) 'KIN energy = ',ekin/anrm
          write(*,*) 'POT energy = ',epot/anrm
          write(*,*) 'QP energy = ',qp/anrm
          write(*,*) 'TOTAL energy = ',(qp+epot)/anrm
          endif
        ener=(qp+epot)/anrm
!---------------------------------------------------------------------------
!----------- examine the initial point distribution ------------------------
        npr=nd
        if(npr.gt.20) npr=20
        if (myid.eq.0) then 
          open(28,file='wf0')
          do i=1,np
            write(28,1000) (x(j,i),j=1,npr)
          enddo
          close(28)
          endif 
!        write(*,*) 'LEAVE INITIALIZE'
!--------------------------------------------------------------------
1000    format(200(e12.5,1x))
1007    format(700(e14.7,1x))
        return
        end
!-------------------------------------------------------------
!-------------------------------------------------------------
        subroutine aver(n,nd,nqd,t,w,s,x,p,myid,numprocs)
        implicit real*8 (a-h,o-z)
!------------------ averages  calculations-------------------------
!        parameter(ndx=250,npx=100000)
        include 'bohm.inc'
        include 'mpif.h'
	real*8 w(n),x(ndx,npx),s(n),p(ndx,npx),xav(ndx)
        common/pot_xyz/ qc(3,ndx),ql,qls,dw,sc,esh,dch,zh,ch,bnd
        common/wf/a0(nqx,nqx),q0(ndx),p0(ndx)
!------------------ normalize ------------------------
c        if(t.lt.1d-10) then
c        write(*,*) t
c         an=0d0
c        do 9 i=1,n
c9        an=an+w(i)
c        do 8 i=1,n
c8        w(i)=w(i)/an
c        endif
!----------------------------------------------------------
        do 10 j=1,2*nd+1
10       xav(j)=0d0
        prob=0d0
        an=0d0
!----------------------------------------------------------
!----------------------------------------------------------
        do 13 i=1,n
         an=an+w(i)
         if(nd.ge.6) then
!c          rch=abs(abs(x(3,i)-x(6,i))-ch)
          rch=abs(x(3,i)-ch)
          if(rch.ge.bnd) prob=prob+w(i)
         endif
         do 11 j=1,nd
          xav(j)=xav(j)+x(j,i)*w(i)
          xav(nd+j)=xav(nd+j)+x(j,i)**2*w(i)
11      continue
         xav(2*nd+1)=xav(2*nd+1)+x(2,i)*x(3,i)*w(i)
13      continue
        do 12 j=1,nd
12       xav(j)=xav(j)/an
        prob=prob/an
c-----------------------------------------------------------
c
c   print        
c
!-----------------------------------------------------------
       if (myid.eq.0) then 
         npr = min(20,nd)
         write(11,1000) t,an,prob,(xav(i),i=1,npr)
         flush(11)
!----------- print 20 QTs for the first quantum atom  -------
         ip = min(n,20)
         write(15,1000) t,(x(1,j),j=1,ip)
         write(16,1000) t,(x(2,j),j=1,ip)
         write(17,1000) t,(x(3,j),j=1,ip)
         flush(15)
         flush(16)
         flush(17)
!----------- print  20 CTs for the first classical atom ------
         if(nd.gt.nqd) then
           n1=nqd+1
           write(25,1000) t,(x(n1,j),j=1,ip)
           write(26,1000) t,(x(n1+1,j),j=1,ip)
           write(27,1000) t,(x(n1+2,j),j=1,ip)
           flush(25)
           flush(26)
           flush(27)
           endif
!------------ if all atoms are quantum print coordinates of the second atom
         if(nd.eq.nqd) then
           n1=4
           write(25,1000) t,(x(n1,j),j=1,ip)
           write(26,1000) t,(x(n1+1,j),j=1,ip)
           write(27,1000) t,(x(n1+2,j),j=1,ip)
           flush(25)
           flush(26)
           flush(27)
           endif
         endif  
1000   format(1000(e12.5,1x))
       return
       end
!-------------------------------------------------------------
!-------------------------------------------------------------
        subroutine prob(n,nd,amt,w,s,x,am,y,cm)
        implicit real*8 (a-h,o-z)
!------------ probability  calculations ----------------------
!        parameter(ndx=250,npx=100000)
        include 'bohm.inc'
        include 'mpif.h'
	real*8 w(n),x(ndx,npx),s(n),y(ndx),cm(ndx),am(ndx)
        real*8 q(ndx,ndx)
c--- average positions of atoms y; bond distances q(i,j)-------
        nat=nd/3
        do 20 i=1,nat
         do 20 j=i+1,nat
20        q(i,j)=0d0
        do 9 i=1,nd
9       y(i)=0d0
        an=0d0
c----------- the trajectory loop -----------------------------
	do 10 i=1,n
         ww=w(i)
         do 11 j=1,nd
11        y(j)=y(j)+x(j,i)*ww
         an=an+ww
c-------------- bond distances -------------------------------
         do 21 j=1,nat
          js=1+(j-1)*3
          do 21 k=j+1,nat
           ks=1+(k-1)*3
           rr=(x(js,i)-x(ks,i))**2+(x(js+1,i)-x(ks+1,i))**2
     &          +(x(js+2,i)-x(ks+2,i))**2
          q(j,k)=q(j,k)+rr*ww
21       continue
10	continue
c--------------- normalize ----------------------------------
        do 22 j=1,nat
         do 22 k=1,nat
22      q(j,k)=sqrt(q(j,k)/an)
        do 12 i=1,nd
12       y(i)=y(i)/an
!ccc------ CofM positions repeated to nd -----------------------
        do 13 j=1,3
         cm(j)=0d0
         do 14 i=1,nat
          is=j+(i-1)*3
14        cm(j)=cm(j)+am(is)*y(is)
         cm(j)=cm(j)/amt
         do 15 i=2,nat
15        cm(j+(i-1)*3)=cm(j)
13      continue
!-------------- bond distances ------------------------------
        write(30,1000) cm(1),cm(2),cm(3),(y(i),i=1,nd)
        write(31,1000) ((q(i,j),j=i+1,nat),i=1,nat)
1000    format(450(e12.5,1x))
	return
	end
!-------------------------------------------------------------
!-------------------------------------------------------------
       subroutine trot(np,nd,iat,t,w,q,p,s,v,dv,u,du,myid,numprocs)
!  dver(0)       : ehrenfest averaged potential over the wavefunction
!  dver(1:ndx)   : ehnrefest averaged forces    over the wavefunction
        implicit real*8 (a-h,o-z)
!	parameter(ndx=250,npx=100000)
        include 'bohm.inc'
        include  'mpif.h'
        real*8 q(ndx,npx),p(ndx,npx),s(npx),u(npx),du(nqx,npx),
     &  w(npx),v(npx),dv(ndx,npx)
        real*8  dver(0:ndx), vlocal(npx), dvlocal(ndx,npx) 
                                                  
        complex*16 im,ov
        integer iat(*)
        parameter(im=(0d0,1d0))
        common/pot_xyz/ qc(3,ndx),ql,qls,dw,sc,esh,dch,zh,ch,bnd
        common/prop/dt,dt2,tau(ndx),tau2(ndx),am(ndx)
     &          ,ipot,ipf,nqd,kmax,kout,krst
        common/traj_cut/qmn,qmx
!--------------- trajectory propagation  --------------------
! half-step
        do  i=1,np
!---<<<<<<<< exclude the fringe QM trajectories <<<<<<<<<<<
c          if(q(3,i).lt.qmn) then
c           w(i)=0d0
c           goto 10
c          endif
c          if(q(3,i).gt.qmx) goto 10
!--------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

         do  k=1,nd             ! dimensions loop
          p(k,i)=p(k,i)-dv(k,i)*dt2
          if(k.le.nqd) p(k,i)=p(k,i)-du(k,i)*dt2
          q(k,i)=q(k,i)+p(k,i)*tau(k)
          s(i)=s(i)+p(k,i)**2*tau2(k)
         enddo

         s(i)=s(i)-(v(i)+u(i))*dt2
        enddo           ! end trajectory loop

!---------------- recompute all forces -----------------
        t = t+dt
c        call qpot(np,nd,t,w,q,s,u,du,qp,myid,numprocs)
        call cqf(np,nd,t,w,q,s,u,du,qp,myid,numprocs)

c
c  local MPI  loop,  distribute potential,forces  evaluation
c      Last task (myid=numprocs-1)  is under-employed  
      
        nwork = (np + numprocs-1) /numprocs   
        ip0 = myid*nwork  + 1
        ip1 = myid*nwork  + nwork   
        if(ip1.gt.np)  ip1 = np     ! the last one has the least work to do

        nat1=nd/3
        !do 24 i=1,np
        do 24 ip=ip0,ip1        
!---<<<<<<<< exclude the fringe QM trajectories <<<<<<<<<<<
          if(q(3,ip).lt.qmn) then
           w(n)=0d0
           goto 24
          endif
          if(q(3,ip).gt.qmx) goto 24
!--------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!         call derivs0(nd,iat,q(1,ip),vlocal(ip),dvlocal(1,ip))      !<-- old (fit)
c         call derivs_ho(nd,nqd,q(1,ip),vlocal(ip),dvlocal(1,ip))      !<-- old (fit)
         call derivs(nat1,iat,q(1,ip),vlocal(ip),dvlocal(1,ip))       !<-- new (dftb) 
         vlocal(ip) = vlocal(ip) -esh                                 !<-- this probably should  be moved  to myid=0
24      continue

c
c  now collect all  potenials and forces from  MPI tasks         
        nlocal = ip1-ip0+1
        call mpi_gather(vlocal(ip0),nlocal,MPI_DOUBLE_PRECISION,
     $        v,nlocal,MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD,ierr)   

        call mpi_gather(dvlocal(1,ip0),ndx*nlocal,MPI_DOUBLE_PRECISION,
     $        dv,ndx*nlocal,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 

        call mpi_barrier(MPI_COMM_WORLD,ierr)
         
        call mpi_Bcast(v,np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        
        call mpi_Bcast(dv,ndx*np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,
     $        ierr)
!--------------------- Ehrenfest averaging <<<<<<<<<<<<<<
!--------------------- added 09/24/2012 SG --------------
        if(ipf.eq.0) then
         do j=0,nd
          dver(j)=0d0
         enddo
         do i=1,np
          dver(0)=dver(0)+w(i)
          do j=nqd+1,nd
           dver(j)=dver(j)+dv(j,i)*w(i)
          enddo
         enddo
         do i=1,np
          do j=nqd+1,nd
           dv(j,i)=dver(j)/dver(0)
          enddo
         enddo
        endif
!------------------------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!-------------------- half-step in p -------------------
        do i=1,np       ! start trajectory loop

         do  k=1,nd     ! dimensions loop
          p(k,i)=p(k,i)-dv(k,i)*dt2
          if(k.le.nqd) p(k,i)=p(k,i)-du(k,i)*dt2
         enddo

         s(i)=s(i)-(v(i)+u(i))*dt2
        enddo           ! end trajectory loop


1000    format(20(e12.5,1x))
1007    format(700(e14.7,1x))
        return
       end
!-------------------------------------------------------------
!-------------------------------------------------------------
        subroutine over(nd,np,x,t,p,s,w,v,dv,u,imax,myid,numprocs)
!------------ compute average energy and autocorrelation -----
        implicit real*8 (a-h,o-z)
        include 'bohm.inc'
        include 'mpif.h'
        parameter(xmin=1.d0,xmax=10.d0,nx=250)
        parameter(pmin=-20.d0,pmax=20.d0,npz=1000)
        real*8 x(ndx,npx),p(ndx,npx),s(npx),w(npx),u(npx),
     &          v(npx),dv(ndx,npx)
        real*8 ek(npx),qs(npx),dis(npx),disz,numap,numrp,
     &          avp(ndx),avp2(ndx)
        real*8 xz(nx),pz(npz),numx(nx),nump(npz),angle(npx)
        real*8 disang(npx),distan(npx),orix(ndx,npx),adsdis(npx)
        complex*16 im,ov
        character*20 filepmv,filewfz,filepz
        parameter(im=(0d0,1d0),BohrtoA=0.52918d0,Hartoev=27.211d0)
	parameter(adsp=5.5d0,borp=11.d0)
        common/prop/dt,dt2,tau(ndx),tau2(ndx),am(ndx)
     &          ,ipot,ipf,nqd,kmax,kout,krst
        common/inixy/ hini(2,npx)
	common/enpt/ qss(npx),vs(npx)
!
! ov=<psi(-t)|psi(t)>; en=<H>; an=<psi|psi>; vav=<V>;qp=<U>
! enq=<K> for _quantum  subsystem; enc=<K> for classical  subsystem
!
        eks =0.0d0
        ov=(0d0,0d0)
        en=0d0
        an=0d0
        vav=0d0
        qp=0d0
        enq=0d0
        enc=0d0

        do 10 i=1,np
         ov=ov+w(i)*exp(2d0*im*s(i))
         an=an+w(i)
         en=en+(v(i)+u(i))*w(i)
         qp=qp+u(i)*w(i)
         vav=vav+v(i)*w(i)
         if(t.eq.0) vs(i)=v(i)
         ek(i)=0d0
         qs(i)=0d0
         do 11 k=1,nd
          if(k.le.nqd) qs(i)=qs(i)+p(k,i)**2/2d0/am(k)
11        ek(i)=ek(i)+p(k,i)**2/2d0/am(k)
         if(t.eq.0) qss(i)=qs(i)
         en=en+ek(i)*w(i)
         enq=enq+qs(i)*w(i)
         enc=enc+(ek(i)-qs(i))*w(i)
10      continue
        en=en/an
        qp=qp/an
        vav=vav/an
        enq=enq/an
        enc=enc/an
        ov=ov/an
!
c--------------- print  ------------------------------------
        if (myid.eq.0) then 
         write(19,1007) t,en,qp,vav,enq,enc
         flush(19)
         write(13,1007) 2d0*t,dreal(ov),imag(ov)
         flush(13)
        endif

1007    format(1200(e14.7,1x))
        return
        end
!-------------------------------------------------------------
!-------------------------------------------------------------
        subroutine basis(nb,nd,q,f)
!------------------------- fitting polynomial basis ----------
        implicit real*8 (a-h,o-z)
        real*8 f(nb),q(nd)
        do 10 i=1,nb
10       f(i)=0d0
!------------------------ linear basis ------------------------
        do 11 i=1,nd
11       f(i)=q(i)
        f(nb)=1d0
        return
        end
!-------------------------------------------------------------
!-------------------------------------------------------------
        subroutine qpot(np,nd,t,w,q,s,v,dv,qp,myid,numprocs)
        implicit real*8 (a-h,o-z)
        include 'bohm.inc'
        include 'mpif.h'
	parameter(eps=1d-5,nd2=nqx+1)
        real*8 q(ndx,npx),v(nqx),dv(nqx,npx),w(npx),s(npx)
        real*8 m2(nd2,nd2),b(nd2,nd2),f(nd2),pf(ndx)
        real*8 m3(nd2,nd2),b3(nd2,nd2)
        common/prop/dt,dt2,tau(ndx),tau2(ndx),am(ndx)
     &          ,ipot,ipf,nqd,kmax,kout,krst
        common/chk/xsv(ndx*npx),psv(ndx,npx),ssv(npx),wsv(npx),ist
!------------------------------------------------------------------
        do 10 i=1,np
         v(i)=0d0
         do 10 k=1,nqd
10      dv(k,i)=0d0
        qp=0d0
        if(ipot.eq.0) return
!------------------ if linear basis ---------------------
         n1=nqd+1
        nqz=nqd/3
         n3=nqz+1
        if (myid.eq.0) then 
          if(t.lt.1d-16) write(*,*) 'NQDim=',nqd,'  NBAS=',n1
          endif
        do 11 j=1,n1
         do 11 i=1,n1
          m2(i,j)=0d0
11      b(i,j)=0d0
!----------------- build up matrices --------------------
        do 214 k=1,np
         ww=w(k)
         call basis(n1,nqd,q(1,k),f)
        do 212 j=1,n1
         do 212 i=1,n1
212     m2(i,j)=m2(i,j)+f(i)*f(j)*ww
214     continue
        do 213 j=1,n1
213      b(j,j)=-m2(n1,n1)/2d0

c----------------3D quantum potential----------------
        if(ipot.eq.1) then

        call dposv('L',n1,nqd,m2,nd2,b,nd2,info)
        if(info.ne.0) then
         write(*,*) 'DPOSV FAILED'
         stop
        endif
        if (myid.eq.0)  then
         if(ist.eq.0) then
          write(14,1000) t,((b(j,i),j=1,n1),i=1,nqd)
         else
          if(t.gt.dt) write(14,1000) t,((b(j,i),j=1,n1),i=1,nqd)
         endif
          flush(14)
        endif
!------------- QUANTUM POTENTIAL ------------
        qp=0.d0
        do 16 n=1,np
!----------- fitted momentum -------------------
         call basis(n1,nqd,q(1,n),f)
         call dotvec(n1,nqd,f,b,pf)
         v(n)=0d0
!         do 222 i=1,nqd
!          do 224 j=1,nqd

c         disz = 0d0
c         disz =abs(q(3,n)-q(6,n))
c        if(n.eq.1) then
c        write(123,*) disz
c        endif

!        if(t.le.600)then
!         do 225 i=nqd,nqd
!          do 226 j=nqd,nqd
!226       dv(i,n)=dv(i,n)-b(i,j)*pf(j)/am(j)
!225       v(n)=v(n)-(pf(i)*pf(i)+b(i,i))/2d0/am(i)
!        else
!         if(disz.le.4.8d0)then
         do 222 i=1,nqd
          do 224 j=1,nqd
224       dv(i,n)=dv(i,n)-b(i,j)*pf(j)/am(j)
222       v(n)=v(n)-(pf(i)*pf(i)+b(i,i))/2d0/am(i)
!         else
!         do 221 i=nqd,nqd
!          do 223 j=nqd,nqd
!223       dv(i,n)=0.d0
!221       v(n)=0.d0
!         endif
!        endif
         qp=qp+v(n)*w(n)
16      continue

        else
c------------quantum potential along z--------------

        do i=1,n3
          do j=1,n3
            b3(i,j)=0d0
          enddo
        enddo

        do i=1,nqz
          do j=1,nqz
            m3(i,j)=m2(i*3,j*3)
          enddo
        enddo

        do i=1,nqz
          m3(i,n3)=m2(i*3,nqz*3+1)
        enddo

        do j=1,nqz
          m3(n3,j)=m2(nqz*3+1,j*3)
        enddo

        m3(n3,n3)=m2(nqz*3+1,nqz*3+1)

        do 203 j=1,n3
203      b3(j,j)=-m2(n1,n1)/2d0
        
!-------------------------------------
        call dposv('L',n3,nqz,m3,nd2,b3,nd2,info)
        if(info.ne.0) then
         write(*,*) 'DPOSV FAILED'
         stop
        endif
        if (myid.eq.0)  then
         if(ist.eq.0) then
          write(14,1000) t,((b3(j,i),j=1,n3),i=1,nqz)
         else
          if(t.gt.dt) write(14,1000) t,((b3(j,i),j=1,n3),i=1,nqz)
         endif
          flush(14)
        endif
!------------- QUANTUM POTENTIAL ------------
	qp=0.d0
	do 17 n=1,np
!----------- fitted momentum -------------------
         call basis(n1,nqd,q(1,n),f)
         call dotvec(n3,nqz,f,b3,pf)
         v(n)=0d0

         do 223 i=3,nqd,3
          do 225 j=3,nqd,3
225       dv(i,n)=dv(i,n)-b3(i/3,j/3)*pf(j/3)/am(j)
223       v(n)=v(n)-(pf(i/3)*pf(i/3)+b3(i/3,i/3))/2d0/am(i)
	 qp=qp+v(n)*w(n)
17	continue
        endif

1000    format(80(e12.5,1x))
1001    format(80(e16.9,',',1x))
        return
        end
c-------------------------------------------------------------
c-------------------------------------------------------------
        subroutine dotvec(n,nd,f,b,p)
c---------- computes dot product of two vectors 
        implicit real*8 (a-h,o-z)
        include 'bohm.inc'
	parameter(eps=1d-5,nd2=nqx+1)
        real*8 f(n),p(nd),b(nd2,nd2)
        common/prop/dt,dt2,tau(ndx),tau2(ndx),am(ndx)
     &          ,ipot,ipf,nqd,kmax,kout,krst
        if(ipot.eq.1) then
        do 10 i=1,nd
         p(i)=0d0
         do 10 j=1,n
10       p(i)=p(i)+f(j)*b(j,i)
        else
        do 11 i=1,nd
         p(i)=0d0
         do j=1,n
           if(j.eq.n) p(i)=p(i)+f((j-1)*3+1)*b(j,i)
           p(i)=p(i)+f(j*3)*b(j,i)
         enddo
11      continue
        endif
        return
        end
!
c-------------------------------------------------------------
c-------------------------------------------------------------
	subroutine determ(n,a,d)
	implicit real*8 (a-h,o-z)
ccc computes a determinant
        include 'bohm.inc'
	dimension a(nqx,nqx),b(nqx,nqx),indx(nqx,nqx)
	do 10 j=1,n
	 do 10 i=1,n
10	  b(i,j)=a(i,j)
	call ludcmp(b,n,nqx,indx,d)
	do 11 j=1,n
11	d=d*b(j,j)
	return
	end
c-------------------------------------------------------------
c-------------------------------------------------------------
        subroutine derivx(nd,x,v,dv)
        implicit real*8 (a-h,o-z)
c-------- system call to an external classical pot -----------
        include 'bohm.inc'
        parameter(btoang=0.529177249d0)
        character*2 asymbol(ndx)
        common/molecule/asymbol,nat
        real*8 x(nd),dv(nd),xa(ndx),dva(ndx)
c------------- convert to angstroms --------------------------
        do 10 i=1,nd
10       xa(i)=x(i)*btoang
c----------- input file for external program -----------------
        open(8,file='IN_pot')
        write(8,*) nd
        do 11 i=1,nd
11       write(8,*) xa(i)
        close(8)
c------------- external call ----------------------------------
        CALL SYSTEM('pot.x')
c-------------- read its  output ------------------------------
        open(12,file='pot.out')
        read(12,*) v
        do i=1,nd
         read(12,*) dva(i)
        enddo
        close(12)
c---------- convert to atomic unit ----------------------------
        do 12 i=1,nd
12       dv(i)=dva(i)*btoang
        return
        end
C=============================================================================
        subroutine derivs2(nd,iatt,x,v,dv)
        implicit real*8 (a-h,o-z)
        real*8 x(3,*),dv(3,*), dedx(3,100)
        include 'bohm.inc'

c------------- harmonic 5x5 lattice   potentials --------------------
        aa =0.18d0   ! harmonic  constant
        dl=1.4d0 !      lattice constant
        q  = 0.1
        qq = q*q

        v=0d0


        Nat = 9
        E = 0.0d0
        do i=1, Nat
          dedx(1,i) =0.0d0
          dedx(2,i) =0.0d0
          dedx(3,i) =0.0d0
          enddo

C  harmonic part
        iat =0     !init counter
        do i =0,2
            do j=0,2
              x0=  real(j)*dl
              y0=  real(i)*dl
              z0=0.0
              iat = iat +1
              xx = x(1,iat)
              yy = x(2,iat)
              zz = x(3,iat)
              dE    =  aa*((xx-x0)**2 + (yy-y0)**2+ (zz-z0)**2)
              E    = E  + dE
              dx = aa*2.0*(xx-x0)
              dy = aa*2.0*(yy-y0)
              dz = aa*2.0*(zz-z0)
              dEdx(1,iat)  = dedx(1,iat)  + dx
              dEdx(2,iat)  = dedx(2,iat)  + dy
              dEdx(3,iat)  = dedx(3,iat)  + dz
              !write(*,123)j,i,xx,yy,zz, dE,dx,dy,dz
              enddo
            enddo
123         format (I6,I6,100(' ',F16.8))
             !stop
c return
        v =E
        do i=1,Nat
           dv(1,i) = dedx(1,i)
           dv(2,i) = dedx(2,i)
           dv(3,i) = dedx(3,i)
!           write(*,100)i,dv(1,i), dv(2,i),dv(3,i)

           enddo
!           write(*,*)'jj:  E = ',  v
!           stop


!          write(*,*)'nat =',nat
           goto 1000
c  now add coulomb
        do i =1,Nat
           do j=1,i-1
            xx = x(1,i) - x(1,j)
            yy = x(2,i) - x(2,j)
            zz = x(3,i) - x(3,j)
            rsq = xx**2+yy**2 +zz**2
            E = E +  qq/sqrt(rsq)
            qqrsq = qq/rsq
            dedx(1,i) =  dedx(1,i) +  xx*qqrsq
            dedx(2,i) =  dedx(2,i) +  yy*qqrsq
            dedx(3,i) =  dedx(3,i) +  zz*qqrsq
            dedx(1,j) =  dedx(1,j) -  xx*qqrsq
            dedx(2,j) =  dedx(2,j) -  yy*qqrsq
            dedx(3,j) =  dedx(3,j) -  zz*qqrsq

!            write(*,*)'ij,rsq:',i,j, rsq
!             write(*,*)'i,j,zz',i,j, zz
            enddo
          enddo
!          stop

c return
1000    continue
        v =E
        do i=1,Nat
           dv(1,i) = dedx(1,i)
           dv(2,i) = dedx(2,i)
           dv(3,i) = dedx(3,i)
           enddo


!         write(*,*)'jj:  E = ',  v
!         write(*,*)'***  x:'
!         do  i=1,Nat
!           write(*,100)i,x(1,i), x(2,i),x(3,i)
!           enddo
!         write(*,*)'*** grad:'
!         do  i=1,Nat
!           write(*,100)i,dv(1,i), dv(2,i),dv(3,i)
!           enddo
!           write(*,*)' derivs2'
!           stop

100     format('jj:', I6, 100(' ' , E16.6))
        end
c-------------------------------------------------------------


C=============================================================================
C  coupled Morse in bond-distances
c
        subroutine derivs1(nd,x,v,dv)
        implicit real*8 (a-h,o-z)
c-------------------------------------------------------------
c------------- define Morse H2 potentials --------------------
        !parameter(d=160d0,z=1.04435d0,dw=174.5067d0,r0=1.4d0)
        parameter(d=0.17429d0,z=1.0435d0,dw=d*z*z,r0=1.4d0)
        real*8 x(3,*),dv(3,*)
        v=0d0
        E =0.0d0
        Nat=nd/3
        do  i=1,nat
          dv(1,i)=0.0d0
          dv(2,i)=0.0d0
          dv(3,i)=0.0d0
          enddo
c--------------------------------------------
        do i=1,Nat
          xi = x(1,i)
          yi = x(2,i)
          zi = x(3,i)
          do j=1,i-1
            xji= x(1,j) - xi
            yji= x(2,j) - yi
            zji= x(3,j) - zi
            rsq = xji*xji +yji*yji +zji*zji
            r = sqrt(rsq)
            ri =1.0/r
            ex = xji*ri
            ey = yji*ri
            ez = zji*ri
            t = exp(-z*(r-r0))
            E   = E + d*(t-1d0)**2
            Ep  = 2d0*d*z*t*(1d0-t)
            Epx = Ep*ex
            Epy = Ep*ey
            Epz = Ep*ez
            dv(1,j) =  dv(1,j)  + Epx
            dv(2,j) =  dv(2,j)  + Epy
            dv(3,j) =  dv(3,j)  + Epz
            dv(1,i) =  dv(1,i)  - Epx
            dv(2,i) =  dv(2,i)  - Epy
            dv(3,i) =  dv(3,i)  - Epz
            enddo
          enddo
        v = E
c--------------------------------------------
        if(v.gt.600d0) then
         v=600d0
        do i=1,nat
         dv(1,i)=0d0
         dv(2,i)=0d0
         dv(3,i)=0d0
        enddo
        endif
        return
        end
c-------------------------------------------------------------
c-------------------------------------------------------------
        subroutine derivs0(nd,iat,x,v,dv)
        implicit real*8 (a-h,o-z)
c-------------------------------------------------------------
c-------- classical potential --------------------------------
        include 'bohm.inc'
c------------- define Morse H2 potentials --------------------
!        parameter(d=0.17429d0,a=1.0435d0,dw=2d0*d*a*a,sc=0d-1*dw)
!        parameter(ql=1.4d0,qls=ql+1e-3)
        real*8 x(3,*),dv(3,*)
        integer iat(*)
        common/prop/dt,dt2,tau(ndx),tau2(ndx),am(ndx)
     &          ,ipot,ipf,nqd,kmax,kout,krst
        common/pot_xyz/ qc(3,ndx),ql,qls,dw,sc,esh,dch,zh,ch,bnd
        nat=nd/3
        v=-esh
        do 11 i=1,nat
         dv(1,i)=0d0
         dv(2,i)=0d0
11       dv(3,i)=0d0
!--------------------------------------------------------------------
!------------ interaction of QM particle 1 with CL particle 2 -------
!-------------- HO oscillator for CH bond ----------------------
c         xi=x(1,1)-qc(1,1)
c         yi=x(2,1)-qc(2,1)
c         zi=x(3,1)-x(3,2)-ch
c         v=v+dch*(xi*xi+yi*yi+zi*zi)/2d0
c         dv(1,1)=dv(1,1)+dch*xi
!!c         dv(1,2)=dv(1,2)-dch*xi
c         dv(2,1)=dv(2,1)+dch*yi
!!c         dv(2,2)=dv(2,2)-dch*yi
c         dv(3,1)=dv(3,1)+dch*zi
c         dv(3,2)=dv(3,2)-dch*zi
!-------------- MORSE oscillator for CH bond ----------------------
!!c         xi=x(1,1)-x(1,2)
!!c         yi=x(2,1)-x(2,2)
         xi=x(1,1)-qc(1,1)
         yi=x(2,1)-qc(2,1)
         zi=x(3,1)-x(3,2)
!!c         ri=exp(-zh*(abs(zi)-ch))
         ri=exp(-zh*(zi-ch))
         v=v+dch*(xi*xi+yi*yi)/2d0+dch*(ri-1d0)**2
         dv(1,1)=dv(1,1)+dch*xi
!!c         dv(1,2)=dv(1,2)-dch*xi
         dv(2,1)=dv(2,1)+dch*yi
!!c         dv(2,2)=dv(2,2)-dch*yi
!!c         dv(3,1)=dv(3,1)+zh*dch*ri*(1d0-ri)*zi/abs(zi)*2d0
!!c         dv(3,2)=dv(3,2)-zh*dch*ri*(1d0-ri)*zi/abs(zi)*2d0
         dv(3,1)=dv(3,1)+zh*dch*ri*(1d0-ri)*2d0
         dv(3,2)=dv(3,2)-zh*dch*ri*(1d0-ri)*2d0
!--------------------------------------------------------------------
!----------------- Coupled harmonic potential -----------------------
        do 12 j=2,nat
         xi=x(1,j)-qc(1,j)
         yi=x(2,j)-qc(2,j)
         zi=x(3,j)-qc(3,j)
         vl=dw*(xi*xi+yi*yi+zi*zi)/2d0
         v=v+vl
         dv(1,j)=dv(1,j)+dw*xi
         dv(2,j)=dv(2,j)+dw*yi
         dv(3,j)=dv(3,j)+dw*zi

         do 13 i=j+1,nat
          qq=(qc(1,i)-qc(1,j))**2+(qc(2,i)-qc(2,j))**2+(qc(3,i)-qc(3,j))**2
          qq=dsqrt(qq)
          if(qq.lt.qls) then
           dx=x(1,i)-x(1,j)
           dy=x(2,i)-x(2,j)
           dz=x(3,i)-x(3,j)
           xx=dsqrt(dx*dx+dy*dy+dz*dz)
           v=v+(xx-ql)**2*sc/2d0
           t50=sc*(xx-ql)/xx
           dv(1,i)=dv(1,i)+t50*dx
           dv(2,i)=dv(2,i)+t50*dy
           dv(3,i)=dv(3,i)+t50*dz
           dv(1,j)=dv(1,j)-t50*dx
           dv(2,j)=dv(2,j)-t50*dy
           dv(3,j)=dv(3,j)-t50*dz
          endif
13       continue
12      continue
c        write(*,*) 'DERS'
c        do i=1,nat
c         write(*,*) dv(1,i),dv(2,i),dv(3,i)
c        enddo
c        stop
!---- jacek
        !call dumpqef(nd,x,v,dv)
C - Sonya uncomment  the line below for scaling tests (1 sec. sleep)
        ! call sleep(1)
        end
c----------        
!     Jacek:
C  derivs0(nd,iat,x,v,dv)        
        subroutine dumpqef(nd,x,v,dv)
        implicit real *8 (a-h,o-z)
        real*8 x(*), dv(*)
111     format('Q=',18(' ',f12.6),' E=',E14.6,' Frc=',18(' ',E14.6))
        write(500,111) (x(i),i=1,nd),v,(dv(i),i=1,nd) 
        return
        end
c-------------------------------------------------------------
c-------------------------------------------------------------
c-------------------------------------------------------------
c-------------------------------------------------------------
        subroutine derivs3(nd,iat,x,v,dv)
        implicit real*8 (a-h,o-z)
c-------------------------------------------------------------
c-------- classical potential --------------------------------
        include 'bohm.inc'
c        parameter(ndx=250,npx=125000000)
c------------- define Morse H2 potentials --------------------
        parameter(d=0.17429d0,a=1.0435d0,dw=2d0*d*a*a,r0=0.0d0,sc=1d-1)
        real*8 x(nd),dv(nd),dvr(ndx),rx(ndx,ndx),u(ndx,ndx)
        integer iat(*)
        common/prop/dt,dt2,tau(ndx),tau2(ndx),am(ndx)
     &          ,ipot,ipf,nqd,kmax,kout,krst
        v=0d0
        do 11 i=1,nd
11       dv(i)=0d0
c----------------- Coupled harmonic potential -----------------------
        do 12 j=1,nd
         do 12 i=1,nd
12       u(i,j)=0d0  
        do 13 j=1,nd
         u(j,j)=dw
         u(j,j+1)=dw*sc
         u(j+1,j)=dw*sc
13      continue
        v=0d0
        do 14 j=1,nd
         dv(j)=0d0
         do 15 i=1,nd
15        dv(j)=dv(j)+u(j,i)*(x(i)-r0)
14       v=v+(x(j)-r0)*dv(j)/2d0
        return
        end
c-----------------------------------------------------------
c
c
c-----------------------------------------------------------
        subroutine derivs(nat,iat,C,Energy,frc)
        implicit real*8 (a-h,o-z)
        real*8  C(3,*), frc(3,*)
        integer iat(*),lmax(100)
        character*70 blabel
        logical lopf,scfhelp
        blabel=' '
        ibohm=1

!         inquire(101,opened=lopf)
!         if (.NOT.LOPF) then
!           open(101,file='run.dftb',status='unknown')
!         else
!           rewind(101)
!           endif

        ntype   = 2     ! max atoms type
!  lmax passed  from dftb_init via common         
!        lmax(1) = 2     ! carbon  
!        lmax(2) = 1     ! hydrogen
        qnel = 0.0
        scfhelp = .true.
        scftol = 1e-6
!        telec = 1000.0
         telec = 6000.0



           call do_dftb(ibohm,nat,iat,C,Frc,Eelec,Emermin,blabel,
     $         ntype, qnel,scfhelp, scftol,telec, V)
         energy =  Emermin

         do i=1,nat
           Frc(1,i) =  - Frc(1,i)
           Frc(2,i) =  - Frc(2,i)
           Frc(3,i) =  - Frc(3,i)
           enddo


          return
          end
c-------------------------------------------------------------
c-------------------------------------------------------------
      SUBROUTINE ludcmp(a,n,np,indx,d)
C------------------------------------------------------------
C   DOUBLE PRECISION VERSION
C------------------------------------------------------------ 
      INTEGER n,np,indx(n),NMAX
      REAL*8 d,a(np,np),TINY
      PARAMETER (NMAX=2500,TINY=1.0e-20)
      INTEGER i,imax,j,k
      REAL*8 aamax,dum,sum,vv(NMAX)
C
      d=1.
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        if (aamax.eq.0.) stop 'singular matrix in ludcmp'
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=TINY
        if(j.ne.n)then
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      END


      SUBROUTINE lubksb(a,n,np,indx,b)
C----------------------------------------------------------
C    DOUBLE PRECISION VERSION
C----------------------------------------------------------
      INTEGER n,np,indx(n)
      REAL*8 a(np,np),b(n)
      INTEGER i,ii,j,ll
      REAL*8 sum
C
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END

      FUNCTION ran1(idum)
C-----------------------------------------------------------------
C    DOUBLE PRECISION VERSION OF THE NUMREC ROUTINE
C    FRANK GROSSMANN, 12.9.1994
C-----------------------------------------------------------------
      INTEGER*4 idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL*8 ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1.d0/IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2d-14,RNMX=1.0d0-EPS)
      INTEGER*4 j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END
!---------------------------------------------------------------
      FUNCTION ran2(idum)
      INTEGER*4 idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL*8 ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER*4 idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END
!---------------------------------------------------------------
      FUNCTION ran4(idum)
      INTEGER*4 idum
      REAL*8 ran4
CU    USES psdes
      INTEGER*4 idums,irword,itemp,jflmsk,jflone,lword
      REAL*8 ftemp
      EQUIVALENCE (itemp,ftemp)
      SAVE idums,jflone,jflmsk
      DATA idums /0/, jflone /Z'3F800000'/, jflmsk /Z'007FFFFF'/
      if(idum.lt.0)then
        idums=-idum
        idum=1
      endif
      irword=idum
      lword=idums
      call psdes(lword,irword)
      itemp=ior(jflone,iand(jflmsk,irword))
      ran4=ftemp-1.d0
      idum=idum+1
      return
      END
!---------------------------------------------------------------
      SUBROUTINE psdes(lword,irword)
      INTEGER*4 irword,lword,NITER
      PARAMETER (NITER=4)
      INTEGER*4 i,ia,ib,iswap,itmph,itmpl,c1(4),c2(4)
      SAVE c1,c2
      DATA c1 /Z'BAA96887',Z'1E17D32C',Z'03BCDC3C',Z'0F33D1B2'/, c2 
     */Z'4B0F3B58',Z'E874F0C3',Z'6955C5A6', Z'55A7CA46'/
      do 11 i=1,NITER
        iswap=irword
        ia=ieor(irword,c1(i))
        itmpl=iand(ia,65535)
        itmph=iand(ishft(ia,-16),65535)
        ib=itmpl**2+not(itmph**2)
        ia=ior(ishft(ib,16),iand(ishft(ib,-16),65535))
        irword=ieor(lword,ieor(c2(i),ia)+itmpl*itmph)
        lword=iswap
11    continue
      return
      END
!---------------------------------------------------------------

      FUNCTION gasdev(idum)
C-----------------------------------------------------
C    DOUBLE PRECISION VERSION OF NUMREC ROUTINE
C    FRANK GROSSMANN, 12.9.1994
C-----------------------------------------------------
      INTEGER*4 idum
      REAL*8 gasdev
CU    USES ran1
      INTEGER*4 iset
      REAL*8 fac,gset,rsq,v1,v2,ran1
	real*4 ran3
      SAVE iset,gset
      DATA iset/0/
      if (iset.eq.0) then
1       v1=2.0d0*ran1(idum)-1.0d0
        v2=2.0d0*ran1(idum)-1.0d0
        rsq=v1**2+v2**2
        if(rsq.ge.1.0d0.or.rsq.eq.0.0d0)goto 1
        fac=dsqrt(-2.0d0*dlog(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      endif
      return
      END

        subroutine monc(Pi,l,y)
        implicit real*8(a-h,o-z)
        real*8 y(4)
        real*4 x(6)
        if (l.eq.1) call sobseq(-1,x)
10      call sobseq(6,x)
        v1=2.d0*x(1)-1.d0
        v2=2.d0*x(2)-1.d0
        u1=2.d0*x(4)-1.d0
        u2=2.d0*x(6)-1.d0
        sv=v1**2+v2**2
        su=u1**2+u2**2
        if (sv.ge.1d0.or.su.ge.1.d0) goto 10
         if (sv.lt.1d-40.or.su.lt.1d-40) then
          ru=0.d0
          rv=0.d0
          goto 11
         end if
        rv=sqrt(-2.d0*log(sv)/sv)
        ru=sqrt(-2.d0*log(su)/su)
11      y(1)=v1*rv
        y(2)=v2*rv
        y(3)=u1*ru
        y(4)=u2*ru
        return
        end


        subroutine sobseq(n,x)
        integer n,maxbit,maxdim
        real x(*)
        parameter (maxbit=30,maxdim=6)
        integer i,im,in,ipp,j,k,l,ip(maxdim),iu(maxdim,maxbit),
     *   iv(maxbit*maxdim),ix(maxdim),mdeg(maxdim)
        real fac
        save ip,mdeg,ix,iv,in,fac
        equivalence (iv,iu)
        data ip /0,1,1,2,1,4/, mdeg/1,2,3,3,4,4/, ix /6*0/
        data iv /6*1,3,1,3,3,1,1,5,7,7,3,3,5,15,11,5,15,13,9,156*0/
        if (n.lt.0) then
          do 11 k=1,maxdim
11          ix(k)=0
          in=0
          if(iv(1).ne.1) return
          fac=1./2**maxbit
          do 15 k=1,maxdim
            do 12 j=1,mdeg(k)
12            iu(k,j)=iu(k,j)*2**(maxbit-j)
            do 14 j=mdeg(k)+1,maxbit
               ipp=ip(k)
                  i=iu(k,j-mdeg(k))
                  i=ieor(i,i/2**mdeg(k))
                  do 13 l=mdeg(k)-1,1,-1
                        if (iand(ipp,1).ne.0) i=ieor(i,iu(k,j-1))
13                      ipp=ipp/2
14                iu(k,j)=i
15      continue
        else
                im=in
                do 16 j=1,maxbit
                  if(iand(im,1).eq.0) goto 1
16                im=im/2
!--- jacek                  
!                pause 'maxbit too small in sobseq'
                write(*,*) 'warning: maxbit too small in sobseq'
1               im=(j-1)*maxdim
                do 17 k=1,min(n,maxdim)
                        ix(k)=ieor(ix(k),iv(im+k))
17                      x(k)=ix(k)*fac
                in=in+1
        end if
        return
        end
!-------------------------------------------------------------
!-------------------------------------------------------------
        subroutine zero_mat(n,a)
        implicit real*8 (a-h,o-z)
        include 'bohm.inc'
!------------------------- zero matrix ------------------------
        real*8 a(nqx,nqx)
        do 11 j=1,n
         do 11 i=1,n
11        a(i,j)=0.d0
        return
        end
!---------------------------------
        subroutine aclear(N,A) 
        real*8 A(N)
        integer i,N
        do i=1,N
           A(i) = 0.0d0
           enddo
        end   
C-------------------------------------------------------------
C-------------------------------------------------------------
!-------------------------------------------------------------
        subroutine mysecond(t)
        implicit none
        real*8 t
        include 'mpif.h'

        t=MPI_WTIME()

        return
        end
c--------------------------------------------------------------
