!-------------------------------------------------------------
!-------------------------------------------------------------
        subroutine set_traj(nd,iat,np,pow,x,p,s,w,u,du,v,dv,
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

        if(myid.eq.0) then

        write(*,*) 'READING TRAJ INFO via set_traj'

1008     format(700(e23.16,1x))
         ic=0
         open(32,file='traj.sav')
         read(32,*,END=101) tc
         do i=1,np
          ic=ic+1
          read(32,1008,END=101) (x(k,i),k=1,nd),w(i)
          read(32,1008,END=101) (p(k,i),k=1,nd),s(i)
         enddo
101      if(ic.ne.np) then
          write(*,*) 'NOT ENOUGH TRAJ DATA PROVIDED'
          call mpi_abort()
         endif
         close(32)

        endif  !<---- myid=0 is done with reading

!---------------------------------------------------------------
c distribute work
        nwork = (np + numprocs-1) /numprocs
        ip0 = myid*nwork  + 1 
        ip1 = myid*nwork  + nwork
        if(ip1.gt.np)  ip1 = np     ! the last one has the least work to do  
        
!---------------- compute classical potential ------------------
        nat1 = nd/3
         call mysecond(ta)
        do 112 ip =ip0,ip1
         !call derivs0(nd,iat,x(1,i),vlocal(i),dvlocal(1,i))
         !call derivs0(nd,iat,x(1,ip),vlocal(ip),dvlocal(1,ip))       !<--- old(fit)
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
        call qpot(np,nd,t,w,x,s,u,du,qp,myid,numprocs)
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
        if(npr.gt.40) npr=40
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
