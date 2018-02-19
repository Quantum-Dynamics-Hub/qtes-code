      subroutine read_input_mpi(nd,np,iat,pow,
     &          myid,numprocs,ist,x,p,s,w,tc,isr)
!-------------------------------------------------------------
!------- read input file: ndx, ndq, npx are from bohm.inc ----
!------- SETS BUFFERS FOR TRAJECTORIES  ----------------------
!------- TRAJs are not generated here, but  in init_traj.f ---
!- IF TRAJ.SAV IS READ, TRAJS ARE NOT PASSED CORRECTLY TO main
!------- NEED TO READ THEM AGAIN via set_traj.f --------------
!-- DEFINE t=tc=0 if ist.eq.0; else tc is from traj.sav ------
!-------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include 'bohm.inc'
      include 'mpif.h'
      integer*4,intent(out):: isr
      character*2 asymbol(ndx),asym0
      real*8 V(10*ndx + nqx*nqx + 30 )
      integer iat(ndx)
        real*8 x(ndx,npx),p(ndx,npx),s(npx),w(npx)
      common/wf/a0(nqx,nqx),q0(ndx),p0(ndx)
      common/prop/dt,dt2,tau(ndx),tau2(ndx),am(ndx)
     &          ,ipot,ipf,nqd,kmax,kout,krst
      common/molecule/asymbol,nat
      common/pot_xyz/ qc(3,ndx),ql,qls,dw,sc,esh,dch,zh,ch,bnd
      common/traj_cut/qmn,qmx
       dl=1.0d0 !      lattice constant

c
c   MPI rank  0 reads  input files:       
c
      if (myid.eq.0) then 
!cc read atomic xyz-positions  <<<<<<<<<<<<<<<<<<<<<<<
        open(8,file='INxyz')
        read(8,*) ql,qls !lattice constant and interaction distance cutoff
        read(8,*) dw,sc  !force consant dw and coupling  in units of dw    
        read(8,*) dch,zh,ch  !force constant and equilib distance for H
        read(8,*) esh   !shift the potential
        read(8,*) bnd   !distance when H is considered bonded to C
        read(8,*) qmn,qmx ! z-range for the first atom
        sc=dw*sc
!------------------------------------------------------------------------
!cc ql,qls,dw,sc,dch,zh can be used in any away (used in derivs0 only) 
!cc ch, esh, bnd are used after call to deriv and in aver,
!can be adjusted to fit the problem 
!------------------------------------------------------------------------
!cc read input file <<<<<<<<<<<<<<<<<<<<<<<
        open(9,file='INdyn')
        read(9,*) nat
!----------- positions of atoms, iat is read  in --------------------------
        bohr2A =  0.52917724d0
        bohr2A =  0.52918d0 
        A2bohr =1.0d0/bohr2A
!---- read geometry from INgeom in Anstrom for Nat atoms ---------
!--- and convert to bohr staight away ----------------------------
        do i=1,nat
         read(8,*) asym0,qc(1,i),qc(2,i),qc(3,i),iat(i)
         qc(1,i) =  qc(1,i) *A2bohr
         qc(2,i) =  qc(2,i) *A2bohr
         qc(3,i) =  qc(3,i) *A2bohr
        enddo
!--------------------------------
        read(9,*) nq
        read(9,*) kmax
        read(9,*) kout
        read(9,*) dt
        read(9,*) ipot,isr
        read(9,*) pow
        read(9,*) ipf
        read(9,*) np,ist
!-------------- total number of dimensions  ---------------------
        nd=nat*3
        nqd=nq*3
        call zero_mat(nqd,a0)
        krst=0
!----------------------------------------------------------------
!----------- atom type, initial WF, mass of each atom -----------
!----- q0 is with respect to equilibrium position from INgeom ---
!----------------------------------------------------------------
        do 129 n=1,nq
         i1=1+(n-1)*3
         i2=i1+1
         i3=i2+1
         read(9,*) asymbol(n),a0(i1,i1),q0(i1),q0(i2),q0(i3),
     &          p0(i1),p0(i2),p0(i3),am(i1)
!------------- displacements from equilbrium are added ----------
!------------ to get Cartesian coordinates in a.u. --------------
         q0(i1)=q0(i1)+qc(1,n)
         q0(i2)=q0(i2)+qc(2,n)
         q0(i3)=q0(i3)+qc(3,n)
c----------------------------------------------------------------
         a0(i2,i2)=a0(i1,i1)
         a0(i3,i3)=a0(i1,i1)
         am(i2)=am(i1)
         am(i3)=am(i1)
!--------- for imaginary time dynamics use real function, p0=0 ------
!         p0(i1)=0d0
!---------------------------------------------------------------------
129     continue
        do 139 n=nq+1,nat
         i1=1+(n-1)*3
         i2=i1+1
         i3=i2+1
         read(9,*) asymbol(n),t50,q0(i1),q0(i2),q0(i3),
     &          p0(i1),p0(i2),p0(i3),am(i1)
!------------- displacements from equilbrium are added ----------------------------
         q0(i1)=q0(i1)+qc(1,n)
         q0(i2)=q0(i2)+qc(2,n)
         q0(i3)=q0(i3)+qc(3,n)
!-------------------------------------------------------------------
         am(i2)=am(i1)
         am(i3)=am(i1)
!--------- for imaginary time dynamics use real function, p0=0 ------
!         p0(i1)=0d0
!---------------------------------------------------------------------
139    continue
       close(9)
       close(8)
!-------------- define time-steps ----------------------------------------
       tc=0d0
       dt=dt/kout
       dt2=dt/2d0
       do 91 i=1,nd
         tau2(i)=dt2/am(i)
         tau(i)=dt/am(i)
91     continue
!-------------- print out setup --------------------------------------
         write(*,*) 'LAGRANGIAN DYNAMICS in REAL TIME'
         if(ipf.eq.0) write(*,*) 'EHRENFEST: AVERAGE OVER QDOF',ipf
         if(ipf.ne.0) write(*,*) 'DO NOT AVERAGE OVER QDOF',ipf
         if(ipot.eq.0) write(*,*) 'NO R-MOMENTUM FIT: QPOT=0'
         if(ipot.eq.1) write(*,*)'LINEAR R-MOMENTUM FIT: LINEAR Q-FORCE'
         if(ipot.eq.2) write(*,*)'LINEAR R_MOMENTUM FIT: LINEAR Q-FORCE
     &                            ALONG Z'
         write(*,*) 'TOTAL NUMBER OF ATOMS ',nat
         write(*,*) 'DIMENSIONALITY ', nd,' QUANTUM DIM=',nqd
         write(*,*)nq,' QUANTUM ATOMS: ', (asymbol(j),j=1,nq)
         write(*,*)'ATOM TYPE: ', (iat(j),j=1,nq)
         write(*,*)nat-nq,' CLASSICAL ATOMS: ',(asymbol(j),j=nq+1,nat)
         write(*,*)'ATOM TYPE: ', (iat(j),j=nq+1,nat)
         write(*,*)'TIME STEPS, SUBSTEPS, dt: ',
     &          kmax,kout,dt
        write(*,*) 'LEAVE MPI INPUT'
c---------- SG read trajectory arrays ---------------
        if(ist.ne.0) then
1008     format(700(e23.16,1x))
         ic=0
         open(32,file='traj.sav')
         read(32,*,END=101) tc
         krst=int(tc/(dt*kout))
         kmax=kmax-tc/(dt*kout)
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
        endif


         endif     !<------ myid=0 done with reading
c
c   now everyone sets points in the buffer         
c

c integer array
        ind      =  1
        inp      =  ind      + 1 
        ikmax    =  inp      + 1
        ikout    =  ikmax    + 1
        iipot    =  ikout    + 1
        iipf     =  iipot    + 1
        inqd     =  iipf     + 1
        inat     =  inqd     + 1

        inq      =  inat     + 1
        iiat     =  inq      + 1

        iasymbol = iiat      + ndx 
        iendInt  = iasymbol  + ndx 
 
c  floats 
        iql      =  iendInt 
        iqls     =  iql      +   1
        idw      =  iqls     +   1
        isc      =  idw      +   1
        idch     =  isc      +   1
        izh      =  idch     +   1
        ich      =  izh      +   1
        iesh     =  ich      +   1
        ibnd     =  iesh     +   1 
        iqmn     =  ibnd     +   1
        iqmx     =  iqmn     +   1

        ipow     = iqmx      +      1
        idt      = ipow      +      1 
        idt2     = idt       +      1
        itau     = idt2      +      1
        itau2    = itau      +      ndx
        iam      = itau2     +      ndx 
        iq0      = iam       +      ndx
        ip0      = iq0       +      ndx   
        iqc      = ip0       +      ndx
        ia0      = iqc       +    3*ndx  
        iend     = ia0       +  nqx*nqx 

        mdv = 10*ndx + nqx*nqx + 30         ! buffer size (see declariotn)
        call tstcor(iend,mdv)
!-------------------------------------------------------------

c
c   now copy data to buffer        
c

      if  (myid.eq.0) then  
c   copy integers & characters to the buffer
        call cp_i2i(nd,      V(ind),      1)
        call cp_i2i(np,      V(inp),      1)
        call cp_i2i(kmax,    V(ikmax),    1)
        call cp_i2i(kout,    V(ikout),    1)
        call cp_i2i(ipot,    V(iipot),    1)
        call cp_i2i(ipf,     V(iipf),     1)
        call cp_i2i(nqd,     V(inqd),     1)
        call cp_i2i(nat,     V(inat),     1)
                                      
        call cp_i2i(nq,      V(inq),      1)
        call cp_i2i(iat,     V(iiat),     ndx)
                                      
        call cp_a2a(asymbol, V(iasymbol), ndx)

c
c   copy floats to the buffer
        call cp_d2d(ql,  V(iql),         1)                     
        call cp_d2d(qls, V(iqls),        1)                    
        call cp_d2d(dw,  V(idw),         1)                     
        call cp_d2d(sc,  V(isc),         1)                     
        call cp_d2d(dch, V(idch),        1)                    
        call cp_d2d(zh,  V(izh),         1)                     
        call cp_d2d(ch,  V(ich),         1)                     
        call cp_d2d(esh, V(iesh),        1)                    
        call cp_d2d(bnd, V(ibnd),        1)                    
        call cp_d2d(qmn, V(iqmn),        1)                    
        call cp_d2d(qmx, V(iqmx),        1)                    

        call cp_d2d(pow,  V(ipow),       1)         
        call cp_d2d(dt,   V(idt),        1)        
        call cp_d2d(dt2,  V(idt2),       1)        
        call cp_d2d(tau,  V(itau),       ndx)      
        call cp_d2d(tau2, V(itau2),      ndx)       
        call cp_d2d(am,   V(iam),        ndx)       
        call cp_d2d(q0,   V(iq0),        ndx)       
        call cp_d2d(p0,   V(ip0),        ndx)       
        call cp_d2d(qc,   V(iqc),      3*ndx)       
        call cp_d2d(a0,   V(ia0),    nqx*nqx)      
        endif


c  distribute it
        call MPI_Bcast(V(1),iend,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,
     $               ierr)


c  ...  and finally copy from   buffer  to the right places after the transfer
      if (myid.gt.0) then
        call cp_i2i(V(ind),      nd,         1)
        call cp_i2i(V(inp),      np,         1)
        call cp_i2i(V(ikmax),    kmax,       1)
        call cp_i2i(V(ikout),    kout,       1)
        call cp_i2i(V(iipot),    ipot,       1)
        call cp_i2i(V(iipf),     ipf,        1)
        call cp_i2i(V(inqd),     nqd,        1)
        call cp_i2i(V(inat),     nat,        1)
                                           
        call cp_i2i(V(inq),      nq,         1)
        call cp_i2i(V(iiat),     iat,      ndx)

        call cp_a2a(V(iasymbol), asymbol,  ndx)

        call cp_d2d(V(iql),      ql,         1)                     
        call cp_d2d(V(iqls),     qls,        1)                    
        call cp_d2d(V(idw),      dw,         1)                     
        call cp_d2d(V(isc),      sc,         1)                     
        call cp_d2d(V(idch),     dch,        1)                    
        call cp_d2d(V(izh),      zh,         1)                     
        call cp_d2d(V(ich),      ch,         1)                     
        call cp_d2d(V(iesh),     esh,        1)                    
        call cp_d2d(V(ibnd),     bnd,        1)                    
        call cp_d2d(V(iqmn),     qmn,        1)                    
        call cp_d2d(V(iqmx),     qmx,        1)                    

        call cp_d2d(V(ipow),     pow,        1)         
        call cp_d2d(V(idt),      dt,         1)        
        call cp_d2d(V(idt2),     dt2,        1)        
        call cp_d2d(V(itau),     tau,        ndx)      
        call cp_d2d(V(itau2),    tau2,       ndx)       
        call cp_d2d(V(iam),      am,         ndx)       
        call cp_d2d(V(iq0),      q0,         ndx)       
        call cp_d2d(V(ip0),      p0,         ndx)       
        call cp_d2d(V(iqc),      qc,       3*ndx)       
        call cp_d2d(V(ia0),      a0,     nqx*nqx)      
        endif

        return

        end

c----------------------
      subroutine cp_a2a(A,B,N)
      character*2 A(N), B(N)
      integer i,N
      do i =1,N
          B(i) = A(i)
          enddo
      end


c----------------------
!      subroutine cp_d2d(A,B,N)
!      real*8 A(N),B(N)
!      integer i,N
!
!      do i =1,N
!          B(i) = A(i)
!          enddo
!      end
!
!c----------------------
!      subroutine cp_i2i(A,B,N)
!      integer A(N),B(N)
!      integer i,N
!
!      do i =1,N
!          B(i) = A(i)
!          enddo
!      end
!c--------------------
!c testwhether enough space allocated      
!      subroutine  tstcor(iend,MDV)
!      integer iend,mdv
!
!      if (iend.gt.mdv)  then
!         write(*,*) 'stop: not enough space.'
!         stop
!         endif
!
!      end

