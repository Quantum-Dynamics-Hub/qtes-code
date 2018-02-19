!--------------------------------------------------
!--------------------------------------------------
! ipot definition: 
!       0:= QP=0; 1:= global LQF for all QDOFs; 
!       2-3:= uncoupled fit for each quantum z-coord
!       2:= linear basis, 3:= quadratic basis; 
!       ELSE: the fit is 
!       coupled within (x,y,z) of a quantum atom 
!       but uncoupled between the atoms; 
!       4:= linear basis; 
!       7:= diagonal quadratic;
!       10:= full quadratic basis
!       basis order f=[x,y,z,x^2,y^2,z^2,xy,xz,yz,1]
!----------------------------------------------------
!----------------------------------------------------

!-----------------------------------------------------
!-----------------------------------------------------
        subroutine ehrnf(nd,nqd,np,x,w,dv)
	implicit real*8 (a-h,o-z)
!------ averages  gradV for CDOFs over QDOFs --------- 
        include 'bohm.inc'
	real*8 x(ndx,npx),dv(ndx,npx),w(npx)

        do k = nqd+1,nd
         t50 = 0d0
         do i=1,np
          t50 = t50 + dv(k,i)*w(i)
         enddo
         do i=1,np
          dv(k,i) = t50
         enddo
        enddo

        return
        end
!-----------------------------------------------------
!-----------------------------------------------------
        
	subroutine derivs_ho(nd,nqd,x,v,dv)
	implicit real*8 (a-h,o-z)
!----------- test potential -------------------------
        include 'bohm.inc'
	real*8 x(ndx),dv(ndx),qc(ndx)
	parameter(om=0.0136d0,om2=1d0)

	v = 0d0
	do i=1,nd
	 qc(i) = 0d0
	 dv(i) = 0d0
	enddo
	qc(3) = 5d0
	qc(6) = 5d0

	do i=1,nqd
	 v = v+om*(x(i)-qc(i))**2/2d0
	 dv(i) = om*(x(i)-qc(i))
	enddo

	do i=nqd+1,nd
	 v = v+om2*(x(i)-qc(i))**2/2d0
	 dv(i) = om2*(x(i)-qc(i))
	enddo

	return
	end



!---------- PUT INTO SEPARATE FILE ------------
c-----------------------------------------------------------
c-------- linear basis to fit gradA/A ----------------------
c-- generates quadratic Qpot and linear quantum force CQF --
c-- atoms are coupled -- 08/31/2016 SG ---------------------
c-----------------------------------------------------------
        subroutine lqf(np,nd,t,w,q,s,v,dv,qp,myid,numprocs)
        implicit real*8 (a-h,o-z)
        include 'bohm.inc'
        parameter(nd2=nqx+1,lwork=200)
        real*8 q(ndx,npx),v(npx),dv(nqx,npx),w(npx),s(npx)
        real*8 m2(nbx,nbx),b(nbx,nbx),f(nbx),d1(nbx)
        real*8 r(ndx,npx),p(ndx,npx)
        real*8 af(nbx,nbx),work(lwork),ferr(nbx),sw(nbx),berr(nbx)
        integer*4 iwork(lwork)
        character*1 equed
        real*8 qcn(ndx)
        common/prop/dt,dt2,tau(ndx),tau2(ndx),am(ndx)
     &          ,ipot,ipf,nqd,kmax,kout,krst
        common/chk/xsv(ndx*npx),psv(ndx,npx),ssv(npx),wsv(npx),ist
!------------------------------------------------------------------
        if(ipot.ne.1) return

        qp = 0d0

        anq = 0d0		!<< define shift for the basis 
        do i=1,nqd
         qcn(i) = 0d0
        enddo

        do 10 i=1,np
         anq = anq+w(i)
         do 10 k=1,nqd
          qcn(k) = qcn(k)+q(k,i)*w(i)
10      dv(k,i)=0d0

        do i=1,nqd
         qcn(i) = qcn(i)/anq
        enddo

        nc = nqd+1	     !<< zero matrices 

        do j=1,nc
         do i=1,nc
          m2(i,j) = 0d0
          b(i,j) = 0d0
         enddo
        enddo

        do n=1,np               !<< begin QT LOOP 

         call basis0(nc,qcn,q(1,n),f)
         do j=1,nc              !<< overlap matrix 
          do i=1,nc
           m2(i,j) = m2(i,j)+f(i)*f(j)*w(n) 
          enddo
         enddo
        enddo           !>> finish QT LOOP

        do j=1,nqd             !<< RHS vector 
         b(j,j) = -0.5d0*m2(nc,nc)
        enddo
	

c        if(myid.eq.0) then
c        do j=1,nc
c         write(777,*) (m2(i,j),i=1,nc)
c        enddo
c        endif

        call dposv('L',nc,nqd,m2,nbx,b,nbx,info)
c        call dgelsy(nc,nc,nqd,m2,nbx,b,nbx,iwork,rcond,irank,
c     &  work,lwork,info)

	if(info.ne.0) then
	 write(*,*) 'LQF problem',info
c,rcond,irank
	 stop
	endif
	
!------------ construct qpot and gradients -----------------
!---------------- analytical gradient ---- <<<<<<<<<<<<<<<<<
	do n=1,np		!<< begin QT LOOP 
         call basis0(nc,qcn,q(1,n),f)

         do j=1,nqd		!<< external loop for QP
          rf = abvec(nc,f,b(1,j))
          v(n) = v(n)-(rf*rf+b(j,j))/2d0/am(j)

	  do i=1,nqd		!<< internal loop for gradQP
	   dv(i,n) = dv(i,n)-rf*b(i,j)/am(j)
	  enddo

	 enddo

         qp = qp + v(n)*w(n)
        enddo                   !>> finish QT LOOP
	

	if(myid.eq.0) then
c	 if(krst.eq.0) then
	  write(14,*) t,(qcn(i),i=1,nqd)	!<< write fitting coefs
     &		,((b(i,j),i=1,nc),j=1,nqd)
	  flush(14)
c	 endif
	endif

        return
        end
        

c-----------------------------------------------------------
         subroutine basis0(nc,q,x,f)
        implicit real*8 (a-h,o-z)
!------- global lienar basis f={x,y,z...,1}-----
        real*8 x(nc),f(nc),q(nc)
        do j=1,nc-1
         f(j) = x(j)-q(j)
        enddo
        f(nc) = 1d0
        return
        end
c-----------------------------------------------------------
c-----------------------------------------------------------

	
c-----------------------------------------------------------
c-------- quadratic basis to fit gradA/A -------------------
c-- generates quartic Qpot and cubic quantum force CQF ----
c-- atom by atom -- no interatomic coupling  08/31/2016 SG -
c-----------------------------------------------------------
c        subroutine cqf(np,nd,t,w,q,p,r,s,v,dv,qp,myid,numprocs)
        subroutine cqf(np,nd,t,w,q,s,v,dv,qp,myid,numprocs)
        implicit real*8 (a-h,o-z)
        include 'bohm.inc'
        parameter(del = 1d-4,nd2=nqx+1,lwork=200)
        real*8 q(ndx,npx),v(npx),dv(nqx,npx),w(npx),s(npx)
        real*8 m2(nbx,nbx),b(nbx,3),f(nbx),d1(nbx)
        real*8 r(ndx,npx),p(ndx,npx)
        real*8 af(nbx,nbx),work(lwork),ferr(nbx),sw(nbx),berr(nbx)
        integer*4 iwork(lwork)
        character*1 equed
        real*8 b3(nbx,3),d(nbx,3),df(nbx,3),qcn(ndx)
        common/prop/dt,dt2,tau(ndx),tau2(ndx),am(ndx)
     &          ,ipot,ipf,nqd,kmax,kout,krst
        common/chk/xsv(ndx*npx),psv(ndx,npx),ssv(npx),wsv(npx),ist
!------------------------------------------------------------------
        qp = 0d0

        anq = 0d0
        do i=1,nqd
         qcn(i) = 0d0
        enddo

        do 10 i=1,np
         v(i)=0d0
         anq = anq+w(i)
         do 10 k=1,nqd
          qcn(k) = qcn(k)+q(k,i)*w(i)
10      dv(k,i)=0d0

        do i=1,nqd
         qcn(i) = qcn(i)/anq
        enddo
        
        if(ipot.lt.1.or.ipot.gt.10) return

        if(ipot.eq.1) then
         call lqf(np,nd,t,w,q,s,v,dv,qp,myid,numprocs)

	 return
        endif

        nq = nqd/3
        nc = ipot

!-------------- polynomial  basis in z-direction only  ----------
!------------------- one atom at a time -------------------------
        if(ipot.eq.2.or.ipot.eq.3) then

!------- nc=3 is quadratic basis f={z,z*z,1} --------------------
!----- degree of polynomial, can be changed to higher  ----------
!-------- nc=2 is linear basis f={z,1} in the z-coordinate only

!----------- optimized one quantum atom at a time ---------------
!---------- loop over the quantum atoms ------------------
        do kkk = 1,nq

         k1 = kkk*3

!---------- zero matrices --------------------------------
         do j=1,nc
          do i=1,nc
           m2(i,j)=0d0
          enddo
          do i=1,3
           b(j,i) = 0d0
          enddo
         enddo
!-------------- build up matrices, d1==df/dz ----
        do k=1,np
         ww=w(k)
         call basis3(nc,qcn(k1),q(k1,k),f,d1)
         do j = 1,nc
          do i = 1,nc
           m2(i,j) = m2(i,j)+f(i)*f(j)*ww
          enddo
          b(j,1)= b(j,1)-d1(j)*ww/2d0
c---------- LSF of clas and nonclas momenta ------
c------ not fully implemented --------------------
c          b(j,2)= b(j,2)+f(j)*r(k1,k)*ww
c          b(j,3)= b(j,3)+f(j)*p(k1,k)*ww
         enddo
        enddo

!----- find optimal coefficients ------------------------
c        nrhs = 3
        nrhs = 1
        call dposv('L',nc,nrhs,m2,nbx,b,nbx,info)
c----------- various solvers ---------------------<<<<<<<<
c        call dposvx('E','L',nc,nrhs,m2,nbx,af,nbx,equed,sw,
c     &  b,nbx,b3,nbx,rcond,ferr,berr,work,iwork,info)
c        call dgelss(nc,nc,nrhs,m2,nbx,b,nbx,sw,rcond,irank,
c     &  work,lwork,info)
c--------------------------------------------------->>>>>>
c        call dgelsy(nc,nc,nrhs,m2,nbx,b,nbx,iwork,rcond,irank,
c     &  work,lwork,info)
c---------------------------------------------
        if(info.ne.0) then
         write(*,*) 'DPOSV FAILED'
         stop
        endif

!------------- QUANTUM POTENTIAL -------------
!--------- kf=1 is LQF (integral by parts) ---
!--------- kf = 2 is r/p fitting ----------------
        kf = 1
!---------- print the fitting coeffs ---------
        if (myid.eq.0)  then
c	  if(krst.eq.0) then
	   kp = 40+kkk
           write(kp,1000) t,qcn(k1),qcn(k1+1),qcn(k1+2),
     &		(b(j,kf),j=1,nc)
           flush(kp)
c	  endif
        endif

        do 16 n=1,np
!----------- fitted momentum -----------------
         z = q(k1,n)
         z0 = qcn(k1)
         call basis3(nc,z0,z,f,d1)
         rf = abvec(nc,b(1,kf),f)
         dr = abvec(nc,b(1,kf),d1)
         v(n) = -(rf*rf+dr)/2d0/am(k1)

!----------- numerical derivative in z -------
         z = z +del
         call basis3(nc,z0,z,f,d1)
         rf = abvec(nc,b(1,kf),f)
         dr = abvec(nc,b(1,kf),d1)
         vp = -(rf*rf+dr)/2d0/am(k1)

         z = z -del*2d0
         call basis3(nc,z0,z,f,d1)
         rf = abvec(nc,b(1,kf),f)
         dr = abvec(nc,b(1,kf),d1)
         vm = -(rf*rf+dr)/2d0/am(k1)

         dv(k1,n) = (vp-vm)/2d0/del


         qp=qp+v(n)*w(n)
16      continue

c------- finish loop over quantum z's -------------
        enddo


!---------------------------------------------------
        else


!-------- ipot > 3 ---------------------------------
!-------- correlated 3D basis up to quadratic ------
!--------  one atom at a time ---------------
!------- f={x,y,z,x*x,y*y,z*z,x*y,x*z,y*z,1} -------
!------ nc = 4 is linear basis ---------------------
!------ nc = 7 is diagonal quadratic ---------------
!------ nc = 10 is full quadratic ------------------

!----- loop over quantum atoms ---------------------
        do kkk=1,nq		!<< LOOP OVER Q-ATOMS


c-------- index for x of the quantum atom kkk=[1,Nq]
         k1 = 3*(kkk-1)+1

         do j=1,nc
          do i=1,nc
           m2(i,j) = 0d0
          enddo
          do i=1,3
           b(j,i) = 0d0
          enddo
         enddo

c------ build up matrices ----------------------------
c-- d(i,1)=df(i)/dx d(i,2)=df(i)/dy,d(i,3)=df(i)/dz --

         do n=1,np 		!<<  QT LOOP 

          call basis10(nc,qcn(k1),q(k1,n),f,d)
          ww = w(n)
          w2 = ww/2d0

          do i=1,nc
c--------- overlap matrices ----------------
           do j=1,nc
            m2(i,j) = m2(i,j)+f(i)*f(j)*ww
           enddo
c---- 3D gradient --------------------------
	   do j=1,3
            b(i,j) = b(i,j)-d(i,j)*w2
	   enddo
	  enddo

         enddo			!>> END QT LOOP
!-------------------------------------
!---- do the fitting -----------------
!-------------------------------------
        nrhs = 3

        call dposv('L',nc,nrhs,m2,nbx,b,nbx,info)
!------ alternate solver ------
c        call dgelsy(nc,nc,nrhs,m2,nbx,b,nbx,iwork,rcond,irank,
c     &  work,lwork,info)

        if(info.ne.0) then
         write(*,*) 'DPOSV FAILED'
         stop
        endif
        if (myid.eq.0)  then
c	  if(krst.eq.0) then
	   kp = 40+kkk
           write(kp,1000) t,qcn(k1),qcn(k1+1),qcn(k1+2),
     &		((b(j,i),j=1,nc),i=1,3)
           flush(kp)
c	  endif
        endif

!------------- QUANTUM POTENTIAL ------------
!------ grad QP is numerical --------------------------

        do n=1,np 	!<< START QT LOOP 
!----------- fitted  noncl momentum -----------

!------ centered difference for the gradient -- 
         x = q(k1,n)
         y = q(k1+1,n)
         z = q(k1+2,n)
         call qp_fit(nc,am(k1),qcn(k1),b,x,y,z,v0)
         v(n) = v(n)+v0

!----------  derivatives in x ------------------
         x = x+del
         call qp_fit(nc,am(k1),qcn(k1),b,x,y,z,vp)
         x = x-del*2d0
         call qp_fit(nc,am(k1),qcn(k1),b,x,y,z,vm)
         dv(k1,n) = (vp-vm)/2d0/del
         x = x+del

!----------  derivatives in y ------------------
         y = y+del
         call qp_fit(nc,am(k1),qcn(k1),b,x,y,z,vp)
         y = y-del*2d0
         call qp_fit(nc,am(k1),qcn(k1),b,x,y,z,vm)
         dv(k1+1,n) = (vp-vm)/2d0/del
         y = y+del

!----------  derivatives in z ------------------
         z = z+del
         call qp_fit(nc,am(k1),qcn(k1),b,x,y,z,vp)
         z = z-del*2d0
         call qp_fit(nc,am(k1),qcn(k1),b,x,y,z,vm)
         dv(k1+2,n) = (vp-vm)/2d0/del
         
         qp = qp+v(n)*w(n)

         enddo           !<< END QT  LOOP 

	enddo           !<< END Q-ATOM LOOP

	endif		!<< finish qpot selection

1000    format(80(e12.5,1x))
1001    format(80(e16.9,',',1x))


        return
        end


!--------------------------------------------------------------
        subroutine qp_fit(nc,am,qc,b,x,y,z,v)
c------------ computes analytical qpot using -----
c------- fitted coefficients  b(nbasis,3) --------
c------- trajectory position q(3) ----------------
        include 'bohm.inc'
        implicit real*8 (a-h,o-z)
        real*8 q(3),b(nbx,3),f(nc),d(nbx,3),am(3),qc(3)

         q(1) = x
         q(2) = y
         q(3) = z

	 v = 0d0
         call basis10(nc,qc,q,f,d)
         do i=1,3
          rf = abvec(nc,f,b(1,i))
          dr = abvec(nc,d(1,i),b(1,i))
          v = v - (rf*rf+dr)/2d0/am(i) 
         enddo

        return
        end

!-------------------------------------------------------------
        function abvec(n,a,b)
!-------------------------------------------------------------
!------  dot-product of two vectors  -------------------------
        implicit real*8 (a-h,o-z)
        real*8 a(n),b(n)

        abvec = 0d0

        do i=1,n
         abvec = abvec+a(i)*b(i)
        enddo

        return
        end
        
!--------------------------------------------------------------
        subroutine basis10(nc,qc,q,f,d)
!----- defines a quadratic basis in 3 Cartesian dimensions ----
!---- options to return at linear or quadratic diagonal -------
        include 'bohm.inc'
        implicit real*8 (a-h,o-z)
        real*8 qc(3),q(3),f(nc),d(nbx,3),x(3)

c------------ shift the basis ----------
	do i=1,3
	 x(i)=q(i)-qc(i)
	enddo

        do j=1,nc
	 f(j) = 0d0
         do i=1,3
          d(j,i) = 0d0
         enddo
        enddo
c------ constant -------------------------------
        f(nc) = 1d0

c------- linear terms --------------------------
        do i=1,3
         f(i) = x(i)
         d(i,i) = 1d0
        enddo
        
        if(nc.eq.4) return

c------ diagonal quadratic terms ----------------
        do i=1,3
         f(i+3) = x(i)**2
         d(i+3,i)=2d0*x(i)
        enddo

        if(nc.eq.7) return

c------ off-diagonal quadratic: xy,xz,yz --------
        f(7) = x(1)*x(2)
        d(7,1) = x(2)
        d(7,2) = x(1)
        f(8) = x(1)*x(3)
        d(8,1) = x(3)
        d(8,3) = x(1)
        f(9) = x(2)*x(3)
        d(9,2) = x(3)
        d(9,3) = x(2)

        return
        end

!---------------------------------------------
        subroutine basis3(nc,q0,q,f,d) 
!---------------------------------------------
!----- 1D polynomial basis of size nc --------
        implicit real*8 (a-h,o-z)
        real*8 f(nc),d(nc)

        do j=1,nc
         d(j)=0d0
        enddo

        f(nc)=1d0
        f(1)=q-q0
        d(1)=1d0

!-----   if the basis is beyond linear do this
        do j=2,nc-1
         f(j)=(q-q0)**j
         d(j)=j*(q-q0)**(j-1)
        enddo

        return
        end
