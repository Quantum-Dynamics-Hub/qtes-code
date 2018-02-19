        subroutine dumpef(myid,aaa, Np,Nd,E,Frc)
        implicit real*8 (a-h,o-z)
        character *10 aaa  
        include 'bohm.inc'
!        parameter (Ndx=250)
        real*8 E(*), Frc(Ndx,*)
        write(1000+myid,*)'*** myid: ',myid, aaa

1100    format(I3,I6' : ',E14.6,' F= ',300(' ',E14.6))          
        do i=1,np 
           write(1000+myid,1100)myid,i,E(i),(Frc(j,i),j=1,Nd)
           enddo
        end   
