module ivariables
implicit none

integer :: j,wn,wn1,wn3,wn4,clmax,ccwn
integer :: ci,cj,kkk,awn,awn1,clmaxind
integer :: cci,ccj,ccwn1,cclmax,cclmaxind 
integer :: mmm1,nnn1,mmmm2,np,iikk,nsum,rcut
integer :: ij,ml,jjn,nkk,nmkl,kjj,ij1,kk1

integer, parameter :: N =48720
real(8), dimension(3,N) :: R      
real(8), dimension(3) :: CellL    
!clmax clmaxind ccwn cs_dist(j,wn) anelem
integer, parameter :: P1A   =621!N
integer, parameter :: P2A   =622!N
integer, parameter :: P3A   =623!CA


integer, parameter :: S1A   =3721
integer, parameter :: S2A   =3722 
integer, parameter :: S3A   =3723
integer, parameter :: S4A   =3724
integer, parameter :: S5A   =3725


integer, parameter :: MG =311!326
integer, parameter :: ZN =1!326

integer, parameter :: Water_per =3
integer, parameter :: SO4_per =5

integer, parameter :: MG_per =1
integer, parameter :: ZN_per =1

integer, parameter :: wt_nnb =15000
integer, parameter :: so4_nnb =620

integer, parameter :: ZMG_nnb =310
integer, parameter :: ZN2_nnb =310
integer, parameter :: wn2 =1
integer, parameter :: iaa=310
integer, parameter :: maxbin2=201

integer,parameter  :: Step=50000

integer :: nbin1,nbin2,nbin3,nbin4,maxbin

real(8) :: r_Q12A(3,ZMG_nnb),r_T12A(3,ZN2_nnb)
real(8) :: r_P1A(3,wt_nnb),r_P2A(3,wt_nnb),r_P3A(3,wt_nnb),r_P4A(3,wt_nnb),r_P5A(3,wt_nnb),r_P6A(3,wt_nnb)
!real(8) :: r_P7A(3,alig),r_P8A(3,alig),r_P9A(3,alig),r_P10A(3,alig),r_P11A(3,alig),r_P12A(3,alig)

real(8) :: r_S1A(3,so4_nnb),r_S2A(3,so4_nnb),r_S3A(3,so4_nnb),r_S4A(3,so4_nnb),r_S5A(3,so4_nnb)


real(8):: nmgw(iaa,maxbin2),nznw(iaa,maxbin2)

real(8):: amg(iaa),azn(iaa),add1,dd,drbin,addmg,addzn

real(8):: nmgS(iaa,maxbin2),nznS(iaa,maxbin2)
real(8):: ntotmgw(iaa),ntotmgS(iaa),ntotznw(iaa),ntotznS(iaa)
real(8):: dnfacmg,dfacmg,facnewmg(iaa)
real(8):: dnfaczn,dfaczn,facnewzn(iaa)

real(8):: apmg(maxbin2),apzn(maxbin2)
real(8):: prfrnmg(iaa,maxbin2),prfrnzn(iaa,maxbin2)


!real(8):: xx1,yy1,zz1,rrr1_sq,rBMCl
real(8):: xxmgw,yymgw,zzmgw,rrrMGW_sq,rMGW
real(8):: xxznw,yyznw,zzznw,rrrZNW_sq,rZNW
real(8):: xxmgso4,yymgso4,zzmgso4,rrrmgso4_sq,rMGSO4
real(8):: xxznso4,yyznso4,zzznso4,rrrznso4_sq,rZNSO4

real(8):: sum,sum1,avgdist,distsq,sum33,sum34,sum35,sum32,sum36

end module ivariables

!##############################

program read_data
use ivariables

implicit none
integer :: Natom, status, nst
integer :: i, k, ii

! ## MOLFILE
integer, parameter :: Nfile = 1
character(len=72), dimension(Nfile) :: Filename   ! filename
integer, dimension(Nfile) :: Ns      ! # of steps
integer, dimension(Nfile) :: Handle
real, dimension(3*N) :: XYZ
real, dimension(6) :: Box
character(len=10), dimension(Nfile) :: InType
     !open(32,file="box-size.dat")
     !open(75,file="output-dist-BMIMCl.dat")
     !open(76,file="Prob-dist-BMIMCl.dat")
     open(32,file='MG-ZN-so4-water-codepref.dat')
     open(33,file='MG-ZN-so4-water-pref-sum.dat')
     !open(81,file="output-t-Runitvec.dat")
     !open(82,file="output-t-Rall.dat")
! ## Filename

    Filename( 1) = '../../PROD-NPT-1us-whole-RUN1.xtc'
   Ns(1) =500

   Handle(:) = -1          
   InType(:) = 'auto'      

   nst = 0 
 
   call f77_molfile_init_ 

   do i = 1 , Nfile

! 
    call f77_molfile_open_read_(Handle(i), Natom, FileName(i), InType(i))

        drbin=0.1
        rcut=20
        maxbin=anint(rcut/drbin) + 1
        do 132 jjn=1,iaa
         do 130 nkk=1,maxbin
           prfrnmg(jjn,nkk)=0.00
           prfrnzn(jjn,nkk)=0.00
           !
 130       continue
 132       continue
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     do j = 1 , Ns(i) 

       status = 1
       call f77_molfile_read_next_(Handle(i),Natom,XYZ(1),Box(1),status)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       do k = 1, N
         ii = (k-1)*3
         R(1,k) = dble( XYZ(ii+1) )
         R(2,k) = dble( XYZ(ii+2) )
         R(3,k) = dble( XYZ(ii+3) )
          !write(*,*) R(1,k),R(2,k),R(3,k),k,Ns(i)
          end do
       
       CellL(1) = dble(Box(1))
       CellL(2) = dble(Box(2))
       CellL(3) = dble(Box(3))
      !write(32,*) CellL(1),CellL(2), CellL(3)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       nst = nst + 1  
!# Analysis  
         do 166 mmm1=1,iaa
         do 135 nnn1=1,maxbin
         nmgw(mmm1,nnn1)=0
         nznw(mmm1,nnn1)=0
         nmgS(mmm1,nnn1)=0
         nznS(mmm1,nnn1)=0         
  135   continue
  166   continue




!       call read_Water
do wn=1, wt_nnb
         r_P1A(1,wn)=R(1,P1A+Water_per*(wn-1))
         r_P1A(2,wn)=R(2,P1A+Water_per*(wn-1))
         r_P1A(3,wn)=R(3,P1A+Water_per*(wn-1))
!
         r_P2A(1,wn)=R(1,P2A+Water_per*(wn-1))
         r_P2A(2,wn)=R(2,P2A+Water_per*(wn-1))
         r_P2A(3,wn)=R(3,P2A+Water_per*(wn-1))
!
         r_P3A(1,wn)=R(1,P3A+Water_per*(wn-1))
         r_P3A(2,wn)=R(2,P3A+Water_per*(wn-1))
         r_P3A(3,wn)=R(3,P3A+Water_per*(wn-1))
!
!  
end do

do wn=1, so4_nnb
         r_S1A(1,wn)=R(1,S1A+SO4_per*(wn-1))
         r_S1A(2,wn)=R(2,S1A+SO4_per*(wn-1))
         r_S1A(3,wn)=R(3,S1A+SO4_per*(wn-1))
!
         r_S2A(1,wn)=R(1,S2A+SO4_per*(wn-1))
         r_S2A(2,wn)=R(2,S2A+SO4_per*(wn-1))
         r_S2A(3,wn)=R(3,S2A+SO4_per*(wn-1))
!
         r_S3A(1,wn)=R(1,S3A+SO4_per*(wn-1))
         r_S3A(2,wn)=R(2,S3A+SO4_per*(wn-1))
         r_S3A(3,wn)=R(3,S3A+SO4_per*(wn-1))
!
         r_S4A(1,wn)=R(1,S4A+SO4_per*(wn-1))
         r_S4A(2,wn)=R(2,S4A+SO4_per*(wn-1))
         r_S4A(3,wn)=R(3,S4A+SO4_per*(wn-1))
!
!
         r_S5A(1,wn)=R(1,S5A+SO4_per*(wn-1))
         r_S5A(2,wn)=R(2,S5A+SO4_per*(wn-1))
         r_S5A(3,wn)=R(3,S5A+SO4_per*(wn-1))
!
end do


do wn=1, ZMG_nnb
         r_Q12A(1,wn)=R(1,MG+MG_per*(wn-1))
         r_Q12A(2,wn)=R(2,MG+MG_per*(wn-1))
         r_Q12A(3,wn)=R(3,MG+MG_per*(wn-1))
!
end do

do wn=1, ZN2_nnb
         r_T12A(1,wn)=R(1,ZN+ZN_per*(wn-1))
         r_T12A(2,wn)=R(2,ZN+ZN_per*(wn-1))
         r_T12A(3,wn)=R(3,ZN+ZN_per*(wn-1))
!
end do
   
!    do wn=1, ZMG_nnb
!!
!          xx1=r_T12A(1,wn)-r_Q12A(1,wn)
!          xx1=xx1-CellL(1)*anint(xx1/CellL(1))
!          yy1=r_T12A(2,wn)-r_Q12A(2,wn)
!          yy1=yy1-CellL(2)*anint(yy1/CellL(2))
!          zz1=r_T12A(3,wn)-r_Q12A(3,wn)
!          zz1=zz1-CellL(3)*anint(zz1/CellL(3))
!          rrr1_sq=(xx1)**2+(yy1)**2+(zz1)**2
!          rBMCl=sqrt(rrr1_sq)
!!
!!                do kkk=1, NNN
!!                if( dl + dd*(kkk-1) < rBMCl .and. rBMCl < dl + dd*kkk ) then
!!                        count_dis(kkk)=count_dis(kkk)+1
!!                end if
!!                end do
!           end do          
!
!-------------------------------------------
     do wn=1,ZMG_nnb 
     do wn1=1,wt_nnb   

         ! xx3=r_T12A(1,ZN2_nnb)-r_P1A(1,wn1)
         ! xx3=xx3-CellL(1)*anint(xx3/CellL(1))
         ! yy3=r_T12A(2,ZN2_nnb)-r_P1A(2,wn1)
         ! yy3=yy3-CellL(2)*anint(yy3/CellL(2))
         ! zz3=r_T12A(3,ZN2_nnb)-r_P1A(3,wn1)
         ! zz3=zz3-CellL(3)*anint(zz3/CellL(3))
         ! rrr3_sq=(xx3)**2+(yy3)**2+(zz3)**2
         ! rZNW=sqrt(rrr3_sq)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          xxmgw=r_Q12A(1,wn)-r_P1A(1,wn1)
          xxmgw=xxmgw-CellL(1)*anint(xxmgw/CellL(1))
          yymgw=r_Q12A(2,wn)-r_P1A(2,wn1)
          yymgw=yymgw-CellL(2)*anint(yymgw/CellL(2))
          zzmgw=r_Q12A(3,wn)-r_P1A(3,wn1)
          zzmgw=zzmgw-CellL(3)*anint(zzmgw/CellL(3))
          rrrMGW_sq=(xxmgw)**2+(yymgw)**2+(zzmgw)**2
          rMGW=sqrt(rrrMGW_sq)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          xxznw=r_T12A(1,wn)-r_P1A(1,wn1)
          xxznw=xxznw-CellL(1)*anint(xxznw/CellL(1))
          yyznw=r_T12A(2,wn)-r_P1A(2,wn1)
          yyznw=yyznw-CellL(2)*anint(yyznw/CellL(2))
          zzznw=r_T12A(3,wn)-r_P1A(3,wn1)
          zzznw=zzznw-CellL(3)*anint(zzznw/CellL(3))
          rrrZNW_sq=(xxznw)**2+(yyznw)**2+(zzznw)**2
          rZNW=sqrt(rrrZNW_sq)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !
!          xx3=r_T12A(1,ZN2_nnb)-r_P1A(1,wn1)
!          xx3=xx3-CellL(1)*anint(xx3/CellL(1))
!          yy3=r_T12A(2,ZN2_nnb)-r_P1A(2,wn1)
!          yy3=yy3-CellL(2)*anint(yy3/CellL(2))
!          zz3=r_T12A(3,ZN2_nnb)-r_P1A(3,wn1)
!          zz3=zz3-CellL(3)*anint(zz3/CellL(3))
!          rrr3_sq=(xx3)**2+(yy3)**2+(zz3)**2
!          rCLW=sqrt(rrr3_sq)
!
          nbin1=anint(rMGW/drbin) + 1
          nbin2=anint(rZNW/drbin) + 1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          if(nbin1 .le. maxbin) then      
!        
         nmgw(wn,nbin1) = nmgw(wn,nbin1) + 1
!        
          endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          
          if(nbin2 .le. maxbin) then
!        
         nznw(wn,nbin2) = nznw(wn,nbin2) + 1
!
          endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
          end do
!-------------------------------------------
     do wn1=1,so4_nnb
!
          xxmgso4=r_Q12A(1,wn)-r_S1A(1,wn1)
          xxmgso4=xxmgso4-CellL(1)*anint(xxmgso4/CellL(1))
          yymgso4=r_Q12A(2,wn)-r_S1A(2,wn1)
          yymgso4=yymgso4-CellL(2)*anint(yymgso4/CellL(2))
          zzmgso4=r_Q12A(3,wn)-r_S1A(3,wn1)
          zzmgso4=zzmgso4-CellL(3)*anint(zzmgso4/CellL(3))
          rrrmgso4_sq=(xxmgso4)**2+(yymgso4)**2+(zzmgso4)**2
          rMGSO4=sqrt(rrrmgso4_sq)
!
          xxznso4=r_T12A(1,wn)-r_S1A(1,wn1)
          xxznso4=xxznso4-CellL(1)*anint(xxznso4/CellL(1))
          yyznso4=r_T12A(2,wn)-r_S1A(2,wn1)
          yyznso4=yyznso4-CellL(2)*anint(yyznso4/CellL(2))
          zzznso4=r_T12A(3,wn)-r_S1A(3,wn1)
          zzznso4=zzznso4-CellL(3)*anint(zzznso4/CellL(3))
          rrrznso4_sq=(xxznso4)**2+(yyznso4)**2+(zzznso4)**2
          rZNSO4=sqrt(rrrznso4_sq)
!
          nbin3=anint(rMGSO4/drbin) + 1
          nbin4=anint(rZNSO4/drbin) + 1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          if(nbin3 .le. maxbin) then
!        
         nmgS(wn,nbin3) = nmgS(wn,nbin3) + 1
!
          endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          
          if(nbin4 .le. maxbin) then
!        
         nznS(wn,nbin4) = nznS(wn,nbin4) + 1
!!!!!!!!!    
          endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          
          end do
          end do
!!!!!!!-----------------------------------
           do 150 kjj=1,iaa
          amg(kjj)=0 
          azn(kjj)=0
          do 140 nmkl=1,maxbin
!          write(21,*) n, nci(n), ncw(n)
           ntotmgw(kjj)=0
           ntotmgS(kjj)=0
           ntotznw(kjj)=0
           ntotznS(kjj)=0           
!           
           nsum=nmkl

         do 162 iikk=1,nsum
         ntotmgw(kjj)=ntotmgw(kjj)+ nmgw(kjj,iikk)
         !
         ntotmgS(kjj)=ntotmgS(kjj) + nmgS(kjj,iikk)
         !
         ntotznw(kjj)=ntotznw(kjj)+ nznw(kjj,iikk)
         !
         ntotznS(kjj)=ntotznS(kjj) + nznS(kjj,iikk)
         !        
 162     continue
! 
        dnfacmg=(so4_nnb-ntotmgS(kjj))
        dfacmg=(wt_nnb-ntotmgw(kjj))
        facnewmg(kjj)=dnfacmg/dfacmg
!
        dnfaczn=(so4_nnb-ntotznS(kjj))
        dfaczn=(wt_nnb-ntotznw(kjj))
        facnewzn(kjj)=dnfaczn/dfaczn
!
       amg(kjj)=ntotmgS(kjj) -(facnewmg(kjj)*ntotmgw(kjj))
       azn(kjj)=ntotznS(kjj) -(facnewzn(kjj)*ntotznw(kjj))
!       
       prfrnmg(kjj,nmkl)=prfrnmg(kjj,nmkl)+ amg(kjj)
       prfrnzn(kjj,nmkl)=prfrnzn(kjj,nmkl)+ azn(kjj)
!
!       avgnion(k,n)=avgnion(k,n)+nci(k,n)
!       avgnh2o(k,n)=avgnh2o(k,n)+ncw(k,n)
!       write(3,*) 'see', n, aaa, prfrn(n)
 140   continue
 150   continue

!!!!-------------------------------------!!!!!!!!!!!!!
!         
        !  end do
!----------------------------------------    
         end do ! frame loop

       do 151 kk1=1,iaa
       do 139 np=1,maxbin
       !
       !avgnion(kk,np)= avgnion(kk,np)/nframes
       !avgnh2o(kk,np)=avgnh2o(kk,np)/nframes

       prfrnmg(kk1,np)=prfrnmg(kk1,np)/Step
       prfrnzn(kk1,np)=prfrnzn(kk1,np)/Step
       !
       dd=(np-1)*drbin
       write(3,*) dd, prfrnmg(kk1,np),prfrnzn(kk1,np)
       !write(2,*) dd, avgnion(kk,np)
       !write(8,*) dd,avgnh2o(kk,np)
 139   continue
 151   continue
       addmg=0.0
       addzn=0.0
       mmmm2=0
       do 171 ij1=1,maxbin
       add1=(ij1-1)*drbin
       apmg(ij1)=0.0
       apzn(ij1)=0.0
!
       do 181 ml=1,iaa
       apmg(ij1)=apmg(ij1)+prfrnmg(ml,ij1)
       apzn(ij1)=apzn(ij1)+prfrnzn(ml,ij1)
!
 181   continue
       apmg(ij1)=apmg(ij1)/iaa
       apzn(ij1)=apzn(ij1)/iaa
       write(32,*) add1,apmg(ij1),apzn(ij1)
       if ((add1.ge.10.0).and.(add1.le.14.0)) then
       addmg=addmg+apmg(ij1)
       addzn=addzn+apzn(ij1)
       mmmm2=mmmm2+1
       !
       end if
       write(33,*) addmg,addzn,apmg(ij1),apzn(ij1)
       !
 171   continue
       write(33,*)"(addmg/mmmm2),(addzn/mmmm2),mmmm2"
       write(33,*) (addmg/mmmm2),(addzn/mmmm2),mmmm2

     call f77_molfile_close_read_(Handle(Nfile),status)

   end do   

end
