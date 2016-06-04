      PROGRAM NANOBINARYREAD

!---------------------------NanoBinaryRead v1 -----------------------------
! purpose: 
! to use subroutine READNF to read data from near-field files written 
! by DDSCAT and then
! postprocess as desired for visualization, etc.
!
! allocatable arrays are passed through modules READNF_ECOM and READNF_BCOM
! other information is passed through READNF argument list
!
! This is like a lite version of DDPOSTPROCESS.
! It just grabs/processes basic data from the E/Bfield binary.
! If no B-field is requested, those fields are filled with zeros.
! We also grab the composition data at the end.
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Added - NanoBioNode, University of Illinois at Urbana-Champaign, 2016
! Note that this file has module dependencies from Flatau & Draine's DDSCAT: 
! readnf.f90,readnf_ecom.f90,readnf_bcom.f90,ddprecision.f90
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      USE DDPRECISION,ONLY: WP
      USE READNF_ECOM,ONLY: CXADIA,CXEINC,CXEPS,CXESCA,CXPOL,ICOMP
      USE READNF_BCOM,ONLY: CXBINC,CXBSCA

      IMPLICIT NONE

      CHARACTER :: CFLENAME*60,CFLPAR*60,CFLPAR_DEFAULT*60,COMMAND(10)*60
      CHARACTER :: CSTAMP*26, strafg*80

      INTEGER ::                                                        &
         IDVOUT,ILINE,IOBIN,IVTR,IX1,IX2,IY1,IY2,IZ1,IZ2,               &
         JA,JX,JY,JZ,K,                                                 &
         NAB,NAT0,NAT3,NCOMP,NRFLDB,NRWORD,NRWORD_NF,NX,NXY,NXYZ,NY,NZ, &
         VERSNUM                                                        !

      REAL(WP) ::                                           &
         AEFF,DPHYS,E2,EINC2,                               &
         NAMBIENT,PI,SNORM,SUMERR2,TINY,                    &
         W1,W2,W3,W4,W5,W6,W7,W8,WAVE,WX,WY,WZ,             &
         XA,XB,XMAX,XMIN,YA,YB,YMAX,YMIN,ZA,ZB,ZETA,ZMAX,ZMIN

      REAL(WP) ::    &
         AK_TF(1:3), &
         SVEC(1:3),  &
         S_INC(1:3), &
         X0(1:3),    &
         XTF(1:3)    !

      REAL(WP),ALLOCATABLE :: &
         S(:,:,:,:)

      COMPLEX(WP) :: CXBX,CXBY,CXBZ,CXERR,CXEX,CXEY,CXEZ
      COMPLEX(WP) ::   &
         CXB(1:3),     &
         CXB0_TF(1:3), &
         CXB_INC(1:3), &
         CXB_SCA(1:3), &
         CXE(1:3),     &
         CXE_INC(1:3), &
         CXE_SCA(1:3), &
         CXE0_TF(1:3), &
         CXP(1:3)      !
 
!VTR note that VTR graphics fields have to be in double precision
      REAL(KIND=8),ALLOCATABLE,DIMENSION(:) :: VTRX,VTRY,VTRZ        ! mesh
      REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:) :: VTR8,VTR8B, VTR8C ! data
      REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:) :: VTRVX,VTRVY,VTRVZ ! data
!VTR supplementary file of problem geometry
      REAL(KIND=8),DIMENSION(1):: XVECT,YVECT,ZVECT
      REAL(KIND=8),DIMENSION(1,1,1):: U,V,W
      COMPLEX(WP):: CXBB(3),CXEE(3)
      INTEGER::IX,IY,IZ

      if (iargc() .ne. 2) then
         print *, 'Usage: myddpostprocess INPUTname SquareOn'
         call exit(1)
      else
         call getarg(1, CFLENAME)
         call getarg(2, strafg)
         read (unit=strafg, fmt=*) IVTR
      endif




!==============================================================================

      CALL READNF(CFLENAME,IDVOUT,CSTAMP,VERSNUM,NRFLDB,AEFF,DPHYS, &
                  NX,NY,NZ,NAT0,X0,XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX,   &
                  NAMBIENT,WAVE,AK_TF,CXE0_TF,CXB0_TF,NCOMP)        !

!==============================================================================
      S_INC(1)=REAL(CXE0_TF(2)*CONJG(CXB0_TF(3))- &
                    CXE0_TF(3)*CONJG(CXB0_TF(2)))
      S_INC(2)=REAL(CXE0_TF(3)*CONJG(CXB0_TF(1))- &
                    CXE0_TF(1)*CONJG(CXB0_TF(3)))
      S_INC(3)=REAL(CXE0_TF(1)*CONJG(CXB0_TF(2))- &
                    CXE0_TF(2)*CONJG(CXB0_TF(1)))
      SNORM=SQRT(S_INC(1)**2+S_INC(2)**2+S_INC(3)**2)
      ALLOCATE(S(1:NX,1:NY,1:NZ,1:3))
      IF(NRFLDB==1)THEN
         DO JZ=1,NZ
            DO JY=1,NY
               DO JX=1,NX
                  CXEX=CXEINC(JX,JY,JZ,1)+CXESCA(JX,JY,JZ,1)
                  CXEY=CXEINC(JX,JY,JZ,2)+CXESCA(JX,JY,JZ,2)
                  CXEZ=CXEINC(JX,JY,JZ,3)+CXESCA(JX,JY,JZ,3)
                  CXBX=CXBINC(JX,JY,JZ,1)+CXBSCA(JX,JY,JZ,1)
                  CXBY=CXBINC(JX,JY,JZ,2)+CXBSCA(JX,JY,JZ,2)
                  CXBZ=CXBINC(JX,JY,JZ,3)+CXBSCA(JX,JY,JZ,3)
                  S(JX,JY,JZ,1)=REAL(CXEY*CONJG(CXBZ)-CXEZ*CONJG(CXBY))/SNORM
                  S(JX,JY,JZ,2)=REAL(CXEZ*CONJG(CXBX)-CXEX*CONJG(CXBZ))/SNORM
                  S(JX,JY,JZ,3)=REAL(CXEX*CONJG(CXBY)-CXEY*CONJG(CXBX))/SNORM
               ENDDO
            ENDDO
         ENDDO
      ENDIF


      ALLOCATE(VTRX(NX),VTRY(NY),VTRZ(NZ))
      ALLOCATE(VTR8(NX,NY,NZ))
      ALLOCATE(VTR8B(NX,NY,NZ))
      ALLOCATE(VTR8C(NX,NY,NZ))

      IF(NRFLDB==1)THEN
         ALLOCATE(VTRVX(NX,NY,NZ))
         ALLOCATE(VTRVY(NX,NY,NZ))
         ALLOCATE(VTRVZ(NZ,NY,NZ))
      ENDIF
       
! calculate x,y,z in dimensional units

      DO JX=1,NX
         VTRX(JX)=(JX+X0(1))*DPHYS
      ENDDO
      DO JY=1,NY
         VTRY(JY)=(JY+X0(2))*DPHYS
      ENDDO
      DO JZ=1,NZ
         VTRZ(JZ)=(JZ+X0(3))*DPHYS
      ENDDO

      open (UNIT=12, FILE='field_datafile.txt', STATUS='UNKNOWN')      
      write (12, FMT=*) "Dimensions of system (x,y,z): ", NX, NY, NZ
      write (12, FMT=*)										  			&
      "          Xcoord      ","     Ycoord      ","     Zcoord      ",	&
      " EField      ",													&
      "     EField-X(Re)   ","  EField-X(Im)    ",			&
      "   EField-Y(Re)","     EField-Y(Im)    ",			&
      "   EField-Z(Re)  ","   EField-Z(Im)    ",			&
       "B On/Off      ", " Bfield         ",			&
      "  BField-X(Re)   ","  BField-X(Im)     ",			&
      "  BField-Y(Re)   ","  BField-Y(Im)     ",			&
      "  BField-Z(Re)   ","  BField-Z(Im) ",			&    
      "  Poynting X      ",									&
      "     Poynting Y      ",									&
      " Poynting Z      "									!
      DO JZ=1,NZ
         DO JY=1,NY
            DO JX=1,NX
               CXEE(1:3) = CXEINC(JX,JY,JZ,1:3)+CXESCA(JX,JY,JZ,1:3)
               IF (IVTR==2) THEN
                  VTR8(JX,JY,JZ)=DOT_PRODUCT(CXEE(1:3),CXEE(1:3))
               ELSE
                  VTR8(JX,JY,JZ)=SQRT(DOT_PRODUCT(CXEE(1:3),CXEE(1:3)))
               ENDIF
               IF(NRFLDB==1)THEN
                  CXBB(1:3)=CXBINC(JX,JY,JZ,1:3)+CXBSCA(JX,JY,JZ,1:3)
                  IF (IVTR==2) THEN
                     VTR8B(JX,JY,JZ)=DOT_PRODUCT(CXBB(1:3),CXBB(1:3))         
                  ELSE
                     VTR8B(JX,JY,JZ)=SQRT(DOT_PRODUCT(CXBB(1:3),CXBB(1:3)))
                  ENDIF
               ENDIF
               
               VTR8C(JX,JY,JZ) = ICOMP(JX,JY,JZ,1)

               write (12, FMT=*)										  &
               real(VTRX(JX),kind=selected_real_kind(6)), 				  &
               real(VTRY(JY),kind=selected_real_kind(6)), 				  &
               real(VTRZ(JZ),kind=selected_real_kind(6)),				  &
               real(VTR8(JX,JY,JZ),kind=selected_real_kind(6)),			  &
               CXEE(1:3),												  &
               NRFLDB,													  &
               real(VTR8B(JX,JY,JZ),kind=selected_real_kind(6)),		  &
               CXBB(1:3),												  &
               S(JX,JY,JZ,1:3)											  !
            ENDDO
         ENDDO
      ENDDO
      close(UNIT=12)

      
      !Write composition data to a separate file for easy handling in nanoDDSCAT
      open (UNIT=12, FILE='composition.txt', STATUS='UNKNOWN')
      write(unit=12,fmt=*) real(vtr8c(1:nx,1:ny,1:nz),kind=selected_real_kind(6))	
      close(UNIT=12)
 
      STOP
    END PROGRAM NANOBINARYREAD
