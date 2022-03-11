!------------------------------------------------------------------------------------
!
!      FILE mod_storm.F
!
!      This file is part of the FUNWAVE-TVD program under the Simplified BSD license
!
!-------------------------------------------------------------------------------------
! 
!    Copyright (c) 2016, FUNWAVE Development Team
!
!    (See http://www.udel.edu/kirby/programs/funwave/funwave.html
!     for Development Team membership)
!
!    All rights reserved.
!
!    FUNWAVE_TVD is free software: you can redistribute it and/or modify
!    it under the terms of the Simplified BSD License as released by
!    the Berkeley Software Distribution (BSD).
!
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions are met:
!
!    1. Redistributions of source code must retain the above copyright notice, this
!       list of conditions and the following disclaimer.
!    2. Redistributions in binary form must reproduce the above copyright notice,
!    this list of conditions and the following disclaimer in the documentation
!    and/or other materials provided with the distribution.
!
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
!    ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
!    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
!    ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
!    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
!    ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
!    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!  
!    The views and conclusions contained in the software and documentation are those
!    of the authors and should not be interpreted as representing official policies,
!    either expressed or implied, of the FreeBSD Project.
!  
!-------------------------------------------------------------------------------------
!
!  METEO is a module to model wind and pressure effects   
!  it can also provide sources for meteo tsunami and landslide 
!
!  HISTORY :
!    05/16/2017  Fengyan Shi
!
!-------------------------------------------------------------------------------------

# if defined (METEO)

MODULE METEO_MODULE
  USE PARAM
  USE GLOBAL,ONLY : Mloc,Nloc,Nghost,Ibeg,Iend,Jbeg,Jend,DX,DY, &
                    H,ETA,ETAx,ETAy,ETAt,HeightMax,MinDepthFrc,ETA0
  USE INPUT_READ
#if defined (PARALLEL)
  USE GLOBAL,ONLY : myid,ier, npx,npy,PX,PY
  USE MPI
# endif
  IMPLICIT NONE
  SAVE

       LOGICAL :: MeteoGausian = .FALSE.
       LOGICAL :: AxisymmetricGaussian = .FALSE.
	   LOGICAL :: LPGaussian = .FALSE.
       LOGICAL :: SlideModel = .FALSE.
       LOGICAL :: WindConstantField = .FALSE.
       LOGICAL :: WindForce = .FALSE.
       LOGICAL :: WindHollandModel = .FALSE.
       LOGICAL :: AirPressure = .FALSE.
       LOGICAL :: WindWaveInteraction = .FALSE.
       CHARACTER(LEN=80) WIND_FILE
       INTEGER :: NumTimeWindData
       REAL(SP), DIMENSION(:,:),ALLOCATABLE :: WindU2D,WindV2D
       REAL(SP),DIMENSION(:),ALLOCATABLE :: TimeWind
       REAL(SP),DIMENSION(:),ALLOCATABLE :: WU,WV
       REAL(SP) :: Cdw,WindCrestPercent
       INTEGER :: icount_winddata = 1   
       INTEGER,DIMENSION(:,:),ALLOCATABLE :: MASK_WIND 

    REAL(SP) :: Xstorm1,Ystorm1,Xstorm2,Ystorm2, &
                Pn_storm1, Pc_storm1,A_storm1,B_storm1, &
                Pn_storm2, Pc_storm2,A_storm2,B_storm2, &
                TimeStorm1,TimeStorm2,ThetaStorm

!   Gausian
    REAL(SP) :: DP_storm2,SigmaX2, SigmaY2, Theta2, &
                DP_storm1,SigmaX1, SigmaY1, Theta1
				
!   Axisymmetric Gaussian
    REAL(SP) :: Radius_storm1,PT_storm1,LC_storm1,LT_storm1,& !PC_storm1 is duplicated
				TP_storm1,PSp_storm1,ExpG_storm1,Decay_storm1,CDis_storm1, &
				Lp_storm,RLp_storm,TLp_storm,Kap_storm,RA_storm
				
!   LP Gaussian
    REAL(SP) :: PSpC_storm1,PSpT_storm1,PAC_storm1,PAT_storm1

    REAL(SP),DIMENSION(:,:),ALLOCATABLE :: StormPressureTotal, &
                                           StormPressureX,StormPressureY
    REAL(SP),DIMENSION(:),ALLOCATABLE ::  Xco,Yco
    LOGICAL :: OUT_METEO = .TRUE.
!    REAL(SP):: PLOT_INTV_STORM,PLOT_COUNT_STORM
! slide
    REAL(SP) :: LengthSlide,WidthSlide,AlphaSlide,BetaSlide,PSlide,epsilon
    LOGICAL :: FirstCall_Met = .TRUE.                


#if defined (PARALLEL)
    REAL(SP) :: myvar
# endif    

CONTAINS
  
! READ STORM

SUBROUTINE METEO_INITIAL
  USE GLOBAL,ONLY : itmp1,itmp2,itmp3,itmp4,SMALL,LARGE,INPUT_FILE_NAME
                    
  USE INPUT_READ
  IMPLICIT NONE
  CHARACTER(LEN=80)::FILE_NAME=' '
  CHARACTER(LEN=80)::TMP_NAME=' '
  INTEGER :: Ifile,ierr
  CHARACTER(LEN=80):: StormName
  CHARACTER(LEN=80) :: WHAT

       ALLOCATE (StormPressureTotal(Mloc,Nloc),&
                StormPressureX(Mloc,Nloc), &
                StormPressureY(Mloc,Nloc))
         StormPressureTotal=ZERO
         StormPressureX=ZERO
         StormPressureY=ZERO

       ALLOCATE(WindU2D(Mloc,Nloc),WindV2D(Mloc,Nloc),MASK_WIND(Mloc,Nloc))
       WindU2D=ZERO
       WindV2D=ZERO
       MASK_WIND=0


! read storm from input.txt
      FILE_NAME=INPUT_FILE_NAME

# if defined (PARALLEL)
      if (myid.eq.0) WRITE(3,*)'                                         '
      if (myid.eq.0) WRITE(3,*)'-------------- METEO INFO -----------------'
# else
      WRITE(3,*)'                                         '
      WRITE(3,*)'-------------- METEO INFO -----------------'   
# endif

! type of meteo
  ! constant wind field
      CALL READ_LOGICAL(WindConstantField,FILE_NAME,'WindConstantField',ierr)
      IF(ierr==1)THEN
        WindConstantField = .FALSE.
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'not constant wind field'
         WRITE(3,'(A40)')'not constant wind field'
      endif
# else
         WRITE(*,'(A40)')'not constant wind field'
         WRITE(3,'(A40)')'not constant wind field'
# endif
       ENDIF

  ! holland wind field
      CALL READ_LOGICAL(WindHollandModel,FILE_NAME,'WindHollandModel',ierr)
      IF(ierr==1)THEN
        WindHollandModel = .FALSE.
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'not WindHollandModel'
         WRITE(3,'(A40)')'not WindHollandModel'
      endif
# else
         WRITE(*,'(A40)')'not WindHollandModel'
         WRITE(3,'(A40)')'not WindHollandModel'
# endif
       ENDIF

  ! meteo gausian 2D
      CALL READ_LOGICAL(MeteoGausian,FILE_NAME,'MeteoGausian',ierr)
      IF(ierr==1)THEN
        MeteoGausian = .FALSE.
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'not use MeteoGausian'
         WRITE(3,'(A40)')'not use MeteoGausian'
      endif
# else
         WRITE(*,'(A40)')'not use MeteoGausian'
         WRITE(3,'(A40)')'not use MeteoGausian'
# endif
       ENDIF
	   
  ! axisymmetric gaussian 2D
      CALL READ_LOGICAL(AxisymmetricGaussian,FILE_NAME,'AxisymmetricGaussian',ierr)
      IF(ierr==1)THEN
        AxisymmetricGaussian = .FALSE.
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'not use AxisymmetricGaussian'
         WRITE(3,'(A40)')'not use AxisymmetricGaussian'
      endif
# else
         WRITE(*,'(A40)')'not use AxisymmetricGaussian'
         WRITE(3,'(A40)')'not use AxisymmetricGaussian'
# endif
       ENDIF
	   
  ! Lynett pressure gaussian 2D
      CALL READ_LOGICAL(LPGaussian,FILE_NAME,'LPGaussian',ierr)
      IF(ierr==1)THEN
        LPGaussian = .FALSE.
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'not use LPGaussian'
         WRITE(3,'(A40)')'not use LPGaussian'
      endif
# else
         WRITE(*,'(A40)')'not use LPGaussian'
         WRITE(3,'(A40)')'not use LPGaussian'
# endif
       ENDIF

  ! slide model
      CALL READ_LOGICAL(SlideModel,FILE_NAME,'SlideModel',ierr)
      IF(ierr==1)THEN
        SlideModel = .FALSE.
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'not SlideModel'
         WRITE(3,'(A40)')'not SlideModel'
      endif
# else
         WRITE(*,'(A40)')'not SlideModel'
         WRITE(3,'(A40)')'not SlideModel'
# endif
       ENDIF

! ------------------ meteo gausian
  IF(MeteoGausian)THEN
         AirPressure = .TRUE.
  ENDIF
  
! ------------------ axisymmetric gaussian
  IF(AxisymmetricGaussian)THEN
         AirPressure = .TRUE.
  ENDIF
  
! ------------------ Lynett pressure gaussian
  IF(LPGaussian)THEN
         AirPressure = .TRUE.
  ENDIF

! -------------------wind and pressure---------------
  IF(WindHollandModel.OR.WindConstantField)THEN

    CALL READ_LOGICAL(WindWaveInteraction,  &
     FILE_NAME,'WindWaveInteraction',ierr)
      IF(ierr==1)THEN
        WindWaveInteraction = .FALSE.
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'no WindWaveInteraction'
         WRITE(3,'(A40)')'no WindWaveInteraction'
      endif
# else
         WRITE(*,'(A40)')'no WindWaveInteraction'
         WRITE(3,'(A40)')'no WindWaveInteraction'
# endif
       ENDIF

! air pressure
      CALL READ_LOGICAL(AirPressure,FILE_NAME,'AirPressure',ierr)

      IF(ierr==1)THEN
        AirPressure = .FALSE.
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'NO air pressure effect'
         WRITE(3,'(A40)')'NO air pressure effect'
      endif
# else
         WRITE(*,'(A40)')'NO air pressure effect'
         WRITE(3,'(A40)')'NO air pressure effect'
# endif
       ENDIF

! wind
      CALL READ_LOGICAL(WindForce,FILE_NAME,'WindForce',ierr)

      IF(ierr==1)THEN
       IF(WindConstantField)THEN
        WindForce = .TRUE.
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'force set windforce=true because use windconstantfield'
         WRITE(3,'(A80)')'force set windforce=true because use windconstantfield'
      endif
# else
         WRITE(*,'(A80)')'force set windforce=true because use windconstantfield'
         WRITE(3,'(A80)')'force set windforce=true because use windconstantfield'
# endif
       ELSE
        WindForce = .FALSE.
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'NO wind effect'
         WRITE(3,'(A40)')'NO wind effect'
      endif
# else
         WRITE(*,'(A40)')'NO wind effect'
         WRITE(3,'(A40)')'NO wind effect'
# endif
        ENDIF ! end windconstantfield
      ENDIF ! end ierr


      IF(WindForce)THEN

      CALL READ_FLOAT(Cdw,FILE_NAME,'Cdw',ierr)
      IF(ierr==1)THEN
        Cdw=0.002_SP
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'Cdw Default:  Cdw=0.002'
         WRITE(3,'(A40)')'Cdw Default:  Cdw=0.002'
      endif
# else
         WRITE(*,'(A40)')'Cdw Default:  Cdw=0.002'
         WRITE(3,'(A40)')'Cdw Default:  Cdw=0.002'
# endif
       ENDIF  ! end ierr default

        CALL READ_FLOAT(WindCrestPercent,FILE_NAME,'WindCrestPercent',ierr)

      IF(ierr==1)THEN
        WindCrestPercent=LARGE
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'WindCrestPercent Default:  LARGE'
         WRITE(3,'(A40)')'WindCrestPercent Default:  LARGE'
      endif
# else
         WRITE(*,'(A40)')'WindCrestPercent Default:  LARGE'
         WRITE(3,'(A40)')'WindCrestPercent Default:  LARGE'
# endif
       ELSE
         IF(.NOT.WindWaveInteraction)THEN
            WindCrestPercent=LARGE
         ENDIF

       ENDIF ! end ierr windcrestpercent default

# if defined (PARALLEL)
      if (myid.eq.0) WRITE(3,'(A7,F12.5)')'Cdw = ', Cdw
      if (myid.eq.0) WRITE(3,'(A21,F12.5)')'WindCrestPercent = ', WindCrestPercent
# else
      WRITE(3,'(A7,F12.5)')'Cdw = ', Cdw
      WRITE(3,'(A21,F12.5)')'WindCrestPercent = ', WindCrestPercent
# endif 

      ENDIF ! end windforce

  ENDIF ! end wind and pressure
! ----------------------------------------------------

  IF(WindHollandModel) THEN

     CALL Holland_Model_Setup

  ENDIF ! end holland model

  IF(MeteoGausian) THEN

     CALL MeteoGausian_Setup

  ENDIF ! end meteogausian
  
  IF(AxisymmetricGaussian) THEN

     CALL AxisymmetricGaussian_Setup

  ENDIF ! end AxisymmetricGaussian
  
  IF(LPGaussian) THEN

     CALL LPGaussian_Setup

  ENDIF ! end LPGaussian

  IF (WindConstantField) THEN

     CALL Contant_Wind_Field_Setup

  ENDIF ! end constant wind field

     IF (SlideModel) THEN
       CALL Slide_Model_Setup
     ENDIF

  CALL READ_LOGICAL(OUT_METEO,FILE_NAME,'OUT_METEO',ierr)

End SUBROUTINE METEO_INITIAL

SUBROUTINE METEO_FORCING
  USE GLOBAL,ONLY : Mloc,Nloc,tmp1,tmp2,SMALL,TIME,ZERO
  USE INPUT_READ
  IMPLICIT NONE
  INTEGER :: Ifile,ierr,I,J

  IF (WindHollandModel) THEN

    CALL Holland_Model_Forcing

  ENDIF ! end holland model

  IF (MeteoGausian) THEN

    CALL MeteoGausian_Forcing

  ENDIF ! end meteogausian
  
  IF (AxisymmetricGaussian) THEN

    CALL AxisymmetricGaussian_Forcing

  ENDIF ! end AxisymmetricGaussian
  
  IF (LPGaussian) THEN

    CALL LPGaussian_Forcing

  ENDIF ! end LPGaussian

  IF(WindConstantField)THEN

    CALL Constant_Wind_Forcing

  ENDIF ! end constand wind field

  IF(SlideModel)THEN
    CALL Slide_Model_Forcing
  ENDIF

END SUBROUTINE METEO_FORCING

SUBROUTINE Holland_Model_Setup
  USE GLOBAL,ONLY : itmp1,itmp2,itmp3,itmp4,SMALL,LARGE, INPUT_FILE_NAME
# if defined (PARALLEL)
  USE GLOBAL,ONLY :   iista,jjsta  !ykchoi Jan/23/2018
# endif
  USE INPUT_READ
  IMPLICIT NONE
  CHARACTER(LEN=80)::FILE_NAME=' '
  CHARACTER(LEN=80)::STORM_FILE =' '
  CHARACTER(LEN=80)::TMP_NAME=' '
  INTEGER :: Ifile,ierr
  CHARACTER(LEN=80):: StormName

! read storm from input.txt
      FILE_NAME=INPUT_FILE_NAME

! storm file
      CALL READ_STRING(STORM_FILE,FILE_NAME,'STORM_FILE',ierr)
# if defined (PARALLEL)
      if (myid.eq.0) WRITE(3,'(A15,A50)')'STORM_FILE:', STORM_FILE
# else
      WRITE(3,'(A15,A50)')'STORM_FILE:', STORM_FILE
# endif


      ALLOCATE (Xco(Mloc),Yco(Nloc))

      
! Xco, and Yco

# if defined (CARTESIAN)
! Cartesian
# if defined (PARALLEL)
![---ykchoi Jan/23/2018
!      Xco(Ibeg) = npx*(Mloc-2*Nghost)*DX
	Xco(Ibeg) = (iista-1)*DX
!---ykchoi Jan/23/2018]

# else
     Xco(Ibeg) = ZERO
# endif
     DO I = Ibeg+1,Mloc
       Xco(I) = Xco(I-1)+DX
     ENDDO
     DO I = Ibeg-1,Ibeg-Nghost,-1
       Xco(I) = Xco(I+1)-DX
     ENDDO

# if defined (PARALLEL)
![---ykchoi Jan/23/2018
!     Yco(Jbeg) = npy*(Nloc-2*Nghost)*DY
      Yco(Jbeg) = (jjsta-1)*DY
!---ykchoi Jan/23/2018]
# else
     Yco(Jbeg) = ZERO
# endif
     DO J = Jbeg+1,Nloc
       Yco(J) = Yco(J-1)+DY
     ENDDO
     DO J = Jbeg-1,Jbeg-Nghost,-1
       Yco(J) = Yco(J+1)-DY
     ENDDO
# endif

    TMP_NAME = TRIM(STORM_FILE)

! check existing

 INQUIRE(FILE=TRIM(TMP_NAME),EXIST=FILE_EXIST)
  IF(.NOT.FILE_EXIST)THEN
# if defined (PARALLEL)
   IF(MYID==0)  &
   WRITE(*,*) TRIM(TMP_NAME), 'CANNOT BE FOUND. STOP'
   CALL MPI_FINALIZE (ier)
   STOP
# else
    WRITE(*,*) TRIM(TMP_NAME), 'CANNOT BE FOUND. STOP'
    STOP
# endif
  ENDIF

! open file
  Ifile=300
  OPEN(Ifile,FILE=TRIM(TMP_NAME))

! read file
         READ(Ifile,*)  ! title
         READ(Ifile,*)  StormName  !  name
         READ(Ifile,*)  ! t,x,y, Pn, Pc, A and B
         READ(Ifile,*)  TimeStorm2,Xstorm2,Ystorm2,  &
                        Pn_storm2,Pc_storm2, A_storm2, B_storm2

         TimeStorm1 = TimeStorm2
         Xstorm1 = Xstorm2
         Ystorm1 = Ystorm2
         Pn_storm1 = Pn_storm2
         Pc_storm1 = Pc_storm2
         A_storm1 = A_storm2
         B_storm1 = B_storm2

# if defined (PARALLEL)
   IF(MYID==0)THEN
   WRITE(3,*) 'Storm Name: ',  TRIM(StormName)
   WRITE(3,*) 'Initial Time, X, Y', TimeStorm2,Xstorm2,Ystorm2
   WRITE(3,*) 'Storm Pn,Pc,A, B: ',Pn_storm2,Pc_storm2, A_storm2, B_storm2
   ENDIF
# else
   WRITE(*,*) 'Storm Name: ',  TRIM(StormName)
   WRITE(*,*) 'Initial Time, X, Y', TimeStorm2,Xstorm2,Ystorm2
   WRITE(*,*) 'Storm Pn,Pc,A, B: ',Pn_storm2,Pc_storm2, A_storm2, B_storm2
   WRITE(3,*) 'Storm Name: ',  TRIM(StormName)
   WRITE(3,*) 'Initial Time, X, Y', TimeStorm2,Xstorm2,Ystorm2
   WRITE(3,*) 'Storm Pn,Pc,A, B: ',Pn_storm2,Pc_storm2, A_storm2, B_storm2
# endif

END SUBROUTINE Holland_Model_Setup

SUBROUTINE Contant_Wind_Field_Setup
  USE GLOBAL,ONLY : itmp1,itmp2,itmp3,itmp4,SMALL,LARGE,INPUT_FILE_NAME
  USE INPUT_READ
  IMPLICIT NONE
  CHARACTER(LEN=80)::FILE_NAME=' '
  CHARACTER(LEN=80)::STORM_FILE =' '
  CHARACTER(LEN=80)::TMP_NAME=' '
  INTEGER :: Ifile,ierr
  CHARACTER(LEN=80):: StormName
  CHARACTER(LEN=80) :: WHAT

! read storm from input.txt
      FILE_NAME=INPUT_FILE_NAME

  CALL READ_STRING(WIND_FILE,FILE_NAME,'CONSTANT_WIND_FILE',ierr)

  INQUIRE(FILE=TRIM(WIND_FILE),EXIST=FILE_EXIST)
    IF(.NOT.FILE_EXIST)THEN
# if defined (PARALLEL)
     IF(MYID==0) &
      WRITE(*,*) TRIM(WIND_FILE), ' specified in input.txt but does not exist,stop'
     CALL MPI_FINALIZE (ier)
     STOP
# else
      WRITE(*,*) TRIM(WIND_FILE), ' specified in input.txt but does not exist,stop'
      STOP
# endif
    ENDIF

       OPEN(2,FILE=TRIM(WIND_FILE))
         READ(2,*)WHAT
         READ(2,*)NumTimeWindData
         ALLOCATE (TimeWind(NumTimeWindData),WU(NumTimeWindData),&
                WV(NumTimeWindData))
         DO I=1,NumTimeWindData
           READ(2,*,END=111)TimeWind(I),WU(I),WV(I)
         ENDDO
111    CONTINUE
       CLOSE(2)

END SUBROUTINE Contant_Wind_Field_Setup

SUBROUTINE Holland_Model_Forcing
  USE GLOBAL,ONLY : Mloc,Nloc,tmp1,tmp2,SMALL,TIME,ZERO,ETAmean
  USE INPUT_READ
  IMPLICIT NONE
  INTEGER :: Ifile,ierr,I,J
  REAL(SP) :: Pn_storm, Pc_storm,A_storm,B_storm,Pw,Vw,Rdis,AngleLocal
  REAL(SP) :: Xstorm, Ystorm
  REAL(SP) :: Celerity,WaveAngle

    StormPressureTotal = ZERO

    IF(TIME>TimeStorm1.AND.TIME>TimeStorm2) THEN

         TimeStorm1=TimeStorm2
         Xstorm1 = Xstorm2
         Ystorm1 = Ystorm2

    Ifile = 300

    READ(Ifile,*,END=120)  TimeStorm2,Xstorm2,Ystorm2,  &
                        Pn_storm2,Pc_storm2, A_storm2, B_storm2

# if defined (PARALLEL)
   IF(MYID==0)THEN
     WRITE(3,*)'T,X,Y = ', TimeStorm2,Xstorm2,Ystorm2
     WRITE(*,*)'T,X,Y = ', TimeStorm2,Xstorm2,Ystorm2
   WRITE(3,*) 'Storm Pn,Pc,A, B: ',Pn_storm2,Pc_storm2, A_storm2, B_storm2
   WRITE(*,*) 'Storm Pn,Pc,A, B: ',Pn_storm2,Pc_storm2, A_storm2, B_storm2
   ENDIF
# else
   WRITE(*,*)'T,X,Y = ', TimeStorm2,Xstorm2,Ystorm2
   WRITE(*,*) 'Storm Pn,Pc,A, B: ',Pn_storm2,Pc_storm2, A_storm2, B_storm2
   WRITE(3,*)'T,X,Y = ', TimeStorm2,Xstorm2,Ystorm2
   WRITE(3,*) 'Storm Pn,Pc,A, B: ',Pn_storm2,Pc_storm2, A_storm2, B_storm2
# endif

    ThetaStorm = ATAN2(Ystorm2-Ystorm1,  &
                              Xstorm2-Xstorm1)


    ENDIF ! end time > timestorm2

! calculate force
    tmp2=ZERO
    tmp1=ZERO

    IF(TIME>TimeStorm1)THEN
      IF(TimeStorm1.EQ.TimeStorm2)THEN
        ! no more data
        tmp2=ZERO
        tmp1=ZERO
      ELSE
      tmp2=(TimeStorm2-TIME) &
            /MAX(SMALL, ABS(TimeStorm2-TimeStorm1))
      tmp1=1.0_SP - tmp2;
      ENDIF  ! no more data?
    ENDIF ! time>time_1

    Xstorm = Xstorm2*tmp1 +Xstorm1*tmp2
    Ystorm = Ystorm2*tmp1 +Ystorm1*tmp2
    Pn_storm = Pn_storm2*tmp1 + Pn_storm1*tmp2
    Pc_storm = Pc_storm2*tmp1 + Pc_storm1*tmp2
    A_storm = A_storm2*tmp1 + A_storm1*tmp2
    B_storm = B_storm2*tmp1 + B_storm1*tmp2

120 CONTINUE  ! no more data for vessel Kves

! sourceX and sourceY

    DO J=1,Nloc
    DO I=1,Mloc
       Rdis=SQRT((Xco(I)-Xstorm)**2+(Yco(J)-Ystorm)**2)/1000.0_SP
       Rdis=MAX(SMALL,Rdis)
          ! Rdis is in km, Pw is in mb
       Pw=Pc_storm +(Pn_storm-Pc_storm)*EXP(-A_storm/Rdis**B_storm)

      IF(AirPressure)THEN
       StormPressureTotal(I,J) = Pw/100.0_SP  ! convert from cm to meter
      ENDIF

       Vw=SQRT(A_storm*B_storm*100.0_SP*ABS(Pn_storm-Pc_storm)  &
           *EXP(-A_storm/Rdis**B_storm) / Rho_air/Rdis**B_storm)

       AngleLocal = ATAN2((Xco(I)-Xstorm),(Yco(J)-Ystorm))

    IF(WindForce)THEN
     IF(WindWaveInteraction)THEN
       tmp3=SQRT(GRAV*MAX(MinDepthFrc,H(I,J)))
       tmp1=MAX(SQRT(ETAx(I,J)*ETAx(I,J)+ETAy(I,J)*ETAy(I,J)),SMALL)
       Celerity=MIN(ABS(ETAt(I,J))/tmp1,SQRT(GRAV*ABS(H(I,J))))
       WaveAngle=ATAN2(ETAy(I,J),ETAx(I,J))

       WindU2D(I,J)=-Vw*COS(AngleLocal)-Celerity*COS(WaveAngle)
       WindV2D(I,J)=Vw*SIN(AngleLocal)-Celerity*SIN(WaveAngle)
     ELSE       
       WindU2D(I,J)=-Vw*COS(AngleLocal)
       WindV2D(I,J)=Vw*SIN(AngleLocal)
     ENDIF

    ENDIF ! end windforce

    ENDDO
    ENDDO


   IF(WindForce)THEN
    MASK_WIND=1
     DO J=1,Nloc
     DO I=1,Mloc
      IF(WindCrestPercent==LARGE)THEN
        ! keep mask_wind = 1
      ELSE
        IF(ETA(I,J)<HeightMax(I,J)*(1-WindCrestPercent)+ETAmean(I,J)) MASK_WIND(I,J)=0
      ENDIF
     ENDDO
     ENDDO
   ENDIF

   IF(AirPressure)THEN
    DO J=Jbeg,Jend
    DO I=Ibeg,Iend
# if defined (CARTESIAN)
       StormPressureX(I,J) = -Grav*H(I,J)*  &
               (StormPressureTotal(I+1,J)-StormPressureTotal(I-1,J))/2.0_SP  &
               /DX
       StormPressureY(I,J) = -Grav*H(I,J)*  &
               (StormPressureTotal(I,J+1)-StormPressureTotal(I,J-1))/2.0_SP  &
               /DY
# endif
    ENDDO
    ENDDO
   ENDIF


END SUBROUTINE Holland_Model_Forcing

SUBROUTINE Constant_Wind_Forcing
  USE GLOBAL,ONLY : Mloc,Nloc,tmp1,tmp2,SMALL,TIME,ZERO,ETAmean
  USE INPUT_READ
  IMPLICIT NONE
  INTEGER :: Ifile,ierr,I,J
  REAL(SP) :: AngleLocal
  REAL(SP) :: Celerity,WaveAngle

    IF(WindForce)THEN
! time series wind
      IF(TIME>TimeWind(icount_winddata).AND. &
               icount_winddata<NumTimeWindData)THEN
        icount_winddata=icount_winddata+1
      ENDIF
      IF(icount_winddata>1)THEN  ! wind start
       IF(TIME>TimeWind(icount_winddata))THEN
        tmp2=ZERO
       ELSE
        tmp2=(TimeWind(icount_winddata)-TIME) &
           /(TimeWind(icount_winddata)-TimeWind(icount_winddata-1))
       ENDIF
       WindU2D=WU(icount_winddata)*(1.0_SP-tmp2)&
                       +WU(icount_winddata-1)*tmp2
       WindV2D=WV(icount_winddata)*(1.0_SP-tmp2)&
                       +WV(icount_winddata-1)*tmp2
      ENDIF

! adjusted by wave celerity based on Chen et al. 2004


    IF(WindWaveInteraction)THEN
     DO J=1,Nloc
     DO I=1,Mloc
       tmp3=SQRT(GRAV*MAX(MinDepthFrc,H(I,J)))
       tmp1=MAX(SQRT(ETAx(I,J)*ETAx(I,J)+ETAy(I,J)*ETAy(I,J)),SMALL)
       Celerity=MIN(ABS(ETAt(I,J))/tmp1,SQRT(GRAV*ABS(H(I,J))))
       AngleLocal=ATAN2(ETAy(I,J),ETAx(I,J))

       WindU2D(I,J)=WindU2D(I,J)-Celerity*COS(AngleLocal)
       WindV2D(I,J)=WindV2D(I,J)-Celerity*SIN(AngleLocal)
     ENDDO
     ENDDO
    ENDIF

     MASK_WIND=1
     DO J=1,Nloc
     DO I=1,Mloc
      IF(WindCrestPercent==LARGE)THEN
        ! keep mask_wind = 1
      ELSE
        IF(ETA(I,J)<HeightMax(I,J)*(1-WindCrestPercent)+ETAmean(I,J))MASK_WIND(I,J)=0
      ENDIF
     ENDDO
     ENDDO

   ENDIF ! end windforce

END SUBROUTINE Constant_Wind_Forcing

SUBROUTINE Slide_Model_Setup
  USE GLOBAL,ONLY : itmp1,itmp2,itmp3,itmp4,SMALL,LARGE,INPUT_FILE_NAME
# if defined (PARALLEL)
  USE GLOBAL,ONLY : iista,jjsta    !ykchoi Jan/23/2018
# endif
  USE INPUT_READ
  IMPLICIT NONE
  CHARACTER(LEN=80)::FILE_NAME=' '
  CHARACTER(LEN=80)::SLIDE_FILE =' '
  CHARACTER(LEN=80)::TMP_NAME=' '
  INTEGER :: Ifile,ierr
  CHARACTER(LEN=80):: SlideName

    epsilon = 0.717
    AirPressure = .TRUE.
    WindForce = .FALSE.

! read slide from input.txt
      FILE_NAME=INPUT_FILE_NAME

! storm file
      CALL READ_STRING(SLIDE_FILE,FILE_NAME,'SLIDE_FILE',ierr)
# if defined (PARALLEL)
      if (myid.eq.0) WRITE(3,'(A15,A50)')'SLIDE_FILE:', SLIDE_FILE
# else
      WRITE(3,'(A15,A50)')'SLIDE_FILE:', SLIDE_FILE
# endif

      ALLOCATE (Xco(Mloc),Yco(Nloc))
      
! Xco, and Yco

# if defined (CARTESIAN)
! Cartesian
# if defined (PARALLEL)
![---ykchoi Jan/23/2018
!     Xco(Ibeg) = npx*(Mloc-2*Nghost)*DX
      Xco(Ibeg) = (iista-1)*DX
!---ykchoi Jan/23/2018]
# else
     Xco(Ibeg) = ZERO
# endif
     DO I = Ibeg+1,Mloc
       Xco(I) = Xco(I-1)+DX
     ENDDO
     DO I = Ibeg-1,Ibeg-Nghost,-1
       Xco(I) = Xco(I+1)-DX
     ENDDO

# if defined (PARALLEL)
![---ykchoi Jan/23/2018
!     Yco(Jbeg) = npy*(Nloc-2*Nghost)*DY
      Yco(Jbeg) = (jjsta-1)*DY
!---ykchoi Jan/23/2018]
# else
     Yco(Jbeg) = ZERO
# endif
     DO J = Jbeg+1,Nloc
       Yco(J) = Yco(J-1)+DY
     ENDDO
     DO J = Jbeg-1,Jbeg-Nghost,-1
       Yco(J) = Yco(J+1)-DY
     ENDDO
# endif

    TMP_NAME = TRIM(SLIDE_FILE)

! check existing

 INQUIRE(FILE=TRIM(TMP_NAME),EXIST=FILE_EXIST)
  IF(.NOT.FILE_EXIST)THEN
# if defined (PARALLEL)
   IF(MYID==0)  &
   WRITE(*,*) TRIM(TMP_NAME), 'SLIDE FILE CANNOT BE FOUND. STOP'
   CALL MPI_FINALIZE (ier)
   STOP
# else
    WRITE(*,*) TRIM(TMP_NAME), 'SLIDE FILE CANNOT BE FOUND. STOP'
    STOP
# endif
  ENDIF

! open file
  Ifile=300
  OPEN(Ifile,FILE=TRIM(TMP_NAME))

! read file
         READ(Ifile,*)  ! title
         READ(Ifile,*)  SlideName  !  name
         READ(Ifile,*)  ! info
         READ(Ifile,*) LengthSlide,WidthSlide,&
                       AlphaSlide,BetaSlide,PSlide
         READ(Ifile,*)  ! t,x,y
         READ(Ifile,*)  TimeStorm2,Xstorm2,Ystorm2
                        
         TimeStorm1 = TimeStorm2
         Xstorm1 = Xstorm2
         Ystorm1 = Ystorm2

# if defined (PARALLEL)
   IF(MYID==0)THEN
   WRITE(3,*) 'Slide Name: ',  TRIM(SlideName)
   WRITE(3,*) 'Initial Time, X, Y', TimeStorm2,Xstorm2,Ystorm2
   WRITE(3,*) 'Slide L W A B P = ', LengthSlide,WidthSlide,&
                       AlphaSlide,BetaSlide,PSlide
   ENDIF
# else
   WRITE(*,*) 'Slide Name: ',  TRIM(SlideName)
   WRITE(*,*) 'Initial Time, X, Y', TimeStorm2,Xstorm2,Ystorm2
   WRITE(*,*) 'Slide L W A B P = ', LengthSlide,WidthSlide,&
                       AlphaSlide,BetaSlide,PSlide
   WRITE(3,*) 'Slide Name: ',  TRIM(SlideName)
   WRITE(3,*) 'Initial Time, X, Y', TimeStorm2,Xstorm2,Ystorm2
   WRITE(3,*) 'Slide L W A B P = ', LengthSlide,WidthSlide,&
                       AlphaSlide,BetaSlide,PSlide
# endif

END SUBROUTINE Slide_Model_Setup

SUBROUTINE Slide_Model_Forcing
  USE GLOBAL,ONLY : Mloc,Nloc,tmp1,tmp2,SMALL,TIME,ZERO
  USE INPUT_READ
  IMPLICIT NONE
  INTEGER :: Ifile,ierr,I,J
  REAL(SP) :: Xstorm, Ystorm, C,kb,kw,sech1,sech2


    StormPressureTotal = ZERO

    IF(TIME>TimeStorm1.AND.TIME>TimeStorm2) THEN

         TimeStorm1=TimeStorm2
         Xstorm1 = Xstorm2
         Ystorm1 = Ystorm2

    Ifile = 300

    READ(Ifile,*,END=120)  TimeStorm2,Xstorm2,Ystorm2

# if defined (PARALLEL)
   IF(MYID==0)THEN
     WRITE(3,*)'T,X,Y = ', TimeStorm2,Xstorm2,Ystorm2
     WRITE(*,*)'T,X,Y = ', TimeStorm2,Xstorm2,Ystorm2
   ENDIF
# else
   WRITE(*,*)'T,X,Y = ', TimeStorm2,Xstorm2,Ystorm2
   WRITE(3,*)'T,X,Y = ', TimeStorm2,Xstorm2,Ystorm2
# endif

    ThetaStorm = ATAN2(Ystorm2-Ystorm1,  &
                              Xstorm2-Xstorm1)

    ENDIF ! end time > timestorm2

! calculate force
    tmp2=ZERO
    tmp1=ZERO

    IF(TIME>TimeStorm1)THEN
      IF(TimeStorm1.EQ.TimeStorm2)THEN
        ! no more data
        tmp2=ZERO
        tmp1=ZERO
      ELSE
      tmp2=(TimeStorm2-TIME) &
            /MAX(SMALL, ABS(TimeStorm2-TimeStorm1))
      tmp1=1.0_SP - tmp2;
      ENDIF  ! no more data?
    ENDIF ! time>time_1

    Xstorm = Xstorm2*tmp1 +Xstorm1*tmp2
    Ystorm = Ystorm2*tmp1 +Ystorm1*tmp2

120 CONTINUE  ! no more data for vessel Kves

! sourceX and sourceY

    C = acosh(1.0_SP/epsilon)
    Kb = 2.0_SP*C/MAX(SMALL,WidthSlide)
    kw = 2.0_SP*C/MAX(SMALL,LengthSlide)

    DO I=1,Mloc
      sech1 = Kw*(Xco(I)-Xstorm)
    DO J=1,Nloc
      sech2=Kb*(Yco(J)-Ystorm)

       StormPressureTotal(I,J)= &
         (Pslide/(1.0_SP-epsilon))*(((1.0_SP/cosh(sech1)) &
           *(1.0_SP/cosh(sech2))) - epsilon)                 
         IF (StormPressureTotal(I,J)<0.0) THEN
                 StormPressureTotal(I,J)=0.0   
         ENDIF
    ENDDO
    ENDDO

    IF(FirstCall_Met)THEN
      Eta = -StormPressureTotal
      Eta0 = -StormPressureTotal
      FirstCall_Met = .FALSE.
!open(10,file='tmp.txt')
!  do j=1,nloc
!     write(10,100)(Eta(I,J),I=1,mloc)
!  enddo
!100 format(500f16.6)
!close(10)
    ENDIF

    DO J=Jbeg,Jend
    DO I=Ibeg,Iend
# if defined (CARTESIAN)
       StormPressureX(I,J) = -Grav*H(I,J)*  &
               (StormPressureTotal(I+1,J)-StormPressureTotal(I-1,J))/2.0_SP  &
               /DX
       StormPressureY(I,J) = -Grav*H(I,J)*  &
               (StormPressureTotal(I,J+1)-StormPressureTotal(I,J-1))/2.0_SP  &
               /DY
# endif
    ENDDO
    ENDDO

END SUBROUTINE Slide_Model_Forcing

SUBROUTINE MeteoGausian_Setup
  USE GLOBAL,ONLY : itmp1,itmp2,itmp3,itmp4,SMALL,LARGE,INPUT_FILE_NAME
# if defined (PARALLEL)
  USE GLOBAL,ONLY : iista,jjsta   !ykchoi Jan/23/2018
# endif
  USE INPUT_READ
  IMPLICIT NONE
  CHARACTER(LEN=80)::FILE_NAME=' '
  CHARACTER(LEN=80)::STORM_FILE =' '
  CHARACTER(LEN=80)::TMP_NAME=' '
  INTEGER :: Ifile,ierr
  CHARACTER(LEN=80):: StormName

! read storm from input.txt
      FILE_NAME=INPUT_FILE_NAME

! storm file
      CALL READ_STRING(STORM_FILE,FILE_NAME,'METEO_GAUSIAN_FILE',ierr)
# if defined (PARALLEL)
      if (myid.eq.0) WRITE(3,'(A15,A50)')'STORM_FILE:', STORM_FILE
# else
      WRITE(3,'(A15,A50)')'STORM_FILE:', STORM_FILE
# endif


      ALLOCATE (Xco(Mloc),Yco(Nloc))

      
! Xco, and Yco

# if defined (CARTESIAN)
! Cartesian
# if defined (PARALLEL)
![---ykchoi Jan/23/2018
!     Xco(Ibeg) = npx*(Mloc-2*Nghost)*DX
     Xco(Ibeg) = (iista-1)*DX
!---ykchoi Jan/23/2018]
# else
     Xco(Ibeg) = ZERO
# endif
     DO I = Ibeg+1,Mloc
       Xco(I) = Xco(I-1)+DX
     ENDDO
     DO I = Ibeg-1,Ibeg-Nghost,-1
       Xco(I) = Xco(I+1)-DX
     ENDDO

# if defined (PARALLEL)
![---ykchoi Jan/23/2018
!     Yco(Jbeg) = npy*(Nloc-2*Nghost)*DY
     Yco(Jbeg) = (jjsta-1)*DY
!---ykchoi Jan/23/2018]
# else
     Yco(Jbeg) = ZERO
# endif
     DO J = Jbeg+1,Nloc
       Yco(J) = Yco(J-1)+DY
     ENDDO
     DO J = Jbeg-1,Jbeg-Nghost,-1
       Yco(J) = Yco(J+1)-DY
     ENDDO
# endif

    TMP_NAME = TRIM(STORM_FILE)

! check existing

 INQUIRE(FILE=TRIM(TMP_NAME),EXIST=FILE_EXIST)
  IF(.NOT.FILE_EXIST)THEN
# if defined (PARALLEL)
   IF(MYID==0)  &
   WRITE(*,*) TRIM(TMP_NAME), 'CANNOT BE FOUND. STOP'
   CALL MPI_FINALIZE (ier)
   STOP
# else
    WRITE(*,*) TRIM(TMP_NAME), 'CANNOT BE FOUND. STOP'
    STOP
# endif
  ENDIF

! open file
  Ifile=300
  OPEN(Ifile,FILE=TRIM(TMP_NAME))

! read file
         READ(Ifile,*)  ! title
         READ(Ifile,*)  StormName  !  name
         READ(Ifile,*)  ! t,x,y, Pn, Pc, A and B
         READ(Ifile,*)  TimeStorm2,Xstorm2,Ystorm2,  &
                        DP_storm2,SigmaX2, SigmaY2, Theta2

         TimeStorm1 = TimeStorm2
         Xstorm1 = Xstorm2
         Ystorm1 = Ystorm2
         DP_storm1 = DP_storm2
         SigmaX1 = SigmaX2
         SigmaY1 = SigmaY2
         Theta1 = Theta2

# if defined (PARALLEL)
   IF(MYID==0)THEN
   WRITE(3,*) 'Storm Name: ',  TRIM(StormName)
   WRITE(3,*) 'Initial Time, X, Y', TimeStorm2,Xstorm2,Ystorm2
   WRITE(3,*) 'dP,SigmaX,SigmaY, Theta: ', DP_storm2,SigmaX2, SigmaY2, Theta2
   ENDIF
# else
   WRITE(*,*) 'Storm Name: ',  TRIM(StormName)
   WRITE(*,*) 'Initial Time, X, Y', TimeStorm2,Xstorm2,Ystorm2
   WRITE(*,*) 'dP,SigmaX,SigmaY, Theta: ', DP_storm2,SigmaX2, SigmaY2, Theta2
   WRITE(3,*) 'Storm Name: ',  TRIM(StormName)
   WRITE(3,*) 'Initial Time, X, Y', TimeStorm2,Xstorm2,Ystorm2
   WRITE(3,*) 'dP,SigmaX,SigmaY, Theta: ', DP_storm2,SigmaX2, SigmaY2, Theta2
# endif

END SUBROUTINE MeteoGausian_Setup

SUBROUTINE MeteoGausian_Forcing
  USE GLOBAL,ONLY : Mloc,Nloc,tmp1,tmp2,SMALL,TIME,ZERO
  USE INPUT_READ
  IMPLICIT NONE
  INTEGER :: Ifile,ierr,I,J
  REAL(SP) :: DP_storm,SigmaX, SigmaY, Theta
  REAL(SP) :: Xstorm, Ystorm
  REAL(SP) :: a,b,c

    StormPressureTotal = ZERO

    IF(TIME>TimeStorm1.AND.TIME>TimeStorm2) THEN

         TimeStorm1=TimeStorm2
         Xstorm1 = Xstorm2
         Ystorm1 = Ystorm2

    Ifile = 300

    READ(Ifile,*,END=120)  TimeStorm2,Xstorm2,Ystorm2,  &
                        DP_storm2,SigmaX2, SigmaY2, Theta2

# if defined (PARALLEL)
   IF(MYID==0)THEN
     WRITE(3,*)'T,X,Y = ', TimeStorm2,Xstorm2,Ystorm2
     WRITE(*,*)'T,X,Y = ', TimeStorm2,Xstorm2,Ystorm2
   WRITE(3,*) 'dP,SigmaX,SigmaY, Theta: ', DP_storm2,SigmaX2, SigmaY2, Theta2
   WRITE(*,*) 'dP,SigmaX,SigmaY, Theta: ', DP_storm2,SigmaX2, SigmaY2, Theta2
   ENDIF
# else
   WRITE(*,*)'T,X,Y = ', TimeStorm2,Xstorm2,Ystorm2
   WRITE(*,*) 'dP,SigmaX,SigmaY, Theta: ', DP_storm2,SigmaX2, SigmaY2, Theta2
   WRITE(3,*)'T,X,Y = ', TimeStorm2,Xstorm2,Ystorm2
   WRITE(3,*) 'dP,SigmaX,SigmaY, Theta: ', DP_storm2,SigmaX2, SigmaY2, Theta2
# endif

    ENDIF ! end time > timestorm2

! calculate force
    tmp2=ZERO
    tmp1=ZERO

    IF(TIME>TimeStorm1)THEN
      IF(TimeStorm1.EQ.TimeStorm2)THEN
        ! no more data
        tmp2=ZERO
        tmp1=ZERO
      ELSE
      tmp2=(TimeStorm2-TIME) &
            /MAX(SMALL, ABS(TimeStorm2-TimeStorm1))
      tmp1=1.0_SP - tmp2;
      ENDIF  ! no more data?
    ENDIF ! time>time_1

    Xstorm = Xstorm2*tmp1 +Xstorm1*tmp2
    Ystorm = Ystorm2*tmp1 +Ystorm1*tmp2
    dP_storm = Dp_storm2*tmp1 + dP_storm1*tmp2
    SigmaX = SigmaX2*tmp1 + SigmaX1*tmp2 
    SigmaY= SigmaY2*tmp1 + SigmaY1*tmp2
    Theta = (Theta2*tmp1 + Theta1*tmp2)*PI/180.0_SP

120 CONTINUE  ! no more data for vessel Kves

# if defined (LINEAR_JUMP)

    DO J=1,Nloc
    DO I=1,Mloc

       a = ATAN2(Yco(J)-Ystorm,Xco(I)-Xstorm)
       b = Theta + PI*0.5_SP  ! angle above (Xstorm,Ystorm)
       c = Theta - PI*0.5_SP  ! angle below (Xstorm,Ystorm)

       IF(a > b .OR. a < c)THEN
         StormPressureTotal(I,J) = dP_Storm*0.01_SP ! convert from cm to meter
       ELSE
         StormPressureTotal(I,J) = 0.0_SP
       ENDIF
    ENDDO
    ENDDO

# else
! sourceX and sourceY

    IF(SigmaX .eq. ZERO .OR. SigmaY .eq. ZERO)THEN
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A40,A40)')'SigmaX or SigmaY:', 'ZERO, STOP'
         WRITE(3,'(A40,A40)')'SigmaX or SigmaY:', 'ZERO, STOP'
      endif
       call MPI_FINALIZE ( ier )
# else
         WRITE(*,'(A40,A40)')'SigmaX or SigmaY:', 'ZERO, STOP'
         WRITE(3,'(A40,A40)')'SigmaX or SigmaY:', 'ZERO, STOP'
# endif
        STOP
     
    ENDIF

    a = (COS(Theta))**2/2.0_SP/SigmaX**2  &
              + (SIN(Theta))**2/2.0_SP/SigmaY**2
    b =-SIN(2.0_SP*Theta)/4.0_SP/SigmaX**2   &
              + SIN(2.0_SP*Theta)/4.0_SP/SigmaY**2
    c = (SIN(Theta))**2/2.0_SP/SigmaX**2  &
              + (COS(Theta))**2/2.0_SP/SigmaY**2

    DO J=1,Nloc
    DO I=1,Mloc

       StormPressureTotal(I,J)=dP_Storm*EXP(-(a*(Xco(I)-Xstorm)**2 &
                         +2.0_SP*b*(Xco(I)-Xstorm)*(Yco(J)-Ystorm)  &
                         +c*(Yco(J)-Ystorm)**2))/100.0_SP

       !   /100.0_SP  ! convert from cm to meter


    ENDDO
    ENDDO

# endif
  ! linear jump

    DO J=Jbeg,Jend
    DO I=Ibeg,Iend
# if defined (CARTESIAN)
       StormPressureX(I,J) = -Grav*H(I,J)*  &
               (StormPressureTotal(I+1,J)-StormPressureTotal(I-1,J))/2.0_SP  &
               /DX
       StormPressureY(I,J) = -Grav*H(I,J)*  &
               (StormPressureTotal(I,J+1)-StormPressureTotal(I,J-1))/2.0_SP  &
               /DY
# endif
    ENDDO
    ENDDO

END SUBROUTINE MeteoGausian_Forcing


SUBROUTINE AxisymmetricGaussian_Setup
  USE GLOBAL,ONLY : itmp1,itmp2,itmp3,itmp4,SMALL,LARGE,INPUT_FILE_NAME,&
                    Lon_West,Lat_South,Dphi,Dtheta
# if defined (PARALLEL)
  USE GLOBAL,ONLY : iista,jjsta   !ykchoi Jan/23/2018
# endif
  USE INPUT_READ
  IMPLICIT NONE
  CHARACTER(LEN=80)::FILE_NAME=' '
  CHARACTER(LEN=80)::EVENT_FILE =' '
  CHARACTER(LEN=80)::TMP_NAME=' '
  INTEGER :: Ifile,ierr
  CHARACTER(LEN=80):: EventName

! read storm from input.txt
      FILE_NAME=INPUT_FILE_NAME

! storm file
      CALL READ_STRING(EVENT_FILE,FILE_NAME,'AXISYMMETRIC_GAUSSIAN_FILE',ierr)
# if defined (PARALLEL)
      if (myid.eq.0) WRITE(3,'(A15,A50)')'EVENT_FILE:', EVENT_FILE
# else
      WRITE(3,'(A15,A50)')'EVENT_FILE:', EVENT_FILE
# endif

    ALLOCATE (Xco(Mloc),Yco(Nloc))
      
! Xco, and Yco

# if defined (CARTESIAN)
! Cartesian
# if defined (PARALLEL)
     Xco(Ibeg) = (iista-1)*DX
# else
     Xco(Ibeg) = ZERO
# endif
     DO I = Ibeg+1,Mloc
       Xco(I) = Xco(I-1)+DX
     ENDDO
     DO I = Ibeg-1,Ibeg-Nghost,-1
       Xco(I) = Xco(I+1)-DX
     ENDDO

# if defined (PARALLEL)
     Yco(Jbeg) = (jjsta-1)*DY
# else
     Yco(Jbeg) = ZERO
# endif
     DO J = Jbeg+1,Nloc
       Yco(J) = Yco(J-1)+DY
     ENDDO
     DO J = Jbeg-1,Jbeg-Nghost,-1
       Yco(J) = Yco(J+1)-DY
     ENDDO
	
# else
! Spherical
# if defined (PARALLEL)
     Xco(Ibeg) = (iista-1)*Dphi + Lon_West
# else
     Xco(Ibeg) = Lon_West
# endif
     DO I = Ibeg+1,Mloc
       Xco(I) = Xco(I-1)+Dphi
     ENDDO
     DO I = Ibeg-1,Ibeg-Nghost,-1
       Xco(I) = Xco(I+1)-Dphi
     ENDDO

# if defined (PARALLEL)
     Yco(Jbeg) = (jjsta-1)*Dtheta + Lat_South
# else
     Yco(Jbeg) = Lat_South
# endif
     DO J = Jbeg+1,Nloc
       Yco(J) = Yco(J-1)+Dtheta
     ENDDO
     DO J = Jbeg-1,Jbeg-Nghost,-1
       Yco(J) = Yco(J+1)-Dtheta
     ENDDO

# endif

    TMP_NAME = TRIM(EVENT_FILE)

! check existing

 INQUIRE(FILE=TRIM(TMP_NAME),EXIST=FILE_EXIST)
  IF(.NOT.FILE_EXIST)THEN
# if defined (PARALLEL)
   IF(MYID==0)  &
   WRITE(*,*) TRIM(TMP_NAME), 'CANNOT BE FOUND. STOP'
   CALL MPI_FINALIZE (ier)
   STOP
# else
    WRITE(*,*) TRIM(TMP_NAME), 'CANNOT BE FOUND. STOP'
    STOP
# endif
  ENDIF

! open file
  Ifile=300
  OPEN(Ifile,FILE=TRIM(TMP_NAME))

! read file
         READ(Ifile,*)  ! title
         READ(Ifile,*)  EventName  !  name
         READ(Ifile,*)  ! t0,t1,x0,y0,R0,Pc0,Pt0
         READ(Ifile,*)  TimeStorm1,TimeStorm2,Xstorm1,Ystorm1, &
		                Radius_storm1,PC_storm1,PT_storm1
         READ(Ifile,*)  ! Lc0,Lt0,Tp,psp,alp,eps,cdis
         READ(Ifile,*)  LC_storm1,LT_storm1,TP_storm1,PSp_storm1, &
		                ExpG_storm1,Decay_storm1,CDis_storm1
! close file
  CLOSE(Ifile)

# if defined (PARALLEL)
   IF(MYID==0)THEN
   WRITE(3,*) 'Event Name: ',  TRIM(EventName)
   WRITE(3,*) 'Initial Time, End Time, X, Y: ', &
               TimeStorm1,TimeStorm2,Xstorm1,Ystorm1
   WRITE(3,*) 'Initial Radius, Crest dP, Trough dP: ', &
               Radius_storm1,PC_storm1,PT_storm1
   WRITE(3,*) 'Crest length, Trough length: ', &
               LC_storm1,LT_storm1
   WRITE(3,*) 'Pulse period & Speed, Exponent for Growth: ', &
               TP_storm1,PSp_storm1,ExpG_storm1
   WRITE(3,*) 'Pulse decay factor, Control distance factor: ', &
               Decay_storm1,CDis_storm1
   ENDIF
# else
   WRITE(*,*) 'Event Name: ',  TRIM(EventName)
   WRITE(*,*) 'Initial Time, End Time, X, Y: ', &
               TimeStorm1,TimeStorm2,Xstorm1,Ystorm1
   WRITE(*,*) 'Initial Radius, Crest dP, Trough dP: ', &
               Radius_storm1,PC_storm1,PT_storm1
   WRITE(*,*) 'Crest length, Trough length: ', &
               LC_storm1,LT_storm1
   WRITE(*,*) 'Pulse period & Speed, Exponent for Growth: ', &
               TP_storm1,PSp_storm1,ExpG_storm1
   WRITE(*,*) 'Pulse decay factor, Control distance factor: ', &
               Decay_storm1,CDis_storm1
   WRITE(3,*) 'Event Name: ',  TRIM(EventName)
   WRITE(3,*) 'Initial Time, End Time, X, Y: ', &
               TimeStorm1,TimeStorm2,Xstorm1,Ystorm1
   WRITE(3,*) 'Initial Radius, Crest dP, Trough dP: ', &
               Radius_storm1,PC_storm1,PT_storm1
   WRITE(3,*) 'Crest length, Trough length: ', &
               LC_storm1,LT_storm1
   WRITE(3,*) 'Pulse period & Speed, Exponent for Growth: ', &
               TP_storm1,PSp_storm1,ExpG_storm1
   WRITE(3,*) 'Pulse decay factor, Control distance factor: ', &
               Decay_storm1,CDis_storm1
# endif

    RA_storm = PT_storm1/PC_storm1
    Kap_storm = 2.0_SP*ACOSH(1.0_SP/SQRT(Decay_storm1))
	! Pulse wavelength in far-field (km) (prop R/R0 increase)
    Lp_storm = TP_storm1*PSp_storm1/1000.0_SP
	! Radius to reach far-field
	RLp_storm = Radius_storm1*(Lp_storm/((0.5_SP + CDis_storm1) &
	            *(LC_storm1 + LT_storm1)))**(1.0_SP/ExpG_storm1)
	! Time to reach far-field period
	TLp_storm = TimeStorm1 + (RLp_storm - Radius_storm1) &
	                         *1000.0_SP/PSp_storm1

END SUBROUTINE AxisymmetricGaussian_Setup


SUBROUTINE AxisymmetricGaussian_Forcing
  USE GLOBAL,ONLY : Mloc,Nloc,tmp1,tmp2,SMALL,TIME,ZERO
  USE INPUT_READ
  IMPLICIT NONE
  INTEGER :: Ifile,ierr,I,J
  REAL(SP) :: Radius_storm,Radius_storm0,PC_storm,LC_storm,LT_storm,&
		      DXi_storm,GC_storm,GT_storm
  REAL(SP) :: Dist_storm ! distance to source center (rxy:m)
  REAL(SP) :: Xia ! Relative Gaussian polar coordinates from crest (km)

    StormPressureTotal = ZERO

    IF(TIME>TimeStorm1.AND.TIME<TimeStorm2) THEN
	  Radius_storm = Radius_storm1 + &
	                 (TIME-TimeStorm1)*PSp_storm1/1000.0_SP
	  PC_storm = PC_storm1*Radius_storm1/Radius_storm
	  IF(TIME<TLp_storm) THEN
	    LC_storm = LC_storm1*(Radius_storm/Radius_storm1)**ExpG_storm1
		LT_storm = LT_storm1*(Radius_storm/Radius_storm1)**ExpG_storm1
	  ELSE
	    Radius_storm0 = Radius_storm1 + &
	                 (TLp_storm-TimeStorm1)*PSp_storm1/1000.0_SP
	    LC_storm = LC_storm1*(Radius_storm0/Radius_storm1)**ExpG_storm1
		LT_storm = LT_storm1*(Radius_storm0/Radius_storm1)**ExpG_storm1
	  ENDIF
	  DXi_storm = CDis_storm1*(LC_storm+LT_storm)
      DO J=1,Nloc
      DO I=1,Mloc
# if defined (CARTESIAN)
        Dist_storm = SQRT((Xco(I)-Xstorm1)**2 &
	                     +(Yco(J)-Ystorm1)**2)
# else
        Dist_storm = 2.0_SP*R_earth*asin(SQRT(  sin( (Yco(J)-Ystorm1) &
	                                                 *pi/360.0_SP)**2 &
	                 + cos(Yco(J)*pi/180.0_SP)*cos(Ystorm1*pi/180.0_SP) &
		  			  *sin( (Xco(I)-Xstorm1)*pi/360.0_SP )**2 ) )
# endif
        Xia = Dist_storm/1000.0_SP - Radius_storm
		GC_storm = Zero; GT_storm = Zero;
		IF(Xia<=0.5*LC_storm.AND.Xia>=-(0.5*LT_storm+DXi_storm)) THEN
          GC_storm = (1.0_SP/cosh( Kap_storm*Xia/LC_storm ))**2
		  GT_storm = (1.0_SP/cosh( Kap_storm &
		                        *(Xia+DXi_storm)/LT_storm ))**2
        ENDIF
		StormPressureTotal(I,J) = PC_storm*(GC_storm-RA_storm*GT_storm)/100.0_SP
      ENDDO
      ENDDO
	
    ENDIF ! end time > timestorm2

    DO J=Jbeg,Jend
    DO I=Ibeg,Iend
# if defined (CARTESIAN)
       StormPressureX(I,J) = -Grav*H(I,J)*  &
               (StormPressureTotal(I+1,J)-StormPressureTotal(I-1,J))/2.0_SP  &
               /DX
       StormPressureY(I,J) = -Grav*H(I,J)*  &
               (StormPressureTotal(I,J+1)-StormPressureTotal(I,J-1))/2.0_SP  &
               /DY
# else
       StormPressureX(I,J) = -Grav*H(I,J)*  &
               (StormPressureTotal(I+1,J)-StormPressureTotal(I-1,J))/2.0_SP  &
               /DX(i,j)
       StormPressureY(I,J) = -Grav*H(I,J)*  &
               (StormPressureTotal(I,J+1)-StormPressureTotal(I,J-1))/2.0_SP  &
               /DY(i,j)
# endif
    ENDDO
    ENDDO

END SUBROUTINE AxisymmetricGaussian_Forcing


SUBROUTINE LPGaussian_Setup
  USE GLOBAL,ONLY : itmp1,itmp2,itmp3,itmp4,SMALL,LARGE,INPUT_FILE_NAME,&
                    Lon_West,Lat_South,Dphi,Dtheta
# if defined (PARALLEL)
  USE GLOBAL,ONLY : iista,jjsta   !ykchoi Jan/23/2018
# endif
  USE INPUT_READ
  IMPLICIT NONE
  CHARACTER(LEN=80)::FILE_NAME=' '
  CHARACTER(LEN=80)::EVENT_FILE =' '
  CHARACTER(LEN=80)::TMP_NAME=' '
  INTEGER :: Ifile,ierr
  CHARACTER(LEN=80):: EventName

! read storm from input.txt
      FILE_NAME=INPUT_FILE_NAME

! storm file
      CALL READ_STRING(EVENT_FILE,FILE_NAME,'LP_GAUSSIAN_FILE',ierr)
# if defined (PARALLEL)
      if (myid.eq.0) WRITE(3,'(A15,A50)')'EVENT_FILE:', EVENT_FILE
# else
      WRITE(3,'(A15,A50)')'EVENT_FILE:', EVENT_FILE
# endif

    ALLOCATE (Xco(Mloc),Yco(Nloc))
      
! Xco, and Yco

# if defined (CARTESIAN)
! Cartesian
# if defined (PARALLEL)
     Xco(Ibeg) = (iista-1)*DX
# else
     Xco(Ibeg) = ZERO
# endif
     DO I = Ibeg+1,Mloc
       Xco(I) = Xco(I-1)+DX
     ENDDO
     DO I = Ibeg-1,Ibeg-Nghost,-1
       Xco(I) = Xco(I+1)-DX
     ENDDO

# if defined (PARALLEL)
     Yco(Jbeg) = (jjsta-1)*DY
# else
     Yco(Jbeg) = ZERO
# endif
     DO J = Jbeg+1,Nloc
       Yco(J) = Yco(J-1)+DY
     ENDDO
     DO J = Jbeg-1,Jbeg-Nghost,-1
       Yco(J) = Yco(J+1)-DY
     ENDDO
	
# else
! Spherical
# if defined (PARALLEL)
     Xco(Ibeg) = (iista-1)*Dphi + Lon_West
# else
     Xco(Ibeg) = Lon_West
# endif
     DO I = Ibeg+1,Mloc
       Xco(I) = Xco(I-1)+Dphi
     ENDDO
     DO I = Ibeg-1,Ibeg-Nghost,-1
       Xco(I) = Xco(I+1)-Dphi
     ENDDO

# if defined (PARALLEL)
     Yco(Jbeg) = (jjsta-1)*Dtheta + Lat_South
# else
     Yco(Jbeg) = Lat_South
# endif
     DO J = Jbeg+1,Nloc
       Yco(J) = Yco(J-1)+Dtheta
     ENDDO
     DO J = Jbeg-1,Jbeg-Nghost,-1
       Yco(J) = Yco(J+1)-Dtheta
     ENDDO

# endif

    TMP_NAME = TRIM(EVENT_FILE)

! check existing

 INQUIRE(FILE=TRIM(TMP_NAME),EXIST=FILE_EXIST)
  IF(.NOT.FILE_EXIST)THEN
# if defined (PARALLEL)
   IF(MYID==0)  &
   WRITE(*,*) TRIM(TMP_NAME), 'CANNOT BE FOUND. STOP'
   CALL MPI_FINALIZE (ier)
   STOP
# else
    WRITE(*,*) TRIM(TMP_NAME), 'CANNOT BE FOUND. STOP'
    STOP
# endif
  ENDIF

! open file
  Ifile=300
  OPEN(Ifile,FILE=TRIM(TMP_NAME))

! read file
         READ(Ifile,*)  ! title
         READ(Ifile,*)  EventName  !  name
         READ(Ifile,*)  ! t0,t1,x0,y0,speed,speed_t,pA_source,trough_pA_source
         READ(Ifile,*)  TimeStorm1,TimeStorm2,Xstorm1,Ystorm1, &
		                PSpC_storm1,PSpT_storm1,PAC_storm1,PAT_storm1
! close file
  CLOSE(Ifile)

# if defined (PARALLEL)
   IF(MYID==0)THEN
   WRITE(3,*) 'Event Name: ',  TRIM(EventName)
   WRITE(3,*) 'Initial Time, End Time, X, Y: ', &
               TimeStorm1,TimeStorm2,Xstorm1,Ystorm1
   WRITE(3,*) 'Crest moving speed, Trough moving speed: ', &
               PSpC_storm1,PSpT_storm1
   WRITE(3,*) 'Crest dP*surface area, Trough dP*surface area: ', &
               PAC_storm1,PAT_storm1
   ENDIF
# else
   WRITE(*,*) 'Event Name: ',  TRIM(EventName)
   WRITE(*,*) 'Initial Time, End Time, X, Y: ', &
               TimeStorm1,TimeStorm2,Xstorm1,Ystorm1
   WRITE(*,*) 'Crest dP*surface area, Trough dP*surface area: ', &
               PSpC_storm1,PSpT_storm1
   WRITE(*,*) 'Crest dP, Trough dP: ', &
               PAC_storm1,PAT_storm1
   WRITE(3,*) 'Event Name: ',  TRIM(EventName)
   WRITE(3,*) 'Initial Time, End Time, X, Y: ', &
               TimeStorm1,TimeStorm2,Xstorm1,Ystorm1
   WRITE(3,*) 'Crest dP*surface area, Trough dP*surface area: ', &
               PSpC_storm1,PSpT_storm1
   WRITE(3,*) 'Crest dP, Trough dP: ', &
               PAC_storm1,PAT_storm1
# endif

END SUBROUTINE LPGaussian_Setup


SUBROUTINE LPGaussian_Forcing
  USE GLOBAL,ONLY : Mloc,Nloc,tmp1,tmp2,SMALL,TIME,ZERO
  USE INPUT_READ
  IMPLICIT NONE
  INTEGER :: Ifile,ierr,I,J
  REAL(SP) :: FDC_storm,FDT_storm ! front distance from center to crest and trough (m)
  REAL(SP) :: WPC_storm,WPT_storm ! Wave period of crest and trough (s)
  REAL(SP) :: WLC_storm,WLT_storm ! Wave length of crest and trough (m)
  REAL(SP) :: SA_storm ! surface area of guassian ring (m^2)
  REAL(SP) :: PC_storm,PT_storm ! max pressure jump at crest and trough
  REAL(SP) :: Dist_storm ! distance to source center (rxy:m)
  REAL(SP) :: LAC_storm,LAT_storm! Relative Gaussian polar coordinates from crest and trough
  REAL(SP) :: PSpC_storm,PSpT_storm 

    StormPressureTotal = ZERO

    IF(TIME>TimeStorm1.AND.TIME<TimeStorm2) THEN
	  PSpC_storm = PSpC_storm1*tanh((TIME-TimeStorm1+300.0_SP)/4000.0_SP)
	  FDC_storm = (TIME-TimeStorm1)*PSpC_storm
	  PSpT_storm = PSpT_storm1*tanh((TIME-TimeStorm1)/7000.0_SP)
	  FDT_storm = (TIME-TimeStorm1)*PSpT_storm
	  WPC_storm = 2300.0_SP*(TIME-TimeStorm1)**(-0.1)
	  WPT_storm = 3300.0_SP*(TIME-TimeStorm1)**(-0.1)
	  WLC_storm = PSpC_storm*WPC_storm
	  WLT_storm = PSpT_storm*WPT_storm
	  SA_storm = 2.0*pi*FDC_storm**0.5
	  PC_storm = MAX(0.0,PAC_storm1/SA_storm*(tanh((FDC_storm)/300000.0_SP))**2.0)
	  SA_storm = 2.0*pi*(FDC_storm+200000.0_SP)**0.91
	  PT_storm = PAT_storm1/SA_storm
      DO J=1,Nloc
      DO I=1,Mloc
# if defined (CARTESIAN)
        Dist_storm = SQRT((Xco(I)-Xstorm1)**2 &
	                     +(Yco(J)-Ystorm1)**2)
# else
        Dist_storm = 2.0_SP*R_earth*asin(SQRT(  sin( (Yco(J)-Ystorm1) &
	                                                 *pi/360.0_SP)**2 &
	                 + cos(Yco(J)*pi/180.0_SP)*cos(Ystorm1*pi/180.0_SP) &
		  			  *sin( (Xco(I)-Xstorm1)*pi/360.0_SP )**2 ) )
# endif
        LAC_storm = (Dist_storm-FDC_storm)/(WLC_storm/1.5)
		LAT_storm = (Dist_storm-FDT_storm)/(WLT_storm/1.5)
		StormPressureTotal(I,J) = (PC_storm*EXP(-LAC_storm**2.0) &
                           		 - PT_storm*EXP(-LAT_storm**2.0))/10000.0_SP
      ENDDO
      ENDDO
	
    ENDIF ! end time > timestorm2

    DO J=Jbeg,Jend
    DO I=Ibeg,Iend
# if defined (CARTESIAN)
       StormPressureX(I,J) = -Grav*H(I,J)*  &
               (StormPressureTotal(I+1,J)-StormPressureTotal(I-1,J))/2.0_SP  &
               /DX
       StormPressureY(I,J) = -Grav*H(I,J)*  &
               (StormPressureTotal(I,J+1)-StormPressureTotal(I,J-1))/2.0_SP  &
               /DY
# else
       StormPressureX(I,J) = -Grav*H(I,J)*  &
               (StormPressureTotal(I+1,J)-StormPressureTotal(I-1,J))/2.0_SP  &
               /DX(i,j)
       StormPressureY(I,J) = -Grav*H(I,J)*  &
               (StormPressureTotal(I,J+1)-StormPressureTotal(I,J-1))/2.0_SP  &
               /DY(i,j)
# endif
    ENDDO
    ENDDO

END SUBROUTINE LPGaussian_Forcing

END MODULE METEO_MODULE

# endif 
! end meteo
