
!
! (C) Copyright 2021 Met Office UK
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

!--


module fft_init_f

  implicit none


  ! factors [dimension: 10]
  INTEGER :: FACTORS(10)

  ! trigonometric coefficient
  REAL, ALLOCATABLE :: TRIGS(:)

  ! increment between successive elements;
  ! e.g. inc=1 for consecutively stored data;
  INTEGER, PARAMETER :: INC = 1

  ! FFT subroutines process data in blocks of 'block_length';
  ! this may have different optimum values for different processes;
  ! 'block_length' can be any positive (non-zero) integer, but values that are
  ! too large or small will result in severe increases in CPU time;
  INTEGER, PARAMETER :: BLOCK_LENGTH = 128


  contains

  subroutine initialize_FFT(N, FACTORS, TRIGS)
  
  !
  ! this is a modified version of the function:
  !   'SET99' [C. Temperton, 1987]
  !
  ! the function is used to initiliaze the procedures for
  ! multiple Fast Fourier Transforms;
  !


  IMPLICIT NONE

  ! subroutine arguments
  INTEGER, INTENT(IN)  :: N
  INTEGER, INTENT(OUT) :: FACTORS(10)
  REAL, ALLOCATABLE, INTENT(OUT) :: TRIGS(:)

  ! local variables
  INTEGER :: I
  INTEGER :: NFAX
  INTEGER :: JFAX(10)
  INTEGER :: LFAX(7)
  INTEGER :: NIL
  INTEGER :: NHL
  INTEGER :: L
  INTEGER :: K
  INTEGER :: NU
  INTEGER :: IFAC
  REAL    :: ANGLE
  REAL    :: DEL

  JFAX = 0
  DATA LFAX/6,8,5,4,3,2,1/

  allocate(TRIGS(3*N))

  ! first initialization
  FACTORS(:) = 0
  TRIGS(:) = 0.0

  !     TRIGS FOR REAL PERIODIC TRANSFORM
  !     ---------------------------------
  DEL=4.0*ASIN(1.0)/N
  NIL=0
  NHL=(N/2)-1
  !CDIR NODEP
  DO K=NIL,NHL
    ANGLE=K*DEL
    TRIGS(2*K+1)=COS(ANGLE)
    TRIGS(2*K+2)=SIN(ANGLE)
  END DO

  !     EXTRA TRIGS FOR (CO)SINE TRANSFORM
  !     ----------------------------------

  DEL=0.5*DEL
  DO K=1,NHL
    ANGLE=K*DEL
    TRIGS(2*N+K)=SIN(ANGLE)
  END DO

  !     EXTRA TRIGS FOR SHIFTED (CO)SINE TRANSFORM
  !     ------------------------------------------
  DEL=0.5*DEL
  DO K=1,N
    ANGLE=K*DEL
    TRIGS(N+K)=SIN(ANGLE)
  END DO

  !     NOW FIND FACTORS OF N
  !     ---------------------
  !
  !     FIND FACTORS OF N (8,6,5,4,3,2; ONLY ONE 8 ALLOWED)
  !     LOOK FOR SIXES FIRST, STORE FACTORS IN DESCENDING ORDER
  NU=N
  IFAC=6
  K=0
  L=1
  DO
    IF (MOD(NU,IFAC) == 0) THEN
      K=K+1
      JFAX(K)=IFAC
      IF (IFAC == 8) THEN
        IF (K /= 1) THEN
          JFAX(1)=8
          JFAX(k)=6
        END IF
      END IF
      NU=NU/IFAC
      IF (NU == 1) EXIT
      IF (IFAC /= 8) CYCLE
    END IF
    L=L+1
    IFAC=LFAX(L)
    IF (IFAC == 1) EXIT
  END DO
  !     NOW REVERSE ORDER OF FACTORS
  NFAX=K
  FACTORS(1)=NFAX
  DO I=1,NFAX
    FACTORS(NFAX+2-I)=JFAX(I)
  END DO


  end subroutine initialize_FFT


end module fft_init_f
