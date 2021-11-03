
!
! (C) Copyright 2021 Met Office UK
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

!--


subroutine fft_gpoint2spectral_f(A, TRIGS, FACTORS, INC, JUMP, N, LOT)

!
! this is a modified version of the function:
!   'FFT991' [C. Temperton, 1987]
!
! https://confluence.ecmwf.int/pages/viewpage.action?pageId=50044238
! https://www.ecmwf.int/en/elibrary/12606-self-sorting-mixed-radix-fast-fourier-transforms
!
! the function is used to perform 'multiple Fast Fourier Transforms' (multiple fast
! real periodic transforms)
!
!
! real transform of length N performed by removing redundant operations
! from complex transform of length N (gridpoint to spetral)
!
! - input:
!   A is the array containing the input data
!   TRIGS is a previously prepared list of trig function values
!   FACTORS is a previously prepared list of factors of N
!   INC is the increment within each data 'vector'
!       (e.g. INC=1 for consecutively stored data)
!   JUMP is the increment between the start of each data vector
!   N is the length of the data vectors
!   LOT is the number of data vectors
!
! - output:
!   A is the array containing the output data; note that input A is overwritten by the results;
!
!
! note that ...
!
! - the length of the data vectors (and the tranforms) N must be an even number
!   that has no other factors except possibly power of 2, 3, and 5
! - the length of the array A is LOT*(N+2);
!   within the array A, each data vector (input) X contains N elements and 2 extra zeros at the end;
! - not suitable for single data vector, i.e. LOT>1
!
!
!                              ***
!
!   Definition of transforms:
!   -------------------------
!
!   A(K)=(1/N)*SUM(J=0,...,N-1)(X(J)*COS(2*J*K*PI/N))
!   B(K)=-(1/N)*SUM(J=0,...,N-1)(X(J)*SIN(2*J*K*PI/N))
!
!   Ordering of data (input):
!       X(0),X(1),X(2),...,X(N-1), 0 , 0 ; (N+2) locations required
!
!   Ordering of coefficients:
!       A(0),B(0),A(1),B(1),A(2),B(2),...,A(N/2),B(N/2)
!       where B(0)=B(N/2)=0 ;              (N+2) locations required
!
!
! note that ...
!
! - the data vector (input) X contains N elements and 2 extra zeros at the end;
!   the (full) length of the data vector X is (N+2)
! - the definition of the transform includes the scale factor (1/N)
!

!--

USE fft_init_f, ONLY: &
    BLOCK_LENGTH

USE fft_pass_proc_f

IMPLICIT NONE

! subroutine arguments
INTEGER, INTENT(IN)  :: INC         ! Increment between elements of data vector
INTEGER, INTENT(IN)  :: JUMP        ! Increment between start of each data vector
INTEGER, INTENT(IN)  :: N           ! Length of data vector in grid-point space
                                    ! without extra zeroes
INTEGER, INTENT(IN)  :: LOT         ! Number of data vectors
INTEGER, INTENT(IN)  :: FACTORS(10) ! List of factors of N
REAL, INTENT(IN)     :: TRIGS(N)    ! Trigonometrical functions
REAL, INTENT (INOUT) :: A( INC*(N+2) + JUMP*(LOT-1) ) ! data array

! workspace; 1 array is required
REAL :: WORK((N+2)*BLOCK_LENGTH) ! general workspace

! local variables
INTEGER :: NFAX       ! Number of factors
INTEGER :: NX         ! N+1 except where N is odd then holds N
INTEGER :: NBLOX      ! Number of blocks LOT is split into
INTEGER :: NB         ! Do loop counter
INTEGER :: ISTART     ! Start address for a block
INTEGER :: NVEX       ! Number of elements in vector
INTEGER :: IA         ! Used to pass ISTART to fft_pass_proc_f
INTEGER :: IX         ! Variable used for addressing
INTEGER :: LA         ! Variable used for addressing
INTEGER :: IGO        ! A control variable
INTEGER :: K          ! Do loop counter
INTEGER :: IFAC       ! Holds current factor
INTEGER :: IERR       ! Holds error status
INTEGER :: I,J,II,IZ,JJ,IBASE,JBASE  ! loop and indexing variables.


! ----------------------------------------------------------------
!    Section 1. Set up information for section 2.
! ----------------------------------------------------------------

WORK = 0.0

! Set number of factors and NX

NFAX = FACTORS(1)
NX = N+1
IF (MOD(N,2) == 1) NX = N

! Calculate number of blocks of BLOCK_LENGTH data vectors are to be split into

NBLOX = 1+(LOT-1)/BLOCK_LENGTH
NVEX = LOT-(NBLOX-1)*BLOCK_LENGTH

! ----------------------------------------------------------------
!    Section 2. Gridpoint to spectral transform
! ----------------------------------------------------------------

ISTART = 1
DO NB=1,NBLOX
  IA = ISTART
  LA = N
  IGO = 1

  DO K=1,NFAX
    IFAC = FACTORS(NFAX+2-K)
    LA = LA/IFAC
    IERR = -1
    IF (IGO == 1) THEN
      CALL fft_qpassm_f(A(IA),             &
                        A(IA+IFAC*LA*INC), &
                        WORK(1),           &
                        WORK(LA+1),        &
                        TRIGS,             &
                        INC,               &
                        1,                 &
                        JUMP,              &
                        NX,                &
                        NVEX,              &
                        N,                 &
                        IFAC,              &
                        LA,                &
                        IERR)
    ELSE
      CALL fft_qpassm_f(WORK(1),         &
                        WORK(IFAC*LA+1), &
                        A(IA),           &
                        A(IA+LA*INC),    &
                        TRIGS,           &
                        1,               &
                        INC,             &
                        NX,              &
                        JUMP,            &
                        NVEX,            &
                        N,               &
                        IFAC,            &
                        LA,              &
                        IERR)
    END IF
    IGO=-IGO
    IA=ISTART+INC
  END DO

  ! If necessary, copy results back to A
  !  ------------------------------------

  IF (MOD(NFAX,2) /= 0) THEN
    IBASE = 1
    JBASE = IA
    DO JJ=1,NVEX
      I = IBASE
      J = JBASE
      DO II=1,N
        A(J) = WORK(I)
        I = I+1
        J = J+INC
      END DO
      IBASE = IBASE+NX
      JBASE = JBASE+JUMP
    END DO
  END IF

  ! Shift A(0) & fill in zero imag parts
  ! ------------------------------------

  IX = ISTART
  DO J=1,NVEX
    A(IX) = A(IX+INC)
    A(IX+INC) = 0.0
    IX = IX+JUMP
  END DO
  IF (MOD(N,2) /= 1) THEN
    IZ = ISTART+(N+1)*INC
    DO J=1,NVEX
      A(IZ) = 0.0
      IZ = IZ+JUMP
    END DO
  END IF
  ISTART = ISTART+NVEX*JUMP
  NVEX = BLOCK_LENGTH

END DO


end subroutine fft_gpoint2spectral_f
