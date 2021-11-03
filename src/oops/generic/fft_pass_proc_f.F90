
!
! (C) Copyright 2021 Met Office UK
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

!--


module fft_pass_proc_f

  implicit none

  contains

  subroutine fft_qpassm_f(A, B, C, D, TRIGS, INC1, INC2, INC3, INC4, LOT, N, IFAC, LA, IERR)

  !
  ! the function is used to carry out data processing -- i.e. to perform one pass
  ! through data as part of multiple FFTs (Temperton method)
  !
  !
  ! A is first real input vector
  !   Equivalence B(1) with A (IFAC*LA*INC1+1)
  ! C is first real output vector
  !   Equivalence D(1) with C(LA*INC2+1)
  ! TRIGS is a precalculated list of sines and cosines
  ! INC1 is the addressing increment for A
  ! INC2 is the addressing increment for C
  ! INC3 is the increment between input vectors A
  ! INC4 is the increment between output vectors C
  ! LOT is the number of vectors
  ! N is the length of the vectors
  ! IFAC is the current factor of N
  ! LA is the product of previous factors
  ! IERR is an error indicator
  !   0 - Pass completed without error
  !   1 - LOT greater than BLOCK_LENGTH
  !   2 - IFAC not catered for
  !   3 - IFAC only catered for if LA=N/IFAC
  !

  !--

  USE fft_init_f, ONLY: &
      BLOCK_LENGTH

  IMPLICIT NONE


  INTEGER, INTENT(IN)  :: N
  INTEGER, INTENT(IN)  :: INC1,INC2,INC3,INC4,LOT,IFAC,LA
  INTEGER, INTENT(OUT) :: IERR
  REAL, INTENT(IN)     :: A(*)
  REAL, INTENT(IN)     :: B(*)
  REAL, INTENT(INOUT)  :: C(*)
  REAL, INTENT(INOUT)  :: D(*)
  REAL, INTENT(IN)     :: TRIGS(N)

  INTEGER :: IH, IG, JF, KF, IF, IE, JE, KE, ID, JD, KD, KC
  INTEGER :: I, IJK, K, KB, IC, JC, IA, IB, JA, JB, L, J, IGO, M
  INTEGER :: IINK, JINK, IJUMP, KSTOP, IBAD, IBASE, JBASE
  REAL    :: C5, S5, ZSIN45, ZSIN72, ZSIN36, B21, B11, A11, A21
  REAL    :: A10, A20, B10, B20, C4, S4, B4, B5, B6, ZQRT5, A5, A6
  REAL    :: A4, B0, C3, S3, A0, SIN45
  REAL    :: ZSIN60, B3, A3, B2, A2, A1, B1, C2, C1, S1
  REAL    :: S2, Z
  REAL, PARAMETER :: SIN36 = 0.587785252292473
  REAL, PARAMETER :: SIN60 = 0.866025403784437
  REAL, PARAMETER :: SIN72 = 0.951056516295154
  REAL, PARAMETER :: QRT5 = 0.559016994374947

  LOGICAL :: FUPDATE

  ! if this flag is set to 'true', set 'IBAD=0' to indicate
  ! that the pass is completed without error
  FUPDATE = .false.

  M = N/IFAC
  IINK = LA*INC1
  JINK = LA*INC2
  IJUMP = (IFAC-1)*IINK
  KSTOP = (N-IFAC)/(2*IFAC)

  ! error: LOT greater than BLOCK_LENGTH
  IBAD = 1

  IF (LOT <= BLOCK_LENGTH) THEN

    IBASE = 0
    JBASE = 0
    IGO = IFAC-1
    ! error: IFAC not catered for (default case)
    IBAD = 2

    SELECT CASE (IGO)

    !--
    CASE (1)
      ! Coding for factor 2

      IA = 1
      IB = IA+IINK
      JA = 1
      JB = JA+(2*M-LA)*INC2

      IF (LA /= M) THEN
        DO L=1,LA
          I = IBASE
          J = JBASE
    !CDIR NODEP
          DO IJK=1,LOT
            C(JA+J) = A(IA+I)+A(IB+I)
            C(JB+J) = A(IA+I)-A(IB+I)
            I = I+INC3
            J = J+INC4
          END DO
          IBASE = IBASE+INC1
          JBASE = JBASE+INC2
        END DO
        JA = JA+JINK
        JINK = 2*JINK
        JB = JB-JINK
        IBASE = IBASE+IJUMP
        IJUMP = 2*IJUMP+IINK
        IF (JA /= JB) THEN
          DO K=LA,KSTOP,LA
            KB = K+K
            C1 = TRIGS(KB+1)
            S1 = TRIGS(KB+2)
            JBASE = 0
            DO L=1,LA
              I = IBASE
              J = JBASE
    !CDIR NODEP
              DO IJK=1,LOT
                C(JA+J) = A(IA+I)+(C1*A(IB+I)+S1*B(IB+I))
                C(JB+J) = A(IA+I)-(C1*A(IB+I)+S1*B(IB+I))
                D(JA+J) = (C1*B(IB+I)-S1*A(IB+I))+B(IA+I)
                D(JB+J) = (C1*B(IB+I)-S1*A(IB+I))-B(IA+I)
                I = I+INC3
                J = J+INC4
              END DO
              IBASE = IBASE+INC1
              JBASE = JBASE+INC2
            END DO
            IBASE = IBASE+IJUMP
            JA = JA+JINK
            JB = JB-JINK
          END DO
          IF (JA > JB) THEN
            FUPDATE = .true.
          END IF
        END IF
        IF (FUPDATE .eqv. .false.) THEN
          JBASE = 0
          DO L=1,LA
            I = IBASE
            J = JBASE
      !CDIR NODEP
            DO IJK=1,LOT
              C(JA+J) = A(IA+I)
              D(JA+J) = -A(IB+I)
              I = I+INC3
              J = J+INC4
            END DO
            IBASE = IBASE+INC1
            JBASE = JBASE+INC2
          END DO
          FUPDATE = .true.
        END IF
      END IF
      IF (FUPDATE .eqv. .false.) THEN
        Z = 1.0/N
        DO L=1,LA
          I = IBASE
          J = JBASE
      !CDIR NODEP
          DO IJK=1,LOT
            C(JA+J) = Z*(A(IA+I)+A(IB+I))
            C(JB+J) = Z*(A(IA+I)-A(IB+I))
            I = I+INC3
            J = J+INC4
          END DO
          IBASE = IBASE+INC1
          JBASE = JBASE+INC2
        END DO
        FUPDATE = .true.
      END IF

    !--
    CASE (2)
      ! coding for factor 3

      IA = 1
      IB = IA+IINK
      IC = IB+IINK
      JA = 1
      JB = JA+(2*M-LA)*INC2
      JC = JB

      IF (LA /= M) THEN
        DO L=1,LA
          I = IBASE
          J = JBASE
    !CDIR NODEP
          DO IJK=1,LOT
            C(JA+J) = A(IA+I)+(A(IB+I)+A(IC+I))
            C(JB+J) = A(IA+I)-0.5*(A(IB+I)+A(IC+I))
            D(JB+J) = SIN60*(A(IC+I)-A(IB+I))
            I = I+INC3
            J = J+INC4
          END DO
          IBASE = IBASE+INC1
          JBASE = JBASE+INC2
        END DO
        JA = JA+JINK
        JINK = 2*JINK
        JB = JB+JINK
        JC = JC-JINK
        IBASE = IBASE+IJUMP
        IJUMP = 2*IJUMP+IINK
        IF (JA /= JC) THEN
          DO K=LA,KSTOP,LA
            KB = K+K
            KC = KB+KB
            C1 = TRIGS(KB+1)
            S1 = TRIGS(KB+2)
            C2 = TRIGS(KC+1)
            S2 = TRIGS(KC+2)
            JBASE = 0
            DO L=1,LA
              I = IBASE
              J = JBASE
    !CDIR NODEP
              DO IJK=1,LOT
                A1 = (C1*A(IB+I)+S1*B(IB+I))+(C2*A(IC+I)+S2*B(IC+I))
                B1 = (C1*B(IB+I)-S1*A(IB+I))+(C2*B(IC+I)-S2*A(IC+I))
                A2 = A(IA+I)-0.5*A1
                B2 = B(IA+I)-0.5*B1
                A3 = SIN60*((C1*A(IB+I)+S1*B(IB+I))-(C2*A(IC+I)+S2*B(IC+I)))
                B3 = SIN60*((C1*B(IB+I)-S1*A(IB+I))-(C2*B(IC+I)-S2*A(IC+I)))
                C(JA+J) = A(IA+I)+A1
                D(JA+J) = B(IA+I)+B1
                C(JB+J) = A2+B3
                D(JB+J) = B2-A3
                C(JC+J) = A2-B3
                D(JC+J) = -(B2+A3)
                I = I+INC3
                J = J+INC4
              END DO
              IBASE = IBASE+INC1
              JBASE = JBASE+INC2
            END DO
            IBASE = IBASE+IJUMP
            JA = JA+JINK
            JB = JB+JINK
            JC = JC-JINK
          END DO
          IF (JA > JC) THEN
            FUPDATE = .true.
          END IF
        END IF
        IF (FUPDATE .eqv. .false.) THEN
          JBASE = 0
          DO L=1,LA
            I = IBASE
            J = JBASE
      !CDIR NODEP
            DO IJK=1,LOT
              C(JA+J) = A(IA+I)+0.5*(A(IB+I)-A(IC+I))
              D(JA+J) = -SIN60*(A(IB+I)+A(IC+I))
              C(JB+J) = A(IA+I)-(A(IB+I)-A(IC+I))
              I = I+INC3
              J = J+INC4
            END DO
            IBASE = IBASE+INC1
            JBASE = JBASE+INC2
          END DO
          FUPDATE = .true.
        END IF
      END IF
      IF (FUPDATE .eqv. .false.) THEN
        Z = 1.0/N
        ZSIN60 = Z*SIN60
        DO L=1,LA
          I = IBASE
          J = JBASE
      !CDIR NODEP
          DO IJK=1,LOT
            C(JA+J) = Z*(A(IA+I)+(A(IB+I)+A(IC+I)))
            C(JB+J) = Z*(A(IA+I)-0.5*(A(IB+I)+A(IC+I)))
            D(JB+J) = ZSIN60*(A(IC+I)-A(IB+I))
            I = I+INC3
            J = J+INC4
          END DO
          IBASE = IBASE+INC1
          JBASE = JBASE+INC2
        END DO
        FUPDATE = .true.
      END IF

    !--
    CASE (3)
      ! Coding for factor 4

      IA = 1
      IB = IA+IINK
      IC = IB+IINK
      ID = IC+IINK
      JA = 1
      JB = JA+(2*M-LA)*INC2
      JC = JB+2*M*INC2
      JD = JB

      IF (LA /= M) THEN
        DO L=1,LA
          I = IBASE
          J = JBASE
    !CDIR NODEP
          DO IJK=1,LOT
            C(JA+J) = (A(IA+I)+A(IC+I))+(A(IB+I)+A(ID+I))
            C(JC+J) = (A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I))
            C(JB+J) = A(IA+I)-A(IC+I)
            D(JB+J) = A(ID+I)-A(IB+I)
            I = I+INC3
            J = J+INC4
          END DO
          IBASE = IBASE+INC1
          JBASE = JBASE+INC2
        END DO
        JA = JA+JINK
        JINK = 2*JINK
        JB = JB+JINK
        JC = JC-JINK
        JD = JD-JINK
        IBASE = IBASE+IJUMP
        IJUMP = 2*IJUMP+IINK
        IF (JB /= JC) THEN
          DO K=LA,KSTOP,LA
            KB = K+K
            KC = KB+KB
            KD = KC+KB
            C1 = TRIGS(KB+1)
            S1 = TRIGS(KB+2)
            C2 = TRIGS(KC+1)
            S2 = TRIGS(KC+2)
            C3 = TRIGS(KD+1)
            S3 = TRIGS(KD+2)
            JBASE = 0
            DO L=1,LA
              I = IBASE
              J = JBASE
    !CDIR NODEP
              DO IJK=1,LOT
                A0 = A(IA+I)+(C2*A(IC+I)+S2*B(IC+I))
                A2 = A(IA+I)-(C2*A(IC+I)+S2*B(IC+I))
                A1 = (C1*A(IB+I)+S1*B(IB+I))+(C3*A(ID+I)+S3*B(ID+I))
                A3 = (C1*A(IB+I)+S1*B(IB+I))-(C3*A(ID+I)+S3*B(ID+I))
                B0 = B(IA+I)+(C2*B(IC+I)-S2*A(IC+I))
                B2 = B(IA+I)-(C2*B(IC+I)-S2*A(IC+I))
                B1 = (C1*B(IB+I)-S1*A(IB+I))+(C3*B(ID+I)-S3*A(ID+I))
                B3 = (C1*B(IB+I)-S1*A(IB+I))-(C3*B(ID+I)-S3*A(ID+I))
                C(JA+J) = A0+A1
                C(JC+J) = A0-A1
                D(JA+J) = B0+B1
                D(JC+J) = B1-B0
                C(JB+J) = A2+B3
                C(JD+J) = A2-B3
                D(JB+J) = B2-A3
                D(JD+J) = -(B2+A3)
                I = I+INC3
                J = J+INC4
              END DO
              IBASE = IBASE+INC1
              JBASE = JBASE+INC2
            END DO
            IBASE = IBASE+IJUMP
            JA = JA+JINK
            JB = JB+JINK
            JC = JC-JINK
            JD = JD-JINK
          END DO
          IF (JB > JC) THEN
            FUPDATE = .true.
          END IF
        END IF
        IF (FUPDATE .eqv. .false.) THEN
          SIN45 = SQRT(0.5)
          JBASE = 0
          DO L=1,LA
            I = IBASE
            J = JBASE
      !CDIR NODEP
            DO IJK=1,LOT
              C(JA+J) = A(IA+I)+SIN45*(A(IB+I)-A(ID+I))
              C(JB+J) = A(IA+I)-SIN45*(A(IB+I)-A(ID+I))
              D(JA+J) = -A(IC+I)-SIN45*(A(IB+I)+A(ID+I))
              D(JB+J) = A(IC+I)-SIN45*(A(IB+I)+A(ID+I))
              I = I+INC3
              J = J+INC4
            END DO
            IBASE = IBASE+INC1
            JBASE = JBASE+INC2
          END DO
          FUPDATE = .true.
        END IF
      END IF
      IF (FUPDATE .eqv. .false.) THEN
        Z = 1.0/N
        DO L=1,LA
          I = IBASE
          J = JBASE
      !CDIR NODEP
          DO IJK=1,LOT
            C(JA+J) = Z*((A(IA+I)+A(IC+I))+(A(IB+I)+A(ID+I)))
            C(JC+J) = Z*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I)))
            C(JB+J) = Z*(A(IA+I)-A(IC+I))
            D(JB+J) = Z*(A(ID+I)-A(IB+I))
            I = I+INC3
            J = J+INC4
          END DO
          IBASE = IBASE+INC1
          JBASE = JBASE+INC2
        END DO
        FUPDATE = .true.
      END IF

    !--
    CASE (4)
      ! Coding for factor 5

      IA = 1
      IB = IA+IINK
      IC = IB+IINK
      ID = IC+IINK
      IE = ID+IINK
      JA = 1
      JB = JA+(2*M-LA)*INC2
      JC = JB+2*M*INC2
      JD = JC
      JE = JB

      IF (LA /= M) THEN
        DO L=1,LA
          I = IBASE
          J = JBASE
    !CDIR NODEP
          DO IJK=1,LOT
            A1 = A(IB+I)+A(IE+I)
            A3 = A(IB+I)-A(IE+I)
            A2 = A(IC+I)+A(ID+I)
            A4 = A(IC+I)-A(ID+I)
            A5 = A(IA+I)-0.25*(A1+A2)
            A6 = QRT5*(A1-A2)
            C(JA+J) = A(IA+I)+(A1+A2)
            C(JB+J) = A5+A6
            C(JC+J) = A5-A6
            D(JB+J) = -SIN72*A3-SIN36*A4
            D(JC+J) = -SIN36*A3+SIN72*A4
            I = I+INC3
            J = J+INC4
          END DO
          IBASE = IBASE+INC1
          JBASE = JBASE+INC2
        END DO
        JA = JA+JINK
        JINK = 2*JINK
        JB = JB+JINK
        JC = JC+JINK
        JD = JD-JINK
        JE = JE-JINK
        IBASE = IBASE+IJUMP
        IJUMP = 2*IJUMP+IINK
        IF (JB /= JD) THEN
          DO K=LA,KSTOP,LA
            KB = K+K
            KC = KB+KB
            KD = KC+KB
            KE = KD+KB
            C1 = TRIGS(KB+1)
            S1 = TRIGS(KB+2)
            C2 = TRIGS(KC+1)
            S2 = TRIGS(KC+2)
            C3 = TRIGS(KD+1)
            S3 = TRIGS(KD+2)
            C4 = TRIGS(KE+1)
            S4 = TRIGS(KE+2)
            JBASE = 0
            DO L=1,LA
              I = IBASE
              J = JBASE
    !CDIR NODEP
              DO IJK=1,LOT
                A1 = (C1*A(IB+I)+S1*B(IB+I))+(C4*A(IE+I)+S4*B(IE+I))
                A3 = (C1*A(IB+I)+S1*B(IB+I))-(C4*A(IE+I)+S4*B(IE+I))
                A2 = (C2*A(IC+I)+S2*B(IC+I))+(C3*A(ID+I)+S3*B(ID+I))
                A4 = (C2*A(IC+I)+S2*B(IC+I))-(C3*A(ID+I)+S3*B(ID+I))
                B1 = (C1*B(IB+I)-S1*A(IB+I))+(C4*B(IE+I)-S4*A(IE+I))
                B3 = (C1*B(IB+I)-S1*A(IB+I))-(C4*B(IE+I)-S4*A(IE+I))
                B2 = (C2*B(IC+I)-S2*A(IC+I))+(C3*B(ID+I)-S3*A(ID+I))
                B4 = (C2*B(IC+I)-S2*A(IC+I))-(C3*B(ID+I)-S3*A(ID+I))
                A5 = A(IA+I)-0.25*(A1+A2)
                A6 = QRT5*(A1-A2)
                B5 = B(IA+I)-0.25*(B1+B2)
                B6 = QRT5*(B1-B2)
                A10 = A5+A6
                A20 = A5-A6
                B10 = B5+B6
                B20 = B5-B6
                A11 = SIN72*B3+SIN36*B4
                A21 = SIN36*B3-SIN72*B4
                B11 = SIN72*A3+SIN36*A4
                B21 = SIN36*A3-SIN72*A4
                C(JA+J) = A(IA+I)+(A1+A2)
                C(JB+J) = A10+A11
                C(JE+J) = A10-A11
                C(JC+J) = A20+A21
                C(JD+J) = A20-A21
                D(JA+J) = B(IA+I)+(B1+B2)
                D(JB+J) = B10-B11
                D(JE+J) = -(B10+B11)
                D(JC+J) = B20-B21
                D(JD+J) = -(B20+B21)
                I = I+INC3
                J = J+INC4
              END DO
              IBASE = IBASE+INC1
              JBASE = JBASE+INC2
            END DO
            IBASE = IBASE+IJUMP
            JA = JA+JINK
            JB = JB+JINK
            JC = JC+JINK
            JD = JD-JINK
            JE = JE-JINK
          END DO
          IF (JB > JD) THEN
            FUPDATE = .true.
          END IF
        END IF
        IF (FUPDATE .eqv. .false.) THEN
          JBASE = 0
          DO L=1,LA
            I = IBASE
            J = JBASE
      !CDIR NODEP
            DO IJK=1,LOT
              A1 = A(IB+I)+A(IE+I)
              A3 = A(IB+I)-A(IE+I)
              A2 = A(IC+I)+A(ID+I)
              A4 = A(IC+I)-A(ID+I)
              A5 = A(IA+I)+0.25*(A3-A4)
              A6 = QRT5*(A3+A4)
              C(JA+J) = A5+A6
              C(JB+J) = A5-A6
              C(JC+J) = A(IA+I)-(A3-A4)
              D(JA+J) = -SIN36*A1-SIN72*A2
              D(JB+J) = -SIN72*A1+SIN36*A2
              I = I+INC3
              J = J+INC4
            END DO
            IBASE = IBASE+INC1
            JBASE = JBASE+INC2
          END DO
          FUPDATE = .true.
        END IF
      END IF
      IF (FUPDATE .eqv. .false.) THEN
        Z = 1.0/N
        ZQRT5 = Z*QRT5
        ZSIN36 = Z*SIN36
        ZSIN72 = Z*SIN72
        DO L=1,LA
          I = IBASE
          J = JBASE
      !CDIR NODEP
          DO IJK=1,LOT
            A1 = A(IB+I)+A(IE+I)
            A3 = A(IB+I)-A(IE+I)
            A2 = A(IC+I)+A(ID+I)
            A4 = A(IC+I)-A(ID+I)
            A5 = Z*(A(IA+I)-0.25*(A1+A2))
            A6 = ZQRT5*(A1-A2)
            C(JA+J) = Z*(A(IA+I)+(A1+A2))
            C(JB+J) = A5+A6
            C(JC+J) = A5-A6
            D(JB+J) = -ZSIN72*A3-ZSIN36*A4
            D(JC+J) = -ZSIN36*A3+ZSIN72*A4
            I = I+INC3
            J = J+INC4
          END DO
          IBASE = IBASE+INC1
          JBASE = JBASE+INC2
        END DO
        FUPDATE = .true.
      END IF

    !--
    CASE (5)
      ! Coding for factor 6

      IA = 1
      IB = IA+IINK
      IC = IB+IINK
      ID = IC+IINK
      IE = ID+IINK
      IF = IE+IINK
      JA = 1
      JB = JA+(2*M-LA)*INC2
      JC = JB+2*M*INC2
      JD = JC+2*M*INC2
      JE = JC
      JF = JB

      IF (LA /= M) THEN
        DO L=1,LA
          I = IBASE
          J = JBASE
    !CDIR NODEP
          DO IJK=1,LOT
            A11 = (A(IC+I)+A(IF+I))+(A(IB+I)+A(IE+I))
            C(JA+J) = (A(IA+I)+A(ID+I))+A11
            C(JC+J) = (A(IA+I)+A(ID+I)-0.5*A11)
            D(JC+J) = SIN60*((A(IC+I)+A(IF+I))-(A(IB+I)+A(IE+I)))
            A11 = (A(IC+I)-A(IF+I))+(A(IE+I)-A(IB+I))
            C(JB+J) = (A(IA+I)-A(ID+I))-0.5*A11
            D(JB+J) = SIN60*((A(IE+I)-A(IB+I))-(A(IC+I)-A(IF+I)))
            C(JD+J) = (A(IA+I)-A(ID+I))+A11
            I = I+INC3
            J = J+INC4
          END DO
          IBASE = IBASE+INC1
          JBASE = JBASE+INC2
        END DO
        JA = JA+JINK
        JINK = 2*JINK
        JB = JB+JINK
        JC = JC+JINK
        JD = JD-JINK
        JE = JE-JINK
        JF = JF-JINK
        IBASE = IBASE+IJUMP
        IJUMP = 2*IJUMP+IINK
        IF (JC /= JD) THEN
          DO K=LA,KSTOP,LA
            KB = K+K
            KC = KB+KB
            KD = KC+KB
            KE = KD+KB
            KF = KE+KB
            C1 = TRIGS(KB+1)
            S1 = TRIGS(KB+2)
            C2 = TRIGS(KC+1)
            S2 = TRIGS(KC+2)
            C3 = TRIGS(KD+1)
            S3 = TRIGS(KD+2)
            C4 = TRIGS(KE+1)
            S4 = TRIGS(KE+2)
            C5 = TRIGS(KF+1)
            S5 = TRIGS(KF+2)
            JBASE = 0
            DO L=1,LA
              I = IBASE
              J = JBASE
    !CDIR NODEP
              DO IJK=1,LOT
                A1 = C1*A(IB+I)+S1*B(IB+I)
                B1 = C1*B(IB+I)-S1*A(IB+I)
                A2 = C2*A(IC+I)+S2*B(IC+I)
                B2 = C2*B(IC+I)-S2*A(IC+I)
                A3 = C3*A(ID+I)+S3*B(ID+I)
                B3 = C3*B(ID+I)-S3*A(ID+I)
                A4 = C4*A(IE+I)+S4*B(IE+I)
                B4 = C4*B(IE+I)-S4*A(IE+I)
                A5 = C5*A(IF+I)+S5*B(IF+I)
                B5 = C5*B(IF+I)-S5*A(IF+I)
                A11 = (A2+A5)+(A1+A4)
                A20 = (A(IA+I)+A3)-0.5*A11
                A21 = SIN60*((A2+A5)-(A1+A4))
                B11 = (B2+B5)+(B1+B4)
                B20 = (B(IA+I)+B3)-0.5*B11
                B21 = SIN60*((B2+B5)-(B1+B4))
                C(JA+J) = (A(IA+I)+A3)+A11
                D(JA+J) = (B(IA+I)+B3)+B11
                C(JC+J) = A20-B21
                D(JC+J) = A21+B20
                C(JE+J) = A20+B21
                D(JE+J) = A21-B20
                A11 = (A2-A5)+(A4-A1)
                A20 = (A(IA+I)-A3)-0.5*A11
                A21 = SIN60*((A4-A1)-(A2-A5))
                B11 = (B5-B2)-(B4-B1)
                B20 = (B3-B(IA+I))-0.5*B11
                B21 = SIN60*((B5-B2)+(B4-B1))
                C(JB+J) = A20-B21
                D(JB+J) = A21-B20
                C(JD+J) = A11+(A(IA+I)-A3)
                D(JD+J) = B11+(B3-B(IA+I))
                C(JF+J) = A20+B21
                D(JF+J) = A21+B20
                I = I+INC3
                J = J+INC4
              END DO
              IBASE = IBASE+INC1
              JBASE = JBASE+INC2
            END DO
            IBASE = IBASE+IJUMP
            JA = JA+JINK
            JB = JB+JINK
            JC = JC+JINK
            JD = JD-JINK
            JE = JE-JINK
            JF = JF-JINK
          END DO
          IF (JC > JD) THEN
            FUPDATE = .true.
          END IF
        END IF
        IF (FUPDATE .eqv. .false.) THEN
          JBASE = 0
          DO L=1,LA
            I = IBASE
            J = JBASE
      !CDIR NODEP
            DO IJK=1,LOT
              C(JA+J) = (A(IA+I)+0.5*(A(IC+I)-A(IE+I)))+ SIN60*(A(IB+I)-A(IF+I))
              D(JA+J) = -(A(ID+I)+0.5*(A(IB+I)+A(IF+I)))-SIN60*(A(IC+I)+A(IE+I))
              C(JB+J) = A(IA+I)-(A(IC+I)-A(IE+I))
              D(JB+J) = A(ID+I)-(A(IB+I)+A(IF+I))
              C(JC+J) = (A(IA+I)+0.5*(A(IC+I)-A(IE+I)))-SIN60*(A(IB+I)-A(IF+I))
              D(JC+J) = -(A(ID+I)+0.5*(A(IB+I)+A(IF+I)))+SIN60*(A(IC+I)+A(IE+I))
              I = I+INC3
              J = J+INC4
            END DO
            IBASE = IBASE+INC1
            JBASE = JBASE+INC2
          END DO
          FUPDATE = .true.
        END IF
      END IF
      IF (FUPDATE .eqv. .false.) THEN
        Z = 1.0/N
        ZSIN60 = Z*SIN60
        DO L=1,LA
          I = IBASE
          J = JBASE
      !CDIR NODEP
          DO IJK=1,LOT
            A11 = (A(IC+I)+A(IF+I))+(A(IB+I)+A(IE+I))
            C(JA+J) = Z*((A(IA+I)+A(ID+I))+A11)
            C(JC+J) = Z*((A(IA+I)+A(ID+I))-0.5*A11)
            D(JC+J) = ZSIN60*((A(IC+I)+A(IF+I))-(A(IB+I)+A(IE+I)))
            A11 = (A(IC+I)-A(IF+I))+(A(IE+I)-A(IB+I))
            C(JB+J) = Z*((A(IA+I)-A(ID+I))-0.5*A11)
            D(JB+J) = ZSIN60*((A(IE+I)-A(IB+I))-(A(IC+I)-A(IF+I)))
            C(JD+J) = Z*((A(IA+I)-A(ID+I))+A11)
            I = I+INC3
            J = J+INC4
          END DO
          IBASE = IBASE+INC1
          JBASE = JBASE+INC2
        END DO
        FUPDATE = .true.
      END IF

    !--
    CASE (6,7)
      ! Coding for factor 8

      IF (LA /= M) THEN
        ! error: IFAC only catered for if LA=N/IFAC
        IBAD = 3
        FUPDATE = .false.
      ELSE
        IA = 1
        IB = IA+IINK
        IC = IB+IINK
        ID = IC+IINK
        IE = ID+IINK
        IF = IE+IINK
        IG = IF+IINK
        IH = IG+IINK
        JA = 1
        JB = JA+LA*INC2
        JC = JB+2*M*INC2
        JD = JC+2*M*INC2
        JE = JD+2*M*INC2
        Z = 1.0/N
        ZSIN45 = Z*SQRT(0.5)

        DO L=1,LA
          I = IBASE
          J = JBASE
      !CDIR NODEP
          DO IJK=1,LOT
            C(JA+J) = Z*(((A(IA+I)+A(IE+I))+(A(IC+I)+A(IG+I)))+ &
                       ((A(ID+I)+A(IH+I))+(A(IB+I)+A(IF+I))))
            C(JE+J) = Z*(((A(IA+I)+A(IE+I))+(A(IC+I)+A(IG+I)))- &
                       ((A(ID+I)+A(IH+I))+(A(IB+I)+A(IF+I))))
            C(JC+J) = Z*((A(IA+I)+A(IE+I))-(A(IC+I)+A(IG+I)))
            D(JC+J) = Z*((A(ID+I)+A(IH+I))-(A(IB+I)+A(IF+I)))
            C(JB+J) = Z*(A(IA+I)-A(IE+I))+ZSIN45*((A(IH+I)-A(ID+I))-(A(IF+I)-A(IB+I)))
            C(JD+J) = Z*(A(IA+I)-A(IE+I))-ZSIN45*((A(IH+I)-A(ID+I))-(A(IF+I)-A(IB+I)))
            D(JB+J) = ZSIN45*((A(IH+I)-A(ID+I))+(A(IF+I)-A(IB+I)))+Z*(A(IG+I)-A(IC+I))
            D(JD+J) = ZSIN45*((A(IH+I)-A(ID+I))+(A(IF+I)-A(IB+I)))-Z*(A(IG+I)-A(IC+I))
            I = I+INC3
            J = J+INC4
          END DO
          IBASE = IBASE+INC1
          JBASE = JBASE+INC2
        END DO
        FUPDATE = .true.
      END IF

    !--
    CASE DEFAULT
      FUPDATE = .false.

    END SELECT

    IF (FUPDATE .eqv. .true.) THEN
      ! final update ...
      ! pass completed without error
      IBAD = 0
    END IF

  END IF

  IERR = IBAD


  end subroutine fft_qpassm_f


  subroutine fft_rpassm_f(A, B, C, D, TRIGS, INC1, INC2, INC3, INC4, LOT, N, IFAC, LA, IERR)

  !
  ! the function is used to carry out data processing -- i.e. to perform one pass
  ! through data as part of multiple FFTs (Temperton method)
  !
  !
  ! A is first real input vector
  !   Equivalence B(1) with A (LA*INC1+1)
  ! C is first real output vector
  !   Equivalence D(1) with C(IFAC*LA*INC2+1)
  ! TRIGS is a precalculated list of sines and cosines
  ! INC1 is the addressing increment for A
  ! INC2 is the addressing increment for C
  ! INC3 is the increment between input vectors A
  ! INC4 is the increment between output vectors C
  ! LOT is the number of vectors
  ! N is the length of the vectors
  ! IFAC is the current factor of N
  ! LA is the product of previous factors
  ! IERR is an error indicator
  !   0 - Pass completed without error
  !   2 - IFAC not catered for
  !   3 - IFAC only catered for if LA=N/IFAC
  !

  !--

  IMPLICIT NONE


  INTEGER, INTENT(IN)  :: N
  INTEGER, INTENT(IN)  :: INC1,INC2,INC3,INC4,LOT,IFAC,LA
  INTEGER, INTENT(OUT) :: IERR
  REAL, INTENT(IN)     :: A(*)
  REAL, INTENT(IN)     :: B(*)
  REAL, INTENT(INOUT)  :: C(*)
  REAL, INTENT(INOUT)  :: D(*)
  REAL, INTENT(IN)     :: TRIGS(N)

  INTEGER :: JF, KF, IF, IE, JE, KE, ID, JD, KD, KC
  INTEGER :: I, IJK, K, KB, IC, JC, IA, IB, JA, JB, L, J, IGO, M
  INTEGER :: IINK, JINK, KSTOP, IBAD, IBASE, JBASE, JUMP, JG, JH
  REAL    :: C5, S5, B21, B11, A11, A21, SSIN60, QQRT5
  REAL    :: A10, A20, B10, B20, C4, S4
  REAL    :: SIN45, C3, S3
  REAL    :: C2, C1, S1
  REAL    :: S2, SSIN36, SSIN72, SSIN45
  REAL, PARAMETER :: SIN36 = 0.587785252292473
  REAL, PARAMETER :: SIN60 = 0.866025403784437
  REAL, PARAMETER :: SIN72 = 0.951056516295154
  REAL, PARAMETER :: QRT5 = 0.559016994374947

  LOGICAL :: FUPDATE

  ! if this flag is set to 'true', set 'IBAD=0' to indicate
  ! that the pass is completed without error
  FUPDATE = .false.

  M = N/IFAC
  IINK = LA*INC1
  JINK = LA*INC2
  JUMP = (IFAC-1)*JINK
  KSTOP = (N-IFAC)/(2*IFAC)

  IBASE = 0
  JBASE = 0
  IGO = IFAC-1
  ! error: IFAC not catered for (default case)
  IBAD = 2

  SELECT CASE (IGO)

  !--
  CASE (1)
    ! Coding for factor 2

    IA = 1
    IB = IA+(2*M-LA)*INC1
    JA = 1
    JB = JA+JINK

    IF (LA /= M) THEN
      DO L=1,LA
        I = IBASE
        J = JBASE
        DO IJK=1,LOT
          C(JA+J) = A(IA+I)+A(IB+I)
          C(JB+J) = A(IA+I)-A(IB+I)
          I = I+INC3
          J = J+INC4
        END DO
        IBASE = IBASE+INC1
        JBASE = JBASE+INC2
      END DO
      IA = IA+IINK
      IINK = 2*IINK
      IB = IB-IINK
      IBASE = 0
      JBASE = JBASE+JUMP
      JUMP = 2*JUMP+JINK
      IF (IA /= IB) THEN
        DO K=LA,KSTOP,LA
          KB = K+K
          C1 = TRIGS(KB+1)
          S1 = TRIGS(KB+2)
          IBASE = 0
          DO L=1,LA
            I = IBASE
            J = JBASE
            DO IJK=1,LOT
              C(JA+J) = A(IA+I)+A(IB+I)
              D(JA+J) = B(IA+I)-B(IB+I)
              C(JB+J) = C1*(A(IA+I)-A(IB+I))-S1*(B(IA+I)+B(IB+I))
              D(JB+J) = S1*(A(IA+I)-A(IB+I))+C1*(B(IA+I)+B(IB+I))
              I = I+INC3
              J = J+INC4
            END DO
            IBASE = IBASE+INC1
            JBASE = JBASE+INC2
          END DO
          IA = IA+IINK
          IB = IB-IINK
          JBASE = JBASE+JUMP
        END DO
        IF (IA > IB) THEN
          FUPDATE = .true.
        END IF
      END IF
      IF (FUPDATE .eqv. .false.) THEN
        IBASE = 0
        DO L=1,LA
          I = IBASE
          J = JBASE
          DO IJK=1,LOT
            C(JA+J) = A(IA+I)
            C(JB+J) = -B(IA+I)
            I = I+INC3
            J = J+INC4
          END DO
          IBASE = IBASE+INC1
          JBASE = JBASE+INC2
        END DO
        FUPDATE = .true.
      END IF
    END IF
    IF (FUPDATE .eqv. .false.) THEN
      DO  L=1,LA
        I = IBASE
        J = JBASE
        DO IJK=1,LOT
          C(JA+J) = 2.0*(A(IA+I)+A(IB+I))
          C(JB+J) = 2.0*(A(IA+I)-A(IB+I))
          I = I+INC3
          J = J+INC4
        END DO
        IBASE = IBASE+INC1
        JBASE = JBASE+INC2
      END DO
      FUPDATE = .true.
    END IF

  !--
  CASE (2)
    ! Coding for factor 3

    IA = 1
    IB = IA+(2*M-LA)*INC1
    IC = IB
    JA = 1
    JB = JA+JINK
    JC = JB+JINK

    IF (LA /= M) THEN
      DO L=1,LA
        I = IBASE
        J = JBASE
        DO IJK=1,LOT
          C(JA+J) = A(IA+I)+A(IB+I)
          C(JB+J) = (A(IA+I)-0.5*A(IB+I))-(SIN60*(B(IB+I)))
          C(JC+J) = (A(IA+I)-0.5*A(IB+I))+(SIN60*(B(IB+I)))
          I = I+INC3
          J = J+INC4
        END DO
        IBASE = IBASE+INC1
        JBASE = JBASE+INC2
      END DO
      IA = IA+IINK
      IINK = 2*IINK
      IB = IB+IINK
      IC = IC-IINK
      JBASE = JBASE+JUMP
      JUMP = 2*JUMP+JINK
      IF (IA /= IC) THEN
        DO K=LA,KSTOP,LA
          KB = K+K
          KC = KB+KB
          C1 = TRIGS(KB+1)
          S1 = TRIGS(KB+2)
          C2 = TRIGS(KC+1)
          S2 = TRIGS(KC+2)
          IBASE = 0
          DO L=1,LA
            I = IBASE
            J = JBASE
            DO IJK=1,LOT
              C(JA+J) = A(IA+I)+(A(IB+I)+A(IC+I))
              D(JA+J) = B(IA+I)+(B(IB+I)-B(IC+I))
              C(JB+J) = C1*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))-(SIN60*(B(IB+I)+B(IC+I)))) &
                       -S1*((B(IA+I)-0.5*(B(IB+I)-B(IC+I)))+(SIN60*(A(IB+I)-A(IC+I))))
              D(JB+J) = S1*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))-(SIN60*(B(IB+I)+B(IC+I)))) &
                       +C1*((B(IA+I)-0.5*(B(IB+I)-B(IC+I)))+(SIN60*(A(IB+I)-A(IC+I))))
              C(JC+J) = C2*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))+(SIN60*(B(IB+I)+B(IC+I)))) &
                       -S2*((B(IA+I)-0.5*(B(IB+I)-B(IC+I)))-(SIN60*(A(IB+I)-A(IC+I))))
              D(JC+J) = S2*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))+(SIN60*(B(IB+I)+B(IC+I)))) &
                       +C2*((B(IA+I)-0.5*(B(IB+I)-B(IC+I)))-(SIN60*(A(IB+I)-A(IC+I))))
              I = I+INC3
              J = J+INC4
            END DO
            IBASE = IBASE+INC1
            JBASE = JBASE+INC2
          END DO
          IA = IA+IINK
          IB = IB+IINK
          IC = IC-IINK
          JBASE = JBASE+JUMP
        END DO
        IF (IA > IC) THEN
          FUPDATE = .true.
        END IF
      END IF
      IF (FUPDATE .eqv. .false.) THEN
        IBASE = 0
        DO L=1,LA
          I = IBASE
          J = JBASE
          DO IJK=1,LOT
            C(JA+J) = A(IA+I)+A(IB+I)
            C(JB+J) = (0.5*A(IA+I)-A(IB+I))-(SIN60*B(IA+I))
            C(JC+J) = -(0.5*A(IA+I)-A(IB+I))-(SIN60*B(IA+I))
            I = I+INC3
            J = J+INC4
          END DO
          IBASE = IBASE+INC1
          JBASE = JBASE+INC2
        END DO
        FUPDATE = .true.
      END IF
    END IF
    IF (FUPDATE .eqv. .false.) THEN
      SSIN60 = 2.0*SIN60
      DO L=1,LA
        I = IBASE
        J = JBASE
        DO IJK=1,LOT
          C(JA+J) = 2.0*(A(IA+I)+A(IB+I))
          C(JB+J) = (2.0*A(IA+I)-A(IB+I))-(SSIN60*B(IB+I))
          C(JC+J) = (2.0*A(IA+I)-A(IB+I))+(SSIN60*B(IB+I))
          I = I+INC3
          J = J+INC4
        END DO
        IBASE = IBASE+INC1
        JBASE = JBASE+INC2
      END DO
      FUPDATE = .true.
    END IF

  !--
  CASE (3)
    ! Coding for factor 4

    IA = 1
    IB = IA+(2*M-LA)*INC1
    IC = IB+2*M*INC1
    ID = IB
    JA = 1
    JB = JA+JINK
    JC = JB+JINK
    JD = JC+JINK

    IF (LA /= M) THEN
      DO L=1,LA
        I = IBASE
        J = JBASE
        DO IJK=1,LOT
          C(JA+J) = (A(IA+I)+A(IC+I))+A(IB+I)
          C(JB+J) = (A(IA+I)-A(IC+I))-B(IB+I)
          C(JC+J) = (A(IA+I)+A(IC+I))-A(IB+I)
          C(JD+J) = (A(IA+I)-A(IC+I))+B(IB+I)
          I = I+INC3
          J = J+INC4
        END DO
        IBASE = IBASE+INC1
        JBASE = JBASE+INC2
      END DO
      IA = IA+IINK
      IINK = 2*IINK
      IB = IB+IINK
      IC = IC-IINK
      ID = ID-IINK
      JBASE = JBASE+JUMP
      JUMP = 2*JUMP+JINK
      IF (IB /= IC) THEN
        DO K=LA,KSTOP,LA
          KB = K+K
          KC = KB+KB
          KD = KC+KB
          C1 = TRIGS(KB+1)
          S1 = TRIGS(KB+2)
          C2 = TRIGS(KC+1)
          S2 = TRIGS(KC+2)
          C3 = TRIGS(KD+1)
          S3 = TRIGS(KD+2)
          IBASE = 0
          DO L=1,LA
            I = IBASE
            J = JBASE
            DO IJK=1,LOT
              C(JA+J) = (A(IA+I)+A(IC+I))+(A(IB+I)+A(ID+I))
              D(JA+J) = (B(IA+I)-B(IC+I))+(B(IB+I)-B(ID+I))
              C(JC+J) = C2*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I))) &
                       -S2*((B(IA+I)-B(IC+I))-(B(IB+I)-B(ID+I)))
              D(JC+J) = S2*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I))) &
                       +C2*((B(IA+I)-B(IC+I))-(B(IB+I)-B(ID+I)))
              C(JB+J) = C1*((A(IA+I)-A(IC+I))-(B(IB+I)+B(ID+I))) &
                       -S1*((B(IA+I)+B(IC+I))+(A(IB+I)-A(ID+I)))
              D(JB+J) = S1*((A(IA+I)-A(IC+I))-(B(IB+I)+B(ID+I))) &
                       +C1*((B(IA+I)+B(IC+I))+(A(IB+I)-A(ID+I)))
              C(JD+J) = C3*((A(IA+I)-A(IC+I))+(B(IB+I)+B(ID+I))) &
                       -S3*((B(IA+I)+B(IC+I))-(A(IB+I)-A(ID+I)))
              D(JD+J) = S3*((A(IA+I)-A(IC+I))+(B(IB+I)+B(ID+I))) &
                       +C3*((B(IA+I)+B(IC+I))-(A(IB+I)-A(ID+I)))
              I = I+INC3
              J = J+INC4
            END DO
            IBASE = IBASE+INC1
            JBASE = JBASE+INC2
          END DO
          IA = IA+IINK
          IB = IB+IINK
          IC = IC-IINK
          ID = ID-IINK
          JBASE = JBASE+JUMP
        END DO
        IF (IB > IC) THEN
          FUPDATE = .true.
        END IF
      END IF
      IF (FUPDATE .eqv. .false.) THEN
        IBASE = 0
        SIN45 = SQRT(0.5)
        DO L=1,LA
          I = IBASE
          J = JBASE
          DO IJK=1,LOT
            C(JA+J) = A(IA+I)+A(IB+I)
            C(JB+J) = SIN45*((A(IA+I)-A(IB+I))-(B(IA+I)+B(IB+I)))
            C(JC+J) = B(IB+I)-B(IA+I)
            C(JD+J) = -SIN45*((A(IA+I)-A(IB+I))+(B(IA+I)+B(IB+I)))
            I = I+INC3
            J = J+INC4
          END DO
          IBASE = IBASE+INC1
          JBASE = JBASE+INC2
        END DO
        FUPDATE = .true.
      END IF
    END IF
    IF (FUPDATE .eqv. .false.) THEN
      DO L=1,LA
        I = IBASE
        J = JBASE
        DO IJK=1,LOT
          C(JA+J) = 2.0*((A(IA+I)+A(IC+I))+A(IB+I))
          C(JB+J) = 2.0*((A(IA+I)-A(IC+I))-B(IB+I))
          C(JC+J) = 2.0*((A(IA+I)+A(IC+I))-A(IB+I))
          C(JD+J) = 2.0*((A(IA+I)-A(IC+I))+B(IB+I))
          I = I+INC3
          J = J+INC4
        END DO
        IBASE = IBASE+INC1
        JBASE = JBASE+INC2
      END DO
      FUPDATE = .true.
    END IF

  !--
  CASE (4)
    ! Coding for factor 5

    IA = 1
    IB = IA+(2*M-LA)*INC1
    IC = IB+2*M*INC1
    ID = IC
    IE = IB
    JA = 1
    JB = JA+JINK
    JC = JB+JINK
    JD = JC+JINK
    JE = JD+JINK

    IF (LA /= M) THEN
      DO L=1,LA
        I = IBASE
        J = JBASE
        DO IJK=1,LOT
          C(JA+J) = A(IA+I)+(A(IB+I)+A(IC+I))
          C(JB+J) = ((A(IA+I)-0.25*(A(IB+I)+A(IC+I)))+QRT5*(A(IB+I)-A(IC+I))) &
              -(SIN72*B(IB+I)+SIN36*B(IC+I))
          C(JC+J) = ((A(IA+I)-0.25*(A(IB+I)+A(IC+I)))-QRT5*(A(IB+I)-A(IC+I))) &
              -(SIN36*B(IB+I)-SIN72*B(IC+I))
          C(JD+J) = ((A(IA+I)-0.25*(A(IB+I)+A(IC+I)))-QRT5*(A(IB+I)-A(IC+I))) &
              +(SIN36*B(IB+I)-SIN72*B(IC+I))
          C(JE+J) = ((A(IA+I)-0.25*(A(IB+I)+A(IC+I)))+QRT5*(A(IB+I)-A(IC+I))) &
              +(SIN72*B(IB+I)+SIN36*B(IC+I))
          I = I+INC3
          J = J+INC4
        END DO
        IBASE = IBASE+INC1
        JBASE = JBASE+INC2
      END DO
      IA = IA+IINK
      IINK = 2*IINK
      IB = IB+IINK
      IC = IC+IINK
      ID = ID-IINK
      IE = IE-IINK
      JBASE = JBASE+JUMP
      JUMP = 2*JUMP+JINK
      IF (IB /= ID) THEN
        DO K=LA,KSTOP,LA
          KB = K+K
          KC = KB+KB
          KD = KC+KB
          KE = KD+KB
          C1 = TRIGS(KB+1)
          S1 = TRIGS(KB+2)
          C2 = TRIGS(KC+1)
          S2 = TRIGS(KC+2)
          C3 = TRIGS(KD+1)
          S3 = TRIGS(KD+2)
          C4 = TRIGS(KE+1)
          S4 = TRIGS(KE+2)
          IBASE = 0
          DO L=1,LA
            I = IBASE
            J = JBASE
            DO IJK=1,LOT

              A10 = (A(IA+I)-0.25*((A(IB+I)+A(IE+I))+(A(IC+I)+A(ID+I)))) &
                  +QRT5*((A(IB+I)+A(IE+I))-(A(IC+I)+A(ID+I)))
              A20 = (A(IA+I)-0.25*((A(IB+I)+A(IE+I))+(A(IC+I)+A(ID+I)))) &
                  -QRT5*((A(IB+I)+A(IE+I))-(A(IC+I)+A(ID+I)))
              B10 = (B(IA+I)-0.25*((B(IB+I)-B(IE+I))+(B(IC+I)-B(ID+I)))) &
                  +QRT5*((B(IB+I)-B(IE+I))-(B(IC+I)-B(ID+I)))
              B20 = (B(IA+I)-0.25*((B(IB+I)-B(IE+I))+(B(IC+I)-B(ID+I)))) &
                  -QRT5*((B(IB+I)-B(IE+I))-(B(IC+I)-B(ID+I)))
              A11 = SIN72*(B(IB+I)+B(IE+I))+SIN36*(B(IC+I)+B(ID+I))
              A21 = SIN36*(B(IB+I)+B(IE+I))-SIN72*(B(IC+I)+B(ID+I))
              B11 = SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))
              B21 = SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))

              C(JA+J) = A(IA+I)+((A(IB+I)+A(IE+I))+(A(IC+I)+A(ID+I)))
              D(JA+J) = B(IA+I)+((B(IB+I)-B(IE+I))+(B(IC+I)-B(ID+I)))
              C(JB+J) = C1*(A10-A11)-S1*(B10+B11)
              D(JB+J) = S1*(A10-A11)+C1*(B10+B11)
              C(JE+J) = C4*(A10+A11)-S4*(B10-B11)
              D(JE+J) = S4*(A10+A11)+C4*(B10-B11)
              C(JC+J) = C2*(A20-A21)-S2*(B20+B21)
              D(JC+J) = S2*(A20-A21)+C2*(B20+B21)
              C(JD+J) = C3*(A20+A21)-S3*(B20-B21)
              D(JD+J) = S3*(A20+A21)+C3*(B20-B21)

              I = I+INC3
              J = J+INC4
            END DO
            IBASE = IBASE+INC1
            JBASE = JBASE+INC2
          END DO
          IA = IA+IINK
          IB = IB+IINK
          IC = IC+IINK
          ID = ID-IINK
          IE = IE-IINK
          JBASE = JBASE+JUMP
        END DO
        IF (IB > ID) THEN
          FUPDATE = .true.
        END IF
      END IF
      IF (FUPDATE .eqv. .false.) THEN
        IBASE = 0
        DO L=1,LA
          I = IBASE
          J = JBASE
          DO IJK=1,LOT
            C(JA+J) = (A(IA+I)+A(IB+I))+A(IC+I)
            C(JB+J) = (QRT5*(A(IA+I)-A(IB+I))+(0.25*(A(IA+I)+A(IB+I))-A(IC+I))) &
                -(SIN36*B(IA+I)+SIN72*B(IB+I))
            C(JE+J) = -(QRT5*(A(IA+I)-A(IB+I))+(0.25*(A(IA+I)+A(IB+I))-A(IC+I))) &
                -(SIN36*B(IA+I)+SIN72*B(IB+I))
            C(JC+J) = (QRT5*(A(IA+I)-A(IB+I))-(0.25*(A(IA+I)+A(IB+I))-A(IC+I))) &
                -(SIN72*B(IA+I)-SIN36*B(IB+I))
            C(JD+J) = -(QRT5*(A(IA+I)-A(IB+I))-(0.25*(A(IA+I)+A(IB+I))-A(IC+I))) &
                -(SIN72*B(IA+I)-SIN36*B(IB+I))
            I = I+INC3
            J = J+INC4
          END DO
          IBASE = IBASE+INC1
          JBASE = JBASE+INC2
        END DO
        FUPDATE = .true.
      END IF
    END IF
    IF (FUPDATE .eqv. .false.) THEN
      QQRT5 = 2.0*QRT5
      SSIN36 = 2.0*SIN36
      SSIN72 = 2.0*SIN72
      DO L=1,LA
        I = IBASE
        J = JBASE
        DO IJK=1,LOT
          C(JA+J) = 2.0*(A(IA+I)+(A(IB+I)+A(IC+I)))
          C(JB+J) = (2.0*(A(IA+I)-0.25*(A(IB+I)+A(IC+I))) &
              +QQRT5*(A(IB+I)-A(IC+I)))-(SSIN72*B(IB+I)+SSIN36*B(IC+I))
          C(JC+J) = (2.0*(A(IA+I)-0.25*(A(IB+I)+A(IC+I))) &
              -QQRT5*(A(IB+I)-A(IC+I)))-(SSIN36*B(IB+I)-SSIN72*B(IC+I))
          C(JD+J) = (2.0*(A(IA+I)-0.25*(A(IB+I)+A(IC+I))) &
              -QQRT5*(A(IB+I)-A(IC+I)))+(SSIN36*B(IB+I)-SSIN72*B(IC+I))
          C(JE+J) = (2.0*(A(IA+I)-0.25*(A(IB+I)+A(IC+I))) &
              +QQRT5*(A(IB+I)-A(IC+I)))+(SSIN72*B(IB+I)+SSIN36*B(IC+I))
          I = I+INC3
          J = J+INC4
        END DO
        IBASE = IBASE+INC1
        JBASE = JBASE+INC2
      END DO
      FUPDATE = .true.
    END IF

  !--
  CASE (5)
    ! Coding for factor 6

    IA = 1
    IB = IA+(2*M-LA)*INC1
    IC = IB+2*M*INC1
    ID = IC+2*M*INC1
    IE = IC
    IF = IB
    JA = 1
    JB = JA+JINK
    JC = JB+JINK
    JD = JC+JINK
    JE = JD+JINK
    JF = JE+JINK

    IF (LA /= M) THEN
      DO L=1,LA
        I = IBASE
        J = JBASE
        DO IJK=1,LOT
          C(JA+J) = (A(IA+I)+A(ID+I))+(A(IB+I)+A(IC+I))
          C(JD+J) = (A(IA+I)-A(ID+I))-(A(IB+I)-A(IC+I))
          C(JB+J) = ((A(IA+I)-A(ID+I))+0.5*(A(IB+I)-A(IC+I))) &
              -(SIN60*(B(IB+I)+B(IC+I)))
          C(JF+J) = ((A(IA+I)-A(ID+I))+0.5*(A(IB+I)-A(IC+I))) &
              +(SIN60*(B(IB+I)+B(IC+I)))
          C(JC+J) = ((A(IA+I)+A(ID+I))-0.5*(A(IB+I)+A(IC+I))) &
              -(SIN60*(B(IB+I)-B(IC+I)))
          C(JE+J) = ((A(IA+I)+A(ID+I))-0.5*(A(IB+I)+A(IC+I))) &
              +(SIN60*(B(IB+I)-B(IC+I)))
          I = I+INC3
          J = J+INC4
        END DO
        IBASE = IBASE+INC1
        JBASE = JBASE+INC2
      END DO
      IA = IA+IINK
      IINK = 2*IINK
      IB = IB+IINK
      IC = IC+IINK
      ID = ID-IINK
      IE = IE-IINK
      IF = IF-IINK
      JBASE = JBASE+JUMP
      JUMP = 2*JUMP+JINK
      IF (IC /= ID) THEN
        DO K=LA,KSTOP,LA
          KB = K+K
          KC = KB+KB
          KD = KC+KB
          KE = KD+KB
          KF = KE+KB
          C1 = TRIGS(KB+1)
          S1 = TRIGS(KB+2)
          C2 = TRIGS(KC+1)
          S2 = TRIGS(KC+2)
          C3 = TRIGS(KD+1)
          S3 = TRIGS(KD+2)
          C4 = TRIGS(KE+1)
          S4 = TRIGS(KE+2)
          C5 = TRIGS(KF+1)
          S5 = TRIGS(KF+2)
          IBASE = 0
          DO L=1,LA
            I = IBASE
            J = JBASE
            DO IJK=1,LOT

              A11 = (A(IE+I)+A(IB+I))+(A(IC+I)+A(IF+I))
              A20 = (A(IA+I)+A(ID+I))-0.5*A11
              A21 = SIN60*((A(IE+I)+A(IB+I))-(A(IC+I)+A(IF+I)))
              B11 = (B(IB+I)-B(IE+I))+(B(IC+I)-B(IF+I))
              B20 = (B(IA+I)-B(ID+I))-0.5*B11
              B21 = SIN60*((B(IB+I)-B(IE+I))-(B(IC+I)-B(IF+I)))

              C(JA+J) = (A(IA+I)+A(ID+I))+A11
              D(JA+J) = (B(IA+I)-B(ID+I))+B11
              C(JC+J) = C2*(A20-B21)-S2*(B20+A21)
              D(JC+J) = S2*(A20-B21)+C2*(B20+A21)
              C(JE+J) = C4*(A20+B21)-S4*(B20-A21)
              D(JE+J) = S4*(A20+B21)+C4*(B20-A21)

              A11 = (A(IE+I)-A(IB+I))+(A(IC+I)-A(IF+I))
              B11 = (B(IE+I)+B(IB+I))-(B(IC+I)+B(IF+I))
              A20 = (A(IA+I)-A(ID+I))-0.5*A11
              A21 = SIN60*((A(IE+I)-A(IB+I))-(A(IC+I)-A(IF+I)))
              B20 = (B(IA+I)+B(ID+I))+0.5*B11
              B21 = SIN60*((B(IE+I)+B(IB+I))+(B(IC+I)+B(IF+I)))

              C(JD+J) = C3*((A(IA+I)-A(ID+I))+A11)-S3*((B(IA+I)+B(ID+I))-B11)
              D(JD+J) = S3*((A(IA+I)-A(ID+I))+A11)+C3*((B(IA+I)+B(ID+I))-B11)
              C(JB+J) = C1*(A20-B21)-S1*(B20-A21)
              D(JB+J) = S1*(A20-B21)+C1*(B20-A21)
              C(JF+J) = C5*(A20+B21)-S5*(B20+A21)
              D(JF+J) = S5*(A20+B21)+C5*(B20+A21)

              I = I+INC3
              J = J+INC4
            END DO
            IBASE = IBASE+INC1
            JBASE = JBASE+INC2
          END DO
          IA = IA+IINK
          IB = IB+IINK
          IC = IC+IINK
          ID = ID-IINK
          IE = IE-IINK
          IF = IF-IINK
          JBASE = JBASE+JUMP
        END DO
        IF (IC > ID) THEN
          FUPDATE = .true.
        END IF
      END IF
      IF (FUPDATE .eqv. .false.) THEN
        IBASE = 0
        DO L=1,LA
          I = IBASE
          J = JBASE
          DO IJK=1,LOT
            C(JA+J) = A(IB+I)+(A(IA+I)+A(IC+I))
            C(JD+J) = B(IB+I)-(B(IA+I)+B(IC+I))
            C(JB+J) = (SIN60*(A(IA+I)-A(IC+I)))-(0.5*(B(IA+I)+B(IC+I))+B(IB+I))
            C(JF+J) = -(SIN60*(A(IA+I)-A(IC+I)))-(0.5*(B(IA+I)+B(IC+I))+B(IB+I))
            C(JC+J) = SIN60*(B(IC+I)-B(IA+I))+(0.5*(A(IA+I)+A(IC+I))-A(IB+I))
            C(JE+J) = SIN60*(B(IC+I)-B(IA+I))-(0.5*(A(IA+I)+A(IC+I))-A(IB+I))
            I = I+INC3
            J = J+INC4
          END DO
          IBASE = IBASE+INC1
          JBASE = JBASE+INC2
        END DO
        FUPDATE = .true.
      END IF
    END IF
    IF (FUPDATE .eqv. .false.) THEN
      SSIN60 = 2.0*SIN60
      DO L=1,LA
        I = IBASE
        J = JBASE
        DO IJK=1,LOT
          C(JA+J) = (2.0*(A(IA+I)+A(ID+I)))+(2.0*(A(IB+I)+A(IC+I)))
          C(JD+J) = (2.0*(A(IA+I)-A(ID+I)))-(2.0*(A(IB+I)-A(IC+I)))
          C(JB+J) = (2.0*(A(IA+I)-A(ID+I))+(A(IB+I)-A(IC+I))) &
              -(SSIN60*(B(IB+I)+B(IC+I)))
          C(JF+J) = (2.0*(A(IA+I)-A(ID+I))+(A(IB+I)-A(IC+I))) &
              +(SSIN60*(B(IB+I)+B(IC+I)))
          C(JC+J) = (2.0*(A(IA+I)+A(ID+I))-(A(IB+I)+A(IC+I))) &
              -(SSIN60*(B(IB+I)-B(IC+I)))
          C(JE+J) = (2.0*(A(IA+I)+A(ID+I))-(A(IB+I)+A(IC+I))) &
              +(SSIN60*(B(IB+I)-B(IC+I)))
          I = I+INC3
          J = J+INC4
        END DO
        IBASE = IBASE+INC1
        JBASE = JBASE+INC2
      END DO
      FUPDATE = .true.
    END IF

  !--
  CASE (6,7)
    ! Coding for factor 8

    IF (LA /= M) THEN
      ! error: IFAC only catered for if LA=N/IFAC
      IBAD = 3
      FUPDATE = .false.
    ELSE  
      IA = 1
      IB = IA+LA*INC1
      IC = IB+2*LA*INC1
      ID = IC+2*LA*INC1
      IE = ID+2*LA*INC1
      JA = 1
      JB = JA+JINK
      JC = JB+JINK
      JD = JC+JINK
      JE = JD+JINK
      JF = JE+JINK
      JG = JF+JINK
      JH = JG+JINK
      SSIN45 = SQRT(2.0)

      DO L=1,LA
        I = IBASE
        J = JBASE
        DO IJK=1,LOT
          C(JA+J) = 2.0*(((A(IA+I)+A(IE+I))+A(IC+I))+(A(IB+I)+A(ID+I)))
          C(JE+J) = 2.0*(((A(IA+I)+A(IE+I))+A(IC+I))-(A(IB+I)+A(ID+I)))
          C(JC+J) = 2.0*(((A(IA+I)+A(IE+I))-A(IC+I))-(B(IB+I)-B(ID+I)))
          C(JG+J) = 2.0*(((A(IA+I)+A(IE+I))-A(IC+I))+(B(IB+I)-B(ID+I)))
          C(JB+J) = 2.0*((A(IA+I)-A(IE+I))-B(IC+I)) &
              +SSIN45*((A(IB+I)-A(ID+I))-(B(IB+I)+B(ID+I)))
          C(JF+J) = 2.0*((A(IA+I)-A(IE+I))-B(IC+I)) &
              -SSIN45*((A(IB+I)-A(ID+I))-(B(IB+I)+B(ID+I)))
          C(JD+J) = 2.0*((A(IA+I)-A(IE+I))+B(IC+I)) &
              -SSIN45*((A(IB+I)-A(ID+I))+(B(IB+I)+B(ID+I)))
          C(JH+J) = 2.0*((A(IA+I)-A(IE+I))+B(IC+I)) &
              +SSIN45*((A(IB+I)-A(ID+I))+(B(IB+I)+B(ID+I)))
          I = I+INC3
          J = J+INC4
        END DO
        IBASE = IBASE+INC1
        JBASE = JBASE+INC2
      END DO
      FUPDATE = .true.
    END IF

  !--
  CASE DEFAULT
    FUPDATE = .false.

  END SELECT

  IF (FUPDATE .eqv. .true.) THEN
    ! final update ...
    ! pass completed without error
    IBAD = 0
  END IF

  IERR = IBAD


  end subroutine fft_rpassm_f


end module fft_pass_proc_f
