C=======================================================================

      SUBROUTINE SDIAG(N,A,D,X,E,IS)

C=======================================================================
C
c     N   order
C     A   matrix to be diagonalized
C     D   eigenvalues
C     X   eigenvectors
C     E   auxiliary field
C     IS = 1  eigenvalues are ordered and major component of X is positiv
C          0  eigenvalues are not ordered
C
C-----------------------------------------------------------------------
C     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(N,N),X(N,N),E(N),D(N)
      DATA TOL,EPS/1.E-32,1.E-10/
      IF (N.EQ.1) THEN
        D(1)=A(1,1)
        X(1,1)=1.
        RETURN
      ENDIF
C
      DO 10 I=1,N
      DO 10 J=1,I
   10 X(I,J)=A(I,J)
C
CCC   HOUSEHOLDER-REDUKTION
      I=N
   15 IF (I-2) 200,20,20
   20 L=I-2
      F=X(I,I-1)
      G=F
      H=0
      IF (L) 31,31,32
   32 DO 30 K=1,L
   30 H=H+X(I,K)*X(I,K)
   31 S=H+F*F
      IF (S-TOL) 33,34,34
   33 H=0
      GOTO 100
   34 IF (H) 100,100,40
   40 L=L+1
      G= SQRT(S)
C      G=DSQRT(S)
      IF (F.GE.0.) G=-G
      H=S-F*G
      HI=1.D0/H
      X(I,I-1)=F-G
      F=0.
      IF (L) 51,51,52
   52 DO 50 J=1,L
      X(J,I)=X(I,J)*HI
      S=0.
      DO 55 K=1,J
   55 S=S+X(J,K)*X(I,K)
      J1=J+1
      IF (L-J1) 57,58,58
   58 DO 59 K=J1,L
   59 S=S+X(K,J)*X(I,K)
   57 E(J)=S*HI
   50 F=F+S*X(J,I)
   51 F=F*HI*.5D0
C
      IF (L) 100,100,62
   62 DO 60 J=1,L
      S=X(I,J)
      E(J)=E(J)-F*S
      P=E(J)
      DO 65 K=1,J
   65 X(J,K)=X(J,K)-S*E(K)-X(I,K)*P
   60 CONTINUE
  100 CONTINUE
      D(I)=H
      E(I-1)=G
      I=I-1
      GOTO 15
CCC   BEREITSTELLEN DER TRANSFORMATIONMATRIX
  200 D(1)=0.
      E(N)=0.
      B=0.
      F=0.
      DO 210 I=1,N
      L=I-1
      IF (D(I).EQ.0.) GOTO 221
      IF (L) 221,221,222
  222 DO 220 J=1,L
      S=0.
      DO 225 K=1,L
  225 S=S+X(I,K)*X(K,J)
      DO 226 K=1,L
  226 X(K,J)=X(K,J)-S*X(K,I)
  220 CONTINUE
  221 D(I)=X(I,I)
      X(I,I)=1
      IF (L) 210,210,232
  232 DO 230 J=1,L
      X(I,J)=0.
  230 X(J,I)=0.
  210 CONTINUE
CCC   DIAGONALISIEREN DER DREIECKSMATRIX
      DO 300 L=1,N
      H=EPS*( ABS(D(L))+ ABS(E(L)))
      IF (H.GT.B) B=H
CCC   TEST FUER SPLITTING
      DO 310 J=L,N
      IF ( ABS(E(J)).LE.B) GOTO 320
  310 CONTINUE
CCC   TEST FUER KONVERGENZ
  320 IF (J.EQ.L) GOTO 300
  340 P=(D(L+1)-D(L))/(2*E(L))
      R= DSQRT(P*P+1.D0)
      PR=P+R
      IF (P.LT.0.) PR=P-R
      H=D(L)-E(L)/PR
      DO 350 I=L,N
  350 D(I)=D(I)-H
      F=F+H
CCC   QR-TRANSFORMATION
      P=D(J)
      C=1.D0
      S=0.
      I=J
  360 I=I-1
      IF (I.LT.L) GOTO 362
      G=C*E(I)
      H=C*P
      IF ( ABS(P)- ABS(E(I))) 363,364,364
  364 C=E(I)/P
      R= DSQRT(C*C+1.D0)
      E(I+1)=S*P*R
      S=C/R
      C=1.D0/R
      GOTO 365
  363 C=P/E(I)
      R= DSQRT(C*C+1.D0)
      E(I+1)=S*E(I)*R
      S=1.D0/R
      C=C/R
  365 P=C*D(I)-S*G
      D(I+1)=H+S*(C*G+S*D(I))
      DO 368 K=1,N
      H=X(K,I+1)
      X(K,I+1)=X(K,I)*S+H*C
  368 X(K,I)=X(K,I)*C-H*S
      GOTO 360
  362 E(L)=S*P
      D(L)=C*P
      IF ( ABS(E(L)).GT.B) GOTO 340
CCC   KONVERGENZ
  300 D(L)=D(L)+F
C
      IF (IS.EQ.0) RETURN
CCC   ORDNEN DER EIGENWERTE
      DO 400 I=1,N
      K=I
      P=D(I)
      J1=I+1
      IF (J1-N) 401,401,400
  401 DO 410 J=J1,N
      IF (D(J).GE.P) GOTO 410
      K=J
      P=D(J)
  410 CONTINUE
  420 IF (K.EQ.I) GOTO 400
      D(K)=D(I)
      D(I)=P
      DO 425 J=1,N
      P=X(J,I)
      X(J,I)=X(J,K)
  425 X(J,K)=P
  400 CONTINUE
C
C     SIGNUM
      DO  71 K=1,N
      S=0.
      DO 72 I=1,N
      H= ABS(X(I,K))
      IF (H.GT.S) THEN
         S=H
         IM=I
      ENDIF
   72 CONTINUE
      IF (X(IM,K).LT.0.) THEN
         DO 73 I=1,N
   73    X(I,K)=-X(I,K)
      ENDIf
   71 CONTINUE
C
      WRITE(*,2000)
 2000 FORMAT(' *** END SDIAG ***')
      RETURN
C-END-SDIAG
      END
