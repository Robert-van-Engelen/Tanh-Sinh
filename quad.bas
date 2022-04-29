' COMBINED TANH-SINH, EXP-SINH, SINH-SINH QUADRATURE
' Estimate integral over open interval with improved Tanh-Sinh, Exp-Sinh and Sinh-Sinh quadratures
' For SHARP PC-12xx to 14xx by Robert van Engelen
' Credit:
'   Michalski & Mosig Tanh-Sinh rule
' See also:
'   https://www.genivia.com/files/qthsh.pdf
'   https://newtonexcelbach.com/2020/10/29/numerical-integration-with-tanh-sinh-quadrature-v-5-0/

' Functions to integrate are defined with label "F1", "F2",... should return Y given X

' VARIABLES
'  A,B     range, -9E99 and 9E99 are -inf and +inf
'  F$      function label to integrate
'  Y       result with error E
'  E       estimated relative error of the result
'  N       levels (up to 6 or 7)
'  C       center
'  D       half distance (Tanh-Sinh) or D=1 (Exp-Sinh)
'  G       sign
'  H       step size h=2^-k
'  K       level counter
'  P,Q,S   quadrature sums
'  T       exp(j*h)
'  L       f(x) point
'  M       f(x) point
'  J,R,U,X scratch

100 "QUAD" E=1E-9,N=6: INPUT "f=F";F$: F$="F"+F$
110 INPUT "a=";A
120 INPUT "b=";B
' init and swap bounds if b<a
130 D=1,G=1,H=1,K=0: IF B<A LET X=A,A=B,B=X,G=-1
140 IF ABS A<9E99 IF ABS B<9E99 GOTO 180
150 IF ABS A<9E99 LET C=A,X=A+D: GOTO 300
160 IF ABS B<9E99 LET C=B,X=B-D,D=-D,G=-G: GOTO 300
170 GOTO 400
' Tanh-Sinh
180 C=(A+B)/2,D=(B-A)/2,X=C: GOSUB F$: S=Y
' outer loop
190 J=1,P=0,L=0,M=0
' inner loop
200 T=EXP(J*H),U=EXP(1/T-T),R=2*U/(1+U)
210 X=A+D*R: IF X>A GOSUB F$: L=Y
220 X=B-D*R: IF X<B GOSUB F$: M=Y
230 Q=(T+1/T)*R/(1+U)*(L+M),P=P+Q,J=J+1+(K>0)
240 IF ABS Q>E*ABS P GOTO 200
' exit inner loop
250 X=S-P,S=S+P,K=K+1
260 IF ABS X>10*E*ABS S IF K<=N LET H=H/2: GOTO 190
' exit, output result (and relative error estimate if >E)
270 Y=G*D*S*H,U=ABS X/(ABS S+E)
280 IF U>E LET E=U: PRINT Y,E: END
290 E=U: PRINT Y: END
' Exp-Sinh
300 GOSUB F$: S=Y
' outer loop
310 J=1,P=0
' inner loop
320 T=EXP(J*H)/2,U=EXP(T-.25/T)
330 X=C+D/U: IF X=C GOTO 370
340 GOSUB F$: L=Y/U,X=C+D*U: GOSUB F$: M=Y*U
350 Q=(T+.25/T)*(L+M),P=P+Q,J=J+1+(K>0)
360 IF ABS Q>E*ABS P GOTO 320
' exit inner loop
370 X=S-P,S=S+P,K=K+1
380 IF ABS X>10*E*ABS S IF K<=N LET H=H/2: GOTO 310
390 GOTO 270
' Sinh-Sinh
400 X=0: GOSUB F$: S=Y
' outer loop
410 J=1,P=0
' inner loop
420 T=EXP(J*H)/2,U=EXP(T-.25/T)/2,R=U-.25/U
430 X=-R: GOSUB F$: L=Y,X=R: GOSUB F$: M=Y
440 Q=(T+.25/T)*(U+.25/U)*(L+M),P=P+Q,J=J+1+(K>0)
450 IF ABS Q>E*ABS P GOTO 420
' exit inner loop
460 X=S-P,S=S+P,K=K+1
470 IF ABS X>10*E*ABS S IF K<=N LET H=H/2: GOTO 410
480 GOTO 270
' For machines with ON ERROR GOTO such as PC-1475, update the following:
' 100 "QUAD" ON ERROR GOTO 490: E=1E-9,N=6: INPUT "f=F";F$: F$="F"+F$
' 210 X=A+D*R,Y=L: IF X>A GOSUB F$: L=Y
' 220 X=B-D*R,Y=M: IF X<B GOSUB F$: M=Y
' 490 IF ERN=2 RETURN

' 500 "F1" Y=4/(1+X*X): RETURN
' 510 "F2" Y=X^-.8: RETURN
' 520 "F3" Y=1/(X*X): RETURN
' 530 "F4" Y=X*EXP(-X*X): RETURN
