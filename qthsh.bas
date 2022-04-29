' TANH-SINH QUADRATURE
' Estimate definite integral over open interval with improved and optimized Tanh-Sinh quadrature
' For SHARP PC-12xx to 14xx by Robert van Engelen
' Credit:
'   Michalski & Mosig Tanh-Sinh rule
' See also:
'   https://www.genivia.com/files/qthsh.pdf
'   https://newtonexcelbach.com/2020/10/29/numerical-integration-with-tanh-sinh-quadrature-v-5-0/

' Functions to integrate are defined with label "F1", "F2",... should return Y given X

' VARIABLES
'  A,B     range
'  F$      function label to integrate
'  Y       result with error E
'  E       estimated relative error of the result
'  N       levels (up to 6 or 7)
'  C       (a+b)/2 center
'  D       (b-a)/2 half distance
'  H       step size h=2^-k
'  K       level counter
'  P,Q,S   quadrature sums
'  T       exp(j*h)
'  L       fp=f(a+x)
'  M       fm=f(b-x)
'  J,R,U,X scratch

100 "QTHSH" E=1E-9,N=6: INPUT "f=F";F$: F$="F"+F$
110 INPUT "a=";A
120 INPUT "b=";B
' init
130 C=(A+B)/2,D=(B-A)/2,X=C: GOSUB F$: S=Y,H=1,K=0
' outer loop
140 J=1,P=0,L=0,M=0
' inner loop
150 T=EXP(J*H),U=EXP(1/T-T),R=2*U/(1+U)
160 X=A+D*R: IF X>A GOSUB F$: L=Y
170 X=B-D*R: IF X<B GOSUB F$: M=Y
180 Q=(T+1/T)*R/(1+U)*(L+M),P=P+Q,J=J+1+(K>0)
190 IF ABS Q>E*ABS P GOTO 150
' exit inner loop
200 X=S-P,S=S+P,K=K+1
210 IF ABS X>10*E*ABS S IF K<=N LET H=H/2: GOTO 140
' exit outer loop, output result (and relative error estimate if >E)
220 Y=D*S*H,U=ABS X/(ABS S+E)
230 IF U>E LET E=U: PRINT Y,E: END
240 E=U: PRINT Y: END
' For machines with ON ERROR GOTO such as PC-1475, update the following:
' 100 "QTHSH" ON ERROR GOTO 290: E=1E-9,N=6: INPUT "f=F";F$: F$="F"+F$
' 160 X=A+D*R,Y=L: IF X>A GOSUB F$: L=Y
' 170 X=B-D*R,Y=M: IF X<B GOSUB F$: M=Y
' 290 IF ERN=2 RETURN

'300 "F1" Y=4/(1+X*X),V=V+1: RETURN
'310 "F2" Y=ATN(1/X),V=V+1: RETURN
'320 "F3" Y=COS X*LN X,V=V+1: RETURN
'330 "F4" Y=X^-.8,V=V+1: RETURN
