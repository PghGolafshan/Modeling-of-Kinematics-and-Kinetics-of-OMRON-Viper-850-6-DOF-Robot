clear 
clc
%% transformation matrix
syms Theta1 Theta2 Theta3 Theta4 Theta5 Theta6
T01 = [cos(Theta1), -sin(Theta1), 0, 0; sin(Theta1),cos(Theta1), 0, 0; 0, 0, 1, 0.205;0, 0, 0, 1];
T12 = [cos(Theta2), -sin(Theta2), 0, 0.075; 0,0, 1, 0; -sin(Theta2), -cos(Theta2), 0, 0; 0, 0, 0, 1];
T23 = [cos(Theta3), -sin(Theta3), 0, 0.385; sin(Theta3),cos(Theta3), 0, 0; 0, 0, 1, 0; 0, 0, 0, 1];
T34 = [cos(Theta4), -sin(Theta4), 0, 0; 0,0, -1, -0.295; sin(Theta4), cos(Theta4), 0, 0; 0, 0, 0, 1];
T45 = [cos(Theta5), -sin(Theta5), 0, 0; 0,0, -1, 0; sin(Theta5), cos(Theta5), 0, 0; 0, 0, 0, 1];
T56 = [cos(Theta6), -sin(Theta6), 0, 0; 0,0, 1, 0; -sin(Theta6), -cos(Theta6), 0, 0; 0, 0, 0, 1];
T6e=[1,0,0,0;0,1,0,0;0,0,1,0.08;0,0,0,1];
T02 = T01 * T12 ;T03 = T01 * T12 * T23 ; 
T04 = T01 * T12 * T23 * T34 ;
T05=  T01 * T12 * T23 * T34 * T45 ;
T06=  T01 * T12 * T23 * T34 * T45 * T56;
T0e=  T01 * T12 * T23 * T34 * T45 * T56 *T6e;
%% velocity & acceleration propagation(Newton-Euler)
% P & R
P01=T01(1:3,4);P12=T12(1:3,4);
P23=T23(1:3,4);P34=T34(1:3,4);
P45=T45(1:3,4);P56=T56(1:3,4);P6e=T6e(1:3,4);
R01=T01(1:3,1:3);R12=T12(1:3,1:3);
R23=T23(1:3,1:3);R34=T34(1:3,1:3);
R45=T45(1:3,1:3);R56=T56(1:3,1:3);R6e=T6e(1:3,1:3);
R10=transpose(R01);R21=transpose(R12);
R32=transpose(R23);R43=transpose(R34);
R54=transpose(R45);R65=transpose(R56);Re6=transpose(R6e);
syms Thetad1 Thetad2 Thetad3 Thetad4 Thetad5 Thetad6 Thetadd1 Thetadd2 Thetadd3 Thetadd4 Thetadd5 Thetadd6
w00=[0;0;0];
w11=R10*w00+[0;0;Thetad1];
w22=R21*w11+[0;0;Thetad2];
w33=R32*w22+[0;0;Thetad3];
w44=R43*w33+[0;0;Thetad4];
w55=R54*w44+[0;0;Thetad5];
w66=R65*w55+[0;0;Thetad6];
wee=R6e*w66+[0;0;0];
v00=[0;0;0];
v11=R10*(v00+cross(w00,P01));
v22=R21*(v11+cross(w11,P12));
v33=R32*(v22+cross(w22,P23));
v44=R43*(v33+cross(w33,P34));
v55=R54*(v44+cross(w44,P45));
v66=R65*(v55+cross(w55,P56));
vee=Re6*(v66+cross(w66,P6e));
R06=R01 * R12 * R23 * R34 * R45 * R56;
R0e=R01 * R12 * R23 * R34 * R45 * R56 *R6e ;
v06= R06 * v66;v0e= R0e * vee; 
w06= R06 * w66;w0e= R0e * wee ;
% w dot
w11_d = simplify(R10*[0;0;0]+R10*(cross([0;0;0],[0;0;Thetad1]))+[0;0;Thetadd1]);
w22_d = simplify(R21*w11+R21*(cross(w11,[0;0;Thetad2]))+[0;0;Thetadd2]);
w33_d = simplify(R32*w22+R32*(cross(w22,[0;0;Thetad3]))+[0;0;Thetadd3]);
w44_d = simplify(R43*w33+R43*(cross(w33,[0;0;Thetad4]))+[0;0;Thetadd4]);
w55_d = simplify(R54*w44+R54*(cross(w44,[0;0;Thetad5]))+[0;0;Thetadd5]);
w66_d = simplify(R65*w55+R65*(cross(w55,[0;0;Thetad6]))+[0;0;Thetadd6]);
% v dot
v11_d=R10*(cross([0;0;0],P01)+cross([0;0;0],cross([0;0;0],P01))+[0;0;9.81]);
v22_d=R21*(cross(w11_d,P12)+cross(w11,cross(w11,P12))+v11_d);
v33_d=R32*(cross(w22_d,P23)+cross(w22,cross(w22,P23))+v22_d);
v44_d=R43*(cross(w33_d,P34)+cross(w33,cross(w33,P34))+v33_d);
v55_d=R54*(cross(w44_d,P45)+cross(w44,cross(w44,P45))+v44_d);
v66_d=R65*(cross(w55_d,P56)+cross(w55,cross(w55,P56))+v55_d);
%v dot center of mass
v1c1_d=cross(w11_d,[0.04;0;0.09])+cross(w11,cross(w11,[0.04;0;0.09]))+v11_d;
v2c2_d=cross(w22_d,[0.16;0;0.12])+cross(w22,cross(w22,[0.16;0;0.12]))+v22_d;
v3c3_d=cross(w33_d,[0.05;0.06;0])+cross(w33,cross(w33,[0.05;0.06;0]))+v33_d;
v4c4_d=cross(w44_d,[0;0;0.1])+cross(w44,cross(w44,[0;0;0.1]))+v44_d;
v5c5_d=cross(w55_d,[0;0;0.02])+cross(w55,cross(w55,[0;0;0.02]))+v55_d;
v6c6_d=cross(w66_d,[0;0;0])+cross(w66,cross(w66,[0;0;0]))+v66_d;
%% mass properties,Inertia,F,f,N,n,taw
m1 = 15.26;%[kg]
m2 = 14.93;%[kg]
m3 = 12.96 ;%[kg]
m4 = 5.14;%[kg]
m5 = 1.11;%[kg]
m6=0;
%F
F11 = v1c1_d * m1;F22 = v2c2_d * m2;
F33 = v3c3_d * m3;F44 = v4c4_d * m4;
F55 = v5c5_d * m5;F66 = v6c6_d * m6;
%inertia
I11 = [0.09,0,-0.03;0,0.1,0;-0.03,0,0.09];
I22 = [0.03,0,0;0,0.32,0;0,0,0.34];
I33 = [0.09,-0.03,0;-0.03,0.06,0;0,0,0.11];
I44 = [0.02,0,0;0,0.02,0;0,0,0.01];
I55 = [0.00137,0,0;0,0.0007,0;0,0,0.0015];
I66=[0.00137,0,0;0,0.0007,0;0,0,0.0015];
%N
N11 = I11 * w11_d + cross(w11,I11*w11);
N22 = I22 * w22_d + cross(w22,I22*w22);
N33 = I33 * w33_d + cross(w33,I33*w33);
N44 = I44 * w44_d + cross(w44,I44*w44);
N55 = I55 * w55_d + cross(w55,I55*w55);
N66 = I66 * w66_d + cross(w66,I66*w66);
%f
f66=0+F66;
f55=R56*f66+F55;
f44=R45*f55+F44;
f33=R34*f44+F33;
f22=R23*f33+F22;
f11=R12*f22+F11;
%n
n66=N66;
n55=N55+R56*n66+cross([0;0;0.02],F55)+cross(-P56,R56*f66);
n44=N44+R45*n55+cross([0;0;0.1],F44)+cross(-P45,R45*f55);
n33=N33+R34*n44+cross([0.05;0.06;0],F33)+cross(-P34,R34*f44);
n22=N22+R23*n33+cross([0.16;0;0.12],F22)+cross(-P23,R23*f33);
n11=N11+R12*n22+cross([0.04;0;0.09],F11)+cross(-P12,R12*f22);
%taw
taw6=transpose(n66)*[0;0;1];
taw5=transpose(n55)*[0;0;1];
taw4=transpose(n44)*[0;0;1];
taw3=transpose(n33)*[0;0;1];
taw2=transpose(n22)*[0;0;1];
taw1=transpose(n11)*[0;0;1];
% taw matrix
taw =[vpa(simplify(taw1),3) ;vpa(simplify(taw2),3) ;vpa(simplify(taw3),3) ;vpa(simplify(taw4),3) ;vpa(simplify(taw5),3) ;vpa(simplify(taw6),3)];
%% G
Thetadd1=0;Thetadd2=0;Thetadd3=0;
Thetadd4=0;Thetadd5=0;Thetadd6=0;
Thetad1=0;Thetad2=0;Thetad3=0;
Thetad4=0;Thetad5=0;Thetad6=0;
g(1)=eval(taw1);g(2)=eval(taw2);g(3)=eval(taw3);
g(4)=eval(taw4);g(5)=eval(taw5);g(6)=eval(taw6);
G=transpose(g);
%% M
%For First column
Thetadd1=1;Thetadd2=0;Thetadd3=0;
Thetadd4=0;Thetadd5=0;Thetadd6=0;
Thetad1=0;Thetad2=0;Thetad3=0;
Thetad4=0;Thetad5=0;Thetad6=0;

M(1,1)=eval(taw1)-G(1);
M(2,1)=eval(taw2)-G(2);
M(3,1)=eval(taw3)-G(3);
M(4,1)=eval(taw4)-G(4);
M(5,1)=eval(taw5)-G(5);
M(6,1)=eval(taw6)-G(6);
% For column 2
Thetadd1=0;Thetadd2=1;Thetadd3=0;
Thetadd4=0;Thetadd5=0;Thetadd6=0;

M(1,2)=eval(taw1)-G(1);
M(2,2)=eval(taw2)-G(2);
M(3,2)=eval(taw3)-G(3);
M(4,2)=eval(taw4)-G(4);
M(5,2)=eval(taw5)-G(5);
M(6,2)=eval(taw6)-G(6);
% For column 3
Thetadd1=0;Thetadd2=0;Thetadd3=1;
Thetadd4=0;Thetadd5=0;Thetadd6=0;

M(1,3)=eval(taw1)-G(1);
M(2,3)=eval(taw2)-G(2);
M(3,3)=eval(taw3)-G(3);
M(4,3)=eval(taw4)-G(4);
M(5,3)=eval(taw5)-G(5);
M(6,3)=eval(taw6)-G(6);

%For column 4
Thetadd1=0;Thetadd2=0;Thetadd3=0;
Thetadd4=1;Thetadd5=0;Thetadd6=0;

M(1,4)=eval(taw1)-G(1);
M(2,4)=eval(taw2)-G(2);
M(3,4)=eval(taw3)-G(3);
M(4,4)=eval(taw4)-G(4);
M(5,4)=eval(taw5)-G(5);
M(6,4)=eval(taw6)-G(6);
% For column 5
Thetadd1=0;Thetadd2=0;Thetadd3=0;
Thetadd4=0;Thetadd5=1;Thetadd6=0;

M(1,5)=eval(taw1)-G(1);
M(2,5)=eval(taw2)-G(2);
M(3,5)=eval(taw3)-G(3);
M(4,5)=eval(taw4)-G(4);
M(5,5)=eval(taw5)-G(5);
M(6,5)=eval(taw6)-G(6);

% For column 6
Thetadd1=0;Thetadd2=0;Thetadd3=0;
Thetadd4=0;Thetadd5=0;Thetadd6=1;

M(1,6)=eval(taw1)-G(1);
M(2,6)=eval(taw2)-G(2);
M(3,6)=eval(taw3)-G(3);
M(4,6)=eval(taw4)-G(4);
M(5,6)=eval(taw5)-G(5);
M(6,6)=eval(taw6)-G(6);

%% C 
%column 1
Thetadd1=0;Thetadd2=0;Thetadd3=0;
Thetadd4=0;Thetadd5=0;Thetadd6=0;
Thetad1=1;Thetad2=0;Thetad3=0;
Thetad4=0;Thetad5=0;Thetad6=0;

C(1,1)=eval(taw1)-G(1);
C(2,1)=eval(taw2)-G(2);
C(3,1)=eval(taw3)-G(3);
C(4,1)=eval(taw4)-G(4);
C(5,1)=eval(taw5)-G(5);
C(6,1)=eval(taw6)-G(6);
%column 2
Thetad1=0;Thetad2=1;Thetad3=0;
Thetad4=0;Thetad5=0;Thetad6=0;

C(1,2)=eval(taw1)-G(1);
C(2,2)=eval(taw2)-G(2);
C(3,2)=eval(taw3)-G(3);
C(4,2)=eval(taw4)-G(4);
C(5,2)=eval(taw5)-G(5);
C(6,2)=eval(taw6)-G(6);
%column 3
Thetad1=0;Thetad2=0;Thetad3=1;
Thetad4=0;Thetad5=0;Thetad6=0;

C(1,3)=eval(taw1)-G(1);
C(2,3)=eval(taw2)-G(2);
C(3,3)=eval(taw3)-G(3);
C(4,3)=eval(taw4)-G(4);
C(5,3)=eval(taw5)-G(5);
C(6,3)=eval(taw6)-G(6);
%column 4
Thetad1=0;Thetad2=0;Thetad3=0;
Thetad4=1;Thetad5=0;Thetad6=0;

C(1,4)=eval(taw1)-G(1);
C(2,4)=eval(taw2)-G(2);
C(3,4)=eval(taw3)-G(3);
C(4,4)=eval(taw4)-G(4);
C(5,4)=eval(taw5)-G(5);
C(6,4)=eval(taw6)-G(6);
%column 5
Thetad1=0;Thetad2=0;Thetad3=0;
Thetad4=0;Thetad5=1;Thetad6=0;

C(1,5)=eval(taw1)-G(1);
C(2,5)=eval(taw2)-G(2);
C(3,5)=eval(taw3)-G(3);
C(4,5)=eval(taw4)-G(4);
C(5,5)=eval(taw5)-G(5);
C(6,5)=eval(taw6)-G(6);
%column 6
Thetad1=0;Thetad2=0;Thetad3=0;
Thetad4=0;Thetad5=0;Thetad6=1;

C(1,6)=eval(taw1)-G(1);
C(2,6)=eval(taw2)-G(2);
C(3,6)=eval(taw3)-G(3);
C(4,6)=eval(taw4)-G(4);
C(5,6)=eval(taw5)-G(5);
C(6,6)=eval(taw6)-G(6);
clear Thetad1 Theatad2 Thetad3 Thetad4 Thetad5 Thetad6
syms Thetad1 Theatad2 Thetad3 Thetad4 Thetad5 Thetad6
Thetad=[Thetad1;Thetad2;Thetad3;Thetad4;Thetad5;Thetad6];
V=C*Thetad;

%% numerical
Theta1=deg2rad(30);Theta2=deg2rad(30);
Theta3=deg2rad(30);Theta4=deg2rad(30);
Theta5=deg2rad(30);Theta6=deg2rad(30);
Thetad1=0.1;Thetad2=0.1;Thetad3=0.1;
Thetad4=0.1;Thetad5=0.1;Thetad6=0.1;
Thetadd1=0.1;Thetadd2=0.1;Thetadd3=0.1;
Thetad4d=0.1;Thetadd5=0.1;Thetadd6=0.1;
M_numerical=eval(M);
C_numerical=eval(C);
G_numerical=eval(G);
Thetadd=[Thetadd1;Thetadd2;Thetadd3;Thetadd4;Thetadd5;Thetadd6];
% %for I66 and m6==0
taw_NewtonEuler_numerical=eval(M*Thetadd + V + G);
taw_answer_symbolic=M*Thetadd + V + G
%% cartesian
%z and o
z1=T01(1:3,3);z2=T02(1:3,3);
z3=T03(1:3,3);z4=T04(1:3,3);
z5=T05(1:3,3);z6=T06(1:3,3);
o1=T01(1:3,4);o2=T02(1:3,4);
o3=T03(1:3,4);o4=T04(1:3,4);
o5=T05(1:3,4);o6=T06(1:3,4);oe=T0e(1:3,4);
%Jacobians
J1=[cross(z1,(oe-o1));z1];J2=[cross(z2,(oe-o1));z2];
J3=[cross(z3,(oe-o3));z3];J4=[cross(z4,(oe-o4));z4];
J5=[cross(z5,(oe-o5));z5];J6=[cross(z6,(oe-o6));z6];
J_alternative=[J1 J2 J3 J4 J5 J6];
%Mx=M*inv(J_alternative)
%Cx=C-M*inv(J_alternative)*diff(J_alternative)*[Thetad1;Thetad2;Thetad3;Thetad4;Thetad5;Thetad6]
