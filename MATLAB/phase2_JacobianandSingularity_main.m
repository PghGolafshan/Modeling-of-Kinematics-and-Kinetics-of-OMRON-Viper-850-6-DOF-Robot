clear
clc
%% Tranformation matrixes
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
%% Alternative
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
J_alternative=[J1 J2 J3 J4 J5 J6]

% Theta1=deg2rad(30);Theta2=deg2rad(30);
% Theta3=deg2rad(30);Theta4=deg2rad(30);
% Theta5=deg2rad(30);Theta6=deg2rad(30);
% J_alternative_numerical=eval(J_alternative)
%% Propagation
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
syms Thetad1 Thetad2 Thetad3 Thetad4 Thetad5 Thetad6
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
R0e=R01 * R12 * R23 * R34 * R45 * R56 *R6e ;
v0e=  R0e * vee; 
w0e= R0e * wee ;
v1(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, Thetad1, Thetad2, Thetad3, Thetad4, Thetad5, Thetad6)=v0e(1,1);
v2(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, Thetad1, Thetad2, Thetad3, Thetad4, Thetad5, Thetad6)=v0e(2,1);
v3(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, Thetad1, Thetad2, Thetad3, Thetad4, Thetad5, Thetad6)=v0e(3,1);
Jv(1,1)= v1(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 1, 0, 0, 0, 0, 0);
Jv(2,1)= v2(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 1, 0, 0, 0, 0, 0);
Jv(3,1)= v3(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 1, 0, 0, 0, 0, 0);
Jv(1,2)= v1(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 1, 0, 0, 0, 0);
Jv(2,2)= v2(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 1, 0, 0, 0, 0);
Jv(3,2)= v3(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 1, 0, 0, 0, 0);
Jv(1,3)= v1(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 0, 1, 0, 0, 0);
Jv(2,3)= v2(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 0, 1, 0, 0, 0);
Jv(3,3)= v3(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 0, 1, 0, 0, 0);
Jv(1,4)= v1(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 0, 0, 1, 0, 0);
Jv(2,4)= v2(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 0, 0, 1, 0, 0);
Jv(3,4)= v3(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 0, 0, 1, 0, 0);
Jv(1,5)= v1(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 0, 0, 0, 1, 0);
Jv(2,5)= v2(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 0, 0, 0, 1, 0);
Jv(3,5)= v3(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 0, 0, 0, 1, 0);
Jv(1,6)= v1(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 0, 0, 0, 0, 1);
Jv(2,6)= v2(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 0, 0, 0, 0, 1);
Jv(3,6)= v3(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 0, 0, 0, 0, 1);
w1(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, Thetad1, Thetad2, Thetad3, Thetad4, Thetad5, Thetad6)=w0e(1,1);
w2(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, Thetad1, Thetad2, Thetad3, Thetad4, Thetad5, Thetad6)=w0e(2,1);
w3(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, Thetad1, Thetad2, Thetad3, Thetad4, Thetad5, Thetad6)=w0e(3,1);
Jw(1,1)= w1(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 1, 0, 0, 0, 0, 0);
Jw(2,1)= w2(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 1, 0, 0, 0, 0, 0);
Jw(3,1)= w3(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 1, 0, 0, 0, 0, 0);
Jw(1,2)= w1(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 1, 0, 0, 0, 0);
Jw(2,2)= w2(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 1, 0, 0, 0, 0);
Jw(3,2)= w3(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 1, 0, 0, 0, 0);
Jw(1,3)= w1(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 0, 1, 0, 0, 0);
Jw(2,3)= w2(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 0, 1, 0, 0, 0);
Jw(3,3)= w3(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 0, 1, 0, 0, 0);
Jw(1,4)= w1(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 0, 0, 1, 0, 0);
Jw(2,4)= w2(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 0, 0, 1, 0, 0);
Jw(3,4)= w3(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 0, 0, 1, 0, 0);
Jw(1,5)= w1(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 0, 0, 0, 1, 0);
Jw(2,5)= w2(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 0, 0, 0, 1, 0);
Jw(3,5)= w3(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 0, 0, 0, 1, 0);
Jw(1,6)= w1(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 0, 0, 0, 0, 1);
Jw(2,6)= w2(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 0, 0, 0, 0, 1);
Jw(3,6)= w3(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 0, 0, 0, 0, 1);
J_propagation=[Jv;Jw]

% Theta1=deg2rad(30);Theta2=deg2rad(30);
% Theta3=deg2rad(30);Theta4=deg2rad(30);
% Theta5=deg2rad(30);Theta6=deg2rad(30);
% J_propagation_numerical=eval(J_propagation)

%% Alternative and Propagation Jacobian compare for same Theta variables

% Theta1=deg2rad(30);Theta2=deg2rad(30);
% Theta3=deg2rad(30);Theta4=deg2rad(30);
% Theta5=deg2rad(30);Theta6=deg2rad(30);
% differenceAlternativePropagation=eval(J_alternative)-eval(J_propagation)

%% determinant of Jacobian matrix & singularity angles
Determinant_J_alternative=simplify(det(J_alternative))
i = 1;
for Theta2=deg2rad(-190):0.05:deg2rad(45)
    for Theta3=deg2rad(-29):0.05:deg2rad(256)
        determinant_t2_t3=((885*sin(2*Theta2 + 2*Theta3))/2 + (6379*cos(Theta2 + Theta3))/2 + (3465*cos(Theta3))/2 + (4543*sin(Theta2))/2 + (5929*cos(Theta2 - Theta3))/2 + (1155*cos(2*Theta2 + Theta3))/2 + (4543*sin(Theta2 + 2*Theta3))/2)/8000000;
        vpa(determinant_t2_t3);
        if (determinant_t2_t3<=0.0001)
            singularity_Data(i,1) = rad2deg(Theta2);
            singularity_Data(i,2) = rad2deg(Theta3);
            i = i+1;
        end
    end
end
writematrix(singularity_Data,'singularity_Data_excel.xlsx','WriteMode','append');
