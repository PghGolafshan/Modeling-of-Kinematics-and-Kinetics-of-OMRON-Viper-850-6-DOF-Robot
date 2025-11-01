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

Theta1=deg2rad(30);Theta2=deg2rad(30);
Theta3=deg2rad(30);Theta4=deg2rad(30);
Theta5=deg2rad(30);Theta6=deg2rad(30);
J_alternative_numerical=eval(J_alternative)