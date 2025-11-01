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
%% Lagrangian
Tc55=[1 0 0 0;0 1 0 0;0 0 1 0.02;0 0 0 1];
Tc44=[1 0 0 0;0 1 0 0;0 0 1 0.1;0 0 0 1];
Tc33=[1 0 0 0.05;0 1 0 0.06;0 0 1 0;0 0 0 1];
Tc22=[1 0 0 0.16;0 1 0 0;0 0 1 0.12;0 0 0 1];
Tc11=[1 0 0 0.04;0 1 0 0;0 0 1 0.09;0 0 0 1];
T0C1=T01*Tc11;
T02 =T01 * T12; 
T0C2=T02*Tc22;
T03=T01 * T12 * T23; 
T0C3=T03*Tc33;
T04 = T01 * T12 * T23 * T34;
T0C4=T04*Tc44;
T05= T01 * T12 * T23 * T34 * T45; 
T0C5=T05*Tc55;
T06=T01 * T12 * T23 * T34 * T45 * T56; 
T0C6=T06;
Z1=T01(1:3,3);Z2=T02(1:3,3);
Z3=T03(1:3,3);Z4=T04(1:3,3);
Z5=T05(1:3,3);Z6=T06(1:3,3);
O1=T01(1:3,4);Oc1=T0C1(1:3,4);
O2=T02(1:3,4);Oc2=T0C2(1:3,4);
O3=T03(1:3,4);Oc3=T0C3(1:3,4);
O4=T04(1:3,4);Oc4=T0C4(1:3,4);
O5=T05(1:3,4);Oc5=T0C5(1:3,4);
O6=T06(1:3,4);Oc6=O6;
% coulmn of jacobians of ee.
J1=[cross(Z1,(O6-O1));Z1];
J2=[cross(Z2,(O6-O2));Z2];
J3=[cross(Z3,(O6-O3));Z3];
J4=[cross(Z4,(O6-O4));Z4];
J5=[cross(Z5,(O6-O5));Z5];
J6=[cross(Z6,(O6-O6));Z6];
J=[J1 J2 J3 J4 J5 J6];
% Jc
zero=[0;0;0;0;0;0];Jc6=J;
Jv6=Jc6(1:3,:);Jw6=Jc6(4:6,:);
Jc51=[cross(Z1,(Oc5-O1));Z1];
Jc52=[cross(Z2,(Oc5-O2));Z2];
Jc53=[cross(Z3,(Oc5-O3));Z3];
Jc54=[cross(Z4,(Oc5-O4));Z4];
Jc55=[cross(Z5,(Oc5-O5));Z5];
Jc5=[Jc51 Jc52 Jc53 Jc54 Jc55 zero];
Jv5=Jc5(1:3,:);Jw5=Jc5(4:6,:);
Jc41=[cross(Z1,(Oc4-O1));Z1];
Jc42=[cross(Z2,(Oc4-O2));Z2];
Jc43=[cross(Z3,(Oc4-O3));Z3];
Jc44=[cross(Z4,(Oc4-O4));Z4];
Jc4=[Jc41 Jc42 Jc43 Jc44 zero zero];
Jv4=Jc4(1:3,:);Jw4=Jc4(4:6,:);
Jc31=[cross(Z1,(Oc3-O1));Z1];
Jc32=[cross(Z2,(Oc3-O2));Z2];
Jc33=[cross(Z3,(Oc3-O3));Z3];
Jc3=[Jc31 Jc32 Jc33 zero zero zero];
Jv3=Jc3(1:3,:);Jw3=Jc3(4:6,:);
Jc21=[cross(Z1,(Oc2-O1));Z1];
Jc22=[cross(Z2,(Oc2-O2));Z2];
Jc2=[Jc21 Jc22 zero zero zero zero];
Jv2=Jc2(1:3,:);Jw2=Jc2(4:6,:);
Jc11=[cross(Z1,(Oc1-O1));Z1];
Jc1=[Jc11 zero zero zero zero zero];
Jv1=Jc1(1:3,:);Jw1=Jc1(4:6,:);
% M
syms m1 m2 m3 m4 m5 m6 I1 I2 I3 I4 I5 I6
R01=T01(1:3,1:3);
R02=T02(1:3,1:3);
R03=T03(1:3,1:3);
R04=T04(1:3,1:3);
R05=T05(1:3,1:3);
R06=T06(1:3,1:3);
M_theta1= m1*transpose(Jv1)*Jv1+transpose(Jw1)*R01*I1*transpose(R01)*Jw1;
M_theta2= m2*transpose(Jv2)*Jv2+transpose(Jw2)*R02*I2*transpose(R02)*Jw2;
M_theta3= m3*transpose(Jv3)*Jv3+transpose(Jw3)*R03*I3*transpose(R03)*Jw3;
M_theta4= m4*transpose(Jv4)*Jv4+transpose(Jw4)*R04*I4*transpose(R04)*Jw4;
M_theta5= m5*transpose(Jv5)*Jv5+transpose(Jw5)*R05*I5*transpose(R05)*Jw5;
M_theta6= m6*transpose(Jv6)*Jv6+transpose(Jw6)*R06*I6*transpose(R06)*Jw6;
M_theta=M_theta1+M_theta2+M_theta3+M_theta4+M_theta5+M_theta6;

% C
% k=1 n1
C_111=0.5*(((diff(M_theta(1,1),Theta1)))+(diff(M_theta(1,1),Theta1))-(diff(M_theta(1,1),Theta1)));
C_121=0.5*(((diff(M_theta(1,2),Theta1)))+(diff(M_theta(1,1),Theta2))-(diff(M_theta(1,2),Theta1)));
C_211=C_121;
C_131=0.5*(((diff(M_theta(1,3),Theta1)))+(diff(M_theta(1,1),Theta3))-(diff(M_theta(1,3),Theta1)));
C_311=C_131;
C_141=0.5*(((diff(M_theta(1,4),Theta1)))+(diff(M_theta(1,1),Theta4))-(diff(M_theta(1,4),Theta1)));
C_411=C_141;
C_151=0.5*(((diff(M_theta(1,5),Theta1)))+(diff(M_theta(1,1),Theta5))-(diff(M_theta(1,5),Theta1)));
C_511=C_151;
C_161=0.5*(((diff(M_theta(1,6),Theta1)))+(diff(M_theta(1,1),Theta6))-(diff(M_theta(1,6),Theta1)));
C_611=C_161;

% k=1 n2
C_221=0.5*(((diff(M_theta(1,2),Theta2)))+(diff(M_theta(1,2),Theta2))-(diff(M_theta(2,2),Theta1)));
C_231=0.5*(((diff(M_theta(1,3),Theta2)))+(diff(M_theta(1,2),Theta3))-(diff(M_theta(2,3),Theta1)));
C_321=C_231;
C_241=0.5*(((diff(M_theta(1,4),Theta2)))+(diff(M_theta(1,2),Theta4))-(diff(M_theta(2,4),Theta1)));
C_421=C_241;
C_251=0.5*(((diff(M_theta(1,5),Theta2)))+(diff(M_theta(1,2),Theta5))-(diff(M_theta(2,5),Theta1)));
C_521=C_251;
C_261=0.5*(((diff(M_theta(1,6),Theta2)))+(diff(M_theta(1,2),Theta6))-(diff(M_theta(2,6),Theta1)));
C_621=C_261;

% k=1 n3
C_331=0.5*(((diff(M_theta(1,3),Theta3)))+(diff(M_theta(1,3),Theta3))-(diff(M_theta(3,3),Theta1)));
C_341=0.5*(((diff(M_theta(1,4),Theta3)))+(diff(M_theta(1,3),Theta4))-(diff(M_theta(3,4),Theta1)));
C_431=C_341;
C_351=0.5*(((diff(M_theta(1,5),Theta3)))+(diff(M_theta(1,3),Theta5))-(diff(M_theta(3,5),Theta1)));
C_531=C_351;
C_361=0.5*(((diff(M_theta(1,6),Theta3)))+(diff(M_theta(1,3),Theta6))-(diff(M_theta(3,6),Theta1)));
C_631=C_361;

% k=1 n4
C_441=0.5*(((diff(M_theta(1,4),Theta4)))+(diff(M_theta(1,4),Theta4))-(diff(M_theta(4,4),Theta1)));
C_451=0.5*(((diff(M_theta(1,5),Theta4)))+(diff(M_theta(1,4),Theta5))-(diff(M_theta(4,5),Theta1)));
C_541=C_451;
C_461=0.5*(((diff(M_theta(1,6),Theta4)))+(diff(M_theta(1,4),Theta6))-(diff(M_theta(4,6),Theta1)));
C_641=C_461;

% k=1 n5
C_551=0.5*(((diff(M_theta(1,5),Theta5)))+(diff(M_theta(1,5),Theta5))-(diff(M_theta(5,5),Theta1)));
C_561=0.5*(((diff(M_theta(1,6),Theta5)))+(diff(M_theta(1,5),Theta6))-(diff(M_theta(5,6),Theta1)));
C_651=C_561;

% k=1 n6
C_661=0.5*(((diff(M_theta(1,6),Theta6)))+(diff(M_theta(1,6),Theta6))-(diff(M_theta(6,6),Theta1)));

% k=2 n1
C_112=0.5*(((diff(M_theta(2,1),Theta1)))+(diff(M_theta(2,1),Theta1))-(diff(M_theta(1,1),Theta2)));
C_122=0.5*(((diff(M_theta(2,2),Theta1)))+(diff(M_theta(2,1),Theta2))-(diff(M_theta(1,2),Theta2)));
C_212=C_122;
C_132=0.5*(((diff(M_theta(2,3),Theta1)))+(diff(M_theta(2,1),Theta3))-(diff(M_theta(1,3),Theta2)));
C_312=C_132;
C_142=0.5*(((diff(M_theta(2,4),Theta1)))+(diff(M_theta(2,1),Theta4))-(diff(M_theta(1,4),Theta2)));
C_412=C_142;
C_152=0.5*(((diff(M_theta(2,5),Theta1)))+(diff(M_theta(2,1),Theta5))-(diff(M_theta(1,5),Theta2)));
C_512=C_152;
C_162=0.5*(((diff(M_theta(2,6),Theta1)))+(diff(M_theta(2,1),Theta6))-(diff(M_theta(1,6),Theta2)));
C_612=C_162;

% k=2 n2
C_222=0.5*(((diff(M_theta(2,2),Theta2)))+(diff(M_theta(2,2),Theta2))-(diff(M_theta(2,2),Theta2)));
C_232=0.5*(((diff(M_theta(2,3),Theta2)))+(diff(M_theta(2,2),Theta3))-(diff(M_theta(2,3),Theta2)));
C_322=C_232;
C_242=0.5*(((diff(M_theta(2,4),Theta2)))+(diff(M_theta(2,2),Theta4))-(diff(M_theta(2,4),Theta2)));
C_422=C_242;
C_252=0.5*(((diff(M_theta(2,5),Theta2)))+(diff(M_theta(2,2),Theta5))-(diff(M_theta(2,5),Theta2)));
C_522=C_252;
C_262=0.5*(((diff(M_theta(2,6),Theta2)))+(diff(M_theta(2,2),Theta6))-(diff(M_theta(2,6),Theta2)));
C_622=C_262;

% k=2 n3
C_332=0.5*(((diff(M_theta(2,3),Theta3)))+(diff(M_theta(2,3),Theta3))-(diff(M_theta(3,3),Theta2)));
C_342=0.5*(((diff(M_theta(2,4),Theta3)))+(diff(M_theta(2,3),Theta4))-(diff(M_theta(3,4),Theta2)));
C_432=C_342;
C_352=0.5*(((diff(M_theta(2,5),Theta3)))+(diff(M_theta(2,3),Theta5))-(diff(M_theta(3,5),Theta2)));
C_532=C_352;
C_362=0.5*(((diff(M_theta(2,6),Theta3)))+(diff(M_theta(2,3),Theta6))-(diff(M_theta(3,6),Theta2)));
C_632=C_362;

% k=2 n4
C_442=0.5*(((diff(M_theta(2,4),Theta4)))+(diff(M_theta(2,4),Theta4))-(diff(M_theta(4,4),Theta2)));
C_452=0.5*(((diff(M_theta(2,5),Theta4)))+(diff(M_theta(2,4),Theta5))-(diff(M_theta(4,5),Theta2)));
C_542=C_452;
C_462=0.5*(((diff(M_theta(2,6),Theta4)))+(diff(M_theta(2,4),Theta6))-(diff(M_theta(4,6),Theta2)));
C_642=C_462;

% k=2 n5
C_552=0.5*(((diff(M_theta(2,5),Theta5)))+(diff(M_theta(2,5),Theta5))-(diff(M_theta(5,5),Theta2)));
C_562=0.5*(((diff(M_theta(2,6),Theta5)))+(diff(M_theta(2,5),Theta6))-(diff(M_theta(5,6),Theta2)));
C_652=C_562;

% k=2 n6
C_662=0.5*(((diff(M_theta(2,6),Theta6)))+(diff(M_theta(2,6),Theta6))-(diff(M_theta(6,6),Theta2)));

% k=3 n1
C_113=0.5*(((diff(M_theta(3,1),Theta1)))+(diff(M_theta(3,1),Theta1))-(diff(M_theta(1,1),Theta3)));
C_123=0.5*(((diff(M_theta(3,2),Theta1)))+(diff(M_theta(3,1),Theta2))-(diff(M_theta(1,2),Theta3)));
C_213=C_123;
C_133=0.5*(((diff(M_theta(3,3),Theta1)))+(diff(M_theta(3,1),Theta3))-(diff(M_theta(1,3),Theta3)));
C_313=C_133;
C_143=0.5*(((diff(M_theta(3,4),Theta1)))+(diff(M_theta(3,1),Theta4))-(diff(M_theta(1,4),Theta3)));
C_413=C_143;
C_153=0.5*(((diff(M_theta(3,5),Theta1)))+(diff(M_theta(3,1),Theta5))-(diff(M_theta(1,5),Theta3)));
C_513=C_153;
C_163=0.5*(((diff(M_theta(3,6),Theta1)))+(diff(M_theta(3,1),Theta6))-(diff(M_theta(1,6),Theta3)));
C_613=C_163;

% k=3 n2
C_223=0.5*(((diff(M_theta(3,2),Theta2)))+(diff(M_theta(3,2),Theta2))-(diff(M_theta(2,2),Theta3)));
C_233=0.5*(((diff(M_theta(3,3),Theta2)))+(diff(M_theta(3,2),Theta3))-(diff(M_theta(2,3),Theta3)));
C_323=C_233;
C_243=0.5*(((diff(M_theta(3,4),Theta2)))+(diff(M_theta(3,2),Theta4))-(diff(M_theta(2,4),Theta3)));
C_423=C_243;
C_253=0.5*(((diff(M_theta(3,5),Theta2)))+(diff(M_theta(3,2),Theta5))-(diff(M_theta(2,5),Theta3)));
C_523=C_253;
C_263=0.5*(((diff(M_theta(3,6),Theta2)))+(diff(M_theta(3,2),Theta6))-(diff(M_theta(2,6),Theta3)));
C_623=C_263;

% k=3 n3
C_333=0.5*(((diff(M_theta(3,3),Theta3)))+(diff(M_theta(3,3),Theta3))-(diff(M_theta(3,3),Theta3)));
C_343=0.5*(((diff(M_theta(3,4),Theta3)))+(diff(M_theta(3,3),Theta4))-(diff(M_theta(3,4),Theta3)));
C_433=C_343;
C_353=0.5*(((diff(M_theta(3,5),Theta3)))+(diff(M_theta(3,3),Theta5))-(diff(M_theta(3,5),Theta3)));
C_533=C_353;
C_363=0.5*(((diff(M_theta(3,6),Theta3)))+(diff(M_theta(3,3),Theta6))-(diff(M_theta(3,6),Theta3)));
C_633=C_363;

% k=3 n4
C_443=0.5*(((diff(M_theta(3,4),Theta4)))+(diff(M_theta(3,4),Theta4))-(diff(M_theta(4,4),Theta3)));
C_453=0.5*(((diff(M_theta(3,5),Theta4)))+(diff(M_theta(3,4),Theta5))-(diff(M_theta(4,5),Theta3)));
C_543=C_453;
C_463=0.5*(((diff(M_theta(3,6),Theta4)))+(diff(M_theta(3,4),Theta6))-(diff(M_theta(4,6),Theta3)));
C_643=C_463;

% k=3 n5
C_553=0.5*(((diff(M_theta(3,5),Theta5)))+(diff(M_theta(3,5),Theta5))-(diff(M_theta(5,5),Theta3)));
C_563=0.5*(((diff(M_theta(3,6),Theta5)))+(diff(M_theta(3,5),Theta6))-(diff(M_theta(5,6),Theta3)));
C_653=C_563;

% k=3 n6
C_663=0.5*(((diff(M_theta(3,6),Theta6)))+(diff(M_theta(3,6),Theta6))-(diff(M_theta(6,6),Theta3)));

% k=4 n1
C_114=0.5*(((diff(M_theta(4,1),Theta1)))+(diff(M_theta(4,1),Theta1))-(diff(M_theta(1,1),Theta4)));
C_124=0.5*(((diff(M_theta(4,2),Theta1)))+(diff(M_theta(4,1),Theta2))-(diff(M_theta(1,2),Theta4)));
C_214=C_124;
C_134=0.5*(((diff(M_theta(4,3),Theta1)))+(diff(M_theta(4,1),Theta3))-(diff(M_theta(1,3),Theta4)));
C_314=C_134;
C_144=0.5*(((diff(M_theta(4,4),Theta1)))+(diff(M_theta(4,1),Theta4))-(diff(M_theta(1,4),Theta4)));
C_414=C_144;
C_154=0.5*(((diff(M_theta(4,5),Theta1)))+(diff(M_theta(4,1),Theta5))-(diff(M_theta(1,5),Theta4)));
C_514=C_154;
C_164=0.5*(((diff(M_theta(4,6),Theta1)))+(diff(M_theta(4,1),Theta6))-(diff(M_theta(1,6),Theta4)));
C_614=C_164;

% k=4 n2
C_224=0.5*(((diff(M_theta(4,2),Theta2)))+(diff(M_theta(4,2),Theta2))-(diff(M_theta(2,2),Theta4)));
C_234=0.5*(((diff(M_theta(4,3),Theta2)))+(diff(M_theta(4,2),Theta3))-(diff(M_theta(2,3),Theta4)));
C_324=C_234;
C_244=0.5*(((diff(M_theta(4,4),Theta2)))+(diff(M_theta(4,2),Theta4))-(diff(M_theta(2,4),Theta4)));
C_424=C_244;
C_254=0.5*(((diff(M_theta(4,5),Theta2)))+(diff(M_theta(4,2),Theta5))-(diff(M_theta(2,5),Theta4)));
C_524=C_254;
C_264=0.5*(((diff(M_theta(4,6),Theta2)))+(diff(M_theta(4,2),Theta6))-(diff(M_theta(2,6),Theta4)));
C_624=C_264;

% k=4 n3
C_334=0.5*(((diff(M_theta(4,3),Theta3)))+(diff(M_theta(4,3),Theta3))-(diff(M_theta(3,3),Theta4)));
C_344=0.5*(((diff(M_theta(4,4),Theta3)))+(diff(M_theta(4,3),Theta4))-(diff(M_theta(3,4),Theta4)));
C_434=C_344;
C_354=0.5*(((diff(M_theta(4,5),Theta3)))+(diff(M_theta(4,3),Theta5))-(diff(M_theta(3,5),Theta4)));
C_534=C_354;
C_364=0.5*(((diff(M_theta(4,6),Theta3)))+(diff(M_theta(4,3),Theta6))-(diff(M_theta(3,6),Theta4)));
C_634=C_364;

% k=4 n4
C_444=0.5*(((diff(M_theta(4,4),Theta4)))+(diff(M_theta(4,4),Theta4))-(diff(M_theta(4,4),Theta4)));
C_454=0.5*(((diff(M_theta(4,5),Theta4)))+(diff(M_theta(4,4),Theta5))-(diff(M_theta(4,5),Theta4)));
C_544=C_454;
C_464=0.5*(((diff(M_theta(4,6),Theta4)))+(diff(M_theta(4,4),Theta6))-(diff(M_theta(4,6),Theta4)));
C_644=C_464;

% k=4 n5
C_554=0.5*(((diff(M_theta(4,5),Theta5)))+(diff(M_theta(4,5),Theta5))-(diff(M_theta(5,5),Theta4)));
C_564=0.5*(((diff(M_theta(4,6),Theta5)))+(diff(M_theta(4,5),Theta6))-(diff(M_theta(5,6),Theta4)));
C_654=C_564;

% k=4 n6
C_664=0.5*(((diff(M_theta(4,6),Theta6)))+(diff(M_theta(4,6),Theta6))-(diff(M_theta(6,6),Theta4)));

% k=5 n1
C_115=0.5*(((diff(M_theta(5,1),Theta1)))+(diff(M_theta(5,1),Theta1))-(diff(M_theta(1,1),Theta5)));
C_125=0.5*(((diff(M_theta(5,2),Theta1)))+(diff(M_theta(5,1),Theta2))-(diff(M_theta(1,2),Theta5)));
C_215=C_125;
C_135=0.5*(((diff(M_theta(5,3),Theta1)))+(diff(M_theta(5,1),Theta3))-(diff(M_theta(1,3),Theta5)));
C_315=C_135;
C_145=0.5*(((diff(M_theta(5,4),Theta1)))+(diff(M_theta(5,1),Theta4))-(diff(M_theta(1,4),Theta5)));
C_415=C_145;
C_155=0.5*(((diff(M_theta(5,5),Theta1)))+(diff(M_theta(5,1),Theta5))-(diff(M_theta(1,5),Theta5)));
C_515=C_155;
C_165=0.5*(((diff(M_theta(5,6),Theta1)))+(diff(M_theta(5,1),Theta6))-(diff(M_theta(1,6),Theta5)));
C_615=C_165;

% k=5 n2
C_225=0.5*(((diff(M_theta(5,2),Theta2)))+(diff(M_theta(5,2),Theta2))-(diff(M_theta(2,2),Theta5)));
C_235=0.5*(((diff(M_theta(5,3),Theta2)))+(diff(M_theta(5,2),Theta3))-(diff(M_theta(2,3),Theta5)));
C_325=C_235;
C_245=0.5*(((diff(M_theta(5,4),Theta2)))+(diff(M_theta(5,2),Theta4))-(diff(M_theta(2,4),Theta5)));
C_425=C_245;
C_255=0.5*(((diff(M_theta(5,5),Theta2)))+(diff(M_theta(5,2),Theta5))-(diff(M_theta(2,5),Theta5)));
C_525=C_255;
C_265=0.5*(((diff(M_theta(5,6),Theta2)))+(diff(M_theta(5,2),Theta6))-(diff(M_theta(2,6),Theta5)));
C_625=C_265;

% k=5 n3

C_335=0.5*(((diff(M_theta(5,3),Theta3)))+(diff(M_theta(5,3),Theta3))-(diff(M_theta(3,3),Theta5)));
C_345=0.5*(((diff(M_theta(5,4),Theta3)))+(diff(M_theta(5,3),Theta4))-(diff(M_theta(3,4),Theta5)));
C_435=C_345;
C_355=0.5*(((diff(M_theta(5,5),Theta3)))+(diff(M_theta(5,3),Theta5))-(diff(M_theta(3,5),Theta5)));
C_535=C_355;
C_365=0.5*(((diff(M_theta(5,6),Theta3)))+(diff(M_theta(5,3),Theta6))-(diff(M_theta(3,6),Theta5)));
C_635=C_365;

% k=5 n4
C_445=0.5*(((diff(M_theta(5,4),Theta4)))+(diff(M_theta(5,4),Theta4))-(diff(M_theta(4,4),Theta5)));
C_455=0.5*(((diff(M_theta(5,5),Theta4)))+(diff(M_theta(5,4),Theta5))-(diff(M_theta(4,5),Theta5)));
C_545=C_455;
C_465=0.5*(((diff(M_theta(5,6),Theta4)))+(diff(M_theta(5,4),Theta6))-(diff(M_theta(4,6),Theta5)));
C_645=C_465;

% k=5 n5
C_555=0.5*(((diff(M_theta(5,5),Theta5)))+(diff(M_theta(5,5),Theta5))-(diff(M_theta(5,5),Theta5)));
C_565=0.5*(((diff(M_theta(5,6),Theta5)))+(diff(M_theta(5,5),Theta6))-(diff(M_theta(5,6),Theta5)));
C_655=C_565;

% k=5 n6
C_665=0.5*(((diff(M_theta(5,6),Theta6)))+(diff(M_theta(5,6),Theta6))-(diff(M_theta(6,6),Theta5)));

% k=6 n1
C_116=0.5*(((diff(M_theta(6,1),Theta1)))+(diff(M_theta(6,1),Theta1))-(diff(M_theta(1,1),Theta6)));
C_126=0.5*(((diff(M_theta(6,2),Theta1)))+(diff(M_theta(6,1),Theta2))-(diff(M_theta(1,2),Theta6)));
C_216=C_126;
C_136=0.5*(((diff(M_theta(6,3),Theta1)))+(diff(M_theta(6,1),Theta3))-(diff(M_theta(1,3),Theta6)));
C_316=C_136;
C_146=0.5*(((diff(M_theta(6,4),Theta1)))+(diff(M_theta(6,1),Theta4))-(diff(M_theta(1,4),Theta6)));
C_416=C_146;
C_156=0.5*(((diff(M_theta(6,5),Theta1)))+(diff(M_theta(6,1),Theta5))-(diff(M_theta(1,5),Theta6)));
C_516=C_156;
C_166=0.5*(((diff(M_theta(6,6),Theta1)))+(diff(M_theta(6,1),Theta6))-(diff(M_theta(1,6),Theta6)));
C_616=C_166;

% k=6 n2
C_226=0.5*(((diff(M_theta(6,2),Theta2)))+(diff(M_theta(6,2),Theta2))-(diff(M_theta(2,2),Theta6)));
C_236=0.5*(((diff(M_theta(6,3),Theta2)))+(diff(M_theta(6,2),Theta3))-(diff(M_theta(2,3),Theta6)));
C_326=C_236;
C_246=0.5*(((diff(M_theta(6,4),Theta2)))+(diff(M_theta(6,2),Theta4))-(diff(M_theta(2,4),Theta6)));
C_426=C_246;
C_256=0.5*(((diff(M_theta(6,5),Theta2)))+(diff(M_theta(6,2),Theta5))-(diff(M_theta(2,5),Theta6)));
C_526=C_256;
C_266=0.5*(((diff(M_theta(6,6),Theta2)))+(diff(M_theta(6,2),Theta6))-(diff(M_theta(2,6),Theta6)));
C_626=C_266;

% k=6 n3
C_336=0.5*(((diff(M_theta(6,3),Theta3)))+(diff(M_theta(6,3),Theta3))-(diff(M_theta(3,3),Theta6)));
C_346=0.5*(((diff(M_theta(6,4),Theta3)))+(diff(M_theta(6,3),Theta4))-(diff(M_theta(3,4),Theta6)));
C_436=C_346;
C_356=0.5*(((diff(M_theta(6,5),Theta3)))+(diff(M_theta(6,3),Theta5))-(diff(M_theta(3,5),Theta6)));
C_536=C_356;
C_366=0.5*(((diff(M_theta(6,6),Theta3)))+(diff(M_theta(6,3),Theta6))-(diff(M_theta(3,6),Theta6)));
C_636=C_366;

% k=6 n4
C_446=0.5*(((diff(M_theta(6,4),Theta4)))+(diff(M_theta(6,4),Theta4))-(diff(M_theta(4,4),Theta6)));
C_456=0.5*(((diff(M_theta(6,5),Theta4)))+(diff(M_theta(6,4),Theta5))-(diff(M_theta(4,5),Theta6)));
C_546=C_456;
C_466=0.5*(((diff(M_theta(6,6),Theta4)))+(diff(M_theta(6,4),Theta6))-(diff(M_theta(4,6),Theta6)));
C_646=C_466;

% k=6 n5
C_556=0.5*(((diff(M_theta(6,5),Theta5)))+(diff(M_theta(6,5),Theta5))-(diff(M_theta(5,5),Theta6)));
C_566=0.5*(((diff(M_theta(6,6),Theta5)))+(diff(M_theta(6,5),Theta6))-(diff(M_theta(5,6),Theta6)));
C_656=C_566;

% k=6 n6
C_666=0.5*(((diff(M_theta(6,6),Theta6)))+(diff(M_theta(6,6),Theta6))-(diff(M_theta(6,6),Theta6)));

% C final 

C(1)=C_111+C_121+C_131+C_141+C_151+C_161+C_211+C_221+C_231+C_241+C_251+C_261+C_311+C_321+C_331+C_341+C_351+C_361+C_411+C_421+C_431+C_441+C_451+C_461+C_511+C_521+C_531+C_541+C_551+C_561+C_611+C_621+C_631+C_641+C_651+C_661;
C(2)=C_112+C_122+C_132+C_142+C_152+C_162+C_212+C_222+C_232+C_242+C_252+C_262+C_312+C_322+C_332+C_342+C_352+C_362+C_412+C_422+C_432+C_442+C_452+C_462+C_512+C_522+C_532+C_542+C_552+C_562+C_612+C_622+C_632+C_642+C_652+C_662;
C(3)=C_113+C_123+C_133+C_143+C_153+C_163+C_213+C_223+C_233+C_243+C_253+C_263+C_313+C_323+C_333+C_343+C_353+C_363+C_413+C_423+C_433+C_443+C_453+C_463+C_513+C_523+C_533+C_543+C_553+C_563+C_613+C_623+C_633+C_643+C_653+C_663;
c(4)=C_114+C_124+C_134+C_144+C_154+C_164+C_214+C_224+C_234+C_244+C_254+C_264+C_314+C_324+C_334+C_344+C_354+C_364+C_414+C_424+C_434+C_444+C_454+C_464+C_514+C_524+C_534+C_544+C_554+C_564+C_614+C_624+C_634+C_644+C_654+C_664;
c(5)=C_115+C_125+C_135+C_145+C_155+C_165+C_215+C_225+C_235+C_245+C_255+C_265+C_315+C_325+C_335+C_345+C_355+C_365+C_415+C_425+C_435+C_445+C_455+C_465+C_515+C_525+C_535+C_545+C_555+C_565+C_615+C_625+C_635+C_645+C_655+C_665;
c(6)=C_116+C_126+C_136+C_146+C_156+C_166+C_216+C_226+C_236+C_246+C_256+C_266+C_316+C_326+C_336+C_346+C_356+C_366+C_416+C_426+C_436+C_446+C_456+C_466+C_516+C_526+C_536+C_546+C_556+C_566+C_616+C_626+C_636+C_646+C_656+C_666;
c=transpose(c);
% G
syms m6 gr
rc1=T0C1(1:3,4);
rc2=T0C2(1:3,4);
rc3=T0C3(1:3,4);
rc4=T0C4(1:3,4);
rc5=T0C5(1:3,4);
rc6=T0C6(1:3,4);
G=[0;0;-9.81];
m1 = 15.26;%[kg]
m2 = 14.93;%[kg]
m3 = 12.96 ;%[kg]
m4 = 5.14;%[kg]
m5 = 1.11;%[kg]
I1 = [0.09,0,-0.03;0,0.1,0;-0.03,0,0.09];
I2 = [0.03,0,0;0,0.32,0;0,0,0.34];
I3 = [0.09,-0.03,0;-0.03,0.06,0;0,0,0.11];
I4 = [0.02,0,0;0,0.02,0;0,0,0.01];
I5 = [0.00137,0,0;0,0.0007,0;0,0,0.0015];
transpose(G);
P1=-m1*transpose(G)*rc1;
P2=-m2*transpose(G)*rc2;
P3=-m3*transpose(G)*rc3;
P4=-m4*transpose(G)*rc4;
P5=-m5*transpose(G)*rc5;
P6=-m6*transpose(G)*rc6;
P = P1 + P2 + P3 + P4 + P5 + P6 ;
% k=1
g1=diff(P,Theta1);
% k=2
g2=diff(P,Theta2);
% k=3
g3=diff(P,Theta3);
% k=4 
g4=diff(P,Theta4);
% k=5 
g5=diff(P,Theta5);
% k=6
g6=diff(P,Theta6);
g=[g1;g2;g3;g4;g5;g6];
% taw 
syms Thetadd1 Thetadd2 Thetadd3 Thetadd4 Thetadd5 Thetadd6
syms Thetad1 Thetad2 Thetad3 Thetad4 Thetad5 Thetad6
taw1=(M_theta(1,1)*Thetadd1+M_theta(1,2)*Thetadd2+M_theta(1,3)*Thetadd3+M_theta(1,4)*Thetadd4+M_theta(1,5)*Thetadd5+M_theta(1,6)*Thetadd6)+(C_111*Thetad1*Thetad1+C_121*Thetad1*Thetad2+C_131*Thetad1*Thetad3+C_141*Thetad1*Thetad4+C_151*Thetad1*Thetad5+C_161*Thetad1*Thetad6+C_211*Thetad2*Thetad1+C_221*Thetad2*Thetad2+C_231*Thetad2*Thetad3+C_241*Thetad2*Thetad4+C_251*Thetad2*Thetad5+C_261*Thetad2*Thetad6+C_311*Thetad3*Thetad1+C_321*Thetad3*Thetad2+C_331*Thetad3*Thetad3+C_341*Thetad3*Thetad4+C_351*Thetad3*Thetad5+C_361*Thetad3*Thetad6+C_411*Thetad4*Thetad1+C_421*Thetad4*Thetad2+C_431*Thetad4*Thetad3+C_441*Thetad4*Thetad4+C_451*Thetad4*Thetad5+C_461*Thetad4*Thetad6+C_511*Thetad5*Thetad1+C_521*Thetad5*Thetad2+C_531*Thetad5*Thetad3+C_541*Thetad5*Thetad4+C_551*Thetad5*Thetad5+C_561*Thetad5*Thetad6+C_611*Thetad6*Thetad1+C_621*Thetad6*Thetad2+C_631*Thetad6*Thetad3+C_641*Thetad6*Thetad4+C_651*Thetad6*Thetad5+C_661*Thetad1*Thetad6)+g1;
taw2=(M_theta(2,1)*Thetadd1+M_theta(2,2)*Thetadd2+M_theta(2,3)*Thetadd3+M_theta(2,4)*Thetadd4+M_theta(2,5)*Thetadd5+M_theta(2,6)*Thetadd6)+(C_112*Thetad1*Thetad1+C_122*Thetad1*Thetad2+C_132*Thetad1*Thetad3+C_142*Thetad1*Thetad4+C_152*Thetad1*Thetad5+C_162*Thetad1*Thetad6+C_212*Thetad2*Thetad1+C_222*Thetad2*Thetad2+C_232*Thetad2*Thetad3+C_242*Thetad2*Thetad4+C_252*Thetad2*Thetad5+C_262*Thetad2*Thetad6+C_312*Thetad3*Thetad1+C_322*Thetad3*Thetad2+C_332*Thetad3*Thetad3+C_342*Thetad3*Thetad4+C_352*Thetad3*Thetad5+C_362*Thetad3*Thetad6+C_412*Thetad4*Thetad1+C_422*Thetad4*Thetad2+C_432*Thetad4*Thetad3+C_442*Thetad4*Thetad4+C_452*Thetad4*Thetad5+C_462*Thetad4*Thetad6+C_512*Thetad5*Thetad1+C_522*Thetad5*Thetad2+C_532*Thetad5*Thetad3+C_542*Thetad5*Thetad4+C_552*Thetad5*Thetad5+C_562*Thetad5*Thetad6+C_612*Thetad6*Thetad1+C_622*Thetad6*Thetad2+C_632*Thetad6*Thetad3+C_642*Thetad6*Thetad4+C_652*Thetad6*Thetad5+C_662*Thetad1*Thetad6)+g2;
taw3=(M_theta(3,1)*Thetadd1+M_theta(3,2)*Thetadd2+M_theta(3,3)*Thetadd3+M_theta(3,4)*Thetadd4+M_theta(3,5)*Thetadd5+M_theta(3,6)*Thetadd6)+(C_113*Thetad1*Thetad1+C_123*Thetad1*Thetad2+C_133*Thetad1*Thetad3+C_143*Thetad1*Thetad4+C_153*Thetad1*Thetad5+C_163*Thetad1*Thetad6+C_213*Thetad2*Thetad1+C_223*Thetad2*Thetad2+C_233*Thetad2*Thetad3+C_243*Thetad2*Thetad4+C_253*Thetad2*Thetad5+C_263*Thetad2*Thetad6+C_313*Thetad3*Thetad1+C_323*Thetad3*Thetad2+C_333*Thetad3*Thetad3+C_343*Thetad3*Thetad4+C_353*Thetad3*Thetad5+C_363*Thetad3*Thetad6+C_413*Thetad4*Thetad1+C_423*Thetad4*Thetad2+C_433*Thetad4*Thetad3+C_443*Thetad4*Thetad4+C_453*Thetad4*Thetad5+C_463*Thetad4*Thetad6+C_513*Thetad5*Thetad1+C_523*Thetad5*Thetad2+C_533*Thetad5*Thetad3+C_543*Thetad5*Thetad4+C_553*Thetad5*Thetad5+C_563*Thetad5*Thetad6+C_613*Thetad6*Thetad1+C_623*Thetad6*Thetad2+C_633*Thetad6*Thetad3+C_643*Thetad6*Thetad4+C_653*Thetad6*Thetad5+C_663*Thetad1*Thetad6)+g3;
taw4=(M_theta(4,1)*Thetadd1+M_theta(4,2)*Thetadd2+M_theta(4,3)*Thetadd3+M_theta(4,4)*Thetadd4+M_theta(4,5)*Thetadd5+M_theta(4,6)*Thetadd6)+(C_114*Thetad1*Thetad1+C_124*Thetad1*Thetad2+C_134*Thetad1*Thetad3+C_144*Thetad1*Thetad4+C_154*Thetad1*Thetad5+C_164*Thetad1*Thetad6+C_214*Thetad2*Thetad1+C_224*Thetad2*Thetad2+C_234*Thetad2*Thetad3+C_244*Thetad2*Thetad4+C_254*Thetad2*Thetad5+C_264*Thetad2*Thetad6+C_314*Thetad3*Thetad1+C_324*Thetad3*Thetad2+C_334*Thetad3*Thetad3+C_344*Thetad3*Thetad4+C_354*Thetad3*Thetad5+C_364*Thetad3*Thetad6+C_414*Thetad4*Thetad1+C_424*Thetad4*Thetad2+C_434*Thetad4*Thetad3+C_444*Thetad4*Thetad4+C_454*Thetad4*Thetad5+C_464*Thetad4*Thetad6+C_514*Thetad5*Thetad1+C_524*Thetad5*Thetad2+C_534*Thetad5*Thetad3+C_544*Thetad5*Thetad4+C_554*Thetad5*Thetad5+C_564*Thetad5*Thetad6+C_614*Thetad6*Thetad1+C_624*Thetad6*Thetad2+C_634*Thetad6*Thetad3+C_644*Thetad6*Thetad4+C_654*Thetad6*Thetad5+C_664*Thetad1*Thetad6)+g4;
taw5=(M_theta(5,1)*Thetadd1+M_theta(5,2)*Thetadd2+M_theta(5,3)*Thetadd3+M_theta(5,4)*Thetadd4+M_theta(5,5)*Thetadd5+M_theta(5,6)*Thetadd6)+(C_115*Thetad1*Thetad1+C_125*Thetad1*Thetad2+C_135*Thetad1*Thetad3+C_145*Thetad1*Thetad4+C_155*Thetad1*Thetad5+C_165*Thetad1*Thetad6+C_215*Thetad2*Thetad1+C_225*Thetad2*Thetad2+C_235*Thetad2*Thetad3+C_245*Thetad2*Thetad4+C_255*Thetad2*Thetad5+C_265*Thetad2*Thetad6+C_315*Thetad3*Thetad1+C_325*Thetad3*Thetad2+C_335*Thetad3*Thetad3+C_345*Thetad3*Thetad4+C_355*Thetad3*Thetad5+C_365*Thetad3*Thetad6+C_415*Thetad4*Thetad1+C_425*Thetad4*Thetad2+C_435*Thetad4*Thetad3+C_445*Thetad4*Thetad4+C_455*Thetad4*Thetad5+C_465*Thetad4*Thetad6+C_515*Thetad5*Thetad1+C_525*Thetad5*Thetad2+C_535*Thetad5*Thetad3+C_545*Thetad5*Thetad4+C_555*Thetad5*Thetad5+C_565*Thetad5*Thetad6+C_615*Thetad6*Thetad1+C_625*Thetad6*Thetad2+C_635*Thetad6*Thetad3+C_645*Thetad6*Thetad4+C_655*Thetad6*Thetad5+C_665*Thetad1*Thetad6)+g5;
taw6=(M_theta(6,1)*Thetadd1+M_theta(6,2)*Thetadd2+M_theta(6,3)*Thetadd3+M_theta(6,4)*Thetadd4+M_theta(6,5)*Thetadd5+M_theta(6,6)*Thetadd6)+(C_116*Thetad1*Thetad1+C_126*Thetad1*Thetad2+C_136*Thetad1*Thetad3+C_146*Thetad1*Thetad4+C_156*Thetad1*Thetad5+C_166*Thetad1*Thetad6+C_216*Thetad2*Thetad1+C_226*Thetad2*Thetad2+C_236*Thetad2*Thetad3+C_246*Thetad2*Thetad4+C_256*Thetad2*Thetad5+C_266*Thetad2*Thetad6+C_316*Thetad3*Thetad1+C_326*Thetad3*Thetad2+C_336*Thetad3*Thetad3+C_346*Thetad3*Thetad4+C_356*Thetad3*Thetad5+C_366*Thetad3*Thetad6+C_416*Thetad4*Thetad1+C_426*Thetad4*Thetad2+C_436*Thetad4*Thetad3+C_446*Thetad4*Thetad4+C_456*Thetad4*Thetad5+C_466*Thetad4*Thetad6+C_516*Thetad5*Thetad1+C_526*Thetad5*Thetad2+C_536*Thetad5*Thetad3+C_546*Thetad5*Thetad4+C_556*Thetad5*Thetad5+C_566*Thetad5*Thetad6+C_616*Thetad6*Thetad1+C_626*Thetad6*Thetad2+C_636*Thetad6*Thetad3+C_646*Thetad6*Thetad4+C_656*Thetad6*Thetad5+C_666*Thetad1*Thetad6)+g6;
taw=taw1 + taw2 + taw3 + taw4 +taw5 + taw6;

%% numerical
Theta1=deg2rad(0);Theta2=deg2rad(0);
Theta3=deg2rad(0);Theta4=deg2rad(0);
Theta5=deg2rad(0);Theta6=deg2rad(0);
Thetad1=0.1;Thetad2=0.1;Thetad3=0.1;
Thetad4=0.1;Thetad5=0.1;Thetad6=0.1;
Thetadd1=0.1;Thetadd2=0.1;Thetadd3=0.1;
Thetadd4=0.1;Thetadd5=0.1;Thetadd6=0.1;
Thetadd=[Thetadd1;Thetadd2;Thetadd3;Thetadd4;Thetadd5;Thetadd6];
%for I66 and m6==0
I6=zeros(3,3);m6=0;
taw_Lagrangian=(M_theta*Thetadd + c + g)
taw_Lagrangian_num=-eval(M_theta*Thetadd + c + g);
