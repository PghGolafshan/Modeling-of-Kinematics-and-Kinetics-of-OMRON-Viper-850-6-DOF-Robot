clear
clc
syms t1 t2 t3 t4 t5 t6;
L(1)=Link([0,0.205,0,0],'modified');
L(2)=Link([0,0,0.075,-pi/2],'modified');
L(3)=Link([0,0,0.385,0],'modified');
L(4)=Link([0,0.295,0,pi/2],'modified');
L(5)=Link([0,0,0,pi/2],'modified');
L(6)=Link([0,0,0,-pi/2],'modified');
L(7)=Link([0,0.08,0,0],'modified');
L(1).qlim=[deg2rad(-170),deg2rad(170)]; 
L(2).qlim=[deg2rad(-190),deg2rad(45)];
L(3).qlim=[deg2rad(-29),deg2rad(256)];
L(4).qlim=[deg2rad(-190),deg2rad(190)];
L(5).qlim=[deg2rad(-120),deg2rad(120)];
L(6).qlim=[deg2rad(-360),deg2rad(360)];

R = SerialLink(L,'name','viper850')
R.teach
% T0e=vpa(R.fkine([t1 t2 t3 t4 t5 t6]));
% digits(3)
T=R.fkine([deg2rad(30) deg2rad(30) deg2rad(30) deg2rad(30) deg2rad(30) deg2rad(30) deg2rad(0)])
eul=tr2eul(T);