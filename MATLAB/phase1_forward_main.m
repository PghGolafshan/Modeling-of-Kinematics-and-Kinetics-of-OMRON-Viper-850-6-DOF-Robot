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
T0e=  T01 * T12 * T23 * T34 * T45 * T56 *T6e
%% user defined joint variant (theta) 
for n=1:6
 theta{n}=input(['theta',num2str(n),': (degrees) \n'],'s');
 theta{n}=str2sym(theta{n});
end
T1 = [cosd(theta{1}), -sind(theta{1}), 0, 0; sind(theta{1}),cosd(theta{1}), 0, 0; 0, 0, 1, 0.205;0, 0, 0, 1];
T2 = [cosd(theta{2}), -sind(theta{2}), 0, 0.075; 0,0, 1, 0; -sind(theta{2}), -cosd(theta{2}), 0, 0; 0, 0, 0, 1];
T3 = [cosd(theta{3}), -sind(theta{3}), 0, 0.385; sind(theta{3}),cosd(theta{3}), 0, 0; 0, 0, 1, 0; 0, 0, 0, 1];
T4 = [cosd(theta{4}), -sind(theta{4}), 0, 0; 0,0, -1, -0.295; sind(theta{4}), cosd(theta{4}), 0, 0; 0, 0, 0, 1];
T5 = [cosd(theta{5}), -sind(theta{5}), 0, 0; 0,0, -1, 0; sind(theta{5}), cosd(theta{5}), 0, 0; 0, 0, 0, 1];
T6 = [cosd(theta{6}), -sind(theta{6}), 0, 0; 0,0, 1, 0  ; -sind(theta{6}), -cosd(theta{6}), 0, 0; 0, 0, 0, 1];
Te=[1,0,0,0;0,1,0,0;0,0,1,0.08;0,0,0,1];
T0_e=  T1 * T2 * T3 * T4 * T5 * T6 * Te;
disp(['Transformation matrix of link 6 to base frame based on given joint variables:']);
eval(T0_e)
%% Euler,Fixed,equivalent Angles and Euler parameters for user defined joint valriables
for i=1:4
    for j=1:4
        r(i,j)=eval(T0_e(i,j));
    end
end
% Z-Y-Z euler angles
beta=atan2(sqrt(r(3,1)^2+r(3,2)^2),r(3,3));
if cos(beta)==0 && sin(beta)==1
    alpha=0;
    gamma=atan2(r(1,2),r(2,2));
elseif cos(beta)==0 && sin(beta)==-1
     alpha=0;
    gamma=-1*atan2(r(1,2),r(2,2));
else
alpha=atan2(r(2,3)/cos(beta),r(1,3)/cos(beta));
gamma=atan2(r(3,2)/sin(beta),-r(3,1)/sin(beta));
end
alpha;
beta;
gamma;
% X-Y-Z fixed angles
beta_fixed=beta;
alpha_fixed=alpha;
gamma_fixed=gamma;
% equivalent Angle-Axis
theta_equivalent=acos((r(1,1)+r(2,2)+r(3,3)-1)/2);
k_equivalent=(1/(2*sin(theta_equivalent)))*[r(3,2)-r(2,3) ; r(1,3)-r(3,1) ; r(2,1)-r(1,2)];
% Euler parameters
epsilon4=0.5*sqrt(1+r(1,1)+r(2,2)+r(3,3));
epsilon1=r(3,2)-r(2,3)/(4*epsilon4);
epsilon2=r(1,3)-r(3,1)/(4*epsilon4);
epsilon3=r(2,1)-r(1,2)/(4*epsilon4);
%% workspace
%2d
% Theta1=0;
% Theta4=0;
% wrist_position=T0e(1:3,4);
% for Theta2=deg2rad(-190):0.1:deg2rad(45) ;
%     for Theta3=deg2rad(-29):0.1:deg2rad(256) ;
%         for Theta5=deg2rad(-120):0.1:deg2rad(120) ;
%             x=gpuArray(eval(wrist_position(1)));
%             z=gpuArray(eval(wrist_position(3)));
%             plot(x,z,'o')
%             hold on
%             grid on  
%         end
%     end
% end
% %% 3d
% wrist_position=T0e(1:3,4);
% for Theta1=deg2rad(-170):0.1:deg2rad(170) ;
%     for Theta2=deg2rad(-190):0.1:deg2rad(45) ;
%         for Theta3=deg2rad(-29):0.1:deg2rad(256) ;
%             x=gpuArray(eval(wrist_position(1)));
%             y=gpuArray(eval(wrist_position(2)));
%             z=gpuArray(eval(wrist_position(3)));
%             plot3(x,y,z,'.')
%             hold on
%             grid on 
%         end
%     end
% end