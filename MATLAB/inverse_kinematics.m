% Define the DH parameters and transformation matrices
syms theta1 theta2 theta3 theta4 theta5 theta6

T01 = [cos(theta1), -sin(theta1), 0, 0; 
       sin(theta1), cos(theta1), 0, 0; 
       0, 0, 1, 0.205; 
       0, 0, 0, 1];

T12 = [cos(theta2), -sin(theta2), 0, 0.075; 
       0, 0, 1, 0; 
       -sin(theta2), -cos(theta2), 0, 0; 
       0, 0, 0, 1];

T23 = [cos(theta3), -sin(theta3), 0, 0.385; 
       sin(theta3), cos(theta3), 0, 0; 
       0, 0, 1, 0; 
       0, 0, 0, 1];

T34 = [cos(theta4), -sin(theta4), 0, 0; 
       0, 0, -1, -0.295; 
       sin(theta4), cos(theta4), 0, 0; 
       0, 0, 0, 1];

T45 = [cos(theta5), -sin(theta5), 0, 0; 
       0, 0, -1, 0; 
       sin(theta5), cos(theta5), 0, 0; 
       0, 0, 0, 1];

T56 = [cos(theta6), -sin(theta6), 0, 0; 
       0, 0, 1, 0; 
       -sin(theta6), -cos(theta6), 0, 0; 
       0, 0, 0, 1];

T6e = [1, 0, 0, 0; 
       0, 1, 0, 0; 
       0, 0, 1, 0.08; 
       0, 0, 0, 1];

% Compute the overall transformation matrix symbolically
T0e = T01 * T12 * T23 * T34 * T45 * T56 * T6e;

% Define the joint angles (30 degrees for each joint) in radians
joint_angles_deg = 30 * ones(1, 6);
joint_angles_rad = deg2rad(joint_angles_deg);

% Substitute the joint angles into the transformation matrix
T0e_subs = subs(T0e, {theta1, theta2, theta3, theta4, theta5, theta6}, num2cell(joint_angles_rad));

% Evaluate the transformation matrix numerically
T0e_eval = double(T0e_subs);

disp('Forward kinematics (30 degrees for all joints):');
disp(T0e_eval);

% Define the desired end-effector transformation matrix based on the forward kinematics result
xd = T0e_eval(1, 4);
yd = T0e_eval(2, 4);
zd = T0e_eval(3, 4);
Rd = T0e_eval(1:3, 1:3);

% Create the desired transformation matrix
Td = [Rd, [xd; yd; zd]; 0, 0, 0, 1];

% Numerical IK solution using nonlinear optimization

% Define the error function for the optimization
error_function = @(angles) sum(ik_error_function(angles, Td).^2);

% Initial guess for the joint angles (radians)
initial_guess = zeros(1, 6);

% Define bounds for the joint angles (in radians)
lb = -pi * ones(1, 6);
ub = pi * ones(1, 6);

% Solve the inverse kinematics problem using fmincon
options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp');
[solution, fval, exitflag] = fmincon(error_function, initial_guess, [], [], [], [], lb, ub, [], options);

% Convert the solution to degrees
solution_deg = rad2deg(solution);

disp('Inverse kinematics solution (radians):');
disp(solution);

disp('Inverse kinematics solution (degrees):');
disp(solution_deg);

function error = ik_error_function(angles, Td)
    % Define the DH parameters and transformation matrices
    theta1 = angles(1);
    theta2 = angles(2);
    theta3 = angles(3);
    theta4 = angles(4);
    theta5 = angles(5);
    theta6 = angles(6);

    T01 = [cos(theta1), -sin(theta1), 0, 0; 
           sin(theta1), cos(theta1), 0, 0; 
           0, 0, 1, 0.205; 
           0, 0, 0, 1];

    T12 = [cos(theta2), -sin(theta2), 0, 0.075; 
           0, 0, 1, 0; 
           -sin(theta2), -cos(theta2), 0, 0; 
           0, 0, 0, 1];

    T23 = [cos(theta3), -sin(theta3), 0, 0.385; 
           sin(theta3), cos(theta3), 0, 0; 
           0, 0, 1, 0; 
           0, 0, 0, 1];

    T34 = [cos(theta4), -sin(theta4), 0, 0; 
           0, 0, -1, -0.295; 
           sin(theta4), cos(theta4), 0, 0; 
           0, 0, 0, 1];

    T45 = [cos(theta5), -sin(theta5), 0, 0; 
           0, 0, -1, 0; 
           sin(theta5), cos(theta5), 0, 0; 
           0, 0, 0, 1];

    T56 = [cos(theta6), -sin(theta6), 0, 0; 
           0, 0, 1, 0; 
           -sin(theta6), -cos(theta6), 0, 0; 
           0, 0, 0, 1];

    T6e = [1, 0, 0, 0; 
           0, 1, 0, 0; 
           0, 0, 1, 0.08; 
           0, 0, 0, 1];

    % Compute the overall transformation matrix numerically
    T0e_num = T01 * T12 * T23 * T34 * T45 * T56 * T6e;

    % Compute the position and orientation error
    position_error = T0e_num(1:3, 4) - Td(1:3, 4);
    orientation_error = T0e_num(1:3, 1:3) - Td(1:3, 1:3);

    % Combine the errors into a single error metric
    error = [position_error; reshape(orientation_error, 9, 1)];
end
