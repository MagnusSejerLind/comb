clc,clear,close all

%%
syms m x y a b l_1 l_2 phi_1 phi_1_dot phi_1_dotdot phi_2 phi_2_dot phi_2_dotdot theta_1 theta_1_dot theta_1_dotdot theta_2 theta_2_dot theta_2_dotdot 
syms k_1 k_2 g_1 g_2 c_1 c_2 d_1 d_2
syms C  % Proportional constant
syms rho U gamma_1 gamma_2 xi

%% Equations

% Kinetic energy (Modified)
T = m*(1/3*l_1*(a^3+b^3)*theta_1_dot^2 + 1/3*l_1^3*(a+b)*phi_1_dot^2 - 1/2*l_1^2*(a^2-b^2)*theta_1_dot*phi_1_dot + ...
    l_1^2*l_2*(a+b)*phi_1_dot^2 + 1/3*l_2^3*(a+b)*phi_2_dot^2 + 1/3*l_2*(a^3+b^3)*theta_2_dot^2 + ...
    l_1*l_2^2*(a+b)*phi_1_dot*phi_2_dot - l_1*l_2*(a^2-b^2)*phi_1_dot*theta_2_dot - 1/2*l_2^2*(a^2-b^2)*phi_2_dot*theta_2_dot   ...
    );

% Potential energy
V = 1/2*( k_1*phi_1^2 + g_1*theta_1^2 + k_2*phi_2^2 + g_2*theta_2^2 );

% Generalized coordinates
q = [phi_1, phi_2, theta_1, theta_2];
q_dot = [phi_1_dot, phi_2_dot, theta_1_dot, theta_2_dot];

% Forces
F_i = [-C*rho/2*U^2*(theta_1-1/U*phi_1_dot*gamma_1)*(a+b)*l_1, ...
-C*rho/2*U^2*(theta_2-1/U*phi_1_dot*gamma_1+phi_2_dot*gamma_2)*(a+b)*l_1];

f_trans = [-sin(theta_1), -sin(theta_1), cos(theta_1)*cos(phi_1);
            -sin(theta_2), -sin(theta_2), cos(theta_2)*cos(phi_2)] ;
for i = 1:length(F_i)
    F_bold(i,:) = F_i(i).*f_trans(i,:);
end

% Position vectors
r_1 = [-xi*cos(theta_1), gamma_1*cos(phi_1), x*1i*sin(theta_1)+gamma_1*sin(phi_1)];
r_2 = [-xi*cos(theta_2), l_1*cos(phi_1)+gamma_2*cos(phi_2), -xi*sin(theta_2)+l_1*sin(phi_1)+gamma_2*sin(phi_2)];

r_bold = r_1;
r_bold(2,:) = r_2;

% Damping
D = 1/2*( c_1*phi_1_dot^2 + d_1*theta_1_dot^2 + c_2*phi_2_dot^2 + d_2*theta_2_dot^2 );

%% Set up the Lagrange equation

% Langrian
L = T - V;

% Derivatives
for j = 1:4
    % Partial derivative of Lagrian with respect to d(q_i)/dt
    dL_dq_dot(j) = diff(L,q_dot(j));


    % Partial derivative of Lagrian with respect to q_i
    dL_dq(j) = diff(L,q(j));


    % Partial derivative of damping with respect to d(q_i)/dt
    dD_dq_dot(j) = diff(D,q_dot(j));

end


% Express generalized forces Q_j
for j = 1:4
    r_bold_dq_k1(j,:) = diff(r_bold(1,:),q(j));
    r_bold_dq_k2(j,:) = diff(r_bold(2,:),q(j));
end
r_bold_dq = [r_bold_dq_k1 ; r_bold_dq_k2];

for j = 1:length(q)
    Q(j,:) = dot(F_bold(1,:)', r_bold_dq_k1(j,:)) + dot(F_bold(2,:)', r_bold_dq_k2(j,:));
end



% Temporal derivative of partial derivative of Lagrian with respect to d(q_i)/dt
%%%%%%%%%%%%% OBS MANUAL %%%%%%%%%%%%%%%
% phi_dot->phi_dotdot - likewise for theta
ddL_dq_dot_dt = [m*(0.6667*l_1^3*phi_1_dotdot*(a + b) - 0.5000*l_1^2*theta_1_dotdot*(a^2 - b^2) + ...
    2*l_1^2*l_2*phi_1_dotdot*(a + b) + l_1*l_2*phi_2_dotdot*(a + b)), m*(0.6667*l_2^3*phi_2_dotdot*(a + b) - ...
    0.5000*l_2^2*theta_2_dotdot*(a^2 - b^2) + l_1*l_2*phi_1_dotdot*(a + b) - ...
    l_1*l_2*theta_2_dotdot*(a^2 - b^2)), m*(0.3333*l_1*(a^3 + b^3) - ...
    0.5000*l_1^2*phi_1_dotdot*(a^2 - b^2)), -m*(0.5000*l_2^2*phi_2_dotdot*(a^2 - b^2) - ...
    0.6667*l_2*theta_2_dotdot*(a^3 + b^3) + l_1*l_2*phi_2_dotdot*(a^2 - b^2))];
%%%%%


%% Lagrange equation

% Q excluded 
LagEq_homo = ddL_dq_dot_dt - dL_dq + dD_dq_dot;


% with Q_j
LagEq = ddL_dq_dot_dt - dL_dq + dD_dq_dot - Q.';

%%

EOM_1 = LagEq(:,1);
EOM_2 = LagEq(:,2);
EOM_3 = LagEq(:,3);
EOM_4 = LagEq(:,4);







