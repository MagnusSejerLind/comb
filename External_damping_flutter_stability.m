clc,clear,close all
set(0,'defaultTextInterpreter','latex');
%%

n = 1000;
u_vec = linspace(0,2.5,n);  % Velocity range

dof = 4;
alpha = zeros(n,2*dof);
omega = zeros(n,2*dof);


for i_case = 1:2
    for i = 1:length(u_vec)
        u = u_vec(i);

        %% Non-dimensional parameters
        l = 1; a = 0.25; b = 0.25; 
        gamma_1 = 0.5; gamma_2 = 0.5; xi = 0.125; 
        chi = 1.0; 
        c_in_1 = 0.1; c_in_2 = c_in_1; d_in_1 = c_in_1; d_in_2 = c_in_1;    % Internal damping
        c_hat = 0.1    % External damping
        
        % Spring constants
        if i_case == 1
            % case 1
            k_1 = 0.1; k_2 = 0.05; 
            g_1 = 1; g_2 = 0.5;
        else
            % case 2
            k_1 = 1; k_2 = 0.5; 
            g_1 = 0.1; g_2 = 0.05;
        end

        %%

        % Mass, damping, stiffness matrices
        M = [2*(1/3+l) l^2 -1/2*(a-b) -l*(a-b);
            l^2 2/3*l^3 0 -1/2*l^2*(a-b);
            -1/2*(a-b) 0 2/3*(a^2+b^2-a*b) 0;
            -l*(a-b) -1/2*l^2*(a-b) 0 2/3*l*(a^2+b^2-a*b)];

        C = [c_in_1 - chi * u * (gamma_1 + l) * (a + b) -chi * u * gamma_2 * l * (a + b) 0 0;
            -chi * u * gamma_2 * l * (a + b) c_in_2 - chi * u * gamma_2 ^ 2 * l * (a + b) 0 0;
            chi * u * gamma_1 * xi * (a + b) 0 d_in_1 0;
            u * chi * xi * l * (a + b) u * chi * xi * l ^ 2 * gamma_2 * (a + b) 0 d_in_2;];

        K = [k_1 0 chi * u ^ 2 * gamma_1 * (a + b) chi * u ^ 2 * l * (a + b);
            0 k_2 0 chi * u ^ 2 * l * gamma_2 * (a + b);
            0 0 g_1 - chi * u ^ 2 * xi * (a + b) 0;
            0 0 0 g_2 - chi * u ^ 2 * xi * l * (a + b)];


        %% Express A matrix
        A = [zeros(dof) eye(dof) ; -M\K -M\C];

        %% Flutter stability
        [Psi,Lambda]=eig(A);    % Eigenvalue problem
        lambda=diag(Lambda);    % Eigenvalues extracted

        alpha(i,:) = sort(real(lambda));    % Real part
        omega(i,:) = sort(imag(lambda));    % Imaginary part



    end

    %%

    figure()
    plot(u_vec,alpha,'k',LineWidth=1.75)
    xlabel('$u$')
    ylabel('$\alpha$',Rotation=360)
    title("Case",i_case)
    grid
    yline(0,'r',LineWidth=1)


    figure()
    plot(u_vec,omega,'k',LineWidth=1.75)
    ylabel('$\omega$',Rotation=360)
    xlabel('$u$')
    title("Case",i_case)
    grid
    yline(0,'r',LineWidth=1)



end