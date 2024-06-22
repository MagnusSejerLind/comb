clc,clear,close all


%%
u_vec = linspace(0,2.5,1000);

for i_case = 1:2
for i = 1:length(u_vec)
    u = u_vec(i);


    %% Values

    l = 1; a = 0.25; b = 0.25; gamma_1 = 0.5; gamma_2 = 0.5; xi = 0.125; chi = 1.0; c_1 = 0.1; c_2 = c_1; d_1 = c_1; d_2 = c_1;

    if i_case == 1
    % case 1
    k_1 = 0.1; k_2 = 0.05; g_1 = 1; g_2 = 0.5;
    else
    % case 2
    k_1 = 1; k_2 = 0.5; g_1 = 0.1; g_2 = 0.05;
    end

    %%

    % Mass, damping, stiffness
    M = [2*(1/3+l) l^2 -1/2*(a-b) -l*(a-b);
        l^2 2/3*l^3 0 -1/2*l^2*(a-b);
        -1/2*(a-b) 0 2/3*(a^2+b^2-a*b) 0;
        -l*(a-b) -1/2*l^2*(a-b) 0 2/3*l*(a^2+b^2-a*b)];


    K = [k_1 0 chi * u ^ 2 * gamma_1 * (a + b) chi * u ^ 2 * l * (a + b);
        0 k_2 0 chi * u ^ 2 * l * gamma_2 * (a + b);
        0 0 g_1 - chi * u ^ 2 * xi * (a + b) 0;
        0 0 0 g_2 - chi * u ^ 2 * xi * l * (a + b)];


    C = [c_1 - chi * u * (gamma_1 + l) * (a + b) -chi * u * gamma_2 * l * (a + b) 0 0;
        -chi * u * gamma_2 * l * (a + b) c_2 - chi * u * gamma_2 ^ 2 * l * (a + b) 0 0;
        chi * u * gamma_1 * xi * (a + b) 0 d_1 0;
        u * chi * xi * l * (a + b) u * chi * xi * l ^ 2 * gamma_2 * (a + b) 0 d_2;];

    %%
    dof = 4;

    % Express A matrix
    A = [zeros(dof) eye(dof) ; -M\K -M\C];


    %% Check for flutter stability
    [Psi,Lambda]=eig(A);
    lambda=diag(Lambda); % Eigenvalues

    alpha(i,:) = sort(real(lambda));
    omega(i,:) = sort(imag(lambda));
    


end

%%

figure()
plot(u_vec,alpha,'k',LineWidth=2)
xlabel('u')
ylabel('\alpha')
grid

figure()
plot(u_vec,omega,'k',LineWidth=2)
ylabel('\omega')
xlabel('u')
grid


end