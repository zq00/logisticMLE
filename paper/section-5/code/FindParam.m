

% gamma = dlmread(strcat(fileloc,'gamma.txt'));

% Find solution to system in three unknowns

alpha0 = 1.5;
tau0 = sqrt(2)*alpha0*gamma;
lambda0 = 2;   

problem.options = optimoptions('fsolve','Display','none','PlotFcn',@optimplotfirstorderopt);
problem.objective = @(x) RootPragya(kappa,gamma,x(1),x(2),x(3));
problem.x0 = [alpha0, tau0, lambda0];
problem.solver = 'fsolve';
x = fsolve(problem);

[alpha_star, tau_star, lambda_star] = deal(x(1),x(2),x(3));
sigma_star = sqrt(tau_star^2 - alpha_star^2 * gamma^2)/sqrt(kappa);
%fprintf(' alpha^star: %2.4f\n', alpha_star)
%fprintf('sigma^star:  %2.4f\n\n', sigma_star)
%fprintf('lambda^star: %2.4f\n', lambda_star)

%write file
filename = strcat(fileloc,'param.txt');
fileID = fopen(filename,'w');
fprintf(fileID,'%f\t %f\t%f\n',[alpha_star, sigma_star, lambda_star]);
fclose(fileID);


function F = RootPragya(kappa,gamma,alpha,tau,lambda)

rho_prime = @(t) 1./(1+exp(-t));
rho_second = @(t) 1./(exp(-t/2) + exp(t/2)).^2;

options = optimset('Display','off');
prox = @(z) fsolve(@(t) lambda*rho_prime(t) + t - z, z,options);

% We now compute the density of Q_1, Q_2 which we express as the density of
% Z_1, Z_2 where Q_1 = gamma Z_1 and Q_2 = tau Z_2
% The correlation between Z_1 and Z_2 is alpha x gamma / tau

rho = alpha*gamma/tau;
density = @(z1,z2) 1/(2*pi*sqrt(1-rho^2)) .* exp(-(z1.^2 - 2*rho*z1.*z2 + z2.^2)/2/(1-rho^2));

% We now compute the density on a fine grid 
step_size = 0.025;
z1 = -8:step_size:8; 
z2 = -8:step_size:8; 
[Z1,Z2] = meshgrid(z1,z2);
D = density(Z1,Z2);

% We compute prox(Q_2) on a fine grid 
proxQ2 = prox(tau*z2);
% We evaluate rho' and rho'' at the prox 
rho_prime_proxQ2 = rho_prime(proxQ2);
rho_second_proxQ2 = rho_second(proxQ2);


% We evaluate the integrands on a fine grid 
make_tensor = @(z1,z2) z2' * z1;
trapezoidal_rule = @(X) step_size^2 * sum(X(:));

RHS1 = 1/kappa * trapezoidal_rule(...
    make_tensor(2*rho_prime(-gamma*z1), (lambda * rho_prime_proxQ2).^2) .* D);

RHS2 = trapezoidal_rule(...
    make_tensor(-2*rho_prime(-gamma*z1) * gamma .* z1, lambda * rho_prime_proxQ2) .* D);

RHS3 = trapezoidal_rule(...
    make_tensor(2*rho_prime(-gamma*z1), 1./(1 + lambda * rho_second_proxQ2)) .* D);

F(1) = tau^2 - alpha^2 * gamma^2 - RHS1;
F(2) = - RHS2;
F(3) = 1 - kappa - RHS3;
end

