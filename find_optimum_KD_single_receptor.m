%% Find Optimum KD
%
% Single Receiver
%
% Inputs c1, c0, mu_noise, sigma_noise, NR
% clearvars

phi1 = 1e17;
phi0 = phi1*1e-2;
mu = 42;
sigma = 0.01;
NR = 120;

mean_noise = exp(mu+sigma^2/2);
median_noise = exp(mu);
mod_noise = exp(mu-sigma^2);

erf_derivative = @(x) -(2 / sqrt(pi)) * exp(-x.^2);
% pb = @(phi,n,KD) (exp(n)+phi)./(exp(n)+phi+(KD));

fun_mean    = @ (phi,n,KD) NR * (((exp(n)+phi))./(((exp(n)+phi))+KD)) .* (1./(sigma.*sqrt(2*pi))) .* exp(-((n-mu).^2)./(2*sigma^2));
fun_squared = @ (phi,n,KD) ((NR*((((exp(n)+phi)))./(((exp(n)+phi))+KD))).^2) .* (1./(sigma.*sqrt(2*pi))) .* exp(-((n-mu).^2)./(2*sigma^2));
fun_var     = @ (phi,n,KD) NR *((KD * (((exp(n)+phi))))./((((exp(n)+phi))+KD).^2)) .* (1./(sigma.*sqrt(2*pi))) .* exp(-((n-mu).^2)./(2*sigma^2));

total_variance = @(phi,n,KD) trapz(n,fun_var(phi,n,KD)) + trapz(n,fun_squared(phi,n,KD))-trapz(n,fun_mean(phi,n,KD))^2;

fun_mean_derivative = @(phi,n,KD) -NR * ((exp(n)+phi)./((exp(n)+phi+KD).^2)) .* (1./(sigma.*sqrt(2*pi))) .* exp(-((n-mu).^2)./(2*sigma^2));
fun_derivative_var1 = @ (phi,n,KD) NR * (((exp(n)+phi).* (exp(n)+phi-KD)) ./ ((exp(n)+phi+KD).^3)) .* (1./(sigma.*sqrt(2*pi))) .* exp(-((n-mu).^2)./(2*sigma^2));
fun_derivative_var2 = @ (phi,n,KD) -2 * (((NR * (exp(n)+phi)).^2) ./((exp(n)+phi+KD).^3)) .* (1./(sigma.*sqrt(2*pi))) .* exp(-((n-mu).^2)./(2*sigma^2));

total_derivative = @(phi,n,KD) trapz(n,fun_derivative_var1(phi,n,KD)) + trapz(n,fun_derivative_var2(phi,n,KD)) - 2 * trapz(n,fun_mean(phi,n,KD)) * trapz(n,fun_mean_derivative(phi,n,KD));

% lambda = NR * 0.5;
KD = sqrt(phi0/phi1)*(phi1);
noise = linspace(log(mod_noise/1000),log(mean_noise*1e6),1e4);
alpha = 0.5;%0.5
tol_KD = 1e-10;
tol_lambda = 1e-2;
KDs(1)   = KD;
KDler(1) = KD;

%% Adjust lambda according to KD
mean_c0 = trapz(noise,fun_mean(phi0,noise,KD));
mean_c1 = trapz(noise,fun_mean(phi1,noise,KD));
var_c0  = total_variance(phi0,noise,KD);
var_c1  = total_variance(phi1,noise,KD);
m0  = mean_c0;
m1  = mean_c1;
s0  = sqrt(2*var_c0); % s = sqrt(2 *var_c0)
s1  = sqrt(2*var_c1);
gamma = var_c1 - var_c0;
lambda = (gamma^-1) * (var_c1 * mean_c0 - var_c0 * mean_c1 + sqrt(var_c1)*sqrt(var_c0) * sqrt((mean_c1-mean_c0)^2 + gamma*log(var_c1/var_c0)));
lambdas(1)=lambda;
% lambda = 600;
prob_error(1) = probability_error_objective_fun(lambda,m0,m1,s0,s1);

for s = 1:10000
    s
    %% KD Gradient Descent
    % Statistics & Gradient
    mean_c0 = trapz(noise,fun_mean(phi0,noise,KD));
    mean_c1 = trapz(noise,fun_mean(phi1,noise,KD));
    var_c0  = total_variance(phi0,noise,KD);
    var_c1  = total_variance(phi1,noise,KD);
    m0  = mean_c0;
    m1  = mean_c1;
    s0  = sqrt(2*var_c0); % s = sqrt(2 *var_c0)
    s1  = sqrt(2*var_c1);
    dm0_dKD = trapz(noise,fun_mean_derivative(phi0,noise,KD));
    dm1_dKD = trapz(noise,fun_mean_derivative(phi1,noise,KD));
    ds0_dKD  = (1/s0)*total_derivative(phi0,noise,KD); % s = sqrt(2 *var_c0)
    ds1_dKD  = (1/s1)*total_derivative(phi1,noise,KD);
    
    % KD Update
    gradient_KD = calculate_derivative(lambda,m0,m1,s0,s1,dm0_dKD,dm1_dKD,ds0_dKD,ds1_dKD);
    KDold = KD;

%     if s==5000
%         alpha = alpha*5;
%     end
    KD = exp(log(KD) -(alpha/(1+log10(s))) * KD * gradient_KD);
    KDler(s+1) = KD;
    deltaKD(s) = KD-KDold;

    % Lambda Update
    gamma = var_c1 - var_c0;
    lambda = (gamma^-1) * (var_c1 * mean_c0 - var_c0 * mean_c1 + sqrt(var_c1)*sqrt(var_c0) * sqrt((mean_c1-mean_c0)^2 + gamma*log(var_c1/var_c0)));
    lambdas(s)=lambda;
    % Prob Error Objective Function
    prob_error(end+1) = probability_error_objective_fun(lambda,m0,m1,s0,s1);
end

% semilogy(1:s+1,prob_error,'LineWidth',2);
% yyaxis right
% semilogy(1:s+1,KDler,'LineWidth',2);
% grid on
% legend('Prob. Error','KD');
% yline(phi0,'LineWidth',2)
% hold on
% yline(phi1,'LineWidth',2)