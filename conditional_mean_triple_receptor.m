function out = conditional_mean_triple_receptor(c,KD,mu,sigma,NR)

N = 1e4;
% c0 = 1e21;
% KD = 1e21;
% mu = 49.04;
% sigma = 2;

KD1 = KD / 1e2;
KD2 = KD;
KD3 = KD * 1e2;

pb  = @ (n) 1/3 * (((exp(n)+exp(c))./(exp(n)+exp(c)+KD1)) + ((exp(n)+exp(c))./(exp(n)+exp(c)+KD2)) + ((exp(n)+exp(c))./(exp(n)+exp(c)+KD3)));
% fun = @ (n) pb(n) .* (1./(n.*sigma.*sqrt(2*pi))) .* exp(-((log(n)-mu).^2)./(2*sigma^2));
fun  = @ (n) NR * pb(n) .* (1./(sigma.*sqrt(2*pi))) .* exp(-((n-mu).^2)./(2*sigma^2));

% noise = logspace(-20,log10(exp(mu+(sigma^2)/2+10*sigma)),N);
noise = linspace(0,100,N);

out = trapz(noise,fun(noise));