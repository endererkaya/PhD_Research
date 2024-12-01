function out = conditional_total_variance_quadruple_receptor(c,KD,mu,sigma,NR)

N = 1e3;
% c = 1e21;
% KD = 1e21;
% mu = 49.04;
% sigma = 2;

KD1 = KD / 1e3;
KD2 = KD / 1e1;
KD3 = KD * 1e1;
KD4 = KD * 1e3;

pb  = @ (n) 1/4 * (((exp(n)+exp(c))./(exp(n)+exp(c)+KD1)) + ((exp(n)+exp(c))./(exp(n)+exp(c)+KD2)) + ((exp(n)+exp(c))./(exp(n)+exp(c)+KD3)) + ((exp(n)+exp(c))./(exp(n)+exp(c)+KD4)));

% fun_mean    = @ (n) pb(n) .* (1./(n.*sigma.*sqrt(2*pi))) .* exp(-((log(n)-mu).^2)./(2*sigma^2));
% fun_squared = @ (n) (pb(n).^2) .* (1./(n.*sigma.*sqrt(2*pi))) .* exp(-((log(n)-mu).^2)./(2*sigma^2));
% fun_var     = @ (n) (pb(n).*(1-pb(n))) .* (1./(n.*sigma.*sqrt(2*pi))) .* exp(-((log(n)-mu).^2)./(2*sigma^2));

% fun_mean    = @ (n) pb(n) .* (1./(sigma.*sqrt(2*pi))) .* exp(-((n-mu).^2)./(2*exp(sigma^2)));
% fun_squared = @ (n) (pb(n).^2) .* (1./(sigma.*sqrt(2*pi))) .* exp(-((n-mu).^2)./(2*exp(sigma^2)));
% fun_var     = @ (n) (pb(n).*(1-pb(n))) .* (1./(sigma.*sqrt(2*pi))) .* exp(-((n-mu).^2)./(2*exp(sigma^2)));

fun_mean    = @ (n) (NR * pb(n)) .* (1./(sigma.*sqrt(2*pi))) .* exp(-((n-mu).^2)./(2*sigma^2));
fun_squared = @ (n) ((NR * pb(n)).^2) .* (1./(sigma.*sqrt(2*pi))) .* exp(-((n-mu).^2)./(2*sigma^2));
fun_var     = @ (n) ((NR * pb(n)).*(1-pb(n))) .* (1./(sigma.*sqrt(2*pi))) .* exp(-((n-mu).^2)./(2*sigma^2));

% noise = logspace(-20,log10(exp(mu+(sigma^2)/2+10*sigma)),N);
% noise = logspace(10,32,N);
noise = linspace(1,100,N);

total_variance = @(n) trapz(n,fun_var(n)) + trapz(n,fun_squared(n))-trapz(n,fun_mean(n))^2;

out = total_variance(noise);