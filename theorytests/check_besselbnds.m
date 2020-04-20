% Barnett see how good Bessel bnds are. 3/21/20

% Siegel UB on J_n(nx) :
figure; x=(1:1e3)/1e3;y=sqrt(1-x.^2);plot(x,log(x.*exp(y)./(1+y)),'-');
hold on; n=1e2; plot(x,log(besselj(n,x*n))/n,'r-');   % check it, large order
hold on; n=1e3; plot(x,log(besselj(n,x*n))/n,'g-');

% z = argument/n in Bessel:
y = @(z) sqrt(1-z.^2);
g = @(z) z.*exp(y(z))./(1+y(z));    % what gets taken to power n in Siegel.
z = logspace(-10,0,11)';
g(z)./z
exp(1)/2  % slope of upper linear bound, empirically.

