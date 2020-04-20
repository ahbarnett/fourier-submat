% check Prop 6: periodized Gaussian DFT pair
N=100; s = 5;
J = -N/2:N/2-1;
f = 0*J;
for n=-10:10, f = f + exp(-0.5*(J+n*N).^2/s^2); end  % approx periodic sum
Ff = 0*J;     % k indices same as J ones
for n=-10:10, Ff = Ff + sqrt(2*pi)*s*exp(-2*(pi*s/N)^2*(J+n*N).^2); end
Ff0 = fftshift(fft(fftshift(f)));
figure; subplot(2,1,1); plot(J,f,'+'); title('f_j');
subplot(2,1,2); plot(J,[Ff; Ff0], '+'); legend('Ff formula', 'Ff numerical');
norm(Ff - Ff0)   % check, good
