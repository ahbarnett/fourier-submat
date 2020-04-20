% check Moitra's Fejer kernel vector (Fejer implies periodic, 2pi-per usually)
% Barnett 4/1/20

% Moitra's proof is pretty loose, esp the application of Poisson summation.

clear
N=100;
l=10;  % initial width (bandwidth of Hhat)
r = 3;   % power
J=-N/2:N/2-1;
ker = @(x) (sin(l*pi*x)./sin(pi.*x)/l).^2;   % note 1-per in x, continuous var.
x = J/N;
H = ker(x).^r; H(x==0) = 1.0;  % rescue dumb nan
figure;
subplot(2,1,1); plot(J,H,'+'); v = axis; v(3)=1e-15;
hold on; plot(J,1./(4*(l*x).^(2*r)),'r-'); axis(v);
set(gca,'ysc','log');
Hhat = fftshift(fft(fftshift(H)));
subplot(2,1,2); semilogy(J,abs(Hhat),'+');   % note, compact support size 2rl
