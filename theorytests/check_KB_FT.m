% Check (deplinthed) Kaiser--Bessel Fourier transform pair formula.
% Barnett 11/2017, updated 4/2020 for dftsubmat paper repo, 2pi.i FT convention.
% needs: hline, vline, ft.

clear
sig = 20;
KB = @(z) besseli(0,sig*sqrt(1-z.^2)) .* (abs(z)<=1);  % K-B
tophat = @(z) (abs(z)<=1);                             % top hat
DKB = @(z) KB(z) - tophat(z);                          % deplinthed

DKBhat0 = @(w) 2*sinc(sqrt(w.^2-sig^2)) - 2*sinc(w);   % in usual FT convention

figure;
subplot(1,2,1); z = -1.1:1e-3:1.1;
h(2)=plot(z,DKB(z),'b-');
hold on; plot([-1 1], DKB([-1 1]), 'b.','markersize',10);
axis tight
v = axis;
xlabel('$x$','interpreter','latex'); ylabel('$f(x)$','interpreter','latex')
text(-1,v(4)*0.9,sprintf('$\\sigma=%g$',sig),'interpreter','latex','color',[0 0 0])
subplot(1,2,2);
wid = 3.0*sig;
k=-wid:1e-2:wid;  % symm
L=1.0;
DKBhat = real(ft(DKB,L,k));    % use quadrature for FT
semilogy(2*pi*k,abs([DKBhat; DKBhat0(k)]),'-');   % overplot theory and numer.
xlabel('$\xi$','interpreter','latex');         % xi = 2pi.k
ylabel('$\hat f(\xi)$','interpreter','latex')
axis tight; v = axis;
vline(2*pi*sig*[-1 1]); hline(2,'b:'); text(2*pi*sig*2,3,'2','color',[0 0 1]);
text(2*pi*sig+1,v(4)*0.7,'$\xi=2\pi\sigma$','interpreter','latex','color',[1 0 0]);

% math check
norm(DKBhat-DKBhat0(k))/norm(DKBhat)
