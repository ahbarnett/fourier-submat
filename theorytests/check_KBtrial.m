% check deplinthed Kaiser-Bessel trial vector has correct DFT formula.
% Barnett 3/30/20, based on check_prop6.m Gaussian case.
% 4/6/20 deplinthed case.
% Needs: sinc, vline

clear
N = 32;             % DFT size
q = 21;              % submat size p*q
p = 15;
s = (pi/2)*q*(1-p/N);   % sigma, the freq param in KB
fprintf('predicted min-sing-val order 1/J_0(sigma) = %.3g\n',1/besseli(0,s))
J = ceil(-N/2):ceil(N/2)-1;   % handles N even or odd
f = real(besseli(0,s*sqrt(1-(2*J/q).^2))) - 1;   % deplinthed (the "-1") KB
f(~(J>=-q/2 & J<q/2)) = 0;                       % give it the cpt supp
%figure; plot(J,f,'+-'); title('f vector')
%f(abs(2*J)==q) = 0.5;       % fix the endpoints - still not enough
Ff = 0*J;     % k indices same as J ones
nt = 3e4;     % # terms: error = O(1/nt^2)   (note 1e5 gets 1e-14 rel digits)
Ffmax = nan(2*nt+1,1);    % measure how big each term is
for m=-nt:nt
  x = pi*q*(J/N+m);     % algebraic progression of freqs
  % this has catastrophic cancellation for large N and m... careful...
  Ffm = q*( sinc(sqrt(x.^2 - s^2)) - sinc(x) ); % this m term in DKB DFT formula
  Ff = Ff + Ffm;
  Ffmax(m+nt+1) = max(abs(Ffm));
end
fprintf('max abs (rel) of last term included (m=%d): %.3g (%.3g)\n',nt,Ffmax(end),Ffmax(end)/max(abs(Ff)))
% note for N odd, fftshift is off by one, hence the hacked circshift...
Ff0 = real(fftshift(fft(circshift(f,ceil(N/2)))));      % real since symmetric
fprintf('max rel err btw theory & num: %.3g\n',norm(Ff-Ff0,inf)/norm(Ff,inf))
% got at N=20, q=7, p=11, nt=1e6: 1e-14 error. I think that proves it!

figure; subplot(3,1,1); semilogy(J,f,'+'); title('trial vector f_j');
axis tight; v=axis; axis([-N/2 N/2 v(3:4)]); xlabel('j');
vline(q/2*[-1 1],'r:','q/2');
subplot(3,1,2); semilogy(J,abs([Ff;Ff0]), '+');
legend('(Ff)_k formula', '(Ff)_k numer.'); title('DFT of trial vector');
axis tight; v=axis; axis([-N/2 N/2 v(3:4)]); xlabel('k');
vline((N-p)/2*[-1 1],'r:','(N-p)/2');
subplot(3,1,3); plot(J,Ff-Ff0, '+-'); title('abs err in Ff: theory-numer.');
axis tight; v=axis; axis([-N/2 N/2 v(3:4)]); xlabel('k');
