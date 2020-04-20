% see if the trial vector in Thm 2 is close to the actual min sing vec.
% Barnett 4/13/20.
clear
N=64;
F = fft(eye(N));
p=N/2, q=N/2-1;   % near diag, where gets nonsharp by factor 1.45 or so.
assert(q<=p,'must be tall or square')    % for trial vec to be in v (right)
[U S V] = svd(F(1:p,1:q));
sig1 = S(1,1); sigmin = S(q,q);
fprintf('N=%d, p=%d, q=%d: cond = %.3g\n',N,p,q,sig1/sigmin)  % cond
vmin = V(:,end); umin = U(:,end);
Ff = F(:,1:q)*vmin;
s = pi/2*(1-p/N)*q;   % our sigma width param
% really, edges can go out to +-(q+1)/2 where hits zero
s = pi/2*(1-p/N)*(q-2);   % slightly optimized sigma width param
J = floor(-q/2+.75):floor(q/2-.25);    % handles both even and odd cases
v0 = besseli(0,s*sqrt(1-(2*J/q).^2)) - 1;  % our DKB trial vec
v0 = v0'/norm(v0);  % unit norm col vec
figure(1); clf
subplot(2,1,1); semilogy([abs(vmin), v0],'+-');
legend('right sing vec v','Thm 2 KB trial vec');
axis tight;
subplot(2,1,2); semilogy(abs(Ff),'+-'); legend('full Ff output vec');
axis tight;
