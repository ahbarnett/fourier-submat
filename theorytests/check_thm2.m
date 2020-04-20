% sanity check main thm 2 in the *non*-asymptotic rate limit. It still applies.
% Barnett 4/13/20

clear
for N=2.^(2:6);   % up to 64, since beyond that, cond>1e16
  F = fft(eye(N));
  K = nan(N); Kt=K;
  for p = 1:N-1                 % sweep left half of submatrix sizes (p,q)
    for q = 1:p
      K(p,q) = cond(F(1:p,1:q));
      Kt(p,q) = (besseli(0,pi/2*(1-p/N)*q)-1) / (2*(sqrt(N/q)+6*sqrt(p*q)));
    end
  end
  fprintf('N=%d\t max cond=%.3g    \tcheck is >=0: %.3g\n',N,max(K(:)),min(log(K(:)) - log(Kt(:))))
  figure; imagesc(1:N,1:N,log(K)-log(Kt)); colorbar; title('log K - log Kthm2');
  xlabel('q'); ylabel('p'); axis tight equal ij;  % (1,0) corner bad by sqrt?
end
% looks good


