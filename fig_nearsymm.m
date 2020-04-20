% figure showing near-symmetry in the cond # of p*q submatrix of DFT,
% wrt p -> N-p and q -> N-q.
% Barnett 3/20/20

clear
if 1      % log cond A plot over (p,q) for several N .....................
Ns = [8 16 32]; nn=numel(Ns);
figure;
for i=1:nn, N = Ns(i); F = fft(eye(N));    % loop over N's
  K = nan(N,N);
  for p=1:N, for q = p:N          % sweep all submatrix sizes (p,q)
      K(p,q) = cond(F(1:p,1:q));
      K(q,p) = K(p,q);   % get this one for free by invariance under adjoint
    end, end
  axes('position',[0.05+(i-1)/3,0,0.28,1]);  % hack my own subplot loc
  imagesc(log10(K));
  xlabel('$q$','interpreter','latex');
  ylabel('$p$','interpreter','latex');
  colormap(goodbw);  %colormap(jet(256));
  colorbar %('southoutside')
  title(sprintf('$\\log_{10}$ cond$(A)$: $\\;N$=%d',N),'interpreter','latex');
  axis equal tight ij;
  text(0.5-N/4,0.5,['(' char(96+i) ')'],'fontsize',12);
  %  print('-depsc2','-loose',sprintf('nearsymm%d.eps',i));  % loose fixes eps bbox
end
set(gcf,'paperposition',[0 0 8 2.4]);
% after dvipdf this fails to display right in evince (acroread fine): :(
%print -depsc2 -loose nearsymm.eps
%print -depsc2 nearsymm.eps
%print -dtiff -r300 -loose nearsymm.tiff
print -dtiff -r300 nearsymm.tiff
system('convert nearsymm.tiff eps3:nearsymm.eps');  % EPS is 2MB, bad
% but seems to compress down when dvipdf done.
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0    % mysterious 5/2 power growth of log on antidiagonal .................
  N = 128;
  F = fft(eye(N));
  K = nan(N-1,N-1);
  for p=1:N-1, for q = p:N-1
      K(p,q) = cond(F(1:p,1:q));
      K(q,p) = K(p,q);   % get this one for free by invariance under adjoint
    end, end
    
  ii = 1:N-1;
  %d = nan*ii; for i=ii, d(i) = K(i,N-i); end   % antidiagonal
  q = round(0.75*N);
  d = K(q,:);       % horiz slice
  mld =max(log(d));
  figure; loglog(ii, log(d), '+-'); hold on;
  %  hold on; semilogy(ii,exp(2.0 * 1.4* min(ii,N-ii).^2.5/N^1.5),'r-'); % 5/2 pow!
g=8/7; loglog(ii,ii.^g*(log(d(N/4))/(N/4)^g),'r-');
g=1.0; loglog(ii,ii.^g*(log(d(N/4))/(N/4)^g),'m-');
g=6/5; loglog(ii,ii.^g*(log(d(N/4))/(N/4)^g),'k--');
%g=1.0; loglog(ii,min(ii,N-ii).^g*(mld/(N/2)^g),'m-');
%g=1.5; loglog(ii,min(ii,N-ii).^g*(mld/(N/2)^g),'g-');
%g=2.7; loglog(ii,min(ii,N-ii).^g*(mld/(N/2)^g),'m-');
%g=2.5; loglog(ii,min(ii,N-ii).^g*(mld/(N/2)^g),'g-');
%hold on; loglog(ii,2.0 * min(ii,N-ii).^2/N,'r-');
  axis tight; axis([3 q+2 1 mld]);
  title(sprintf('$\\log_{10}$ cond $A$: $\\;N$=%d',N),'interpreter','latex');
  % ** text (d)
end

% Results: 1<=q<=N/2 slice at p=N/2:  cond ~ exp( c q^g)
% power g is not stable as N grows:
% N=128: best power is 8/7, fitted through N/4
% N=96: best power is 7/6, fitted through N/4
% N = 64:     6/5
% N=48:  5/4

% if pow tends to 1, good, that matches expected asymptotics. We are seeing
% pre-asymp behavior.

% not really universal power,. Also it dep on the q slice, eg q=0.75N different

