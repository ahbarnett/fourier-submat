% check proposition on sing vals and cond # of complementary subblocks
% 3/20/20

clear
if 0     % single example, detailed check of proof steps
N = 60; p = 20; q = 15;
fprintf('sqrt N = %g\n',sqrt(N))
F = fft(eye(N));
A = F(1:p,1:q);      % contiguous submatrix
B = F(1:p,q+1:end); C = F(p+1:end,1:q); D = F(p+1:end,q+1:end);  % other blocks
sA = svd(A); sB = svd(B); sC = svd(C); sD = svd(D);
fprintf('A:\t%.16g\t%.16g\n',max(sA),min(sA))
fprintf('B:\t%.16g\t%.16g\n',max(sB),min(sB))
fprintf('C:\t%.16g\t%.16g\n',max(sC),min(sC))
fprintf('D:\t%.16g\t%.16g\n',max(sD),min(sD))

max(sA)^2 + min(sC)^2 - N   % check steps in pf
min(sA)^2 + max(sC)^2 - N

max(sC)^2 + min(sD)^2 - N
min(sC)^2 + max(sD)^2 - N


cond(A)
cond(D)

'check bnds:'
cond(D)/cond(A) - (1-min(sC)^2/N)^-.5  % good
return
end

N = 65;           % sweep cond over (p,q) plane
F = fft(eye(N));
K = nan(N-1,N-1);
for p=1:N-1, for q = p:N-1
    K(p,q) = cond(F(1:p,1:q));
    K(q,p) = K(p,q);   % get this one for free by invariance under adjoint
  end, end
K
figure;
subplot(2,2,1); imagesc(K); colorbar; axis equal tight;
xlabel('q'); ylabel('p');
title(sprintf('N=%d: cond submatrices of F',N))
subplot(2,2,2); imagesc(K-flipud(fliplr(K))); colorbar; axis equal tight;
xlabel('q'); ylabel('p');
title(sprintf('cond diff under inversion symm'))
subplot(2,2,3); imagesc(log(K)); colorbar; axis equal tight;
xlabel('q'); ylabel('p');
title(sprintf('N=%d: log cond submatrices of F',N))
subplot(2,2,4); imagesc(log(K)-flipud(fliplr(log(K)))); colorbar; axis equal tight;
xlabel('q'); ylabel('p');
title(sprintf('log cond diff under inversion symm'))

if 1
ii =1:N-1;
d = nan*ii; for i=ii, d(i) = K(i,N-i); end 
figure; semilogy(ii, d, '+-'); title('skew diagonal slice');
hold on; semilogy(ii,exp(2.0 * 1.4* min(ii,N-ii).^2.5/N^1.5),'r-'); % 5/2 pow!
end

if 0   % fit a model
  figure;
end

% demo a bad case for paper
%F=fft(eye(50)); p=5; q=6; cond(F(1:p,1:q)), cond(F(p+1:end,q+1:end))

% check statement in symm chat: sigma_min of A
N = 60; fprintf('sqrt N = %g\n',sqrt(N))
F = fft(eye(N)); A = F(2:N,2:N); v=ones(N-1,1); A*v, norm(A*v)/norm(v)
