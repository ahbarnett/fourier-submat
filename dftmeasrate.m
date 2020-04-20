function [r status] = dftmeasrate(n0,p0,q0,o)
% DFTMEASRATE  Measure empirical exponential growth rate of a submatrix shape
%
% r = dftmeasrate(n,p,q) fits condition number growth of a family of submatrices
%  of DFT matrices that share a common shape alpha := p/n and beta := q/n.
%  The fit is to the form exp(rho.N), where rho depends on the shape
%  (alpha, beta).
%
% [r status] = ...  also returns status = 0: if r good, else r unreliable
%                                1: even j=1 exceeds max cond (r is a LB)
%                                2: j=2 exceeds max cond, so just j=1 (")
%                                3: matrix predicted SVD size got too big
%
% [...] = dftmeasrate(...,opts) controls: opts.debug = 0,1,2
%
% Note: Internally, only N that is multiple of n are used, to get
%  (alpha,beta) exactly correct.
%
% Barnett 4/10/20

if nargin==0, test_dftmeasrate; return; end
if nargin<4, o=[]; end
if ~isfield(o,'debug'), o.debug=0; end

status = 0;           % ok
lcf = @(j) logcondmult(n0,p0,q0,j,o.debug);   % func handle log(cond(jth mult size A))

%if n0==50 & p0==9 & q0==9, r = (lcf(2)-lcf(1))/n0; return; end   % bypass hack

lcmax = log(1e13);    % max log cond # we can reliably measure (<1e16)
if n0>=50, rat = 1.2; box = 0.1; lcmaxx=log(1e16);  % corner zoom neff=50
else, rat = 5/4; box=0.25; lcmaxx=log(4e14); end    % main plot, neff=30

if p0/q0>1/rat & p0/q0<rat & p0/n0>box & p0/n0<1-box & q0/n0>box & q0/n0<1-box, lcmax = lcmaxx; end    % helps near diag, careful to adjust!

if o.debug>1, figure; end
stop = false; j = 1; lc = lcf(j);
if lc>lcmax, r=lc/n0; status=1; return; end  % only gives lower bnd for r
lcold = -inf;          % just to make sure lc>lcold 1st time around...
while ~stop & lc<=lcmax & lc>lcold   % find largest j st lcm(j) < lcmax
  lcold = lc;      % keep old val
  j=2*j;
  if max(p0,q0)*min(p0,q0)^2*j^3>5e8, stop = true; status=3;   % SVD too big
  else, lc = lcf(j);
  end
end
j0 = j/2; lc0 = lcold; j1 = j; lc1 = lc;  % brackets
while j1-j0>1
  jm = round((j1+j0)/2);   % midpoint
  lcm = lcf(jm);
  if lcm>lcmax             % take lower bracket
    j1=jm; lc1=lcm;
  else                     % upper bracket
    j0=jm; lc0=lcm;
  end
end
% done; j0 is below max, j1 above.
% use rate from this and half of it...
if j0==1, r = lc0/n0; status=2;        % only one valid size, not enough
  if o.debug, fprintf('(j,p,q,cond) used only: (%d,%d,%d,%.3g)\n',j0,j0*p0,j0*q0,exp(lc0)), end
  if o.debug>1, plot([0 j0], [0 lc0], '-'); end
else
  jm = floor(j0/2);   % half the size
  lcm = lcf(jm);
  r = (lc0-lcm)/(n0*(j0-jm));
  if o.debug, fprintf('(j,p,q,cond) used: (%d,%d,%d,%.3g) and (%d,%d,%d,%.3g)\n',jm,jm*p0,jm*q0,exp(lcm),j0,j0*p0,j0*q0,exp(lc0)), end
  if o.debug>1, plot([jm j0], [lcm lc0], '-'); end
end
if o.debug>1, hline(lcmax); end

function lc = logcondmult(n0,p0,q0,j,debug)
% returns log of cond # for j*p0 by j*q0 submatrix of the j*n0 DFT matrix
p=j*p0; q=j*q0; N=j*n0;
A = exp(2i*pi*(1:p)'*(1:q)/N);   % build p*q submat A, via outer prod
lc = log(cond(A));
if debug>1, fprintf('j=%d:\tlc=%.3g\n',j,lc), plot(j,lc,'+'); axis tight; hold on; v=axis; v([1 3])=0; axis(v); end


%%%%%%%%%%%%%%%
function test_dftmeasrate
% simple tests...
if 0
o.debug = 1;
r1 = dftmeasrate(10,1,1,o)
r2 = dftmeasrate(20,2,2,o)
r1-r2

o.debug=2; [r st] = dftmeasrate(30,15,15,o)
%o.debug=2; r1 = dftmeasrate(30,29,3,o)
%r1 = dftmeasrate(30,16,5,o)

stop

o.debug=2;
r1 = dftmeasrate(50,12,12,o)
r2 = dftmeasrate(25,6,6,o)
r1-r2
stop
end

% test looping over (alpha,beta)...
na = 200;   % 1 / step size in alpha or beta
as = (1:20)/na;
r = nan(numel(as)); st=r; r2=r; st2=r;   % data arrays
for a=1:numel(as), p0=as(a)*na;
  for b=1:a, q0=as(b)*na;
    [r(a,b) st(a,b)] = dftmeasrate(na,p0,q0);
    if mod(p0,2)==0 & mod(q0,2)==0
      [r2(a,b) st2(a,b)] = dftmeasrate(na/2,p0/2,q0/2); end
  end
end
figure;
subplot(2,2,1); imagesc(as,as,r); colorbar; axis equal tight ij; c = caxis;
xlabel('\beta'); ylabel('\alpha'); title('r rate');
subplot(2,2,2); imagesc(as,as,r2-r); colorbar; axis equal tight ij; 
title('r-r2 (err estimate of r)'), xlabel('\beta'); ylabel('\alpha');
fprintf('max diff r vs r2: %.3g\n',max(max(abs(r-r2))))
subplot(2,2,3); imagesc(as,as,st); colorbar; axis equal tight ij; 
title('status of each r call'), xlabel('\beta'); ylabel('\alpha');
subplot(2,2,4); imagesc(as,as,st2); colorbar; axis equal tight ij; 
title('status of each r2 call'), xlabel('\beta'); ylabel('\alpha');

