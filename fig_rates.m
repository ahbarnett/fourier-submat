% figure showing actual and theoretical exponential cond(submatrix) rates rho.
% Barnett 2/28/30 - 3/20/20. Restart 4/9/20.
% Needs: goodbw

clear
ox = .07; oy = .04; sx = .4; sy = .4;  % layout needed for all subplots

if 1 % MAIN SQUARE (0,1)^2 PLOT ..............................................
if 1  % regenerate data - 20 sec or so.
n = 30; as = (0:n)/n;  % sample set in a and in b
g = nan(numel(as)); gt=g; st=g;
for a=2:n, alpha=as(a)
  for b=2:a, beta=as(b);
    gt(a,b) =  pi/2*(min(alpha,beta)-alpha*beta);  % KB theorem rate, Thm 2
    %gt(a,b) =  beta*log(cot(pi*alpha/4));    % prolate Slepian V asym
    [g(a,b) st(a,b)] = dftmeasrate(n,alpha*n,beta*n);  % new rate meas
    g(b,a) = g(a,b); gt(b,a) = gt(a,b);
  end
end
g(isnan(g)) = 0.0;    % fill the zero-size results w/ zero rate
%save rates2.mat n as g st gt
else, load rates2.mat; end
a=find(as==0.5); fprintf('ratio at (.5,.5): %.5g\n',g(a,a)/gt(a,a))
da=as(2)-as(1); rc=1.2;   % cutoff ratio to count area frac (avg w & w/o edges)
af = (sum(g(:)./gt(:)<rc)+sum(sum(g(2:n,2:n)./gt(2:n,2:n)<rc)))/2*da^2;
fprintf('area frac below %g is %.3g\n',rc,af)
% figure(2); imagesc(as,as,st); colorbar;   % check status of g meas

figure(1); clf
v = [0:0.05:0.6]; vl = 0.1:0.1:0.5;  % show contour lines, & to label, (a-b)
axes('position',[ox oy+.5 sx sy]);
[C,h] = contourf(as,as,g,v);
clabel(C,h,vl,'labelspacing',600);
ylabel('$\alpha=p/N$','interpreter','latex');
xlabel('$\beta=q/N$','interpreter','latex');
axis([0 1 0 1]); axis equal tight ij; colorbar, c=caxis;
title('(a)  \quad measured asymptotic rate $\tilde\rho(\alpha,\beta)$','interpreter','latex');

axes('position',[ox+.5 oy+.5 sx sy]);
[C,h] = contourf(as,as,gt,v);
clabel(C,h,vl,'labelspacing',600);
ylabel('$\alpha=p/N$','interpreter','latex');
xlabel('$\beta=q/N$','interpreter','latex');
axis equal tight ij; caxis(c); colorbar;
title('(b) \quad  rate $\rho(\alpha,\beta)$ from Theorem 2','interpreter','latex');
colormap(goodbw)

v = [1:0.05:1.7]; vl = [1:.1:1.6];  % show contour lines, & to label, (c)
rat = g./gt; rat(gt==0)=nan;   % ratio, kill the edges
axes('position',[ox oy sx sy]);
[C,h] = contourf(as,as,rat,v); caxis([min(v),max(v)]);
clabel(C,h,vl,'labelspacing',600);
ylabel('$\alpha=p/N$','interpreter','latex');
xlabel('$\beta=q/N$','interpreter','latex');
axis equal tight ij; colorbar;
title('(c) \quad  ratio $\tilde\rho/\rho$, with Theorem 2','interpreter','latex');

end   % done with (0,1)^2 main plot

% CORNER DETAIL.....................................
amax = 4/(pi*exp(1));   % from Siegel-Chebyshev bnd a la ONeil-Rokhlin '07
astar = 0.117019018358476;  % crossover of rates
if 1    % regenerate data (fast)
%neff = 100; n = 15; as = (0:n)/neff;  % sample set in a and in b
neff = 50; n = 9; as = (0:n)/neff;  % sample set in a and in b
%neff = 100; n = 20; as = (0:n)/neff;  % sample set in a and in b
%neff = 50; n = 10; as = (0:n)/neff;  % sample set in a and in b
g = nan(numel(as)); gt=g; gc=g; st=g;
for a=2:n+1, alpha=as(a);
  for b=2:a, beta=as(b);
    gt(a,b) =  pi/2*(min(alpha,beta)-alpha*beta);  % KB theorem rate, Thm 2
    gc(a,b) = -log(max(alpha,beta)/amax)*min(alpha,beta); % Thm 3 (Siegel-Cheby)
    [g(a,b) st(a,b)] = dftmeasrate(neff,alpha*neff,beta*neff);  % new rate meas
    g(b,a) = g(a,b); gt(b,a) = gt(a,b); gc(b,a) = gc(a,b);
  end
end
g(isnan(g)) = 0.0;    % fill the zero-size results w/ zero rate
g(n+1,n+1) = (g(n,n)+dftmeasrate(5,1,1))/2;   % hack to approx (9/50,9/50) !
% figure; plot(as,diag(g),'+-');  % check diag curve is smooth
save ratescorner.mat n neff as g gt gc st
else, load ratescorner.mat; end

if 0   % draft of corner of fig as its own investigative fig. IIIIIIIIIIII
figure; subplot(2,2,1); imagesc(as,as,g);
v = [0:0.05:0.4]; vl = v;
[C,h] = contourf(as,as,g,v); clabel(C,h,vl,'labelspacing',300);
axis(n/neff*[0 1 0 1]); axis equal ij;
caxis([0 .4]); colorbar
subplot(2,2,2); imagesc(as,as,max(gt,gc));
[C,h] = contourf(as,as,max(gt,gc),v); clabel(C,h,vl,'labelspacing',300);
axis(n/neff*[0 1 0 1]); axis equal ij;
caxis([0 .4]); colorbar
hold on; plot(astar*[0 1 1],astar*[1 1 0],'--','color',.5*[1 1 1],'linewidth',2);
title('max gc,gt')
rat = g./max(gt,gc); 
subplot(2,2,3); imagesc(as,as,rat); v=[1.4:0.1:2]; vl=v;
[C,h] = contourf(as,as,rat,v); clabel(C,h,vl,'labelspacing',100);
axis(n/neff*[0 1 0 1]); axis equal ij;
colorbar; title('ratio')
hold on; plot(astar*[0 1 1],astar*[1 1 0],'--','color',.5*[1 1 1],'linewidth',2);
%subplot(2,2,4); imagesc(as,as,st); colorbar; title('g status');
colormap(goodbw)
end                                                    % IIIIIIIIIIIIII

axes('position',[ox+.5 oy+.2 .19 .19]);      % (d1)
v = [0:0.05:0.4]; vl = v;
[C,h] = contourf(as,as,g,v); clabel(C,h,vl,'labelspacing',300);
axis(n/neff*[0 1 0 1]); axis equal ij;
caxis([0 .4]); %colorbar
ylabel('$\alpha=p/N$','interpreter','latex');
%text(.03,-.02,'(d)  \quad  corner zooms of (a--c), also using Theorem $3$','interpreter','latex','fontsize',11);
text(-0.02,-.02,'corner zooms: $\;$ (d) $\;$ $\tilde\rho$','interpreter','latex','fontsize',11);
axes('position',[ox+.72 oy+.2 .19 .19]);      % (d2)
[C,h] = contourf(as,as,max(gt,gc),v); clabel(C,h,vl,'labelspacing',300);
axis(n/neff*[0 1 0 1]); axis equal ij;
caxis([0 .4]); %colorbar
hold on; plot(astar*[0 1 1],astar*[1 1 0],'--','color',0*[1 1 1],'linewidth',1);
xlabel('$\beta=q/N$','interpreter','latex');
text(-.03,astar,'$\alpha_\ast$','interpreter','latex','fontsize',12);
text(.02,-.02,'(e) $\; \rho$ (Thm.~2, 3)','interpreter','latex','fontsize',11);
rat = g./max(gt,gc);
axes('position',[ox+.515 oy-0.03 .23 .19]);   % (d3)
v=[1.4:0.1:2]; vl=v;
[C,h] = contourf(as,as,rat,v); clabel(C,h,vl,'labelspacing',100);
axis(n/neff*[0 1 0 1]); axis equal ij; colorbar
hold on; plot(astar*[0 1 1],astar*[1 1 0],'--','color',0*[1 1 1],'linewidth',1);
text(-.03,astar,'$\alpha_\ast$','interpreter','latex','fontsize',12);
text(0.26,0.1,'(f) $\;$  ratio $\tilde\rho/\rho$','interpreter','latex','fontsize',11);
text(0.29,0.125,'(Thm.~2, 3)','interpreter','latex','fontsize',11);
colormap(goodbw)

% output:
set(gcf,'paperposition',[0 0 10 8]); print -depsc2 rates.eps


if 0
axes('position',[ox+.5 oy sx sy]);
v = [1:0.05:1.7]; vl = [1:.1:1.6];
rat = g./gt; rat(gt==0)=nan;   % ratio, kill the edges
[C,h] = contourf(as,as,rat,v); caxis([min(v),max(v)]);
clabel(C,h,vl,'labelspacing',600);
ylabel('$\alpha=p/N$','interpreter','latex');
xlabel('$\beta=q/N$','interpreter','latex');
axis equal tight ij; colorbar;
end

% check (0,0) limit of ratio, for sharpness of Thm.3 in corner...
%figure; plot(as,diag(rat),'+-'); axis tight; v=axis; v([1 3])=[0 1.2]; axis(v);

