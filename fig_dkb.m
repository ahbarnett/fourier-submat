% make deplinthed KB fig for paper. Needs: sinc, arrow, vline
% Barnett 4/9/20

clear;

if 1 % ............. deplinthed Kaiser-Bessel pics:

  a = 0.7;   % alpha
q = 7;
s = (pi/2)*(1-a)*q;   % sigma
t = (q/2)*linspace(-1,1,1e3);    % t grid (all inside supp)
kb = real(besseli(0,s*sqrt(1-(2*t/q).^2)));
kb(abs(t)>q/2) = 0;
f = real(besseli(0,s*sqrt(1-(2*t/q).^2))) - 1;   % deplinthed (the "-1") KB
f(abs(t)>q/2) = 0;
om = 2*linspace(-1,1,1e3);   % freq grid
x = pi*q*om; fhat = q*( sinc(sqrt(x.^2 - s^2)) - sinc(x) );
om0 = 0.17;     % offset of example algebraic progression
oms = om0 + (-5:5);   % example summation grid "s"
xs = pi*q*oms; fhats = q*( sinc(sqrt(xs.^2 - s^2)) - sinc(xs) );

figure(1); clf; subplot(2,1,1);
plot(q/2*[-1 1 -1 1],[1 1 0 0],'g.','markersize',10); hold on;
h1 = plot(t,kb,'g-'); h2 = plot(t,f,'b-', 'linewidth',1);
flatbit = 1.5; plot((q/2)*[-flatbit 1; -1 flatbit],zeros(2,2),'b-','linewidth',2); % 2 segs
axis tight; xlabel('$t$','interpreter','latex'); ylabel('$f(t)$','interpreter','latex');
set(gca,'xtick',[-5 -q/2 0 q/2 5],'xticklabel',{'-5','-q/2','0','q/2','5'});
v = axis; v(4) = v(4)*1.03; axis(v);
h=legend([h1 h2], {'Kaiser-Bessel','deplinthed KB'});   % move the legend...
hp=get(h,'position'); hp(3:4)=hp(3:4)+[-.03,.03]; set(h,'position',hp);
arrow([0 0],[0 max(f)],'tipangle',10);
text(0.15,max(f)/2,'$I_0(\sigma)-1$','color',[0 0 0],'interpreter','latex');
text(-5,5.5,'(a)','fontsize',12);

subplot(2,1,2);
plot(om,fhat,'-'); hold on; axis tight;
xlabel('$\omega = x/\pi q$','interpreter','latex');
ylabel('$\hat f(\omega)$','interpreter','latex'); v = axis;
plot(oms,fhats,'k.','markersize',15);
v(3:4) = v(3:4) + v(4)*0.03*[-1 1]; axis(v);
vline(s/(pi*q)*[-1 1],'r--');
text(s/(pi*q)+0.05,max(fhat)*0.87,'$\pm \sigma/\pi q$','color',[1 0 0],'interpreter','latex');
for m=-2:1
  text(om0+m-0.15+0.2*(m==0), fhats(oms==om0+m)+max(fhat)/8, sprintf('$m=%d$',m),'color',[0 0 0],'interpreter','latex');
end
arrow([om0-2,0.4*max(fhat)], [om0-1,0.4*max(fhat)], 'ends','both','tipangle',10);
text(om0-1.8,0.5*max(fhat),'$a/\pi q = 1$','interpreter','latex');
text(-1.9,19,'(b)','fontsize',12);

set(gcf,'paperposition',[0 0 4.5 4.5]);
print -depsc2 dkbpair.eps
end

if 1 % ........... warped sinc pics:

s = 25;
x0 = .2*s^2;
lw = 0.1;
figure(2); clf;

x = linspace(0,x0,1e3);                % grid
si = sinc(x);            % funcs to show
wsi = (1./(x>=s)).*sinc(sqrt(x.^2-s^2));   % warped sinc
y0 = 2/s; v = [0 x0 -y0 y0];    % start panels
xsiz = 0.7; xsiz2 = 0.15;   % widths of main vs tail plots in [0,1] rel units
p = [.1 .7 xsiz .25]; axes('position',p);
plot(x,si, '-','linewidth',lw); axis(v); set(gca,'xticklabel',{});
tx = 45; ty = 0.06; fs = 12;
text(tx,ty,'(c) \quad sinc $x$','fontsize',fs,'interpreter','latex');
p(2) = p(2)-.3; axes('position',p);
plot(x,wsi, '-','linewidth',lw); axis(v); set(gca,'xticklabel',{});
vline(s,'r--'); text(0.8*s,0.03,'$\sigma$','color',[1 0 0],'fontsize',fs,'interpreter','latex');
text(tx,ty,'(d) \quad sinc $\sqrt{x^2-\sigma^2}$','fontsize',fs,'interpreter','latex');
p(2) = p(2)-.3; axes('position',p);
plot(x,(x>=s).*(wsi-si), '-','linewidth',lw); axis(v);
vline(s,'r--'); text(0.8*s,0.03,'$\sigma$','color',[1 0 0],'fontsize',fs,'interpreter','latex');
text(tx,ty,'(e)  sinc $\sqrt{x^2-\sigma^2}\; - $ sinc $x$','fontsize',fs,'interpreter','latex');
xlabel('$x$','interpreter','latex');

x = 10*s^2 + linspace(0,x0*(xsiz2/xsiz),1e3);   % tail panels, same abs x scale
si = sinc(x);
wsi = sinc(sqrt(x.^2-s^2));
y0 = 1.5/min(x); v = [min(x) max(x) -y0 y0];
p = [.15+xsiz .7 xsiz2 .25]; axes('position',p);
plot(x,si, '-','linewidth',lw); axis(v); set(gca,'xticklabel',{});
p(2) = p(2)-.3; axes('position',p);
plot(x,wsi, '-','linewidth',lw); axis(v); set(gca,'xticklabel',{});
p(2) = p(2)-.3; axes('position',p);
plot(x,wsi-si, '-','linewidth',lw); axis(v);
xlabel('$x$','interpreter','latex');

set(gcf,'paperposition',[0 0 4.5 4.5]);
print -depsc2 sincwarp.eps
end
