clear all; close all;

k = linspace(0,6*pi,1000);
c = 1;
g = 9.81;
figure(1)
hold on;
omega = @(k,l) sqrt(c*g*k.^2./(k.^2 + l.^2));
for n = 0:100
    l = n*pi;
    plot(k,omega(k,l))
end
xlabel('k')
ylabel('w')

bb = 1; % Entropy coeff
gg = 1; % total mixing ratio coeff 
c1p = 1.1; % above interface entropy
c1m = 0.9; % below interface entropy
c2p = 1.1; % above interface total mixing ratio
c2m = 0.9; % below interface total mixing ratio
a1 = 1; a2 = 1; % liquid water mixing ratio coeffs
sb_z = 1;
rTb_z = 1;
cp = c1p*sb_z + c2p*rTb_z;
cm = c1m*sb_z + c2m*rTb_z;
h = 0.5;
H = 1;
alphap = @(k,omega) k.*sqrt(g*cp./omega.^2 - 1);
alpham = @(k,omega) k.*sqrt(g*cm./omega.^2 - 1);

F = @(k,w) alpham(k,w).*sin(alpham(k,w)*h).*cos(alphap(k,w)*(h-H)).*...
    (g*(bb*(c1p - c1m)+gg*(c2p - c2m))./(a1*bb+a2*gg)*(a1*sb_z + a2*rTb_z)+1) -...
    alphap(k,w).*sin(alphap(k,w)*(h-H)).*cos(alpham(k,w)*h);

%F = @(k,w) alpham(k,w).*sin(alpham(k,w)*h).*cos(alphap(k,w)*(h-H)).*(g*(cp-cm)./(g*cm-w.^2)+1) ...
%        - alphap(k,w).*sin(alphap(k,w)*(h-H)).*cos(alpham(k,w)*h); 

%wtmp = linspace(sqrt(g/2)-1,sqrt(g/2)+1,1000);
wtmp = linspace(0,5,1000);
%wtmp = linspace(0,3.5,1000);
ktmp = linspace(0,6*pi,1000);
figure(2)
[Ktmp,Wtmp] = meshgrid(ktmp,wtmp);
contour(Ktmp,Wtmp,F(Ktmp,Wtmp),[0 0])
xlabel('k')
ylabel('omega')

%% Finding particular omegas
k0 = 2*pi;
omegaAvg = omega(k0,pi);
figure(1)
hold on
plot(k0,omegaAvg,'x')

figure(3)
wtmp = linspace(sqrt(g/2)-1,sqrt(g/2)+1,1000);
plot(wtmp,F(k0,wtmp));
ylabel('F(k_0,w)')
xlabel('w')

omegaSplt = fzero(@(w) F(k0,w),omegaAvg)

figure(3)
hold on
plot(omegaSplt,0,'x')

figure(2)
hold on
plot(k0,omegaSplt,'x')

figure(4)
surf(Ktmp,Wtmp,F(Ktmp,Wtmp),'Edgecolor','none')
xlabel('k')
ylabel('w')