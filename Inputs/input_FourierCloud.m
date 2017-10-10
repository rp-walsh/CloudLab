%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% - Input file for WaveTank.m using periodic BCs in x and homogenous
%   Neumann in y for the pressure inversion.
% - This input file contains ONLY the problem set up and NOT the exact
%   solution i.e. this input file is IS NOT well suited for a convergence
%   study. 
%
% This input file MUST specify atleast the following:
%                 - uInit(x,z) = initial data for u
%                 - wInit(x,z) = initial data for w
%                 - sInit(x,z) = initial data for s 
%                 - rho(x,z,s) = density as a function of position and
%                                entropy 
%                 - PBCT(x,t) = Top Neumann BC for pressure
%                 - PBCB(x,t) = Top Neumann BC for pressure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define solution parameters
g = 9.81;
bb = 1; gg = 1; % Entropy/total mixing ratio coeff 
c1p = 2; c1m = 1; % above/below interface entropy
c2p = 2; c2m = 1; % above/below interface total mixing ratio
H = 1; h = 1/2; % height of the domain/interface
n = 2; k = n*pi; % fourier/horiz number (n = even)
a1 = 1; a2 = 1; % liquid water mixing ratio coeffs
sb_z = 1; rTb_z = 1; % Background entropy/total mixing ratio coeff
cp = c1p*sb_z + c2p*rTb_z; % Effective cp
cm = c1m*sb_z + c2m*rTb_z; % Effective cm

%% Define verticle wave strcture
alphap = @(k,omega) k.*sqrt(g*cp./omega.^2 - 1);
alpham = @(k,omega) k.*sqrt(g*cm./omega.^2 - 1);

%% Define dispersion relation
F = @(k,w) alpham(k,w).*sin(alpham(k,w)*h).*cos(alphap(k,w)*(h-H)).*...
    (g*(bb*(c1p - c1m)+gg*(c2p - c2m))./(a1*bb+a2*gg)*(a1*sb_z + a2*rTb_z)+1) -...
    alphap(k,w).*sin(alphap(k,w)*(h-H)).*cos(alpham(k,w)*h);

omega = fzero(@(w) F(k,w),sqrt(g/2));

%% Value determining solution strcture (>0 = waves)
%SolnStructP = g*cp/omega^2
%SolnStructM = g*cm/omega^2

%% Define functions
alphap = alphap(k,omega);
alpham = alpham(k,omega);
%eps = 0.1;
%eps = .01;
eps = .05;
Ap = eps*cos(alpham*h); Am = eps*cos(alphap*(h-H));
sHatM = @(x,z,t ) -sb_z*alpham*Am/(cm*g-omega^2).*sin(alpham*z).*sin(k*x-omega*t);
rTHatM = @(x,z,t) -rTb_z*alpham*Am/(cm*g-omega^2).*sin(alpham*z).*sin(k*x-omega*t);
Gmm = @(x,t) h - (a1*sHatM(x,h,t) + a2*rTHatM(x,h,t))./(a1*bb + a2^gg);

Pm = @(x,z,t) g*(c1m*bb + c2m*gg)*(z-h).^2./2 + Am*cos(alpham*z).*sin(k*x-omega*t); 
Pp = @(x,z,t) g*(c1p*bb + c2p*gg)*(z-h).^2./2 + Ap*cos(alphap*(z-H)).*sin(k*x-omega*t); 
PExact =@(x,z,t) Pm(x,z,t).*((Gmm(x,t)-z)>=0) + Pp(x,z,t).*((z-Gmm(x,t))>0);
%PBkg = @(x,z,t) bbeta*g*cm*(z-h).^2./2.*((Gmm(x,t)-z)>=0) + bbeta*g*cp*(z-h).^2./2.*((z-Gmm(x,t))>0);
PBkg = @(x,z,rl) g*(c1m*bb+c2m*gg)*(z-h).^2./2.*(rl<0) + g*(c1p*bb+c2p*gg)*(z-h).^2./2.*(rl>=0);

%Pmx = @(x,z,t) k*Am.*cos(alpham*z).*cos(k*x-omega*t);
%Ppx = @(x,z,t) k*Ap.*cos(alphap*(z-H)).*cos(k*x-omega*t);
%PxExact = @(x,z,t) Pmx(x,z,t).*((Gmm(x,t)-z)>=0) + Ppx(x,z,t).*((z-Gmm(x,t))>0);

%Pmz = @(x,z,t) bbeta*g*cm*(z-h) - alpham*Am*sin(alpham*z).*sin(k*x-omega*t); 
%Ppz = @(x,z,t) bbeta*g*cp*(z-h) - alphap*Ap*sin(alphap*(z-H)).*sin(k*x-omega*t); 
%PzExact = @(x,z,t) Pmz(x,z,t).*((Gmm(x,t)-z)>=0) + Ppz(x,z,t).*((z-Gmm(x,t))>0);

sBkg = @(x,z) bb*(z-h);
Sm = @(x,z,t) - sb_z*alpham*Am/(cm*g-omega^2).*sin(alpham*z).*sin(k*x-omega*t);
Sp = @(x,z,t) - sb_z*alphap*Ap/(cp*g-omega^2).*sin(alphap*(z-H)).*sin(k*x-omega*t); 
SExact =@(x,z,t) Sm(x,z,t).*((Gmm(x,t)-z)>=0) + Sp(x,z,t).*((z-Gmm(x,t))>0);
sInit = @(x,z) SExact(x,z,0);

rTBkg = @(x,z) gg*(z-h);
rTm = @(x,z,t) - rTb_z*alpham*Am/(cm*g-omega^2).*sin(alpham*z).*sin(k*x-omega*t);
rTp = @(x,z,t) - rTb_z*alphap*Ap/(cp*g-omega^2).*sin(alphap*(z-H)).*sin(k*x-omega*t); 
rTExact =@(x,z,t) rTm(x,z,t).*((Gmm(x,t)-z)>=0) + rTp(x,z,t).*((z-Gmm(x,t))>0);
rTInit = @(x,z) rTExact(x,z,0);

rhoExact =@(x,z,t) (-c1m*Sm(x,z,t)-c2m*rTm(x,z,t)).*((Gmm(x,t)-z)>=0) +...
                   (-c1p*Sp(x,z,t)-c2p*rTp(x,z,t)).*((z-Gmm(x,t))>0);
rho =@(s,rT) (-c1p*s-c2p*rT).*((a1*s + a2*rT)>=0) + (-c1m*s-c2m*rT).*((a1*s + a2*rT)<0);

%drhodz2 = @(x,z,t) -cm*bbeta + cm*alpham.^2*Am/(cm*g-omega^2).*cos(alpham*z).*sin(k*x-omega*t);%minus
%drhodz1 = @(x,z,t) -cp*bbeta + cp*alphap.^2*Ap/(cp*g-omega^2).*cos(alphap*(z-H)).*sin(k*x-omega*t);%plus

Wm = @(x,z,t) -omega/(g*cm-omega^2)*Am*alpham.*sin(alpham*z).*cos(k*x-omega*t);
Wp = @(x,z,t) -omega/(g*cp-omega^2)*Ap*alphap.*sin(alphap*(z-H)).*cos(k*x-omega*t);
wExact =@(x,z,t) Wm(x,z,t).*((Gmm(x,t)-z)>=0) + Wp(x,z,t).*((z-Gmm(x,t))>0);
wInit = @(x,z) wExact(x,z,0);

Um = @(x,z,t) k/omega*Am*cos(alpham*z).*sin(k*x-omega*t); 
Up = @(x,z,t) k/omega*Ap*cos(alphap*(z-H)).*sin(k*x-omega*t); 
uExact =@(x,z,t) Um(x,z,t).*((Gmm(x,t)-z)>=0) + Up(x,z,t).*((z-Gmm(x,t))>0);
uInit = @(x,z) uExact(x,z,0);

%PBCT = @(x,t) (H-h)*bbeta*g*cp*ones(size(x)); 
%PBCB = @(x,t) -h*bbeta*g*cm*ones(size(x)); 

PBCT = @(x,t) zeros(size(x));
PBCB = @(x,t) zeros(size(x));

wBkg = @(x,z) zeros(size(x));
uBkg = @(x,z) zeros(size(x));

PBkg = @(x,z,rl) g*(c1m*bb + c2m*gg)*(z-h).^2./2.*(rl<0) +...
       g*(c1p*bb + c2p*gg)*(z-h).^2./2.*(rl>=0);
PBkgBCT = @(x,t) (H-h)*g*(c1p*bb + c2p*gg)*ones(size(x)); 
PBkgBCB = @(x,t) -h*g*(c1m*bb + c2m*gg)*ones(size(x));