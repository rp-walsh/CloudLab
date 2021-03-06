clear all; close all;

addpath(genpath('..'));
input_FourierCloud;
L = 1; H = 1;
x = linspace(0,1,1000);
[x,z] = meshgrid(x);

time = 0;

wtmp = linspace(0,5,1000);
ktmp = linspace(0,6*pi,1000);
figure(1)
[Ktmp,Wtmp] = meshgrid(ktmp,wtmp);
contour(Ktmp,Wtmp,F(Ktmp,Wtmp),[0 0])
xlabel('k')
ylabel('omega')

for j = time
    fig = figure(3);
    
    subplot(2,2,1)
    pcolor(x,z,uExact(x,z,j))
    title(['Solution u at time = ' num2str(j)])
    shading flat; axis equal; colorbar; axis([0 L 0 H]);
    
    subplot(2,2,2)
    surf(x,z,wExact(x,z,j))
    view([0,90])
    title(['Solution w at time = ' num2str(j)])
    shading flat; axis equal; colorbar; axis([0 L 0 H]);
    hold on
    plot(x(1,:),Gmm(x(1,:),j),'r-');
    hold off
    axis square

    subplot(2,2,3)
    pcolor(x,z,rhoExact(x,z,j))
    title(['Solution rho at time = ' num2str(j)])
    shading flat; axis equal; colorbar; axis([0 L 0 H]);
    
    subplot(2,2,4)
    pcolor(x,z,PExact(x,z,j))
    title(['Solution P at time = ' num2str(j)])
    shading flat; axis equal; colorbar; axis([0 L 0 H]);
    drawnow
    pause(0.1)

    W = PDist(x,z,j);
    figure(4)
    plot(z(:,225),W(:,225))
    title('P')
    hold on
    plot(Gmm(x(1,225),j),0,'ro')

    W = rhoExact(x,z,j);
    figure(5)
    plot(z(:,225),W(:,225))
    title('rho')
    hold on
    plot(Gmm(x(1,225),j),0,'ro')

    W = wExact(x,z,j);
    figure(6)
    plot(z(:,225),W(:,225))
    title('w')
    hold on
    plot(Gmm(x(1,225),j),0,'ro')

    W = uExact(x,z,j);
    figure(7)
    plot(z(:,225),W(:,225))
    title('u')
    hold on
    plot(Gmm(x(1,225),j),0,'ro')

    W = PzExact(x,z,j);
    figure(8)
    plot(z(:,225),W(:,225))
    title('P_z')
    hold on
    plot(Gmm(x(1,225),j),0,'ro')
end
