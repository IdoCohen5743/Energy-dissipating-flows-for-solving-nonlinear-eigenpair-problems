close all;
clear all;
clc;
rng(41001)

dt = 0.01;
dTheta = 4*pi/126;
N_std = 0.2;
theta = [-2*pi:dTheta:2*pi];
beta = 0.5;
alpha = 12*beta^2;
c = 4*beta^2;


lambda = 2;
u0 = alpha*(sech(beta*(theta)*sqrt(lambda))).^2;

u =  u0+0.1+N_std*randn(size(u0));

% u = alpha*(sech(beta*(theta-2*pi))).^2 + alpha*(sech(beta*(theta+2*pi)*5)).^2;
uini = u;
invdt2 = 1/(dTheta^2);
lapir = [1,-2,1,zeros(1,length(u)-3)];
lapic = [1,zeros(1,length(u)-3)];
lapi = toeplitz(lapic,lapir)*invdt2;
lapi = [[-1/dTheta^2,1/dTheta^2,zeros(1,length(u)-2)];...
        lapi;...
        [zeros(1,length(u)-2),1/dTheta^2,-1/dTheta^2]];
% lapi = [[-2/dTheta^2,1/dTheta^2,zeros(1,length(u)-2)];...
%         lapi;...
%         [zeros(1,length(u)-2),1/dTheta^2,-2/dTheta^2]];
    
lapi = -lapi';
initJ = -sum((u0*lapi).*u0);
[u] = zeroQ(c,u,dTheta,lapi);

iii = 0;
Ju = 0;


while true
    iii = iii+1;

    
    Tu =u*lapi;
	Qu = -c*u+(u.^2)/2;
    
    %% cohen-gilboa
    q = Qu;
    q = q/norm(q);
    p = Tu;
    p = p/norm(p);
	
    s =  sign(sum(q.*p));
    ut =  s*q  -p;
    
    Angle(iii) = acos(abs(Tu*Qu')/norm(Tu)/norm(Qu))*180/pi;

    lambda = (Tu*u')/(Qu*u');
%     lambda = sign(sum(Tu.*Qu))*norm(Tu)/norm(Qu);
    Ju(end+1) = 0.5*(Tu*u');
    dJu(iii) = Ju(end) - Ju(end-1);
    
%     deltaU = dt*ut;


    e = -Tu+lambda*Qu;
    normE(iii) = norm(e);
    deltaU = 0.5*(dTheta^2)*e/2;
    if(~isnan(deltaU))
    else
        oooo = 0;
    end
    u = u + deltaU;  
    
    if(~isnan(u))
    else
        oooo = 0;
    end
    normDu(iii) = (norm(deltaU));
    [u] = zeroQ(c,u,dTheta,lapi);
    if(~isnan(u))
    else
        oooo = 0;
    end
   
    
    if(Angle(end)>179.9999||Angle(end)<0.0001)
        break;
    end
    if(mod(iii,1000)==0)
    figure(26);

        subplot(3,3,1);plot(Qu);title('Qu');
        subplot(3,3,2);plot(p);title('p');
        subplot(3,3,3);plot(u0);hold on;plot(uini);hold on;plot(u);hold off;title('u evolution');
        
        subplot(3,3,4);plot(normE);title('norm e');
        legend(['e = ',num2str(normE(end)),' th = 5.0000e-05'])
        subplot(3,3,5);plot(Angle);title('Angle');
        legend(['Angle = ',num2str(Angle(end))]);
        subplot(3,3,7);semilogy(normDu);title('normDu');
        subplot(3,3,8);plot(deltaU);title('deltaU');
%         legend(['meanQ = ',num2str(meanQ(end))])
        subplot(3,3,9);semilogy(abs(dJu));
        drawnow;
    end
    
   
end
iii

         figure(26);

        subplot(3,3,1);plot(Qu);title('Qu');
        subplot(3,3,2);plot(p);title('p');
        subplot(3,3,3);plot(u0);hold on;plot(uini);hold on;plot(u);hold off;title('u evolution');
        
        subplot(3,3,4);plot(normE);title('norm e');
        legend(['e = ',num2str(normE(end)),' th = 5.0000e-05'])
        subplot(3,3,5);plot(Angle);title('Angle');
        legend(['Angle = ',num2str(Angle(end))]);
        subplot(3,3,7);semilogy(normDu);title('normDu');
        subplot(3,3,8);plot(deltaU);title('deltaU');
%         subplot(3,3,8);plot(meanQ);title('meanQ');
%         legend(['meanQ = ',num2str(meanQ(end))])
        %subplot(3,3,9);plot(lambda);title('\lambda');
        subplot(3,3,9);semilogy(abs(dJu));
        drawnow;
figure()
plot(Angle)

h = figure();plot(theta,u0);hold on;plot(theta,uini);hold on;plot(theta,u,'k');hold off;title('Result');
h.Children.Title.Interpreter = 'Latex';
h.Children.Title.FontSize = 40;

h.Children.XTick = [-7:2:7]*pi/4;
h.Children.XTickLabel = [{'$-\frac{7\pi}{4}$'},{'$-\frac{5\pi}{4}$'},{'$-\frac{3\pi}{4}$'},...
    {'$-\frac{\pi}{4}$'},{'$\frac{\pi}{4}$'},{'$\frac{3\pi}{4}$'},{'$\frac{5\pi}{4}$'},{'$\frac{7\pi}{4}$'}];
h.Children.FontSize = 38;
h.Children.TickLabelInterpreter = 'Latex';
h.Children.XLim = [-2*pi,2*pi];
h.Children.YLim = [-1,4];

h.Children.Children(1).LineWidth = 4;
h.Children.Children(2).LineWidth = 4;

h.Children.Children(3).LineWidth = 4;

legend('$u^*$','Noisy','Ours')
h.Children(1).Interpreter = 'Latex';

g=figure();
plot(theta,Qu/norm(Qu)*s,'r','LineWidth',4);hold on;plot(theta,p/norm(p),'b:','LineWidth',4);hold off;
title(['$s\frac{Q\left(u^*\right)}{||Q\left(u^*\right)||}$ $-\frac{\Delta u^*}{||\Delta u^*||}$'])


g.Children.Title.Interpreter = 'Latex';
g.Children.Title.FontSize = 40;

g.Children.XTick = [-7:2:7]*pi/4;
g.Children.XTickLabel = [{'$-\frac{7\pi}{4}$'},{'$-\frac{5\pi}{4}$'},{'$-\frac{3\pi}{4}$'},...
    {'$-\frac{\pi}{4}$'},{'$\frac{\pi}{4}$'},{'$\frac{3\pi}{4}$'},{'$\frac{5\pi}{4}$'},{'$\frac{7\pi}{4}$'}];
g.Children.FontSize = 38;
g.Children.TickLabelInterpreter = 'Latex';

g.Children.XLim = [-2*pi,2*pi];
g.Children.YLim = [-0.15,0.3];

g.Children.Children(1).LineWidth = 4;
g.Children.Children(2).LineWidth = 4;
legend('$s\frac{Q\left(u^*\right)}{||Q\left(u^*\right)||}$','$-\frac{\Delta u^*}{||\Delta u^*||}$')
g.Children(1).Interpreter = 'Latex';

k = figure();
plot(theta,Tu./Qu);hold on;plot(theta,lambda*ones(size(u0)),'--');hold off; title('-$\Delta u^*$ and $Q\left(u^*\right)$ Quotient')
k.Children.Title.Interpreter = 'Latex';
% k.Children.Title.FontSize = 40;

k.Children.XLim = [-2*pi,2*pi];
k.Children.XTick = [-7:2:7]*pi/4;
k.Children.XTickLabel = [{'$-\frac{7\pi}{4}$'},{'$-\frac{5\pi}{4}$'},{'$-\frac{3\pi}{4}$'},...
    {'$-\frac{\pi}{4}$'},{'$\frac{\pi}{4}$'},{'$\frac{3\pi}{4}$'},{'$\frac{5\pi}{4}$'},{'$\frac{7\pi}{4}$'}];
k.Children.FontSize = 60;
k.Children.TickLabelInterpreter = 'Latex';

k.Children.YTick = [1.5:0.02:1.6];
k.Children.YLim = [1.5,1.6];
k.Children.Children(1).LineWidth = 4;
k.Children.Children(2).LineWidth = 4;
legend('$\frac{-\Delta u^*}{Q\left(u^*\right)}$','$\lambda^*$')
k.Children(1).Interpreter = 'Latex';
k.Children(1).FontSize = 60;
grid on



g=figure();
plot(theta,Qu/norm(Qu)*s,'k','LineWidth',4);hold on;plot(theta,p/norm(p),'--y','LineWidth',4);hold off;
%title(['$s\frac{Q\left(u^*\right)}{||Q\left(u^*\right)||}$, $-\frac{\Delta u^*}{||\Delta u^*||}$'])


g.Children.Title.Interpreter = 'Latex';


g.Children.XTick = [-7:2:7]*pi/4;
g.Children.XTickLabel = [{'$-\frac{7\pi}{4}$'},{'$-\frac{5\pi}{4}$'},{'$-\frac{3\pi}{4}$'},...
    {'$-\frac{\pi}{4}$'},{'$\frac{\pi}{4}$'},{'$\frac{3\pi}{4}$'},{'$\frac{5\pi}{4}$'},{'$\frac{7\pi}{4}$'}];
g.Children.FontSize = 60;
g.Children.TickLabelInterpreter = 'Latex';
% g.Children.Title.FontSize = 60;
g.Children.XLim = [-2*pi,2*pi];
g.Children.YLim = [-0.15,0.3];

g.Children.Children(1).LineWidth = 4;
g.Children.Children(2).LineWidth = 4;
legend('$s\frac{Q\left(u^*\right)}{||Q\left(u^*\right)||}$','$-\frac{\Delta u^*}{||\Delta u^*||}$')
g.Children(1).Interpreter = 'Latex';
grid on;

k = axes('Position',[.7 .7 .2 .2]);
box on
plot(theta,Tu./Qu);hold on;plot(theta,lambda*ones(size(u0)),'--');hold off; title('-$\Delta u^*$ and $Q(u^*)$ Quotient')
k.Title.Interpreter = 'Latex';
% k.Children.Title.FontSize = 40;

k.XLim = [-2*pi,2*pi];
k.XTick = [-2:1:2]*pi;
% k.XTickLabel = [{'$-\frac{7\pi}{4}$'},{'$-\frac{5\pi}{4}$'},{'$-\frac{3\pi}{4}$'},...
%     {'$-\frac{\pi}{4}$'},{'$\frac{\pi}{4}$'},{'$\frac{3\pi}{4}$'},{'$\frac{5\pi}{4}$'},{'$\frac{7\pi}{4}$'}];
k.XTickLabel = [{'$-2\pi$'},{'$-\pi$'},{'$0$'},{'$\pi$'},{'$2\pi$'}];
k.FontSize = 30;
k.TickLabelInterpreter = 'Latex';

k.YTick = [1.5:0.05:1.6];
k.YLim = [1.5,1.6];
k.Children(1).LineWidth = 4;
k.Children(2).LineWidth = 4;
kl = legend('$-\frac{\Delta u^*}{Q(u^*)}$','$\lambda^*$');
kl.Interpreter = 'latex';
% k.Children(1).Interpreter = 'Latex';
% k.Children(1).FontSize = 60;
grid on
