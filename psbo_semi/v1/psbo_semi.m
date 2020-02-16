function psbo_semi
close all


q=18;               % HO frequency in Hz
% q=0;
omega=2*pi*q;      % HO angular frequency Hz*rad

TB0=17E-3;          % bloch oscillation period in s
fB0=1/TB0;           % bloch oscillation frequency in Hz
Tf=200E-3;          % final time
U0=4;               % lattice depth
phi0=0;
U1=1;

% Define physical constants
amu=1.66E-27;       % 1 AMU
m=7*amu;            % Lithium mass
lambda=1064E-9;     % Wavelength of light
h=6.626E-34;        % Planck's Const.
kL=2*pi/lambda;     % Wave Vector
hb=h/(2*pi);        % Reduced planck's constant
c=3E8;              % Speed of light
Er=hb^2*kL^2/(2*m); % Recoil Energy
vR=hb*kL/m;         % recoil velocity
d=lambda/2;         % Lattice site distance
fR=Er/h;            % recoil frequency in Hz

% Simulation parameters
F0=h*fB0/d;                               % Initial force in N
frequencyRamp = 0;  %yes or no
dfdt = 151;
%% Calculate and plot dispersions

% Load bandwidth figure
hF1=open('../bands.fig');

% Get the axes for the bandwidth
ax=hF1.Children;
axes(ax(2))
yyaxis left
pS=ax(2).Children(2);

% Get data for bandwidth
X=pS.XData;
Y=pS.YData;
close(hF1);

% Replot the bandwidth
hF2=figure(2);
hF2.Color='w';
hF2.Position(3:4)=[400 300];
ax=axes;
set(ax,'box','on','linewidth',1,'fontsize',14,'fontname','times');
hold on
plot(X,Y,'k-','linewidth',2);
xlabel('lattice depth ($E_R$)','interpreter','latex');
ylabel('ground bandwidth ($E_R$)','interpreter','latex');

% Construct the bandwidth function
BW=@(U) interp1(X,Y,U);

% Construct dEdK
Ek=@(U,k) -cos(pi*k)*0.5*BW(U);
dEdK=@(U,K) sin(pi*K)*0.5*BW(U)*pi;

% Check the band structure
K=linspace(-1,1,500);
hF3=figure(3);
hF3.Color='w';
hF3.Position(3:4)=[400 300];
axes;
set(gca,'box','on','linewidth',1,'fontsize',14,'fontname','times');
xlabel('quasimomentum ($\hbar k_L$)','interpreter','latex');
hold on

yyaxis left
ylabel('energy ($E_R$)','interpreter','latex');
plot(K,Ek(U0,K),'-','linewidth',2);

yyaxis right
ylabel('$\partial E/\partial k~(E_R/\hbar k_L)$','interpreter','latex');
plot(K,dEdK(U0,K),'-','linewidth',2);

text(.98,.02,['$U_0=' num2str(U0) 'E_R$'],'interpreter','latex',...
    'fontsize',14,'units','normalized','horizontalalignment','right',...
    'verticalalignment','bottom');

L=(fR*BW(U0)/fB0);

%% Time Evolve Static Lattice
% Time evolvs a classical particle in a harmonic potential
% with a custom dispersion.

% momentum is in recoil momentum in hbar*k_L = m*vR
% position is in lattice spacing d

% Phase space time evolution operator d(x,p)/dt
% dx/dt=dE/dp=dE/dk=(E_r/hbar*k_L)*dE/dK=vR/2*dE/dK
% dp/dt=-dU/dx=F-m*omega^2*x

% Then parametrize x and p in terms of d and hbar*k_L
% dX/dt = (1/d)(vR/2*dE/dK)
% dP/dt = (F-m*omega^2*x)/(hbar*k_L) = (F-m*omega^2*X*d)/(m*vR)
% P = [x/d; p/(m*vR)]

% Time evolution operator
dPdT=@(t,P) [....
    (1/d)*(vR/2)*dEdK(U0,P(2));
    (F0-m*omega^2*P(1)*d)/(m*vR)];

% time evolution that sucker
Tlim=[0; Tf];   % time limits
Pini=[0; 0];    % initial phase space vector (x=0,k=0);
opt=odeset('MaxStep',1E-4); % max step size in seconds
[T,Y]=ode45(@(t,V) dPdT(t,V),Tlim,Pini, opt);


hF4=figure(4);
hF4.Color='w';
hF4.Position(3:4)=[600 400];
subplot(211);
set(gca,'box','on','linewidth',1,'fontsize',14,...
    'fontname','times');
xlabel('time (ms)','interpreter','latex');
yyaxis left
plot(T*1e3,Y(:,1),'-','linewidth',2);
ylabel('position ($d$)','interpreter','latex');

yyaxis right
plot(T*1e3,mod(Y(:,2)+1,2)-1,'-','linewidth',.5);
ylabel('quasimomentum ($\hbar k_L$)','interpreter','latex');

str=['$f_B=' num2str(round(fB0,1)) '~\mathrm{Hz}$;~'  ...
    '$T_B=' num2str(round(1E3*TB0,1)) '~\mathrm{ms}$;~'  ...
    '$U_0=' num2str(round(U0,1)) 'E_R$;~'  ...
    '$\Delta E=' num2str(round(BW(U0),2)) 'E_R$;~'  ...
    '$2L_{\mathrm{WS}}=' num2str(round(L)) 'd$;~' ...
    '$q=' num2str(round(q,1)) '~\mathrm{Hz}$'];

ax=gca;
ax.Units='pixels';
ax.Position(4)=ax.Position(4)-20;


text(0,ax.Position(4)+10,str,'units','pixels','interpreter','latex',...
    'fontsize',8,...
    'verticalalignment','bottom',...
    'horizontalalignment','left');

subplot(212)
plot([0 max(T)]*1e3,[1 1]*U0,'k-','linewidth',1);
set(gca,'box','on','linewidth',1,'fontsize',14,...
    'fontname','times');
xlabel('time (ms)','interpreter','latex');

ylabel('lattice depth ($E_R$)','interpreter','latex');


%% Time evolve with modulated lattice depth
% Time evolvs a classical particle in a harmonic potential
% with a custom dispersion.

% momentum is in recoil momentum in hbar*k_L = m*vR
% position is in lattice spacing d

% Phase space time evolution operator d(x,p)/dt
% dx/dt=dE/dp=dE/dk=(E_r/hbar*k_L)*dE/dK=vR/2*dE/dK
% dp/dt=-dU/dx=F-m*omega^2*x

% Then parametrize x and p in terms of d and hbar*k_L
% dX/dt = (1/d)(vR/2*dE/dK)
% dP/dt = (F-m*omega^2*x)/(hbar*k_L) = (F-m*omega^2*X*d)/(m*vR)
% P = [x/d; p/(m*vR)]

% Time evolution operator

if frequencyRamp == 0
    Ufunc=@(t) U0+U1*sin(2*pi*fB0*t+phi0);  %no frequency ramp
else
    Ufunc=@(t) U0+U1*sin(2*pi*(fB0+dfdt*t).*t+phi0);  %frequency ramp
end

dPdT=@(t,P) [....
    (1/d)*(vR/2)*dEdK(Ufunc(t),P(2));
    (F0-m*omega^2*P(1)*d)/(m*vR)];

% time evolution that sucker
Tlim=[0; Tf];   % time limits
Pini=[0; 0];    % initial phase space vector (x=0,k=0);
opt=odeset('MaxStep',1E-4); % max step size in seconds
[T,Y]=ode45(@(t,V) dPdT(t,V),Tlim,Pini, opt);


hF5=figure(5);
hF5.Color='w';
hF5.Position(3:4)=[600 400];


subplot(211);
set(gca,'box','on','linewidth',1,'fontsize',14,...
    'fontname','times');
xlabel('time (ms)','interpreter','latex');
yyaxis left
plot(T*1e3,Y(:,1),'-','linewidth',2);
ylabel('position ($d$)','interpreter','latex');

yyaxis right
plot(T*1e3,mod(Y(:,2)+1,2)-1,'-','linewidth',.5);
ylabel('quasimomentum ($\hbar k_L$)','interpreter','latex');

str=['$f_B=' num2str(round(fB0,1)) '~\mathrm{Hz}$;~'  ...
    '$T_B=' num2str(round(1E3*TB0,1)) '~\mathrm{ms}$;~'  ...
    '$U_0=' num2str(round(U0,1)) 'E_R$;~'  ...
    '$\Delta E=' num2str(round(BW(U0),2)) 'E_R$;~'  ...
    '$U_1=' num2str(round(U1,1)) 'E_R$;~' ...
    '$\phi=' num2str(round(phi0/pi,1)) '\pi$;~' ...
    '$2L_{\mathrm{WS}}=' num2str(round(L)) 'd$;~' ...
    '$q=' num2str(round(q,1)) '~\mathrm{Hz}$'];

ax=gca;
ax.Units='pixels';
ax.Position(4)=ax.Position(4)-20;



text(0,ax.Position(4)+10,str,'units','pixels','interpreter','latex',...
    'fontsize',8,...
    'verticalalignment','bottom',...
    'horizontalalignment','left');

subplot(212)

plot(T*1e3,Ufunc(T),'k-','linewidth',1);
set(gca,'box','on','linewidth',1,'fontsize',14,...
    'fontname','times');
xlabel('time (ms)','interpreter','latex');

ylabel('lattice depth ($E_R$)','interpreter','latex');


end

