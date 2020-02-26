close all; clear;
% Semiclassical simulation of noninteracting particles in an amplitude
% modulated optical cosine lattice with a weak overall confining harmonic 
% potential. Computes trajectories of initial conditions spanning a finite
% region of phase space for a minimum uncertainty wavepacket using
% semiclassical equations of motions. Trajectories of different initial
% conditions are weighted with a normal distribution to calculate first and
% second moments of the cloud's position over time.

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

%% Initialize Parameters HERE
q=15.5;              % HO frequency in Hz
omega=2*pi*q;      % HO angular frequency Hz*rad

TB0=18E-3;      % bloch oscillation period in s
fB0=1/TB0;         % bloch oscillation frequency in Hz
fD = 53.56;        % set drive frequency
Tf=200E-3;         % final time
U0=4.345;            % lattice depth in Er
phi0=pi/2;       % phase of sine wave drive
U1=1.02;            % modulation half amplitude in Er
F0=h*fB0/d;        % Initial force in N

% Simulation parameters
frequencyRamp = 0;  %yes or no
modulation = 1; %0 is static lattice, 1 is modulated lattice
dfdt = 115;
sizex = 51;      % number of points in position space
sizep = 51;      % number of points in momentum space
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

%% Compensate for initial half BO for experimental comparison

J = Er*BW(U0+U1*sin(phi0))/4;   %tunneling energy in Joules
lws = 2*J/F0;                   %localization length in meters
F = F0 - m*omega^2*lws;         %actual force at start of the drive
%% Time evolve with modulated lattice depth
% Time evolves a classical ensemble of particles in a harmonic potential
% with a custom dispersion.

% momentum is in recoil momentum in hbar*k_L = m*vR
% position is in lattice spacing d

% math for differential equations
% Phase space time evolution operator d(x,p)/dt
% dx/dt=dE/dp=dE/dk=(E_r/hbar*k_L)*dE/dK=vR/2*dE/dK
% dp/dt=-dU/dx=F-m*omega^2*x

% Then parametrize x and p in terms of d and hbar*k_L
% dX/dt = (1/d)(vR/2*dE/dK)
% dP/dt = (F-m*omega^2*x)/(hbar*k_L) = (F-m*omega^2*X*d)/(m*vR)
% P = [x/d; p/(m*vR)]


% Time evolution operator
% Ufunc is the modulated lattice depth
if frequencyRamp == 0
    Ufunc=@(t) U0+U1*sin(2*pi*fD*t+phi0);  %no frequency ramp
else
    Ufunc=@(t) U0+U1*sin(2*pi*(fD+dfdt*t).*t+phi0);  %frequency ramp
end
%differential equations for advancing position and momentum
if modulation == 1
    dPdT=@(t,P) [....
        (1/d)*(vR/2)*dEdK(Ufunc(t),P(2));
        (F-m*omega^2*P(1)*d)/(m*vR)];
else
    dPdT=@(t,P) [....
        (1/d)*(vR/2)*dEdK(U0,P(2));
        (F-m*omega^2*P(1)*d)/(m*vR)];

end

% Classically time evolve a gaussian distribution of initial phase space conditions
% centered at x=0 and p=p0 when the modulation starts. A box in phase
% space is integrated over and a probability weight is associated with each
% initial condition. The solutions are interpolated for identical time vectors
% so that the moments can be computed.
xsigma = 60; %position width
psigma = pi/(2*xsigma); %minimum uncertainty
p0 = 1; %if starting modulation not at the center of the Brillouin zone
xwidth = 3*xsigma; pwidth = 3*psigma; %integrate phase space to 3 sigma point
xinit = linspace(-xwidth,xwidth,sizex);
pinit = linspace(-pwidth,pwidth,sizep)+p0;
Tlim=[0; Tf];   % time limits for ode solver
opt=odeset('MaxStep',min([TB0 1/fD])/10); % max step size in seconds
tvec = linspace(0,Tf,1000); %time vector for interpolated points
xvec = cell(sizex,sizep); %position vector for interpolated points
pvec = cell(sizex,sizep); %momentum vector for interpolated points
phasespaceprob = cell(sizex,sizep); %gaussian weighted array
disp('Intergrating classical equations of motion');
tic
parfor ix=1:sizex*sizep
    [ii,jj]=ind2sub([sizex sizep],ix);
    [T,Y]=ode45(@(t,V) dPdT(t,V),Tlim,[xinit(ii);pinit(jj)],opt); %numerically integrate
    xvec{ix} = interp1(T,Y(:,1),tvec); %interpolate for identical length vectors
    pvec{ix} = mod(interp1(T,Y(:,2),tvec)+1,2)-1;
    phasespaceprob{ix} = exp(-(xinit(ii)/xsigma)^2/2)*exp(-((pinit(jj)-p0)/psigma)^2/2); %assign probability
end
xvec = cell2mat(cellfun(@(x)reshape(x,1,1,[]),xvec,'un',0));
pvec = cell2mat(cellfun(@(x)reshape(x,1,1,[]),pvec,'un',0));
phasespaceprob=cell2mat(phasespaceprob);
toc
normalization = sum(sum(phasespaceprob,1)); %normalize probability
phasespaceprob = phasespaceprob / normalization;

%% Calculate the first and second moments of the distribution over time and
%  interpolate a probability distribution over time for a mixOD
x_avg = zeros(length(tvec),1); p_avg = x_avg;  %initialize moment vectors
xspread_avg = x_avg; pspread_avg = x_avg;
maxx = 1.2*max(xvec,[],'all'); minx = 1.2*min(xvec,[],'all'); %set mixOD array
xpoints = minx:1:maxx;
P= reshape(phasespaceprob,1,sizex*sizep); %turn probability matrix into array
OD = zeros(length(xpoints),length(tvec));  %initialize OD
for ii=1:length(tvec)
    % weight trajectories
    weightedx = xvec(:,:,ii).*phasespaceprob;
    weightedxsquared = weightedx.*xvec(:,:,ii);
    weightedp = pvec(:,:,ii).*phasespaceprob;
    weightedpsquared = weightedp.*pvec(:,:,ii);
    
    %compute moments <x> and sqrt(<x^2>-<x>^2)
    x_avg(ii) = sum(sum(weightedx,1));
    p_avg(ii) = sum(sum(weightedp,1));
    xspread_avg(ii) = sqrt(sum(sum(weightedxsquared,1)) - x_avg(ii)^2);
    pspread_avg(ii) = sqrt(sum(sum(weightedpsquared,1)) - p_avg(ii)^2);
    
    %organize probabilities and interpolate to created normalized OD
    x = xvec(:,:,ii);
    x = reshape(x,1,sizex*sizep); %turn matrix into 1D array
    binsize = 5.2;  %pixel size on camera
    xbins = min(x):binsize:max(x);
    probs = zeros(1,length(xbins));
    for jj=1:length(xbins)
        probs(jj) = sum(P(abs(x-xbins(jj))<binsize/2)); %sum probabilities in each xbin
    end
    if ii~=1
        od = interp1(xbins,probs,xpoints,'makima',0); %interpolate using modified akima spline, extrapolated points set to 0
        od = smoothdata(od,'sgolay',15); %smooth over interpolation oscillations
        od(od<0) = 0; %get rid of negative od from smoothing
    else
        od = normpdf(xpoints,0,xsigma); %since values not unique just make gaussian at t=0
    end
    OD(:,ii) = od/sum(od);  %normalize
    %plotting
%     figure(9);
%     clf;
%     plot(xpoints,od);
%     hold on;
%     plot(xbins,probs,'b.');
%     ylim([0 .05]); 
%     hold off;
%     drawnow;
%     del = .01;
%     filename = 'interpolate.gif';
%     im = frame2im(getframe(gcf));
%     [imind,cm] = rgb2ind(im,256);
%     if ii == 1
%         imwrite(imind,cm,filename,'gif','Loopcount',Inf,'DelayTime',del);
%     else
%         imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',del);
%     end
end

%% Plot OD
hF6=figure(6);
clf;
hF6.Color='w';
hF6.Position(1:2) = [400 100];
hF6.Position(3:4)=[1000 500];
ax = gca;
surf(tvec*1E3,xpoints,OD);
colormap(jet);
caxis([0 .01]);
shading flat;
view(2);
ylim([minx maxx]);
c=colorbar;
% set(gca,'ColorScale','log');
xlim([0 200]);
xlabel('Time (ms)');
ylabel('Position ($d$)','Interpreter','latex');
saveas(gcf,'psbo_semi_ODmix.fig');
saveas(gcf,'psbo_semi_ODmix.png');

%% Plotting 
% Plots for first and second moments
hF5=figure(5);
clf;
hF5.Color='w';
hF5.Position(1:2) = [100 100];
hF5.Position(3:4)=[600 600];
subplot(311);
set(gca,'box','on','linewidth',1,'fontsize',14,...
    'fontname','times');
yyaxis left
plot(tvec*1e3,x_avg,'-','linewidth',2);
ylabel('$\langle x \rangle \, (d)$','interpreter','latex');
% ylim([-250 250]);
hold on;
yyaxis right
plot(tvec*1e3,mod(p_avg+1,2)-1,'-','linewidth',.5);
ylabel('$\langle p \rangle \, (\hbar k_L)$','interpreter','latex');
ylim([-1 1]);

str=['$f_B=' num2str(round(fB0,1)) '~\mathrm{Hz}$;~'  ...
    '$T_B=' num2str(round(1E3*TB0,1)) '~\mathrm{ms}$;~'  ...
    '$U_0=' num2str(round(U0,1)) 'E_R$;~'  ...
    '$\phi_0=' num2str(round(phi0*180/pi)) '^{\circ}$;~'  ...
    '$2L_{\mathrm{WS}}=' num2str(round(L)) 'd$;~' ...
    '$q=' num2str(round(q,1)) '~\mathrm{Hz}$'];

ax=gca;
ax.Units='pixels';
% ax.Position(4)=ax.Position(4)-20;


text(-40,ax.Position(4)+10,str,'units','pixels','interpreter','latex',...
    'fontsize',12,...
    'verticalalignment','bottom',...
    'horizontalalignment','left');

% plot position spread over time
subplot(3,1,2)
ax = gca;
set(ax,'box','on','linewidth',1,'fontsize',14,...
    'fontname','times');
hold on;
yyaxis left;
plot(tvec*1e3,xspread_avg,'-','linewidth',2);
ylabel('$\sigma_x \, (d)$','interpreter','latex');
ylim([0 250]);

yyaxis right
plot(tvec*1e3,pspread_avg,'-','linewidth',.5);
ylabel('$\sigma_p (\hbar k_L)$','interpreter','latex');
ylim([0 1]);

% plot lattice depth
subplot(313)
set(gca,'box','on','linewidth',1,'fontsize',14,...
    'fontname','times');
hold on;
if modulation == 1
    latticedepth = Ufunc(tvec);
    plot(tvec*1e3,latticedepth,'k-','linewidth',1);
else
    plot([0 max(tvec)]*1e3,[1 1]*U0,'k-','linewidth',1);
end
xlabel('time (ms)','interpreter','latex');
ylabel('$V_0 (E_R)$','interpreter','latex');

figure(73);
clf;
plot(tvec*1e3,xspread_avg.*pspread_avg*2/pi,'Linewidth',1.5);
xlabel('Time (ms)')
ylabel('Uncertainty ($\hbar/2$)','interpreter','latex');
ylim([0 2]);
xlim([0 200]);

saveas(gcf,'psbo_semi_sim.fig');
saveas(gcf,'psbo_semi_sim.png');
save('data.mat','tvec','x_avg','xspread_avg','latticedepth','U0','U1','phi0',...
    'TB0','fD','dfdt','xpoints','OD');