% blohOscillate.m
% Author      : CJ Fujiwara
% Last Edited : 2019/10/05
% 
%
% This code is meant to simulate bloch oscillations in a slowly modulated band
% structure

close all;
clear all;
doSave=1;
filename='BO.gif';

%% experiment parameters

% Specify fundamental constants and other fixed parameters
amu=1.66E-27;       % atomic mass unit in kg
mLi=7*amu;          % mass of lithium in kg
h=6.626E-34;        % planck's constant in SI (Joules*second)
hbar = h/(2*pi);
d=532E-9;            % lattice spacing in m
kL = pi/d;

%%%%% Specification of initial wavefunction %%%%%
% Spatial extent of wavefunction
% Remember that a lattice spacing in these units is pi
s=60*pi;   %initial width in lattice site * pi

% initial eigenstate
k=1;           % initial quasimomentum
n=1;            % initial band index

%%%%% lattice parameters %%%%%

% lattice depth
depth=4.345;        % the lattice depth in Er
fR=25.127*1E3;   % recoil frequency in Hz
wR=2*pi*fR;     % angular frequency of lattice
Er = h*fR;

% force as expressed in bloch frequency
% fB=29.1;        % bloch frequency in Hz
% fB=150;           % Aproximately our fastest bloch frequency
tB0 = 18E-3;          % bloch period in s
fB0 = 1/tB0;
F0 = h*fB0/d;

% lattice drive
fD = 53.56;                % drive frequency in Hz
wD = 2*pi*fD;
theta = pi/2;         % initial drive phase in rad
A = 1.02/depth;            % modulation amplitude

% Quadratic portion of magnetic potential 
fho=15.5;         % harmonic oscillator frequency in Hz
who = 2*pi*fho;

%%%%% Convert experimental parameters into lattice recoils

% Convert the experimental values into fundamental units of Er
Q=(mLi*who^2*d^2)/(2*Er);         % curvature frequency scaled in Er
W = wD/wR;      % drive frequency in Er

% Display the simulation parameters
disp(' ');
disp('blochOscillate.m');
disp(' ');
disp('Simulating Bloch oscillations for the following parameters...');
disp(['lattice depth             : ' num2str(round(depth,1)) ' Er']);
disp(['bloch frequency           : ' num2str(round(fB0,1)) ' Hz']);
disp(['curvature                 : ' num2str(fho) ' Hz']);
disp(['initial band              : ' num2str(n) ]);
disp(['initial quasimomentum     : ' num2str(k)  ' recoils']);

%% theory calculations
% Calcuate and display relevant theoretical values

% Spectrum at zone center
H1=makeHmatrix(0,depth*(1+A*sin(theta)));                      
[~,eCen]=eig(full(H1));       
eCen=diag(eCen);

% spectrum at zone edge
H2=makeHmatrix(1,depth*(1+A*sin(theta)));                         
[~,eEdge]=eig(full(H2));                             
eEdge=diag(eEdge);

% Find the bandwidths
BW=abs(eEdge-eCen);
% Compensate for initial half BO
J = Er*BW(1)/4;   %tunneling energy in Joules
lws = 2*J/F0;                   %localization length in meters
F = F0 - mLi*who^2*lws;         %actual force at start of the drive
fB = F*d/h;
F = fB/fR;
% lattice potential
v_mod=@(phi,t) -cos(2*phi)*(depth/2)*(1+A*sin(W*t+theta))-...       
    phi/pi*F+Q/pi^2 *phi.^2;

%% Simulation Parameters
% Define the parameters of the simulation
% Utype='RK4';

% Define position mesh
nSites=[-400 600];              % lattice site vector
nPerSite=5;                    % points per site (sets max momentum)
X=(linspace(nSites(1),nSites(2),...
    nPerSite*range(nSites))*pi)';
dX=X(2)-X(1);                   % Separation between positions

% Define time vectors. Time is in phase angle of hbar*omega_R.
% T=abs(2*pi/F);                  % Bloch period
tau=120E-3*wR;             

dtau=.4*dX^2;                   % 0.4 for RK4 stability (basically need 
%that dT/dX^2 be small for the kinetic energy matrix to have reasonable
%values that don't cause divergeneces, smaller is better)
tVec=0:dtau:tau;                % Time vector
dtau=tVec(2)-tVec(1);           % Time spacing

disp(' ');
disp(['number of sites           : ' num2str(nSites)]);
disp(['points per site           : ' num2str(nPerSite)]);
disp(['max time                  : ' num2str(tau/(pi/F)) ' bloch periods']);
disp(['dT/dX^2                   : ' num2str(dtau/dX^2)]);
disp(' ');
%% Make wavefunction over all space
fprintf('Initializing wavefuction...');
maxBand=20;                                     % Max band index   
H=makeHmatrix(k,depth);                         % Static hamiltonian
[v,b]=eig(full(H));                             % Eigenstates
lambda=b*ones(length(b),1);                     % Eigenvalues
lambda=lambda(1:maxBand);                       % Look at subspace
un=makeBlochState(v(:,n));                      % Create u_{n,k}

% Wavefunction is u_(n,k)*exp(ikx) with the gaussian overlap
Y=un(X).*exp(-X.^2/(4*s^2)).*exp(1i*k*X);
Y=Y/(sqrt(sum(conj(Y).*Y)));

disp('done');
%% Make Hamiltonians and time Evolution Operator
disp(' ');
fprintf('creating matix operators...');

% identity operator
Imat=speye(length(X));           
% kinetic energy operator
Tmat=gallery('tridiag',...          
    length(X),-1,2,-1)/dX^2;

disp('done');
disp(' ');
%% Initialize figure for live update
fprintf('initializing figure for live update...');

hF2=figure(2);
clf
set(hF2,'color','w');
hF2.Position(3)=600;
hF2.Position(4)=400;
hF2.Position(2)=50;

% subplot 1 : real space wavefunction
hAx1=subplot(221);
set(gca,'xgrid', 'off', 'ygrid', 'off',...
    'box','on','YTickLabel',{},'FontSize',12,'xdir','normal',...
    'fontname','times');
ylabel('$\Psi(\phi)$','interpreter','latex');
xlabel('position (sites)','interpreter','latex');
xlim([min(X/pi) max(X/pi)]);
yMax=max([max(real(Y)) max(imag(Y))]);
ylim([-1.2*yMax 1.2*yMax]);
hold on

pPsiR=plot(X/pi,real(Y),'.-');      % real part
pPsiC=plot(X/pi,imag(Y),'k.-');     % imaginary part


% subplot 2 : real space density
hAx2=subplot(222);
set(gca,'xgrid', 'off', 'ygrid', 'off',...
    'box','on','YTickLabel',{},'FontSize',12,'xdir','normal');
xlim([min(X/pi) max(X/pi)]);
xlabel('position (sites)','interpreter','latex');
ylim([0 1.1*max(real(conj(Y).*Y))]);
ylabel('$|\Psi(\phi)|^2$','interpreter','latex');
xlim(nSites);

hold on
pPsiD=plot(X/pi,Y.*Y,'.-');     % probability density

% subplot 3 : momentum space density
[kV,P]=computeFFT(X,Y);
hAx4=subplot(223);
set(gca,'xgrid', 'off', 'ygrid', 'off',...
    'box','on','YTickLabel',{},'FontSize',12,'xdir','normal');
xlabel('momentum ($\hbar k_L$)','interpreter','latex')
ylabel('|\Psi(k)|^2');
ylim([0 1.2*max(P)]);
xlim([min(kV) max(kV)]*pi);
hold on

pK=plot(kV*pi,P,'.-');      % momentum space density

drawnow

 
disp('done');

%% Time Evolve
disp(' ');
disp('beginning time evolution!!!!');

nStp=0;
inc=.001;
tt=[];

YMat=Y;
disp(length(tVec))
tic
for kk=1:length(tVec)
    t=tVec(kk);     % the time
    % time-dependent potential
    Vmat=sparse(1:length(X),1:length(X),...
    v_mod(X,t),...
    length(X),length(X));
    % make Hamiltonian
    Hmat=sparse(Tmat+Vmat);

    % make propagator
    U=makeRK4(Hmat,dtau,X);  %4th order Runge-Kutta method
%     U=makeCK(Hmat,dtau,X);   %2nd order midpoint Crank Nicholson method
    
    Y=U*Y;          % apply the TDSE

    % Check to update the plot
    if kk/length(tVec)>=nStp*inc
        nStp=nStp+1;
        tt(end+1)=tVec(kk);
        
        % calculate observables
        [~,Pk]=computeFFT(X,Y);        
        Px=real(conj(Y).*Y);
        YMat(1:length(Y),nStp)=Px;
        
        % update each plot
        pPsiR.YData=real(Y);
        pPsiC.YData=imag(Y);            
        pPsiD.YData=Px;
        pK.YData=Pk;     
        
        drawnow;
        
        if doSave        
%             filename='.gif';
            frame = getframe(hF2);        
            im = frame2im(frame);
            [A,map] = rgb2ind(im,256);    
             if kk == 1
                imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1.5);
                kk=0;
             else            
                    imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',.1);
            end
        end
    end
end
toc
%% Compute moments
t = tt'/wR;   %time vector in s
mean=X'/pi*YMat;   %average position in lattice sites
spread = sqrt((X'/pi).^2 * YMat - mean.^2);
%% Plotting of the Mean Position
hF3=figure(12);
clf
hF3.Position(3:4)=[700 300];
hF3.Color='w';
hF3.Name='center of mass motion';

subplot(131)
imagesc(t*1E3,X/pi,YMat)
set(gca,'fontname','times','fontsize',12,'box','on',...
    'linewidth',1);
colormap(flip(bone))
drawnow

ylabel('position (sites)');
xlabel('Time (ms)','interpreter','latex')
hold on
caxis([0 5E-4]);

subplot(132)
cla
plot(t*1E3,mean,'color','k','linewidth',2);
set(gca,'fontname','times','fontsize',12,'box','on',...
    'linewidth',1);

ylabel('Mean Position ($d$)','interpreter','latex');
xlabel('Time (ms)','interpreter','latex')
hold on

xlim([min(t) max(t)]*1E3);
ylim([min(mean)*1.2 max(mean)*1.2]);

caxis([0 5E-4]);

subplot(133);
plot(t*1E3,spread,'color','k','linewidth',2);
set(gca,'fontname','times','fontsize',12,'box','on','linewidth',1);
ylabel('Spread ($d$)','interpreter','latex');
xlabel('Time (ms)','interpreter','latex');
saveas(gcf,'moments.fig');
saveas(gcf,'moments.png');
save('data.mat','t','mean','spread','depth','A','theta',...
    'tB0','fD');