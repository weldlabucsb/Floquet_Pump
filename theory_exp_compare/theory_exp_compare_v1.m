% This codes runs comparison between the experimental data and semiclassical simulation.
% Put the corresponding atomdata.mat from the experiment and data.mat from
% psbo_semi into a single folder. The plots the theoretical and
% experimentally calculated first and second moments on top of each other.

%% Initialize
oldfolder = cd;
dirname = 'DNA';  %put name of folder with the atomdata.mat (from exp) and data.mat (from psbo_semi) 
cd(dirname);
atomd = 'atomdata.mat';
simd = 'data.mat';
load(atomd);
load(simd);
pixconv = 2.6E-6;  %pixel size to meters
d = 532E-9;  %lattice spacing in meters
x = x_avg;  %simulation calculations
t = tvec;
sigma = xspread_avg;
mean = []; std=[]; time=[];

%% Use if gaussian fits are good (for consistency probably don't use this section)

% for ii = 1:length(atomdata)
% mean = [mean -(atomdata(1).cloudCenter_y-atomdata(ii).cloudCenter_y)/d]; %mean position
% % meanpix = (atomdata.cloudCenter_y / pixconv)-atomdata.ROI(1); %mean pixel position from top
% stdx = [std atomdata(ii).cloudSD_y/d]; %width
% time = [time (atomdata(ii).vars.dummy)];
% end

%% Integrate OD to compute moments for nongaussian clouds
f=figure(7)
clf;
ax = axes;
for ii = 1:length(atomdata)
    thisODy = atomdata(ii).summedODy;
    smoothedODy = smoothdata(thisODy,'sgolay',25); %smooth data
    normalizedODy = smoothedODy / sum(smoothedODy);
    normalizedODy(normalizedODy<.001)=0; %remove negative OD, can set the threshold
    normalizedODy = normalizedODy / sum(normalizedODy); %normalize
    hold off;
    pos = 1:length(thisODy);
    thismean = sum(normalizedODy.*pos);  %compute first moment
    mean = [mean thismean];
    delta = pos-thismean;
    thisstd = sqrt(sum((normalizedODy.*(delta.^2)))); %compute second moment
    std = [std thisstd];
    time = [time (atomdata(ii).vars.dummy)];
    
    % Plot distribution processing
    plot(thisODy,'Linewidth',1.5);
    hold on;
    plot(smoothedODy,'Linewidth',1.5);
    plot(normalizedODy * sum(smoothedODy),'Linewidth',1.5);
    hold off;
    ylabel('OD');
    xlabel('Position (au)');
    ylim([-500 3000]);
    xlim([0 length(thisODy)]);
    legend('Raw','Filtered','Noise Floor Removed','Location','NW');
    drawnow;
        del = .25;
        filename = 'smoothing.gif';
        frame = getframe(gcf);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if ii == 1;
            imwrite(imind,cm,filename,'gif','Loopcount',Inf,'DelayTime',del);
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',del);
        end
end
mean = (mean - mean(1))*pixconv/d;
std = std*pixconv/d;

%% Plotting of data vs semiclassical simulation
f = figure(10);
clf;
ax1 = axes; ax1.Position=[.12 .57 .8 .35];
set(gca,'box','on','linewidth',1.5,'fontsize',16,...
    'fontname','times');
plot(tvec*1E3,x_avg,'-b','Linewidth',1.5);
hold on;
plot(time,mean,'.k','MarkerSize',15);
ylabel('Mean Position ($d$)','Interpreter','Latex','fontsize',12);
% title('No frequency Ramp, DNA, $\varphi_0=3\pi/2$','Interpreter','Latex','fontsize',12)

ax2 = axes; ax2.Position=[.12 .1 .8 .4];
set(gca,'box','on','linewidth',1.5,'fontsize',16,...
    'fontname','times');
plot(tvec*1E3,xspread_avg,'-b','Linewidth',1.5);
hold on;
plot(time,std,'.k','MarkerSize',15);
ylabel('Cloud Width ($d$)','Interpreter','Latex','fontsize',12);
xlabel('Time (ms)','Interpreter','Latex','fontsize',12);
ylim([0 250]);


cd(oldfolder);