close all
%clear all
Vdrop = 4.4e-4; %ml

Ndrop=1000;
inocula = logspace(0,3,10); %range of inocula to simulate
dt = 1/60;
timeSpan = 0:dt:30;

Nthresh = 1.6e8; %cell/ml. Threshold to calculate the lag time like in the experiments.

stdNoise = 0.88; %variance of the noise 0.87 with lambda and 2 sigma for the calib
stdNoiseTitle = stdNoise;
MM = 6.8246; %mean value of the Exp distribution of lag for logn
VV = 1.3322^2 - stdNoise^2; %std of the corrected distribution of the experimental lag that gollows a logn
%calculate the new param of the logn to remove the noise
mu = log(MM^2 / sqrt(MM^2 + VV));
s = sqrt(log(VV/MM^2 + 1));


mg = 0.8430; %average grate %1/h from inoculum 1
vg = 0.02; %variance grate from inoculum 1

%% simulation of the simple growth of bacteria in the droplets. Measure of the lag for eqch droplet.

tmes = zeros(Ndrop,length(inocula));

clear timeSeries tmes
for inoc = 1:length(inocula)  %loop over the inoculum
    for i = 1: Ndrop %loop over the drop
        
        r = round(poissrnd(inocula(inoc))); %draw a random inoculum according to the poisson distribution.
        
        timeSeries = ones(r,length(timeSpan))*nan; %timeSerie of the growth for each bacteria in this drop
        clear tlag
        tlag = lognrnd(mu,s,r,1);  %draw randomly the cell-lag according to the corrected cell-lag distribution measured.
        
        grate = normrnd(mg,vg,r,1); %draw randomly the cell growth rate according to the growth rate distribution measured.
        
        timeSeries(:,:)= exp(grate.*(timeSpan-tlag)); % proceed to the exponential growth of every bacteria of this drop.
        timeSeries(timeSeries<1)=1; % the timeseries must start at 1 before the division of the bacteria
        
        %calculation of the cell concentration in this drop along the
        %time by suming the timeSerie of its bacteria
        l = size(timeSeries);
        if l(1)==1
            totDrop = timeSeries/Vdrop;
        else
            totDrop = nansum(timeSeries)/Vdrop;
        end
        
        %measure the lag of the droplet by finding the time at which the
        %cell concentration gets above Nth like in the experiments
        if r~=0
            tau = timeSpan(find(totDrop>Nthresh,1,'first'));
            if isempty(tau)
                tmes(i,inoc) = nan;
            else
                tmes(i,inoc) = timeSpan(find(totDrop>Nthresh,1,'first'))-log(Nthresh*Vdrop/r)/nanmean(grate);
            end
        else
            tmes(i,inoc) = nan;
        end
    end
end

for i = 1:length(inocula)
    tstat(i) = nanmean(tmes(:,i));
    tstsdstat(i) = nanstd(tmes(:,i));
end

%% synchronistation demonstration of the synchronisation effect
clear tmes tsync tstsdsync
tmes = zeros(Ndrop,length(inocula));
for inoc = 1:length(inocula)
    
    for i = 1: Ndrop
        
        r = round(poissrnd(inocula(inoc)));  %draw a random inoculum according to the poisson distribution.
        
        if r~=0
            timeSeries = ones(r,length(timeSpan))*nan; %timeSerie of the growth for each bacteria in this drop
            
            clear tlag
            %lognormal
            tlag = lognrnd(mu,s,r,1);  %draw randomly the cell-lag according to the corrected cell-lag distribution measured.
            
            tlag = ones(r,1) .* (min(tlag)); %set all the cell-lags to the minimal cell-lag of the leader
            
            grate = normrnd(mg,vg,r,1); %draw randomly the cell growth rate according to the growth rate distribution measured.
            
            
            timeSeries(:,:)= exp(grate.*(timeSpan-tlag)); % proceed to the exponential growth of every bacteria of this drop.
            timeSeries(timeSeries<1)=1; % the timeseries must start at 1 before the division of the bacteria
            
            %calculation of the cell concentration in this drop along the
            %time by suming the timeSerie of its bacteria
            l = size(timeSeries);
            if l(1)==1
                totDrop = timeSeries/Vdrop;
            else
                totDrop = nansum(timeSeries)/Vdrop;
            end
            
            %measure the lag of the droplet by finding the time at which the
            %cell concentration gets above Nth like in the experiments
            tau = timeSpan(find(totDrop>Nthresh,1,'first'));
            if isempty(tau)
                tmes(i,inoc) = nan;
            else
                tmes(i,inoc) = timeSpan(find(totDrop>Nthresh,1,'first'))-log(Nthresh*Vdrop/r)/nanmean(grate);
            end
        else
            tmes(i,inoc) = nan;
        end
    end
end

for i = 1:length(inocula)
    tsync(i) = nanmean(tmes(:,i));
    tstdsync(i) = nanstd(tmes(:,i));
end

%% synchronistation without noise

stdNoise = 0; %variance of the noise

MM = 6.8246; %mean value of the Exp distribution of lag for logn
VV = 1.3322 - stdNoise^2; %std of the Exp distribution of lag for logn
%calculate the new param of the logn to remove the noise
mu = log(MM^2 / sqrt(MM^2 + VV));
s = sqrt(log(VV/MM^2 + 1));

clear tmes tsyncNoNoise tstsdsyncNoNoise
tmes = zeros(Ndrop,length(inocula));
for inoc = 1:length(inocula)
    
    for i = 1: Ndrop
        
        r = round(poissrnd(inocula(inoc)));  %draw a random inoculum according to the poisson distribution.
        
        if r~=0
            timeSeries = ones(r,length(timeSpan))*nan; %timeSerie of the growth for each bacteria in this drop
            
            clear tlag
            %lognormal
            tlag = lognrnd(mu,s,r,1);  %draw randomly the cell-lag according to the corrected cell-lag distribution measured.
            
            tlag = ones(r,1) .* (min(tlag)); %set all the cell-lags to the minimal cell-lag of the leader
            
            grate = normrnd(mg,vg,r,1); %draw randomly the cell growth rate according to the growth rate distribution measured.
            
            
            timeSeries(:,:)= exp(grate.*(timeSpan-tlag)); % proceed to the exponential growth of every bacteria of this drop.
            timeSeries(timeSeries<1)=1; % the timeseries must start at 1 before the division of the bacteria
            
            %calculation of the cell concentration in this drop along the
            %time by suming the timeSerie of its bacteria
            l = size(timeSeries);
            if l(1)==1
                totDrop = timeSeries/Vdrop;
            else
                totDrop = nansum(timeSeries)/Vdrop;
            end
            
            %measure the lag of the droplet by finding the time at which the
            %cell concentration gets above Nth like in the experiments
            tau = timeSpan(find(totDrop>Nthresh,1,'first'));
            if isempty(tau)
                tmes(i,inoc) = nan;
            else
                tmes(i,inoc) = timeSpan(find(totDrop>Nthresh,1,'first'))-log(Nthresh*Vdrop/r)/nanmean(grate);
            end
        else
            tmes(i,inoc) = nan;
        end
    end
end

for i = 1:length(inocula)
    tsyncNoNoise(i) = nanmean(tmes(:,i));
    tstdsyncNoNoise(i) = nanstd(tmes(:,i));
end

%% plot the curves of stat effect vs sync effect need to run the cell above first.

%experimental points
lag =  [ 6.564471  6.018853  5.392048  4.729353  4.375776  4.259628; ...
    6.349308  5.901676  5.470768  4.871843  4.669578  4.680528; ...
    6.377465  5.937269  5.553422  5.204688  4.640315  4.353197]';

slag = [...
    0.843766  0.700523 0.491303 0.519735  0.203902  0.199777;...
    1.445141  0.974685  0.682811  0.469918  0.363288  0.316424; ...
    1.120196  0.956412  0.731424  0.720511  0.455629  0.434152 ]';

N0 = repmat([1 4 16 64 256 1024], size(lag,2), 1)';


%plot the figure

figure('Renderer', 'painters', 'Position', [10 10 600 600])
hold on

alpha = 0.1;
yRd = tsync; % your mean vector;
x = log(inocula);
std_dev = tstdsync;
curve1 = yRd + std_dev;
curve2 = yRd - std_dev;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'b','FaceAlpha',alpha);
plot(x, yRd, '-b', 'LineWidth', 3,'MarkerSize',20,'DisplayName', 'Leader');

y = tstat; % your mean vector;
x = log(inocula);
std_dev = tstsdstat;
curve1 = y + std_dev;
curve2 = y - std_dev;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'b','FaceAlpha',alpha);
plot(x, y, '-.b', 'LineWidth', 3,'MarkerSize',20,'DisplayName', 'Stat');


y = tsyncNoNoise; % your mean vector;
x = log(inocula);
std_dev = tstdsyncNoNoise;
curve1 = y + std_dev;
curve2 = y - std_dev;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'k','FaceAlpha',alpha);
plot(x, y, ':k', 'LineWidth', 3,'MarkerSize',20,'DisplayName','No Noise');

clr = {[0.6350 0.0780 0.1840] [0.4660 0.8 0.1880] 'k'};
mkr = {'d' 'x' '*'};
for i = 1:size(lag,2)
    h = errorbar(log(N0(:,i)), lag(:,i), slag(:,i),mkr{i}, 'MarkerSize',20, 'color', clr{i}, 'LineWidth',3, 'CapSize', 40);
end
hold off


xlim([-0.2 log(1500)])
ylim([3 8])
xticks( log(N0(:,1)));
xticklabels(N0(:,1));
xlabel('inoculum (cell/drop)');
ylabel('mean lag time (h)');
title(['N0=1 noise std=' num2str(stdNoiseTitle) ])

box('on')
grid('on')
set(gca,'FontName','Helvetica')
set(0, 'defaultAxesFontSize', 40);
set(gca,'LineWidth',2) 
set(gca,'GridAlpha', 0.1)
p=0.3;
scale=.6;
set(gca,'position',[p,p,scale,scale])

