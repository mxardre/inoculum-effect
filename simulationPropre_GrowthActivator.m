close all
%clear all
Vdrop = 4.4e-4; %ml

Ndrop=1000;
inocula = logspace(0,3,10); %range of inocula to simulate
dt = 1/60; %heures
timeSpan = 0:dt:30; %heures

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

spanThActv = logspace(-5,2,11);

rActv = mg/log(2); %the rate of production of growth activator is the inverse of the doubling time of cells.

%% synchronistation demonstration of the synchronisation effect
clear lagPop  stdLagPop   decsyncPop leadsyncPop

tmes = zeros(Ndrop,length(inocula));
k=0;
for thActv = spanThActv
    k = k +1;
    clear tmes tsync tstsdsync decsync dec lead leadsync
    
    dec = ones(Ndrop,length(inocula))*nan;
    lead = ones(Ndrop,length(inocula))*nan;
    tmes = ones(Ndrop,length(inocula))*nan;
    
    for inoc = 1:length(inocula)
        
        for i = 1: Ndrop
            
            r = round(poissrnd(inocula(inoc)));  %draw a random inoculum according to the poisson distribution.
            if r~=0
                timeSeries = ones(r,length(timeSpan))*nan; %timeSerie of the growth for each bacteria in this drop
                timeSeriesActv = zeros(r,length(timeSpan)); %timeSerie of the growth activator for each bacteria in this drop
                clear tlag
                %lognormal
                tlag = lognrnd(mu,s,r,1);  %draw randomly the cell-lag according to the corrected cell-lag distribution measured.
                grate = normrnd(mg,vg,r,1); %draw randomly the cell growth rate according to the growth rate distribution measured.
                
                
                timeSeriesActv = exp(grate.*(timeSpan-tlag))-1; %production of molecule is linear with number of cells so it follows the exp growth
                timeSeriesActv(timeSeriesActv<0) = 0; % molecule contration cannot be negative.
                timeSeriesActv = (timeSpan-tlag).*rActv.*timeSeriesActv; %multiplication by time and production rate.
                
                %calculate the total concentration of molecule produced by all
                %the cell in the drop. Need of condition for drop with one cell
                %(no need to sum over cells)
                if r>1
                    actv = sum(timeSeriesActv);
                else
                    actv = timeSeriesActv;
                end
                
                %find the time at which the concentration of molecule gets
                %above a given threshold.
                tActv = timeSpan(find(actv>thActv,1,'first'));
                
                dec(i,inoc) = tActv - min(tlag); % difference of lag time of the leader cell and the lag time due to production of molecule.
                lead(i,inoc) = sum(tlag<=tActv);
                
                
                tlagActv = tlag;
                tlagActv(tlagActv>tActv)=tActv; %every cells lag time end when the molecule gets above the threshold.
                
                
                timeSeries(:,:)= exp(grate.*(timeSpan-tlagActv)); % proceed to the exponential growth of every bacteria of this drop.
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
        decsync(i) = nanmean(dec(:,i));
        leadsync(i) = nanmean(lead(:,i));
    end
    
    lagPop(k,:)=tsync;
    stdLagPop(k,:)=tstdsync;
    decsyncPop(k,:)=decsync;
    leadsyncPop(k,:)=leadsync;
end

%%
close all



Y = repmat(log10(spanThActv)',1,10);
X = repmat(log10(inocula),11,1);
Z = lagPop;

%population lag time vs threshold and inoculum
figure('Renderer', 'painters', 'Position', [10 10 900 900]),
surf(X,Y,Z);

%title('population lag time')
xlabel('inoculum')
xticks(log10(round(inocula)))
xticklabels(round(inocula))

ylabel('threshold (au)')
yticks(log10(spanThActv(1:2:end)))
yticklabels(num2str(spanThActv(1:2:end)','%1.0e'))
zlabel('population lag time (h)')
c = colorbar;
set(gca,'FontSize',30)
c.Location='northoutside';
view(35.207879105520632,39.388548057259705)

% population lag time minus leader cell lag time
figure('Renderer', 'painters', 'Position', [10 10 900 900]),
Y = repmat(log10(spanThActv)',1,10);
X = repmat(log10(inocula),11,1);
decsyncPop(decsyncPop<=dt)=0;
Z = decsyncPop;

surf(X,Y,Z);

%title('time difference between cell leader lag time and lag time of population')
xlabel('inoculum')
xticks(log10(round(inocula)))
xticklabels(round(inocula))

ylabel('threshold (au)')
yticks(log10(spanThActv(1:2:end)))
yticklabels(num2str(spanThActv(1:2:end)','%1.0e'))
zlabel('lag pop - lag leader cell (h)')
c = colorbar;
set(gca,'FontSize',30)
caxis([0, 2]);
c.Location='northoutside';
view(35.207879105520632,39.388548057259705)

%number of cells that multiply before synchro
figure('Renderer', 'painters', 'Position', [10 10 900 900]),
Y = repmat(log10(spanThActv)',1,10);
X = repmat(log10(inocula),11,1);
Z = leadsyncPop;

surf(X,Y,Z);

%title('number of cell leaders')
xlabel('inoculum')
xticks(log10(round(inocula)))
xticklabels(round(inocula))

ylabel('threshold (au)')
yticks(log10(spanThActv(1:2:end)))
yticklabels(num2str(spanThActv(1:2:end)','%1.0e'))
zlabel('number of leader cells')
set(gca,'FontSize',30)
c = colorbar;
caxis([1 5]);
c.Limits = [1 5];
c.Ticks = [1 2 3 4 5 6 7 8];
c.Location='northoutside';
view(35.207879105520632,39.388548057259705)
