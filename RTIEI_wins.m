clear all; close all;
display('Importing...');
%load('/Users/joshsalvi/Documents/Lab/Lab/Simulation Data/ONHFishJosh/Analysis/Fp7-kp2-start1end2e5-xfish1.0Noise.mat');
load('/Users/joshsalvi/Documents/Lab/Lab/Simulation Data/Sinusoids/noisysinewave22movingavg.mat');

wins = 1:10:150;                                    % window range (for minimum time between events), pts
nup = -3:0.2:3;
ndown = -3:0.2:3;

display('Find up events...');
up_events(1,1:length(nup),1:length(wins))=0;    %initialize
up_ind(1:length(movingavgX),1:length(nup),1:length(wins))=0;
for l = 1:length(wins)
for j = 1:length(nup)
    Mvec_up = zeros(1,length(Xvec));
for i = 1:length(movingavgX)
    if i <= win
        clear a;
        a = Xvec(1:round(win/2)+i-1) > movingavgX(i) + nup(j)*movingstdX(i);            % left edge
        Mvec_up(1:round(win/2)+i-1) = Mvec_up(1:round(win/2)+i-1)+a(:)';
    else if i < length(movingavgX) - win
            clear a;
            a = Xvec(i-floor(win/2):round(win/2)+i-1) > movingavgX(i) + nup(j)*movingstdX(i);           % middle
            Mvec_up(i-floor(win/2):round(win/2)+i-1) = Mvec_up(i-floor(win/2):round(win/2)+i-1)+a(:)';
        else
            clear a;
            a = Xvec(i-floor(win/2):end) > movingavgX(i) + nup(j)*movingstdX(i);          % right edge
            Mvec_up(i-floor(win/2):end) = Mvec_up(i-floor(win/2):end)+a(:)';
        end
    end
end
Mvec_up(Mvec_up > 0) = 1;   % reduce multiple counts to 1
up_ind(:,j,l) = Mvec_up;
clear updiffend updiffstart updiff updiffind
updiffend = find(diff(up_ind(:,j,l))==-1);           % concatenate nearby points, stepping through window sizes
updiffstart = find(diff(up_ind(:,j,l))==1);

if isempty(updiffend)==0 || isempty(updiffstart)==0
if updiffend(1)<updiffstart(1)
    updiffend(1)=[];
end
for k = 2:length(updiffend)
    if abs(updiffend(k) - updiffstart(k-1)) < wins(l)
        if updiffend(k)>updiffstart(k-1)
            up_ind(updiffstart(k-1):updiffend(k),j,l) = 1;
        else
            up_ind(updiffend(k-1):updiffstart(k),j,l) = 1;
        end
    end
end
end
%}
up_events(:,j,l) = floor(sum(abs(diff(up_ind(:,j,l))))/2);
end
if wins(l) == wins(round(length(wins)/2))
    display('50%');
end
end
events_up_diff = diff(up_events);


display('Find down events...');
down_events(1,1:length(ndown),1:length(wins))=0;
down_ind(1:length(movingavgX),1:length(ndown),1:length(wins))=0;
for l = 1:length(wins)
for j = 1:length(ndown)
    Mvec_down = zeros(1,length(Xvec));
for i = 1:length(movingavgX)
    if i <= win
        clear a;
        a = Xvec(1:round(win/2)+i-1) < movingavgX(i) - ndown(j)*movingstdX(i);            % left edge
        Mvec_down(1:round(win/2)+i-1) = Mvec_down(1:round(win/2)+i-1)+a(:)';
    else if i < length(movingavgX) - win
            clear a;
            a = Xvec(i-floor(win/2):round(win/2)+i-1) < movingavgX(i) - ndown(j)*movingstdX(i);           % middle
            Mvec_down(i-floor(win/2):round(win/2)+i-1) = Mvec_down(i-floor(win/2):round(win/2)+i-1)+a(:)';
        else
            clear a;
            a = Xvec(i-floor(win/2):end) < movingavgX(i) - ndown(j)*movingstdX(i);          % right edge
            Mvec_down(i-floor(win/2):end) = Mvec_down(i-floor(win/2):end)+a(:)';
        end
    end
end
Mvec_down(Mvec_down > 0) = 1;   % reduce multiple counts to 1
down_ind(:,j,l) = Mvec_down;
clear downdiffend downdiffstart downdiff downdiffind
downdiffend = find(diff(down_ind(:,j,l))==-1);           % concatenate nearby points, stepping through window sizes
downdiffstart = find(diff(down_ind(:,j,l))==1);

if isempty(downdiffend)==0 || isempty(downdiffstart)==0
if downdiffend(1)<downdiffstart(1)
    downdiffend(1)=[];
end
for k = 2:length(downdiffend)
    if abs(downdiffend(k) - downdiffstart(k-1)) < wins(l)
        if downdiffend(k)>downdiffstart(k-1)
            down_ind(downdiffstart(k-1):downdiffend(k),j,l) = 1;
        else
            down_ind(downdiffend(k):downdiffstart(k-1),j,l) = 1;
        end
    end
end
end
%}
down_events(:,j,l) = floor(sum(abs(diff(down_ind(:,j,l))))/2);
end
if wins(l) == wins(round(length(wins)/2))
    display('50%');
end
end
events_down_diff = diff(down_events);


for l = 1:length(wins)
    clear nupconv_hold ndownconv_hold
uppeaks{l} = findpeaks(up_events(:,:,l));
downpeaks{l} = findpeaks(down_events(:,:,l));
if length(uppeaks{l}) > 1
    nupconv_hold = find(abs(diff(up_events(1,uppeaks{l}(end-1):uppeaks{l}(end),l)))==min(abs(diff(up_events(1,uppeaks{l}(end-1):uppeaks{l}(end),l)))))+uppeaks{l}(end-1);
else
    nupconv_hold = uppeaks{l}(1);
end
if length(downpeaks{l}) > 1
    ndownconv_hold = find(abs(diff(down_events(1,downpeaks{l}(end-1):downpeaks{l}(end),l)))==min(abs(diff(down_events(1,downpeaks{l}(end-1):downpeaks{l}(end),l)))))+downpeaks{l}(end-1);
else
    ndownconv_hold = downpeaks{l}(1);
end
nupconv(l) = nupconv_hold(1);
ndownconv(l) = ndownconv_hold(1);
end


%uppeakswins = findpeaks((diff(diag(squeeze(up_events(1,nupconv(),:))'))))-1;
%downpeakswins = findpeaks((diff(diag(squeeze(down_events(1,ndownconv(),:))'))))-1;
%uppeakswins = find(diff(diff((diag(squeeze(up_events(1,nupconv(),:))'))))==max(diff(diff((diag(squeeze(up_events(1,nupconv(),:))'))))));
%downpeakswins = find(((diag(squeeze(down_events(1,ndownconv(),:))')))==min(((diag(squeeze(down_events(1,ndownconv(),:))')))));
%uppeakswins = find(((diag(squeeze(up_events(1,nupconv(),:))')))==min(((diag(squeeze(up_events(1,nupconv(),:))')))));
%downpeakswins = find(((diag(squeeze(down_events(1,ndownconv(),:))')))==min(((diag(squeeze(down_events(1,ndownconv(),:))')))));
%uppeakswins = find(sign(diff(diag(squeeze(up_events(1,nupconv(),:))')))==0);
%downpeakswins = find(sign(diff(diag(squeeze(down_events(1,nupconv(),:))')))==0);

%nupconvwins = uppeakswins(1);
%ndownconvwins = downpeakswins(1);

upstd = std(diag(squeeze(up_events(1,nupconv(),length(wins)/2:end))'));
downstd = std(diag(squeeze(down_events(1,ndownconv(),length(wins)/2:end))'));
upmean = mean(diag(squeeze(up_events(1,nupconv(),length(wins)/2:end))'));
downmean = mean(diag(squeeze(down_events(1,ndownconv(),length(wins)/2:end))'));
nupconvrange = find(diag(squeeze(up_events(1,nupconv(),:))') > upmean-upstd & diag(squeeze(up_events(1,nupconv(),:))') < upmean+upstd);
ndownconvrange = find(diag(squeeze(down_events(1,ndownconv(),:))') > downmean-downstd & diag(squeeze(down_events(1,ndownconv(),:))') < downmean+downstd);
uppeakswins = find(((diag(squeeze(up_events(1,nupconv(),nupconvrange))')))==min(((diag(squeeze(up_events(1,nupconv(),nupconvrange))')))));
downpeakswins = find(((diag(squeeze(down_events(1,ndownconv(),ndownconvrange))')))==min(((diag(squeeze(down_events(1,ndownconv(),ndownconvrange))')))));

nupconvwins=nupconvrange(uppeakswins(1));
ndownconvwins=ndownconvrange(downpeakswins(1));


figure(2);
subplot(3,2,1);plot(nup,up_events(:,:,nupconvwins));title('Up Events'); ylabel('Events'); xlabel('nup');
subplot(3,2,2);plot(ndown,down_events(:,:,ndownconvwins));title('Down Events'); ylabel('Events'); xlabel('ndown');
subplot(3,2,3);plot(nup(2:end),events_up_diff(:,:,nupconvwins)); ylabel('dEvents'); xlabel('nup');
subplot(3,2,4);plot(ndown(2:end),events_down_diff(:,:,nupconvwins)); ylabel('dEvents'); xlabel('nup');
subplot(3,2,1);hold on; scatter(nup(nupconv(nupconvwins)),up_events(:,nupconv(nupconvwins),nupconvwins),'ro');scatter(nup(uppeaks{nupconvwins}),up_events(:,uppeaks{nupconvwins},nupconvwins),'k.');
subplot(3,2,2);hold on; scatter(ndown(ndownconv(ndownconvwins)),down_events(:,ndownconv(ndownconvwins),ndownconvwins),'ro');scatter(ndown(downpeaks{ndownconvwins}),down_events(:,downpeaks{ndownconvwins},ndownconvwins),'k.');
subplot(3,2,5);plot(wins,squeeze(up_events(:,nupconv(nupconvwins),:)));hold on; scatter(wins(nupconvwins),up_events(1,nupconv(nupconvwins),nupconvwins),'ro');title('up windows');ylabel('events');xlabel('window size(pts)');
subplot(3,2,6);plot(wins,squeeze(down_events(:,ndownconv(ndownconvwins),:)));hold on; scatter(wins(ndownconvwins),down_events(1,ndownconv(ndownconvwins),ndownconvwins),'ro'); title('down windows');ylabel('events');xlabel('window size(pts)');
%}

%{
% fit exponential and find time to pass Ntau time constants
Ntau=1;
[fitup, gofup] = fit(wins',diag(squeeze(up_events(1,nupconv(),:))'),'exp1');
[fitdown, gofdown] = fit(wins',diag(squeeze(down_events(1,ndownconv(),:))'),'exp1');
nupconvfit = round(Ntau/abs(fitup.b));
ndownconvfit = round(Ntau/abs(fitdown.b));
nupconvwins = findnearest(nupconvfit,wins);
ndownconvwins = findnearest(ndownconvfit,wins);
nupconvwins = nupconvwins(1);
ndownconvwins = ndownconvwins(1);
if isempty(nupconvwins) == 1
    nupconvwins = 1;
end
if isempty(ndownconvwins) == 1
    ndowconvwins = 1;
end


figure(2);
subplot(3,2,1);plot(nup,up_events(:,:,nupconvwins));title('Up Events'); ylabel('Events'); xlabel('nup');
subplot(3,2,2);plot(ndown,down_events(:,:,ndownconvwins));title('Down Events'); ylabel('Events'); xlabel('ndown');
subplot(3,2,3);plot(nup(2:end),events_up_diff(:,:,nupconvwins)); ylabel('dEvents'); xlabel('nup');
subplot(3,2,4);plot(ndown(2:end),events_down_diff(:,:,nupconvwins)); ylabel('dEvents'); xlabel('nup');
subplot(3,2,1);hold on; scatter(nup(nupconv(nupconvwins)),up_events(:,nupconv(nupconvwins),nupconvwins),'ro');scatter(nup(uppeaks{nupconvwins}),up_events(:,uppeaks{nupconvwins},nupconvwins),'k.');
subplot(3,2,2);hold on; scatter(ndown(ndownconv(ndownconvwins)),down_events(:,ndownconv(ndownconvwins),ndownconvwins),'ro');scatter(ndown(downpeaks{ndownconvwins}),down_events(:,downpeaks{ndownconvwins},ndownconvwins),'k.');
subplot(3,2,5);plot(wins,squeeze(up_events(:,nupconv(nupconvwins),:)));hold on; scatter(wins(nupconvwins),up_events(1,nupconv(nupconvwins),nupconvwins),'ro');;ylabel('events');xlabel('window size(pts)');
subplot(3,2,6);plot(wins,squeeze(down_events(:,ndownconv(ndownconvwins),:)));hold on; scatter(wins(ndownconvwins),down_events(1,ndownconv(ndownconvwins),ndownconvwins),'ro'); title('down windows');ylabel('events');xlabel('window size(pts)');
subplot(3,2,5);hold on;plot(fitup,'fit');title(sprintf('%s %s%s','up windows','exp fit R2= ',num2str(gofup.rsquare)))
subplot(3,2,6);hold on;plot(fitdown,'fit');title(sprintf('%s %s%s','down windows','exp fit R2= ',num2str(gofdown.rsquare)))
%}

figure(3);
plot(time,Xvec);hold on; ylabel('displacement');
plot(time(down_ind(:,ndownconv(ndownconvwins),ndownconvwins)==1),Xvec(down_ind(:,ndownconv(ndownconvwins),ndownconvwins)==1),'r.');
plot(time(up_ind(:,nupconv(nupconvwins),nupconvwins)==1),Xvec(up_ind(:,nupconv(nupconvwins),nupconvwins)==1),'g.');
title(sprintf('%s %s%s %s%s %s%s%s%s','Spikes ','nup= ',num2str(nup(nupconv(nupconvwins))),'ndown= ',num2str(ndown(ndownconv(ndownconvwins))),'eventsup/eventsdown= ',num2str(up_events(1,nupconv(nupconvwins),nupconvwins)),'/',num2str(down_events(1,ndownconv(ndownconvwins),ndownconvwins))));


% Find residence times, buffer times, etc....
display('Find residence times...');
upindex = up_ind(:,nupconv(nupconvwins),nupconvwins);
downindex = down_ind(:,ndownconv(ndownconvwins),ndownconvwins);          % indices for up and down positions


upinddiff = diff(upindex);
downinddiff = diff(downindex);     
spikeindup_start = find(upinddiff==1);
spikeindup_end = find(upinddiff==-1);

if spikeindup_start(1) > spikeindup_end(1)
    spikeindup_end(1) = [];
end

spikeinddown_start = find(downinddiff==1);
spikeinddown_end = find(downinddiff==-1);

if spikeinddown_start(1) > spikeinddown_end(1)
    spikeinddown_end(1) = [];
end


% Calculate residence times for the up position and for the down position
for i = 1:length(spikeindup_start)-1
    RT_up(i) = (time(spikeindup_end(i))-time(spikeindup_start(i))) * 1e3;         % convert to ms
end
for i = 1:length(spikeinddown_start)-1
    RT_down(i) = (time(spikeinddown_end(i))-time(spikeinddown_start(i))) * 1e3;   % convert to ms
end



% Calculate inter-event intervals
for i = 1:length(spikeindup_start)-1
    IEI_up(i) = (time(spikeindup_start(i+1))-time(spikeindup_end(i))) * 1e3;         % convert to ms
end
for i = 1:length(spikeinddown_start)-1
    IEI_down(i) = (time(spikeinddown_start(i+1))-time(spikeinddown_end(i))) * 1e3;   % convert to ms
end

figure(4);
binsizeup = 2*iqr(RT_up)*length(RT_up)^(-1/3);      % freedman-diaconis rule
binsizedown = 2*iqr(RT_down)*length(RT_down)^(-1/3);
nbinsup = round((max(RT_up) - min(RT_up))/binsizeup);
nbinsdown = round((max(RT_down) - min(RT_down))/binsizedown);
if nbinsup > 1e4                    % maximum limit on bins to reduce memory use
    nbinsup = 1e4;
else if isempty(nbinsup) ==1 | isnan(nbinsup)
        nbinsup=2;
    end    
end
if nbinsdown > 1e4
    nbinsdown = 1e4;
else if isempty(nbinsdown) ==1 | isnan(nbinsdown)
        nbinsdown=2;
    end
end
subplot(2,2,1); 
hist(RT_up,nbinsup); title('residence time - UP');
subplot(2,2,2);
hist(RT_down,nbinsdown); title('residence time - DOWN');
clear binsizeup binsizedown nbinsup nbinsdown
binsizeup = 2*iqr(IEI_up)*length(IEI_up)^(-1/3);      % freedman-diaconis rule
binsizedown = 2*iqr(IEI_down)*length(IEI_down)^(-1/3);
nbinsup = round((max(IEI_up) - min(IEI_up))/binsizeup);
nbinsdown = round((max(IEI_down) - min(IEI_down))/binsizedown);
if nbinsup > 1e4                    % maximum limit on bins to reduce memory use
    nbinsup = 1e4;
else if isempty(nbinsup) ==1 | isnan(nbinsup)
        nbinsup=2;
    end    
end
if nbinsdown > 1e4
    nbinsdown = 1e4;
else if isempty(nbinsdown) ==1 | isnan(nbinsdown)
        nbinsdown=2;
    end
end
subplot(2,2,3); 
hist(IEI_up,nbinsup); title('inter-event interval - UP');
subplot(2,2,4);
hist(IEI_down,nbinsdown); title('inter-event interval - DOWN');
    


display('Find buffer region and recalculate...');
bufferind = upindex(:)+downindex(:);     % Buffer index is the index for all those points not included in the up or down excursions (0 is in buffer, 1 is not)
bufferinddiff = diff(bufferind);
buffer_start = find(bufferinddiff==1);
buffer_end = find(bufferinddiff==-1);

if buffer_start(1) > buffer_end(1)
    buffer_end(1) = [];
end

% Calculate residence times within the buffer
for i = 1:length(buffer_start)-1
    RT_buffer(i) = (time(buffer_end(i))-time(buffer_start(i))) * 1e3;         % convert to ms
end
% Calculate inter-buffer times
for i = 1:length(buffer_start)-1
    IEI_buffer(i) = (time(buffer_start(i+1))-time(buffer_end(i))) * 1e3;         % convert to ms
end
figure(5);
binsizebuffer = 2*iqr(RT_buffer)*length(RT_buffer)^(-1/3);      % freedman-diaconis rule
nbinsbuffer = round((max(RT_buffer) - min(RT_buffer))/binsizebuffer);
subplot(1,2,1); 
hist(RT_buffer,nbinsbuffer); title('residence time - buffer');
binsizebuffer = 2*iqr(IEI_buffer)*length(IEI_buffer)^(-1/3);      % freedman-diaconis rule
nbinsbuffer = round((max(IEI_buffer) - min(IEI_buffer))/binsizebuffer);
if nbinsbuffer > 1e4                    % maximum limit on bins to reduce memory use
    nbinsbuffer = 1e4;
end

subplot(1,2,2); 
hist(IEI_buffer,nbinsbuffer); title('inter-event interval - buffer');

display('Saving...');
%save('/Users/joshsalvi/Documents/Lab/Lab/Simulation Data/ONHFishJosh/Analysis/Fp7-kp2-start1end2e5-xfish1.0Noise-movingavg_concwins2.mat');
save('/Users/joshsalvi/Documents/Lab/Lab/Simulation Data/Sinusoids/noisysinewave22movingavg_concwins2_exp.mat');

display('Finished.');
