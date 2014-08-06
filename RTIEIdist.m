clear all; close all;
display('Importing...');
load('/Users/joshsalvi/Documents/Lab/Lab/Simulation Data/ONHFishJosh/Analysis/Fp4-kp2-xfish1.0Noisemovingavg_1_1e6.mat');

wins = 1:5:30;                                    % window range (for minimum time between events), pts
nup = -3:0.2:3;
ndown = -3:0.2:3;

up_ind = zeros(length(Xvec),length(nup));
down_ind = zeros(length(Xvec),length(ndown));

display('Find up events...');
for l = 1:length(wins)
for j = 1:length(nup)
    Mvec_up = zeros(1,length(Xvec));
for i = 1:length(movingavgX)
    if i <= win
        clear a;
        a = Xvec(i:round(win/2)+i-1) > movingavgX(i) + nup(j)*movingstdX(i);            % left edge
        Mvec_up(i:round(win/2)+i-1) = Mvec_up(i:round(win/2)+i-1)+a(:)';
    else if i < length(movingavgX) - win
            clear a;
            a = Xvec(i-floor(win/2):round(win/2)+i-1) > movingavgX(i) + nup(j)*movingstdX(i);           % middle
            Mvec_up(i-floor(win/2):round(win/2)+i-1) = Mvec_up(i-floor(win/2):round(win/2)+i-1)+a(:)';
        else
            clear a;
            a = Xvec(i-floor(win/2):i) > movingavgX(i) + nup(j)*movingstdX(i);          % right edge
            Mvec_up(i-floor(win/2):i) = Mvec_up(i-floor(win/2):i)+a(:)';
        end
    end
end
Mvec_up(Mvec_up > 0) = 1;   % reduce multiple counts to 1
up_ind(:,j,l) = Mvec_up;
up_ind(:,j,l) = ceil(smooth(up_ind(:,j,l),wins(l)));        % remove points below window threshold (i.e., concatenate spikes)
up_events(:,j,l) = floor(sum(abs(diff(up_ind(:,j,l))))/2);
end
if wins(l) == wins(round(length(wins)/2))
    display('50%');
end
end
events_up_diff = diff(up_events);


display('Find down events...');
for l = 1:length(wins)
for j = 1:length(ndown)
    Mvec_down = zeros(1,length(Xvec));
for i = 1:length(movingavgX)
    if i <= win
        clear a;
        a = Xvec(i:round(win/2)+i-1) < movingavgX(i) - ndown(j)*movingstdX(i);            % left edge
        Mvec_down(i:round(win/2)+i-1) = Mvec_down(i:round(win/2)+i-1)+a(:)';
    else if i < length(movingavgX) - win
            clear a;
            a = Xvec(i-floor(win/2):round(win/2)+i-1) < movingavgX(i) - ndown(j)*movingstdX(i);           % middle
            Mvec_down(i-floor(win/2):round(win/2)+i-1) = Mvec_down(i-floor(win/2):round(win/2)+i-1)+a(:)';
        else
            clear a;
            a = Xvec(i-floor(win/2):i) < movingavgX(i) - ndown(j)*movingstdX(i);          % right edge
            Mvec_down(i-floor(win/2):i) = Mvec_down(i-floor(win/2):i)+a(:)';
        end
    end
end
Mvec_down(Mvec_down > 0) = 1;   % reduce multiple counts to 1
down_ind(:,j,l) = Mvec_down;
down_ind(:,j,l) = ceil(smooth(down_ind(:,j,l),wins(l)));        % remove points below window threshold (i.e., concatenate spikes)
down_events(:,j,l) = floor(sum(abs(diff(down_ind(:,j,l))))/2);
end
if wins(l) == wins(round(length(wins)/2))
    display('50%');
end
end
events_down_diff = diff(down_events);


for l = 1:length(wins)
uppeaks{l} = findpeaks(up_events(:,:,l));
downpeaks{l} = findpeaks(down_events(:,:,l));
nupconv{l} = uppeaks{l}(end);
ndownconv{l} = downpeaks{l}(end);

if length(nupconv{l})>1
    nupconv{l} = nupconv{l}(1);
end
if length(ndownconv{l})>1
    ndownconv{l} = ndownconv{l}(1);
end
end
nupconv = cell2num(nupconv);
ndownconv = cell2num(ndownconv);

uppeakswins = findpeaks((diff(diag(squeeze(up_events(1,nupconv(),:))'))))-1;
downpeakswins = findpeaks((diff(diag(squeeze(down_events(1,ndownconv(),:))'))))-1;
%uppeakswins = find(((diag(squeeze(up_events(1,nupconv(),:))')))==min(((diag(squeeze(up_events(1,nupconv(),:))')))));
%downpeakswins = find(((diag(squeeze(down_events(1,ndownconv(),:))')))==min(((diag(squeeze(down_events(1,ndownconv(),:))')))));
%uppeakswins = find(sign(diff(diag(squeeze(up_events(1,nupconv(),:))')))==0);
%downpeakswins = find(sign(diff(diag(squeeze(down_events(1,nupconv(),:))')))==0);
nupconvwins = uppeakswins(1);
ndownconvwins = downpeakswins(1);

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
end
if nbinsdown > 1e4
    nbinsdown = 1e4;
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
end
if nbinsdown > 1e4
    nbinsdown = 1e4;
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
save('/Users/joshsalvi/Documents/Lab/Lab/Simulation Data/ONHFishJosh/Analysis/Fp4-kp2-xfish1.0NoiseRTIEI_1_1e6.mat');

display('Finished.');
