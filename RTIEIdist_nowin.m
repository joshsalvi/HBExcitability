clear all; close all;
display('Importing...');
load('/Users/joshsalvi/Documents/Lab/Lab/Simulation Data/ONHFishJosh/Analysis/Fp3-kp3-xfish1.0Noisemovingavg.mat');
%load('/Users/joshsalvi/Documents/Lab/Lab/Simulation Data/Sinusoids/noisysinewave44movingavg.mat');

nup = -3:0.2:3;
ndown = -3:0.2:3;

up_ind = zeros(length(Xvec),length(nup));
down_ind = zeros(length(Xvec),length(ndown));

display('Find up events...');
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
up_ind(:,j) = Mvec_up;
up_events(:,j) = floor(sum(abs(diff(up_ind(:,j))))/2);
end
events_up_diff = diff(up_events);


display('Find down events...');
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
down_ind(:,j) = Mvec_down;
down_events(:,j) = floor(sum(abs(diff(down_ind(:,j))))/2);
end
events_down_diff = diff(down_events);


% Find local minima
uppeaks = findpeaks(up_events);
downpeaks = findpeaks(down_events);
if length(uppeaks)>1
    nupconv = find(up_events(uppeaks(1):uppeaks(2))==min(up_events(uppeaks(1):uppeaks(2))))+uppeaks(1)-1;
else
    nupconv = uppeaks(1);
end
if length(downpeaks)>1
    ndownconv = find(down_events(downpeaks(1):downpeaks(2))==min(down_events(downpeaks(1):downpeaks(2))))+downpeaks(1)-1;
else
    ndownconv = downpeaks(1);
end

if length(nupconv)>1
    nupconv = nupconv(1);
end
if length(ndownconv)>1
    ndownconv = ndownconv(1);
end


figure(2);
subplot(2,2,1);plot(nup,up_events);ylabel('Events');xlabel('nup');title('Up Events');
subplot(2,2,2);plot(ndown,down_events);ylabel('Events');xlabel('ndown');title('Down Events');
subplot(2,2,3);plot(nup(2:end),diff(up_events));ylabel('dEvents');xlabel('nup');title('Up Events, Diff');
subplot(2,2,4);plot(ndown(2:end),diff(down_events));ylabel('dEvents');xlabel('nup');title('Up Events, Diff');
subplot(2,2,1);hold on; scatter(nup(nupconv),up_events(nupconv),'ro');scatter(nup(uppeaks),up_events(uppeaks),'k.');
subplot(2,2,2);hold on; scatter(ndown(ndownconv),down_events(ndownconv),'ro');scatter(ndown(downpeaks),down_events(downpeaks),'k.');


figure(3);
plot(time,Xvec);hold on; ylabel('displacement');
plot(time(down_ind(:,ndownconv)==1),Xvec(down_ind(:,ndownconv)==1),'r.');
plot(time(up_ind(:,nupconv)==1),Xvec(up_ind(:,nupconv)==1),'g.');
title(sprintf('%s %s%s %s%s %s%s%s%s','Spikes ','nup= ',num2str(nup(nupconv)),'ndown= ',num2str(ndown(ndownconv)),'eventsup/eventsdown= ',num2str(up_events(nupconv)),'/',num2str(down_events(ndownconv))));


% Find residence times, buffer times, etc....
display('Find residence times...');
upindex = up_ind(:,nupconv);
downindex = down_ind(:,ndownconv);          % indices for up and down positions
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
%save('/Users/joshsalvi/Documents/Lab/Lab/Simulation Data/ONHFishJosh/Analysis/noisysinewave44RTDIEI_nowin.mat');
save('/Users/joshsalvi/Documents/Lab/Lab/Simulation Data/ONHFishJosh/Analysis/Fp3-kp3-xfish1.0NoiseRTIEI_nowin.mat');


display('Finished.');
