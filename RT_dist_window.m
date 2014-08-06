% This function uses a moving filter as opposed to detrending the data

clear all; close all;

display('Importing...');

%load('/Users/joshsalvi/Downloads/4.mat');
%load('/Users/joshsalvi/Documents/Lab/Lab/Simulation Data/ONHFishJosh/xfish1.0Noise.mat');
load('/Users/joshsalvi/Documents/Lab/Lab/Simulation Data/ONHFishJosh/noisysinewave.mat');

Fs = 10e3;                          % Choose the sampling frequency (note that this won't work properly for simulation data)
% Select operating points

% FORCE AND STIFFNESS


Ord_k = unique(k_rand);
Ord_F = unique(F_rand);

Fp = 1;                             % SELECT OPERATING POINT
kp = 1;

[Np,Nt] = ind2sub(size(F_rand),find(k_rand==Ord_k(kp) & F_rand==Ord_F(Fp)));

startl = 1;                      % how much time would you like to use?
endl = 1e4;
Xvec = Xd(startl:endl,Np,Nt);
clear Xd;

% FORCE ONLY
%{
Fi = 3;
ind = Fi;
Xvec = Xd(:,ind);
clear Xd;
%}
display('Moving average and standard deviation...');
time = linspace(0,length(Xvec)/Fs,length(Xvec));    % time in seconds

figure(1); subplot(3,1,1); plot(time,Xvec); title(sprintf('%s%s %s%s','F= ',num2str(F_rand(Np,Nt)),' k= ',num2str(k_rand(Np,Nt))));ylabel('displacement')
% Moving window, running average using filter
win = 2e3+1;                                        % CHOOSE WINDOW, must be ODD
wins = 2:10:200;                                    % CHOOSE SMOOTHING WINDOW RANGE (if used)
movingavgX = filter(ones(1,win)/win,1,Xvec);
subplot(3,1,2); plot(time,movingavgX); title(sprintf('%s %s','window = ',num2str(win))); ylabel('moving average');
movingstdX = movingstd(Xvec,win,'central');
subplot(3,1,3); plot(time,movingstdX); title(sprintf('%s %s','window = ',num2str(win))); ylabel('moving standard deviation');


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
end
events_down_diff = diff(down_events);


%{
figure(2);
subplot(2,2,1);plot(nup,up_events);title('Up Events'); ylabel('Events'); xlabel('nup');
subplot(2,2,2);plot(ndown,down_events);title('Down Events'); ylabel('Events'); xlabel('ndown');
subplot(2,2,3);plot(nup(2:end),events_up_diff); ylabel('dEvents'); xlabel('nup');
subplot(2,2,4);plot(ndown(2:end),events_down_diff); ylabel('dEvents'); xlabel('nup');
%}

% METHOD A - Find local maxima
%{
nupconv = find(diff(sign(events_up_diff))==-2)+1;
ndownconv = find(diff(sign(events_down_diff))==-2)+1;
if length(nupconv) > 1                                      % in case there is more than one point found
    nupconv = find(up_events == max(up_events(nupconv)));
end
if length(ndownconv) > 1
    ndownconv = find(down_events == max(down_events(ndownconv)));
end
if length(nupconv) > 1                                      % if there's still a problem
    nupconv = nupconv(1);
end
if length(ndownconv) > 1
    ndownconv = ndownconv(1);
end
figure(2);
subplot(2,2,1);hold on; scatter(nup(nupconv),up_events(nupconv),'ro');
subplot(2,2,2);hold on; scatter(ndown(ndownconv),down_events(ndownconv),'ro');
%}
        
% METHOD B - Find local minima, method 1
%{
display('Find local minima...');
nupconv = find(diff(sign(events_up_diff))==2)+1;
ndownconv = find(diff(sign(events_down_diff))==2)+1;
if isempty(nupconv) == 1                                % If no local minimum can be found, use a local maximum
    nupconv = find(diff(sign(events_up_diff))==-2)+1;
end
if isempty(ndownconv) == 1                                % If no local minimum can be found, use a local maximum
    ndownconv = find(diff(sign(events_down_diff))==-2)+1;
end

if length(nupconv) > 1                                      % in case there is more than one point found
    nupconv = find(up_events == min(up_events(nupconv)));
end
if length(ndownconv) > 1
    ndownconv = find(down_events == min(down_events(ndownconv)));
end
if length(nupconv) > 1                                      % if there's still a problem
    nupconv = nupconv(1);
end
if length(ndownconv) > 1
    ndownconv = ndownconv(1);
end
figure(2);
subplot(2,2,1);hold on; scatter(nup(nupconv),up_events(nupconv),'ro');
subplot(2,2,2);hold on; scatter(ndown(ndownconv),down_events(ndownconv),'ro');
%}


% METHOD C - Find local minima, method 2
%{
updiffmin = find(events_up_diff==min(events_up_diff));
updiffmax = find(events_up_diff==max(events_up_diff));
downdiffmin = find(events_down_diff==min(events_down_diff));
downdiffmax = find(events_down_diff==max(events_down_diff));

if updiffmin<updiffmax
    nupconv = (find(up_events(updiffmin:updiffmax)==min(up_events(updiffmin:updiffmax)))) + updiffmin - 1;
else
    nupconv = (find(up_events(updiffmax:updiffmin)==min(up_events(updiffmax:updiffmin)))) + updiffmax - 1;
end

if downdiffmin<downdiffmax
    ndownconv = (find(down_events(downdiffmin:downdiffmax)==min(down_events(downdiffmin:downdiffmax)))) + downdiffmin - 1;
else
    ndownconv = (find(down_events(downdiffmax:downdiffmin)==min(down_events(downdiffmax:downdiffmin)))) + downdiffmax - 1;
end

if length(nupconv) > 1                                      % in case there is more than one point found
    nupconv = find(up_events == min(up_events(nupconv)));
end
if length(ndownconv) > 1
    ndownconv = find(down_events == min(down_events(ndownconv)));
end

figure(2);
subplot(2,2,1);hold on; scatter(nup(nupconv),up_events(nupconv),'ro');scatter(nup(updiffmin),up_events(updiffmin),'k.');scatter(nup(updiffmax),up_events(updiffmax),'k.');
subplot(2,2,2);hold on; scatter(ndown(ndownconv),down_events(ndownconv),'ro');scatter(ndown(downdiffmin),down_events(downdiffmin),'k.');scatter(ndown(downdiffmax),down_events(downdiffmax),'k.');
subplot(2,2,3);hold on; scatter(nup(nupconv),events_up_diff(nupconv),'ro');scatter(nup(updiffmin),events_up_diff(updiffmin),'k.');scatter(nup(updiffmax),events_up_diff(updiffmax),'k.');
subplot(2,2,4);hold on; scatter(ndown(ndownconv),events_down_diff(ndownconv),'ro');scatter(ndown(downdiffmin),events_down_diff(downdiffmin),'k.');scatter(ndown(downdiffmax),events_down_diff(downdiffmax),'k.');
%}

% METHOD D - Find local minima, method 3
%{
uppeaks = findpeaks(up_events);
downpeaks = findpeaks(down_events);
nupconv = find(up_events(uppeaks(1):uppeaks(2))==min(up_events(uppeaks(1):uppeaks(2))))+uppeaks(1)-1;
ndownconv = find(down_events(downpeaks(1):downpeaks(2))==min(down_events(downpeaks(1):downpeaks(2))))+downpeaks(1)-1;

if length(nupconv)>1
    nupconv = nupconv(1);
end
if length(ndownconv)>1
    ndownconv = ndownconv(1);
end


figure(2);
subplot(2,2,1);hold on; scatter(nup(nupconv),up_events(nupconv),'ro');scatter(nup(uppeaks),up_events(uppeaks),'k.');
subplot(2,2,2);hold on; scatter(ndown(ndownconv),down_events(ndownconv),'ro');scatter(ndown(downpeaks),down_events(downpeaks),'k.');
%}


% METHOD E - Find local maxima, method 2

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

%uppeakswins = findpeaks((diff(diag(squeeze(up_events(1,nupconv(),:))'))));
%downpeakswins = findpeaks((diff(diag(squeeze(down_events(1,ndownconv(),:))'))));
uppeakswins = find(((diag(squeeze(up_events(1,nupconv(),:))')))==min(((diag(squeeze(up_events(1,nupconv(),:))')))));
downpeakswins = find(((diag(squeeze(down_events(1,ndownconv(),:))')))==min(((diag(squeeze(down_events(1,ndownconv(),:))')))));
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
if nbinsup > 1e4                    % maximum limit on bins to reduce memory
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
if nbinsup > 1e4                    % maximum limit on bins to reduce memory
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
if nbinsbuffer > 1e4                    % maximum limit on bins to reduce memory
    nbinsbuffer = 1e4;
end

subplot(1,2,2); 
hist(IEI_buffer,nbinsbuffer); title('inter-event interval - buffer');


display('Finished.');

%% SUPPLEMENTAL SCRIPTS (not in use)

% Concatenate residence times, method 2
display('Concatenate residence times...')
winrt =6;                            % CHOOSE MINIMUM BETWEEN EVENTS (ms)
for i = 1:length(RT_up)
    if i > 1
    if RT_up(i) < winrt
        RT_up_conc(i) = RT_up(i)+RT_up_conc(i-1);    % convert to ms
        RT_up_conc(i-1) = 0;
    else
        RT_up_conc(i) = RT_up(i);
    end
    else
        RT_up_conc(i) = RT_up(i);
    end
end
RT_up_conc(RT_up_conc==0)=[];

for i = 1:length(RT_down)
    if i > 1
    if RT_down(i) < winrt
        RT_down_conc(i) = RT_down(i)+RT_down_conc(i-1);    % convert to ms
        RT_down_conc(i-1) = 0;
    else
        RT_down_conc(i) = RT_down(i);
    end
    else
        RT_down_conc(i) = RT_down(i);
    end
end
RT_down_conc(RT_down_conc==0)=[];

display('Concatenate inter-event interval times...')
                            % CHOOSE MINIMUM BETWEEN EVENTS (ms)
for i = 1:length(IEI_up)
    if i > 1
    if IEI_up(i) < winrt
        IEI_up_conc(i) = IEI_up(i)+IEI_up_conc(i-1);    % convert to ms
        IEI_up_conc(i-1) = 0;
    else
        IEI_up_conc(i) = IEI_up(i);
    end
    else
        IEI_up_conc(i) = IEI_up(i);
    end
end
IEI_up_conc(IEI_up_conc==0)=[];

for i = 1:length(IEI_down)
    if i > 1
    if IEI_down(i) < winrt
        IEI_down_conc(i) = IEI_down(i)+IEI_down_conc(i-1);    % convert to ms
        IEI_down_conc(i-1) = 0;
    else
        IEI_down_conc(i) = IEI_down(i);
    end
    else
        IEI_down_conc(i) = IEI_down(i);
    end
end
IEI_down_conc(IEI_down_conc==0)=[];


figure(6);
binsizeup = 2*iqr(RT_up_conc)*length(RT_up_conc)^(-1/3);      % freedman-diaconis rule
binsizedown = 2*iqr(RT_down_conc)*length(RT_down_conc)^(-1/3);
nbinsup = round((max(RT_up_conc) - min(RT_up_conc))/binsizeup);
nbinsdown = round((max(RT_down_conc) - min(RT_down_conc))/binsizedown);
subplot(2,2,1); 
hist(RT_up_conc,nbinsup); title(sprintf('%s %s%s%s%s','residence time - UP, concatenated','countup/countdown= ',num2str(length(RT_up_conc)),', ',num2str(length(RT_down_conc))));
subplot(2,2,2);
hist(RT_down_conc,nbinsdown); title('residence time - DOWN, concatenated');
clear binsizeup binsizedown nbinsup nbinsdown
binsizeup = 2*iqr(IEI_up_conc)*length(IEI_up_conc)^(-1/3);      % freedman-diaconis rule
binsizedown = 2*iqr(IEI_down_conc)*length(IEI_down)^(-1/3);
nbinsup = round((max(IEI_up_conc) - min(IEI_up_conc))/binsizeup);
nbinsdown = round((max(IEI_down_conc) - min(IEI_down_conc))/binsizedown);
subplot(2,2,3); 
hist(IEI_up_conc,nbinsup); title(sprintf('%s %s%s%s%s','inter-event interval - UP, concatenated','countup/countdown= ',num2str(length(IEI_up_conc)),', ',num2str(length(IEI_down_conc))));
subplot(2,2,4);
hist(IEI_down_conc,nbinsdown); title('inter-event interval - DOWN, concatenated');
