clear all; close all;

display('Importing...');

%load('/Users/joshsalvi/Downloads/4.mat');
load('/Users/joshsalvi/Documents/Lab/Lab/Simulation Data/ONHFishJosh/xfish1.0Noise.mat');

Fs = 10e3;                          % Choose the sampling frequency (note that this won't work properly for simulation data)
% Select operating points

% FORCE AND STIFFNESS

Np = 4;
Nt = 4;
Xvec = Xd(:,Np,Nt);
clear Xd;

% FORCE ONLY
%{
Fi = 3;
ind = Fi;
Xvec = Xd(:,ind);
clear Xd;
%}

time = linspace(0,length(Xvec)/Fs,length(Xvec));    % time in seconds

figure(1); subplot(2,1,1); plot(time,Xvec); title(sprintf('%s%s %s%s','F= ',num2str(F_rand(Np,Nt)),' k= ',num2str(k_rand(Np,Nt))));
% Detrend data
win = 8e3;          % CHOOSE WINDOW
Xvec = Xvec - smooth(Xvec,win);
subplot(2,1,2); plot(time,Xvec); title(sprintf('%s %s','window = ',num2str(win)));
clear win;


display('Convergence of standard deviations...');
% Find number of events

% Convergence for nup and ndown
nup = 0:0.05:5;            % Initialize nup and ndown
ndown = 0:0.05:5;
wintime = 0.001;             % CHOOSE WINDOW TIME IN SECONDS
win = findnearest(time,wintime); 

% OPTION 1 - No stepping
mux = mean(Xvec);
stdx = std(Xvec);
display('UP');
for i = 1:length(nup);
    Xvecind_up(:,i) = Xvec >= mux + nup(i)*stdx;
    Xvecind_up(:,i) = ceil(smooth(Xvecind_up(:,i),win));        % remove points below window threshold (i.e., concatenate spikes)
    events_up(i) = floor(sum(abs(diff(Xvecind_up(:,i))))/2);
end
events_up_diff = diff(events_up);
display('DOWN');
for j = 1:length(ndown)
    Xvecind_down(:,j) = Xvec <= mux - ndown(j)*stdx;
    Xvecind_down(:,j) = ceil(smooth(Xvecind_down(:,j),win));      % remove points below window threshold (i.e., concatenate spikes)
    events_down(j) = floor(sum(abs(diff(Xvecind_down(:,j))))/2);
end 
events_down_diff = diff(events_down);

figure(2);
subplot(2,2,1);plot(nup,events_up);title('Up Events'); ylabel('Events'); xlabel('nup');
subplot(2,2,2);plot(ndown,events_down);title('Down Events'); ylabel('Events'); xlabel('ndown');
subplot(2,2,3);plot(nup(2:end),events_up_diff); ylabel('dEvents'); xlabel('nup');
subplot(2,2,4);plot(ndown(2:end),events_down_diff); ylabel('dEvents'); xlabel('nup');

% METHOD A - Find local maxima
%{
nupconv = find(diff(sign(events_up_diff))==-2)+1;
ndownconv = find(diff(sign(events_down_diff))==-2)+1;
if length(nupconv) > 1                                      % in case there is more than one point found
    nupconv = find(events_up == min(events_up(nupconv)));
end
if length(ndownconv) > 1
    ndownconv = find(events_down == min(events_down(ndownconv)));
end
figure(2);
subplot(2,2,1);hold on; scatter(nup(nupconv),events_up(nupconv),'ro');
subplot(2,2,2);hold on; scatter(ndown(ndownconv),events_down(ndownconv),'ro');
%}
    
% METHOD B - Find local minima, method 1
nupconv = find(diff(sign(events_up_diff))==2)+1;
ndownconv = find(diff(sign(events_down_diff))==2)+1;
if isempty(nupconv) == 1                                % If no local minimum can be found, use a local maximum
    nupconv = find(diff(sign(events_up_diff))==-2)+1;
end
if isempty(ndownconv) == 1                                % If no local minimum can be found, use a local maximum
    ndownconv = find(diff(sign(events_down_diff))==-2)+1;
end

if length(nupconv) > 1                                      % in case there is more than one point found
    nupconv = find(events_up == min(events_up(nupconv)));
end
if length(ndownconv) > 1
    ndownconv = find(events_down == min(events_down(ndownconv)));
end

figure(2);
subplot(2,2,1);hold on; scatter(nup(nupconv),events_up(nupconv),'ro');
subplot(2,2,2);hold on; scatter(ndown(ndownconv),events_down(ndownconv),'ro');


% METHOD C - Find local minima, method 2
%{
updiffmin = find(events_up_diff==min(events_up_diff));
updiffmax = find(events_up_diff==max(events_up_diff));
downdiffmin = find(events_down_diff==min(events_down_diff));
downdiffmax = find(events_down_diff==max(events_down_diff));

if updiffmin<updiffmax
    nupconv = (find(events_up(updiffmin:updiffmax)==min(events_up(updiffmin:updiffmax)))) + updiffmin - 1;
else
    nupconv = (find(events_up(updiffmax:updiffmin)==min(events_up(updiffmax:updiffmin)))) + updiffmax - 1;
end

if downdiffmin<downdiffmax
    ndownconv = (find(events_down(downdiffmin:downdiffmax)==min(events_down(downdiffmin:downdiffmax)))) + downdiffmin - 1;
else
    ndownconv = (find(events_down(downdiffmax:downdiffmin)==min(events_down(downdiffmax:downdiffmin)))) + downdiffmax - 1;
end

if length(nupconv) > 1                                      % in case there is more than one point found
    nupconv = find(events_up == min(events_up(nupconv)));
end
if length(ndownconv) > 1
    ndownconv = find(events_down == min(events_down(ndownconv)));
end

figure(2);
subplot(2,2,1);hold on; scatter(nup(nupconv),events_up(nupconv),'ro');scatter(nup(updiffmin),events_up(updiffmin),'k.');scatter(nup(updiffmax),events_up(updiffmax),'k.');
subplot(2,2,2);hold on; scatter(ndown(ndownconv),events_down(ndownconv),'ro');scatter(ndown(downdiffmin),events_down(downdiffmin),'k.');scatter(ndown(downdiffmax),events_down(downdiffmax),'k.');
subplot(2,2,3);hold on; scatter(nup(nupconv),events_up_diff(nupconv),'ro');scatter(nup(updiffmin),events_up_diff(updiffmin),'k.');scatter(nup(updiffmax),events_up_diff(updiffmax),'k.');
subplot(2,2,4);hold on; scatter(ndown(ndownconv),events_down_diff(ndownconv),'ro');scatter(ndown(downdiffmin),events_down_diff(downdiffmin),'k.');scatter(ndown(downdiffmax),events_down_diff(downdiffmax),'k.');
%}


figure(3);
plot(time,Xvec);hold on; 
plot(time(Xvecind_down(:,ndownconv)),Xvec(Xvecind_down(:,ndownconv)),'r.');
plot(time(Xvecind_up(:,nupconv)),Xvec(Xvecind_up(:,nupconv)),'g.');
title(sprintf('%s %s%s %s%s %s%s%s%s','Spikes ','nup= ',num2str(nup(nupconv)),'ndown= ',num2str(ndown(ndownconv)),'eventsup/eventsdown= ',num2str(events_up(nupconv)),'/',num2str(events_down(ndownconv))));


% OPTION 2 - Moving window
%{
wintime = 10;     % CHOOSE WINDOW TIME IN SECONDS
win = findnearest(time,wintime);
overlap = 0;
for i = 1:length(Xvec)-win
    clear Xvec2
    Xvec2 = Xvec(i:i+win-1);
    mux = mean(Xvec2);
    stdx = std(Xvec2);
    
    for j = 1:length(nup)
        Xvec2ind_up(:,i,j) = Xvec2 > mux + nup(j)*stdx;
    end
    
    for k = 1:length(ndown)
        Xvec2ind_down(:,i,j) = Xvec2 > mux - ndown(j)*stdx;
    end
end
%}
    
display('Find residence times...');
% Find the residence time distribution
upindex = Xvecind_up(:,nupconv);
downindex = Xvecind_down(:,ndownconv);          % indices for up and down positions
upinddiff = diff(upindex);
downinddiff = diff(downindex);                  % first derivative (+1 is start of excursion, -1 is end of excursion)
spikeindup_start = find(upinddiff==1);
spikeindup_end = find(upinddiff==-1);
%{
if spikeindup_start(1) > spikeindup_end(1)
    spikeindup_start(1) == [];
end
%}
spikeinddown_start = find(downinddiff==1);
spikeinddown_end = find(downinddiff==-1);
%{
if spikeinddown_start(1) > spikeinddown_end(1)
    spikeinddown_start(1) == [];
end
%}

% Calculate residence times for the up position and for the down position
for i = 1:length(spikeindup_start)
    RT_up(i) = (time(spikeindup_end(i))-time(spikeindup_start(i))) * 1e3;         % convert to ms
end
for i = 1:length(spikeinddown_start)
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
subplot(2,2,1); 
hist(RT_up,nbinsup); title('residence time - UP');
subplot(2,2,2);
hist(RT_down,nbinsdown); title('residence time - DOWN');
clear binsizeup binsizedown nbinsup nbinsdown
binsizeup = 2*iqr(IEI_up)*length(IEI_up)^(-1/3);      % freedman-diaconis rule
binsizedown = 2*iqr(IEI_down)*length(IEI_down)^(-1/3);
nbinsup = round((max(IEI_up) - min(IEI_up))/binsizeup);
nbinsdown = round((max(IEI_down) - min(IEI_down))/binsizedown);
subplot(2,2,3); 
hist(RT_up,nbinsup); title('inter-event interval - UP');
subplot(2,2,4);
hist(RT_down,nbinsdown); title('inter-event interval - DOWN');
    

display('Find buffer region and recalculate...');
bufferind = upindex(:)+downindex(:);     % Buffer index is the index for all those points not included in the up or down excursions (0 is in buffer, 1 is not)
bufferinddiff = diff(bufferind);
buffer_start = find(bufferinddiff==1);
buffer_end = find(bufferinddiff==-1);

% Calculate residence times within the buffer
for i = 1:length(spikeindup_start)
    RT_buffer(i) = (time(buffer_end(i))-time(buffer_start(i))) * 1e3;         % convert to ms
end
% Calculate inter-buffer times
for i = 1:length(spikeindup_start)-1
    IEI_buffer(i) = (time(buffer_start(i+1))-time(buffer_end(i))) * 1e3;         % convert to ms
end
figure(5);
binsizebuffer = 2*iqr(RT_buffer)*length(RT_buffer)^(-1/3);      % freedman-diaconis rule
nbinsbuffer = round((max(RT_buffer) - min(RT_buffer))/binsizebuffer);
subplot(1,2,1); 
hist(RT_buffer,nbinsbuffer); title('residence time - buffer');
binsizebuffer = 2*iqr(IEI_buffer)*length(IEI_buffer)^(-1/3);      % freedman-diaconis rule
nbinsbuffer = round((max(IEI_buffer) - min(IEI_buffer))/binsizebuffer);
subplot(1,2,2); 
hist(IEI_buffer,nbinsbuffer); title('inter-event interval - buffer');
