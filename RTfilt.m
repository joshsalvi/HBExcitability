function [uptimes, downtimes, Xvecfilt_detrended,RTupexp_p,RTupgauss_p,RTdownexp_p,RTdowngauss_p,RTupexp_stat,RTdownexp_stat,RTupgauss_stat,RTdowngauss_stat] = RTfilt(dir,Fp,kp,win1,freq1,Fs,sopt,out,mintime,altsol,figures)
% Residence Time Calculations
% This function defines the residence times for a signal by first finding
% the moving average, then applies a high-pass filter at an input
% frequency. The function will output plots of the original time signal,
% moving average, moving standard deviation, and filtered/detrended signal
% for analysis. A second plot then calculates the histogram from the
% Freedman-Diaconis rule and calculates the statistics for both exponential
% and normal distributions.
%
% [uptimes downtimes, Xvecfilt_detrended,RTupexp_p,RTupgauss_p,RTdownexp_p,RTdowngauss_p,RTupexp_stat,RTdownexp_stat,RTupgauss_stat,RTdowngauss_stat] = RTfilt(dir,Fp,kp,win1,freq1,Fs,sopt,out,mintime,altsol,figures)
%
%  dir  :    directory of time signal
%  Fp   :    force index
%  kp   :    stiffness index
%  win1 :    window for moving average and standard deviation (in sec)
%  freq1:    frequency of high-pass filter (in Hz)
%  Fs   :    scan rate (Hz)
%  sopt :    save signal? (1=yes, 0=no)
%  out  :    saved output file with all of the parameters
%  mintime : minimum time allowed for residence times
%  altsol  : run an alternative solution? (1=yes, 0=no)
%  figures : display figures? (1=yes, 0=no)
%
%      out not required if sopt==0; figures,mintime not required,
%
%     RTfilt('/dir/sine.mat',1,1,3,800,1,'/dir/sineout.mat');
%

warning off;
%display('Importing...');

%load('/Users/joshsalvi/Downloads/4.mat');
%load('/Users/joshsalvi/Documents/Lab/Lab/Simulation Data/ONHFishJosh/xfish1.0Noise.mat');
%load('/Users/joshsalvi/Documents/Lab/Lab/Simulation Data/Sinusoids/noisysinewave.mat');
if isempty(dir) == 0
load(dir);
end

%Fs = 10e3;                          % Choose the sampling frequency (note that this won't work properly for simulation data)
% Select operating points

% FORCE AND STIFFNESS


Ord_k = unique(k_rand);
Ord_F = unique(F_rand);

%Fp = 1;                             % SELECT OPERATING POINT
%kp = 1;

[Np,Nt] = ind2sub(size(F_rand),find(k_rand==Ord_k(kp) & F_rand==Ord_F(Fp)));

startl = 1;                      % how much time would you like to use?
endl = length(Xd);
Xvec = Xd(startl:endl,Np,Nt);

Xvec = Xvec(1:10:end);             % downsample (optional)

clear Xd;

% FORCE ONLY
%{
Fi = 3;
ind = Fi;
Xvec = Xd(:,ind);
clear Xd;
%}
%display('Moving average and standard deviation...');
time = linspace(0,length(Xvec)/Fs,length(Xvec));    % time in seconds

if figures==1
figure; subplot(4,1,1); plot(time,Xvec); title(sprintf('%s%s %s%s','F= ',num2str(F_rand(Np,Nt)),' k= ',num2str(k_rand(Np,Nt))));ylabel('%displacement')
end
% Moving window, running average using filter
win = round(win1*Fs);                                     % CHOOSE WINDOW, must be ODD
if rem(floor(win),2) == 1 %odd
    win = floor(win);
else %even
    win = floor(win)+1; 
end
filtfreq2 = freq1;                                 % CHOOSE FILTER WINDOW FOR HIGH-PASS FILTERING (in Hz)
win2 = round(Fs/filtfreq2);

movingavgX = filter(ones(1,win)/win,1,Xvec);
Xvecfilt = medfilt1(Xvec,win2);
% Use only last 80% of signal
Xvecfilt_detrended = Xvecfilt(round(length(Xvecfilt)*0.2):end) - movingavgX(round(length(Xvecfilt)*0.2):end);
time = time(round(length(Xvecfilt)*0.2:end));

if figures==1
    subplot(4,1,1);subplot(4,1,1); hold on;plot(time,Xvec(round(length(Xvecfilt)*0.2:end)),'r'); title(sprintf('%s%s %s%s','F= ',num2str(F_rand(Np,Nt)),' k= ',num2str(k_rand(Np,Nt))));ylabel('%displacement')
subplot(4,1,2); plot(time,movingavgX(round(length(Xvecfilt)*0.2:end))); title(sprintf('%s %s','window = ',num2str(win))); ylabel('moving average');
movingstdX = movingstd(Xvec,win,'central');
subplot(4,1,3); plot(time,movingstdX(round(length(Xvecfilt)*0.2:end))); title(sprintf('%s %s','window = ',num2str(win))); ylabel('moving standard deviation');
subplot(4,1,4); plot(time,Xvecfilt_detrended(1:length(time)),'r'); ylabel('filtered time signal'); title('filtered and detrended signal');
end

% Find up and down positions
Mup = zeros(1,length(Xvecfilt_detrended));          % initialize variables
Mdown = zeros(1,length(Xvecfilt_detrended));
sup = Xvecfilt_detrended > 0;                 % find up/down positions, those that are >0 or <0
sdown = Xvecfilt_detrended < 0;
Mup(sup) = 1;
Mdown(sdown) = 1;

updiffend = find(diff(Mup)==-1);           % find up times
updiffstart = find(diff(Mup)==1);

if isempty(updiffend)==0 && isempty(updiffstart)==0
if updiffend(1)<updiffstart(1)
    updiffend(1)=[];
end
for k = 2:length(updiffend)
        if updiffend(k)>updiffstart(k-1)
            Mup(updiffstart(k-1):updiffend(k)) = 1;
        else
            Mup(updiffend(k-1):updiffstart(k)) = 1;
        end
end
end
up_events = min([length(updiffstart) length(updiffend)]);



for i = 1:up_events
    uptimes(i) = time(updiffend(i)) - time(updiffstart(i));
end

downdiffend = find(diff(Mdown)==-1);         % find down times
downdiffstart = find(diff(Mdown)==1);

if isempty(downdiffend)==0 && isempty(downdiffstart)==0
if downdiffend(1)<downdiffstart(1)
    downdiffend(1)=[];
end
for k = 2:length(downdiffend)
        if downdiffend(k)>downdiffstart(k-1)
            Mdown(downdiffstart(k-1):downdiffend(k)) = 1;
        else
            Mdown(downdiffend(k-1):downdiffstart(k)) = 1;
        end
end
end
down_events = min([length(downdiffstart) length(downdiffend)]);



for i = 1:down_events
    downtimes(i) = time(downdiffend(i)) - time(downdiffstart(i));
end

binsizeup = 2*iqr(uptimes)*length(uptimes)^(-1/3);      % freedman-diaconis rule
binsizedown = 2*iqr(downtimes)*length(downtimes)^(-1/3);
nbinsup = round((max(uptimes) - min(uptimes))/binsizeup);
nbinsdown = round((max(downtimes) - min(downtimes))/binsizedown);
if nbinsup > 1e4                    % maximum limit on bins to reduce memory use
    nbinsup = 1e4;
else if isempty(nbinsup) ==1 || isnan(nbinsup)
        nbinsup=2;
    else if nbinsup < 4
            nbinsup = 4;
        end
    end    
end
if nbinsdown > 1e4
    nbinsdown = 1e4;
else if isempty(nbinsdown) ==1 || isnan(nbinsdown)
        nbinsdown=2;
    else if nbinsdown < 4
            nbinsdown = 4;
        end
    end
end

if length(uptimes) >=4 && length(downtimes)>=4
%disp('Running statistics...');
% Run statistics
[~, RTupexp_p, RTupexp_stat] = lillietest(uptimes,'Distr','exp');
[~, RTupgauss_p, RTupgauss_stat] = lillietest(uptimes);
[~, RTdownexp_p, RTdownexp_stat] = lillietest(downtimes,'Distr','exp');
[~, RTdowngauss_p, RTdowngauss_stat] = lillietest(downtimes);
[~, RTupgauss_p(2), RTupgauss_stat(2)] = kstest(uptimes);
[~, RTdowngauss_p(2), RTdowngauss_stat(2)] = kstest(downtimes);
[~, RTupgauss_p(3), RTupgauss_stat(3)] = jbtest(uptimes);
[~, RTdowngauss_p(3), RTdowngauss_stat(3)] = jbtest(downtimes);

if figures==1
figure;
subplot(2,2,1);histfit(uptimes,nbinsup); title('RT up');
subplot(2,2,2);set(gca, 'visible', 'off');text(0,0.6,sprintf('Mean = %.5f\nStd. Dev. = %.5f\nMedian = %.5f\nMin = %.5f\nMax = %.5f\nLillie(exp) = %.1e\nLillie(gauss) = %.1e\nKS = %.1e\nJB = %.1e\nNumber of counts = %s',mean(uptimes),std(uptimes),median(uptimes),min(uptimes),max(uptimes),RTupexp_stat,RTupgauss_stat(1),RTupgauss_stat(2),RTupgauss_stat(3),num2str(up_events)));
subplot(2,2,3);histfit(downtimes,nbinsdown); title('RT down');
subplot(2,2,4);set(gca, 'visible', 'off');text(0,0.6,sprintf('Mean = %.5f\nStd. Dev. = %.5f\nMedian = %.5f\nMin = %.5f\nMax = %.5f\nLillie(exp) = %.1e\nLillie(gauss) = %.1e\nKS = %.1e\nJB = %.1e\nNumber of counts = %s',mean(downtimes),std(downtimes),median(downtimes),min(downtimes),max(downtimes),RTdownexp_stat,RTdowngauss_stat(1),RTdowngauss_stat(2),RTdowngauss_stat(3),num2str(down_events)));
end
end

if isempty(mintime) == 0
uptimes(uptimes<mintime) = [];
downtimes(downtimes<mintime) = [];
up_events = length(uptimes);
down_events = length(downtimes);

% Run statistics
if length(uptimes) >=4 && length(downtimes)>=4
[~, RTupexp_p, RTupexp_stat] = lillietest(uptimes,'Distr','exp');
[~, RTupgauss_p, RTupgauss_stat] = lillietest(uptimes);
[~, RTdownexp_p, RTdownexp_stat] = lillietest(downtimes,'Distr','exp');
[~, RTdowngauss_p, RTdowngauss_stat] = lillietest(downtimes);
[~, RTupgauss_p(2), RTupgauss_stat(2)] = kstest(uptimes);
[~, RTdowngauss_p(2), RTdowngauss_stat(2)] = kstest(downtimes);
[~, RTupgauss_p(3), RTupgauss_stat(3)] = jbtest(uptimes);
[~, RTdowngauss_p(3), RTdowngauss_stat(3)] = jbtest(downtimes);

if figures==1
figure;
subplot(2,2,1);histfit(uptimes,nbinsup); title(sprintf('%s%s','Minimum time = ',num2str(mintime)));ylabel('RT up');
subplot(2,2,2);set(gca, 'visible', 'off');text(0,0.6,sprintf('Mean = %.5f\nStd. Dev. = %.5f\nMedian = %.5f\nMin = %.5f\nMax = %.5f\nLillie(exp) = %.1e\nLillie(gauss) = %.1e\nKS = %.1e\nJB = %.1e\nNumber of counts = %s',mean(uptimes),std(uptimes),median(uptimes),min(uptimes),max(uptimes),RTupexp_stat,RTupgauss_stat(1),RTupgauss_stat(2),RTupgauss_stat(3),num2str(up_events)));
subplot(2,2,3);histfit(downtimes,nbinsdown); title(sprintf('%s%s','Minimum time = ',num2str(mintime)));ylabel('RT down');
subplot(2,2,4);set(gca, 'visible', 'off');text(0,0.6,sprintf('Mean = %.5f\nStd. Dev. = %.5f\nMedian = %.5f\nMin = %.5f\nMax = %.5f\nLillie(exp) = %.1e\nLillie(gauss) = %.1e\nKS = %.1e\nJB = %.1e\nNumber of counts = %s',mean(downtimes),std(downtimes),median(downtimes),min(downtimes),max(downtimes),RTdownexp_stat,RTdowngauss_stat(1),RTdowngauss_stat(2),RTdowngauss_stat(3),num2str(down_events)));
end

else
    %disp('Minimum number of observations not met for statistical test');
end
end

% Alternative way to solve the problem - Recenters the detrended signal
% using the min/max within the middle 50% of the signal.
if altsol ==1
%disp('Alternate solution...');
rangelow = round(0.25*length(Xvecfilt_detrended));
rangehigh = round(0.75*length(Xvecfilt_detrended));
Xvecfilt_detrended = Xvecfilt_detrended - mean([max(Xvecfilt_detrended(rangelow:rangehigh)) min(Xvecfilt_detrended(rangelow:rangehigh))]);

if figures==1
figure;
subplot(4,1,2); plot(time,movingavgX(round(length(Xvecfilt)*0.2:end))); title(sprintf('%s %s','window = ',num2str(win))); ylabel('moving average');
movingstdX = movingstd(Xvec,win,'central');
subplot(4,1,3); plot(time,movingstdX(round(length(Xvecfilt)*0.2:end))); title(sprintf('%s %s','window = ',num2str(win))); ylabel('moving standard deviation');
subplot(4,1,4); plot(time,Xvecfilt_detrended(1:length(time))); ylabel('filtered time signal'); title('filtered and detrended signal');
end

% Find up and down positiosn
Mup = zeros(1,length(Xvecfilt_detrended));          % initialize variables
Mdown = zeros(1,length(Xvecfilt_detrended));
sup = Xvecfilt_detrended > 0;                 % find up/down positions
sdown = Xvecfilt_detrended < 0;
Mup(sup) = 1;
Mdown(sdown) = 1;

updiffend = find(diff(Mup)==-1);           % find up times
updiffstart = find(diff(Mup)==1);

if isempty(updiffend)==0 && isempty(updiffstart)==0
if updiffend(1)<updiffstart(1)
    updiffend(1)=[];
end
for k = 2:length(updiffend)
        if updiffend(k)>updiffstart(k-1)
            Mup(updiffstart(k-1):updiffend(k)) = 1;
        else
            Mup(updiffend(k-1):updiffstart(k)) = 1;
        end
end
end
up_events = min([length(updiffstart) length(updiffend)]);



for i = 1:up_events
    uptimes(i) = time(updiffend(i)) - time(updiffstart(i));
end

downdiffend = find(diff(Mdown)==-1);         % find down times
downdiffstart = find(diff(Mdown)==1);

if isempty(downdiffend)==0 && isempty(downdiffstart)==0
if downdiffend(1)<downdiffstart(1)
    downdiffend(1)=[];
end
for k = 2:length(downdiffend)
        if downdiffend(k)>downdiffstart(k-1)
            Mdown(downdiffstart(k-1):downdiffend(k)) = 1;
        else
            Mdown(downdiffend(k-1):downdiffstart(k)) = 1;
        end
end
end
down_events = min([length(downdiffstart) length(downdiffend)]);



for i = 1:down_events
    downtimes(i) = time(downdiffend(i)) - time(downdiffstart(i));
end

binsizeup = 2*iqr(uptimes)*length(uptimes)^(-1/3);      % freedman-diaconis rule
binsizedown = 2*iqr(downtimes)*length(downtimes)^(-1/3);
nbinsup = round((max(uptimes) - min(uptimes))/binsizeup);
nbinsdown = round((max(downtimes) - min(downtimes))/binsizedown);
if nbinsup > 1e4                    % maximum limit on bins to reduce memory use
    nbinsup = 1e4;
else if isempty(nbinsup) ==1 || isnan(nbinsup)
        nbinsup=2;
    else if nbinsup < 4
            nbinsup = 4;
        end
    end    
end
if nbinsdown > 1e4
    nbinsdown = 1e4;
else if isempty(nbinsdown) ==1 || isnan(nbinsdown)
        nbinsdown=2;
    else if nbinsdown < 4
            nbinsdown = 4;
        end
    end
end

%disp('Running statistics...');
% Run statistics
[~, RTupexp_p, RTupexp_stat] = lillietest(uptimes,'Distr','exp');
[~, RTupgauss_p, RTupgauss_stat] = lillietest(uptimes);
[~, RTdownexp_p, RTdownexp_stat] = lillietest(downtimes,'Distr','exp');
[~, RTdowngauss_p, RTdowngauss_stat] = lillietest(downtimes);
[~, RTupgauss_p(2), RTupgauss_stat(2)] = kstest(uptimes);
[~, RTdowngauss_p(2), RTdowngauss_stat(2)] = kstest(downtimes);
[~, RTupgauss_p(3), RTupgauss_stat(3)] = jbtest(uptimes);
[~, RTdowngauss_p(3), RTdowngauss_stat(3)] = jbtest(downtimes);

if figures==1
figure; title('BEFORE MINIMUM TIME IMPLEMENTED');
subplot(2,2,1);histfit(uptimes,nbinsup);title('BEFORE MINIMUM TIME IMPLEMENTED'); ylabel('RT up');
subplot(2,2,2);set(gca, 'visible', 'off');text(0,0.6,sprintf('Mean = %.5f\nStd. Dev. = %.5f\nMedian = %.5f\nMin = %.5f\nMax = %.5f\nLillie(exp) = %.1e\nLillie(gauss) = %.1e\nKS = %.1e\nJB = %.1e\nNumber of counts = %s',mean(uptimes),std(uptimes),median(uptimes),min(uptimes),max(uptimes),RTupexp_stat,RTupgauss_stat(1),RTupgauss_stat(2),RTupgauss_stat(3),num2str(up_events)));
subplot(2,2,3);histfit(downtimes,nbinsdown);title('BEFORE MINIMUM TIME IMPLEMENTED');ylabel('RT down');
subplot(2,2,4);set(gca, 'visible', 'off');text(0,0.6,sprintf('Mean = %.5f\nStd. Dev. = %.5f\nMedian = %.5f\nMin = %.5f\nMax = %.5f\nLillie(exp) = %.1e\nLillie(gauss) = %.1e\nKS = %.1e\nJB = %.1e\nNumber of counts = %s',mean(downtimes),std(downtimes),median(downtimes),min(downtimes),max(downtimes),RTdownexp_stat,RTdowngauss_stat(1),RTdowngauss_stat(2),RTdowngauss_stat(3),num2str(down_events)));
end

if isempty(mintime) == 0
uptimes(uptimes<mintime) = [];
downtimes(downtimes<mintime) = [];
up_events = length(uptimes);
down_events = length(downtimes);

if length(uptimes)>=4 && length(downtimes)>=4
% Run statistics
[~, RTupexp_p, RTupexp_stat] = lillietest(uptimes,'Distr','exp');
[~, RTupgauss_p, RTupgauss_stat] = lillietest(uptimes);
[~, RTdownexp_p, RTdownexp_stat] = lillietest(downtimes,'Distr','exp');
[~, RTdowngauss_p, RTdowngauss_stat] = lillietest(downtimes);
[~, RTupgauss_p(2), RTupgauss_stat(2)] = kstest(uptimes);
[~, RTdowngauss_p(2), RTdowngauss_stat(2)] = kstest(downtimes);
[~, RTupgauss_p(3), RTupgauss_stat(3)] = jbtest(uptimes);
[~, RTdowngauss_p(3), RTdowngauss_stat(3)] = jbtest(downtimes);

if figures==1
figure;title(sprintf('%s%s','Minimum time = ',num2str(mintime)));
subplot(2,2,1);histfit(uptimes,nbinsup);title(sprintf('%s%s','Minimum time = ',num2str(mintime))); ylabel('RT up');
subplot(2,2,2);set(gca, 'visible', 'off');text(0,0.6,sprintf('Mean = %.5f\nStd. Dev. = %.5f\nMedian = %.5f\nMin = %.5f\nMax = %.5f\nLillie(exp) = %.1e\nLillie(gauss) = %.1e\nKS = %.1e\nJB = %.1e\nNumber of counts = %s',mean(uptimes),std(uptimes),median(uptimes),min(uptimes),max(uptimes),RTupexp_stat,RTupgauss_stat(1),RTupgauss_stat(2),RTupgauss_stat(3),num2str(up_events)));
subplot(2,2,3);histfit(downtimes,nbinsdown); title(sprintf('%s%s','Minimum time = ',num2str(mintime)));ylabel('RT down');
subplot(2,2,4);set(gca, 'visible', 'off');text(0,0.6,sprintf('Mean = %.5f\nStd. Dev. = %.5f\nMedian = %.5f\nMin = %.5f\nMax = %.5f\nLillie(exp) = %.1e\nLillie(gauss) = %.1e\nKS = %.1e\nJB = %.1e\nNumber of counts = %s',mean(downtimes),std(downtimes),median(downtimes),min(downtimes),max(downtimes),RTdownexp_stat,RTdowngauss_stat(1),RTdowngauss_stat(2),RTdowngauss_stat(3),num2str(down_events)));
end

else
    %disp('Minimum number of observations not met for statistical test');
end
end
end

if sopt==1
%display('Saving...');
%save('/Users/joshsalvi/Documents/Lab/Lab/Simulation Data/Sinusoids/whitenoisemovingavg.mat','Xvec','Fs','time','win','movingavgX','movingstdX','startl','endl','Fp','kp','F_rand','k_rand');
%save('/Users/joshsalvi/Documents/Lab/Lab/Simulation Data/ONHFishJosh/Analysis/Fp1-kp1-start1end4e5-xfish1.0Noise.mat','Xvec','Fs','time','win','movingavgX','movingstdX','startl','endl','Fp','kp','F_rand','k_rand');
save(out);
end

%display('Finished.');
warning on;
end
