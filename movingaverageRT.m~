% Define moving window

clear all; close all;

display('Importing...');

%load('/Users/joshsalvi/Downloads/4.mat');
%load('/Users/joshsalvi/Documents/Lab/Lab/Simulation Data/ONHFishJosh/xfish1.0Noise.mat');
%load('/Users/joshsalvi/Documents/Lab/Lab/Simulation Data/Sinusoids/noisysinewave.mat');
load('/Users/joshsalvi/Documents/Lab/Lab/Simulation Data/Sinusoids/brownsine.mat');

Fs = 10e3;                          % Choose the sampling frequency (note that this won't work properly for simulation data)
% Select operating points

% FORCE AND STIFFNESS


Ord_k = unique(k_rand);
Ord_F = unique(F_rand);

Fp = 1;                             % SELECT OPERATING POINT
kp = 1;

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
display('Moving average and standard deviation...');
time = linspace(0,length(Xvec)/Fs,length(Xvec));    % time in seconds

figure(1); subplot(3,1,1); plot(time,Xvec); title(sprintf('%s%s %s%s','F= ',num2str(F_rand(Np,Nt)),' k= ',num2str(k_rand(Np,Nt))));ylabel('displacement')
% Moving window, running average using filter
win = 51;                                        % CHOOSE WINDOW, must be ODD
movingavgX = filter(ones(1,win)/win,1,Xvec);
subplot(3,1,2); plot(time,movingavgX); title(sprintf('%s %s','window = ',num2str(win))); ylabel('moving average');
movingstdX = movingstd(Xvec,win,'central');
subplot(3,1,3); plot(time,movingstdX); title(sprintf('%s %s','window = ',num2str(win))); ylabel('moving standard deviation');

display('Saving...');
%save('/Users/joshsalvi/Documents/Lab/Lab/Simulation Data/Sinusoids/whitenoisemovingavg.mat','Xvec','Fs','time','win','movingavgX','movingstdX','startl','endl','Fp','kp','F_rand','k_rand');
%save('/Users/joshsalvi/Documents/Lab/Lab/Simulation Data/ONHFishJosh/Analysis/Fp1-kp1-start1end4e5-xfish1.0Noise.mat','Xvec','Fs','time','win','movingavgX','movingstdX','startl','endl','Fp','kp','F_rand','k_rand');
display('Finished.');
