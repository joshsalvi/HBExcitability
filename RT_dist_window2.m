% This function uses a moving filter as opposed to detrending the data, v.2

clear all; close all;

display('Importing...');

%load('/Users/joshsalvi/Downloads/4.mat');
load('/Users/joshsalvi/Documents/Lab/Lab/Simulation Data/ONHFishJosh/xfish1.0Noise.mat');

Fs = 10e3;                          % Choose the sampling frequency (note that this won't work properly for simulation data)
% Select operating points

% FORCE AND STIFFNESS


Ord_k = unique(k_rand);
Ord_F = unique(F_rand);

Fp = 4;                             % SELECT OPERATING POINT
kp = 3;

[Np,Nt] = ind2sub(size(F_rand),find(k_rand==Ord_k(kp) & F_rand==Ord_F(Fp)));

startl = 30e4;                      % how much time would you like to use?
endl = 35e4;
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
wins = 100;                                         % CHOOSE SMOOTHING WINDOW (if used)
movingavgX = filter(ones(1,win)/win,1,Xvec);
subplot(3,1,2); plot(time,movingavgX); title(sprintf('%s %s','window = ',num2str(win))); ylabel('moving average');
movingstdX = movingstd(Xvec,win,'central');
subplot(3,1,3); plot(time,movingstdX); title(sprintf('%s %s','window = ',num2str(win))); ylabel('moving standard deviation');


n = -3:0.2:3;

up_ind = zeros(length(Xvec),length(n));



display('Find events...');
for j = 1:length(n)
    Mvec = zeros(1,length(Xvec));
for i = 1:length(movingavgX)
    if i <= win
        clear a;
        a = Xvec(i:round(win/2)+i-1) > movingavgX(i) + n(j)*movingstdX(i);            % left edge
        Mvec(i:round(win/2)+i-1) = Mvec(i:round(win/2)+i-1)+a(:)';
    else if i < length(movingavgX) - win
            clear a;
            a = Xvec(i-floor(win/2):round(win/2)+i-1) > movingavgX(i) + n(j)*movingstdX(i);           % middle
            Mvec(i-floor(win/2):round(win/2)+i-1) = Mvec(i-floor(win/2):round(win/2)+i-1)+a(:)';
        else
            clear a;
            a = Xvec(i-floor(win/2):i) > movingavgX(i) + n(j)*movingstdX(i);          % right edge
            Mvec(i-floor(win/2):i) = Mvec(i-floor(win/2):i)+a(:)';
        end
    end
end
Mvec(Mvec > 0) = 1;   % reduce multiple counts to 1
up_ind(:,j) = Mvec;
up_ind(:,j) = ceil(smooth(up_ind(:,j),wins));        % remove points below window threshold (i.e., concatenate spikes)
events(:,j) = floor(sum(abs(diff(up_ind(:,j))))/2);
end
events_diff = diff(events);
