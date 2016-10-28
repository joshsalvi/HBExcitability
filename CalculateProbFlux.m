function [VLsum, pval, VLsum_H0, b, a] = CalculateProbFlux(x,varargin)
%
% Calculates the probability fluxes for a two-dimensional probability
% distribution.
%
% [VLsum, pval, VLsum_H0, b, a] = CalculateProbFlux(x, Nboot, plotyn)
%
% INPUTS
% ------
%   x:          input signal
%   Nboot:      number of boostraps (default = 1)
%   plotyn:     plot? (1 = yes; default = 0)
%
% OUTPUTS
% -------
%   VLsum:      vector length summation for x
%   pval:       one-sided p-value for VLsum
%   VLsum_H0:   bootstrapped vector lengths for normally-distributed values
%   (b,a):      kernel density estimate used for calculation of p-value
%
% Written by:
%   Joshua D. Salvi
%   jsalvi@rockefeller.edu
%

if nargin == 1
    Nboot = 1;
    plotyn = 0;
elseif nargin == 2
    Nboot = varargin{1};
    plotyn = 0;
elseif nargin == 3
    Nboot = varargin{1};
    plotyn = varargin{2};
end

% Length of x
L = length(x);

% Hilbert transform
xh = hilbert(x);xi=(imag(xh));

% 2D phase space
Nbins = freedmandiaconis(x)*freedmandiaconis(xi)/2;
if iscolumn(x) == 0
    [bw dens mx my]=kde2d([x' xi'],Nbins);
else
    [bw dens mx my]=kde2d([x xi],Nbins);
end

% Initialization
fluxx = zeros(length(mx),length(my));
fluxy = zeros(length(mx),length(my));
xycount = zeros(length(mx),length(my));

disp('Calculating flux for X...')

for j = 1:L
    % Calulate fluxes
    
    xind0 = find(x(j)>=mx(1,:));xind0 = xind0(end);
    xval(j) = xind0;
    
    yind0 = find(xi(j)>=my(:,1));yind0 = yind0(end);
    yval(j) = yind0;
    
    if j > 1
        xycount(xval(j-1),yval(j-1)) = xycount(xval(j-1),yval(j-1)) + 1;
        if yval(j) - yval(j-1) > 0
            fluxx(xval(j-1),yval(j-1)) = fluxx(xval(j-1),yval(j-1)) + 1;
        else
            fluxx(xval(j-1),yval(j-1)) = fluxx(xval(j-1),yval(j-1)) - 1;
        end
        if xval(j) - xval(j-1) > 0
            fluxy(xval(j-1),yval(j-1)) = fluxy(xval(j-1),yval(j-1)) + 1;
        else
            fluxy(xval(j-1),yval(j-1)) = fluxy(xval(j-1),yval(j-1)) - 1;
        end
    end
    
end

% Summary statistics
fluxx = fluxx.*xycount;fluxy = fluxy.*xycount;  % Normalize by number of counts
phase0 = atan2(fluxy,fluxx);
veclength = sqrt(fluxx.^2+fluxy.^2)./(size(fluxx,1)*size(fluxx,2)); % Calculate vector lengths
VLsum = sum(sum(veclength));    % Summation of vector lengths

% Plot the distribution and vectors
if plotyn == 1
    figure;
    subplot(1,2,1);
    pcolor(mx(1,:),my(:,1),dens); shading interp;load jetnew.mat;colormap(bone);
    hold on;
    [xn yn] = meshgrid(mx(1,:),my(:,1)');
    quiver(xn,yn,fluxx,fluxy,'LineStyle','-','AutoScaleFactor',5,'Color','r');
end


% BOOTSTRAPPING

if Nboot > 0
    
disp('Bootstrapping...')

for m = 1:Nboot
    
    x0 = randn(1,L);
    x0h = hilbert(x0);x0i=(imag(x0h));
    Nbins = freedmandiaconis(x0)*freedmandiaconis(x0i);
    [bw dens mx my]=kde2d([x0' x0i'],Nbins);
    
    fluxx = zeros(length(mx),length(my));
    fluxy = zeros(length(mx),length(my));
    xycount = zeros(length(mx),length(my));
    xval = zeros(1,L);
    yval = zeros(1,L);
    xs = zeros(1,L);
    
    for j = 1:L
        
        xind0 = find(x0(j)>=mx(1,:));xind0 = xind0(end);
        xval(j) = xind0;

        yind0 = find(x0i(j)>=my(:,1));yind0 = yind0(end);
        yval(j) = yind0;

        if j > 1
            xycount(xval(j-1),yval(j-1)) = xycount(xval(j-1),yval(j-1)) + 1;
            if yval(j) - yval(j-1) > 0
            fluxx(xval(j-1),yval(j-1)) = fluxx(xval(j-1),yval(j-1)) + 1;
            else
                fluxx(xval(j-1),yval(j-1)) = fluxx(xval(j-1),yval(j-1)) - 1;
            end
            if xval(j) - xval(j-1) > 0
                fluxy(xval(j-1),yval(j-1)) = fluxy(xval(j-1),yval(j-1)) + 1;
            else
                fluxy(xval(j-1),yval(j-1)) = fluxy(xval(j-1),yval(j-1)) - 1;
            end
        end

    end
    
    % Summary statistics
    fluxx = fluxx.*xycount;fluxy = fluxy.*xycount;
    phase0 = atan2(fluxy,fluxx);
    veclength = sqrt(fluxx.^2+fluxy.^2)./(size(fluxx,1)*size(fluxx,2));
    VLsum_H0(m) = sum(sum(veclength));
    
    % Status indicator
    if mod(m,round(Nboot/10)) == 0
        disp([num2str(m/Nboot*100) '% complete.']);
    end
end

if plotyn == 1
    subplot(1,2,2);
    pcolor(mx(1,:),my(:,1),dens); shading interp;load jetnew.mat;colormap(gray);
    hold on;
    [xn yn] = meshgrid(mx(1,:),my(:,1)');
    quiver(xn,yn,fluxx,fluxy,'LineStyle','-','AutoScaleFactor',5,'Color','r');
end

% Kernel density estimate for null distribution
[a, b] = ksdensity(VLsum_H0,min(VLsum_H0):1e-5:max(VLsum_H0));
a = a./sum(a);

% Calculate p-value
q = findnearest(b,VLsum);
pval = 1 - sum(a(1:q));

else
    
    a=NaN;b=NaN;pval=NaN;VLsum_H0=NaN;
    
end

end
