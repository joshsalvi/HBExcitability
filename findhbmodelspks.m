function [Fspk, kc, X, numspk, pks, trs, ampl] = findhbmodelspks(Fmin,Fmax,kc,noiselevel,thresh,dF,tvec)
%
% This function calculates the border of the region over which a hair
% bundle exhibits spontaneous motion.
%
% [Fspk, kc, X, numspk, pks, trs, ampl] = findhbmodelspks(Fmin,Fmax,kc,noiselevel,thresh,dF,tvec)
%
% INPUTS:
% Fmin, Fmax: minimum and maximum force value to be analyzed
% kc: stiffness
% noiselevel: noise level (std of the white noise vector)
% thresh: peak-detection threshold
% dF: minimum distance between force values that must be reached upon
% finding the border region (the final resolution)
% tvec: time vector
%
% OUTPUTS:
% Fspk: force value at which the number of spikes first reaches zero
% kc: input stiffness
% X: time series at the border region
% numspk: number of spikes at (Fspk,kc)
% pks/trs: array of peaks and troughs (as indices)
% ampl: pk-tr amplitude at (Fspk,kc)
%
% REQUIREMENTS:
% This function calls hbtoymodel().
%
%
%
% Joshua D. Salvi
% jsalvi@rockefeller.edu
%
%


Fc = (Fmin + Fmax) / 2;                                     % starting force
Fdiff = Fc - Fmin;                                          % starting difference between forces
Fmax_old = Fmax; Fmin_old = Fmin;

dispyn = 1;                                                 % display progress (1=yes)

% Loop through values of force until a minimum difference, df, between the
% subsequent forces is reached.
if dispyn == 1
    disp('Calculating...');
end

while Fdiff > dF
    
    if dispyn == 1
        disp(['dF = ' num2str(Fdiff) '  Fc = ' num2str(Fc)]);
    end
    
    clear Xsto X pk tr
    
    [~, Xsto] = hbtoymodel(Fc,kc,noiselevel,0,[0 0],tvec);      % generate simulation data
    X = Xsto(1,round(size(Xsto,2)/2):end);                  % extract last 50% of signal
    [pk, tr] = PTDetect(X,thresh);                           % find peaks and troughs
    
    Fc_old = Fc;
    
    if isempty(pk) == 0 || isempty(tr) == 0                 % if peaks exist, increase the force
        Fmin_old = Fc;
        Fc = (Fc + Fmax_old) / 2;
        Fdiff = Fmax_old - Fc;
    elseif isempty(pk) == 1 || isempty(tr) == 1             % if no peaks, decrease the force
        Fmax_old = Fc;
        Fc = (Fc + Fmin_old) / 2;
        Fdiff = Fc - Fmin_old;
    end  
    
end
if dispyn == 1
    disp('Final check...');
end

clear X Xsto pk tr

% Check which time series possess a spike
[~, Xsto] = hbtoymodel(Fc,kc,noiselevel,0,[0 0],tvec);          % Force 1
X1 = Xsto(1,round(size(Xsto,2)/2):end);
[pk1, tr1] = PTDetect(X1,thresh);

[~, Xsto] = hbtoymodel(Fc_old,kc,noiselevel,0,[0 0],tvec);      % Force 2
X2 = Xsto(1,round(size(Xsto,2)/2):end);
[pk2, tr2] = PTDetect(X2,thresh);

if (isempty(pk1) == 0 || isempty(tr1) == 0) && (isempty(pk2) == 1 || isempty(tr2) == 1)
    Fspk = Fc;
    numspk = max([length(tr1) length(pk1)]);
    X = X1;
    pks = pk1;trs=tr1;
elseif (isempty(pk1) == 1 || isempty(tr1) == 1) && (isempty(pk2) == 0 || isempty(tr2) == 0)
    Fspk = Fc_old;
    numspk = max([length(tr2) length(pk2)]);
    X = X2;
    pks = pk2;trs=tr2;
elseif (isempty(pk1) == 0 || isempty(tr1) == 0) && (isempty(pk2) == 0 || isempty(tr2) == 0)
    [~, Xsto] = hbtoymodel(Fmax_old,kc,noiselevel,0,[0 0],tvec);      % Force 2
    X2 = Xsto(1,round(size(Xsto,2)/2):end);
    [pk2, tr2] = PTDetect(X2,thresh);
    if isempty(pk2) == 0 || isempty(tr2) == 0
        Fspk = Fmax_old;
        numspk = max([length(tr2) length(pk2)]);
        X = X2;
        pks = pk2;trs=tr2;
    else
        Fspk = (Fc + Fc_old) / 2;
        [~, Xsto] = hbtoymodel(Fspk,kc,noiselevel,0,[0 0],tvec);
        X = Xsto(1,round(size(Xsto,2)/2):end);
        [pk2, tr2] = PTDetect(X2,thresh);
        numspk = max([length(tr2) length(pk2)]);
        pks = pk2;trs=tr2;
    end
elseif (isempty(pk1) == 1 || isempty(tr1) == 1) && (isempty(pk2) == 1 || isempty(tr2) == 1)
    [~, Xsto] = hbtoymodel(Fmin_old,kc,noiselevel,0,[0 0],tvec);      % Force 2
    X2 = Xsto(1,round(size(Xsto,2)/2):end);
    [pk2, tr2] = PTDetect(X2,thresh);
    if isempty(pk2) == 0 || isempty(tr2) == 0
        Fspk = Fmin_old;
        numspk = max([length(tr2) length(pk2)]);
        X = X2;
        pks=pk2;trs=tr2;
    else
        disp('No peaks found.');
        Fspk = NaN;
        numspk = NaN;
        X = NaN;
    end
    return;
else
    disp('Error.');
    return;
end
%}
if isempty(pks)==0 && isempty(trs)==0
    ampl = abs(X(pks(1)) - X(trs(1)));
elseif isempty(pks)==0 && isempty(trs)==1
    ampl = abs(X(pks(1)) - min(X));
elseif isempty(pks)==1 && isempty(trs)==0
    ampl = abs(max(X) - X(trs(1)));
else
    ampl = abs(max(X) - min(X));
end
if exist('ampl')==0
    ampl=NaN;
end
if isempty(ampl)==1
    ampl=NaN;
end

if dispyn == 1
    disp('Finished.');
end


end
