% morris_lecar_gillespie.m
% solves the morris-lecar equations with stochastic K+ channels
% based on Bard Ermentrout's example of M-L model generating Hopf
% bifurcation
% http://www.math.pitt.edu/~bard/bardware/meth3/odes/mlfig2.ode
% call from MATLAB prompt as follows: output = morris_lecar_gilliespie
% output is a cell array containing vectors of time, membrane potential
% from stochastic simulations, and membrane potential from deterministic
% simulations

function output = morris_lecar_gillespie

Nkvect = [1 10 100 1000]; % number of K+ channels

output = cell(length(Nkvect),1); % output cell array
figure
hold
for i = 1:length(Nkvect)
    [t,v,vdet] = morris_lecar_gillespie_iteration(Nkvect(i));
    size(t)
    output{i,1} = [t' v' vdet'];
    plot(t,v+(length(Nkvect)+1-i)*120)
end
% add deterministic plot from the last simulation
plot(t,vdet,'k')
tscale = [1000 1100 1100];
vscale = [-70 -70 30];
plot(tscale,vscale,'k')
axis off

end % main function



% call sequence for each iteration
% input: Nk, the number of K+ channels
% output: t, the sequence of event times
% output: v, the voltage vector for stochastic simulation
% output: vdet, the voltage vector for deterministic simulation
function [t, v, vdet] = morris_lecar_gillespie_iteration(Nk)

% initial conditions
v0=-60.855;
w0=0.014915;
t0 = 0;

% other parameters
iapp=230;vk=-84;vl=-60;vca=120;
gk=8;gl=2;gca=4.4;c=20;
tstart = 100; % time to turn on stimulus
tend = 5000; % time to end simulation
dtmin = 0.1; % minimum step size

i = 1;  % counter, updated with each stochastic transition
t(i) = t0; % first reported time

% initial conditions of dynamic and stochastic variables
No = round(winf(v0)*Nk); % set number of open channels at equilibrium value
Nc = Nk - No; % number of closed channels
N_open(i) = No; 
N_closed(i) = Nc;
v(i) = v0; 
vdet(i) = v0;
w(i) = w0;
wdet(i) = w0;

% seed random number generator to unique value each run
rand('twister', sum(100*clock));

% begin stochastic simulation
while (t(i) < tend)
    i = i+1;  % i-1 = counter for number of state transitions
    r = rand(2,1);  % generate r(1) and r(2), uniformly distributed over [0,1)
    lm = Nc*aw(v(i-1)) + No*bw(v(i-1)); % effective rate constant of next transition
    dt = 1/lm*log(1/r(1));
    if (dt > dtmin)
        dt = dtmin; % impose a minimum step size with no change in No or Nc
    elseif (No*bw(v(i-1)) > r(2)*lm) % channel closing
        No = No - 1;
        Nc = Nc + 1;
    else
        No = No + 1;
        Nc = Nc - 1;
    end
    t(i) = t(i-1) + dt; %update time vector; note that time is not equally spaced
    N_open(i)=No; 
    N_closed(i)=Nc;
    gkstoch = gk*No/Nk; % stochastic K+ conductance
    % integrate current-balance ODE using forward Euler -- stochastic
    % version
    v(i) = v(i-1) + dt/c*(iapp*(t(i-1)>tstart) - gca*minf(v(i-1))*(v(i-1)-vca) - gkstoch*(v(i-1)-vk) - gl*(v(i-1)-vl));
    % integrate current-balance ODE using forward Euler -- deterministic
    % version    
    vdet(i) = vdet(i-1) + dt/c*( iapp*(t(i-1)>tstart) - gca*minf(vdet(i-1))*(vdet(i-1)-vca) - gk*wdet(i-1)*(vdet(i-1)-vk) - gl*(vdet(i-1)-vl) );
    wdet(i) = wdet(i-1) + dt*(winf(vdet(i-1))-wdet(i-1))/tauw(vdet(i-1));        
end

end

function winf = winf(v)
v3=2;v4=30;
winf=.5*(1+tanh((v-v3)/v4));
end

function minf = minf(v)
v1=-1.2;v2=18;
minf = .5*(1+tanh((v-v1)/v2));
end

function tauw = tauw(v)
v3=2;v4=30;phi=0.04;
tauw = 1/phi/cosh((v-v3)/(2*v4)); % we include phi in tauw to make it easy to do stochastic calculations
end

% we'll use winf and tauw to get the rate constants

function aw = aw(v)
aw = winf(v)/tauw(v);
end

function bw = bw(v)
bw = (1-winf(v))/tauw(v);
end
