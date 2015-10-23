
import numpy as np
import scipy
import matcompat

# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass

clear(all)
plt.close(all)
display('Importing...')
#%load('/Users/joshsalvi/Downloads/4.mat');
np.load('/Users/joshsalvi/Documents/Lab/Lab/Simulation Data/ONHFishJosh/xfish1.0Noise.mat')
Fs = 10e3
#% Choose the sampling frequency (note that this won't work properly for simulation data)
#% Select operating points
#% FORCE AND STIFFNESS
Np = 4.
Nt = 4.
Xvec = Xd(:, Np, Nt)
clear(Xd)
#% FORCE ONLY
#%{
Fi = 3.
ind = Fi
Xvec = Xd(:, ind)
clear(Xd)
#%}
time = np.linspace(0., matdiv(length(Xvec), Fs), length(Xvec))
#% time in seconds
plt.figure(1.)
plt.subplot(2., 1., 1.)
plt.plot(time, Xvec)
plt.title(sprintf('%s%s %s%s', 'F= ', num2str(F_rand(Np, Nt)), ' k= ', num2str(k_rand(Np, Nt))))
#% Detrend data
win = 8e3
#% CHOOSE WINDOW
Xvec = Xvec-smooth(Xvec, win)
plt.subplot(2., 1., 2.)
plt.plot(time, Xvec)
plt.title(sprintf('%s %s', 'window = ', num2str(win)))
clear(win)
display('Convergence of standard deviations...')
#% Find number of events
#% Convergence for nup and ndown
nup = np.arange(0., 5.05, 0.05)
#% Initialize nup and ndown
ndown = np.arange(0., 5.05, 0.05)
wintime = 0.001
#% CHOOSE WINDOW TIME IN SECONDS
win = findnearest(time, wintime)
#% OPTION 1 - No stepping
mux = np.mean(Xvec)
stdx = np.std(Xvec)
display('UP')
for i in np.arange(1., (length(nup))+1):
    Xvecind_up[:,int(i)-1] = Xvec >= mux+np.dot(nup[int(i)-1], stdx)
    Xvecind_up[:,int(i)-1] = np.ceil(smooth(Xvecind_up[:,int(i)-1], win))
    #% remove points below window threshold (i.e., concatenate spikes)
    events_up[int(i)-1] = np.floor((np.sum(np.abs(np.diff(Xvecind_up[:,int(i)-1])))/2.))
    
events_up_diff = np.diff(events_up)
display('DOWN')
for j in np.arange(1., (length(ndown))+1):
    Xvecind_down[:,int(j)-1] = Xvec<=mux-np.dot(ndown[int(j)-1], stdx)
    Xvecind_down[:,int(j)-1] = np.ceil(smooth(Xvecind_down[:,int(j)-1], win))
    #% remove points below window threshold (i.e., concatenate spikes)
    events_down[int(j)-1] = np.floor((np.sum(np.abs(np.diff(Xvecind_down[:,int(j)-1])))/2.))
    
events_down_diff = np.diff(events_down)
plt.figure(2.)
plt.subplot(2., 2., 1.)
plt.plot(nup, events_up)
plt.title('Up Events')
plt.ylabel('Events')
plt.xlabel('nup')
plt.subplot(2., 2., 2.)
plt.plot(ndown, events_down)
plt.title('Down Events')
plt.ylabel('Events')
plt.xlabel('ndown')
plt.subplot(2., 2., 3.)
plt.plot(nup[1:], events_up_diff)
plt.ylabel('dEvents')
plt.xlabel('nup')
plt.subplot(2., 2., 4.)
plt.plot(ndown[1:], events_down_diff)
plt.ylabel('dEvents')
plt.xlabel('nup')
#% METHOD A - Find local maxima
#%{
nupconv = nonzero((np.diff(np.sign(events_up_diff)) == -2.))+1.
ndownconv = nonzero((np.diff(np.sign(events_down_diff)) == -2.))+1.
if length(nupconv) > 1.:
    #% in case there is more than one point found
nupconv = nonzero((events_up == matcompat.max(events_up[int(nupconv)-1])))

if length(ndownconv) > 1.:
    ndownconv = nonzero((events_down == matcompat.max(events_down[int(ndownconv)-1])))


plt.figure(2.)
plt.subplot(2., 2., 1.)
plt.hold(on)
plt.scatter(nup[int(nupconv)-1], events_up[int(nupconv)-1], 'ro')
plt.subplot(2., 2., 2.)
plt.hold(on)
plt.scatter(ndown[int(ndownconv)-1], events_down[int(ndownconv)-1], 'ro')
#%}
#% METHOD B - Find local minima, method 1
nupconv = nonzero((np.diff(np.sign(events_up_diff)) == 2.))+1.
ndownconv = nonzero((np.diff(np.sign(events_down_diff)) == 2.))+1.
if isempty(nupconv) == 1.:
    #% If no local minimum can be found, use a local maximum
nupconv = nonzero((np.diff(np.sign(events_up_diff)) == -2.))+1.

if isempty(ndownconv) == 1.:
    #% If no local minimum can be found, use a local maximum
ndownconv = nonzero((np.diff(np.sign(events_down_diff)) == -2.))+1.

if length(nupconv) > 1.:
    #% in case there is more than one point found
nupconv = nonzero((events_up == matcompat.max(events_up[int(nupconv)-1])))

if length(ndownconv) > 1.:
    ndownconv = nonzero((events_down == matcompat.max(events_down[int(ndownconv)-1])))


plt.figure(2.)
plt.subplot(2., 2., 1.)
plt.hold(on)
plt.scatter(nup[int(nupconv)-1], events_up[int(nupconv)-1], 'ro')
plt.subplot(2., 2., 2.)
plt.hold(on)
plt.scatter(ndown[int(ndownconv)-1], events_down[int(ndownconv)-1], 'ro')
#% METHOD C - Find local minima, method 2
#%{
updiffmin = nonzero((events_up_diff == matcompat.max(events_up_diff)))
updiffmax = nonzero((events_up_diff == matcompat.max(events_up_diff)))
downdiffmin = nonzero((events_down_diff == matcompat.max(events_down_diff)))
downdiffmax = nonzero((events_down_diff == matcompat.max(events_down_diff)))
if updiffmin<updiffmax:
    nupconv = nonzero((events_up[int(updiffmin)-1:updiffmax] == matcompat.max(events_up[int(updiffmin)-1:updiffmax])))+updiffmin-1.
else:
    nupconv = nonzero((events_up[int(updiffmax)-1:updiffmin] == matcompat.max(events_up[int(updiffmax)-1:updiffmin])))+updiffmax-1.
    

if downdiffmin<downdiffmax:
    ndownconv = nonzero((events_down[int(downdiffmin)-1:downdiffmax] == matcompat.max(events_down[int(downdiffmin)-1:downdiffmax])))+downdiffmin-1.
else:
    ndownconv = nonzero((events_down[int(downdiffmax)-1:downdiffmin] == matcompat.max(events_down[int(downdiffmax)-1:downdiffmin])))+downdiffmax-1.
    

if length(nupconv) > 1.:
    #% in case there is more than one point found
nupconv = nonzero((events_up == matcompat.max(events_up[int(nupconv)-1])))

if length(ndownconv) > 1.:
    ndownconv = nonzero((events_down == matcompat.max(events_down[int(ndownconv)-1])))


plt.figure(2.)
plt.subplot(2., 2., 1.)
plt.hold(on)
plt.scatter(nup[int(nupconv)-1], events_up[int(nupconv)-1], 'ro')
plt.scatter(nup[int(updiffmin)-1], events_up[int(updiffmin)-1], 'k.')
plt.scatter(nup[int(updiffmax)-1], events_up[int(updiffmax)-1], 'k.')
plt.subplot(2., 2., 2.)
plt.hold(on)
plt.scatter(ndown[int(ndownconv)-1], events_down[int(ndownconv)-1], 'ro')
plt.scatter(ndown[int(downdiffmin)-1], events_down[int(downdiffmin)-1], 'k.')
plt.scatter(ndown[int(downdiffmax)-1], events_down[int(downdiffmax)-1], 'k.')
plt.subplot(2., 2., 3.)
plt.hold(on)
plt.scatter(nup[int(nupconv)-1], events_up_diff[int(nupconv)-1], 'ro')
plt.scatter(nup[int(updiffmin)-1], events_up_diff[int(updiffmin)-1], 'k.')
plt.scatter(nup[int(updiffmax)-1], events_up_diff[int(updiffmax)-1], 'k.')
plt.subplot(2., 2., 4.)
plt.hold(on)
plt.scatter(ndown[int(ndownconv)-1], events_down_diff[int(ndownconv)-1], 'ro')
plt.scatter(ndown[int(downdiffmin)-1], events_down_diff[int(downdiffmin)-1], 'k.')
plt.scatter(ndown[int(downdiffmax)-1], events_down_diff[int(downdiffmax)-1], 'k.')
#%}
plt.figure(3.)
plt.plot(time, Xvec)
plt.hold(on)
plt.plot(time[Xvecind_down[:,int(ndownconv)-1]], Xvec[Xvecind_down[:,int(ndownconv)-1]], 'r.')
plt.plot(time[Xvecind_up[:,int(nupconv)-1]], Xvec[Xvecind_up[:,int(nupconv)-1]], 'g.')
plt.title(sprintf('%s %s%s %s%s %s%s%s%s', 'Spikes ', 'nup= ', num2str(nup[int(nupconv)-1]), 'ndown= ', num2str(ndown[int(ndownconv)-1]), 'eventsup/eventsdown= ', num2str(events_up[int(nupconv)-1]), '/', num2str(events_down[int(ndownconv)-1])))
#% OPTION 2 - Moving window
#%{
wintime = 10.
#% CHOOSE WINDOW TIME IN SECONDS
win = findnearest(time, wintime)
overlap = 0.
for i in np.arange(1., (length(Xvec)-win)+1):
    clear(Xvec2)
    Xvec2 = Xvec[int(i)-1:i+win-1.]
    mux = np.mean(Xvec2)
    stdx = np.std(Xvec2)
    for j in np.arange(1., (length(nup))+1):
        Xvec2ind_up[:,int(i)-1,int(j)-1] = Xvec2 > mux+np.dot(nup[int(j)-1], stdx)
        
    for k in np.arange(1., (length(ndown))+1):
        Xvec2ind_down[:,int(i)-1,int(j)-1] = Xvec2 > mux-np.dot(ndown[int(j)-1], stdx)
        
    
#%}
display('Find residence times...')
#% Find the residence time distribution
upindex = Xvecind_up[:,int(nupconv)-1]
downindex = Xvecind_down[:,int(ndownconv)-1]
#% indices for up and down positions
upinddiff = np.diff(upindex)
downinddiff = np.diff(downindex)
#% first derivative (+1 is start of excursion, -1 is end of excursion)
spikeindup_start = nonzero((upinddiff == 1.))
spikeindup_end = nonzero((upinddiff == -1.))
#%{
if spikeindup_start[0] > spikeindup_end[0]:
    spikeindup_start[0] == np.array([])


#%}
spikeinddown_start = nonzero((downinddiff == 1.))
spikeinddown_end = nonzero((downinddiff == -1.))
#%{
if spikeinddown_start[0] > spikeinddown_end[0]:
    spikeinddown_start[0] == np.array([])


#%}
#% Calculate residence times for the up position and for the down position
for i in np.arange(1., (length(spikeindup_start))+1):
    RT_up[int(i)-1] = np.dot(time[int(spikeindup_end[int(i)-1])-1]-time[int(spikeindup_start[int(i)-1])-1], 1e3)
    #% convert to ms
    
for i in np.arange(1., (length(spikeinddown_start))+1):
    RT_down[int(i)-1] = np.dot(time[int(spikeinddown_end[int(i)-1])-1]-time[int(spikeinddown_start[int(i)-1])-1], 1e3)
    #% convert to ms
    
#% Calculate inter-event intervals
for i in np.arange(1., (length(spikeindup_start)-1.)+1):
    IEI_up[int(i)-1] = np.dot(time[int(spikeindup_start[int((i+1.))-1])-1]-time[int(spikeindup_end[int(i)-1])-1], 1e3)
    #% convert to ms
    
for i in np.arange(1., (length(spikeinddown_start)-1.)+1):
    IEI_down[int(i)-1] = np.dot(time[int(spikeinddown_start[int((i+1.))-1])-1]-time[int(spikeinddown_end[int(i)-1])-1], 1e3)
    #% convert to ms
    
plt.figure(4.)
binsizeup = np.dot(2.*iqr(RT_up), matixpower(length(RT_up), -1./3.))
#% freedman-diaconis rule
binsizedown = np.dot(2.*iqr(RT_down), matixpower(length(RT_down), -1./3.))
nbinsup = np.round(matdiv(matcompat.max(RT_up)-matcompat.max(RT_up), binsizeup))
nbinsdown = np.round(matdiv(matcompat.max(RT_down)-matcompat.max(RT_down), binsizedown))
plt.subplot(2., 2., 1.)
plt.hist(RT_up, nbinsup)
plt.title('residence time - UP')
plt.subplot(2., 2., 2.)
plt.hist(RT_down, nbinsdown)
plt.title('residence time - DOWN')
clear(binsizeupbinsizedownnbinsupnbinsdown)
binsizeup = np.dot(2.*iqr(IEI_up), matixpower(length(IEI_up), -1./3.))
#% freedman-diaconis rule
binsizedown = np.dot(2.*iqr(IEI_down), matixpower(length(IEI_down), -1./3.))
nbinsup = np.round(matdiv(matcompat.max(IEI_up)-matcompat.max(IEI_up), binsizeup))
nbinsdown = np.round(matdiv(matcompat.max(IEI_down)-matcompat.max(IEI_down), binsizedown))
plt.subplot(2., 2., 3.)
plt.hist(RT_up, nbinsup)
plt.title('inter-event interval - UP')
plt.subplot(2., 2., 4.)
plt.hist(RT_down, nbinsdown)
plt.title('inter-event interval - DOWN')
display('Find buffer region and recalculate...')
bufferind = upindex.flatten(1)+downindex.flatten(1)
#% Buffer index is the index for all those points not included in the up or down excursions (0 is in buffer, 1 is not)
bufferinddiff = np.diff(bufferind)
buffer_start = nonzero((bufferinddiff == 1.))
buffer_end = nonzero((bufferinddiff == -1.))
#% Calculate residence times within the buffer
for i in np.arange(1., (length(spikeindup_start))+1):
    RT_buffer[int(i)-1] = np.dot(time[int(buffer_end[int(i)-1])-1]-time[int(buffer_start[int(i)-1])-1], 1e3)
    #% convert to ms
    
#% Calculate inter-buffer times
for i in np.arange(1., (length(spikeindup_start)-1.)+1):
    IEI_buffer[int(i)-1] = np.dot(time[int(buffer_start[int((i+1.))-1])-1]-time[int(buffer_end[int(i)-1])-1], 1e3)
    #% convert to ms
    
plt.figure(5.)
binsizebuffer = np.dot(2.*iqr(RT_buffer), matixpower(length(RT_buffer), -1./3.))
#% freedman-diaconis rule
nbinsbuffer = np.round(matdiv(matcompat.max(RT_buffer)-matcompat.max(RT_buffer), binsizebuffer))
plt.subplot(1., 2., 1.)
plt.hist(RT_buffer, nbinsbuffer)
plt.title('residence time - buffer')
binsizebuffer = np.dot(2.*iqr(IEI_buffer), matixpower(length(IEI_buffer), -1./3.))
#% freedman-diaconis rule
nbinsbuffer = np.round(matdiv(matcompat.max(IEI_buffer)-matcompat.max(IEI_buffer), binsizebuffer))
plt.subplot(1., 2., 2.)
plt.hist(IEI_buffer, nbinsbuffer)
plt.title('inter-event interval - buffer')