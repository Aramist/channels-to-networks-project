function  [ph,edges,r] = psth(times, binsize, fs, ntrials, triallen, plotflag, varargin)
% PSTH Computes the peri-stimulus time histogram from spike times.
% The routine plots the trial averaged spike rate as a function of time.
% H = PSTH(TIMES, BINSIZE, FS,NTRIALS,TRIALLEN)
% H = PSTH(TIMES, BINSIZE, FS,NTRIALS,TRIALLEN ,AXESHANDLE)
% TIMES - spike times (samples)
% BINSIZE - binwidth (ms) - 1000ms
% FS - sampling rate (hz)
% NTRIALS - number of trials
% TRIALLEN - length of a trial (samples)
% PLOTFLAG - 0/1 (no plot/plot)
% H - plot handle
%
% An example:
% %spike times can be specified in continuous time 
% %here we have 3 trials and a trial length of 1000 samples
% t = [10, 250, 900, 1300, 1600, 2405, 2900];
%
% %the same spike times can also be specified per trial
% t2 =[10, 250, 900, 300, 600, 405, 900];
% r = psth(t,10,1000,3,1000) ;
% r2 = psth(t2,10,1000,3,1000);
%
% Author: Rajiv Narayan
% askrajiv@gmail.com
% Boston University, Boston, MA

h_color ='k';
nin=nargin;

narginchk(6,7);

switch(nin)
 
 case 6 %no plot handle
%   figure;
%   h=gca;
  
 case 7
  if(ishandle(varargin{1}))
    h=varargin{1};
  else
    error('Invalid Plot handle');
  end

end

%Compute PSTH        
lastBin = binsize * ceil((triallen-1)*(1000/(fs*binsize)));
edges = 0 : binsize : lastBin;	
x = (mod(times-1,triallen)+1)*(1000/fs);
r = (histc(x,edges)*1000) / (ntrials*binsize);
ph=0;
% ph=bar(edges(1:end-1),r(1:end-1),'histc');

if plotflag==1
%Plot histogram
figure;
h=gca;
axes(h);
ph=bar(edges(1:end-1),r(1:end-1),'histc');
set(ph,'edgecolor',h_color,'facecolor',h_color);
end
end


