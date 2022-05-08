function [mfr]=mfrwindow(matrix,start,stop,windowtime,nspk)

%% Computing instantaenous mean firing rate (with windows size (stop-start))

% Arguments
%matrix: Spike train
%start: Start of the window
%stop: Stop of the window
%windowtime: Time of window size given
%nspk: Spiklet (singlet or doublet or triplet)

% Output
%mfr: mean firing rate

%%
%Created on 2016
%@author: Vignayanandam R. Muddapu (CNS@IIT-Madras)

%%
% mfr=((sum(matrix(:,start:stop),2)/windowtime)/nspk);
sum_over_cols = matrix(:, stop + 1) - matrix(:, start + 1);
mfr = (sum_over_cols / windowtime) / nspk;
end