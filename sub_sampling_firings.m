function [reduced_len]=sub_sampling_firings(fr_cache,Nneur,Ttime,dt)

%% Subsampling of spike time data

% Arguments
%linear_S: Spike times (linear_S=[times,number ID])
%Nneur: Number of neurons
%Ttime: Simulation time
%dt: TIme step of simulation

% Output
%temptime2: Subsampled input data

%%
%Created on 2016
%@author: Vignayanandam R. Muddapu (CNS@IIT-Madras)

%%
Ntime=(Ttime)*dt;
% Basically reduces the time-bin resolution to 1ms
reduced_len = zeros(Nneur, Ntime);
% The input is now expected to be the spike hist cache matrix. It has shape
% (n_neurons, n_timesteps) and is cumulative in the spiking activity.


for te=1:Ntime
    start = 1 + 10 * (te-1) + 1;
    stop = 10 + 10 * (te-1) + 1;
    % will index neurons for which at least one spike occurred between
    % start and stop
    condition = (fr_cache(:, stop) - fr_cache(:, start)) > 0;
    reduced_len(condition, te) = 1;
    reduced_len(~condition, te) = 0;
end

