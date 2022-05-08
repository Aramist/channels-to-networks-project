function [output]=sub_sampling(input,dt)

%% Subsampling of time series data

% Arguments
%input: time series data

% Output
%output: Subsampled input data

%%
%Created on 2016
%@author: Vignayanandam R. Muddapu (CNS@IIT-Madras)

%%
[mm,nn]=size(input);
Ntime=(nn)*dt;

output=zeros(mm,Ntime);
jump = cast(1/dt, 'int32');

for te=1:Ntime
    start = 1 + jump * (te - 1);
    stop = jump + jump * (te - 1);
    output(:,te)=mean(input(:,start:stop), 2);
end

end