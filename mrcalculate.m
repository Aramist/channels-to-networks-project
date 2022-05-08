function [Rvalue,Ravg]=mrcalculate(spike_times,Nneur,Ntime)

%% Computing synchrony value across time

% Arguments
%linear_S: Spike times (linear_S=[times,number ID])
%Nneur: Number of neurons
%Ntime: Simulation time

% Output
%Rvalue: Synchrony value across time
%Ravg: Average synchrony value

% References
%Pinsky PF, Rinzel J (1995) Synchrony measures for biological neural networks. Biol Cybern 73:129–137.

%%
%Created on 2016
%@author: Vignayanandam R. Muddapu (CNS@IIT-Madras)

%%
phi=3000*ones(Nneur,Ntime-1);

for neur=1:Nneur
    temptime = find(spike_times(neur, :));
    
    for j=(1:numel(temptime) - 1)
        for i=temptime(j):1:temptime(j+1)-1
            phi(neur,i)=(2*pi*(i-temptime(j)))/(temptime(j+1)-temptime(j));
        end
    end
end
a=sqrt(-1);
tempM=sum(phi)/numel(phi);
M=exp(a*tempM);
Rvalue=((sum(exp(a*phi))/neur))./M;
Ravg=sum(abs(Rvalue))/numel(Rvalue);

end