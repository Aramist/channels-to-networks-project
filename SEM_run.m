%% Spiking Excitotoxicity Model (SEM)
%Script for running the simulation

%%
%Created on 2016
%@author: Vignayanandam R. Muddapu (CNS@IIT-Madras)

%References
%Muddapu VR, Mandali A, Chakravarthy VS, Ramaswamy S (2019) A computational model of loss of dopaminergic cells in Parkinson’s disease
%due to glutamate-induced excitotoxicity. Front Neural Circuits 13:11
%Available at: https://www.frontiersin.org/articles/10.3389/fncir.2019.00011/abstract [Accessed February 25, 2019].

%%
clc;clear;close all;
tic
time=clock; time=fix(time);
curdate=time(3);curmonth=time(2); hour=time(4);minute=time(5);sec=time(6); % Tracking time

dur=50000; % Duration of simulation in milliseconds

randinit=1; % 0-No random initializatiion, 1-Random initializatiion
nolat=1; % 0-Lateral connections off, 1-Lateral connections on
pd=0; % 0-PD condition on, 1-PD condition off
lstsn=1; % 0-STN-SNc connection off, 1-STN-SNc connection on

wsg=1; % Connection strength from STN to GPe
wgs=20; % Connection strength from GPe to STN

sthrsnc=11; % Stress thershold
sthrsnc1=deci2str(sthrsnc);

% Small, 64 neurons: 11.229s
% Smallish, 100 neurons: 10.870s
% Medium, 144 neurons: 12.444s
% Large, 400 neurons


% High thresh:
% 64 neurons: 10.575s
% 100 neurons: 14.171
% 144 neurons: 16.730s
% 256 neurons: 31.107s

wstsn=1; % Connection strength from STN to SNc
wstsn1=deci2str(wstsn);

filename0=strcat( ...
    'SEM_Wstsn',num2str(wstsn1), ...
    '_Sthr',num2str(sthrsnc1), ...
    '_pd', num2str(pd), ...
    '_',num2str(dur),'msec_', ...
    num2str(curdate),'_',num2str(curmonth), ...
    '_', num2str(hour), '_',num2str(minute), '_', ...
    num2str(sec));

[snc_firings,stn_firings,gpe_firings,DA,srnd,simtime]=SEM_small(dur,filename0,randinit,nolat,pd,lstsn,wsg,wgs,sthrsnc,wstsn);

toc