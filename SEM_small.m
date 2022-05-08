function [snc_subsampled,stn_subsampled,gpe_subsampled,DA2,srnd,simtime]=SEM_small(dur,filename0,randinit,nolat,pd,lstsn,wsg1,wgs1,sthrsnc,wstsn)

%% Excitotoxicity Izhikevich spiking (Integrate and fire) network model

% Arguments
%x: input (decimal number)

% Output
%ff: string

%%
%Created on 2016
%@author: Vignayanandam R. Muddapu (CNS@IIT-Madras)


%% Flags
tic
time1=clock; % Tracking simulation

LFPcal=0; % LFP recording on-1, off-0

% PLOTS
mfr=1; % Computing mean firing rate 0-no, 1-yes
rastplot_GPe_STN=1;
rastplot_SNc_STN=1;
LFPplot=0;
spk_STN=1;
spk_GPe=1;
spk_SNc=1;

% Saving data
savedata=1;

%% Network Properties

disp('Initialization of the simulation ...')

% Number of neurons
%STN
nSTN=32; % (nSTNxnSTN network size)
Mstn=nSTN;
Nstn=nSTN;
Pstn=Mstn*Nstn;

%SNc
nSNc=8; % (nSNcxnSNc network size)
Msnc=nSNc;
Nsnc=nSNc;
Psnc=Msnc*Nsnc;

%GPe
nGPe=32; % (nGPexnGPe network size)
Mgpe=nGPe;
Ngpe=nGPe;
Pgpe=Mgpe*Ngpe;

% Seeding the random number
srnd = rng; % Saving randomly generated numbers
% rng(srnd); % Restoring the set of random numbers generated

%STN_SNc projections (no. of STN (projXproj) to no. (1) of SNc)
proj=4;
nproj=proj*proj;
idx = randsample(1:Pstn,Pstn);

% Neuron properties
%STN
astn=0.005; % (1/ms)
bstn=0.265; % (1/mV)
cstn=-65; % (mV)
dstn=1.5;
vpeak_stn=30;
% taustn = 0.1;

%SNc
asnc=0.0025; % (1/ms)
bsnc=0.2; % (1/mV)
csnc=-55; % (mV)
dsnc=2;
vpeak_snc=30; %(mV)
% tausnc=0.1;

%GPe
agpe=0.1; % (1/ms)
bgpe=0.2; % (1/mV)
cgpe=-65; % (mV)
dgpe=2;
vpeak_gpe=30;

% Time parameters
tau=0.1;
tausnc=tau;
taustn=tau;
taugpe=tau;
dt=tau;
tspan=dt:tau:dur;
Ttime=numel(tspan);

% Membrane capacitances
Cstn = 1; %(microF)
Csnc = 1; %(microF)
Cgpe = 1; %(microF)

% V, U initialization
if(randinit==1)
    % V, U initialization
    %STN
    Vstn = -62.5.*(rand(Mstn,Nstn)-0.5.*ones(Mstn,Nstn));
    Ustn = ((-15)-(-5)).*rand(Mstn,Nstn) + (-5);
    
    %SNc
    Vsnc = -59.47.*(rand(Msnc,Nsnc)-0.5.*ones(Msnc,Nsnc));
    Usnc = ((-15)-(-5)).*rand(Msnc,Nsnc) + (-5);
    
    %GPe
    Vgpe = -53.67.*(rand(Mgpe,Ngpe)-0.5.*ones(Mgpe,Ngpe));
    Ugpe = ((-15)-(-5)).*rand(Mgpe,Ngpe) + (-5);
    
elseif(randinit==0)
    Vstn = -65.*ones(Mstn,Nstn);
    Ustn = zeros(size(Vstn));
    % Vstns=Vstn;Ustns=Ustn;
    
    Vsnc = -59.47.*ones(Msnc,Nsnc);
    Usnc = zeros(size(Vsnc));
    % Vsncs=Vsnc;Usncs=Usnc;
    
    Vgpe = -53.67.*ones(Mgpe,Ngpe);
    Ugpe = zeros(size(Vgpe));
    % Vgpes=Vgpe;Ugpes=Ugpe;
    
end

% Currents initilization
%STN
Istn=3*ones(size(Vstn));
stn_zeros = zeros(size(Vstn));
stncurr_spk = stn_zeros;

%SNc
Isnc=8*ones(size(Vsnc));
snc_zeros = zeros(size(Vsnc));
snccurr_spk = snc_zeros;

%GPe
Igpe=4.25*ones(size(Vgpe));
gpe_zeros = zeros(size(Vgpe));
gpecurr_spk = gpe_zeros;

% psp variable initilization
%STN
h_nmdastn=stn_zeros;
h_ampastn=stn_zeros;
h_gs = stn_zeros;

%SNc
h_nmdasnc=zeros(Msnc,Nsnc);
h_ampasnc=zeros(Msnc,Nsnc);
hlat_gaba_snc=snc_zeros;

%GPe
h_nmdagpe=gpe_zeros;
h_ampagpe=gpe_zeros;
hlat_gaba_gpe=gpe_zeros;

% decay constants(ms)
taunmda=160;
tauampa=6;
taugaba=4;

% dt/T in PSP
lam_nmda = dt/taunmda;
lam_ampa = dt/tauampa;
lam_gaba= dt/taugaba;

mg0=1; % magnesium conc.

% RMP of receptors
Enmda = 0;
Eampa = 0;
Egaba = -60;

% Initial DA conc.
DA2=0.1*ones(1,Ttime);
Avg_tda=0.1*ones(1,Ttime);

% Effect of DA on post-synaptic currents
cd2=0.1;CD2=0.1;

% STN-SNc connections
wstnsnc_matrix=normrnd(wstsn,0.1,Msnc,Nsnc); % Weights are randomly assigned

% Parameters for calculating mean firing rate
windowsize=1000; %in msec
binsize=windowsize/dt;
windowtime=windowsize/1000; %in sec
nspk=1; % Type of spikes (singlet or doublet)
kkkk=10000/dt; % Check for stress after this time (in msec)
s_thrsnc=zeros(Psnc,Ttime);
taufsnc=5000;
lam_thrsnc = dt/taufsnc;
nc=zeros(1,Ttime);

%% Recording data

% Recording membrane potential of selected neurons
dVstn1=zeros(1,Ttime);dVstn2=zeros(1,Ttime);dVstn3=zeros(1,Ttime);dVstn4=zeros(1,Ttime);
dVsnc1=zeros(1,Ttime);dVsnc2=zeros(1,Ttime);dVsnc3=zeros(1,Ttime);dVsnc4=zeros(1,Ttime);
dVgpe1=zeros(1,Ttime);dVgpe2=zeros(1,Ttime);dVgpe3=zeros(1,Ttime);dVgpe4=zeros(1,Ttime);

if(LFPcal==1)
    % LFP
    % STN
    LFP_STN_exc=zeros(1,Ttime);
    LFP_STN_inh=zeros(1,Ttime);
    LFP_STN_tot=zeros(1,Ttime);
    
    % SNc
    LFP_SNc_exc=zeros(1,Ttime);
    LFP_SNc_inh=zeros(1,Ttime);
    LFP_SNc_tot=zeros(1,Ttime);
    
    % GPe
    LFP_GPe_exc=zeros(1,Ttime);
    LFP_GPe_inh=zeros(1,Ttime);
    LFP_GPe_tot=zeros(1,Ttime);
    
end

% Recording the mean firing rates
mfrstn=zeros(Pstn,Ttime);
mfrstn_cache = zeros(Pstn, Ttime + 1);
mfrsnc=zeros(Psnc,Ttime);
mfrsnc_cache = zeros(Psnc, Ttime + 1);
mfrgpe_cache = zeros(Pgpe, Ttime + 1);

% Others
rndpick_array=[];

%% SIMULATION
disp('Simulation started ...')

for k = 1:Ttime
    
    %-----------------------DA modulating lateral STN connections------------------------%
    
    %No Dopamine condition
    %     da=0;
    %     DA(k)=0;
    
    %Dopamine condition
    DA=Avg_tda(k);
    
    %Computing lateral connections
    da=0.5;
    wda_snc=1;wda_stn=1;wda_gpe=1;wda_snclat=1;
    wda_snc = (1-(cd2*DA))*wda_snc;
    wlatsnc = nolat.*weightcal_snc(DA);
    wlatstn = nolat.*weightcal_stn(DA*pd);
    wlatgpe = nolat.*weightcal_gpe(da);
    
    %Computing DA effect on post-synaptic current
    wsg=((1-CD2*da))*wsg1;
    wgs=((1-CD2*DA))*wgs1;
    
    %----------------------------------------STN-----------------------------------------%
    %---------------------------Input from stn to stn(laterals)--------------------------%
    
    % psp variable
    h_nmdastn = (1-lam_nmda).* h_nmdastn + lam_nmda.*stncurr_spk; %psp nmda stn lat
    h_ampastn = (1-lam_ampa).* h_ampastn + lam_ampa.*stncurr_spk; %psp ampa stn lat
    
    tmplat_nmda_stn = wda_stn.*h_nmdastn.*(Enmda - Vstn);
    tmplat_ampa_stn = wda_stn.*h_ampastn.*(Eampa - Vstn);
    
    
    Ilat_nmda_stn = conv2(tmplat_nmda_stn, wlatstn, 'same'); % lat nmda stn
    Ilat_ampa_stn = conv2(tmplat_ampa_stn, wlatstn, 'same');% lat ampa stn
    
    B = 1./(1 + (mg0/3.57).*exp(-0.062.*Vstn));
    
    %----------------------------------------Input from gpe to stn---------------------%
    % psp variable
    h_gs = (1-lam_gaba).* h_gs + lam_gaba.*gpecurr_spk;% input from gpe to nmda stn
    % gaba current
    I_gs = wgs.*h_gs.*(Egaba - Vstn);
    
    % total current stn recieves
    I_ss=B.*Ilat_nmda_stn + Ilat_ampa_stn;
    Itmpstn =  I_ss + I_gs; % total currents
    % Itmpstn =  I_gs; % w/o laterals total currents
    
    if(LFPcal==1)
        % LFP
        LFP_STN_exc(k)=sum(sum(I_ss));
        LFP_STN_inh(k)=sum(sum(I_gs));
        LFP_STN_tot(k)=LFP_STN_exc(k)+LFP_STN_inh(k);
    end
    
    % V,U updated
    dvstn = taustn.*(((0.04.*Vstn.*Vstn)+5.*Vstn - Ustn + Itmpstn +140 + Istn)./Cstn);%+wgn(Mstn,Nstn,lnoise_st,imp_st);
    dustn = taustn.*astn.*(bstn.*(Vstn) - Ustn);
    Vstn_nxt = Vstn + dvstn;
    Ustn_nxt = Ustn + dustn;
    
    inds = find(Vstn_nxt > vpeak_stn);
    
    Vstn_nxt(inds) = cstn.*ones(size(inds));
    Ustn_nxt(inds) = Ustn(inds) + dstn.*ones(size(inds));
    stncurr_spk = stn_zeros;
    stncurr_spk(inds) = ones(size(inds));
    
    Vstn = Vstn_nxt;
    Ustn = Ustn_nxt;
    
    %----------------------------------------GPe-----------------------------------------%
    %--------------------- Input from stn to gpe----------------------------------%
    % psp variable
    h_nmdagpe = (1-lam_nmda).* h_nmdagpe + lam_nmda.*stncurr_spk;
    h_ampagpe = (1-lam_ampa).* h_ampagpe + lam_ampa.*stncurr_spk;
    
    % nmda and ampa currents
    Inmdagpe = wsg.*h_nmdagpe.*(Enmda - Vgpe);
    Iampagpe = wsg.*h_ampagpe.*(Eampa - Vgpe);
    
    %-------------------------- Input from gpe to gpe(laterals)-------------------%
    % psp variable
    hlat_gaba_gpe = (1-lam_gaba).* hlat_gaba_gpe + lam_gaba.*gpecurr_spk;
    tmplat_gaba_gpe = wda_gpe.*hlat_gaba_gpe.*(Egaba - Vgpe);
    
    % lateral currents
    Ilatgpe = conv2(tmplat_gaba_gpe, wlatgpe, 'same');
    
    B = 1./(1 + (mg0/3.57).*exp(-0.062.*Vgpe));
    I_sg=B.*Inmdagpe+Iampagpe;
    Itmpgpe = Ilatgpe + I_sg;
    
    if(LFPcal==1)
        %LFP
        LFP_GPe_exc(k)=sum(sum(I_sg));
        LFP_GPe_inh(k)=sum(sum(Ilatgpe));
        LFP_GPe_tot(k)=LFP_GPe_exc(k)+LFP_GPe_inh(k);
    end
    
    % V ,U updated
    dvgpe = taugpe.*((0.04.*Vgpe.*Vgpe)+5.*Vgpe+140 - Ugpe + Itmpgpe+Igpe)./Cgpe;
    dugpe = taugpe.*agpe.*(bgpe.*(Vgpe) - Ugpe);
    Vgpe_nxt = Vgpe + dvgpe;
    Ugpe_nxt = Ugpe + dugpe;
    
    inds = find(Vgpe_nxt > vpeak_gpe);
    
    Vgpe_nxt(inds) = cgpe.*ones(size(inds));
    Ugpe_nxt(inds) = Ugpe(inds) + dgpe.*ones(size(inds));
    gpecurr_spk = gpe_zeros;
    gpecurr_spk(inds) = ones(size(inds));
    
    Vgpe = Vgpe_nxt;
    Ugpe = Ugpe_nxt;
    
    %----------------------------------------SNc-----------------------------------------%
    %------------------------------STN projections to SNc--------------------------------%
    
    % PATTERN -> Random picking of STN (4 X 4) projects to single SNc
    start=1;stop=nproj;
    totspk=zeros(1,Psnc);
    for ll=1:Psnc
        % stn_snccurr1=zeros(proj,proj);
        stn_snccurr1=reshape(stncurr_spk(idx(start:stop)),proj,proj);
        totspk(ll)=sum(sum(stn_snccurr1));
        start=start+nproj;stop=stop+nproj;
    end
    stn_snccurr=reshape(totspk,Msnc,Nsnc);
    
    h_nmdasnc = (1-lam_nmda).* h_nmdasnc + lam_nmda.*stn_snccurr; %psp nmda snc
    h_ampasnc = (1-lam_ampa).* h_ampasnc + lam_ampa.*stn_snccurr;
    
    tmp_nmda_snc = h_nmdasnc.*(Enmda - Vsnc);
    tmp_ampa_snc = h_ampasnc.*(Eampa - Vsnc);
    
    I_nmda_snc = lstsn.*wda_snc.*wstnsnc_matrix.*tmp_nmda_snc;
    I_ampa_snc = lstsn.*wda_snc.*wstnsnc_matrix.*tmp_ampa_snc;
    
    B = 1./(1 + (mg0/3.57).*exp(-0.062.*Vsnc));
    
    Istnsnc=B.*I_nmda_snc + I_ampa_snc;
    
    %---------------------------Input from snc to snc(laterals)--------------------------%
    % psp variable
    hlat_gaba_snc = (1-lam_gaba).* hlat_gaba_snc + lam_gaba.*snccurr_spk;
    tmplat_gaba_snc = wda_snclat.*hlat_gaba_snc.*(Egaba - Vsnc);
    
    % lateral currents
    Ilatsnc = conv2(tmplat_gaba_snc, wlatsnc, 'same');
    
    % total current snc recieves
    Itmpsnc = Istnsnc + Ilatsnc; % total currents
    
    if(LFPcal==1)
        %LFP
        LFP_SNc_exc(k)=sum(sum(Istnsnc));
        LFP_SNc_inh(k)=sum(sum(Ilatsnc));
        LFP_SNc_tot(k)=LFP_SNc_exc(k)+LFP_SNc_inh(k);
    end
    
    % V,U updated
    dvsnc = tausnc.*(((0.04.*Vsnc.*Vsnc)+5.*Vsnc - Usnc + Itmpsnc +140 + Isnc)./Csnc);%+wgn(Msnc,Nsnc,lnoise_sn,imp_sn);
    dusnc = tausnc.*asnc.*(bsnc.*(Vsnc) - Usnc);
    Vsnc_nxt = Vsnc + dvsnc;
    Usnc_nxt = Usnc + dusnc;
    
    inds = find(Vsnc_nxt > vpeak_snc);
    
    Vsnc_nxt(inds) = csnc.*ones(size(inds));
    Usnc_nxt(inds) = Usnc(inds) + dsnc.*ones(size(inds));
    snccurr_spk = snc_zeros;
    snccurr_spk(inds) = ones(size(inds));
    
    % DA signal from SNc (spatial integration and normalized)
    Avg_tda(k+1)=sum(snccurr_spk, 'all')/(Msnc*Nsnc);
    DA2(k)=Avg_tda(k);
    
    % Computing mean firing rates
    if mfr==1
        
        % Saving spike train for each trial (STN)
        tmps1 = reshape(stncurr_spk, Pstn, 1);
        mfrstn_cache(:, k+1) = mfrstn_cache(:, k) + tmps1;
        
        % Saving spike train for each trial (SNc)
        tmps1 = reshape(snccurr_spk, Psnc, 1);
        mfrsnc_cache(:, k+1) = mfrsnc_cache(:, k) + tmps1;

        mfrgpe_cache(:, k+1) = mfrgpe_cache(:, k) + reshape(gpecurr_spk, Pgpe, 1);
        
        % Calculating Mean firing rate with fixed sliding window
        start=k-(binsize);
        stop=k;
        if start<1
            start=1;
        end
        
        mfrstn(:,k)=mfrwindow(mfrstn_cache,start,stop,windowtime,nspk);
        mfrsnc(:,k)=mfrwindow(mfrsnc_cache,start,stop,windowtime,nspk);
        mfrsncs=mfrsnc(:,k);
        
        % Calculating stress variable as a function of mean firing rate
        s_thrsnc(:,k+1) = (1-lam_thrsnc).*s_thrsnc(:,k) + lam_thrsnc.*mfrsncs;
        
        scel=[];
        if k>kkkk
            scel=find(s_thrsnc(:,k+1)>sthrsnc);
            if isempty(scel)==0
                rndpick_array=[rndpick_array scel'];
            end
        end
        
        rndpick_array=unique(rndpick_array);
        
        Vsnc_nxt(rndpick_array) = 0*Vsnc_nxt(rndpick_array);
        Usnc_nxt(rndpick_array) = 0*Usnc_nxt(rndpick_array);
        
        nc(k)=numel(rndpick_array);
        if (k > 1) && (nc(k) >= Psnc / 2) && (nc(k-1) < Psnc/2)
            disp(['TD50: ', num2str(k*dt), 'ms.'])
            break
        end
        
    end
    
    Vsnc = Vsnc_nxt;
    Usnc = Usnc_nxt;
    
    
    dVstn1(k)=Vstn_nxt(17,17);
    dVstn2(k)=Vstn_nxt(17,12);
    dVstn3(k)=Vstn_nxt(12,17);
    dVstn4(k)=Vstn_nxt(5,5);
    
    dVgpe1(k)=Vgpe_nxt(17,17);
    dVgpe2(k)=Vgpe_nxt(17,12);
    dVgpe3(k)=Vgpe_nxt(12,17);
    dVgpe4(k)=Vgpe_nxt(5,5);
    
    dVsnc1(k)=Vsnc_nxt(1,1);
    dVsnc2(k)=Vsnc_nxt(3,3);
    dVsnc3(k)=Vsnc_nxt(6,6);
    dVsnc4(k)=Vsnc_nxt(8,8);
    
    if mod(k, 1000) == 0
        disp(['Simulation progress: ', num2str(k*dt), 'ms.']);
    end
    
end
disp('Simulation ended ...')
toc

%% Processing data
if false
    tic
    disp('Subsampling started ...')
    % Sub-sampling
    [stn_subsampled]=sub_sampling_firings(mfrstn_cache,Pstn,Ttime,dt);
    [gpe_subsampled]=sub_sampling_firings(mfrgpe_cache,Pgpe,Ttime,dt);
    [snc_subsampled]=sub_sampling_firings(mfrsnc_cache,Psnc,Ttime,dt);
    [mfrstn_subsampled]=sub_sampling(mfrstn,dt);
    [mfrsnc_subsampled]=sub_sampling(mfrsnc,dt);
    s_thrsnc(:,end)=[];
    [s_thrsnc_subsampled]=sub_sampling(s_thrsnc,dt);
    disp('Subsampling ended ...')
    toc
    
    clear mfrstn
    clear mfrsnc
    clear mfrstn_cache
    clear mfrgpe_cache
    clear mfrsnc_cache
    clear s_thrsnc
    
    %% Saving data without analysis
    if savedata==1
        if mfr==1
            kid=((Psnc)-nc);
            save(filename0,'stn_subsampled','gpe_subsampled','snc_subsampled','DA2','kid','mfrstn_subsampled','mfrsnc_subsampled','s_thrsnc_subsampled','srnd')
        else
            save(filename0,'stn_subsampled','gpe_subsampled','snc_subsampled','DA2','srnd');
        end
    end
    tic
    
    %% Analysis
    
    disp('Analysis started ...')
    [Rvalstn,Ravgstn]=mrcalculate(stn_subsampled,Pstn,Ttime*dt);
    [Rvalgpe,Ravggpe]=mrcalculate(gpe_subsampled,Pgpe,Ttime*dt);
    [Rvalsnc,Ravgsnc]=mrcalculate(snc_subsampled,Psnc,Ttime*dt);
    [Rvalstngpe,Ravgstngpe]=mrintercalculate(stn_subsampled,gpe_subsampled,Pstn,Ttime*dt);
    
    %Frequency
    % Individual neurons
    fstn1=sum((dVstn1)>15)/(Ttime.*dt.*1e-3);
    fstn2=sum((dVstn2)>15)/(Ttime.*dt.*1e-3);
    fstn3=sum((dVstn3)>15)/(Ttime.*dt.*1e-3);
    fstn4=sum((dVstn4)>15)/(Ttime.*dt.*1e-3);
    
    fgpe1=sum((dVgpe1)>15)/(Ttime.*dt.*1e-3);
    fgpe2=sum((dVgpe2)>15)/(Ttime.*dt.*1e-3);
    fgpe3=sum((dVgpe3)>15)/(Ttime.*dt.*1e-3);
    fgpe4=sum((dVgpe4)>15)/(Ttime.*dt.*1e-3);
    
    fsnc1=sum((dVsnc1)>15)/(Ttime.*dt.*1e-3);
    fsnc2=sum((dVsnc2)>15)/(Ttime.*dt.*1e-3);
    fsnc3=sum((dVsnc3)>15)/(Ttime.*dt.*1e-3);
    fsnc4=sum((dVsnc4)>15)/(Ttime.*dt.*1e-3);
    
    % Whole network
    base=Pstn/2;
    stnfrequency=sum(stn_subsampled, 'all')/(2*base.*dt.*Ttime.*1e-3);
    
    base=Pgpe/2;
    gpefrequency=sum(gpe_subsampled, 'all')/(2*base.*dt.*Ttime.*1e-3);
    
    base=Psnc/2;
    sncfrequency=sum(snc_subsampled, 'all')/(2*base.*dt.*Ttime.*1e-3)/2;
    disp('Analysis ended ...')
    
    simtime=[];
    time2=clock;
    simtime(1)=time1(1);simtime(2)=time1(2);simtime(3)=time1(3);simtime(4)=time1(4);simtime(5)=time1(5);
    simtime(6)=time2(1);simtime(7)=time2(2);simtime(8)=time2(3);simtime(9)=time2(4);simtime(10)=time2(5);
    if savedata==1
        if mfr==1
            save(filename0,'stn_subsampled','gpe_subsampled','snc_subsampled','DA2', ...
                'kid','mfrstn_subsampled','mfrsnc_subsampled','s_thrsnc_subsampled','Rvalstn','Ravgstn', ...
                'Rvalsnc','Ravgsnc','Rvalgpe','Ravggpe','Rvalstngpe', ...
                'Ravgstngpe','srnd','simtime', 'dVstn1', 'dVstn2', 'dVstn3', ...
                'dVstn4', 'dVsnc1', 'dVsnc2', 'dVsnc3', 'dVsnc4', 'dVgpe1', ...
                'dVgpe2', 'dVgpe3', 'dVgpe4');
        else
            save(filename0,'stn_subsampled','gpe_subsampled','snc_subsampled','DA2', ...
                'Rvalstn','Ravgstn','Rvalsnc','Ravgsnc','Rvalgpe','Ravggpe', ...
                'Rvalstngpe','Ravgstngpe','srnd','simtime', 'dVstn1', 'dVstn2', ...
                'dVstn3', 'dVstn4', 'dVsnc1', 'dVsnc2', 'dVsnc3', 'dVsnc4', ...
                'dVgpe1', 'dVgpe2', 'dVgpe3', 'dVgpe4');
        end
    end
    toc
end