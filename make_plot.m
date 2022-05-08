clear all;
load('/Users/aramis/Documents/MATLAB/small_pd_off_test.mat','stn_subsampled','gpe_subsampled','snc_subsampled','DA2', ...
    'kid','mfrstn_subsampled','mfrsnc_subsampled','s_thrsnc_subsampled','Rvalstn','Ravgstn', ...
    'Rvalsnc','Ravgsnc','Rvalgpe','Ravggpe','Rvalstngpe', ...
    'Ravgstngpe','srnd','simtime', 'dVstn1', 'dVstn2', 'dVstn3', ...
    'dVstn4', 'dVsnc1', 'dVsnc2', 'dVsnc3', 'dVsnc4', 'dVgpe1', ...
    'dVgpe2', 'dVgpe3', 'dVgpe4');

filename0='plot_small';

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

tau=0.1;
dt=tau;
dur=25000;
tspan=dt:tau:dur;
Ttime=numel(tspan);
nc=zeros(1,Ttime);

mfr=1; % Computing mean firing rate 0-no, 1-yes
rastplot_GPe_STN=0;
rastplot_SNc_STN=0;
LFPplot=0;
spk_STN=0;
spk_GPe=0;
spk_SNc=0;


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



base=Pstn/2;
stnfrequency=sum(stn_subsampled, 'all')/(2*base.*dt.*Ttime.*1e-3);

base=Pgpe/2;
gpefrequency=sum(gpe_subsampled, 'all')/(2*base.*dt.*Ttime.*1e-3);

base=Psnc/2;
sncfrequency=sum(snc_subsampled, 'all')/(2*base.*dt.*Ttime.*1e-3)/2;

n_plotted_neurons=64

if(rastplot_GPe_STN==1)
    sec=0.001;
    sizz=10;
    fig0=figure(110);
    set(fig0, 'Position', [5, 50, 1920, 955]);
    % set(fig1, 'Position', [50, 100, 2000, 920]);
    subplot(511)
    set(gca,'fontsize',sizz);
    x_axis = (1:Ttime*dt) * sec;
    used_neurons = randsample(Pstn, n_plotted_neurons);
    for neuron=(1:n_plotted_neurons)
        plot(x_axis(stn_subsampled(used_neurons(neuron), :)==1), neuron, '.b', 'MarkerSize', 2);
        hold on;
    end
    % plot(sec*(stn_firings2(:,1)),stn_firings2(:,2),'.','MarkerSize',12);
    xlim([0 sec*dt*Ttime]);
    ylim([0 n_plotted_neurons]);
    tit1=strcat('STN overall network freq = ',num2str(stnfrequency),'Hz');
    title(tit1)
    %     title('STN firing')
    % xlabel('Time (msec)')
    ylabel('# of neurons')
    subplot(512)
    set(gca,'fontsize',sizz);
    plot(sec*(1:numel(Rvalstn)),abs(Rvalstn))
    xlim([0 sec*dt*Ttime]);
    ylim([0 1.2]);
    title('STN Rsyn')
    refline([0 Ravgstn]);
    fh=strcat(num2str(Ravgstn));legend(fh);
    % xlabel('Time (msec)')
    ylabel('R syn')
    subplot(513)
    set(gca,'fontsize',sizz);
    used_neurons = randsample(Pgpe, n_plotted_neurons);
    for neuron=(1:n_plotted_neurons)
        plot(x_axis(gpe_subsampled(used_neurons(neuron), :)==1), neuron, '.b', 'MarkerSize', 2);
        hold on;
    end
    % plot(sec*(gpe_firings2(:,1)),gpe_firings2(:,2),'.','MarkerSize',12);
    xlim([0 sec*dt*Ttime]);
    ylim([0 n_plotted_neurons]);
    tit1=strcat('GPe overall network freq = ',num2str(gpefrequency),'Hz');
    title(tit1)
    %     title('GPe firing')
    % xlabel('Time (msec)')
    ylabel('# of neurons')
    subplot(514)
    set(gca,'fontsize',sizz);
    plot(sec*(1:numel(Rvalgpe)),abs(Rvalgpe))
    xlim([0 sec*dt*Ttime]);
    ylim([0 1.2]);
    title('GPe Rsyn')
    refline([0 Ravggpe]);
    fh1=strcat(num2str(Ravggpe));legend(fh1);
    % xlabel('Time (sec)')
    ylabel('R syn')
    subplot(515)
    set(gca,'fontsize',sizz);
    plot(sec*(1:numel(Rvalstngpe)),abs(Rvalstngpe));
    xlim([0 sec*dt*Ttime]);
    ylim([0 1.2]);
    refline([0 Ravgstngpe]);
    fh3=strcat(num2str(Ravgstngpe));legend(fh3);
    title('STN-GPe Rsyn')
    xlabel('Time (sec)')
    ylabel('R syn')
    f0=strcat('GPe_STN_',filename0);
    saveas(fig0,f0,'png');
    %     saveas(fig0,filename0,'png');
    clear fig0
end

if(rastplot_SNc_STN==1)
    sec=0.001;
    sizz=10;
    fig1=figure(111);
    set(fig1, 'Position', [5, 50, 1920, 955]);
    % set(fig1, 'Position', [50, 100, 2000, 920]);
    subplot(411)
    set(gca,'fontsize',sizz);
    x_axis = (1:Ttime*dt) * sec;
    used_neurons = randsample(Pstn, n_plotted_neurons);
    for neuron=(1:n_plotted_neurons)
        plot(x_axis(stn_subsampled(used_neurons(neuron), :)==1), neuron, '.b', 'MarkerSize', 2);
        hold on;
    end
    % plot(sec*(stn_firings2(:,1)),stn_firings2(:,2),'.','MarkerSize',12);
    xlim([0 sec*dt*Ttime]);
    ylim([0 n_plotted_neurons]);
    tit1=strcat('STN overall network freq = ',num2str(stnfrequency),'Hz');
    title(tit1)
    % xlabel('Time (msec)')
    ylabel('# of neurons')
    subplot(412)
    set(gca,'fontsize',sizz);
    plot(sec*(1:numel(Rvalstn)),abs(Rvalstn))
    xlim([0 sec*dt*Ttime]);
    ylim([0 1.2]);
    title('STN Rsyn')
    % xlabel('Time (msec)')
    ylabel('R syn')
    refline([0 Ravgstn]);
    fh=strcat(num2str(Ravgstn));legend(fh);
    subplot(413)
    set(gca,'fontsize',sizz);
    used_neurons = randsample(Psnc, n_plotted_neurons);
    for neuron=(1:n_plotted_neurons)
        plot(x_axis(snc_subsampled(used_neurons(neuron), :)==1), neuron, '.b', 'MarkerSize', 2);
        hold on;
    end
    % plot(sec*(snc_firings2(:,1)),snc_firings2(:,2),'.','MarkerSize',12);
    xlim([0 sec*dt*Ttime]);
    ylim([0 n_plotted_neurons]);
    tit2=strcat('SNc overall network freq = ',num2str(sncfrequency),'Hz');
    title(tit2)
    % xlabel('Time (msec)')
    ylabel('# of neurons')
    subplot(414)
    set(gca,'fontsize',sizz);
    plot(sec*(1:numel(Rvalsnc)),abs(Rvalsnc))
    xlim([0 sec*dt*Ttime]);
    ylim([0 1.2]);
    title('SNc Rsyn')
    xlabel('Time (sec)')
    ylabel('R syn')
    refline([0 Ravgsnc]);
    fh2=strcat(num2str(Ravgsnc));legend(fh2);
    f1=strcat('SNc_STN_',filename0);
    saveas(fig1,f1,'png');
    %     saveas(fig3,filename0,'png');
    clear fig1
end

if(LFPplot==1)
    sec=0.001;
    sizz=10;
    fig2=figure(112);
    set(fig2, 'Position', [5, 50, 1920, 955]);
    % set(fig1, 'Position', [50, 100, 2000, 920]);
    subplot(311)
    set(gca,'fontsize',sizz);
    plot(sec*dt*(1:Ttime),LFP_STN_exc,'b','DisplayName','STN_{exc}');hold on;
    plot(sec*dt*(1:Ttime),LFP_STN_inh,'r','DisplayName','STN_{inh}');
    plot(sec*dt*(1:Ttime),LFP_STN_tot,'k','DisplayName','STN_{tot}');hold off;
    xlim([0 sec*dt*Ttime]);
    %     ylim([-5 500]);
    title('STN LFP')
    % xlabel('Time (msec)')
    legend('show');
    %     ylabel('# of neurons')
    subplot(312)
    set(gca,'fontsize',sizz);
    plot(sec*dt*(1:Ttime),LFP_GPe_exc,'b','DisplayName','GPe_{exc}');hold on;
    plot(sec*dt*(1:Ttime),LFP_GPe_inh,'r','DisplayName','GPe_{inh}');
    plot(sec*dt*(1:Ttime),LFP_GPe_tot,'k','DisplayName','GPe_{tot}');hold off;
    xlim([0 sec*dt*Ttime]);
    %     ylim([-5 500]);
    title('GPe LFP')
    % xlabel('Time (msec)')
    legend('show');
    %     ylabel('# of neurons')
    subplot(313)
    set(gca,'fontsize',sizz);
    plot(sec*dt*(1:Ttime),LFP_SNc_exc,'b','DisplayName','SNc_{exc}');hold on;
    plot(sec*dt*(1:Ttime),LFP_SNc_inh,'r','DisplayName','SNc_{inh}');
    plot(sec*dt*(1:Ttime),LFP_SNc_tot,'k','DisplayName','SNc_{tot}');hold off;
    xlim([0 sec*dt*Ttime]);
    %     ylim([-5 500]);
    title('SNc LFP')
    % xlabel('Time (msec)')
    legend('show');
    %     ylabel('# of neurons')
    f2=strcat('LFP_',filename0);
    saveas(fig2,f2,'png');
    clear fig2
end

if(spk_STN==1)
    sec=0.001;
    fig3=figure(113);
    set(fig3, 'Position', [5, 50, 1920, 955]);
    main1=strcat('STN overall network freq = ',num2str(stnfrequency),'Hz');
    sgtitle(main1)
    subplot(411)
    plot(sec*dt*(1:numel(dVstn1)),dVstn1,'linewidth',1)
    ylabel('Voltage (mV)')
    tit=strcat('STN (17,17) (STN freq.= ',num2str(fstn1),' Hz)');
    title(tit)
    xlim([0 sec*dt*Ttime]);
    ylim([-90 50]);
    %     xlabel('Time (sec)')
    subplot(412)
    plot(sec*dt*(1:numel(dVstn2)),dVstn2,'linewidth',1)
    ylabel('Voltage (mV)')
    tit=strcat('STN (17,12) (STN freq.= ',num2str(fstn2),' Hz)');
    title(tit)
    xlim([0 sec*dt*Ttime]);
    ylim([-90 50]);
    %     xlabel('Time (sec)')
    subplot(413)
    plot(sec*dt*(1:numel(dVstn3)),dVstn3,'linewidth',1)
    ylabel('Voltage (mV)')
    tit=strcat('STN (12,17) (STN freq.= ',num2str(fstn3),' Hz)');
    title(tit)
    xlim([0 sec*dt*Ttime]);
    ylim([-90 50]);
    %     xlabel('Time (sec)')
    subplot(414)
    plot(sec*dt*(1:numel(dVstn4)),dVstn4,'linewidth',1)
    ylabel('Voltage (mV)')
    tit=strcat('STN (5,5) (STN freq.= ',num2str(fstn4),' Hz)');
    title(tit)
    xlim([0 sec*dt*Ttime]);
    ylim([-90 50]);
    xlabel('Time (sec)')
    f3=strcat('zV_STN_',filename0);
    saveas(fig3,f3,'png');
end

if(spk_GPe==1)
    sec=0.001;
    fig4=figure(114);
    set(fig4, 'Position', [5, 50, 1920, 955]);
    main1=strcat('GPe overall network freq = ',num2str(gpefrequency),'Hz');
    sgtitle(main1)
    subplot(411)
    plot(sec*dt*(1:numel(dVgpe1)),dVgpe1,'linewidth',1)
    ylabel('Voltage (mV)')
    tit=strcat('GPe (17,17) (GPe freq.= ',num2str(fgpe1),' Hz)');
    title(tit)
    xlim([0 sec*dt*Ttime]);
    ylim([-90 50]);
    %     xlabel('Time (sec)')
    subplot(412)
    plot(sec*dt*(1:numel(dVgpe2)),dVgpe2,'linewidth',1)
    ylabel('Voltage (mV)')
    tit=strcat('GPe (17,12) (GPe freq.= ',num2str(fgpe2),' Hz)');
    title(tit)
    xlim([0 sec*dt*Ttime]);
    ylim([-90 50]);
    %     xlabel('Time (sec)')
    subplot(413)
    plot(sec*dt*(1:numel(dVgpe3)),dVgpe3,'linewidth',1)
    ylabel('Voltage (mV)')
    tit=strcat('GPe (12,17) (GPe freq.= ',num2str(fgpe3),' Hz)');
    title(tit)
    xlim([0 sec*dt*Ttime]);
    ylim([-90 50]);
    %     xlabel('Time (sec)')
    subplot(414)
    plot(sec*dt*(1:numel(dVgpe4)),dVgpe4,'linewidth',1)
    ylabel('Voltage (mV)')
    tit=strcat('GPe (5,5) (GPe freq.= ',num2str(fgpe4),' Hz)');
    title(tit)
    xlim([0 sec*dt*Ttime]);
    ylim([-90 50]);
    xlabel('Time (sec)')
    f4=strcat('zV_GPe_',filename0);
    saveas(fig4,f4,'png');
end

if(spk_SNc==1)
    sec=0.001;
    fig5=figure(115);
    set(fig5, 'Position', [5, 50, 1920, 955]);
    main1=strcat('SNc overall network freq = ',num2str(sncfrequency),'Hz');
    sgtitle(main1)
    subplot(411)
    plot(sec*dt*(1:numel(dVsnc1)),dVsnc1,'linewidth',1)
    ylabel('Voltage (mV)')
    tit=strcat('SNc (1,1) (SNc freq.= ',num2str(fsnc1),' Hz)');
    title(tit)
    xlim([0 sec*dt*Ttime]);
    ylim([-90 50]);
    %     xlabel('Time (sec)')
    subplot(412)
    plot(sec*dt*(1:numel(dVsnc2)),dVsnc2,'linewidth',1)
    ylabel('Voltage (mV)')
    tit=strcat('SNc (3,3) (SNc freq.= ',num2str(fsnc2),' Hz)');
    title(tit)
    xlim([0 sec*dt*Ttime]);
    ylim([-90 50]);
    %     xlabel('Time (sec)')
    subplot(413)
    plot(sec*dt*(1:numel(dVsnc3)),dVsnc3,'linewidth',1)
    ylabel('Voltage (mV)')
    tit=strcat('SNc (6,6) (SNc freq.= ',num2str(fsnc3),' Hz)');
    title(tit)
    xlim([0 sec*dt*Ttime]);
    ylim([-90 50]);
    %     xlabel('Time (sec)')
    subplot(414)
    plot(sec*dt*(1:numel(dVsnc4)),dVsnc4,'linewidth',1)
    ylabel('Voltage (mV)')
    tit=strcat('SNc (8,8) (SNc freq.= ',num2str(fsnc4),' Hz)');
    title(tit)
    xlim([0 sec*dt*Ttime]);
    ylim([-90 50]);
    xlabel('Time (sec)')
    f5=strcat('zV_SNc_',filename0);
    saveas(fig5,f5,'png');
end

if mfr==1
    plotflag=0;
    if(plotflag == 1)
        kid=((Psnc)-nc);
        sec=0.001;
        sizz=10;
        fig1=figure(2);
        set(fig1, 'Position', [5, 50, 1920, 955]);
        % set(fig1, 'Position', [50, 100, 2000, 920]);
        subplot(511)
        set(gca,'fontsize',sizz);
        plot(sec*(stn_firings2(:,1)),stn_firings2(:,2),'.','MarkerSize',12);
        xlim([0 sec*dt*Ttime]);
        ylim([0 Pstn]);
        title('STN firing')
        % xlabel('Time (msec)')
        ylabel('# of neurons')
        subplot(512)
        set(gca,'fontsize',sizz);
        plot(sec*(1:numel(Rvalstn)),abs(Rvalstn))
        xlim([0 sec*dt*Ttime]);
        ylim([0 1.2]);
        title('STN Rsyn')
        % xlabel('Time (msec)')
        ylabel('R syn')
        subplot(513)
        set(gca,'fontsize',sizz);
        plot(sec*dt*(snc_firings(:,1)),snc_firings(:,2),'.','MarkerSize',12);
        xlim([0 sec*dt*Ttime]);
        ylim([0 Ttime]);
        title('SNc firing')
        % xlabel('Time (msec)')
        ylabel('# of neurons')
        subplot(514)
        set(gca,'fontsize',sizz);
        plot(sec*dt*(1:numel(Rvalsnc)),abs(Rvalsnc))
        xlim([0 sec*dt*Ttime]);
        ylim([0 1.2]);
        title('SNc Rsyn')
        % xlabel('Time (msec)')
        ylabel('R syn')
        subplot(515)
        set(gca,'fontsize',sizz);
        plot(sec*dt*(1:Ttime),kid)
        xlim([0 sec*dt*Ttime]);
        ylim([0 Psnc]);
        title('SNc cell death')
        xlabel('Time (sec)')
        ylabel('# of cells')
        f01=strcat('Org1_',filename0);
        saveas(fig1,f01,'png');
        clear fig1
    end
end

if mfr==1
    plotflag2=1;
    if(plotflag2==1)
        kid=((Psnc)-nc);
        sec=0.001;
        sizz=10;
        fig2=figure(3);
        set(fig2, 'Position', [5, 50, 1920, 955]);
        % set(fig1, 'Position', [50, 100, 2000, 920]);
        subplot(511)
        set(gca,'fontsize',sizz);
        xx=linspace(0,sec*dt*Ttime,dt*Ttime);
        yy=linspace(1,Pstn,Pstn);
        % plot(sec*dt*(stn_firings2(:,1)),stn_firings2(:,2),'.','MarkerSize',12);
        % xlim([0 sec*dt*Ttime]);
        % ylim([0 Pstn]);
        imagesc(xx,yy,mfrstn_subsampled)
        colorbar
        % xlim([0 Ttime]);
        title('STN firing')
        % xlabel('Time (msec)')
        ylabel('# of neurons')
        subplot(512)
        set(gca,'fontsize',sizz);
        plot(sec*(1:numel(Rvalstn)),abs(Rvalstn))
        xlim([0 sec*dt*Ttime]);
        ylim([0 1.2]);
        title('STN Rsyn')
        % xlabel('Time (msec)')
        ylabel('R syn')
        subplot(513)
        set(gca,'fontsize',sizz);
        % plot(sec*dt*(snc_firings(:,1)),snc_firings(:,2),'.','MarkerSize',12);
        % xlim([0 sec*dt*Ttime]);
        % ylim([0 Ttime]);
        xxx=linspace(0,sec*dt*Ttime,dt*Ttime);
        yyy=linspace(1,Psnc,Psnc);
        imagesc(xxx,yyy,mfrsnc_subsampled);
        colorbar
        % xlim([0 sec*dt*Ttime]);
        title('SNc firing')
        % xlabel('Time (msec)')
        ylabel('# of neurons')
        subplot(514)
        set(gca,'fontsize',sizz);
        plot(sec*(1:numel(Rvalsnc)),abs(Rvalsnc))
        xlim([0 sec*dt*Ttime]);
        ylim([0 1.2]);
        title('SNc Rsyn')
        % xlabel('Time (msec)')
        ylabel('R syn')
        subplot(515)
        set(gca,'fontsize',sizz);
        plot(sec*dt*(1:Ttime),kid)
        xlim([0 sec*dt*Ttime]);
        ylim([0 Psnc+5]);
        title('SNc cell death')
        xlabel('Time (sec)')
        ylabel('# of cells')
        f02=strcat('mfr_',filename0);
        saveas(fig2,f02,'png');
        clear fig2
    end
end
