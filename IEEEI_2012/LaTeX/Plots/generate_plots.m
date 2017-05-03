%% simulation 1 EDC plots
clear all; close all;
sim_par = 5;
sim_error = -33;
if (sim_error == -33)
    err_ind = 1;
else
    err_ind = 2;
end
fs = 16000;
Lw = 0.05*fs;
Lw_all = [0.01 0.02 0.03 0.04 0.05]*fs;
Lw_ind = find(Lw_all==Lw);

% Load the processed results
str_load = 'Processed_RPMINT_EDC_PESQ_c';
load(fullfile('Processed Results',str_load));
str_load = strcat('Processed_RPMINT_auto_',num2str(sim_par));
load(fullfile('Results', 'Automatic delta',str_load));

% Load the real acoustic system and calculate EDC
str_load = strcat('Simulation_parameters_',num2str(sim_par));
load(fullfile('Acoustic Systems',str_load));
edc_h = 10*log10(edc(h(:,1)));

c_edcauto = edc_autorpmint(:,err_ind,Lw_ind);
c_edcopt = edc_optrpmint(:,sim_par,err_ind,Lw_ind);

% Plotting EDC of P-MINT opt, P-MINT auto
t = linspace(0,size(edc_optrpmint,1)/fs*1000,size(edc_optrpmint,1));

plot(t,c_edcopt,t(1:120:end),c_edcopt(1:120:end),'xb','LineWidth',1.5, 'MarkerSize',8)
hold all
plot(t,c_edcauto, 'k','LineWidth',1.5);
plot(t(1:length(edc_h)),edc_h,'g',[t(1:120:length(edc_h)) t(length(edc_h) - 20) t(length(edc_h) - 13) t(length(edc_h) - 3) t(length(edc_h))],[edc_h(1:120:end); edc_h(end-20); edc_h(end-13); edc_h(end-3); edc_h(end)],'dg','Linewidth',1.5,'MarkerSize',8);
xlabel('Time [ms]')
ylabel('EDC [dB]')
axis([0 150 -60 0])
grid on
[legend_h,object_h,plot_h,text_strings] = legend('Regularized P-MINT with \delta_{\rm opt}','Regularized P-MINT with \delta_{\rm auto}','{\bf{h}}_1','Location','SouthWest')
str_save = strcat('EDC_optauto_',num2str(sim_par),'_Cm_',num2str(sim_error),'_Ld_',num2str(Lw))
saveas(gca,str_save,'png')
saveas(gca,str_save,'fig')

%% simulation 1: PESQ plots

clear all; close all
sim_par = 5;
sim_error = -33;
if (sim_error == -33)
    err_ind = 1;
else
    err_ind = 2;
end
fs = 16000;
Lw_all = [0.01 0.02 0.03 0.04 0.05]*fs;

% Generate PESQ score for each reverberant signal
clean = wavread('clean_speech.wav');
c_sys = sim_par;
% Load the exact acoustic system
str_load = strcat('Simulation_parameters_',num2str(c_sys));
load(fullfile('Acoustic Systems',str_load));
% Calculate the reverberant signal
rev = fftfilt(h(:,1),clean);
for k = 1:length(Lw_all)
    ref_speech(:,k) = fftfilt(h(1:Lw_all(k),1),clean);
    pesq_rev(k) = pesq(ref_speech(:,k),rev);
end


% Load results
str_load = 'Processed_RPMINT_EDC_PESQ_c';
load(fullfile('Processed Results',str_load));
str_load = strcat('Processed_RPMINT_auto_',num2str(sim_par));
load(fullfile('Results', 'Automatic delta',str_load));


for i = 1:length(Lw_all)
    c_pesqrpmint(i) = pesq_optrpmint(sim_par,err_ind,i);
end
c_pesqautorpmint = pesq_autorpmint(err_ind,:);


PESQ1 = [c_pesqrpmint' c_pesqautorpmint' pesq_rev'];

subplot(1,2,1)
h_bar = bar(Lw_all*1000/fs,PESQ1);
set(h_bar(1),'facecolor',[0 0 0])
set(h_bar(2),'facecolor',[1 1 1])
set(h_bar(3),'facecolor',[0.75 0.75 0.75])
grid on;
ylabel('PESQ Score');
xlabel('Desired Window Length [ms]')
axis([5 55 1 5])

sim_error = [-33:1:-15];
ind = [1 4 9 14 19];
c_error = sim_error(ind);
fs = 16000;
Lw = 0.05*fs;
str_load = strcat('Processed_RPMINT_EDC_PESQ_C_severalCm_simpar_',num2str(sim_par));
load(fullfile('Processed Results',str_load));
str_load = strcat('Processed_RPMINT_auto_Ld_',num2str(Lw),'_simpar_',num2str(sim_par));
load(fullfile('Results', 'Automatic delta',str_load));

c_pesqopt = pesq_optrpmint(ind);
c_pesqauto = pesq_autorpmint(ind);


pesq_rev = ones(5,1)*pesq_rev(end);

PESQ2 = [c_pesqopt' c_pesqauto' pesq_rev];


subplot(1,2,2)
h_bar = bar(Lw_all/fs*1000,PESQ2);
set(h_bar(1),'facecolor',[0 0 0])
set(h_bar(2),'facecolor',[1 1 1])
set(h_bar(3),'facecolor',[0.75 0.75 0.75])
grid on;
legend('Regularized P-MINT with \delta_{\rm opt}','Regularized P-MINT with \delta_{\rm auto}','Reverberant','Location',[0.25, 0.25, .25, .25])
%legend('\delta_{\rm opt}','\delta_{\rm auto}','{\bf{h}}_1','Location',[0.25, 0.25, .25, .25])
ylabel('PESQ Score');
xlabel('Normalized Channel Mismatch E_m [dB]')
%title('L_d = 0.05 f_s (50 ms)')
set(gca,'XTickLabel',{'-33', '-30', '-25', '-20','-15'})
axis([5 55 1 5])
legend('Regularized P-MINT with \delta_{\rm opt}','Regularized P-MINT with \delta_{\rm auto}','Reverberant','Location',[0.25, 0.25, .25, .25])

str_save = strcat('PESQ_optauto_',num2str(sim_par))
saveas(gca,str_save,'png')
saveas(gca,str_save,'fig')

%% simulation 2 EDC plot

clear all; close all
sim_par = 5;
sim_error = -33;
if (sim_error == -33)
    err_ind = 1;
else
    err_ind = 2;
end
fs = 16000;
Lw = 0.05*fs;
Lw_all = [0.01 0.02 0.03 0.04 0.05]*fs;
Lw_ind = find(Lw_all==Lw);
c_Lw = num2str(Lw/fs*1000);

% Load the processed results
str_load = 'Processed_MINT_EDC_PESQ_C';
load(fullfile('Processed Results',str_load));
str_load = 'Processed_RMCLS_EDC_PESQ_C';
load(fullfile('Processed Results',str_load));
str_load = strcat('Processed_CS_minnorm_',num2str(sim_par));
load(fullfile('Processed Results',str_load));
str_load = strcat('Processed_RPMINT_auto_',num2str(sim_par));
load(fullfile('Results', 'Automatic delta',str_load));
c_edcauto = edc_autorpmint(:,err_ind,Lw_ind);

% % Generate the time axis
t = linspace(0,size(c_edcauto,1)/fs*1000,size(c_edcauto,1));
% 
% Load the exact acoustic system
str_load = strcat('Simulation_parameters_',num2str(sim_par));
load(fullfile('Acoustic Systems',str_load));
% Calculate its edc
edc_h = 10*log10(edc(h(:,1)));
c_edcmint = edc_mint(:,sim_par,err_ind);
c_edcrmcls = edc_rmcls(:,sim_par,err_ind,Lw_ind);
c_edccs = edc_cs(:,err_ind,Lw_ind);



plot(t,c_edcmint,t(1:120:end),c_edcmint(1:120:end),'ob','LineWidth',1.5)
hold all
plot(t,c_edccs,'r--','Linewidth',2);
plot(t,c_edcrmcls,'c','Linewidth',3);
plot(t,c_edcauto, 'k','LineWidth',1.5);
plot(t(1:length(edc_h)),edc_h,'g',[t(1:120:length(edc_h)) t(length(edc_h) - 20) t(length(edc_h) - 13) t(length(edc_h) - 3) t(length(edc_h))],[edc_h(1:120:end); edc_h(end-20); edc_h(end-13); edc_h(end-3); edc_h(end)],'dg','Linewidth',1.5,'MarkerSize',8);
grid on
xlabel('Time [ms]')
ylabel('EDC [dB]')
[legend_h,object_h,plot_h,text_strings] = legend('MINT', 'CS', 'RMCLS','Regularized P-MINT with \delta_{\rm auto}','{\bf{h}}_1','Location','SouthWest')
axis([0 150 -60 0])
legendhelp
str_save = strcat('EDC_all_sys_',num2str(sim_par),'_Cm_',num2str(sim_error),'_Ld_',num2str(Lw))
saveas(gca,str_save,'png')

%% simulation 2 PESQ plot

clear all; close all
sim_par = [5];
sim_error = [-33 -15];
c_err = -33;
ind = find(sim_error == c_err);
fs = 16000;
Lw = [0.01 0.02 0.03 0.04 0.05]*fs;

% Load the processed results
str_load = 'Processed_MINT_EDC_PESQ_C';
load(fullfile('Processed Results',str_load));
str_load = 'Processed_RMCLS_EDC_PESQ_C';
load(fullfile('Processed Results',str_load));
str_load = strcat('Processed_CS_minnorm_',num2str(sim_par));
load(fullfile('Processed Results',str_load));
str_load = strcat('Processed_RPMINT_auto_',num2str(sim_par));
load(fullfile('Results', 'Automatic delta',str_load));

clean = wavread('clean_speech.wav');

% Generate PESQ score for each reverberant signal

c_sys = sim_par;
% Load the exact acoustic system
str_load = strcat('Simulation_parameters_',num2str(c_sys));
load(fullfile('Acoustic Systems',str_load));
% Calculate the reverberant signal
rev = fftfilt(h(:,1),clean);
for k = 1:length(Lw)
    ref_speech(:,k) = fftfilt(h(1:Lw(k),1),clean);
    pesq_rev(k) = pesq(ref_speech(:,k),rev);
end

pesq_all = [squeeze(pesq_mint(c_sys,ind,:)) pesq_cs(ind,:)' squeeze(pesq_rmcls(c_sys,ind,:)) pesq_autorpmint(ind,:)' pesq_rev'];
subplot(1,2,1)
h_bar = bar(Lw*1000/fs,pesq_all);
set(h_bar(1),'facecolor',[0.5 0.5 0.5])
set(h_bar(2),'facecolor',[1 1 1])
set(h_bar(3),'facecolor',[0.25 0.25 0.25])
set(h_bar(4),'facecolor',[0 0 0])
set(h_bar(5),'facecolor',[0.75 0.75 0.75])
grid on;
ylabel('PESQ Score');
xlabel('Desired Window Length [ms]')
axis([5 55 1 5])

sim_error = [-33:1:-15];
ind = [1 4 9 14 19];
c_error = sim_error(ind);
fs = 16000;
Lw_s = 0.05*fs;

str_load = strcat('Processed_RPMINT_auto_Ld_',num2str(Lw_s),'_simpar_',num2str(sim_par));
load(fullfile('Results', 'Automatic delta',str_load));

c_pesqauto = pesq_autorpmint(ind);
pesq_rev = ones(5,1)*pesq_rev(end);

% Load the processed results
str_load = 'Processed_MINTieeei_PESQ';
load(fullfile('Processed Results',str_load));
str_load = 'Processed_RMCLSieeei_PESQ';
load(fullfile('Processed Results',str_load));
str_load = strcat('Processed_CSieeei_PESQ');
load(fullfile('Processed Results',str_load));

pesq_all = [pesq_mint' pesq_cs' pesq_rmcls' c_pesqauto' pesq_rev];

subplot(1,2,2)
h_bar = bar(Lw*1000/fs,pesq_all);
set(h_bar(1),'facecolor',[0.5 0.5 0.5])
set(h_bar(2),'facecolor',[1 1 1])
set(h_bar(3),'facecolor',[0.25 0.25 0.25])
set(h_bar(4),'facecolor',[0 0 0])
set(h_bar(5),'facecolor',[0.75 0.75 0.75])
grid on;
ylabel('PESQ Score');
xlabel('Normalized Channel Mismatch E_m [dB]')
%title('L_d = 0.05 f_s (50 ms)')
set(gca,'XTickLabel',{'-33', '-30', '-25', '-20','-15'})
axis([5 55 1 5])
legend('MINT','CS','RMCLS','Regularized P-MINT with \delta_{\rm auto}','Reverberant Signal','Location','NorthWest')
