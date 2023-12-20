%% load spike 
load 'Cond';
load 'SpikeActBZ1';
spike_1_1 = mu0;
data.units{1}= spike_1_1 ;

load 'SpikeActBZ2';
spike_2_1 = mu0;
data.units{2}= spike_2_1 ;

load 'SpikeActBZ3';
spike_3_1 = mu0;
data.units{3}= spike_3_1 ;

load 'SpikeActBZ4';
spike_4_1 = mu0;
spike_4_2 = mu_all;
data.units{4}= spike_4_1 ;
data.units{5}= spike_4_2;


load 'SpikeActBZ5'
spike_5_1 = mu0;
spike_5_2 = mu_all;
data.units{6}= spike_5_1 ;
data.units{7}= spike_5_2 ;


load 'SpikeActBZ6';
spike_6_1 = mu0;
spike_6_2 = mu_all;
data.units{8}= spike_6_1 ;
data.units{9}= spike_6_2;


load 'SpikeActBZ7';
spike_7_1 = mu0;
spike_7_2 = mu_all;
data.units{10}= spike_7_1 ;
data.units{11}= spike_7_2 ;

load 'SpikeActBZ8';
spike_8_1 = mu0;
spike_8_2 = mu_all;
data.units{12}= spike_8_1 ;
data.units{13}= spike_8_2 ;

load 'SpikeActBZ9';
spike_9_1 = mu0;
data.units{14}= spike_9_1 ;

load 'SpikeActBZ10';
spike_10_1 = mu0;
spike_10_2 = mu_all;
data.units{15}= spike_10_1 ;
data.units{16}= spike_10_2;

load 'SpikeActBZ11';
spike_11_1 = mu0;
data.units{17}= spike_11_1 ;

load 'SpikeActBZ12';
spike_12_1 = mu0;
spike_12_2 = mu_all;
data.units{18}= spike_12_1 ;
data.units{19}= spike_12_2 ;

load 'SpikeActBZ13';
spike_13_1 = mu0;
spike_13_2 = mu_all;
data.units{20}= spike_13_1 ;
data.units{21}= spike_13_2 ;

load 'SpikeActBZ14';
spike_14_1 = mu0;
spike_14_2 = mu_all;
data.units{22}= spike_14_1 ;
data.units{23}= spike_14_2;

load 'SpikeActBZ15';
spike_15_1 = mu0;
spike_15_2 = mu_all;
data.units{24}= spike_15_1 ;
data.units{25}= spike_15_2 ;

load 'SpikeActBZ16';
spike_16_1 = mu0;
spike_16_2 = mu_all;
data.units{26}= spike_16_1 ;
data.units{27}= spike_16_2;
clear 'mu0', clear 'mu_all', clear 'su',

%% raster plot:
grp_th = Cond; grp_th(Cond > 8) = Cond(Cond > 8)-8; grp_r = []; grp_r(Cond < 9) = 1; grp_r(Cond > 8) = 2;

% spike_4_1, spike_7_1 , spike_12_2, spike_16_2 are chosen

% >>>> raster plot:

time_samples = 1:8000;
time_stamps = time_samples -1000; 
 
for i= 1:16
    ind_h = find(Cond ==i );
    trial_num = 1:size(ind_h, 2); 

    var_h1_neuron1 = [];
    var_h1_neuron1 = spike_4_1(ind_h, :);
    var_h1_neuron2 = [];
    var_h1_neuron2 = spike_7_1(ind_h, :);

    var_h1_neuron3 = [];
    var_h1_neuron3 = spike_12_2(ind_h, :);

    var_h1_neuron4 = [];
    var_h1_neuron4 = spike_16_2(ind_h, :);
    
    i_ = num2str(i);
    
    figure,
    subplot(221);
    imagesc(time_stamps, trial_num,var_h1_neuron1)
    title_ = (['condition',i_,' ','for neuron1']);
    colormap(flipud(colormap('gray'))); % erase flipud : flip upside down: white to black . black to white
    title(title_);
    xlabel('time');
    ylabel('trial number')

    subplot(222);
    imagesc(time_stamps, trial_num,var_h1_neuron2)
    colormap(flipud(colormap('gray')));
    title_ = (['condition',i_,' ','for neuron2']);
    title (title_);
    xlabel('time');
    ylabel('trial number');

    subplot(223);
    imagesc(time_stamps, trial_num,var_h1_neuron3)
    colormap(flipud(colormap('gray')));
    title_ = (['condition',i_,' ','for neuron3']);
    title(title_);
    xlabel('time');
    ylabel('trial number');

    subplot(224);
    imagesc(time_stamps, trial_num,var_h1_neuron4)
    colormap(flipud(colormap('gray')));
    title_ = (['condition',i_,' ','for neuron4']);
    title(title_);
    xlabel('time');
    ylabel('trial number');  
end


%% PSTH

win = 100; 
psth  = @(x) ndass_smooth(1000*mean(x,1), win); % @ define a function with x as input. mean on first dimention of x(trials) and a window size of win: 100ms so we use 1000 to transfer it from ms to s 1000*1/1000= 1 always remain 1. if you want to change smoothess do it with win. try win = 1. win is sd of gaussian kernel   
t_h = time_stamps;


for i=1:16
    figure;

    
    ind_h = find(Cond == i);
    trial_num = 1 : size(ind_h,2);

    subplot(2, 2, 1);
    hold on
    var_h_neuron1 = [];
    var_h_neuron1 = spike_4_1(ind_h, :);
    plot(t_h, psth(var_h_neuron1), 'r')
    % draw a line for stim onset
    line([0 0], ylim, 'color', 'r')                                                 % Insert Onset line
    xlabel('Time from sample onset (sec.)');
    ylabel('Firing rate (Hz)');
    title_ = (['condition',i_,' ','for neuron1']);
    title (title_);
    
    subplot(2, 2, 2);
    hold on
    var_h_neuron2 = [];
    var_h_neuron2 = spike_7_1(ind_h, :);
    plot(t_h, psth(var_h_neuron2), 'r')
    % draw a line for stim onset
    line([0 0], ylim, 'color', 'r')                                                 % Insert Onset line
    xlabel('Time from sample onset (sec.)');
    ylabel('Firing rate (Hz)');
    title_ = (['condition',i_,' ','for neuron2']);
    title (title_);
    
    subplot(2, 2, 3);
    hold on
    var_h_neuron3 = [];
    var_h_neuron3 = spike_12_2(ind_h, :);
    plot(t_h, psth(var_h_neuron3), 'r')
    % draw a line for stim onset
    line([0 0], ylim, 'color', 'r')                                                 % Insert Onset line
    xlabel('Time from sample onset (sec.)');
    ylabel('Firing rate (Hz)');
    title_ = (['condition',i_,' ','for neuron3']);
    title (title_);
    
    subplot(2, 2, 4);
    hold on
    var_h_neuron4 = [];
    var_h_neuron4 = spike_16_2(ind_h, :);
    plot(t_h, psth(var_h_neuron4), 'r')
    % draw a line for stim onset
    line([0 0], ylim, 'color', 'r')                                                 % Insert Onset line
    xlabel('Time from sample onset (sec.)');
    ylabel('Firing rate (Hz)');
    title_ = (['condition',i_,' ','for neuron4']);
    title (title_);

end


%% mutual information and ROC: fundemental process 

win_  = 60;
psth  = @(x) ndass_smooth(1000*mean(x, 1), win_);
T_st  =1;
T_end = 8008;


units = data.units;
psth_resp = [];

for tri = 1 : size(Cond, 1) % for each trial 
    for sui = 1 : 27
        psth_resp(sui, tri, :) = psth(units{sui}(tri, :)); % not units(1)
    end
end

% >> normalization

ind_b = [1 1000];
for  sui = 1 : size(units, 2)
    var_h=[squeeze(psth_resp(sui, :, :))]; % sui: one neuron.  neuron dimension trial dimension time dimension >> squeeze remove the first D  
    var_b=mean2(var_h(:, ind_b(1):ind_b(2)));
    var_max=max(max(var_h));
    fef_psth_resp_norm(sui, :, :)= (squeeze(psth_resp(sui, :, :)) - var_b)/(var_max - var_b);
end

%% part one 1000:1500 interval 


win_  = 100;
step_ = 5;
win_h = [1000: step_: 1400 ; 1100: step_ :1500]';

IN = [1 2 3 4 5 6 7 8];
OUT=[9 10 11 12 13 14 15 16];

condition_num = input('enter the targeted angle you want to calculate MI for')
comparition_num = input('enter another number 1:8 to compare with the targeted IN data')

ind_h_in = find(Cond == IN(condition_num));
ind_h_in_compare = find(Cond == IN(comparition_num));
ind_h_out = find(Cond ==OUT(condition_num));

mi_in_out = [];
mi_in_in = [];

for sui = 1 : size(units,2) % first for units 
    for ti = 1 : size(win_h,1) % second for time
        t1 = win_h(ti,1);
        t2 = win_h(ti,2);
        
        pref   = nanmean(fef_psth_resp_norm(sui, ind_h_in, t1:t2), 3)';
        npref  = nanmean(fef_psth_resp_norm(sui, ind_h_out, t1:t2), 3)';
       mi_in_out(sui,ti) = ndass_mi([pref; npref],...
        [ones(length(pref), 1); 2*ones(length(npref), 1)], 20, 5);    
         
        npref  = nanmean(fef_psth_resp_norm(sui, ind_h_in_compare, t1:t2), 3)';
         mi_in_in(sui,ti) = ndass_mi([pref; npref],...
        [ones(length(pref), 1); 2*ones(length(npref), 1)], 20, 5);    
     
    end
end

t_h = nanmean(win_h,2)';


figure %for plot : will keep the second(within times) dimension and will mean on the first dimension then calculate standard error of mean and draw it(first dimension: neurons)
hold on
ndass_niceplot(mi_in_out, t_h, 1, 1, 0, 0) 
ndass_niceplot(mi_in_in, t_h, 1, 0, 0, 1)
xlabel('Time from sample onset (sec.)');
ylabel('Mutual information  (a.u.)');
set(gca, 'fontsize', 14, 'fontweight', 'bold');
xlim([999 1501])

legend('MI in&out ','mean MI in&out' , 'MI in&in data' ,'mean MI in&in')
% text(800,0.0,'Chance level','fontsize',18,'fontweight','bold','Color','k');

figure
subplot(121);
hist(mi_in_out);
subplot(122);
hist(mi_in_in);
%% p-value : p-value function returnes p-value 

%% part two 2700:3200 interval 


win_  = 100;
step_ = 5;
win_h = [2700: step_: 3100 ; 2800: step_ :3200]';


condition_num = input('enter the targeted angle you want to calculate MI for');
comparition_num = input('enter another number 1:8 to compare with the targeted IN data');

ind_h_in = find(Cond == IN(condition_num));
ind_h_in_compare = find(Cond == IN(comparition_num));
ind_h_out = find(Cond ==OUT(condition_num));

mi_in_out = [];
mi_in_in = [];

for sui = 1 : size(units,2) % first for units 
    for ti = 1 : size(win_h,1) % second for time
        t1 = win_h(ti,1);
        t2 = win_h(ti,2);
        
        pref   = nanmean(fef_psth_resp_norm(sui, ind_h_in, t1:t2), 3)';
        npref  = nanmean(fef_psth_resp_norm(sui, ind_h_out, t1:t2), 3)';
       mi_in_out(sui,ti) = ndass_mi([pref; npref],...
        [ones(length(pref), 1); 2*ones(length(npref), 1)], 20, 5);    
         
        npref  = nanmean(fef_psth_resp_norm(sui, ind_h_in_compare, t1:t2), 3)';
         mi_in_in(sui,ti) = ndass_mi([pref; npref],...
        [ones(length(pref), 1); 2*ones(length(npref), 1)], 20, 5);    
     
    end
end

t_h = nanmean(win_h,2)';


figure 
hold on
ndass_niceplot(mi_in_out, t_h, 1, 1, 0, 0) 
ndass_niceplot(mi_in_in, t_h, 1, 0, 0, 1)
xlabel('Time from sample onset (sec.)');
ylabel('Mutual information  (a.u.)');
set(gca, 'fontsize', 14, 'fontweight', 'bold');
xlim([2690 3220])

legend('MI in&out ','mean MI in&out' , 'MI in&in data' ,'mean MI in&in')
% text(800,0.0,'Chance level','fontsize',18,'fontweight','bold','Color','k');

figure
subplot(121);
hist(mi_in_out);
subplot(122);
hist(mi_in_in);

%% ROC:
win_  = 60; 
T_st  = 1;
T_end = 8000;

% >>>> roc
condition_num = input('enter the targeted angle you want to calculate roc for');
comparition_num = input('enter another number 1:8 to compare with the targeted IN data');

ind_h_in = find(Cond == IN(condition_num));
ind_h_in_compare = find(Cond == IN(comparition_num));
ind_h_out = find(Cond ==OUT(condition_num));

win_  = 10;
step_ = 5;
win_h = [1000: step_: 1400 ; 1100: step_ :1500]';

roc_in_out = [];
roc_in_in = [];

for sui = 1 : size(units, 2)
    for ti = 1 : size(win_h, 1)
        t1 = win_h(ti, 1);
        t2 = win_h(ti, 2);
        
        pref  = nanmean(fef_psth_resp_norm(sui, ind_h_in, t1:t2), 3)'; % fef_.. our function 
        npref = nanmean(fef_psth_resp_norm(sui, ind_h_out, t1:t2), 3)';
        roc_in_out(sui,ti) = ndass_roc(pref, npref);
                
        npref = nanmean(fef_psth_resp_norm(sui, ind_h_in_compare, t1:t2), 3)';
        roc_in_in(sui,ti) = ndass_roc(pref, npref);
         
    end
end
t_h = nanmean(win_h, 2)';


figure
ndass_niceplot(roc_in_out, t_h, 1, 1, 0, 0) 
ndass_niceplot(roc_in_in, t_h, 1, 0, 0, 1)
xlabel('Time from sample onset (sec.)');
ylabel('ROC  (a.u.)');
set(gca, 'fontsize', 14, 'fontweight', 'bold');
xlim([990 1510 ])

legend('roc in&out ','mean roc in&out' , 'roc in&in data' ,'roc MI in&in')
% text(800,0.0,'Chance level','fontsize',18,'fontweight','bold','Color','k');

figure
subplot(121);
hist(roc_in_out);
subplot(122);
hist(roc_in_in);


%% ROC:
win_  = 60; 
T_st  = 1;
T_end = 8000;

% >>>> roc
condition_num = input('enter the targeted angle you want to calculate roc for');
comparition_num = input('enter another number 1:8 to compare with the targeted IN data');

ind_h_in = find(Cond == IN(condition_num));
ind_h_in_compare = find(Cond == IN(comparition_num));
ind_h_out = find(Cond ==OUT(condition_num));

win_  = 10;
step_ = 5;
win_h = [2700: step_: 3100 ; 2800: step_ :3200]';

roc_in_out = [];
roc_in_in = [];

for sui = 1 : size(units, 2)
    for ti = 1 : size(win_h, 1)
        t1 = win_h(ti, 1);
        t2 = win_h(ti, 2);
        
        pref  = nanmean(fef_psth_resp_norm(sui, ind_h_in, t1:t2), 3)'; % fef_.. our function 
        npref = nanmean(fef_psth_resp_norm(sui, ind_h_out, t1:t2), 3)';
        roc_in_out(sui,ti) = ndass_roc(pref, npref);
                
        npref = nanmean(fef_psth_resp_norm(sui, ind_h_in_compare, t1:t2), 3)';
        roc_in_in(sui,ti) = ndass_roc(pref, npref);
         
    end
end
t_h = nanmean(win_h, 2)';


figure 
hold on
ndass_niceplot(roc_in_out, t_h, 1, 1, 0, 0) 
ndass_niceplot(roc_in_in, t_h, 1, 0, 0, 1)
xlabel('Time from sample onset (sec.)');
ylabel('ROC  (a.u.)');
set(gca, 'fontsize', 14, 'fontweight', 'bold');
xlim([2690 3300 ])

legend('roc in&out ','mean roc in&out' , 'roc in&in data' ,'roc MI in&in')
% text(800,0.0,'Chance level','fontsize',18,'fontweight','bold','Color','k');

figure
subplot(121);
hist(roc_in_out);
subplot(122);
hist(roc_in_in);

%%  SVM
win_  = 100;
step_ = 5;
win_h = [1: step_: 7900; win_+1: step_: 7900+win_]';

preformance = [];

same_angle_num = input('enter the targeted angle you want to recognize');
different_angle_num = input('enter another number angle  1:8 ');

ind_in_angle = ismember(Cond,IN(same_angle_num ));
ind_out_another_angle = ismember(Cond,OUT(different_angle_num));


for ti = 1 : size(win_h, 1)
    t1 = win_h(ti, 1);
    t2 = win_h(ti, 2);

    group1 = nanmean(fef_psth_resp_norm(:,ind_in_angle , t1:t2), 3)';
    group2 = nanmean(fef_psth_resp_norm(:, ind_out_another_angle, t1:t2), 3)';

    svm_res = ndass_svm([group1; group2], [ones(length(group1),1); 2*ones(length(group2),1)],0.7, 10) ;
    preformance(:,ti) = ((svm_res.pt)-1/2)/(1-1/2);

end

t_h = nanmean(win_h, 2)';
t_h = t_h - 1000;

figure
hold on
ndass_niceplot(preformance, t_h, 1, 1, 0, 0)
line([0 0], ylim, 'color', 'r')                                                 % Insert Onset line
xlabel('Time from sample onset (sec.)');
ylabel('Performance (a.u)');
title('SVM')
line([0 0], ylim, 'Color', 'k');
line(xlim, [0.0 0.0], 'Color', 'k');
set(gca, 'fontsize', 14, 'fontweight', 'bold');
xlim([-1000 7000]);
ylim([-0.3 1]);


%% load lfp 

load ('Cond');

load 'LFP1.mat';
LFP1 = LFP;
load ('LFP2.mat');
LFP2 = LFP;
load ('LFP3.mat');
LFP3 = LFP;
load ('LFP4.mat');
LFP4 = LFP;
load ('LFP5.mat');
LFP5 = LFP;
load ('LFP6.mat');
LFP6 = LFP;
load ('LFP7.mat');
LFP7 = LFP;
load ('LFP8.mat');
LFP8 = LFP;
load ('LFP9.mat');
LFP9 = LFP;
load ('LFP10.mat');
LFP10 = LFP;
load ('LFP11.mat');
LFP11 = LFP;
load ('LFP12.mat');
LFP12 = LFP;
load ('LFP13.mat');
LFP13 = LFP;
load ('LFP14.mat');
LFP14 = LFP;
load ('LFP15.mat');
LFP15 = LFP;
load ('LFP16.mat');
LFP16 = LFP;
clear LFP

%% preprocessing

T_st=1;T_end=8000;

% >>>> notch filter 
LFP1 = ndass_nothfilter(LFP1,50);
LFP2 = ndass_nothfilter(LFP2,50);
LFP3 = ndass_nothfilter(LFP3,50);
LFP4 = ndass_nothfilter(LFP4,50);
LFP5 = ndass_nothfilter(LFP5,50);
LFP6 = ndass_nothfilter(LFP6,50);
LFP7 = ndass_nothfilter(LFP7,50);
LFP8 = ndass_nothfilter(LFP8,50);
LFP9 = ndass_nothfilter(LFP9,50);
LFP10 = ndass_nothfilter(LFP10,50);
LFP11 = ndass_nothfilter(LFP11,50);
LFP12 = ndass_nothfilter(LFP12,50);
LFP13 = ndass_nothfilter(LFP13,50);
LFP14 = ndass_nothfilter(LFP14,50);
LFP15 = ndass_nothfilter(LFP15,50);
LFP16 = ndass_nothfilter(LFP16,50);

% isequal(LFP,LFP10) >>> returns false >>> notch filter is done 

% >>>> artifact remove 
% sigma = 2 
% commented the plot. if observation is required uncomment 

artif_lfp1 = ndass_rmartifact(LFP1,1);
artif_lfp2 = ndass_rmartifact(LFP2,1);
artif_lfp3 = ndass_rmartifact(LFP3,1);
artif_lfp4 = ndass_rmartifact(LFP4,1);
artif_lfp5 = ndass_rmartifact(LFP5,1);
artif_lfp6 = ndass_rmartifact(LFP6,1);
artif_lfp7 = ndass_rmartifact(LFP7,1);
artif_lfp8 = ndass_rmartifact(LFP8,1);
artif_lfp9 = ndass_rmartifact(LFP9,1);
artif_lfp10 = ndass_rmartifact(LFP10,1);
artif_lfp11 = ndass_rmartifact(LFP11,1);
artif_lfp12 = ndass_rmartifact(LFP12,1);
artif_lfp13 = ndass_rmartifact(LFP13,1);
artif_lfp14 = ndass_rmartifact(LFP14,1);
artif_lfp15 = ndass_rmartifact(LFP15,1);
artif_lfp16 = ndass_rmartifact(LFP16,1);


% >>>> normalization: z score

normalizer = @(x) (x - nanmean(x(:)))/(nanstd(x(:)));
LFP1 = normalizer(LFP1);
LFP2 = normalizer(LFP2);
LFP3 = normalizer(LFP3);
LFP4 = normalizer(LFP4);
LFP5 = normalizer(LFP5);
LFP6 = normalizer(LFP6);
LFP7 = normalizer(LFP7);
LFP8 = normalizer(LFP8);
LFP9 = normalizer(LFP9);
LFP10 = normalizer(LFP10);
LFP11 = normalizer(LFP11);
LFP12 = normalizer(LFP12);
LFP13 = normalizer(LFP13 );
LFP14 = normalizer(LFP14);
LFP15 = normalizer(LFP15);
LFP16 = normalizer(LFP16);

% we chose electrode one lfp  but it can be done for all lfps 
ind_h1 = ismember(Cond,1)&(~artif_lfp1);
ind_h2 = ismember(Cond,2)&(~artif_lfp1);
ind_h3 = ismember(Cond,3)&(~artif_lfp1);
ind_h4 = ismember(Cond,4)&(~artif_lfp1);
ind_h5 = ismember(Cond,5)&(~artif_lfp1);
ind_h6 = ismember(Cond,6)&(~artif_lfp1);
ind_h7 = ismember(Cond,7)&(~artif_lfp1);
ind_h8 = ismember(Cond,8)&(~artif_lfp1);
ind_h9 = ismember(Cond,9)&(~artif_lfp1);
ind_h10 = ismember(Cond,10)&(~artif_lfp1);
ind_h11 = ismember(Cond,11)&(~artif_lfp1);
ind_h12 = ismember(Cond,12)&(~artif_lfp1);
ind_h13 = ismember(Cond,13)&(~artif_lfp1);
ind_h14 = ismember(Cond,14)&(~artif_lfp1);
ind_h15 = ismember(Cond,15)&(~artif_lfp1);
ind_h16 = ismember(Cond,16)&(~artif_lfp1);

% check if the 26 artifacts are subtracted 
% i= sum(ind_h1) +  sum(ind_h2)+ sum(ind_h3)+ sum(ind_h4)+ sum(ind_h5)+ sum(ind_h6)+ sum(ind_h7)+ sum(ind_h8)+ sum(ind_h9)+ sum(ind_h10)+ sum(ind_h11)+ sum(ind_h12)+ sum(ind_h13)+ sum(ind_h14)+ sum(ind_h15)+ sum(ind_h16)

%% PSTH using Multitaper 1000:1500 interval 

% Settings Chronux
params.Fs=1000; % sampling frequency
params.fpass=[1 100]; % band of frequencies to be kept
params.tapers=[2 3];%[3 5]; % taper parameters
params.pad=2; % pad factor for fft
params.err=[0 0.05];
params.trialave=0;
movingwin=[0.2 0.050];%movingwin=[0.5 0.05];

[S,f]=mtspectrumc(LFP1(:,1000:1500)',params);
S = S';
[Sb,f]=mtspectrumc(LFP1(:,1:2:1000)',params);
Sb= Sb';


% relative power
figure;

hold on
plot(f,nanmean(S(ind_h1,:)./ Sb(ind_h1,:)),'r') % condition one
plot(f,nanmean(S(ind_h2,:)./ Sb(ind_h2,:)),'b')
plot(f,nanmean(S(ind_h3,:)./ Sb(ind_h3,:)),'g')
plot(f,nanmean(S(ind_h4,:)./ Sb(ind_h4,:)),'color', [0.9290 0.6940 0.1250])
plot(f,nanmean(S(ind_h5,:)./ Sb(ind_h5,:)),'m')
plot(f,nanmean(S(ind_h6,:)./ Sb(ind_h6,:)),'color',[0.3010 0.7450 0.9330])
plot(f,nanmean(S(ind_h7,:)./ Sb(ind_h7,:)),'k')
plot(f,nanmean(S(ind_h8,:)./ Sb(ind_h8,:)),'y')
title ('Baseline normalized')
xlabel('Frequency');
ylabel('Power (a.u.)');
set(gca, 'fontsize', 14, 'fontweight', 'bold');
legend()

%% PSTH using Multitaper 2700:3200 interval 

% Settings Chronux
params.Fs=1000; % sampling frequency
params.fpass=[1 100]; % band of frequencies to be kept
params.tapers=[2 3];%[3 5]; % taper parameters
params.pad=2; % pad factor for fft
params.err=[0 0.05];
params.trialave=0;
movingwin=[0.2 0.050];%movingwin=[0.5 0.05];

[S,f]=mtspectrumc(LFP1(:,2700:3200)',params);
S = S';
[Sb,f]=mtspectrumc(LFP1(:,1:2:1000)',params);
Sb= Sb';


% relative power
figure;

hold on
plot(f,nanmean(S(ind_h1,:)./ Sb(ind_h1,:)),'r') % condition one
plot(f,nanmean(S(ind_h2,:)./ Sb(ind_h2,:)),'b')
plot(f,nanmean(S(ind_h3,:)./ Sb(ind_h3,:)),'g')
plot(f,nanmean(S(ind_h4,:)./ Sb(ind_h4,:)),'color', [0.9290 0.6940 0.1250])
plot(f,nanmean(S(ind_h5,:)./ Sb(ind_h5,:)),'m')
plot(f,nanmean(S(ind_h6,:)./ Sb(ind_h6,:)),'color',[0.3010 0.7450 0.9330])
plot(f,nanmean(S(ind_h7,:)./ Sb(ind_h7,:)),'k')
plot(f,nanmean(S(ind_h8,:)./ Sb(ind_h8,:)),'y')
title ('Baseline normalized')
xlabel('Frequency');
ylabel('Power (a.u.)');
set(gca, 'fontsize', 14, 'fontweight', 'bold');
legend()

%% Calculation of phase amplitude coupling by Tort et al 2010
% Specific band  phase amplitude coupling 

artif  = (artif_lfp1 | artif_lfp2);
fs = 1000;
t1_delay = 1000; t2_delay = 1500;
fl = [4:1:8];
fh = [30:3:60];
f = [fl , fh];
[analytic_sig1, f] = ndass_wavelet(LFP1, f, fs);
[analytic_sig2, f] = ndass_wavelet(LFP2, f, fs);

N_shuff = 10; % 
ind_t = t1_delay:t2_delay;
sig_fef_ph = angle(analytic_sig1(:,ismember(f,fl),:));
sig_it_ph  = angle(analytic_sig2(:,ismember(f,fl),:));
sig_fef_po = abs(analytic_sig1(:,ismember(f,fh),:));
sig_it_po  = abs(analytic_sig2(:,ismember(f,fh),:));

for ci=1:2
    
    if ci == 1 vv = [1 2 3 4 5 6 7 8];else  vv = [9 10 11 12 13 14 15 16];end
    ind_h = []; ind_h = ismember(Cond ,vv)&(~artif);
    
    for f1= 1:length(fh)        
        for f2= 1:length(fl)
         
            tic1 = tic;
            % open ndass_pac_mi..
            [pac_it_phase_fef_pow_raw(ci,f1,f2) ,pac_it_phase_fef_pow_n(ci,f1,f2),pac_it_phase_fef_pow_shuff_m(ci,f1,f2),~] ...
                = ndass_pac_mi_shuff_corrected(reshape(squeeze(sig_it_ph(ind_h,f2,ind_t))',[],1),...
                reshape(squeeze(sig_fef_po(ind_h,f1,ind_t))',[],1),N_shuff);
            
            [pac_fef_phase_it_pow_raw(ci,f1,f2) ,pac_fef_phase_it_pow_n(ci,f1,f2),pac_fef_phase_it_pow_shuff_m(ci,f1,f2),~] ...
                = ndass_pac_mi_shuff_corrected(reshape(squeeze(sig_fef_ph(ind_h,f2,ind_t))',[],1),...
                reshape(squeeze(sig_it_po(ind_h,f1,ind_t))',[],1),N_shuff);
            
            [pac_fef_phase_fef_pow_raw(ci,f1,f2),pac_fef_phase_fef_pow_n(ci,f1,f2),pac_fef_phase_fef_pow_shuff_m(ci,f1,f2),~] ...
                = ndass_pac_mi_shuff_corrected(reshape(squeeze(sig_fef_ph(ind_h,f2,ind_t))',[],1),...
                reshape(squeeze(sig_fef_po(ind_h,f1,ind_t))',[],1),N_shuff);
            
            [pac_it_phase_it_pow_raw(ci,f1,f2)  ,pac_it_phase_it_pow_n(ci,f1,f2), pac_it_phase_it_pow_shuff_m(ci,f1,f2),~]   ...
                = ndass_pac_mi_shuff_corrected(reshape(squeeze(sig_it_ph(ind_h,f2,ind_t))',[],1),...
                reshape(squeeze(sig_it_po(ind_h,f1,ind_t))',[],1),N_shuff);
            
            toc(tic1)
        end
    end        
end

% >>>> Preview PAC shuffel corrected 

% result : not a significant result no effect 
a=5e-4;b=-5e-4; % 3 * 10^-5 do it with different numbers replace 3 with 4 
for si = 1 : 2
    if si == 1
        figure('name','PAC for In')
    else
        figure('name','PAC for Out')        
    end
    for ci =1:8
        
        switch ci
            case 1
                val_h =  pac_fef_phase_it_pow_n; title_ = 'electrode 2 Phase electrode 1 Power';
            case 2
                val_h =  pac_it_phase_fef_pow_n; title_ = 'electrode 1 Phase electrode 2 Power';
            case 3
                val_h =  pac_fef_phase_fef_pow_n; title_ = 'electrode 2 Phase electrode 2  Power';
            case 4
                val_h =  pac_it_phase_it_pow_n; title_ = 'electrode 1 Phase electrode 1 Power';
        end
        
        ax=subplot(2,4,ci);
        var_h=(squeeze((val_h(si,:,:))));
        h=pcolor(fl,fh,var_h);
        h.EdgeColor = 'none';
        colormap(jet); clb= colorbar;
        ax.YDir='normal';
        ax.YTick=[4 15 35 100];
         caxis([b a])
        title([title_] );
        xlabel('Phase Freq. (Hz)')
        ylabel('Power Freq. (Hz)')        
        
    end
end



%% field_field  coherency and spike_field cohernecy vector summation rate matched

 %4-8 hertz , 30-60
 % >>>  Wavelet
 % electrode 2 and electrode 9 is chosen  
fs = 1000;
ind_t1 = 1000:1500
t  = [1:8000]/1000; % in second
f1_1  = [1:4]; 
f2_1 = [30:60];
i = input('give the lfp index you want to calculate phase-spike coherence');
lfp = ['LFP',num2str(i)];
lfp = eval(lfp);
[analytic_sig1, f1_1] = ndass_wavelet(lfp, f1_1, fs);
[analytic_sig2, f2_1] = ndass_wavelet(lfp, f2_1, fs);
spike = ['spike_',num2str(i),'_1'];
spike =eval(spike);
phase_lfp_tfmap_f1 = angle(analytic_sig1); % phase of lfp 
phase_lfp_tfmap_f2 = angle(analytic_sig2);
rate_spike = (spike(:,ind_t)); 


%n_trials=size(rate_fef,1); % number of trials 

for fi=1:length(f1_1) % f [2:30]
    
    val_fef_it_=[];val_it_it_=[];val_diff_it_=[];
    phase_lfp_f1 =squeeze(phase_lfp_tfmap_f1(:,fi,ind_t)); % har dafe zamano dare triala ro dare  ferekanso behesh mide va bad bod ferekanso hazf mikone
    for ci=1:2
        if ci<2; condi_h=[1];else condi_h=[9]; end % in and out and same object in both 
        ind_h = [];
        ind_h = ismember(Cond,condi_h)&(~artif);
        
     
        phase_lfp_h = []; phase_lfp_h =phase_lfp_f1(ind_h,:); % phase
        rate_spike_h = [];rate_spike_h = rate_spike(ind_h,:); % spike 
        
        N_sam_sp = 30; % number of sample windowing. 30 be bala behtare 
        
        %  ~ : dont save the second output
     
        [spl_lfp_spike_f1(ci,fi),~] = ndass_spl_percondition(phase_lfp_h,rate_spike_h,N_sam_sp);
 
        
    end
end



for fi=1:length(f2_1) % f [2:30]
    
    val_fef_it_=[];val_it_it_=[];val_diff_it_=[];
    phase_lfp_f2 =squeeze(phase_lfp_tfmap_f2(:,fi,ind_t)); % har dafe zamano dare triala ro dare  ferekanso behesh mide va bad bod ferekanso hazf mikone
    for ci=1:2
        if ci<2; condi_h=[1];else condi_h=[9]; end % in and out and same object in both 
        ind_h = [];
        ind_h = ismember(Cond,condi_h)&(~artif);
        
     
        phase_lfp_h = []; phase_lfp_h =phase_lfp_f2(ind_h,:); % phase
        rate_spike_h = [];rate_spike_h = rate_spike(ind_h,:); % spike 
        
        N_sam_sp = 30; 
        
       
        [spl_lfp_spike_f2(ci,fi),~] = ndass_spl_percondition(phase_lfp_h,rate_spike_h,N_sam_sp);
 
        
    end
end

%% Preview SPL
figure( 'Name','Spike Field Coherency')
subplot(121);
hold on
var_h = [];
var_h = spl_lfp_spike_f1(1,:);
plot(f1_1,var_h,'g')

var_h = [];
var_h = spl_lfp_spike_f1(2,:);
plot(f1_1,var_h,'b')
title ('SPL of FEF spike  with FEF lfp ')
xlabel('Frequency (Hz)');
ylabel('Coherency (a.u.)');
set(gca, 'fontsize', 14, 'fontweight', 'bold');
legend('In','Out')

subplot(122);
hold on
var_h = [];
var_h = spl_lfp_spike_f2(1,:);
plot(f2_1,var_h,'g')

var_h = [];
var_h = spl_lfp_spike_f2(2,:);
plot(f2_1,var_h,'b')
title ('SPL of FEF spike  with FEF lfp ')
xlabel('Frequency (Hz)');
ylabel('Coherency (a.u.)');
set(gca, 'fontsize', 14, 'fontweight', 'bold');
legend('In','Out')



%% spike field 2700=3200 interval

% field_field  coherency and spike_field cohernecy vector summation rate matched

 %4-8 hertz , 30-60
 % >>>  Wavelet
 % electrode 2 and electrode 9 is chosen  
fs = 1000;
ind_t1 = 2700:3200
t  = [1:8000]/1000; % in second
f1_1  = [1:4]; % for cognitive tasks this frequency interval mostly include spikes
f2_1 = [30:60];
i = input('give the lfp index you want to calculate phase-spike coherence');
lfp = ['LFP',num2str(i)];
lfp = eval(lfp);
[analytic_sig1, f1_1] = ndass_wavelet(lfp, f1_1, fs);
[analytic_sig2, f2_1] = ndass_wavelet(lfp, f2_1, fs);
spike = ['spike_',num2str(i),'_1'];
spike =eval(spike);
phase_lfp_tfmap_f1 = angle(analytic_sig1); % phase of lfp 
phase_lfp_tfmap_f2 = angle(analytic_sig2);
rate_spike = (spike(:,ind_t)); 


%n_trials=size(rate_fef,1); % number of trials 

for fi=1:length(f1_1) % f [2:30]
    
    val_fef_it_=[];val_it_it_=[];val_diff_it_=[];
    phase_lfp_f1 =squeeze(phase_lfp_tfmap_f1(:,fi,ind_t)); % har dafe zamano dare triala ro dare  ferekanso behesh mide va bad bod ferekanso hazf mikone
    for ci=1:2
        if ci<2; condi_h=[1];else condi_h=[9]; end % in and out and same object in both 
        ind_h = [];
        ind_h = ismember(Cond,condi_h)&(~artif);
        
     
        phase_lfp_h = []; phase_lfp_h =phase_lfp_f1(ind_h,:); % phase
        rate_spike_h = [];rate_spike_h = rate_spike(ind_h,:); % spike 
        
        N_sam_sp = 30; % number of sample windowing. 30 be bala behtare 
        
        % the below function has two outputs the ~ say dont save the second
       % open the function 
        [spl_lfp_spike_f1(ci,fi),~] = ndass_spl_percondition(phase_lfp_h,rate_spike_h,N_sam_sp);
 
        
    end
end



for fi=1:length(f2_1) % f [2:30]
    
    val_fef_it_=[];val_it_it_=[];val_diff_it_=[];
    phase_lfp_f2 =squeeze(phase_lfp_tfmap_f2(:,fi,ind_t)); % har dafe zamano dare triala ro dare  ferekanso behesh mide va bad bod ferekanso hazf mikone
    for ci=1:2
        if ci<2; condi_h=[1];else condi_h=[9]; end % in and out and same object in both 
        ind_h = [];
        ind_h = ismember(Cond,condi_h)&(~artif);
        
     
        phase_lfp_h = []; phase_lfp_h =phase_lfp_f2(ind_h,:); % phase
        rate_spike_h = [];rate_spike_h = rate_spike(ind_h,:); % spike 
        
        N_sam_sp = 30; 
        
    
        [spl_lfp_spike_f2(ci,fi),~] = ndass_spl_percondition(phase_lfp_h,rate_spike_h,N_sam_sp);
 
        
    end
end

%% Preview SPL
figure( 'Name','Spike Field Coherency')
subplot(121);
hold on
var_h = [];
var_h = spl_lfp_spike_f1(1,:);
plot(f1_1,var_h,'g')

var_h = [];
var_h = spl_lfp_spike_f1(2,:);
plot(f1_1,var_h,'b')
title ('SPL of FEF spike  with FEF lfp ')
xlabel('Frequency (Hz)');
ylabel('Coherency (a.u.)');
set(gca, 'fontsize', 14, 'fontweight', 'bold');
legend('In','Out')

subplot(122);
hold on
var_h = [];
var_h = spl_lfp_spike_f2(1,:);
plot(f2_1,var_h,'g')

var_h = [];
var_h = spl_lfp_spike_f2(2,:);
plot(f2_1,var_h,'b')
title ('SPL of FEF spike  with FEF lfp ')
xlabel('Frequency (Hz)');
ylabel('Coherency (a.u.)');
set(gca, 'fontsize', 14, 'fontweight', 'bold');
legend('In','Out')


