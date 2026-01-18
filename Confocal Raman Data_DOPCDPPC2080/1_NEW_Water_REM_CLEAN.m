
% This script is used to remove water from lipid film scan! It will create
% data for next script. 

%%
clc 
clear all
close all 

set(0, 'DefaultAxesFontSize', 14);

%cd 
%loading files 
files = dir('*.mat');

for i=1:length(files)
    load(files(i).name, '-mat');
end

%% Data 

name_scan = Struct_shortscanGspatGspecwrongspatLine026SpecData1 ; 
name_water = Struct_waterSpectrum042SpecData1 ; 

data_scan_raw = name_scan.data;
Water_spec_raw = name_water.data; 


%% Line Scan select

i_scan =  1;                % start scan
f_scan = size(data_scan_raw,1) ; % final scan

data_scan = data_scan_raw(i_scan:f_scan,:);
x_nm = name_scan.axisscale{2,1};
laser = 531.98 ; %nm
x_cm1 = 10^7 .*(1./laser - 1./x_nm); %Raman shift cm-1

%% Remove BG substract and correct 

i_bg = 599; f_bg = 684; %BG

data_all = [Water_spec_raw; data_scan]; 

for cn = 1:1:size(data_all,1);
    BG(1,1) = mean(data_all(cn,i_bg:f_bg));
    %M(1,cn) = BG(1,1);
    data_all_BG(cn,:) = data_all(cn,:) - BG(1,1); % BG substraction
end 

%BG correction 
% Define finalx positions
x1 = f_bg;
x2 = size(data_all,2);
%dx = x2-x1+1; %pluse 1 becuse you have water in ALL data. 
x_vals = (x1:x2)';
y_BASE = zeros(1600,1);

for cz = 1:1:size(data_all,1);
    y1 = mean(data_all_BG(cz,x1-85:x1), 2);
    y2 = mean(data_all_BG(cz,x2-60:x2), 2);
    y2_all(cz,1) = y2;
    k = (y2 - y1) / (x2 - x1); 
    m = y1 - k * x1;
    
    y_BASE(x1:end) = k * x_vals + m; %replacing part of y_base vector.
    y_BASE(1:x1-1) = 0; % making sure that rest is zero.
    
    data_all_BG(cz,:) = data_all_BG(cz,:) - y_BASE';
end


figure 
plot(y2_all,'ko-','linewidth',1);
title('all Y2 values , 1 = water spec');

%% Normalizide to total integral 

i_tot = 150; f_tot = 1550; % Tot normalizaton 

for cn = 1:1:size(data_all,1);
    totI_integral(cn,1) = trapz(data_all_BG(cn,i_tot:f_tot)); %taking out integral
    Data_all_totInorm(cn,:) = data_all_BG(cn,:)./totI_integral(cn,1); %normalizing to CN integral
end 

%%

Water_BG = data_all_BG(1,:)'; %BG substracted

Data_raw = data_all_BG(2:end,:)';
Data_BG = data_all_BG(2:end,:)'; %BG substracted
Data_scan = Data_all_totInorm(2:end,:)'; %Normalized to totI integral

ScanV = 1:1:size(Data_scan,2); 

%% REPALCEMENT of wrong data and smootening 

%Data_scan(:,58) = Data_scan(:,57); %replace contaminated SCAN nr 58.
%Data_scan(:,61) = Data_scan(:,62); %replace contaminated SCAN nr 58.
%Data_scan(:,76) = Data_scan(:,77); %replace contaminated SCAN nr 58.

SCAN_smooth = smoothdata(Data_scan,1,"movmean",15); 

%% CHECK ALL SPEC

spec_init = 10; %spect_fin = end;

fig1 = figure; %1
plot(x_cm1,Data_raw(:,spec_init:end),'-');
grid minor
set(gca, 'ColorOrder', cool(size(Data_scan,2)));
%legend('Data BG+Norm')
title('RAW Scan BG corrected');
ylabel('Normalized CCD cts');
%ylim([-0.001 0.025]);
xlim([x_cm1(1) x_cm1(end)]);
%print(fig1,'RAWScanSpce', '-dsvg');

fig2 = figure; %1
plot(x_cm1,Data_scan(:,spec_init:end),'-');
grid minor
set(gca, 'ColorOrder', cool(size(Data_scan,2)));
%legend('Data BG+Norm')
title('Scan Normalized to Tot Int');
ylabel('Normalized CCD cts');
ylim([-0.001 0.025]);
xlim([x_cm1(1) x_cm1(end)]);
%print(fig2,'BGnormScanSpce', '-dsvg');


fig3 = figure; %1
plot(x_cm1,SCAN_smooth(:,spec_init:end),'-');
grid minor
set(gca, 'ColorOrder', cool(size(Data_scan,2)));
%legend('Data BG+Norm')
title('Scan Normalized to Tot Int');
ylabel('Normalized CCD cts');
ylim([-0.001 0.025]);
xlim([x_cm1(1) x_cm1(end)]);
%print(fig3,'BGnormScanSpceSmooth', '-dsvg');


fig4 = figure; %1
plot(ScanV,totI_integral(2:end),'ok-');
grid minor
%set(gca, 'ColorOrder', cool(size(Data_scan,2)));
%legend('Data BG+Norm')
title('Scan Normalized to Tot Int');
ylabel('Normalized CCD cts');
%ylim([-0.001 0.025]);
%xlim([x_cm1(1) x_cm1(end)]);


%% FIT and Remove WATER
i_water = 1230 ; f_water = 1550; 
fit_vector = i_water:1:f_water ; 

%C_film = [DOPC_fit_film, DPPC_fit_film]; C_film = double(C_film); %format neede for lsqlin fit
W_fit = Water_BG(i_water:f_water,1); W_fit = double(W_fit); 

D_all = [];
D_fit_all =[]; 
D = [];
mean_D=[];
TSS = [];
TSS_all = [];
RSS = [];
RSS_all =[]; 
R2_value = []; 
R2_all_film = []; 

cm = 1;
for cm = 1:1:size(Data_scan,2)
    D = SCAN_smooth(i_water:f_water,cm); %Data set to fit

    %D(exapt_X,:) = [];
    D_all(:,cm) = D(:,1);
    D = double(D);
    lin = lsqlin(W_fit, D,[],[],[],[],[0],[inf]);
    
    w1_lsq_film_raw(cm,1) = lin(1,1); %Coeficents
    
    D_fited = lin(1,1) .* W_fit;
    D_fit_all(:,cm) = D_fited(:,1); %save fitting 
    
    %R2 Values
    mean_D = mean(D); 
    TSS = sum((D - mean_D).^2);
    TSS_all(cm,1) = TSS(1,1);
    RSS = sum((D - D_fited).^2);
    RSS_all(cm,1) = RSS;
    R2_value = 1 -(RSS/TSS);
    R2_all_film(cm,1) = R2_value;   
end 

%% FITTING VALUES 

spec = 76;

figure; %2
plot(fit_vector,D_all(:,spec),'-','linewidth',1);
hold on 
plot(fit_vector,D_fit_all(:,spec),'-','linewidth',1);
hold off 
grid minor
%legend('Data BG+Norm')
legend('Raw','FIT')
ylabel('Normalized CCD cts');
xlabel('INDEX [nr.]');
%ylim([0 6]);

figure; %2
plot(ScanV,w1_lsq_film_raw,'ko-','linewidth',1);
grid minor
%legend('Data BG+Norm')
legend('W1')
ylabel('W1');
xlabel('INDEX [scan nr.]');
%ylim([0 0.1]);

figure; %2
plot(ScanV,R2_all_film,'ob-','linewidth',1);
hold off 
grid minor
legend('R2 Water FIT')
ylabel('R2 coefficint');
xlabel('INDEX [scan nr.]');
ylim([0 1]);

%% COMPARE WITH WATER REF

WaterScanNr = 100;
Coef = w1_lsq_film_raw(WaterScanNr,1);
Norm_WaterRefSpec = Coef.* Water_BG; 

figure; %1
%plot(x_cm1,Data_scan(:,70),'k-','linewidth',1.5);
hold on 
plot(x_cm1,Data_scan(:,WaterScanNr),'o-','linewidth',1);
plot(x_cm1,SCAN_smooth(:,WaterScanNr),'k-','linewidth',2);
plot(x_cm1,Norm_WaterRefSpec,'r-','linewidth',1.5);
hold off 
grid minor
%set(gca, 'ColorOrder', cool(size(Data_scan,2)));
%legend('Data BG+Norm')
legend('Scan Spec','Fitted WATER')
title('Scan smoothing');
ylabel('Normalized CCD cts');
xlabel('INDEX [nr.]');
%ylim([0 6]);

%% GET ALL WATER FITTED and Remove water contribution

WATER_ToRemove = w1_lsq_film_raw' .* Water_BG; 
SCAN_WRemoved = Data_scan - WATER_ToRemove; 

figure; %1
plot(x_cm1,WATER_ToRemove,'-','linewidth',1.5); 
hold on 
grid minor
set(gca, 'ColorOrder', cool(size(Data_scan,2)));
%legend('Data BG+Norm')
legend('Scan')
title('Watter fitted');
ylabel('Normalized CCD cts');
xlabel('INDEX [nr.]');
%ylim([0 6]);

figure; %1
plot(x_cm1,SCAN_WRemoved(:,60:end),'-','linewidth',1.5); 
grid minor
set(gca, 'ColorOrder', cool(size(Data_scan,2)));
%legend('Data BG+Norm')
legend('Scan')
title('Scan removed water');
ylabel('Normalized CCD cts');
xlabel('INDEX [nr.]');
%ylim([0 6]);

%% Compare without Removed water 

Spec_Comp = 100;

figure; %1
plot(x_cm1,Data_scan(:,Spec_Comp),'-','linewidth',1.5);
hold on 
plot(x_cm1,SCAN_WRemoved(:,Spec_Comp),'-','linewidth',1.5);
hold off 
grid minor
%set(gca, 'ColorOrder', cool(size(Data_scan,2)));
%legend('Data BG+Norm')
legend('RAW','WATER Removed')
title('Scan smoothing');
ylabel('Normalized CCD cts');
xlabel('INDEX [nr.]');
%ylim([0 6]);


%% SAVE DATA 

%save('Data_scan.mat','Data_scan'); %<----- open
save('SCAN_WRemoved.mat','SCAN_WRemoved'); %<----- open
save('WATER_ToRemove.mat','WATER_ToRemove'); %<----- open
%save('SCAN_Smooth.mat','SCAN_Smooth'); %<----- open
save('x_cm1.mat','x_cm1'); %<----- open

