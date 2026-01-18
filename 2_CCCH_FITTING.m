
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

%% Files 

DOPC_REFscan = DOPC_Scan_WRemoved'; 
DPPC_REFscan = DPPC_Scan_WRemoved'; 

%name_scan = Struct_shortGspatPspeclineLine036SpecData1 ; 
%scan_raw = name_scan.data;

scan_raw = SCAN_WRemoved'; 

%Scan_Data_NormTotI = Data_scan'; % Extracted data from the capillary normalized to totInt
%Scan_WRemoved = SCAN_WRemoved'; % Data after water removal 
%Scan_WRemovedSmoth = SCAN_Smooth'; % Data after water removal and smothing. 

%% Scan select and repalce 

i_scan =  1;                % start scan
f_scan = size(scan_raw,1) ; % final scan

%{
Scan_Data = Scan_Data_NormTotI(i_scan:f_scan,:);
Scan_WRem = Scan_WRemoved(i_scan:f_scan,:);
Scan_WRemSmoth = Scan_WRemovedSmoth(i_scan:f_scan,:);
%}

%% SELECT REF 

% 0 value of film was set to 0
DOPC_nr_film = 55 ;%40;
DOPC_nr_bulk = 95 ;%80;

% 0 value of film was set to 0.
DPPC_nr_film = 52; %30;
DPPC_nr_bulk = 71; %65;

DOPC_RefSpec_film = DOPC_REFscan(DOPC_nr_film,:);
DOPC_RefSpec_bulk = DOPC_REFscan(DOPC_nr_bulk,:);

DPPC_RefSpec_film = DPPC_REFscan(DPPC_nr_film,:);
DPPC_RefSpec_bulk = DPPC_REFscan(DPPC_nr_bulk,:);

Data_all = [DOPC_RefSpec_film; DOPC_RefSpec_bulk; DPPC_RefSpec_film; DPPC_RefSpec_bulk; scan_raw]; %ProcessedData]; %; data_scan]; 

%% Remove BG substract and correct 
i_bg = 599; f_bg = 684; %BG

for cn = 1:1:size(Data_all,1);
    BG(1,1) = mean(Data_all(cn,i_bg:f_bg));
    %M(1,cn) = BG(1,1);
    data_all_BG(cn,:) = Data_all(cn,:) - BG(1,1); % BG substraction
end 

%BG correction 
% Define finalx positions
x1 = f_bg;
x2 = size(Data_all,2);
%dx = x2-x1+1; %pluse 1 becuse you have water in ALL data. 
x_vals = (x1:x2)';
y_BASE = zeros(1600,1);

for cz = 1:1:size(Data_all,1);
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

i_CN = 166; f_CN = 181; % C-N normalizaton 

for cn = 1:1:size(Data_all,1);
    
    data_all_CNI(cn,1) = trapz(data_all_BG(cn,i_CN:f_CN)); %taking out integral
    Data_all(cn,:) = data_all_BG(cn,:)./data_all_CNI(cn,1); %normalizing to CN integral

end 

%%

%Separating data not to cause confusion
DOPC_BG_film = data_all_BG(1,:)'; %BG substracted
DOPC_ref_film = Data_all(1,:)'; %Normalized to CN integral
DOPC_totI_film = data_all_CNI(1,:)'; %Integral of spectrum 

DOPC_BG_bulk = data_all_BG(2,:)'; %BG substracted
DOPC_ref_bulk = Data_all(2,:)'; %Normalized to CN integral
DOPC_totI_bulk = data_all_CNI(2,:)'; %Integral of spectrum 

DPPC_BG_film = data_all_BG(3,:)'; %BG substracted
DPPC_ref_film = Data_all(3,:)'; %Normalized to CN integral
DPPC_totI_film = data_all_CNI(3,:)'; %Integral of spectrum 

DPPC_BG_bulk = data_all_BG(4,:)'; %BG substracted
DPPC_ref_bulk = Data_all(4,:)'; %Normalized to CN integral
DPPC_totI_bulk = data_all_CNI(4,:)'; %Integral of spectrum 

%Capillary scan data
Data_BG = data_all_BG(5:end,:)'; %BG substracted
Data_scan = Data_all(5:end,:)'; %Normalized to CN integral
Data_Int_CN = data_all_CNI(5:end,:)'; %Integral of spectrum 

ScanV = 1:1:size(Data_Int_CN,2); 

%% Controll 

figure; %1
plot(x_cm1,Data_scan(:,60:end),'-','linewidth',2);
grid minor
set(gca, 'ColorOrder', cool(size(Data_scan,2)));
%legend('Data BG+Norm')
title('Scan spectrum Normalized to C-N Int');
ylabel('Normalized CCD cts');
xlabel('INDEX [nr.]');
ylim([0 6]);

fInt = figure; %2
plot(ScanV,Data_Int_CN,'-ko','linewidth',1.5);
grid minor;

xlim([0 ScanV(end)]);
%xlim([0 110]);
xticks(0:10:ScanV(end));

%xlabel('Scan Index [nr.]','FontSize', 16);
xlabel('Distance from the capillary tip [\mum]','FontSize', 16);
ylabel('Integral of C-N peak');
%titel('C-N peak integral along the capillary')
yLabelHandle = get(gca, 'YLabel');
set(yLabelHandle, 'Position', get(yLabelHandle, 'Position') + [-1 0 0]);

%print(fInt,'CN_INT_1090', '-dpdf', '-r300','-bestfit'); 
%print(fInt,'CN_INT_1090', '-dsvg');


%% LSQ Fitting Film 
i_fit = 430 ; f_fit = 550; %C=C fitting 
%i_fit = 1051 ; f_fit = 1143; %select fitting area in spectrum
%i_fit = 1051 ; f_fit = 1183; %select fitting area in spectrum

fit_vector = i_fit:1:f_fit; 

%separating references for fitting
DOPC_fit_film = DOPC_ref_film(i_fit:f_fit,1);
DOPC_fit_bulk = DOPC_ref_bulk(i_fit:f_fit,1);

DPPC_fit_film = DPPC_ref_film(i_fit:f_fit,1);
DPPC_fit_bulk = DPPC_ref_bulk(i_fit:f_fit,1);

%exapt_X = 65:90; %=490:520
%DOPC_fit_film(exapt_X,:) = []; 
%DPPC_fit_film(exapt_X,:) = [];

C_film = [DOPC_fit_film, DPPC_fit_film]; C_film = double(C_film); %format neede for lsqlin fit
C_bulk = [DOPC_fit_bulk, DPPC_fit_bulk]; C_bulk = double(C_bulk);

cm = 1;
D_all = [];
D_fit_all =[]; 

for cm = 1:size(Data_scan,2)
    D = Data_scan(i_fit:f_fit,cm); %Data set to fit
    %D(exapt_X,:) = [];
    D_all(:,cm) = D(:,1);
    D = double(D);
    %lin = lsqlin(C, D,[],[],[],[],[0 0],[1 1]); %gives same results as
    lin = lsqlin(C_film, D,[],[],[],[],[0 0],[inf inf]);
    
    s1_lsq_film_raw(cm,1) = lin(1,1); %Coeficents
    s2_lsq_film_raw(cm,1) = lin(2,1); %Coeficents
    
    D_fited = (lin(1,1) .* DOPC_fit_film) + (lin(2,1) .* DPPC_fit_film);
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


%% To controll LSQ fitting film 

%{
scan_nr = 76; %<----- spectra to see

figure
plot(fit_vector,D_all(:,scan_nr),'ko','linewidth',1)
%plot(D_all(:,scan_nr),'ko','linewidth',1)
hold on 
plot(fit_vector,D_fit_all(:,scan_nr),'r-','linewidth',1)
%plot(D_fit_all(:,scan_nr),'r-','linewidth',1)
hold off
grid minor 
legend('Data','LSQ fitt')
ylabel('Normalzied CCD cts');
xlabel('index [nr.]')
xlim([i_fit f_fit])
title(sprintf('LSQ fitting area of scan %d%d',scan_nr));
%}


%% LSQ fitting Bulk 

cm = 1;
for cm = 1:size(Data_scan,2);
    %DOPC bulk fitted
    D = Data_scan(i_fit:f_fit,cm); %Data set to fit
    D_all_bulk(:,cm) = D(:,1);
    
    D = double(D);
    %lin_bulk = lsqlin(C_bulk, D,[],[],[1 1],1,[0 0],[1 1]); %constrained
    lin_bulk = lsqlin(C_bulk, D,[],[],[],[],[0 0],[inf inf]); %constrained
    
    s1_lsq_bulk_raw(cm,1) = lin_bulk(1,1); %Coeficents
    s2_lsq_bulk_raw(cm,1) = lin_bulk(2,1); %Coeficents
    
    D_fited_bulk = (lin_bulk(1,1) .* DOPC_fit_bulk) + (lin_bulk(2,1) .* DPPC_fit_bulk);
    D_fit_all_bulk(:,cm) = D_fited_bulk; %save fitting
    
    %R2 Values
    mean_D = mean(D); 
    TSS = sum((D - mean_D).^2);
    TSS_all_bulk(cm,1) = TSS(1,1);
    RSS = sum((D - D_fited_bulk).^2);
    RSS_all_bulk(cm,1) = RSS;
    R2_value = 1 -(RSS/TSS);
    R2_all_bulk(cm,1) = R2_value;
end 


%% To controll LSQ fitting bulk 
%{
scan_nr = 76; %<----- spectra to see

fCC = figure;
plot(x_cm1(fit_vector),D_all_bulk(:,scan_nr),'ko','linewidth',1)
%plot(D_all(:,scan_nr),'ko','linewidth',1)
hold on 
plot(x_cm1(fit_vector),D_fit_all_bulk(:,scan_nr),'r-','linewidth',1)
%plot(D_fit_all(:,scan_nr),'r-','linewidth',1)
hold off

grid minor 
legend('Data points','LSQ fitting')
ylabel('Normalzied CCD cts');
xlabel('Raman shift [1/cm]')
%xlim([i_fit f_fit])
ylim([0 0.7])
xlim([x_cm1(i_fit) x_cm1(f_fit)])
title(sprintf('Scan nr. %d%d',scan_nr));

%print(fCC,'C=Cfit_1090_40', '-dpdf', '-r300','-bestfit'); 
%print(fCC,'C=Cfit_1090_40', '-dsvg');
%}


%% LSQlin fitting C-H area film 
%i_fitCH = 1051 ; f_fitCH = 1143; %select fitting area in spectrum
i_fitCH = 1055 ; f_fitCH = 1155; %select fitting area in spectrum

fit_vectorCH = i_fitCH:1:f_fitCH ; 

%separating references for fitting
DOPC_fitCH_film = DOPC_ref_film(i_fitCH:f_fitCH,1);
DPPC_fitCH_film = DPPC_ref_film(i_fitCH:f_fitCH,1);

C_filmCH = [DOPC_fitCH_film, DPPC_fitCH_film]; C_filmCH = double(C_filmCH); %format neede for lsqlin fit

cm = 1;
D_CH = [];

for cm = 1:size(Data_scan,2)
    D_CH = Data_scan(i_fitCH:f_fitCH,cm); %Data set to fit
    D_allCH(:,cm) = D_CH(:,1);
    D_CH = double(D_CH);
    
    linCH = lsqlin(C_filmCH, D_CH,[],[],[],[],[0 0],[inf inf]);
    
    s1_lsq_CH(cm,1) = linCH(1,1); %Coeficents
    s2_lsq_CH(cm,1) = linCH(2,1); 
    
    D_fitedCH = (linCH(1,1) .* DOPC_fitCH_film) + (linCH(2,1) .* DPPC_fitCH_film);
    D_fitCH_all(:,cm) = D_fitedCH(:,1); %save fitting 
    
    %R2 Values
    mean_D_CH = mean(D_CH); 
    TSS_CH = sum((D_CH - mean_D_CH).^2);
    TSS_all_CH(cm,1) = TSS_CH(1,1);
    RSS_CH = sum((D_CH - D_fitedCH).^2);
    RSS_all_CH(cm,1) = RSS_CH;
    R2_value_CH = 1 -(RSS_CH/TSS_CH);
    R2_all_CH(cm,1) = R2_value_CH;   
end 

% To controll LSQ fitting bulk 

%{
scan_nr = 58; %<----- spectra to see

fCH = figure;
plot(x_cm1(fit_vectorCH),D_allCH(:,scan_nr),'ko','linewidth',1)
hold on 
plot(x_cm1(fit_vectorCH),D_fitCH_all(:,scan_nr),'r-','linewidth',1)
hold off
grid minor 
legend('Data points','LSQ fitting')
ylabel('Normalzied CCD cts');
xlabel('Raman shift [1/cm]')
%xlabel('Index')
%xlim([i_fitCH f_fitCH])
xlim([x_cm1(i_fitCH) x_cm1(f_fitCH)])
title(sprintf('Scan nr. %d%d',scan_nr));

%print(fCH,'CHfit_1090_58', '-dpdf', '-r300','-bestfit'); 
%print(fCH,'CHfit_1090_58', '-dsvg');
%}


%% PULK fitt CH 

DOPC_fitCH_bulk = DOPC_ref_bulk(i_fitCH:f_fitCH,1);
DPPC_fitCH_bulk = DPPC_ref_bulk(i_fitCH:f_fitCH,1);

C_filmCH_bulk = [DOPC_fitCH_bulk, DPPC_fitCH_bulk]; C_filmCH_bulk = double(C_filmCH_bulk); %format neede for lsqlin fit

for cm = 1:size(Data_scan,2)
    D_CH = Data_scan(i_fitCH:f_fitCH,cm); %Data set to fit
    D_allCH(:,cm) = D_CH(:,1);
    D_CH = double(D_CH);
    
    linCH_bulk = lsqlin(C_filmCH_bulk, D_CH,[],[],[],[],[0 0],[inf inf]);
    
    s1_lsq_CH_bulk(cm,1) = linCH_bulk(1,1); %Coeficents
    s2_lsq_CH_bulk(cm,1) = linCH_bulk(2,1); 
    
    D_fitedCH_bulk = (linCH_bulk(1,1) .* DOPC_fitCH_bulk) + (linCH_bulk(2,1) .* DPPC_fitCH_bulk);
    D_fitCH_all_bulk(:,cm) = D_fitedCH_bulk(:,1); %save fitting 
    
    %R2 Values
    mean_D_CH = mean(D_CH); 
    TSS_CH = sum((D_CH - mean_D_CH).^2);
    TSS_all_CH(cm,1) = TSS_CH(1,1);
    RSS_CH = sum((D_CH - D_fitedCH_bulk).^2);
    RSS_all_CH(cm,1) = RSS_CH;
    R2_value_CH = 1 -(RSS_CH/TSS_CH);
    R2_all_CH_bulk(cm,1) = R2_value_CH;   
end 



%% Plot REF you use to FIT 

figCC = figure; %B
plot(x_cm1, DOPC_ref_film,'-','color',"#0072BD",'linewidth',2);
hold on;
plot(x_cm1, DOPC_ref_bulk,'o-','color',"#0072BD",'linewidth',2);
plot(x_cm1, DPPC_ref_film,'-','color',"#D95319",'linewidth',2);
plot(x_cm1, DPPC_ref_bulk,'o-','color',"#D95319",'linewidth',2);
%plot([x_cm1(i_fit) x_cm1(f_fit)], [0 0],'k-','linewidth',0.5);
grid minor;
hold off;
xlim([x_cm1(i_fit) x_cm1(f_fit)]);
ylim([0 1]);
%legend('DOPC nr.55','DOPC nr.92','DPPC nr.50','DPPC nr.31');
ylabel('Normalzied CCD cts');
xlabel('Raman shift [1/cm]');
%print(figCC,sprintf('RefFit_1400-1700'), '-dpdf', '-r300','-bestfit'); 
%print(figCC,sprintf('RefFit_1400-1700'), '-dsvg');

fig = figure; %C
plot(x_cm1, DOPC_ref_film,'-','color',"#0072BD",'linewidth',2);
hold on;
plot(x_cm1, DOPC_ref_bulk,'o-','color',"#0072BD",'linewidth',2);
plot(x_cm1, DPPC_ref_film,'-','color',"#D95319",'linewidth',2);
plot(x_cm1, DPPC_ref_bulk,'o-','color',"#D95319",'linewidth',2);
%plot([x_cm1(i_fitCH) x_cm1(f_fitCH)], [0 0],'k-','linewidth',0.5);
grid minor;
hold off;
xlim([x_cm1(i_fitCH) x_cm1(f_fitCH)]);
ylim([0 6]);
%legend('DOPC film','DOPC bulk','DPPC film','DPPC bulk');
%legend('DOPC nr.55','DOPC nr.92','DPPC nr.50','DPPC nr.31');
ylabel('Normalzied CCD cts');
xlabel('Raman shift [1/cm]');
%title('Refernce spectrum PHASE fitting')
%lgd = legend;
%lgd.Position = [0.68, 0.75, 0.2, 0.15];

%print(fig,sprintf('RefFit_2800-3000'), '-dpdf', '-r300','-bestfit'); 
%print(fig,sprintf('RefFit_2800-3000'), '-dsvg');

%% WT Correction

wt_dopc = 0.9905; %ml/g
wt_dppc = 0.941; %ml/g
r_wt = wt_dopc/wt_dppc; 

s1_lsq_film = s1_lsq_film_raw .* r_wt;
s2_lsq_film = s2_lsq_film_raw ./ r_wt;

s1_lsq_bulk = s1_lsq_bulk_raw .* r_wt; 
s2_lsq_bulk = s2_lsq_bulk_raw ./ r_wt; 

s1_AVR = (s1_lsq_film + s1_lsq_bulk)./2;
s2_AVR = (s2_lsq_film + s2_lsq_bulk)./2;

%Coeficients
coef1_sumN_film = s1_lsq_film(:,1)./(s1_lsq_film(:,1)+s2_lsq_film(:,1));
coef1_sumN_film = coef1_sumN_film .*100; %Coef corrected for wt%
coef2_sumN_film = 100 - coef1_sumN_film; 
%coef2_sumN_film = s2_lsq_film./(s1_lsq_film(:,1)+s2_lsq_film(:,1));
%coef2_sumN_film = coef2_sumN_film.*100;

coef1_sumN_bulk = s1_lsq_bulk(:,1)./(s1_lsq_bulk(:,1)+s2_lsq_bulk(:,1));
coef1_sumN_bulk = coef1_sumN_bulk.*100; %Coef corrected for wt%
coef2_sumN_bulk = 100 - coef1_sumN_bulk; 

%coef2_sumN_bulk = s2_lsq_bulk./(s1_lsq_bulk(:,1)+s2_lsq_bulk(:,1));
%coef2_sumN_bulk = coef2_sumN_bulk.*100;

%AVREAGE
Av_c1 = (coef1_sumN_film + coef1_sumN_bulk)./2; 
Av_c2 = (coef2_sumN_film + coef2_sumN_bulk)./2; 


%PHASE FITTING 
coef_CH1 = s1_lsq_CH(:,1) ./ (s1_lsq_CH(:,1) + s2_lsq_CH(:,1)); 
coef_CH1 = coef_CH1.*100;
coef_CH2 = s2_lsq_CH(:,1) ./ (s1_lsq_CH(:,1) + s2_lsq_CH(:,1));
coef_CH2 = coef_CH2.*100;

coef_CH1_bulk = s1_lsq_CH_bulk(:,1) ./ (s1_lsq_CH_bulk(:,1) + s2_lsq_CH_bulk(:,1)); 
coef_CH1_bulk = coef_CH1_bulk.*100;
coef_CH2_bulk = s2_lsq_CH_bulk(:,1) ./ (s1_lsq_CH_bulk(:,1) + s2_lsq_CH_bulk(:,1));
coef_CH2_bulk = coef_CH2_bulk.*100;


Phase1_AVR = (coef_CH1 + coef_CH1_bulk)./2;
Phase2_AVR = (coef_CH2 + coef_CH2_bulk)./2;

%% SAVING DATA 


% PUT in to table
Table_coef_all(:,1) = coef1_sumN_film;
Table_coef_all(:,2) = coef2_sumN_film;
Table_coef_all(:,3) = coef1_sumN_bulk;
Table_coef_all(:,4) = coef2_sumN_bulk;
Table_coef_all(:,5) = Av_c1;
Table_coef_all(:,6) = Av_c2;
Table_coef_all(:,7) = i_scan:1:f_scan;
Table_coef_all(:,8) = coef_CH1;
Table_coef_all(:,9) = coef_CH1_bulk;
Table_coef_all(:,10) = Phase1_AVR;
Table_coef_all(:,11) = Data_Int_CN;

Table_allcoef1090_New = Table_coef_all; 
save('Table_allcoef1090_New.mat','Table_allcoef1090_New'); %<----- open


%{
Table_shortR2(:,1) = R2_all_film;
Table_shortR2(:,2) = R2_all_bulk;
Table_shortR2(:,3) = R2_all_CH;
Table_shortR2(:,4) = i_scan:1:f_scan;

save('Table_shortR2.mat','Table_shortR2'); %<----- open
%}

%% PLOTTING DATA ----------------------

% Make X vector  
X_initial_film = 20;
step = 20;
X_final_film = size(coef1_sumN_film,1); 

ScanV = []; 
ScanV = 1:1:X_final_film; 
ScanV = ScanV' ;
ScanV_corr = ScanV - X_initial_film; 

% Plottig lipid ratio 
fig = figure;
hold on ;
p1 = plot(ScanV_corr,coef1_sumN_film(:,1),'-','color',[0.7 0.7 0.7],'linewidth',1);
p2 = plot(ScanV_corr,Av_c1(:,1),'r-','linewidth',2);
p3 = plot(ScanV_corr,coef1_sumN_bulk(:,1),'-','color',[0.5 0.5 0.5],'linewidth',1);
%p4 = plot(smooth_vector(:,rep_i:rep_f)-X_initial_film, Av_c1_smooth(rep_i:rep_f,:),'-','color',[0.6 0.6 0.6],'linewidth',2);

p4 = plot(ScanV_corr,coef2_sumN_film(:,1),'-','color',[0.7 0.7 0.7],'linewidth',1);
p5 = plot(ScanV_corr,Av_c2(:,1),'b-','linewidth',2);
p6 = plot(ScanV_corr,coef2_sumN_bulk(:,1),'-','color',[0.5 0.5 0.5],'linewidth',1);
%p7 = plot(smooth_vector(:,rep_i:rep_f)-X_initial_film, Av_c2_smooth(rep_i:rep_f,:),'-','color',[0.6 0.6 0.6],'linewidth',2);

hold off ;

grid on ;
set(gca,'ycolor','b');
ylim([0 100]);

%xlim([-1*X_initial_film ScanV_corr(end)]);
xlim([-10 ScanV_corr(end)]);
xticks(ScanV_corr(1)-1:10:ScanV_corr(end));

lgd1 = legend([p2,p5],'DOPC Avr.','DPPC Avr.','location','northoutside','Orientation','horizontal');
ylabel('Lipid ratio [wt%]','FontSize', 16);
xlabel('Distance from the capillary tip [\mum]','FontSize', 16);

%print(fig,sprintf('Lipidfit_1090_LLSS'), '-dpdf', '-r300','-bestfit'); 
%print(fig,sprintf('Lipidfit_1090'), '-dsvg');


% Plottig lipid ratio 
fig = figure;
plot(ScanV_corr,coef_CH1_bulk(:,1),'k-.','linewidth',1);
hold on ;
p2 = plot(ScanV_corr,Phase1_AVR(:,1),'ko-','linewidth',1.5);
p3 = plot(ScanV_corr,coef_CH1(:,1),'k--','linewidth',1);
hold off ;

grid on ;
%set(gca,'ycolor','b');
ylim([0 100]);
yticks(0:25:100);

%xlim([-1*X_initial_film ScanV_corr(end)]);
xlim([-10 ScanV_corr(end)]);
xticks(ScanV_corr(1)-1:10:ScanV_corr(end));

lgd1 = legend([p2],'DOPC Avr.','PHASE.','location','northoutside','Orientation','horizontal');
ylabel('Lipid Phase L\alpha vs L\beta [%]','FontSize', 16);
xlabel('Distance from the capillary tip [\mum]','FontSize', 16);
%print(fig,sprintf('CHfit_1090'), '-dsvg');


% Plot All r2 values ####################################################
fr = figure;
p1 = plot(ScanV_corr,R2_all_film,'o-','color',[0.2 0.6 0.2],'linewidth',1.5);
hold on 
p2 = plot(ScanV_corr,R2_all_bulk,'*-','color',[0.2 0.6 0.2],'linewidth',1.5);
p3 = plot(ScanV_corr,R2_all_CH,'ok-','linewidth',1.5);
p3 = plot(ScanV_corr,R2_all_CH_bulk,'*k-','linewidth',1.5);

%xlim([1 size(R2_all_CH,1)]);
hold off 
grid on ;
legend('Ref. Film','Ref. Bulk','PHASE film','PHASE bulk','location','southeast');
%legend([p1 p3],'CC Area fitt','CH Area fitt','location','southeast');

%xlabel('Distance from the capillary tip [\mum]','FontSize', 16);
xlabel('Distance from the capillary tip [\mum]','FontSize', 16);
ylabel('R^2 coefficient')

ylim([0.5 1]);
xlim([-10 ScanV_corr(end)]);
xticks(ScanV_corr(1)-1:10:ScanV_corr(end));

yLabelHandle = get(gca, 'YLabel');
set(yLabelHandle, 'Position', get(yLabelHandle, 'Position') + [-1 0 0]);

%print(fr,'R2_fit_1090_LLSS', '-dpdf', '-r300','-bestfit'); 
%print(fr,'R2_fit_1090', '-dsvg');

%% Smoothening of and saving data 
%{
%You cannot smooth RAW for CC and CH fitting data becuse phase will be affected.
%Smooth coeficient cureve insted

rep_i = 1; rep_f = 117;
Av_c1_smooth = smoothdata(Av_c1,1,"movmean",10); 
Av_c1_Rep = Av_c1(:,:);
Av_c1_Rep(rep_i:rep_f,:) = Av_c1_smooth(rep_i:rep_f,:);


fig = figure; %C
plot(ScanV_corr,Av_c1(:,1),'-','color',[0.7 0.7 0.7],'linewidth',2);
hold on 
plot(ScanV_corr, Av_c1_smooth,'r-','linewidth',2);
plot(ScanV_corr, Av_c1_Rep,'b-','linewidth',2);
grid minor;
hold off;
ylabel('C1');
xlabel('[Spec nr]');
%title('Smooth data')
legend('RAW','Smooth 10','Combined')
xlim([-20 ScanV_corr(end)]);

%}

%% 
%{
% SAVE in to table
Table_allcoef_Smooth(:,1) = Av_c1;
Table_allcoef_Smooth(:,2) = Av_c1_smooth;
Table_allcoef_Smooth(:,3) = Av_c1_Rep;
Table_allcoef_Smooth(:,5) = ScanV_corr;
Table_allcoef_Smooth(:,5) = ScanV;
save('Table_allcoef_Smooth.mat','Table_allcoef_Smooth'); %<----- open

%}

%% OLD smooting 
%{
smooth_window = 5;
smooth_vector = 1:1:size(Av_c1,1) - smooth_window ;

Av_c1_smooth = []; 
for dr = 1:1:size(smooth_vector,2)
    Av_c1_smooth(dr,1) = mean(Av_c1(dr:dr+smooth_window,1),1);
    Av_c2_smooth(dr,1) = mean(Av_c2(dr:dr+smooth_window,1),1);

end

%Replace A1 part with smoothed
rep_i = 1; rep_f = 117;
Av_c1_Rep = Av_c1(:,:);
Av_c1_Rep(rep_i:rep_f,:) = Av_c1_smooth(rep_i:rep_f,:);
Av_c2_smooth = 100 - Av_c1_smooth(:,1);

%
fig = figure; %C
plot(ScanV_corr,Av_c1(:,1),'-','color',[0.6 0.6 0.6],'linewidth',2);
hold on 
plot(smooth_vector-X_initial_film, Av_c1_smooth,'r-','linewidth',2);
plot(ScanV_corr, Av_c1_Rep,'b-','linewidth',2);
grid minor;
hold off;
ylabel('C1');
xlabel('[Spec nr]');
%title('Smooth data')
xlim([-20 ScanV_corr(end)]);

%}

%% PLOT CC AND CH 
%{
X_init = -10; %CUT img 
X_final = 110; 
xTicks = 0:10:120;
yTicks = 0:50:100;

%Define relative axis 
film_end = 90 ; %<--------- Select final film position
film_ini = 0 ; %<--------- Select initial position 
size_ofFilm = film_end - film_ini; 

tickVx2 = linspace(film_ini,film_end,6); %<----- SCALE
tickslable_conv = tickVx2 - film_ini; %DONT USE
tickslable_rel = (tickslable_conv.*100)./tickslable_conv(end); 


fig = figure;
fig.Position = [200 200 600 300];
ax1 = axes(fig);

%Left side 
yyaxis right
set(gca,'ycolor','k') ;
p10 = plot(ScanV_corr,coef_CH1(:,1),'ko-','linewidth',1.5);

xlim([X_init X_final]);
xticks(xTicks);
ylim([0 100]);
yticks(yTicks);

%yticklabels({'0:1','2:8','4:6','6:4','8:2','1:0'});
yticklabels({'0:1','1:1','1:0'});
set(gca,'TickDir','out');
%set(gca,'yticklabel',[]);

ylabel('Phase   L\alpha:L\beta','FontSize', 18);
yLabelHandle = get(gca, 'YLabel');
set(yLabelHandle, 'Position', get(yLabelHandle, 'Position') + [1 0 0]);


%Right side 
yyaxis left ;
set(gca,'ycolor','r');

p2 = plot(ScanV_corr,Av_c1_Rep(:,1),'r-','linewidth',2);

grid on;
set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
%set(gca,'TickDir','out');

xlim([X_init X_final]);
xticks(xTicks);
ylim([0 100]);
yticks(0:20:100);

ylabel('DOPC/DOPC+DPPC [wt%]','FontSize', 16);
yLabelHandle = get(gca, 'YLabel');
set(yLabelHandle, 'Position', get(yLabelHandle, 'Position') + [0 0 0]);

xlabel('Distance from the capillary tip [\mum]','FontSize', 18);
xLabelHandle = get(gca, 'XLabel');
set(xLabelHandle, 'Position', get(xLabelHandle, 'Position') + [0 0 0]);

%legend([p5 p10],'Smooth','PHASE fit')

ax1.Box = 'off';
xline(ax1,ax1.XLim(2));
%{
%Relative Top AXIS --------

ax2 = axes('Position', get(ax1, 'Position'),'Color', 'none');
set(ax2, 'XAxisLocation', 'top','YAxisLocation','left');
%set(ax2,'Xtick',[],'Ytick',[]) % To hide ax2 ticks
set(gca,'XMinorTick','on')
set(gca,'TickDir','out');
ax2.YAxis.Visible = 'off';

step_rel = tickVx2(2) - tickVx2(1);
ax2.XAxis.MinorTickValues = tickVx2(1):step_rel/2:tickVx2(end);

set(ax2,'XLim',[[X_init X_final]]...
    ,'xtick',tickVx2,'xticklabels',string(tickslable_rel));

%xlabel(ax2,'[\mum]','FontSize', 12);
%xLabelHandle = get(gca, 'XLabel');
%set(xLabelHandle, 'Position', get(xLabelHandle, 'Position') + [0 0.02 0]);

%}
%Size=[Left  Bottom  Width   Height]
%Size=[0.125  0.185  0.76  0.70]; %[0.1300 0.1100 0.7750 0.8150] (default, ) 
Size=[0.125  0.185  0.76  0.70]; %[0.1300 0.1100 0.7750 0.8150] (default, ) 
set(ax1,'InnerPosition',Size);
%set(ax2,'InnerPosition',Size);

%print(fig,sprintf('LipidPhase_1090_Fig2_REL_A'), '-dpdf', '-r300','-bestfit'); 
print(fig,sprintf('LipidPhase_1090_Fig2_REL_A'), '-dsvg'); 
%}

%% INTEGRALs 

%{
%X_initial_film = 60; %ScanV_corr is shifted by 69 bit not the data points. 
region_Ai = X_initial_film; 
region_Af = ScanV(end); 
Integral_DOPCA = trapz(Av_c1_Rep(region_Ai:region_Af,1));

bas_fini = 30; %How far from the back! 
Baseline = mean(Av_c1_Rep(end-bas_fini:end,1));

Integral_Baseline = Baseline.* (region_Af-region_Ai); 
Integral_Total = 100.* (region_Af-region_Ai); 

IntRatio_10wt = Integral_Baseline.*100./Integral_Total;
IntRatio_fitted = Integral_DOPCA.*100./Integral_Total;

fig = figure;
hold on ;
%p1 = plot(ScanV(:,1), Av_c1_Rep(:,1),'ko-','linewidth',2);
p2 = plot(ScanV(region_Ai:region_Af,1), Av_c1_Rep(region_Ai:region_Af,1),'ro-','linewidth',2);
p3 = plot([ScanV(end-bas_fini) ScanV(end)],[Baseline Baseline],'b-','linewidth',3);
hold off ;
legend([p3,p2],sprintf('Baseline %.3g wt%%',Baseline),sprintf('Integral %.3g %%',IntRatio_fitted),'location','northwest','FontSize', 16);
grid on ;
ylabel('DOPC/(DOPC+DPPC) [wt%]','FontSize', 16);
xlabel('[Spec nr]');
ylim([0 100]);
xlim([region_Ai region_Af]);

%print(fig,'Integral_1090', '-dsvg');
%}

%{
Baseline = 11;
DOPC_0 = Av_c1; % - Baseline;

region_Ai = 60; 
region_Af = 133;
region_Bi = 133;
region_Bf = 190; 

DOPC_regionA = DOPC_0(region_Ai:region_Af,1);
DOPC_regionB = DOPC_0(region_Bi:region_Bf,1);

Integ_DOPCA = trapz(DOPC_regionA);
Integ_DOPCB = trapz(DOPC_regionB);

Ratio = Integ_DOPCA*100./(Integ_DOPCA+Integ_DOPCB);

fig_0 = figure;
hold on ;
p2 = plot(region_Ai:1:region_Af,DOPC_regionA,'ro-','linewidth',2);
p3 = plot(region_Bi:1:region_Bf,DOPC_regionB,'b-','linewidth',2);
p3 = plot([region_Ai region_Bf],[Baseline Baseline],'k-','linewidth',3);
hold off ;

legend(sprintf('A %.5g',Integ_DOPCA),sprintf('B %.5g',Integ_DOPCB),'location','northwest','FontSize', 16);
grid on ;
set(gca,'ycolor','b');
xlim([-20 ScanV_corr(end)]);

xlim([region_Ai region_Bf]);
%xticks(Vector_X);
%xticklabels(Vector_Xle); 

ylabel('Lipid ratio [%]','FontSize', 16);
xlabel('Distance from the capillary tip [\mum]','FontSize', 16);

%print(fig_0,'Integral_1090', '-dpdf', '-r300','-bestfit'); 
%print(fig_0,'Integral_1090', '-dsvg');

%}

%% Taking out integral under DOPC curve
%{

zero_int = 60; 
final_int = 190 ;%ScanV(end);
Vector_X = zero_int:10:final_int;
Vector_Xle = Vector_X - zero_int; 

DOPC_vector = Av_c1(zero_int:final_int,1);
Integ_DOPCcurve = trapz(DOPC_vector);

totaArea = (final_int - zero_int)*100; 
Int_Prop = Integ_DOPCcurve.* 100./totaArea;

fig = figure;
hold on ;
p2 = plot(zero_int:1:final_int,DOPC_vector,'r-','linewidth',2);
hold off ;

legend(sprintf('Integral %.3g %%',Int_Prop),'location','northwest','FontSize', 16);
grid on ;
set(gca,'ycolor','b');
ylim([0 100]);

xlim([zero_int final_int]);
xticks(Vector_X);
xticklabels(Vector_Xle); 

ylabel('Lipid ratio [%]','FontSize', 16);
xlabel('Distance from the capillary tip [\mum]','FontSize', 16);

%print(fig,'Integral_1090', '-dpdf', '-r300','-bestfit'); 
%print(fig,'Integral_1090', '-dsvg');

%}
%}
