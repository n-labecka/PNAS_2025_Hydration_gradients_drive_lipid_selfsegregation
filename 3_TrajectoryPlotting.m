clc 
clear all
close all 

addpath '/Users/labecka/Documents/3.Lab Work/0.LU_DOPCDPPC studdy/Raman data/Final Fitting IMG/2022_DOPCDPPC_1090/TdiagramFunctions'
warning off MATLAB:griddata:DuplicateDataPoints

set(0, 'DefaultAxesFontSize', 14);

%% Files 

%cd 
%loading files 
files = dir('*.mat');

for i=1:length(files)
    load(files(i).name, '-mat');
end

%% TAKE RAW data for water fitting

name_scan = Struct_shortscanGspatGspecwrongspatLine026SpecData1; %RAW spectra becuse you wantTotI normalized
data_scan_RAW = name_scan.data;

name_water = Struct_waterSpectrum042SpecData1 ; 
Water_spec = name_water.data; 

%Scan select data 
i_scan = Table_allcoef1090_New(1,7); %70; % start scan
f_scan = Table_allcoef1090_New(end,7); %190; %size(data_scan,1) ; % final scan
data_all_raw = data_scan_RAW(i_scan:f_scan,:);

x_nm = name_scan.axisscale{2,1};
laser = 531.98 ; %nm
x_cm1 = 10^7 .*(1./laser - 1./x_nm); %Raman shift cm-1


%%
i_bg = 594; f_bg = 680; %BG

data_all = [Water_spec; data_all_raw]; 

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
%%
%data_all_BG(40:70, 688:732) = data_all_BG(40:70, 750:794);

data_all_BG(15:17,200:210) = data_all_BG(30:32,200:210);
%data_all_BG(15:17,:) = data_all_BG(19:,:);

data_all_BG(15:30,86:104) = data_all_BG(35:50,86:104);


%% Check water in raw spectra 

i_water = 1200 ; f_water = 1550; 
%Water_Int_vector = i_water:f_water; 
%i_fitCH = 1055 ; f_fitCH = 1155; %select fitting area in spectrum

for ck = 1:size(data_all_BG,1);
    Water_Integral_BGControll(ck,1) = trapz(data_all_BG(ck,i_water:f_water));    
end

Water_Integral_BGControll = Water_Integral_BGControll(2:end,:);
figure 
plot(Water_Integral_BGControll,'linewidth',1.5)
ylim([-10 50]);
%xlim([60 110]);

Av_c1_smooth1 = smoothdata(data_all_BG,2,"movmean",15); 

figure 
plot(x_cm1,data_all_BG(:,:));
set(gca, 'ColorOrder', jet);

figure 
plot(x_cm1,Av_c1_smooth1(:,:));
set(gca, 'ColorOrder', jet);

%% Normalizide to total integral 
i_tot = 150; f_tot = 1550; % Tot normalizaton 
%i_tot2 = 1046; f_tot2 = 1550; % Tot normalizaton 

i_CN = 166; f_CN = 181; % C-N normalizaton 

for cn = 1:1:size(data_all,1);
    dataValue_CNI(cn,1) = trapz(data_all_BG(cn,i_CN:f_CN)); %taking out integral

    dataValue_totI(cn,1) = trapz(data_all_BG(cn,i_tot:f_tot)); %taking out integral
    data_all_Norm(cn,:) = data_all_BG(cn,:)./dataValue_totI(cn,1); %normalizing to CN integral
end 

%{
figure  
plot(dataValue_totI,'-ro','linewidth',1.5);
hold on 
plot(dataValue_totI2+50,'-bo','linewidth',1.5);
hold off
xlim([50 180])
%}

%% Normalization Correction 

i_tot = 150; f_tot = 1550; % Tot normalizaton 
value = 19+1; 

Norm_value = trapz(data_all_BG(value,i_tot:f_tot)); 
dataValue_MIX = [];
dataValue_totI_MIX_A = [];
dataValue_totI_MIX_B = [];

for cn = 1:1:value; 
    dataValue_totI_MIX_A(cn,:) = data_all_BG(cn,:)./Norm_value; %normalizing to CN integral
end 

Norm_MIX_A = repmat(Norm_value,value,1);

for cm = value+1:1:size(data_all,1);
    
    dataValue_totI_B = trapz(data_all_BG(cm,i_tot:f_tot)); %taking out integral
    Norm_MIX_B(cm-value,:) = dataValue_totI_B; 
    dataValue_totI_MIX_B(cm-value,:) = data_all_BG(cm,:)./dataValue_totI_B; %normalizing to CN integral

end 

dataValue_MIX = [dataValue_totI_MIX_A;dataValue_totI_MIX_B]; 
%dataValue_MIX = dataValue_MIX(2:end,1);
Norm_MIX = [Norm_MIX_A;Norm_MIX_B];
Norm_MIX = Norm_MIX(2:end,1);


%cct_Norm = Norm_MIX(value-1)

%%
%Water_BG = data_all_BG(1,:)'; %BG substracted
Water_Norm = data_all_Norm(1,:)';%Normalized to tot integral
Data_CNIvalue = dataValue_CNI(2:end,:)';

Data_BG = data_all_BG(2:end,:)'; %BG substracted
Data_totIvalue = dataValue_totI(2:end,:)';

Data_scan = data_all_Norm(2:end,:)'; %Normalized to tot integral
Data_scan_MIX = dataValue_MIX(2:end,:)'; %Normalized to tot integral

%ScanV = 1:1:size(Data_scan,2); 
ScanV = linspace(0,140,170);
%%
figure
plot(x_cm1, Data_scan(:,:),'-','linewidth',1.5);
set(gca, 'ColorOrder', jet);

figure
plot(x_cm1, Data_BG(:,:),'-','linewidth',1.5);
set(gca, 'ColorOrder', jet);

figure
plot(x_cm1, Data_scan_MIX(:,15:120),'-','linewidth',1.5);
set(gca, 'ColorOrder', jet);

FR = figure;
FR.Position = [300 300 600 340]; 
plot(x_cm1, Data_scan_MIX(:,15:120),'-','linewidth',1.5);
%plot(x_cm1, Data_scan_MIX(:,180:-1:60),'-','linewidth',1.5);
%set(gca, 'ColorOrder', jet);
set(gca, 'ColorOrder', cool(size(Data_scan_MIX(:,15:120),2)));
%set(gca, 'ColorOrder', jet);
grid on 
xlim([x_cm1(1) x_cm1(end)]);
%ylim([x_cm1(1) x_cm1(end)]);

c = colorbar;
colormap(cool(size(Data_scan_MIX(:,15:120),2)));
ticks_lableBar = 15:15:120; %size(Data_scan_BG,2);
ticks_lableBar = ticks_lableBar';
ticks_scaleBar = linspace(0,1,size(ticks_lableBar,1));
set(c, 'Ticks', ticks_scaleBar);
set(c, 'TickLabels', num2str(ticks_lableBar));

%print(FR,'Raman_1090', '-dsvg');


%% CHECK Norm Factor 

Norm_MIX_AA = Norm_MIX'; 

Norm_F = []; 
for sm = 1:1:170;
    Norm_F(1,sm) = (Data_totIvalue(1,sm) + Norm_MIX_AA(1,sm))./2; 
end 

fC = figure;
fC.Position = [300 300 600 340]; 
plot(ScanV,Data_totIvalue,'o-','color',[0 0.6 0],'linewidth',1.5);
hold on 
%p2 = plot(Norm_MIX,'ko-','linewidth',1.5);
p3 = plot(ScanV,Norm_F,'ko-','linewidth',1.5);

hold off
xlim([0 120]);
grid on 
legend([p3],'Norm. factor','location','northwest');
set(gca,'TickDir','out');
set(gca,'XMinorTick','on');
ylabel('CCD cts');
xlabel('[μm]');
%print(fC,'Normalization factor', '-dsvg');


%%

pos_1 = 20; %cct_1 = Data_totIvalue(pos_1);
pos_2 = 19;
pos_3 = 15; 
pos_4 = 8; 

pos_norm = value-1; 

fA = figure;
%plot(x_cm1, Data_scan(:,pos_1),'-','color',[0.000, 0.447, 0.741],'linewidth',1.5);
hold on
%plot(x_cm1, Data_scan_MIX(:,pos_re),'-','color',[0.000, 0.00, 0.000],'linewidth',1);
p4 = plot(x_cm1, Data_scan(:,pos_4),'-','linewidth',1.5);
p3 = plot(x_cm1, Data_scan(:,pos_3),'-','linewidth',1.5);
p2 = plot(x_cm1, Data_scan(:,pos_2),'-','linewidth',1.5);
p1 = plot(x_cm1, Data_scan(:,pos_1),'-','linewidth',1.5);
hold off
%set(gca, 'ColorOrder', hsv(8).*0.9);
legend([p1 p2 p3 p4], sprintf('nr. %d',pos_1),sprintf('nr. %d',pos_2),sprintf('nr. %d',pos_3),sprintf('nr. %d',pos_4));
xlim([2700 3800]);
ylim([-0.001 0.025]);
title('Normalized to total integral (cct)');
%print(fA,'Normalized to total integral', '-dsvg');


fB = figure;
%plot(x_cm1, Data_BG(:,pos_norm),'-','color',[0.850, 0.325, 0.098],'linewidth',1.5);
hold on 
%plot(x_cm1, Data_BG(:,pos_re),'-','color',[0.000, 0.447, 0.741],'linewidth',1.5);
p4 = plot(x_cm1, Data_BG(:,pos_4),'-','linewidth',1.5);
p3 = plot(x_cm1, Data_BG(:,pos_3),'-','linewidth',1.5);
p2 = plot(x_cm1, Data_BG(:,pos_2),'-','linewidth',1.5);
p1 = plot(x_cm1, Data_BG(:,pos_1),'-','linewidth',1.5);
hold off
legend([p1 p2 p3 p4], sprintf('nr. %d',pos_1),sprintf('nr. %d',pos_2),sprintf('nr. %d',pos_3),sprintf('nr. %d',pos_4));
hold off
xlim([2700 3800]);
ylim([-1 40]);
title('Raw-data');
%print(fB,'Raw Data', '-dsvg');


%% Taking out water integral and translating

i_water = 1200 ; f_water = 1550; 
%Water_Int_vector = i_water:f_water; 
%i_fitCH = 1055 ; f_fitCH = 1155; %select fitting area in spectrum

for ck = 1:size(Data_scan,2)
    Water_Integral(ck,1) = trapz(Data_scan(i_water:f_water,ck));    
    Water_Integral_MIX(ck,1) = trapz(Data_scan_MIX(i_water:f_water,ck));

end

REF_integral = trapz(Water_Norm(i_water:f_water,1));
Water_Integral_NormToREF = Water_Integral.*100./REF_integral; 
Water_Integral_NormToREF(Water_Integral_NormToREF <= 0) = eps;

Water_Integral_MIX = Water_Integral_MIX.*100./REF_integral; 
Water_Integral_MIX(Water_Integral_MIX <= 0) = eps;

%% CalibCurve Fitting 
%{
a1 = 36.93;
b1 = 0.2748; 
c1 = -30.66;
%}

a1 = 39.52  ;
b1 = 0.2629  ;
c1 = -32.9 ;

Y_fitted_water = a1.*Water_Integral_NormToREF.^b1 +c1; 
Y_fitted_water(Y_fitted_water < 0) = 0;
Y_fitted_water(Y_fitted_water > 100) = 100;

Y_fitted_water_MIX = a1.*Water_Integral_MIX.^b1 +c1; 
Y_fitted_water_MIX(Water_Integral_MIX < 0) = 0;
Y_fitted_water_MIX(Water_Integral_MIX > 100) = 100;


%% Raw Water profile 

set(0, 'DefaultAxesFontSize', 14);
X_Initial = 0; %60;
X_vector = i_scan:1:f_scan; %size(Data_totWaterI,1);
X_vector = X_vector' ; 
X_ticks = 0:10:size(Water_Integral,1);
X_ticks_corr = X_ticks - X_Initial; 


f4 = figure;
f4.Position = [300 300 600 300]; 

plot(X_vector,Water_Integral_NormToREF,'o-','color',[0.000, 0.447, 0.741],'linewidth',1);
hold on 
plot(X_vector,Water_Integral_MIX,'o-','color',[0.850, 0.325, 0.098],'linewidth',1);
hold off
grid minor;
xlim([X_vector(1) X_vector(end)]);
ylim([0 100]);
%legend('a.*exp(b.*x) + c.*exp(d.*x)','Data points','location','southeast')
title('Water profile raw');
ylabel('Water peak Integral [% of total]');
xlabel('Spectrum [nr.]');
%print(f4,'WaterInt_RAW','-dpng');


% Fitted water profile 
f5 = figure;
f5.Position = [300 300 600 300]; 
%ax1 = axes;
plot(X_vector(:,1),Y_fitted_water(:,1),'o-','color',[0.000, 0.447, 0.741],'linewidth',1.5);
hold on 
plot(X_vector,Y_fitted_water_MIX,'o-','color',[0.850, 0.325, 0.098],'linewidth',1);
hold off
ylim([0 100]);
yticks(0:20:100);

%xlim([40 X_vector(end)]);
xlim([X_vector(1) X_vector(end)]);
xticks(X_ticks);
xticklabels(X_ticks_corr);

grid on
%set(gca,'XMinorTick','on')
%set(gca,'YMinorTick','on')

%legend('Film','Bulk','location','northwest','Orientation','horizontal')
ylabel('Water content [wt%]','FontSize', 16);
xlabel('Distance from the capillary tip [\mum]','FontSize', 16);

%print(f5,'WaterCalib_RelScal', '-dpdf', '-r300','-bestfit'); 
%print(f5,'WaterIntegral_translated', '-dsvg');

%% 
%{
Max_water = max(Y_fitted_water(60:value-1,1)); 
Min_water = min(Y_fitted_water(60:value-1,1)); 

Max_water_MIX = max(Y_fitted_water_MIX(60:value-1,1));
Min_water_MIX = min(Y_fitted_water_MIX(60:value-1,1));

MAXMINvalues = [Max_water;Min_water;Max_water_MIX;Min_water_MIX];
Min_ABS = min(MAXMINvalues,[],"all");
Max_ABS = max(MAXMINvalues,[],"all");
ABS_abs_avr = (Max_ABS+Min_ABS)./2; 

%Y_fitted_water_REP = Y_fitted_water; 
%Rep_Vect = repmat(ABS_abs_avr,value-60,1); 
%Y_fitted_water_REP(60:value-1,1) = Rep_Vect(:,1);
%}

%% 
for cd = 1:1:size(Y_fitted_water,1);
    AV_VEC_water(cd,1) = (Y_fitted_water(cd,1)+Y_fitted_water_MIX(cd,1))./2; 
end 

Untill = 20;
Smoothing_value = 20;
AV_VEC_water_smooth = smoothdata(AV_VEC_water,1,"movmean",Smoothing_value); 
AV_VEC_water_smooth(Untill:end,1) = AV_VEC_water(Untill:end,1); 

%% Fitted water profile 
f5 = figure;
f5.Position = [300 300 600 300]; 
%ax1 = axes;
plot(X_vector(:,1),Y_fitted_water(:,1),'-','color',[0.000, 0.447, 0.741],'linewidth',1);
hold on 
plot(X_vector,Y_fitted_water_MIX,'-','color',[0.850, 0.325, 0.098],'linewidth',1);
%plot(X_vector,Y_fitted_water_MIX,'o','color',[0.000, 0.447, 0.741],'linewidth',1.5);
%plot(X_vector,Y_fitted_water_REP,'o-','color',[0.5, 0.5, 0.5],'linewidth',1);
plot(X_vector,AV_VEC_water,'-','color',[0, 0, 0],'linewidth',1.5);
%hold on 
%plot(X_vector,AV_VEC_water_smooth,'o-','color',[1, 0, 0],'linewidth',1);

hold off
ylim([0 100]);
yticks(0:20:100);

xlim([15 X_vector(end)]);
%xlim([X_vector(1) X_vector(end)]);
%xlim([60 X_vector(end)]);
xticks(X_ticks);
xticklabels(X_ticks_corr);

grid on
%set(gca,'XMinorTick','on')
%set(gca,'YMinorTick','on')

%legend('Film','Bulk','location','northwest','Orientation','horizontal')
ylabel('Water content [wt%]','FontSize', 16);
xlabel('Distance from the capillary tip [\mum]','FontSize', 16);

%print(f5,sprintf('LipidPhase_1090_Fig2_AVR'), '-dsvg'); 



% Fitted water profile 
f5 = figure;
f5.Position = [300 300 600 300]; 
%ax1 = axes;
plot(X_vector(:,1),Y_fitted_water(:,1),'o','color',[0.000, 0.447, 0.741],'linewidth',1);
hold on 
plot(X_vector,Y_fitted_water_MIX,'o','color',[0.850, 0.325, 0.098],'linewidth',1);
%plot(X_vector,Y_fitted_water_MIX,'o','color',[0.000, 0.447, 0.741],'linewidth',1.5);
%plot(X_vector,Y_fitted_water_REP,'o-','color',[0.5, 0.5, 0.5],'linewidth',1);
%plot(X_vector,AV_VEC_water,'-','color',[0, 0, 0],'linewidth',1.5);
%hold on 
plot(X_vector,AV_VEC_water_smooth,'o-','color',[0, 0, 0],'linewidth',1);

hold off
ylim([0 100]);
yticks(0:20:100);

xlim([15 X_vector(end)]);
%xlim([X_vector(1) X_vector(end)]);
%xlim([60 X_vector(end)]);
xticks(X_ticks);
xticklabels(X_ticks_corr);

grid on
%set(gca,'XMinorTick','on')
%set(gca,'YMinorTick','on')

%legend('Film','Bulk','location','northwest','Orientation','horizontal')
ylabel('Water content [wt%]','FontSize', 16);
xlabel('Distance from the capillary tip [\mum]','FontSize', 16);

%print(f5,sprintf('WATERAavrSmoothed_1090_'), '-dsvg'); 


%% CC fitting from the previouse script - Tenery Diagram
% HERE you need to know what tabel position coresponds to what! 

% In the tables:  
%{ 
Table_coef_all(:,1) = coef1_sumN_film;
Table_coef_all(:,2) = coef2_sumN_film;
Table_coef_all(:,3) = coef1_sumN_bulk;
Table_coef_all(:,4) = coef2_sumN_bulk;
Table_coef_all(:,5) = Av_c1;
Table_coef_all(:,6) = Av_c2;
Table_coef_all(:,7) = 1:1:190;
Table_coef_all(:,8) = coef_CH1;
Table_coef_all(:,9) = coef_CH1_bulk;
Table_coef_all(:,10) = Phase1_AVR;


Table_allcoef_Smooth(:,1) = Av_c1;
Table_allcoef_Smooth(:,2) = Av_c1_smooth;
Table_allcoef_Smooth(:,3) = Av_c1_Rep;
Table_allcoef_Smooth(:,5) = ScanV_corr;
Table_allcoef_Smooth(:,5) = ScanV;
%}

Av_c1 = Table_allcoef1090_New(:,5)./100; 
Av_c2 = Table_allcoef1090_New(:,6)./100; 
Av_PHASE_1 = Table_allcoef1090_New(:,10); 

%Av_c1_Smooth = Table_allcoef_Smooth(:,1)./100;
%Vector_Smooth = Table_allcoef_Smooth(:,2);


%Untill = 115;
%Smoothing_value = 20;

%replacing smooth data
%rep_i = 1; rep_f = 115;
Av_c1_smooth = smoothdata(Av_c1,1,"movmean",Smoothing_value); 
Av_c1_Rep = Av_c1_smooth;  
Av_c1_Rep(Untill:end,1) = Av_c1(Untill:end,1);
Av_c2_Rep = 1 - Av_c1_Rep; 

C1C2_Table = [];
C1C2_Table(:,1) = Av_c1_Rep(:,1);
C1C2_Table(:,2) = Av_c2_Rep(:,1);

%Av_c1_Rep = Av_c1(:,:);
%Av_c2_Rep = 1 - Av_c1_Rep; 


%% SAVE WATER PROFILE 
% PUT in to table

%{
Table_Water(:,1) = Av_c1_Rep;
Table_Water(:,2) = AV_VEC_water_smooth;
Table_Water(:,3) = Data_totIvalue;
Table_Water(:,4) = Data_CNIvalue;
Table_Water(:,5) = Water_Integral_MIX;
Table_Water(:,6) = Y_fitted_water_MIX;
Table_Short_WaterProfile_1090 = Table_Water; 
save('Table_Short_WaterProfile_1090.mat','Table_Short_WaterProfile_1090'); %<----- open
%}

%% Paper FIGURES  - Both Axis 2X  
set(0, 'DefaultAxesFontSize', 18);

X_init = 15;
X_final = 170-20; %Scan_Vector(end); %110; 
xTicks = X_init:20:X_final;
xTicks_lable = xTicks - X_init; 
yTicks = 0:50:100;

%Define relative axis 
film_ini = 15 ; %<--------- Select initial position 
film_end = film_ini+99 ; %<--------- Select final film position
size_ofFilm = film_end - film_ini; 

tickVx2 = linspace(film_ini,film_end,6); %<----- SCALE
tickslable_conv = tickVx2 - film_ini; %DONT USE
tickslable_rel = (tickslable_conv.*100)./tickslable_conv(end); 


fig = figure;
%fig.Position = [200 200 600 350];
fig.Position = [400 400 600 310];
ax1 = axes(fig);

%Left side 
yyaxis right
set(gca,'ycolor','k') ;

p10 = plot(ScanV(1,X_init+4:X_final),Av_PHASE_1(X_init+4:X_final,1),'ko-','linewidth',1.5);

xlim([X_init X_final]);
xticks(xTicks);
xticks(xTicks);

ylim([0 100]);
yticks(yTicks);

%yticklabels({'0:1','2:8','4:6','6:4','8:2','1:0'});
%yticklabels({'0:1','1:1','1:0'});
%yticklabels({'L\beta','1:1','L\alpha'});
set(gca,'TickDir','out');
set(gca,'yticklabel',[]);

ylabel('Phase Ratio','FontSize', 18);
yLabelHandle = get(gca, 'YLabel');
set(yLabelHandle, 'Position', get(yLabelHandle, 'Position') + [1 0 0]);


%Right side 
yyaxis left ;
set(gca,'ycolor','r');

hold on 
%p1 = plot(ScanV,Av_c1(:,1).*100,'-','linewidth',1.5); %,'color',[0.7 0.7 0.7],'linewidth',1.5);
p1 = plot(ScanV(1,X_init+4:X_final),Av_c1_Rep(X_init+4:X_final,1).*100,'ro-','linewidth',1.5); %,'color',[0.7 0.7 0.7],'linewidth',1.5);

hold off ;

grid on;

ax = gca; % Get handle to current axes.
%ax.XColor = 'r'; % Red
%ax.YColor = 'b'; % Blue
ax.GridColor = [0, 0, 0];
ax.MinorGridColor = [0, 0, 0];

set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
%set(gca,'TickDir','out');

xlim([X_init 125]);
xticks(xTicks);
xticklabels(xTicks_lable);

ylim([0 100]);
yticks(0:20:100);

ylabel('DOPC/DOPC+DPPC [wt%]','FontSize', 16);
yLabelHandle = get(gca, 'YLabel');
set(yLabelHandle, 'Position', get(yLabelHandle, 'Position') + [0 0 0]);

xlabel('Distance [\mum]','FontSize', 18);
xLabelHandle = get(gca, 'XLabel');
set(xLabelHandle, 'Position', get(xLabelHandle, 'Position') + [0 0 0]);

%legend([p5 p10],'Smooth','PHASE fit')

ax1.Box = 'off';
xline(ax1,ax1.XLim(2));


%Relative Top AXIS --------
%{
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

set(gca,'clipping','off');
%print(fig,sprintf('LipidPhase_1090_Fig2_REL_A'), '-dpdf', '-r300','-bestfit'); 
%print(fig,sprintf('LipidPhase_1090_Fig2_REL_A'), '-dsvg'); 



% WATER ------------------------------------------------- 

f4 = figure;
%f4.Position = [200 200 600 350]; % for figure 2
%f4.Position = [300 300 500 300]; % for appendix
f4.Position = [400 400 600 310];

ax1 = axes(f4);

plot(ScanV,AV_VEC_water_smooth,'ob-','color',"#3342FF",'linewidth',1.5);
xlim([X_init 125]);
xticks(xTicks);
xticklabels(xTicks_lable);

ylim([0 100]);
yticks(0:20:100);

grid on;
ax = gca; % Get handle to current axes.
%ax.XColor = 'r'; % Red
%ax.YColor = 'b'; % Blue
ax.GridColor = [0, 0, 0];
ax.MinorGridColor = [0, 0, 0];

set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
set(gca,'TickDir','out');

%title('Water profile')
ylabel('Water content [wt%]','FontSize', 18);
xlabel('Distance [\mum]','FontSize', 18);

set(gca,'ycolor','b');

ax1.Box = 'off';
xline(ax1,ax1.XLim(2))

%Relative Top AXIS --------
%{
ax2 = axes('Position', get(ax1, 'Position'),'Color', 'none');
set(ax2, 'XAxisLocation', 'top','YAxisLocation','left');
set(ax2,'Xtick',[],'Ytick',[]) % To hide ax2 ticks
set(gca,'XMinorTick','on')
set(gca,'TickDir','out');
ax2.YAxis.Visible = 'off';

step_rel = tickVx2(2) - tickVx2(1);
ax2.XAxis.MinorTickValues = tickVx2(1):step_rel/2:tickVx2(end);

set(ax2,'XLim',[[X_init X_final]]...
    ,'xtick',tickVx2,'xticklabels',string(tickslable_rel));

%}

Size=[0.125  0.185  0.76  0.70]; %[0.1300 0.1100 0.7750 0.8150] (default, ) 
set(ax1,'InnerPosition',Size);
%set(ax2,'InnerPosition',Size);

%print(f4,sprintf('Water_1090_fig_REL_A'), '-dpdf', '-r300','-bestfit'); 
%print(f4,sprintf('Water_1090_fig_REL_A'), '-dsvg'); 

%%
%{
Table_2080_20(:,1) = ScanV;
Table_2080_20(:,2) = AV_VEC_water_smooth;
Table_2080_20(:,3) = Av_c1_Rep;
Table_2080_20(:,4) = Av_c2_Rep;
Table_2080_20(:,5) = Av_PHASE_1; 

Table_2080_short_20 = Table_2080_20; 
save('Table_2080_short_20.mat','Table_2080_short_20'); %<----- open
writematrix(Table_2080_short_20, 'Table_2080_short_20.xlsx');
%}

%% TERNARY PHASE DIAGRAM

% Graph Parameters
W(:,1) = AV_VEC_water_smooth(:,1)./100; 
c1_water = W(:,1); % water content
G(:,1) = 1 - W(:,1); %lipid fraction.

B_AV = [];
%%B_AV(:,1) = G(:,1).* Av_c1(:,1);
%%B_AV(:,2) = G(:,1).* Av_c2(:,1);
B_AV(:,1) = G(:,1).* Av_c1_Rep(:,1);
B_AV(:,2) = G(:,1).* Av_c2_Rep(:,1);
%B_AV(:,1) = G(:,1).* SCAN_smooth_CC(:,1);
%B_AV(:,2) = G(:,1).* SCAN_smooth_CC_2(:,1);
B_AV(:,3) = W(:,1);
B_AV(:,4) = W(:,1) + B_AV(:,1) + B_AV(:,2);

B_AV(:,5) = ScanV; 

B_AV_2080_short = B_AV; 
%save('TrajectoryPath_2080_short.mat','B_AV_2080_short'); %<----- open
%writematrix(B_AV_2080_short, 'TrajectoryPath_2080_short.xlsx');


c2_AV = B_AV(:,1);
Xt_AV = 0.5-c1_water.*cos(pi/3)+c2_AV./2;
Yt_AV = 0.866-c1_water.*sin(pi/3)-c2_AV.*cot(pi/6)./2;

%% Plot CC vs watergradient 
%AV_VEC_water_smooth_c = AV_VEC_water_smooth./100;
%Theo_left = 1 - AV_VEC_water_smooth_c; 
%Theo_left = Theo_left;
wtA = 0.8;
wtB = 1-wtA; 

vector_wtA = G(:,1).*0 + wtA;  
vector_wtB = 1 -vector_wtA;

Theo_c1 = G(:,1).*wtA; 
Theo_c2 = G(:,1).*wtB; 
Sum = W(:,1) + Theo_c1 + Theo_c2; 

Theo_c1_wt = Theo_c1./(Theo_c1+Theo_c2);
Theo_c2_wt = Theo_c2./(Theo_c2+Theo_c1);

%AV_VEC_water_smooth
Av_c2_Rep = 1 - Av_c1_Rep; 

B_AV1_wt = B_AV(:,1)./(B_AV(:,1)+B_AV(:,2));
B_AV2_wt = 1- B_AV1_wt; 

Values2080 = []; 
Values2080(:,1) = W(:,1);
Values2080(:,2) = B_AV(:,1); 
Values2080(:,3) =  B_AV(:,2); 
Values2080(:,4) = Theo_c1(:,1); 
Values2080(:,5) =  Theo_c2(:,1); 

figure
plot(AV_VEC_water_smooth(15:150,1),B_AV(15:150,1).*100,'o-','linewidth',1.5);
%plot(AV_VEC_water_smooth,B_AV1_wt.*100,'o-','linewidth',1.5);
hold on 
plot(AV_VEC_water_smooth(15:150,1),B_AV(15:150,2).*100,'o-','linewidth',1.5);
%plot(AV_VEC_water_smooth,B_AV2_wt.*100,'o-','linewidth',1.5);
plot(AV_VEC_water_smooth,Theo_c1.*100 ,'ro-','linewidth',1.5);
plot(AV_VEC_water_smooth,Theo_c2.*100,'ko-','linewidth',1.5);
hold off
legend('DOPC','DPPC');
xlabel('Water gradinet [wt%]','FontSize', 18);
%xlabel('Distance [\mum]','FontSize', 18);

%{

figure
plot(AV_VEC_water_smooth(15:170,1),Av_c1_Rep(15:170,1).*100,'ro-','linewidth',1.5);
hold on 
plot(AV_VEC_water_smooth(15:170,1),vector_wtB(15:170,1).*100 ,'o-','color',[0.6 0.6 0.6],'linewidth',1.5);
hold off
legend('DOPC','DPPC');
xlabel('Water gradinet [wt%]','FontSize', 18);
%xlabel('Distance [\mum]','FontSize', 18);

%}
%%

X_init_w = 13;
X_final_w = 170;


fig = figure;
%fig.Position = [200 200 600 350];
fig.Position = [400 400 620 320];
ax1 = axes(fig);

%Left side 
%yyaxis right
%set(gca,'ycolor','k') ;

%p10 = plot(AV_VEC_water_smooth(X_init_w:X_final_w,1),Av_PHASE_1(X_init_w:X_final_w,1),'ko-','linewidth',1.5);

%xlim([X_init X_final]);
%xticks(xTicks);
%xticks(xTicks);

%ylim([0 100]);
%yticks(yTicks);

%yticklabels({'0:1','2:8','4:6','6:4','8:2','1:0'});
%yticklabels({'0:1','1:1','1:0'});
%yticklabels({'L\beta','1:1','L\alpha'});
%set(gca,'TickDir','out');
%set(gca,'yticklabel',[]);

%ylabel('Phase Ratio','FontSize', 18);
%yLabelHandle = get(gca, 'YLabel');
%set(yLabelHandle, 'Position', get(yLabelHandle, 'Position') + [1 0 0]);

%Right side 
%yyaxis left ;
set(gca,'ycolor','r');

hold on 
p3 = plot(AV_VEC_water_smooth(13:170,1),vector_wtB(13:170,1).*100 ,'o-','color',[0.5 0.5 0.5],'linewidth',1.5,'markersize',8);
%p3 = plot(AV_VEC_water_smooth(15:5:170,1),vector_wtB(15:5:170,1).*100 ,'--','color',[0.6 0.6 0.6],'linewidth',2.5);
p1 = plot(AV_VEC_water_smooth(X_init_w:145,1),Av_c1_Rep(X_init_w:145,1).*100,'ro-','linewidth',1.5,'markersize',8); %,'color',[0.7 0.7 0.7],'linewidth',1.5);

hold off ;

grid on;

ax = gca; % Get handle to current axes.
%ax.XColor = 'r'; % Red
%ax.YColor = 'b'; % Blue
ax.GridColor = [0, 0, 0];
ax.MinorGridColor = [0, 0, 0];

set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
%set(gca,'TickDir','out');

%xlim([X_init 125]);
%xticks(xTicks);
%xticklabels(xTicks_lable);

ylim([0 100]);
yticks(0:20:100);

xlim([10 56]);

set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
set(gca,'TickDir','out');

ylabel('DOPC/DOPC+DPPC [wt%]','FontSize', 16);
yLabelHandle = get(gca, 'YLabel');
set(yLabelHandle, 'Position', get(yLabelHandle, 'Position') + [0 0 0]);

%set(gca,'xcolor','b') ;
legend([p1, p3],'Serpentine trajectory','Linear trajectory','location','northwest');
%legend([p1, p3],'','','location','northwest');

xlabel('Water content [wt%]','FontSize', 18);
xLabelHandle = get(gca, 'XLabel');
set(xLabelHandle, 'Position', get(xLabelHandle, 'Position') + [0 0 0]);

Size=[0.125  0.185  0.76  0.70]; %[0.1300 0.1100 0.7750 0.8150] (default, ) 
set(ax1,'InnerPosition',Size);
%set(ax2,'InnerPosition',Size);

set(gca,'clipping','off');
%print(fig,sprintf('LipidPhase_1090_Fig2_REL_A'), '-dpdf', '-r300','-bestfit'); 
%print(fig,sprintf('LipidPhase_2080_Fig2_'), '-dsvg'); 


%% 
%{
Table_TrajectoryLines(:,1) = AV_VEC_water_smooth;
Table_TrajectoryLines(:,2) = vector_wtB; 
Table_TrajectoryLines(:,3) = Av_c1_Rep; 

save('Table_TrajectoryLines_2080.mat','Table_TrajectoryLines'); %<----- open
writematrix(Table_TrajectoryLines, 'Table_TrajectoryLines_2080.xlsx');
%}

%% ----------------------

%{

fig = figure;
%fig.Position = [200 200 600 350];
fig.Position = [400 400 600 310];
ax1 = axes(fig);


%set(gca,'TickDir','out');
%set(gca,'yticklabel',[]);

%set(gca,'ycolor','r');

hold on 
pA = plot(AV_VEC_water_smooth(13:170,1),vector_wtB(13:170,1).*100 ,'o-','color',[0.6 0.6 0.6],'linewidth',1.5);
%p2 = plot(AV_VEC_water_smooth(88:94,1),Av_c1_Rep(88:94,1).*100,'o-','color',[0 0.3 0.6],'linewidth',1.5); %,'color',[0.7 0.7 0.7],'linewidth',1.5);

%{
p1 = plot(AV_VEC_water_smooth(X_init_w:88,1),Av_c1_Rep(X_init_w:88,1).*100,'^-','color',[0 0.6 0.2],'linewidth',1.5); %,'color',[0.7 0.7 0.7],'linewidth',1.5);
p2 = plot(AV_VEC_water_smooth(88:93,1),Av_c1_Rep(88:93,1).*100,'o-','color',[0 0.6 0.2],'linewidth',1.5); %,'color',[0.7 0.7 0.7],'linewidth',1.5);
p3 = plot(AV_VEC_water_smooth(93:107,1),Av_c1_Rep(93:107,1).*100,'*-','color',[0 0.6 0.2],'linewidth',1.5); %,'color',[0.7 0.7 0.7],'linewidth',1.5);
p4 = plot(AV_VEC_water_smooth(108:145,1),Av_c1_Rep(108:145,1).*100,'o-','color',[0 0.6 0.2],'linewidth',1.5); %,'color',[0.7 0.7 0.7],'linewidth',1.5);
%}
p1 = plot(AV_VEC_water_smooth(X_init_w:88,1),Av_c1_Rep(X_init_w:88,1).*100,'k^-','linewidth',1.5); %,'color',[0.7 0.7 0.7],'linewidth',1.5);
p2 = plot(AV_VEC_water_smooth(88:93,1),Av_c1_Rep(88:93,1).*100,'ko-','linewidth',1.5); %,'color',[0.7 0.7 0.7],'linewidth',1.5);
p3 = plot(AV_VEC_water_smooth(93:107,1),Av_c1_Rep(93:107,1).*100,'k*-','linewidth',1.5); %,'color',[0.7 0.7 0.7],'linewidth',1.5);
p4 = plot(AV_VEC_water_smooth(108:133,1),Av_c1_Rep(108:133,1).*100,'ko-','linewidth',1.5); %,'color',[0.7 0.7 0.7],'linewidth',1.5);

%p1 = plot(AV_VEC_water_smooth(X_init_w:88,1),Av_c1_Rep(X_init_w:88,1).*100,'o-','color',[0 0.6 0.2],'linewidth',1.5); %,'color',[0.7 0.7 0.7],'linewidth',1.5);
%p2 = plot(AV_VEC_water_smooth(88:94,1),Av_c1_Rep(88:94,1).*100,'o-','color',[0 0.3 0.6],'linewidth',1.5); %,'color',[0.7 0.7 0.7],'linewidth',1.5);
%p3 = plot(AV_VEC_water_smooth(94:108,1),Av_c1_Rep(94:108,1).*100,'ro-','linewidth',1.5); %,'color',[0.7 0.7 0.7],'linewidth',1.5);
%p4 = plot(AV_VEC_water_smooth(108:135,1),Av_c1_Rep(108:135,1).*100,'o-','color',[0 0.3 0.6],'linewidth',1.5); %,'color',[0.7 0.7 0.7],'linewidth',1.5);

hold off ;

grid on;

ax = gca; % Get handle to current axes.
%ax.XColor = 'r'; % Red
%ax.YColor = 'b'; % Blue
ax.GridColor = [0, 0, 0];
ax.MinorGridColor = [0, 0, 0];

set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
%set(gca,'TickDir','out');

%xlim([X_init 125]);
%xticks(xTicks);
%xticklabels(xTicks_lable);

ylim([0 100]);
yticks(0:20:100);

xlim([10 56]);

ylabel('DOPC/DOPC+DPPC [wt%]','FontSize', 16);
yLabelHandle = get(gca, 'YLabel');
set(yLabelHandle, 'Position', get(yLabelHandle, 'Position') + [0 0 0]);

%set(gca,'xcolor','b') ;

xlabel('Water content [wt%]','FontSize', 18);
xLabelHandle = get(gca, 'XLabel');
set(xLabelHandle, 'Position', get(xLabelHandle, 'Position') + [0 0 0]);

Size=[0.125  0.185  0.76  0.70]; %[0.1300 0.1100 0.7750 0.8150] (default, ) 
set(ax1,'InnerPosition',Size);
%set(ax2,'InnerPosition',Size);

set(gca,'clipping','off');
%print(fig,sprintf('LipidPhase_1090_Fig2_REL_A'), '-dpdf', '-r300','-bestfit'); 
%print(fig,sprintf('LipidPhase_2080_Fig2_nolable'), '-dsvg'); 

%}
%% 

%EnergyLine = []; 

%EnergyLine
%save('EnergyLine_2080.mat','EnergyLine'); %<----- open

load('EnergyLine_2080.mat'); 

%%

Gtot = sum(EnergyLine(:,2),1); 

fig = figure;
%fig.Position = [200 200 600 350];
fig.Position = [400 400 620 320];

ax1 = axes(fig);

hold on 
p2 = plot([10 56],[0 0],'--','color',[0.5 0.5 0.5],'linewidth',3); %,'color',[0.7 0.7 0.7],'linewidth',1.5);

p2 = plot(EnergyLine(:,1).*100,EnergyLine(:,2),'o-','color',[0 0.5 0.1],'linewidth',1.5,'markersize',8); %,'color',[0.7 0.7 0.7],'linewidth',1.5);

hold off ;

grid on;

ax = gca; % Get handle to current axes.
%ax.XColor = 'r'; % Red
%ax.YColor = 'b'; % Blue
ax.GridColor = [0, 0, 0];
ax.MinorGridColor = [0, 0, 0];

set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
set(gca,'TickDir','out');

%xlim([X_init 125]);
%xticks(xTicks);
%xticklabels(xTicks_lable);

ylim([-15*10^-3 5*10^-3]);
%yticks(0:20:100);

xlim([10 56]);
%legend([p2],sprintf('SUM %.3g',Gtot),'location','southeast','FontSize', 16);
%legend([p2],'','location','southeast','FontSize', 16);

ylabel('ΔG [J/m^2]','FontSize', 16);
yLabelHandle = get(gca, 'YLabel');
set(yLabelHandle, 'Position', get(yLabelHandle, 'Position') + [0 0 0]);

%set(gca,'xcolor','b') ;

xlabel('Water content [wt%]','FontSize', 18);
xLabelHandle = get(gca, 'XLabel');
set(xLabelHandle, 'Position', get(xLabelHandle, 'Position') + [0 0 0]);

Size=[0.125  0.185  0.76  0.70]; %[0.1300 0.1100 0.7750 0.8150] (default, ) 
set(ax1,'InnerPosition',Size);
%set(ax2,'InnerPosition',Size);

set(gca,'clipping','off');
%print(fig,sprintf('LipidPhase_1090_Fig2_REL_A'), '-dpdf', '-r300','-bestfit'); 
%print(fig,sprintf('LipidPhase_2080_Fig2_ENERGY_noLabel'), '-dsvg'); 

%% 
%{
Table_FreeEnergyLine(:,1) = EnergyLine(:,1).*100;
Table_FreeEnergyLine(:,2) = EnergyLine(:,2); 

save('Table_FreeEnergyLine_2080.mat','Table_FreeEnergyLine'); %<----- open
writematrix(Table_FreeEnergyLine, 'Table_FreeEnergyLine_2080.xlsx');
%}

%% Repulsion and attracion terms 

fig = figure;
%fig.Position = [200 200 600 350];
fig.Position = [400 400 600 310];
ax1 = axes(fig);

hold on 
%p2 = plot([10 56],[0 0],'k-','linewidth',1.5); %,'color',[0.7 0.7 0.7],'linewidth',1.5);

p2 = plot(EnergyLine(:,1).*100,EnergyLine(:,3),'*-','linewidth',1.5); %,'color',[0.7 0.7 0.7],'linewidth',1.5);
p2 = plot(EnergyLine(:,1).*100,EnergyLine(:,4).*-1,'*-','linewidth',1.5); %,'color',[0.7 0.7 0.7],'linewidth',1.5);
p2 = plot(EnergyLine(:,1).*100,EnergyLine(:,5),'o-','linewidth',1.5); %,'color',[0.7 0.7 0.7],'linewidth',1.5);
p2 = plot(EnergyLine(:,1).*100,EnergyLine(:,6).*-1,'o-','linewidth',1.5); %,'color',[0.7 0.7 0.7],'linewidth',1.5);
hold off ;

grid on;

ax = gca; % Get handle to current axes.
%ax.XColor = 'r'; % Red
%ax.YColor = 'b'; % Blue
ax.GridColor = [0, 0, 0];
ax.MinorGridColor = [0, 0, 0];

set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
set(gca,'TickDir','out');

%xlim([X_init 125]);
%xticks(xTicks);
%xticklabels(xTicks_lable);

%ylim([-15*10^-3 5*10^-3]);
%yticks(0:20:100);

xlim([10 56]);
%legend([p2],sprintf('SUM %.3g',Gtot),'location','southeast','FontSize', 16);
%legend([p2],'','location','southeast','FontSize', 16);
legend

%ylabel('ΔG','FontSize', 16);
yLabelHandle = get(gca, 'YLabel');
set(yLabelHandle, 'Position', get(yLabelHandle, 'Position') + [0 0 0]);

%set(gca,'xcolor','b') ;

xlabel('Water content [wt%]','FontSize', 18);
xLabelHandle = get(gca, 'XLabel');
set(xLabelHandle, 'Position', get(xLabelHandle, 'Position') + [0 0 0]);

Size=[0.125  0.185  0.76  0.70]; %[0.1300 0.1100 0.7750 0.8150] (default, ) 
set(ax1,'InnerPosition',Size);
%set(ax2,'InnerPosition',Size);

set(gca,'clipping','off');
%print(fig,sprintf('LipidPhase_1090_Fig2_REL_A'), '-dpdf', '-r300','-bestfit'); 
%print(fig,sprintf('LipidPhase_2080_Fig2_ENERGY_'), '-dsvg'); 


%% Color line

% Select SCANS in the path of phase diagram 
Cbarline_st = 13; Cbarline_fi = 159; %SpecSize;
%Cbarline_st = 120; Cbarline_fi = 135; %SpecSize;

%###################
%film_st = 95; film_fi = 111; %SpecSize;
film_st = 13; film_fi = 133; %SpecSize;

%###################

linelenght = film_fi - film_st +1;
line_cVector = 1:1:linelenght;
line_cNr = film_st:1:film_fi; 
%{
LineCmap = jet(linelenght).*[0.95 0.95 0.95]; %color maps
%LineCmap = turbo(linelenght); 
LineCmap = flip(LineCmap); 
%}

%
Barlenght = Cbarline_fi - Cbarline_st +1;
LineCmap = jet(Barlenght).*[0.95 0.95 0.95]; %color maps
LineCmap = flip(LineCmap); %to reverse color order of the bar.
Cc_init = film_st-Cbarline_st+1;
Cc_finit = Cc_init + linelenght -1; 
%}

%Set the colorbar 
Bar_size = linelenght/10; % Divident for the colorbar
ticks_scaleBar = linspace(0,1,Bar_size+1); % +1 
ticks_lableBar = linspace(film_st-film_st,film_fi-film_st,size(ticks_scaleBar,2))'; % if
%you want to set film_st as 0 point. 
%ticks_lableBar = linspace(film_st,film_fi,size(ticks_scaleBar,2))';

%cut data to put in T-diagram
Xt_line = Xt_AV(film_st:film_fi,1);
Yt_line = Yt_AV(film_st:film_fi,1);
Yt_line(end,1) = NaN; % so the Patch plot dont make a loop.
C_line = line_cVector; %create color vector 

fig2 = figure; %-- Plot the tenery system dotted line
[h,hg,htick]=terplot;
hold on 

%line_start = 20; line_end = 80;
%C_line1 = [1:1:51];
%p = patch(Xt_line(line_start:line_end,1),Yt_line(line_start:line_end,1),C_line(1,line_start:line_end),'FaceColor','None','EdgeColor','interp','Marker','o','Markersize',10,'linewidth',2);
p = patch(Xt_line(),Yt_line(),C_line(),'FaceColor','None','EdgeColor','interp','Marker','o','Markersize',10,'linewidth',2);

Xt_line = double(Xt_line);
Yt_line = double(Yt_line);
%text(Xt_line+0.015,Yt_line, num2str(line_cNr'),'FontWeight','bold','Fontsize',12);

p.MarkerFaceColor = 'flat';
p.MarkerEdgeColor = 'k'; %'flat';
hold off ;
%colormap(LineCmap());
colormap(LineCmap(Cc_init:Cc_finit,:));

%print(fig2,sprintf('T_diagram_1090_Cdot_SMooth_WATER'), '-dpdf', '-r300','-bestfit'); 
%print(fig2,sprintf('T_Cdot_60170_smoothed'), '-dsvg'); 
%print(fig2,sprintf('T_diagram_10120'), '-dsvg'); 

%{
fig2 = figure; %-- Plot the tenery system line 
[h,hg,htick]=terplot;
hold on ;

p = patch(Xt_line,Yt_line,C_line,'FaceColor','None','EdgeColor','flat','linewidth',6);

ax = gca; 
ax.Position = ax.Position - [0.07 0 0 0];
hold off ;
%colormap(LineCmap(1:film_fi-film_st,:));
colormap(LineCmap(Cc_init:Cc_finit,:));

cb = colorbar;
pos = get(cb,'Position');
set(cb,'Position',pos + [0.13,0.05,0,0]);
%}
%%

FR = figure;
FR.Position = [300 300 600 340]; 
plot(x_cm1, Data_scan_MIX(:,13:146),'-','linewidth',1.5);
%plot(x_cm1, Data_scan_MIX(:,180:-1:60),'-','linewidth',1.5);
%set(gca, 'ColorOrder', jet);
%set(gca, 'ColorOrder', cool(size(Data_scan_MIX(:,15:120),2)));
set(gca, 'ColorOrder', LineCmap(Cc_init:Cc_finit,:));

%set(gca, 'ColorOrder', jet);
grid on 
xlim([x_cm1(1) x_cm1(end)]);
%ylim([x_cm1(1) x_cm1(end)]);

c = colorbar;
colormap(LineCmap(Cc_init:Cc_finit,:));
ticks_lableBar = 10:10:120; %size(Data_scan_BG,2);
ticks_lableBar = ticks_lableBar';
ticks_scaleBar = linspace(0,1,size(ticks_lableBar,1));
set(c, 'Ticks', ticks_scaleBar);
set(c, 'TickLabels', num2str(ticks_lableBar));
set(c,'direction','reverse');

%print(FR,'Raman_2080', '-dsvg');

%%

CbarFig = figure; %Color bar figure 
set(0, 'DefaultAxesFontSize', 18);
ax = axes;
colormap(LineCmap);
c = colorbar(ax);
ax.Visible = 'off';
%colorbar('direction','reverse')

%ticks_scaleBar = linspace(0,1,17);
%ticks_lableBar = linspace(1,17,17)';
set(c, 'Ticks', ticks_scaleBar,'TickDir','both');
set(c, 'TickLabels', num2str(ticks_lableBar));
%ylabel(c, 'Distance from the capillary tip [\mum]','Fontsize',18);
set(c,'direction','reverse');
%print(CbarFig,sprintf('Cbar_1090_0_130_fliped'), '-dpdf', '-r300','-bestfit'); 
%print(CbarFig,sprintf('Cbar_2080_15135'), '-dsvg'); 
%print(CbarFig,sprintf('Cbar_2080_15125'), '-dsvg'); 


%%
CbarFig = figure; %Color bar figure 
set(0, 'DefaultAxesFontSize', 18);
ax = axes;
colormap(LineCmap(Cc_init:Cc_finit,:));
c = colorbar(ax);
ax.Visible = 'off';
%colorbar('direction','reverse')

%Cc_init = film_st-Cbarline_st+1;
%Cc_finit = Cc_init + linelenght -1;

ticks_lableBar = linspace(Cc_init-1,Cc_finit-1,size(ticks_scaleBar,2))'; % if

set(c, 'Ticks', ticks_scaleBar,'TickDir','both');
set(c, 'TickLabels', num2str(ticks_lableBar));
ylabel(c, 'Distance from the capillary tip [\mum]','Fontsize',18);
set(c,'direction','reverse');
%print(CbarFig,sprintf('Cbar_1090_0_110_fCUT'), '-dpdf', '-r300','-bestfit'); 
%print(CbarFig,sprintf('Cbar_1090_60180'), '-dsvg'); 

%}
%%
%{

%################### Plot SCANS in the path
film_st = 60; %film_fi = end_spec;
%###################

fig2 = figure; %-- Plot the tenery system
[h,hg,htick]=terplot;
hold on 

%ErrorLines
%pt2 = plot(Xt_film(i_BulkR(1):film_fi,1),Yt_film(i_BulkR(1):film_fi,1),'k-.','linewidth',1);
%pt4 = plot(Xt_bulk(i_BulkR(1):film_fi,1),Yt_bulk(i_BulkR(1):film_fi,1),'k-.','linewidth',1);

%pt1 = plot(Xt_film(film_st:i_BulkR(1),1),Yt_film(film_st:i_BulkR(1),1),'-.','color',[0 0.5 0],'linewidth',1);
%pt3 = plot(Xt_bulk(film_st:i_BulkR(1),1),Yt_bulk(film_st:i_BulkR(1),1),'-.','color',[0 0.5 0],'linewidth',1);

%Mainline
%pt5 = plot(Xt_AV(film_st:i_BulkR(1),1),Yt_AV(film_st:i_BulkR(1),1),'ro-','linewidth',1.5);
%pt6 = plot(Xt_AV(i_BulkR(1)-1:film_fi,1),Yt_AV(i_BulkR(1)-1:film_fi,1),'o-','color',[0.5 0.5 0.5],'linewidth',2);
%pt5 = plot(Xt_AV(film_st:i_BulkR(1)-1,1),Yt_AV(film_st:i_BulkR(1)-1,1),'o-','color',[0 0.5 0],'linewidth',2);

%pt6 = plot(Xt_AV(film_st:end,1),Yt_AV(film_st:end,1),'o-','color',[0 0.5 0],'linewidth',2);
pt6 = plot(Xt_AV(film_st:end,1),Yt_AV(film_st:end,1),'o-','color',[0.2 0.6 0.2],'linewidth',2);

hold off 
%hlabels=terlabel('H2O [wt%]','DOPC [?%]','DPPC [?%]');
legend([pt6],'from 0 to 190')

% Set the renderer to 'painters'
%set(gcf, 'Renderer', 'painters');

% Set the Alpha property to 'color'
%set(gcf, 'Color', 'none');
 
%print(fig2,sprintf('T_diagram_1090'), '-dpng'); 
%print(fig2,sprintf('T_diagram_1090_negative'), '-dpdf', '-r300','-bestfit'); 
%print(fig2,sprintf('T_diagram_1090_negative'), '-dsvg'); 
%}

%% 

Av_c1_line = Av_c1_Rep.*100; 
%Av_c1_short_rep(Inital_short:end,1)

%X_initial_film = 60; %ScanV_corr is shifted by 69 bit not the data points.
X_initial_film = 13;

region_Ai = X_initial_film; 
region_Af = 170; 
Integral_DOPCA = trapz(Av_c1_line(region_Ai:region_Af,1));
Integral_DOPCA_mean = mean(Av_c1_line(region_Ai:region_Af,1));

bas_fini = 30; %Avreage how far from the back! 
Baseline = mean(Av_c1_line(end-bas_fini:end,1));

Integral_Baseline = Baseline.* (region_Af-region_Ai); 
Integral_Total = 100.* (region_Af-region_Ai); 

IntRatio_10wt = Integral_Baseline.*100./Integral_Total;
IntRatio_fitted = Integral_DOPCA.*100./Integral_Total;

fig = figure;
fig.Position = [200 200 600 350]; % for figure 2
hold on ;
%p1 = plot(ScanV, Av_c1_line,'-','Color', [0.5 0.5 0.5], 'LineWidth', 1);
p2 = plot(ScanV(1,region_Ai:region_Af), Av_c1_line(region_Ai:region_Af,1),'ro-','linewidth',2);
p3 = plot([ScanV(1,end-bas_fini) ScanV(1,end)],[Baseline Baseline],'-','color',[0.5 0.5 0.5],'linewidth',3);
p3 = plot([ScanV(1,1) ScanV(1,end)],[Baseline Baseline],'--','color',[0.5 0.5 0.5],'linewidth',3);

hold off ;
legend([p3,p2],sprintf('Baseline %.3g wt%%',Baseline),sprintf('Gradient %.3g wt%%%',IntRatio_fitted),'location','northwest','FontSize', 16);
grid on ;
ylabel('DOPC/(DOPC+DPPC) [wt%]','FontSize', 16);
xlabel('[μm]');
ylim([0 100]);
xlim([region_Ai 140]);

%print(fig,'Integral_2080', '-dsvg');


%{
X_initial_film = 15;
Av_c1_RepPC = Av_c1_Rep.*100; 

region_Ai = X_initial_film; 
region_Af = 170; %ScanV_short(end); 
Integral_DOPCA = trapz(Av_c1_RepPC(region_Ai:region_Af,1));

bas_fini = 20; %Avreage how far from the back! 
Baseline = mean(Av_c1_RepPC(end-bas_fini:end,1));

Integral_Baseline = Baseline.* (region_Af-region_Ai); 
Integral_Total = 100.* (region_Af-region_Ai); 

IntRatio_10wt = Integral_Baseline.*100./Integral_Total;
IntRatio_fitted = Integral_DOPCA.*100./Integral_Total;

fig = figure;
hold on ;
p1 =plot(ScanV, Av_c1_RepPC,'-','Color', [0.5 0.5 0.5], 'LineWidth', 1);
p2 = plot(ScanV(1,region_Ai:region_Af), Av_c1_RepPC(region_Ai:region_Af,1),'ro-','linewidth',2);
p3 = plot([ScanV(end-bas_fini) ScanV(end)],[Baseline Baseline],'b-','linewidth',3);
hold off ;
legend([p3,p2],sprintf('Baseline %.3g wt%%',Baseline),sprintf('Integral %.3g %%',IntRatio_fitted),'location','northwest','FontSize', 16);
grid on ;
ylabel('DOPC/(DOPC+DPPC) [wt%]','FontSize', 16);
xlabel('Spec nr.');
ylim([0 100]);
xlim([region_Ai region_Af-30]);

%print(fig,'Integral_2080', '-dsvg');

ScanV1 = ScanV';
%}


%% THESIS FIG

%% Paper FIGURES  - Both Axis 2X  
set(0, 'DefaultAxesFontSize', 18);

X_init = 15;
X_final = 170-20; %Scan_Vector(end); %110; 
xTicks = X_init:20:X_final;
xTicks_lable = xTicks - X_init; 
yTicks = 0:50:100;

%Define relative axis 
film_ini = 15 ; %<--------- Select initial position 
film_end = film_ini+99 ; %<--------- Select final film position
size_ofFilm = film_end - film_ini; 

tickVx2 = linspace(film_ini,film_end,6); %<----- SCALE
tickslable_conv = tickVx2 - film_ini; %DONT USE
tickslable_rel = (tickslable_conv.*100)./tickslable_conv(end); 


fig = figure;
%fig.Position = [200 200 600 350];
fig.Position = [400 400 600 310];
ax1 = axes(fig);


%Right side 
%yyaxis left ;
set(gca,'ycolor','r');

hold on 
%p1 = plot(ScanV,Av_c1(:,1).*100,'-','linewidth',1.5); %,'color',[0.7 0.7 0.7],'linewidth',1.5);
p1 = plot(ScanV(1,X_init+4:X_final),Av_c1_Rep(X_init+4:X_final,1).*100,'ro-','linewidth',1.5); %,'color',[0.7 0.7 0.7],'linewidth',1.5);

hold off ;

grid on;

ax = gca; % Get handle to current axes.
%ax.XColor = 'r'; % Red
%ax.YColor = 'b'; % Blue
ax.GridColor = [0, 0, 0];
ax.MinorGridColor = [0, 0, 0];

set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
%set(gca,'TickDir','out');

xlim([X_init 125]);
xticks(xTicks);
xticklabels(xTicks_lable);

ylim([0 100]);
yticks(0:20:100);

ylabel('DOPC/DOPC+DPPC [wt%]','FontSize', 16);
yLabelHandle = get(gca, 'YLabel');
set(yLabelHandle, 'Position', get(yLabelHandle, 'Position') + [0 0 0]);

xlabel('Distance [\mum]','FontSize', 18);
xLabelHandle = get(gca, 'XLabel');
set(xLabelHandle, 'Position', get(xLabelHandle, 'Position') + [0 0 0]);

%legend([p5 p10],'Smooth','PHASE fit')

ax1.Box = 'off';
xline(ax1,ax1.XLim(2));


%Relative Top AXIS --------
%{
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

set(gca,'clipping','off');
%print(fig,sprintf('LipidPhase_1090_Fig2_REL_A'), '-dpdf', '-r300','-bestfit'); 
%print(fig,sprintf('THESIS_LipidPhase_1090_Fig2_REL_A'), '-dsvg'); 



% WATER ------------------------------------------------- 

f4 = figure;
%f4.Position = [200 200 600 350]; % for figure 2
%f4.Position = [300 300 500 300]; % for appendix
f4.Position = [400 400 600 310];

ax1 = axes(f4);

plot(ScanV,AV_VEC_water_smooth,'ob-','color',"#3342FF",'linewidth',1.5);
xlim([X_init 125]);
xticks(xTicks);
xticklabels(xTicks_lable);

ylim([0 100]);
yticks(0:20:100);

grid on;
ax = gca; % Get handle to current axes.
%ax.XColor = 'r'; % Red
%ax.YColor = 'b'; % Blue
ax.GridColor = [0, 0, 0];
ax.MinorGridColor = [0, 0, 0];

set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
set(gca,'TickDir','out');

%title('Water profile')
ylabel('Water content [wt%]','FontSize', 18);
xlabel('Distance [\mum]','FontSize', 18);

set(gca,'ycolor','b');

ax1.Box = 'off';
xline(ax1,ax1.XLim(2))

%Relative Top AXIS --------
%{
ax2 = axes('Position', get(ax1, 'Position'),'Color', 'none');
set(ax2, 'XAxisLocation', 'top','YAxisLocation','left');
set(ax2,'Xtick',[],'Ytick',[]) % To hide ax2 ticks
set(gca,'XMinorTick','on')
set(gca,'TickDir','out');
ax2.YAxis.Visible = 'off';

step_rel = tickVx2(2) - tickVx2(1);
ax2.XAxis.MinorTickValues = tickVx2(1):step_rel/2:tickVx2(end);

set(ax2,'XLim',[[X_init X_final]]...
    ,'xtick',tickVx2,'xticklabels',string(tickslable_rel));

%}

Size=[0.125  0.185  0.76  0.70]; %[0.1300 0.1100 0.7750 0.8150] (default, ) 
set(ax1,'InnerPosition',Size);
%set(ax2,'InnerPosition',Size);

%print(f4,sprintf('Water_1090_fig_REL_A'), '-dpdf', '-r300','-bestfit'); 
%print(f4,sprintf('THESIS_Water_1090_fig_REL_A'), '-dsvg'); 

