clear all
close all
clc

set(0, 'DefaultAxesFontSize', 14)

%%
scan_file = 00347; 
n_acquisitions = 6; %#<----------------- CHANGE HERE 

%path = sprintf('/Users/labecka/Documents/BeamData/PSI 2023/cSAXS/analysis/eiger9m_integ/S%05d_eiger_integ.h5', scan_file);
path = sprintf('DOPC_S00347_eiger_integ.h5');
h5disp(path); I = h5read(path,'/I_all'); q = h5read(path,'/q'); %norm = h5read(path,'/norm');

I_raw = I; 
I_raw = reshape(I_raw, size(I_raw, 1), size(I_raw, 2), n_acquisitions, []);
I_raw = sum(I_raw, 3);
I_raw(:,1) = [];
dim = size(I_raw);
q = q';

I_raw = log(I_raw);

scan_end = dim(1,2);
X = 1:1:scan_end;
[qq,XX]=meshgrid(X,q);


%% BG substraction 
i_bg = 575; f_bg = 600; %Position of BG substraction

BG(:,:) = mean(I_raw(i_bg:f_bg,:),1);
I_BG(:,:) = I_raw(:,:) - BG(:,:); 
I_flip = flip(I_BG,2);

fig = figure; %
plot(q,I_raw(:,:),'-','linewidth',1);
set(gca, 'ColorOrder', parula(size(I_raw,2)));
grid on 
ylabel('log Intensity [AU]');
xlabel('Scattering vector q [1/nm]');
colorbar
xlim([0.05 4]);
%print(fig,sprintf('Raw_Scan%d%d',scan_file), '-dpng'); 

fig = figure; %
plot(q,I_flip(:,:),'-','linewidth',1);
hold on
plot([q(i_bg,1) q(f_bg,1)],[-0.4 -0.4],'-k','linewidth',2);
hold off
%set(gca, 'ColorOrder', parula(size(I_flip,2)));
grid on 
ylabel('log Intensity [AU]');
xlabel('Scattering vector q [1/nm]');
%ax = gca;
%ax.FontSize = 12; 
ylim([-0.5 8]);
xlim([0.05 4]);

% Adjust the position of the Y-axis label (increase the second value for more space)
yLabelHandle = get(gca, 'YLabel');
set(yLabelHandle, 'Position', get(yLabelHandle, 'Position') + [-0.1 0 0]);
%print(fig,sprintf('BGCorr_Scan%d%d',scan_file), '-dpdf'); 


%% 

%{
A = q;
B = I_raw(:,:);
DOPC_SAXS_347 = [A, B];

save('DOPC_SAXS_347.mat','DOPC_SAXS_347'); %<----- open
writematrix(DOPC_SAXS_347, 'DOPC_SAXS_347.xlsx');
%}

%%  Create 2D map raw data

% Create figure
figure1 = figure;
%figure1.Position = [150 150 700 500]; 

% Create axes
axes1 = axes('Parent',figure1,'CLim',[0 6]);
grid(axes1,'on');
hold(axes1,'on');

% Create axis labels
%zlabel('I ');
ylabel('q [1/nm]');
xlabel('X [spec nr.]');
title(sprintf('spectra BG corrected %d%d',scan_file))

%set(gca,'Xdir','reverse') 

% Create surf
surf(X,q,I_flip,'Parent',axes1,'FaceColor','interp',...
    'EdgeColor','none');

%Colormap
colormap(jet)
colorbar('peer',axes1);

ylim(axes1,[0.5 4]);
xlim(axes1,[X(1) X(end)]);
%xlim(axes1,[X(1) 20]);

%print(figure1,sprintf('ImapBG_Scan%d%d',scan_file), '-dpng'); 



%% Create 2D and 3D map manipulated
%set(0, 'DefaultAxesFontSize', 14)

%X-axis manipulation Nikol
%X = [1:1:scan_end];

initial_X = 4;
final_X = 29; 
stepsize = 8 ;  %um i did lenght/number of steps? 

XX = (X - initial_X) .* stepsize ; %um scale for x-axis 
i_Map = XX(initial_X-2);
f_Map = XX(final_X);
x_spaceing = 0:2*stepsize:f_Map;

figure2 = figure;
%figure2.Position = [500 500 590 450]; 

axes2 = axes('Parent',figure2,'CLim',[0 5]);
ticks_scaleBar = 0:1:5;
hold(axes2,'on');
grid(axes2,'on');

%set(gca,'Xdir','reverse') 

% Create surf
surf(XX,q,I_flip,'Parent',axes2,'FaceColor','interp',...
    'EdgeColor','none');

colormap(jet)
colorbar('peer',axes2);

%X and Y-limits
ylim(axes2,[0.5 4]);
xlim(axes2,[i_Map f_Map]); %xlim(axes2,[XX(1) XX(end)]);
xticks(x_spaceing)

set(gca,'TickDir','both')

h = colorbar;
set(h, 'Ticks', ticks_scaleBar);
%colormap(h, newCmap_nr)
%ylabel(h, 'log Intensity','FontSize', 14);
%h.Label.Position = [0 3 0];

%title('Spatially resoleved SAXS')
%zlabel('I');
ylabel('Scattering vector q [1/nm]','FontSize', 16);
xlabel('Distance from the capillary tip [\mum]','FontSize', 16);

% Adjust the position of the Y-axis label (increase the second value for more space)
yLabelHandle = get(gca, 'YLabel');
set(yLabelHandle, 'Position', get(yLabelHandle, 'Position') + [-0.8 0 0]);

xLabelHandle = get(gca, 'XLabel');
set(xLabelHandle, 'Position', get(xLabelHandle, 'Position') + [0 -.005 0]);

%print(figure2,sprintf('ImapR_Scan%d%d',scan_file), '-dpng'); 
%print(figure2,sprintf('ImapR_Scan%d%d',scan_file), '-dpdf', '-r300','-bestfit'); 
%print(figure2,sprintf('ImapR_Scan%d%d',scan_file), '-dsvg'); 

%% %% NIKOL PEAK ANALYSIS

%To save folder with graphs
%path = pwd; %'/Users/labecka/Documents/3.Lab Work/0.Lu_DOPCDPPC studdy/Scattering Data/10_90' %pwd; % or mention your path 
%myfolder = 'PeakAnalysis525' ;   % new folder to save data
%folder = mkdir([path,filesep,myfolder]) ; %make folder in path
%path  = [path,filesep,myfolder] ;

%Pick graph nr. 
%pos_i = initial_X - 1;
%pos_f = final_X;

pos_i = initial_X - 1;
pos_f = 5;
step = 1; 

SpecNr_vector = pos_i:step:pos_f;

% working with I_baseline shifted Spacing 
D_all= []; I_all = [];

for nr = pos_i:step:pos_f %<----- OBS its spec nr. not um position

    [pks_bs,locs_bs,w,p] = findpeaks(I_flip(100:end-100,nr),q(100:end-100,1),'MinPeakWidth',0.0009,'MinPeakProminence',0.45);

    d_spac = (2.*pi.*10)./locs_bs;
    
    q_spac_all(1:length(locs_bs),nr) = locs_bs(:,1); %extract q-values  of peaks
    D_all(1:length(d_spac),nr) = d_spac(:,1); %exctracting d-spacing of peaks
    I_all(1:length(pks_bs),nr) = pks_bs(:,1); %extracting intensity of peaks
  
    
    fig = figure;
    plot(q, I_flip(:,nr),'linewidth', 1);
    hold on
    plot(locs_bs(:,1),pks_bs(:,1),'o','MarkerSize',12,'LineWidth',2)
    text(locs_bs(:,1)+.08,pks_bs(:,1)+.02, num2str(d_spac(:,1),4))
    grid on 
    hold off
    
    legend( sprintf('%d \\mum\t', XX(nr)),'FontSize', 16);

    xlim([0.05 3.5])
    ylim([-0.5 8])

    xlabel('Scattering vector q [1/nm]','FontSize', 16);
    ylabel('log Intensity [AU]','FontSize', 16);
    yLabelHandle = get(gca, 'YLabel');
    set(yLabelHandle, 'Position', get(yLabelHandle, 'Position') + [-0.05 0 0]);

    %print(fig ,sprintf('SAXSSpectra_nr%d%d',nr), '-dpdf', '-r300','-bestfit'); 
    %print(fig ,sprintf('SAXSSpectra_nr%d%d',nr), '-dsvg');
    
    %To save png figures in separate folder
    %temp=[path,filesep,'fig',num2str(nr),'.png'];
    %saveas(gca,temp);
    
end 

%To save png figures in separate folder
%{
for k = [pos_i:step:pos_f]
    figure(k);
    temp=[path,filesep,'fig',num2str(k),'.png'];
    saveas(gca,temp);
end
%}
%}

%% Relative peak relations

D_all_1 = D_all(1,pos_i:step:pos_f); 
D_all_2 = D_all(2,pos_i:step:pos_f);
D_all_3 = D_all(3,pos_i:step:pos_f);

I_all_1 = I_all(1,pos_i:step:pos_f); 
I_all_2 = I_all(2,pos_i:step:pos_f);

PosIn_V = XX(SpecNr_vector);

f = figure;
plot(SpecNr_vector, D_all_1,'kx','MarkerSize',12,'linewidth',2);
%plot(1:1:size(D_all_1,2), D_all_1,'kx','MarkerSize',12,'linewidth',2)
hold on
plot(SpecNr_vector, D_all_2,'bx','MarkerSize',12,'linewidth',2);
%plot(1:1:size(D_all_2,2), D_all_2,'bx','MarkerSize',12,'linewidth',2)
hold off
xlim([SpecNr_vector(1) SpecNr_vector(end)]);
xticks(SpecNr_vector)
xlabel('X [spec nr.]');
ylabel('d-spacing [Å]');
grid minor 
%print(f ,sprintf('dALL_Scan%d%d',scan_file), '-dpdf', '-r300','-bestfit'); 

f=figure;
plot(SpecNr_vector, I_all_1,'ro','MarkerSize',12,'linewidth',2);
hold on 
plot(SpecNr_vector, I_all_2,'bo','MarkerSize',12,'linewidth',2);
hold off
grid minor 
xlim([SpecNr_vector(1) SpecNr_vector(end)]);
xticks(SpecNr_vector);
%ylim([0 6])
legend('I_1','I_2');
xlabel('X [spec nr.]');
ylabel('log Intensity');
%print(f ,sprintf('IALL_Scan%d%d',scan_file), '-dpdf', '-r300','-bestfit'); 

%% D-spaing plotting 

set(0, 'DefaultAxesFontSize', 14);

yD_space = 48:2:64;
newTics = [PosIn_V(1):stepsize:PosIn_V(end)];

dfig = figure;
plot(PosIn_V(1,2:end), D_all_1(1,2:end),'rx','MarkerSize',12,'linewidth',2);
hold on
plot(PosIn_V(1,end-6:end), D_all_1(1,end-6:end),'kx','MarkerSize',12,'linewidth',2);
plot(PosIn_V(1,end-6:end-2), D_all_2(1,end-6:end-2),'rx','MarkerSize',12,'linewidth',2);
%plot(PosIn_V(1,end-4:end), D_all_2(1,end-4:end),'kx','MarkerSize',12,'linewidth',2)

hold off
grid on
legend('L\alpha','L\alpha','FontSize', 16)

%xticks(PosIn_V)
xticks(x_spaceing)
yticks(yD_space)
xlim([PosIn_V(1) PosIn_V(end)])
ylim([yD_space(1) yD_space(end)])

xlabel('Distance from the capillary tip [\mum]','FontSize', 16);
ylabel('d-spacing [Å]','FontSize', 16);

yLabelHandle = get(gca, 'YLabel');
set(yLabelHandle, 'Position', get(yLabelHandle, 'Position') + [-5 0 0]);

xLabelHandle = get(gca, 'XLabel');
set(xLabelHandle, 'Position', get(xLabelHandle, 'Position') + [0 -.005 0]);

%print(dfig ,sprintf('d_Scan%d%d',scan_file), '-dpdf', '-r300','-bestfit'); 
%print(dfig ,sprintf('d_Scan%d%d',scan_file), '-dsvg');

%% Intensity plotting

yI_space = 0:2:10;

Ifig = figure;
plot(PosIn_V(1,2:end), I_all_1(1,2:end),'ro','MarkerSize',12,'linewidth',2)
hold on
plot(PosIn_V(1,end-6:end), I_all_1(1,end-6:end),'ko','MarkerSize',12,'linewidth',2)
plot(PosIn_V(1,end-6:end-2), I_all_2(1,end-6:end-2),'ro','MarkerSize',12,'linewidth',2)
hold off
grid on 
legend('L\alpha','L\alpha','FontSize', 16)

xlim([PosIn_V(1) PosIn_V(end)])
xticks(x_spaceing)

ylim([yI_space(1) yI_space(end)])
yticks(yI_space)

xlabel('Distance from the capillary tip [\mum]','FontSize', 16);
ylabel('log Intensity','FontSize', 16);

yLabelHandle = get(gca, 'YLabel');
set(yLabelHandle, 'Position', get(yLabelHandle, 'Position') + [-0.5 0 0]);

xLabelHandle = get(gca, 'XLabel');
set(xLabelHandle, 'Position', get(xLabelHandle, 'Position') + [0 -.005 0]);

%print(Ifig ,sprintf('I_Scan%d%d',scan_file), '-dpdf', '-r300','-bestfit'); 
%print(Ifig ,sprintf('I_Scan%d%d',scan_file), '-dsvg');

%% SEMILOG intensity data

yIlog_space = 10.^(yI_space);

I_all_1log = 10.^(I_all_1);
I_all_2log = 10.^(I_all_2);

Ifiglog = figure;
semilogy(PosIn_V(1,2:end), I_all_1log(1,2:end),'ro','MarkerSize',12,'linewidth',2)
hold on
semilogy(PosIn_V(1,end-6:end), I_all_1log(1,end-6:end),'ko','MarkerSize',12,'linewidth',2)
semilogy(PosIn_V(1,end-6:end-2), I_all_2log(1,end-6:end-2),'ro','MarkerSize',12,'linewidth',2)
hold off
grid on 
legend('L\alpha','L\alpha','FontSize', 16)

xlim([PosIn_V(1) PosIn_V(end)])
xticks(x_spaceing)

ylim([10.^0 10.^10])


xlabel('Distance from the capillary tip [\mum]','FontSize', 16);
ylabel('Peak Intensity','FontSize', 16);

yLabelHandle = get(gca, 'YLabel');
set(yLabelHandle, 'Position', get(yLabelHandle, 'Position') + [-0.5 0 0]);

xLabelHandle = get(gca, 'XLabel');
set(xLabelHandle, 'Position', get(xLabelHandle, 'Position') + [0 -.005 0]);

%print(Ifiglog ,sprintf('Ilog_Scan%d%d',scan_file), '-dpdf', '-r300','-bestfit'); 
%print(Ifiglog ,sprintf('Ilog_Scan%d%d',scan_file), '-dsvg');

%% DATA TABLE

PeakDataTable(:,1) = SpecNr_vector(:,:)';
PeakDataTable(:,2) = PosIn_V(:,:)';
PeakDataTable(:,3) = D_all_1(:,:)';
PeakDataTable(:,4) = D_all_2(:,:)';
PeakDataTable(:,5) = I_all_1(:,:)';
PeakDataTable(:,6) = I_all_2(:,:)';

%}