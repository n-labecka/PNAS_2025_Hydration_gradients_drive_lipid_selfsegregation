clear all
close all
clc

set(0, 'DefaultAxesFontSize', 14)

scan_file = 00347; 
n_acquisitions = 6; %#<----------------- CHANGE HERE 

%path = sprintf('/Users/labecka/Documents/BeamData/PSI 2023/cSAXS/analysis/pilatus300k_integ/S%05d_pilatus_2_integ.h5', scan_file);
path = sprintf('DOPC_S00347_pilatus_2_integ.h5', scan_file);
h5disp(path); I = h5read(path,'/I_all'); q = h5read(path,'/q'); %norm = h5read(path,'/norm');

waxs_raw = I; 
waxs_raw = reshape(waxs_raw, size(waxs_raw, 1), size(waxs_raw, 2), n_acquisitions, []);
waxs_raw = sum(waxs_raw, 3);
waxs_raw(:,1) = [];

dim = size(waxs_raw);
scan_end = dim(1,2);
X = 1:1:scan_end;

q = q';
waxs_raw(:,:) = log(waxs_raw(:,:));

%% Data Processing 

%MASKING detector edges  
M1_i = 331; M1_f = 335;
M2_i = 668; M2_f = 673;

for cn = 1:1:scan_end  
    waxs_raw(M1_i:M1_f,cn) = mean(waxs_raw(336:340,cn));
    waxs_raw(M2_i:M2_f,cn) = mean(waxs_raw(650:660,cn)); 
end 

%BG Corecctions 
%i_bg = 550; f_bg = 590; %Position of BG substraction
i_bg = 265; f_bg = 285; %Position of BG substraction

BG(:,:) = mean(waxs_raw(i_bg:f_bg,:),1);
I_waxs(:,:) = waxs_raw(:,:) - BG(:,:); 
I_waxs = flip(I_waxs,2);

%%
fig = figure; %
plot(q,waxs_raw(:,:),'-','linewidth',1)
set(gca, 'ColorOrder', parula(size(waxs_raw,2)));
grid on 
ylabel('log Intensity [AU]');
xlabel('Scattering vector q [1/nm]');
colorbar
xlim([8 24])
%print(fig,sprintf('Raw_Scan%d%d',scan_file), '-dpng'); 

fig = figure;
plot(q,I_waxs(:,:));
hold on
plot([q(i_bg,1) q(f_bg,1)],[-0.4 -0.4],'-k','linewidth',2)
hold off
set(gca, 'ColorOrder', parula(size(I_waxs,2)));
grid on
xlim([6 24])
ylabel('log Intensity [AU]');
xlabel('Scattering vector q [1/nm]');
colorbar
%print(fig, sprintf('WAXS_all_nr%d%d',scan_file), '-dpng'); 


%%
%{
A = q;
B = waxs_raw(:,:);
DOPC_WAXS_347 = [A, B];

save('DOPC_WAXS_347.mat','DOPC_WAXS_347'); %<----- open
writematrix(DOPC_WAXS_347, 'DOPC_WAXS_347.xlsx');
%}

%%
% Create figure
figure1 = figure;
%figure1.Position = [150 150 700 500]; 

% Create axes
axes1 = axes('Parent',figure1,'CLim',[0 1]);
grid(axes1,'on');
hold(axes1,'on');

% Create axis labels
%zlabel('I ');
ylabel('q [1/nm]');
xlabel('X [spec nr.]');
title(sprintf('spectra BG corrected %d%d',scan_file))

%set(gca,'Xdir','reverse') 

% Create surf
surf(X,q,I_waxs,'Parent',axes1,'FaceColor','interp',...
    'EdgeColor','none');

%Colormap
colormap(jet)
colorbar('peer',axes1);

ylim(axes1,[10 22]);
xlim(axes1,[X(1) X(end)]);
%xlim(axes1,[X(1) 20]);

%print(figure1,sprintf('ImapBG_Scan%d%d',scan_file), '-dpng'); 
%% MAPS
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

axes2 = axes('Parent',figure2,'CLim',[0 1]);
hold(axes2,'on');
grid(axes2,'on');

%set(gca,'Xdir','reverse') 

% Create surf
surf(XX,q,I_waxs,'Parent',axes2,'FaceColor','interp',...
    'EdgeColor','none');

colormap(jet)
colorbar('peer',axes2);

%X and Y-limits
ylim(axes2,[10 22]);
xlim(axes2,[i_Map f_Map]); %xlim(axes2,[XX(1) XX(end)]);
xticks(x_spaceing)

set(gca,'TickDir','both')

h = colorbar;
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

%print(figure2,sprintf('WAXSmap_nr%d%d',scan_file), '-dpng'); 
%print(figure2,sprintf('WAXSmap_nr%d%d',scan_file), '-dpdf', '-r300','-bestfit'); 
%print(figure2,sprintf('WAXSmap_nr%d%d',scan_file), '-dsvg'); 

%%
%Pick graph nr. 
pos_i = initial_X - 1;
pos_f = final_X;
step = 1; 

Spec_vector = pos_i:step:pos_f;

% working with I_baseline shifted Spacing 
D_all= []; I_all = [];
for nr = 14 %[pos_i:step:pos_f]
  
    fig = figure(100);
    plot(q, I_waxs(:,nr),'linewidth', 1);
    hold on
    
    [pks_bs,locs_bs,w,p] = findpeaks(I_waxs(1:end,nr),q(1:end,1),'MinPeakWidth',0.0009,'MinPeakProminence',0.45);
    
    d_spac = (2.*pi.*10)./locs_bs;
    legend( sprintf('%d \\mum\t', XX(nr)),'FontSize', 16);

    plot(locs_bs,pks_bs,'o','MarkerSize',12)
    text(locs_bs+.4,pks_bs+.01, num2str(d_spac,4))
    
    ylim([-1 1.5])
    xlim([8 24]);
    xlabel('q [Å^{-1}]');
    ylabel('I [AU]');
    
    grid on 
    %hold off
    
    xlabel('Scattering vector q [1/nm]','FontSize', 16);
    ylabel('log Intensity [AU]','FontSize', 16);
    yLabelHandle = get(gca, 'YLabel');
    set(yLabelHandle, 'Position', get(yLabelHandle, 'Position') + [-0.05 0 0]);
    
    D_all(1:length(d_spac),nr) = d_spac(:,1); %exctracting d-spacing of peaks
    I_all(1:length(pks_bs),nr) = pks_bs(:,1); %extracting intensity of peaks
    
    %print(fig ,sprintf('WAXSpectra_nr%d%d',nr), '-dpdf', '-r300','-bestfit'); 
    %print(fig ,sprintf('WAXSSpectra_nr%d%d',nr), '-dsvg');
    
    %To save png figures in separate folder
    %temp=[path,filesep,'fig',num2str(nr),'.png'];
    %saveas(gca,temp);
end 

%To save png figures in separate folder
%{
for k = [pos_i:steps:pos_f]
    figure(k);
    temp=[path,filesep,'fig',num2str(k),'.png'];
    saveas(gca,temp);
end
%}
%}