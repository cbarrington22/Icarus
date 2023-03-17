
% Plotting results for NOVAC analysis 
outDir = '/Users/charlotteb/Documents/Chapter 2/NOVAC spectra/Results/';
outDirFig = '/Users/charlotteb/Documents/Chapter 2/NOVAC spectra/Results/Figures/';

% Loads results files form linear model 
% AW3
inDirLM_AW3 = '/Users/charlotteb/Documents/Chapter 2/NOVAC spectra/AW3'; % Later add ReRun; % MODEL VERSION 1 IS AW3 cropped frequencies
% 1417 
inFileLM_AW3_1417 = 'ReferenceSpectra_AW3cF_mEst_NOVAC_Masaya_1417.mat'; 
inLM_AW3_1417 = fullfile(inDirLM_AW3, inFileLM_AW3_1417); 
% mEst 
mEstLM_AW3_1417 = load(inLM_AW3_1417,'Results_mEst_real', 'Results_stdErr_real', 'elevAn', 'Results_res_real'); 
mEstLM_AW3_mEst_1417 = mEstLM_AW3_1417.Results_mEst_real;  
mEstLM_AW3_stdErr_1417 = mEstLM_AW3_1417.Results_stdErr_real;  
mEstLM_AW3_angle_1417 = mEstLM_AW3_1417.elevAn; 
Res_AW3_1417 = mEstLM_AW3_1417.Results_res_real;
% 1941
inFileLM_AW3_1941 = 'ReferenceSpectra_AW3cF_mEst_NOVAC_Masaya_1941.mat'; 
inLM_AW3_1941 = fullfile(inDirLM_AW3, inFileLM_AW3_1941); 
% mEst 
mEstLM_AW3_1941 = load(inLM_AW3_1941,'Results_mEst_real', 'Results_stdErr_real', 'elevAn', 'Results_res_real'); 
mEstLM_AW3_mEst_1941 = mEstLM_AW3_1941.Results_mEst_real;  
mEstLM_AW3_stdErr_1941 = mEstLM_AW3_1941.Results_stdErr_real;  
mEstLM_AW3_angle_1941 = mEstLM_AW3_1941.elevAn;  
Res_AW3_1941 = mEstLM_AW3_1941.Results_res_real;

% AW2
inDirLM_AW2 = '/Users/charlotteb/Documents/Chapter 2/NOVAC spectra/AW2'; % Later add ReRun; % MODEL VERSION 1 IS AW2 cropped frequencies
% 1417 
inFileLM_AW2_1417 = 'ReferenceSpectra_AW2cF_mEst_NOVAC_Masaya_1417.mat'; 
inLM_AW2_1417 = fullfile(inDirLM_AW2, inFileLM_AW2_1417); 
% mEst 
mEstLM_AW2_1417 = load(inLM_AW2_1417,'Results_mEst_real', 'Results_stdErr_real', 'elevAn'); 
mEstLM_AW2_mEst_1417 = mEstLM_AW2_1417.Results_mEst_real;  
mEstLM_AW2_stdErr_1417 = mEstLM_AW2_1417.Results_stdErr_real;  
mEstLM_AW2_angle_1417 = mEstLM_AW2_1417.elevAn;  
% 1941
inFileLM_AW2_1941 = 'ReferenceSpectra_AW2cF_mEst_NOVAC_Masaya_1941.mat'; 
inLM_AW2_1941 = fullfile(inDirLM_AW2, inFileLM_AW2_1941); 
% mEst 
mEstLM_AW2_1941 = load(inLM_AW2_1941,'Results_mEst_real', 'Results_stdErr_real', 'elevAn'); 
mEstLM_AW2_mEst_1941 = mEstLM_AW2_1941.Results_mEst_real;  
mEstLM_AW2_stdErr_1941 = mEstLM_AW2_1941.Results_stdErr_real;  
mEstLM_AW2_angle_1941 = mEstLM_AW2_1941.elevAn;  

% AW1
inDirLM_AW1 = '/Users/charlotteb/Documents/Chapter 2/NOVAC spectra/AW1'; % Later add ReRun; % MODEL VERSION 1 IS AW1 cropped frequencies
% 1417 
inFileLM_AW1_1417 = 'ReferenceSpectra_AW1cF_mEst_NOVAC_Masaya_1417.mat'; 
inLM_AW1_1417 = fullfile(inDirLM_AW1, inFileLM_AW1_1417); 
% mEst 
mEstLM_AW1_1417 = load(inLM_AW1_1417,'Results_mEst_real', 'Results_stdErr_real', 'elevAn'); 
mEstLM_AW1_mEst_1417 = mEstLM_AW1_1417.Results_mEst_real;  
mEstLM_AW1_stdErr_1417 = mEstLM_AW1_1417.Results_stdErr_real;  
mEstLM_AW1_angle_1417 = mEstLM_AW1_1417.elevAn;  
% 1941
inFileLM_AW1_1941 = 'ReferenceSpectra_AW1cF_mEst_NOVAC_Masaya_1941.mat'; 
inLM_AW1_1941 = fullfile(inDirLM_AW1, inFileLM_AW1_1941); 
% mEst 
mEstLM_AW1_1941 = load(inLM_AW1_1941,'Results_mEst_real', 'Results_stdErr_real', 'elevAn'); 
mEstLM_AW1_mEst_1941 = mEstLM_AW1_1941.Results_mEst_real;  
mEstLM_AW1_stdErr_1941 = mEstLM_AW1_1941.Results_stdErr_real;  
mEstLM_AW1_angle_1941 = mEstLM_AW1_1941.elevAn;  

% Variables for plotting 
markeredgecolor = 'black'; 
markerfacecolorD =  'white';
markerfacecolor_AW1 = [0.11 0.24 0.81]; % AW1 
markerfacecolor_AW2 = [0.91 0.05 0.05]; % AW2
markerfacecolor_AW3 = [0.09 0.54 0.03]; % AW3
markersize = 4; 
capsize = 0; 
errorbarlinewidthV = 0.5; 
linLim_AW1 = 5e17; 
linLim_AW2 = 5.2e18;
linLim_AW3 = 1.2e18; 
indAng = (7:49); 
% indAng(:,1) = 1; 
% indAng(:,2:44) = indAngx; 
angle_mEst_AW1_1417 = mEstLM_AW1_mEst_1417(2,3:end);
angle_mEst_AW2_1417 = mEstLM_AW2_mEst_1417(2,3:end);
angle_mEst_AW1_1941 = mEstLM_AW1_mEst_1941(2,3:end);
angle_mEst_AW2_1941 = mEstLM_AW2_mEst_1941(2,3:end);
angle_mEst_AW3_1417 = mEstLM_AW3_mEst_1417(2,3:end);
angle_mEst_AW3_1941 = mEstLM_AW3_mEst_1941(2,3:end);
% minX = -70; 
% Axis limits 
% minY = min(cellConD)-1e18;
% maxY = 8.2e18; % max(cellConDconv)+3.5e18;

% FIGURE 1 - AW2 with stdErr only 
figure
% AW3 
% 1417 
subplot(1,2,1)
m = errorbar(mEstLM_AW3_angle_1417, mEstLM_AW3_mEst_1417(2,:), mEstLM_AW3_stdErr_1417(2,:), mEstLM_AW3_stdErr_1417(2,:), 'vertical', 'LineStyle', 'none', 'Marker', '^', 'Color', markerfacecolor_AW3, 'linewidth', errorbarlinewidthV); hold on; 
m.CapSize = capsize; 
m.MarkerFaceColor = markerfacecolor_AW3; 
m.MarkerEdgeColor = markeredgecolor; 
title('14:17 UTC'); hold on; 
xlim([-95 95]);
ylim([-9.9e16 1.5e18]);
ylabel('SCD SO_2 (molec/cm^2)')
set(gca,'YMinorTick','on')
yline(linLim_AW3, ':k')
Fig = gca; Fig.FontSize = 12; set(gcf,'color','w'); Fig.Box = 'on';
% 1941
subplot(1,2,2)
m = errorbar(mEstLM_AW3_angle_1941, mEstLM_AW3_mEst_1941(2,:), mEstLM_AW3_stdErr_1941(2,:), mEstLM_AW3_stdErr_1941(2,:), 'vertical', 'LineStyle', 'none', 'Marker', '^', 'Color', markerfacecolor_AW3, 'linewidth', errorbarlinewidthV); hold on; 
m.CapSize = capsize; 
m.MarkerFaceColor = markerfacecolor_AW3; 
m.MarkerEdgeColor = markeredgecolor; 
title('19:41 UTC'); hold on; 
xlim([-95 95]);
ylim([-9.9e16 1.5e18]);
yline(linLim_AW3, ':k')
set(gca,'YMinorTick','on')
xlabel('Scan elevation angle (°)') 
Fig = gca; Fig.FontSize = 12; set(gcf,'color','w'); Fig.Box = 'on';
%% Make any adjustments 
% SAVES FIGURE
% fnOut = 'Fig1_AW3';
% fname = fullfile(outDirFig, fnOut);
% saveas(gca, fname)

% FIGURE 2 - AW2 and AW1 mEstSO2 with stdErr only 
figure
% AW2 
% 1417 
subplot(2,2,1)
m = errorbar(mEstLM_AW2_angle_1417(:,indAng), mEstLM_AW2_mEst_1417(2,indAng), mEstLM_AW2_stdErr_1417(2,indAng), mEstLM_AW2_stdErr_1417(2,indAng), 'vertical', 'LineStyle', 'none', 'Marker', '^', 'Color', markerfacecolor_AW2, 'linewidth', errorbarlinewidthV); hold on; 
m.CapSize = capsize; 
m.MarkerFaceColor = markerfacecolor_AW2; 
m.MarkerEdgeColor = markeredgecolor; 
title('14:17 UTC'); hold on; 
xlim([-80 80]);
ylim([-4.4e18 -1.9e18]);
% ylabel('SCD SO_2 (molec/cm^2)')
set(gca,'YMinorTick','on')
yline(min(angle_mEst_AW2_1417), '-k')
Fig = gca; Fig.FontSize = 12; set(gcf,'color','w'); Fig.Box = 'on';
% 1941
subplot(2,2,2)
m = errorbar(mEstLM_AW2_angle_1941(:,indAng), mEstLM_AW2_mEst_1941(2,indAng), mEstLM_AW2_stdErr_1941(2,indAng), mEstLM_AW2_stdErr_1941(2,indAng), 'vertical', 'LineStyle', 'none', 'Marker', '^', 'Color', markerfacecolor_AW2, 'linewidth', errorbarlinewidthV); hold on; 
m.CapSize = capsize; 
m.MarkerFaceColor = markerfacecolor_AW2; 
m.MarkerEdgeColor = markeredgecolor; 
title('19:41 UTC'); hold on; 
xlim([-80 80]);
ylim([-4.4e18 -1.9e18]);
% ylabel('SCD SO_2 (molec/cm^2)')
set(gca,'YMinorTick','on')
yline(min(angle_mEst_AW2_1941), '-k')
% xlabel('Scan elevation angle (°)') 
Fig = gca; Fig.FontSize = 12; set(gcf,'color','w'); Fig.Box = 'on';
% AW1
subplot(2,2,3)
m = errorbar(mEstLM_AW1_angle_1417(:,indAng), mEstLM_AW1_mEst_1417(2,indAng), mEstLM_AW1_stdErr_1417(2,indAng), mEstLM_AW1_stdErr_1417(2,indAng), 'vertical', 'LineStyle', 'none', 'Marker', '^', 'Color', markerfacecolor_AW1, 'linewidth', errorbarlinewidthV); hold on; 
m.CapSize = capsize; 
m.MarkerFaceColor = markerfacecolor_AW1; 
m.MarkerEdgeColor = markeredgecolor; 
% title('14e17 UTC'); hold on; 
xlim([-80 80]);
ylim([1e17 12.5e17]);
ylabel('SCD SO_2 (molec/cm^2)')
set(gca,'YMinorTick','on')
yline(min(angle_mEst_AW1_1417), '-k')
Fig = gca; Fig.FontSize = 12; set(gcf,'color','w'); Fig.Box = 'on';
% 1941
subplot(2,2,4)
m = errorbar(mEstLM_AW1_angle_1941(:,indAng), mEstLM_AW1_mEst_1941(2,indAng), mEstLM_AW1_stdErr_1941(2,indAng), mEstLM_AW1_stdErr_1941(2,indAng), 'vertical', 'LineStyle', 'none', 'Marker', '^', 'Color', markerfacecolor_AW1, 'linewidth', errorbarlinewidthV); hold on; 
m.CapSize = capsize; 
m.MarkerFaceColor = markerfacecolor_AW1; 
m.MarkerEdgeColor = markeredgecolor; 
% title('19e41 UTC'); hold on; 
xlim([-80 80]);
ylim([1e17 12.5e17]);
% ylabel('SCD SO_2 (molec/cm^2)')
set(gca,'YMinorTick','on')
yline(min(angle_mEst_AW1_1941), '-k')
xlabel('Scan elevation angle (°)') 
Fig = gca; Fig.FontSize = 12; set(gcf,'color','w'); Fig.Box = 'on';
%% Make any adjustments 
% SAVES FIGURE
% fnOut = 'Fig4_AW2_AW1';
% fname = fullfile(outDirFig, fnOut);
% saveas(gca, fname)

% FIGURE 3 - Residuals AW3 1419 
figure
% AW3 
% 1417 
res = plot(mEstLM_AW3_angle_1417, Res_AW3_1417, '+'); hold on; 
res.Color = markerfacecolor_AW3; 
res = plot(mEstLM_AW3_angle_1941, Res_AW3_1941, 'x'); hold on; 
res.Color = markerfacecolor_AW3; 
% title('14:17 UTC'); 
xlim([-95 95]);
% ylim([-9.9e16 1.5e18]);
% xlim([min(mEstLM_AW3_angle_1417) max((mEstLM_AW3_angle_1417))]);
xlim([-95 95]);
xline(-75, '-k'); 
xline(75, '-k'); 
ylabel('Total error (E)')
Fig = gca; Fig.FontSize = 14; set(gcf,'color','w'); Fig.Box = 'on';
xlabel('Scan elevation angle (°)') 
%% Make any adjustments 
% SAVES FIGURE
% fnOut = 'Fig3_AW3_totalError';
% fname = fullfile(outDirFig, fnOut);
% saveas(gca, fname)

% Offset applied 
% Offset data - using min 
% % 1417 
mEstLM_AW1_mEst_1417_offset = mEstLM_AW1_mEst_1417(2,indAng)-(abs(min(angle_mEst_AW1_1417)));  
mEstLM_AW2_mEst_1417_offset = mEstLM_AW2_mEst_1417(2,indAng)+(abs(min(angle_mEst_AW2_1417)));  
mEstLM_AW3_mEst_1417_offset = mEstLM_AW3_mEst_1417(2,indAng)+(abs(min(angle_mEst_AW3_1417))); 
% 1941 
mEstLM_AW1_mEst_1941_offset = mEstLM_AW1_mEst_1941(2,indAng)-(abs(min(angle_mEst_AW1_1941)));  
mEstLM_AW2_mEst_1941_offset = mEstLM_AW2_mEst_1941(2,indAng)+(abs(min(angle_mEst_AW2_1941)));  
mEstLM_AW3_mEst_1941_offset = mEstLM_AW3_mEst_1941(2,indAng)+(abs(min(angle_mEst_AW3_1941))); 
% % 1417 
% mEstLM_AW1_mEst_1417_offset = mEstLM_AW1_mEst_1417(2, :)-(abs(min(mEstLM_AW1_mEst_1417(2, :))));  
% mEstLM_AW2_mEst_1417_offset = mEstLM_AW2_mEst_1417(2, :)+(abs(min(mEstLM_AW2_mEst_1417(2, :))));  
% mEstLM_AW3_mEst_1417_offset = mEstLM_AW3_mEst_1417(2, :)+(abs(min(mEstLM_AW3_mEst_1417(2, :)))); 
% % 1941 
% mEstLM_AW1_mEst_1941_offset = mEstLM_AW1_mEst_1941(2, :)-(abs(min(mEstLM_AW1_mEst_1941(2, :))));  
% mEstLM_AW2_mEst_1941_offset = mEstLM_AW2_mEst_1941(2, :)+(abs(min(mEstLM_AW2_mEst_1941(2, :))));  
% mEstLM_AW3_mEst_1941_offset = mEstLM_AW3_mEst_1941(2, :)+(abs(min(mEstLM_AW3_mEst_1941(2, :)))); 
% FIGURE 
figure
% All AW 
% 1417 
subplot(1,2,1)
m = errorbar(mEstLM_AW1_angle_1417(:,indAng), mEstLM_AW1_mEst_1417_offset, mEstLM_AW1_stdErr_1417(2,indAng), mEstLM_AW1_stdErr_1417(2,indAng), 'vertical', 'LineStyle', '-', 'Marker', 'd', 'Color', markerfacecolor_AW1, 'linewidth', errorbarlinewidthV); hold on; 
m.CapSize = capsize; 
m.MarkerFaceColor = markerfacecolor_AW1; 
m.MarkerEdgeColor = markeredgecolor; 
m = errorbar(mEstLM_AW2_angle_1417(:,indAng), mEstLM_AW2_mEst_1417_offset, mEstLM_AW2_stdErr_1417(2,indAng), mEstLM_AW2_stdErr_1417(2,indAng), 'vertical', 'LineStyle', '-', 'Marker', 'd', 'Color', markerfacecolor_AW2, 'linewidth', errorbarlinewidthV); hold on; 
m.CapSize = capsize; 
m.MarkerFaceColor = markerfacecolor_AW2; 
m.MarkerEdgeColor = markeredgecolor; 
m = errorbar(mEstLM_AW3_angle_1417(:,indAng), mEstLM_AW3_mEst_1417_offset, mEstLM_AW3_stdErr_1417(2,indAng), mEstLM_AW3_stdErr_1417(2,indAng), 'vertical', 'LineStyle', '-', 'Marker', 'd', 'Color', markerfacecolor_AW3, 'linewidth', errorbarlinewidthV); hold on; 
m.CapSize = capsize; 
m.MarkerFaceColor = markerfacecolor_AW3; 
m.MarkerEdgeColor = markeredgecolor; 
title('14:17 UTC'); hold on; 
xlim([-80 80]);
ylim([0 2.15e18]);
ylabel('SCD SO_2 (molec/cm^2)')
set(gca,'YMinorTick','on')
yl = yline(linLim_AW1-(abs(min(angle_mEst_AW1_1417))), ':'); yl.Color = markerfacecolor_AW1;
yl = yline(linLim_AW2+(abs(min(angle_mEst_AW2_1417))), ':'); yl.Color = markerfacecolor_AW2;
yl = yline(linLim_AW3+(abs(min(angle_mEst_AW3_1417))), ':'); yl.Color = markerfacecolor_AW3;
% yline(linLim_AW3, ':k')
Fig = gca; Fig.FontSize = 12; set(gcf,'color','w'); Fig.Box = 'on';
% 1941  
subplot(1,2,2)
m = errorbar(mEstLM_AW1_angle_1941(:,indAng), mEstLM_AW1_mEst_1941_offset, mEstLM_AW1_stdErr_1941(2,indAng), mEstLM_AW1_stdErr_1941(2,indAng), 'vertical', 'LineStyle', '-', 'Marker', 'd', 'Color', markerfacecolor_AW1, 'linewidth', errorbarlinewidthV); hold on; 
m.CapSize = capsize; 
m.MarkerFaceColor = markerfacecolor_AW1; 
m.MarkerEdgeColor = markeredgecolor; 
m = errorbar(mEstLM_AW2_angle_1941(:,indAng), mEstLM_AW2_mEst_1941_offset, mEstLM_AW2_stdErr_1941(2,indAng), mEstLM_AW2_stdErr_1941(2,indAng), 'vertical', 'LineStyle', '-', 'Marker', 'd', 'Color', markerfacecolor_AW2, 'linewidth', errorbarlinewidthV); hold on; 
m.CapSize = capsize; 
m.MarkerFaceColor = markerfacecolor_AW2; 
m.MarkerEdgeColor = markeredgecolor; 
m = errorbar(mEstLM_AW3_angle_1941(:,indAng), mEstLM_AW3_mEst_1941_offset, mEstLM_AW3_stdErr_1941(2,indAng), mEstLM_AW3_stdErr_1941(2,indAng), 'vertical', 'LineStyle', '-', 'Marker', 'd', 'Color', markerfacecolor_AW3, 'linewidth', errorbarlinewidthV); hold on; 
m.CapSize = capsize; 
m.MarkerFaceColor = markerfacecolor_AW3; 
m.MarkerEdgeColor = markeredgecolor; 
title('19:41 UTC'); hold on; 
xlim([-80 80]);
ylim([0 2.15e18]);
% ylabel('SCD SO_2 (molec/cm^2)')
xlabel('Scan elevation angle (°)') 
set(gca,'YMinorTick','on')
yl = yline(linLim_AW1-(abs(min(angle_mEst_AW1_1941))), ':'); yl.Color = markerfacecolor_AW1;
yl = yline(linLim_AW2+(abs(min(angle_mEst_AW2_1941))), ':'); yl.Color = markerfacecolor_AW2;
yl = yline(linLim_AW3+(abs(min(angle_mEst_AW3_1941))), ':'); yl.Color = markerfacecolor_AW3;
% yline(linLim_AW3, ':k')
Fig = gca; Fig.FontSize = 12; set(gcf,'color','w'); Fig.Box = 'on';
%% Make any adjustments 
% SAVES FIGURE
% fnOut = 'Fig5_AW3_AW2_AW1';
% fname = fullfile(outDirFig, fnOut);
% saveas(gca, fname)

% % Offset applied 
% % Offset data - using mean  
% evInd1 = 1;
% evInd2 = 29; 
% % 1417 
% mEstLM_AW1_mEst_1417_offset = mEstLM_AW1_mEst_1417(2,indAng)-(abs(mean(angle_mEst_AW1_1417(evInd1:evInd2))));  
% mEstLM_AW2_mEst_1417_offset = mEstLM_AW2_mEst_1417(2,indAng)+(abs(mean(angle_mEst_AW2_1417(evInd1:evInd2))));  
% mEstLM_AW3_mEst_1417_offset = mEstLM_AW3_mEst_1417(2,indAng)+(abs(mean(angle_mEst_AW3_1417(evInd1:evInd2))));  
% % 1941 
% evInd1 = 1;
% evInd2 = 25; 
% mEstLM_AW1_mEst_1941_offset = mEstLM_AW1_mEst_1941(2,indAng)-(abs(mean(angle_mEst_AW1_1941(evInd1:evInd2))));  
% mEstLM_AW2_mEst_1941_offset = mEstLM_AW2_mEst_1941(2,indAng)+(abs(mean(angle_mEst_AW2_1941(evInd1:evInd2))));  
% mEstLM_AW3_mEst_1941_offset = mEstLM_AW3_mEst_1941(2,indAng)+(abs(mean(angle_mEst_AW3_1941(evInd1:evInd2))));  
% % FIGURE 
% figure
% % All AW 
% % 1417 
% subplot(1,2,1)
% m = errorbar(mEstLM_AW1_angle_1417(:,indAng), mEstLM_AW1_mEst_1417_offset, mEstLM_AW1_stdErr_1417(2,indAng), mEstLM_AW1_stdErr_1417(2,indAng), 'vertical', 'LineStyle', '-', 'Marker', '^', 'Color', markerfacecolor_AW1, 'linewidth', errorbarlinewidthV); hold on; 
% m.CapSize = capsize; 
% m.MarkerFaceColor = markerfacecolor_AW1; 
% m.MarkerEdgeColor = markeredgecolor; 
% m = errorbar(mEstLM_AW2_angle_1417(:,indAng), mEstLM_AW2_mEst_1417_offset, mEstLM_AW2_stdErr_1417(2,indAng), mEstLM_AW2_stdErr_1417(2,indAng), 'vertical', 'LineStyle', '-', 'Marker', '^', 'Color', markerfacecolor_AW2, 'linewidth', errorbarlinewidthV); hold on; 
% m.CapSize = capsize; 
% m.MarkerFaceColor = markerfacecolor_AW2; 
% m.MarkerEdgeColor = markeredgecolor; 
% m = errorbar(mEstLM_AW3_angle_1417(:,indAng), mEstLM_AW3_mEst_1417_offset, mEstLM_AW3_stdErr_1417(2,indAng), mEstLM_AW3_stdErr_1417(2,indAng), 'vertical', 'LineStyle', '-', 'Marker', '^', 'Color', markerfacecolor_AW3, 'linewidth', errorbarlinewidthV); hold on; 
% m.CapSize = capsize; 
% m.MarkerFaceColor = markerfacecolor_AW3; 
% m.MarkerEdgeColor = markeredgecolor; 
% title('14:17 UTC'); hold on; 
% xlim([-80 80]);
% ylim([-0.3e18 2.3e18]);
% ylabel('SCD SO_2 (molec/cm^2)')
% set(gca,'YMinorTick','on')
% % yline(linLim_AW3, ':k')
% Fig = gca; Fig.FontSize = 12; set(gcf,'color','w'); Fig.Box = 'on';
% % 1941  
% subplot(1,2,2)
% m = errorbar(mEstLM_AW1_angle_1941(:,indAng), mEstLM_AW1_mEst_1941_offset, mEstLM_AW1_stdErr_1941(2,indAng), mEstLM_AW1_stdErr_1941(2,indAng), 'vertical', 'LineStyle', '-', 'Marker', '^', 'Color', markerfacecolor_AW1, 'linewidth', errorbarlinewidthV); hold on; 
% m.CapSize = capsize; 
% m.MarkerFaceColor = markerfacecolor_AW1; 
% m.MarkerEdgeColor = markeredgecolor; 
% m = errorbar(mEstLM_AW2_angle_1941(:,indAng), mEstLM_AW2_mEst_1941_offset, mEstLM_AW2_stdErr_1941(2,indAng), mEstLM_AW2_stdErr_1941(2,indAng), 'vertical', 'LineStyle', '-', 'Marker', '^', 'Color', markerfacecolor_AW2, 'linewidth', errorbarlinewidthV); hold on; 
% m.CapSize = capsize; 
% m.MarkerFaceColor = markerfacecolor_AW2; 
% m.MarkerEdgeColor = markeredgecolor; 
% m = errorbar(mEstLM_AW3_angle_1941(:,indAng), mEstLM_AW3_mEst_1941_offset, mEstLM_AW3_stdErr_1941(2,indAng), mEstLM_AW3_stdErr_1941(2,indAng), 'vertical', 'LineStyle', '-', 'Marker', '^', 'Color', markerfacecolor_AW3, 'linewidth', errorbarlinewidthV); hold on; 
% m.CapSize = capsize; 
% m.MarkerFaceColor = markerfacecolor_AW3; 
% m.MarkerEdgeColor = markeredgecolor; 
% title('19:41 UTC'); hold on; 
% xlim([-80 80]);
% ylim([-2e17 1.4e18]);
% % ylabel('SCD SO_2 (molec/cm^2)')
% xlabel('Scan elevation angle (°)') 
% set(gca,'YMinorTick','on')
% % yline(linLim_AW3, ':k')
% Fig = gca; Fig.FontSize = 12; set(gcf,'color','w'); Fig.Box = 'on';
% %% Make any adjustments 
% % SAVES FIGURE
% fnOut = 'Final_offset_mean';
% fname = fullfile(outDirFig, fnOut);
% saveas(gca, fname)

% figure
% % All AW 
% % 1417 
% subplot(1,2,1)
% 
%    [mEstLM_AW1_mEst_1417_offset(:,indAng); mEstLM_AW1_mEst_1417_offset(:,indAng); mEstLM_AW1_mEst_1417_offset(:,indAng)]
% 
% %% MUST REDEFINE OFFSET 
% one = [];
% two = [];
% three = [];
% for eA = 7:49
%     if mEstLM_AW1_angle_1417(eA) <= 30  % Index 1 to 30 plot AW1 
%        m = errorbar(mEstLM_AW1_angle_1417(:,eA), mEstLM_AW1_mEst_1417_offset(:,eA-6), mEstLM_AW1_stdErr_1417(2,eA), mEstLM_AW1_stdErr_1417(2,eA), 'vertical', 'LineStyle', '-', 'Marker', 'd', 'Color', markerfacecolor_AW1, 'linewidth', errorbarlinewidthV); hold on; 
%        m.CapSize = capsize; 
%        m.MarkerFaceColor = markerfacecolor_AW1; 
%        m.MarkerEdgeColor = markeredgecolor; 
%        one = [one, eA];
%     elseif mEstLM_AW1_angle_1417(eA) > 30 && eA < 45 % Index 1 to 30 plot AW1 
%        m = errorbar(mEstLM_AW3_angle_1417(:,eA), mEstLM_AW3_mEst_1417_offset(:,eA-6), mEstLM_AW3_stdErr_1417(2,eA), mEstLM_AW3_stdErr_1417(2,eA), 'vertical', 'LineStyle', '-', 'Marker', 'd', 'Color', markerfacecolor_AW3, 'linewidth', errorbarlinewidthV); hold on; 
%        m.MarkerFaceColor = markerfacecolor_AW3; 
%        m.MarkerEdgeColor = markeredgecolor; 
%        m.CapSize = capsize; 
%        two = [two, eA];
%     elseif mEstLM_AW1_angle_1417(eA) >= 45 % Index 1 to 30 plot AW1 
%        m = errorbar(mEstLM_AW2_angle_1417(:,eA), mEstLM_AW2_mEst_1417_offset(:,eA-6), mEstLM_AW2_stdErr_1417(2,eA), mEstLM_AW2_stdErr_1417(2,eA), 'vertical', 'LineStyle', '-', 'Marker', 'd', 'Color', markerfacecolor_AW2, 'linewidth', errorbarlinewidthV); hold on; 
%        m.CapSize = capsize; 
%        m.MarkerFaceColor = markerfacecolor_AW2; 
%        m.MarkerEdgeColor = markeredgecolor; 
%        m.CapSize = capsize; 
%        three = [three, eA];
%     else 
%     end 
% end 
% title('14:17 UTC'); hold on; 
% xlim([-80 80]);
% ylim([0 2.15e18]);
% ylabel('SCD SO_2 (molec/cm^2)')
% set(gca,'YMinorTick','on')
yl = yline(linLim_AW1-(abs(min(angle_mEst_AW1_1417))), ':'); yl.Color = markerfacecolor_AW1;
yl = yline(linLim_AW2+(abs(min(angle_mEst_AW2_1417))), ':'); yl.Color = markerfacecolor_AW2;
yl = yline(linLim_AW3+(abs(min(angle_mEst_AW3_1417))), ':'); yl.Color = markerfacecolor_AW3;
% % yline(linLim_AW3, ':k')
% Fig = gca; Fig.FontSize = 12; set(gcf,'color','w'); Fig.Box = 'on';
% newVrY = [mEstLM_AW1_mEst_1417_offset(:,one-6), mEstLM_AW2_mEst_1417_offset(:,two-6), mEstLM_AW3_mEst_1417_offset(:,three-6)]; 
% newVrX = [mEstLM_AW1_angle_1417(:,one), mEstLM_AW2_angle_1417(:,two), mEstLM_AW3_angle_1417(:,three)]; 
% % newVrErr = [mEstLM_AW1_stdErr_1417(:,one), mEstLM_AW2_stdErr_1417(:,two), mEstLM_AW3_stdErr_1417(:,three)]; 
% plot(newVrX, newVrY, '-k'); 
% 
% % 1941  
% subplot(1,2,2)
% for eA = 7:49
%     if mEstLM_AW1_angle_1941(eA) <= 0  % Index 1 to 30 plot AW1 
%        m = errorbar(mEstLM_AW1_angle_1941(:,eA), mEstLM_AW1_mEst_1941_offset(:,eA-6), mEstLM_AW1_stdErr_1941(2,eA), mEstLM_AW1_stdErr_1941(2,eA), 'vertical', 'LineStyle', '-', 'Marker', 'd', 'Color', markerfacecolor_AW1, 'linewidth', errorbarlinewidthV); hold on; 
%        m.CapSize = capsize; 
%        m.MarkerFaceColor = markerfacecolor_AW1; 
%        m.MarkerEdgeColor = markeredgecolor; 
%        m.CapSize = capsize; 
%     elseif mEstLM_AW1_angle_1941(eA) > 0 && eA < 75 % Index 1 to 30 plot AW1 
%        m = errorbar(mEstLM_AW3_angle_1941(:,eA), mEstLM_AW3_mEst_1941_offset(:,eA-6), mEstLM_AW3_stdErr_1941(2,eA), mEstLM_AW3_stdErr_1941(2,eA), 'vertical', 'LineStyle', '-', 'Marker', 'd', 'Color', markerfacecolor_AW3, 'linewidth', errorbarlinewidthV); hold on; 
%        m.MarkerFaceColor = markerfacecolor_AW3; 
%        m.MarkerEdgeColor = markeredgecolor; 
%        m.CapSize = capsize; 
% %     elseif mEstLM_AW1_angle_1941(eA) >= 45 % Index 1 to 30 plot AW1 
% %        m = errorbar(mEstLM_AW2_angle_1941(:,eA), mEstLM_AW2_mEst_1941_offset(:,eA-6), mEstLM_AW2_stdErr_1941(2,eA), mEstLM_AW2_stdErr_1941(2,eA), 'vertical', 'LineStyle', '-', 'Marker', 'd', 'Color', markerfacecolor_AW2, 'linewidth', errorbarlinewidthV); hold on; 
% %        m.CapSize = capsize; 
% %        m.MarkerFaceColor = markerfacecolor_AW2; 
% %        m.MarkerEdgeColor = markeredgecolor; 
%     else 
%     end 
% end 
% title('19:41 UTC'); hold on; 
% xlim([-80 80]);
% ylim([0 2.15e18]);
% % ylabel('SCD SO_2 (molec/cm^2)')
% xlabel('Scan elevation angle (°)') 
% set(gca,'YMinorTick','on')
% % yl = yline(linLim_AW1-(abs(min(angle_mEst_AW1_1941))), ':'); yl.Color = markerfacecolor_AW1;
% % yl = yline(linLim_AW2+(abs(min(angle_mEst_AW2_1941))), ':'); yl.Color = markerfacecolor_AW2;
% % yl = yline(linLim_AW3+(abs(min(angle_mEst_AW3_1941))), ':'); yl.Color = markerfacecolor_AW3;
% % yline(linLim_AW3, ':k')
% Fig = gca; Fig.FontSize = 12; set(gcf,'color','w'); Fig.Box = 'on';
% %% Make any adjustments 
% % SAVES FIGURE
% % fnOut = 'Fig6_combined';
% % fname = fullfile(outDirFig, fnOut);
% % saveas(gca, fname)
% 
% 
% 
% 
% 
% 
% 
% 
% 
% Loads DOAS data 
addpath '/Users/charlotteb/Documents/Chapter 2/NOVAC spectra/'
inDirD = '/Users/charlotteb/Documents/Chapter 2/NOVAC spectra/'; 
inFileD_1417 = 'mDOASResults_2022-09-27T11-37-44_D2J2375_140319_1417_0.csv'; 
inFileD_1941 = 'mDOASResults_2022-09-27T11-40-34_D2J2375_140319_1941_0.csv'; 
hlinesD = 1; 
% 1417 
inD = fullfile(inDirD, inFileD_1417); 
% Loads all data in results file
[SpectrumID_1417, FitCoeff_6101_1417,	FitCoeffError_6101_1417,	Shift_6101_1417,	Squeeze_6101_1417, FitCoeff_Ring_1417, FitCoeffError_Ring_1417,	Shift_Ring_1417,...
    Squeeze_Ring_1417, FitCoeff_FLMS195681_O3_Voigt_223K_HP500_PPMM_1417,	FitCoeffError_FLMS195681_O3_Voigt_223K_HP500_PPMM_1417, Shift_FLMS195681_O3_Voigt_223K_HP500_PPMM_1417,...
    Squeeze_FLMS195681_O3_Voigt_223K_HP500_PPMM_1417, FitCoeff_FLMS195681_SO2_Bogumil_293K_HP500_MOLECCM2_1417,...
    FitCoeffError_FLMS195681_SO2_Bogumil_293K_HP500_PPMM_1417, nSO2D_1417, Shift_FLMS195681_SO2_Bogumil_293K_HP500_PPMM_1417,	Squeeze_FLMS195681_SO2_Bogumil_293K_HP500_PPMM_1417,...
    WavelengthRangeLow_1417,	WavelengthRangeHigh_1417, FitChi2_1417, TimeStamp_1417, Latitude_1417, Longitude_1417, Altitude_1417,	Speed_1417, Course_1417, GPSWarnCode_1417,	GPSQuality_1417,	ElevationAngle_1417,	AzimuthAngle_1417,...
    ExposureTime_1417, Exposures_1417, MaxIntensity_1417, MaxIntensityFitRange_1417, FileName_1417, SpectrometerType_1417, SpectrometerSerialNumber_1417,...
    SpectrometerChannel_1417, Remark_1417] = textread(inD,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%s%s%s%s%f%s', 'headerlines', hlinesD, 'delimiter', ',');  
% 1941 
inD = fullfile(inDirD, inFileD_1941); 
[SpectrumID_1941, FitCoeff_6101_1941,	FitCoeffError_6101_1941,	Shift_6101_1941,	Squeeze_6101_1941, FitCoeff_Ring_1941, FitCoeffError_Ring_1941,	Shift_Ring_1941,...
    Squeeze_Ring_1941, FitCoeff_FLMS195681_O3_Voigt_223K_HP500_PPMM_1941,	FitCoeffError_FLMS195681_O3_Voigt_223K_HP500_PPMM_1941, Shift_FLMS195681_O3_Voigt_223K_HP500_PPMM_1941,...
    Squeeze_FLMS195681_O3_Voigt_223K_HP500_PPMM_1941, FitCoeff_FLMS195681_SO2_Bogumil_293K_HP500_MOLECCM2_1941,...
    FitCoeffError_FLMS195681_SO2_Bogumil_293K_HP500_PPMM_1941, nSO2D_1941, Shift_FLMS195681_SO2_Bogumil_293K_HP500_PPMM_1941,	Squeeze_FLMS195681_SO2_Bogumil_293K_HP500_PPMM_1941,...
    WavelengthRangeLow_1941,	WavelengthRangeHigh_1941, FitChi2_1941, TimeStamp_1941, Latitude_1941, Longitude_1941, Altitude_1941,	Speed_1941, Course_1941, GPSWarnCode_1941,	GPSQuality_1941,	ElevationAngle_1941,	AzimuthAngle_1941,...
    ExposureTime_1941, Exposures_1941, MaxIntensity_1941, MaxIntensityFitRange_1941, FileName_1941, SpectrometerType_1941, SpectrometerSerialNumber_1941,...
    SpectrometerChannel_1941, Remark_1941] = textread(inD,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%s%s%s%s%f%s', 'headerlines', hlinesD, 'delimiter', ',');  
 
SO2D_1417 = FitCoeff_FLMS195681_SO2_Bogumil_293K_HP500_MOLECCM2_1417;  
SO2ErrD_1417 = FitCoeffError_FLMS195681_SO2_Bogumil_293K_HP500_PPMM_1417; 

SO2D_1941 = FitCoeff_FLMS195681_SO2_Bogumil_293K_HP500_MOLECCM2_1941;  
SO2ErrD_1941 = FitCoeffError_FLMS195681_SO2_Bogumil_293K_HP500_PPMM_1941; 

% markerfacecolorD =  'white';
% % % FIGURE 1 - AW2 with stdErr only 
% % figure
% % % AW3 
% % % 1417 
subplot(1,2,1)
% % % m = errorbar(mEstLM_AW3_angle_1417(:,indAng), mEstLM_AW3_mEst_1417(2,indAng), mEstLM_AW3_stdErr_1417(2,indAng), mEstLM_AW3_stdErr_1417(2,indAng), 'vertical', 'LineStyle', '-', 'Marker', '^', 'Color', markerfacecolor_AW3, 'linewidth', errorbarlinewidthV); hold on; 
% % % m.CapSize = capsize; 
% % % m.MarkerFaceColor = markerfacecolor_AW3; 
% % % m.MarkerEdgeColor = markeredgecolor; 
% % % % DOAS 
% d = errorbar(ElevationAngle_1417(3:end,:), nSO2D_1417(3:end,:), SO2ErrD_1417(3:end,:), SO2ErrD_1417(3:end,:), 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);
% d.MarkerFaceColor = markerfacecolorD; 
% d.MarkerEdgeColor = markeredgecolor; 
% d.CapSize = capsize; 
% % 
% d = errorbar(ElevationAngle_1417(32:34,:), nSO2D_1417(32:34,:), SO2ErrD_1417(32:34,:), SO2ErrD_1417(32:34,:), 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);
% d.MarkerFaceColor = markerfacecolorD; 
% d.MarkerEdgeColor = markeredgecolor; 
% d.CapSize = capsize; 
% d.MarkerSize = 8; 
% % 
% % 
% d = errorbar(ElevationAngle_1417(39:41,:), nSO2D_1417(39:41,:), SO2ErrD_1417(39:41,:), SO2ErrD_1417(39:41,:), 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);
% d.MarkerFaceColor = markerfacecolorD; 
% d.MarkerEdgeColor = markeredgecolor; 
% d.CapSize = capsize; 
% d.MarkerSize = 8; 
% % % % 
% d = errorbar(ElevationAngle_1417(42:end,:), nSO2D_1417(42:end,:), SO2ErrD_1417(42:end,:), SO2ErrD_1417(42:end,:), 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);
% d.MarkerFaceColor = markerfacecolorD; 
% d.MarkerEdgeColor = markeredgecolor; 
% d.CapSize = capsize; 
% d.MarkerSize = 10; 

% % m = errorbar(mEstLM_AW3_angle_1417(:,indAng), mEstLM_AW3_mEst_1417(2,indAng), mEstLM_AW3_stdErr_1417(2,indAng), mEstLM_AW3_stdErr_1417(2,indAng), 'vertical', 'LineStyle', '-', 'Marker', '^', 'Color', markerfacecolor_AW3, 'linewidth', errorbarlinewidthV); hold on; 
% % m.CapSize = capsize; 
% % m.MarkerFaceColor = markerfacecolor_AW3; 
% % m.MarkerEdgeColor = markeredgecolor; 
% % % DOAS 
d = errorbar(ElevationAngle_1417(5:47,:), nSO2D_1417(5:47,:), SO2ErrD_1417(5:47,:), SO2ErrD_1417(5:47,:), 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);
d.MarkerFaceColor = markerfacecolorD; 
d.MarkerEdgeColor = markeredgecolor; 
d.CapSize = capsize; 
% 
d = errorbar(ElevationAngle_1417(32:34,:), nSO2D_1417(32:34,:), SO2ErrD_1417(32:34,:), SO2ErrD_1417(32:34,:), 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);
d.MarkerFaceColor = markerfacecolorD; 
d.MarkerEdgeColor = markeredgecolor; 
d.CapSize = capsize; 
d.MarkerSize = 8; 
% 
% 
d = errorbar(ElevationAngle_1417(39:41,:), nSO2D_1417(39:41,:), SO2ErrD_1417(39:41,:), SO2ErrD_1417(39:41,:), 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);
d.MarkerFaceColor = markerfacecolorD; 
d.MarkerEdgeColor = markeredgecolor; 
d.CapSize = capsize; 
d.MarkerSize = 8; 
% % % 
d = errorbar(ElevationAngle_1417(42:47,:), nSO2D_1417(42:47,:), SO2ErrD_1417(42:47,:), SO2ErrD_1417(42:47,:), 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);
d.MarkerFaceColor = markerfacecolorD; 
d.MarkerEdgeColor = markeredgecolor; 
d.CapSize = capsize; 
d.MarkerSize = 10; 

% % title('14:17 UTC'); hold on; 
% % xlim([-95 95]);
% % ylim([-9.9e16 2.8e18]);
% % ylabel('SCD SO_2 (molec/cm^2)')
% % set(gca,'YMinorTick','on')
% % % yl = yline(linLim_AW1, ':'); yl.Color = markerfacecolor_AW1;
% % % yl = yline(linLim_AW2, ':'); yl.Color = markerfacecolor_AW2;
% % yl = yline(linLim_AW3, ':'); yl.Color = markerfacecolor_AW3;
% % Fig = gca; Fig.FontSize = 12; set(gcf,'color','w'); Fig.Box = 'on';
% % % 1941
subplot(1,2,2)
% % m = errorbar(mEstLM_AW3_angle_1941(:,indAng), mEstLM_AW3_mEst_1941(2,indAng), mEstLM_AW3_stdErr_1941(2,indAng), mEstLM_AW3_stdErr_1941(2,indAng), 'vertical', 'LineStyle', '-', 'Marker', '^', 'Color', markerfacecolor_AW3, 'linewidth', errorbarlinewidthV); hold on; 
% % m.CapSize = capsize; 
% % m.MarkerFaceColor = markerfacecolor_AW3; 
% % m.MarkerEdgeColor = markeredgecolor; 
% % % DOAS 
d = errorbar(ElevationAngle_1941(5:47,:), nSO2D_1941(5:47,:), SO2ErrD_1941(5:47,:), SO2ErrD_1941(5:47,:), 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);
d.MarkerFaceColor = markerfacecolorD; 
d.MarkerEdgeColor = markeredgecolor; 
d.CapSize = capsize; 
% % title('19:41 UTC'); hold on; 
% % xlim([-95 95]);
% % ylim([-9.9e16 2.8e18]);
% % % yl = yline(linLim_AW1, ':'); yl.Color = markerfacecolor_AW1;
% % % yl = yline(linLim_AW2, ':'); yl.Color = markerfacecolor_AW2;
% % yl = yline(linLim_AW3, ':'); yl.Color = markerfacecolor_AW3;
% % set(gca,'YMinorTick','on')
% % xlabel('Scan elevation angle (°)') 
% % Fig = gca; Fig.FontSize = 12; set(gcf,'color','w'); Fig.Box = 'on';
% % %% Make any adjustments 
% % % SAVES FIGURE
% % fnOut = 'Fig25_AW3_DOASnormSO2';
% % fname = fullfile(outDirFig, fnOut);
% % saveas(gca, fname)
% % 
% % % FIGURE F1 FOR APPENDIX 
% % % FIGURE 1 - AW2 with stdErr only 
% % figure
% % % AW3 
% % % 1417 
% % subplot(1,2,1)
% % m = errorbar(mEstLM_AW3_angle_1417(:,indAng), mEstLM_AW3_mEst_1417(2,indAng), mEstLM_AW3_stdErr_1417(2,indAng), mEstLM_AW3_stdErr_1417(2,indAng), 'vertical', 'LineStyle', '-', 'Marker', '^', 'Color', markerfacecolor_AW3, 'linewidth', errorbarlinewidthV); hold on; 
% % m.CapSize = capsize; 
% % m.MarkerFaceColor = markerfacecolor_AW3; 
% % m.MarkerEdgeColor = markeredgecolor; 
% % % DOAS 
% % d = errorbar(ElevationAngle_1417(3:end,:), SO2D_1417(3:end,:), SO2ErrD_1417(3:end,:), SO2ErrD_1417(3:end,:), 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);
% % d.MarkerFaceColor = markerfacecolorD; 
% % d.MarkerEdgeColor = markeredgecolor; 
% % d.CapSize = capsize; 
% % % 
% % d = errorbar(ElevationAngle_1417(32:34,:), SO2D_1417(32:34,:), SO2ErrD_1417(32:34,:), SO2ErrD_1417(32:34,:), 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);
% % d.MarkerFaceColor = markerfacecolorD; 
% % d.MarkerEdgeColor = markeredgecolor; 
% % d.CapSize = capsize; 
% % d.MarkerSize = 8; 
% % % 
% % % 
% % d = errorbar(ElevationAngle_1417(39:41,:), SO2D_1417(39:41,:), SO2ErrD_1417(39:41,:), SO2ErrD_1417(39:41,:), 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);
% % d.MarkerFaceColor = markerfacecolorD; 
% % d.MarkerEdgeColor = markeredgecolor; 
% % d.CapSize = capsize; 
% % d.MarkerSize = 8; 
% % % 
% % d = errorbar(ElevationAngle_1417(42:end,:), SO2D_1417(42:end,:), SO2ErrD_1417(42:end,:), SO2ErrD_1417(42:end,:), 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);
% % d.MarkerFaceColor = markerfacecolorD; 
% % d.MarkerEdgeColor = markeredgecolor; 
% % d.CapSize = capsize; 
% % d.MarkerSize = 10; 
% % title('14:17 UTC'); hold on; 
% % xlim([-95 95]);
% % ylim([-9.9e16 2.8e18]);
% % ylabel('SCD SO_2 (molec/cm^2)')
% % set(gca,'YMinorTick','on')
% % % yl = yline(linLim_AW1, ':'); yl.Color = markerfacecolor_AW1;
% % % yl = yline(linLim_AW2, ':'); yl.Color = markerfacecolor_AW2;
% % yl = yline(linLim_AW3, ':'); yl.Color = markerfacecolor_AW3;
% % Fig = gca; Fig.FontSize = 12; set(gcf,'color','w'); Fig.Box = 'on';
% % % 1941
% % subplot(1,2,2)
% % m = errorbar(mEstLM_AW3_angle_1941(:,indAng), mEstLM_AW3_mEst_1941(2,indAng), mEstLM_AW3_stdErr_1941(2,indAng), mEstLM_AW3_stdErr_1941(2,indAng), 'vertical', 'LineStyle', '-', 'Marker', '^', 'Color', markerfacecolor_AW3, 'linewidth', errorbarlinewidthV); hold on; 
% % m.CapSize = capsize; 
% % m.MarkerFaceColor = markerfacecolor_AW3; 
% % m.MarkerEdgeColor = markeredgecolor; 
% % % DOAS 
% % d = errorbar(ElevationAngle_1941(3:end,:), SO2D_1941(3:end,:), SO2ErrD_1941(3:end,:), SO2ErrD_1941(3:end,:), 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);
% % d.MarkerFaceColor = markerfacecolorD; 
% % d.MarkerEdgeColor = markeredgecolor; 
% % d.CapSize = capsize; 
% % title('19:41 UTC'); hold on; 
% % xlim([-95 95]);
% % ylim([-9.9e16 2.8e18]);
% % % yl = yline(linLim_AW1, ':'); yl.Color = markerfacecolor_AW1;
% % % yl = yline(linLim_AW2, ':'); yl.Color = markerfacecolor_AW2;
% % yl = yline(linLim_AW3, ':'); yl.Color = markerfacecolor_AW3;
% % set(gca,'YMinorTick','on')
% % xlabel('Scan elevation angle (°)') 
% % Fig = gca; Fig.FontSize = 12; set(gcf,'color','w'); Fig.Box = 'on';
% % %% Make any adjustments 
% % % SAVES FIGURE
% % fnOut = 'FigF1_AW3_DOASfitcoeSO2';
% % fname = fullfile(outDirFig, fnOut);
% % saveas(gca, fname)
%%%%%%%%%% 
%% Loads DOAS data 
addpath '/Users/charlotteb/Documents/Chapter 2/NOVAC spectra/'
inDirD = '/Users/charlotteb/Documents/Chapter 2/NOVAC spectra/'; 
inFileD_1417 = 'mDOASResults_2022-09-27T11-37-44_D2J2375_140319_1417_0.csv'; 
inFileD_1941 = 'mDOASResults_2022-09-27T11-40-34_D2J2375_140319_1941_0.csv'; 
hlinesD = 1; 
% 1417 
inD = fullfile(inDirD, inFileD_1417); 
% Loads all data in results file
[SpectrumID_1417, FitCoeff_6101_1417,	FitCoeffError_6101_1417,	Shift_6101_1417,	Squeeze_6101_1417, FitCoeff_Ring_1417, FitCoeffError_Ring_1417,	Shift_Ring_1417,...
    Squeeze_Ring_1417, FitCoeff_FLMS195681_O3_Voigt_223K_HP500_PPMM_1417,	FitCoeffError_FLMS195681_O3_Voigt_223K_HP500_PPMM_1417, Shift_FLMS195681_O3_Voigt_223K_HP500_PPMM_1417,...
    Squeeze_FLMS195681_O3_Voigt_223K_HP500_PPMM_1417, FitCoeff_FLMS195681_SO2_Bogumil_293K_HP500_MOLECCM2_1417,...
    FitCoeffError_FLMS195681_SO2_Bogumil_293K_HP500_PPMM_1417, nSO2D_1417, Shift_FLMS195681_SO2_Bogumil_293K_HP500_PPMM_1417,	Squeeze_FLMS195681_SO2_Bogumil_293K_HP500_PPMM_1417,...
    WavelengthRangeLow_1417,	WavelengthRangeHigh_1417, FitChi2_1417, TimeStamp_1417, Latitude_1417, Longitude_1417, Altitude_1417,	Speed_1417, Course_1417, GPSWarnCode_1417,	GPSQuality_1417,	ElevationAngle_1417,	AzimuthAngle_1417,...
    ExposureTime_1417, Exposures_1417, MaxIntensity_1417, MaxIntensityFitRange_1417, FileName_1417, SpectrometerType_1417, SpectrometerSerialNumber_1417,...
    SpectrometerChannel_1417, Remark_1417] = textread(inD,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%s%s%s%s%f%s', 'headerlines', hlinesD, 'delimiter', ',');  
% 1941 
inD = fullfile(inDirD, inFileD_1941); 
[SpectrumID_1941, FitCoeff_6101_1941,	FitCoeffError_6101_1941,	Shift_6101_1941,	Squeeze_6101_1941, FitCoeff_Ring_1941, FitCoeffError_Ring_1941,	Shift_Ring_1941,...
    Squeeze_Ring_1941, FitCoeff_FLMS195681_O3_Voigt_223K_HP500_PPMM_1941,	FitCoeffError_FLMS195681_O3_Voigt_223K_HP500_PPMM_1941, Shift_FLMS195681_O3_Voigt_223K_HP500_PPMM_1941,...
    Squeeze_FLMS195681_O3_Voigt_223K_HP500_PPMM_1941, FitCoeff_FLMS195681_SO2_Bogumil_293K_HP500_MOLECCM2_1941,...
    FitCoeffError_FLMS195681_SO2_Bogumil_293K_HP500_PPMM_1941, nSO2D_1941, Shift_FLMS195681_SO2_Bogumil_293K_HP500_PPMM_1941,	Squeeze_FLMS195681_SO2_Bogumil_293K_HP500_PPMM_1941,...
    WavelengthRangeLow_1941,	WavelengthRangeHigh_1941, FitChi2_1941, TimeStamp_1941, Latitude_1941, Longitude_1941, Altitude_1941,	Speed_1941, Course_1941, GPSWarnCode_1941,	GPSQuality_1941,	ElevationAngle_1941,	AzimuthAngle_1941,...
    ExposureTime_1941, Exposures_1941, MaxIntensity_1941, MaxIntensityFitRange_1941, FileName_1941, SpectrometerType_1941, SpectrometerSerialNumber_1941,...
    SpectrometerChannel_1941, Remark_1941] = textread(inD,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%s%s%s%s%f%s', 'headerlines', hlinesD, 'delimiter', ',');  
 
SO2D_1417 = FitCoeff_FLMS195681_SO2_Bogumil_293K_HP500_MOLECCM2_1417;  
SO2ErrD_1417 = FitCoeffError_FLMS195681_SO2_Bogumil_293K_HP500_PPMM_1417; 

SO2D_1941 = FitCoeff_FLMS195681_SO2_Bogumil_293K_HP500_MOLECCM2_1941;  
SO2ErrD_1941 = FitCoeffError_FLMS195681_SO2_Bogumil_293K_HP500_PPMM_1941; 





% FIGURE 1 - AW2 with stdErr only 
figure
% AW3 
% 1417 
subplot(2,2,1)
yyaxis left 
ylabel('mEst_I_0'); hold on 
ylim([min(mEstLM_AW3_mEst_1417(1,5:end)) max(mEstLM_AW3_mEst_1417(1,5:end))]);
m = errorbar(mEstLM_AW3_angle_1417(:,1), mEstLM_AW3_mEst_1417(1,1), mEstLM_AW3_stdErr_1417(1,1), mEstLM_AW3_stdErr_1417(1,1), 'vertical', 'LineStyle', 'none', 'Marker', '^', 'Color', markerfacecolor_AW3, 'linewidth', errorbarlinewidthV); hold on; 
m.CapSize = capsize; 
m.MarkerFaceColor = markerfacecolor_AW3; 
m.MarkerEdgeColor = markeredgecolor; 
m = errorbar(mEstLM_AW3_angle_1417(:,5:end), mEstLM_AW3_mEst_1417(1,5:end), mEstLM_AW3_stdErr_1417(1,5:end), mEstLM_AW3_stdErr_1417(1,5:end), 'vertical', 'LineStyle', 'none', 'Marker', '^', 'Color', markerfacecolor_AW3, 'linewidth', errorbarlinewidthV); hold on; 
m.CapSize = capsize; 
m.MarkerFaceColor = markerfacecolor_AW3; 
m.MarkerEdgeColor = markeredgecolor; 
yyaxis right 
ylabel('I_0 fit coefficient'); hold on 
d = errorbar(ElevationAngle_1417(5:47,:), FitCoeff_6101_1417(5:47,:), FitCoeffError_6101_1417(5:47,:), FitCoeffError_6101_1417(5:47,:), 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV); hold on 
d.MarkerFaceColor = markerfacecolorD; 
d.MarkerEdgeColor = markeredgecolor; 
d.CapSize = capsize; 
d = errorbar(ElevationAngle_1417(32:34,:), FitCoeff_6101_1417(32:34,:), FitCoeffError_6101_1417(32:34,:), FitCoeffError_6101_1417(32:34,:), 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV); hold on 
d.MarkerFaceColor = markerfacecolorD; 
d.MarkerEdgeColor = markeredgecolor; 
d.CapSize = capsize; 
d.MarkerSize = 8; 
d = errorbar(ElevationAngle_1417(39:41,:), FitCoeff_6101_1417(39:41,:), FitCoeffError_6101_1417(39:41,:), FitCoeffError_6101_1417(39:41,:), 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);
d.MarkerFaceColor = markerfacecolorD; 
d.MarkerEdgeColor = markeredgecolor; 
d.CapSize = capsize; 
d.MarkerSize = 8; 
d = errorbar(ElevationAngle_1417(42:47,:), FitCoeff_6101_1417(42:47,:), FitCoeffError_6101_1417(42:47,:), FitCoeffError_6101_1417(42:47,:), 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);
d.MarkerFaceColor = markerfacecolorD; 
d.MarkerEdgeColor = markeredgecolor; 
d.CapSize = capsize; 
d.MarkerSize = 10; 
ylim([min(FitCoeff_6101_1417) max(FitCoeff_6101_1417)]);
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
% xlim([-95 95]);
xlim([-78 78]);
% Insert formatting 
Fig = gca; Fig.FontSize = 12; set(gcf,'color','w'); Fig.Box = 'on';
subplot(2,2,2)
yyaxis left 
ylabel('mEst_O_3 (molec/cm^2)'); hold on 
m = errorbar(mEstLM_AW3_angle_1417(:,1), mEstLM_AW3_mEst_1417(3,1), mEstLM_AW3_stdErr_1417(3,1), mEstLM_AW3_stdErr_1417(3,1), 'vertical', 'LineStyle', 'none', 'Marker', '^', 'Color', markerfacecolor_AW3, 'linewidth', errorbarlinewidthV); hold on; 
m.CapSize = capsize; 
m.MarkerFaceColor = markerfacecolor_AW3; 
m.MarkerEdgeColor = markeredgecolor; 
m = errorbar(mEstLM_AW3_angle_1417(:,5:end), mEstLM_AW3_mEst_1417(3,5:end), mEstLM_AW3_stdErr_1417(3,5:end), mEstLM_AW3_stdErr_1417(3,5:end), 'vertical', 'LineStyle', 'none', 'Marker', '^', 'Color', markerfacecolor_AW3, 'linewidth', errorbarlinewidthV); hold on; 
m.CapSize = capsize; 
m.MarkerFaceColor = markerfacecolor_AW3; 
m.MarkerEdgeColor = markeredgecolor; 
ylim([min(mEstLM_AW3_mEst_1417(3,5:end)) max(mEstLM_AW3_mEst_1417(3,5:end))]);

yyaxis right 
ylabel('O_3 fit coefficient (molec/cm^2)'); hold on 
d = errorbar(ElevationAngle_1417(5:47,:), FitCoeff_FLMS195681_O3_Voigt_223K_HP500_PPMM_1417(5:47,:), FitCoeffError_FLMS195681_O3_Voigt_223K_HP500_PPMM_1417(5:47,:), FitCoeffError_FLMS195681_O3_Voigt_223K_HP500_PPMM_1417(5:47,:), 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);hold on 
d.MarkerFaceColor = markerfacecolorD; 
d.MarkerEdgeColor = markeredgecolor; 
d.CapSize = capsize; 
d = errorbar(ElevationAngle_1417(32:34,:), FitCoeff_FLMS195681_O3_Voigt_223K_HP500_PPMM_1417(32:34,:), FitCoeffError_FLMS195681_O3_Voigt_223K_HP500_PPMM_1417(32:34,:), FitCoeffError_FLMS195681_O3_Voigt_223K_HP500_PPMM_1417(32:34,:), 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);
d.MarkerFaceColor = markerfacecolorD; 
d.MarkerEdgeColor = markeredgecolor; 
d.CapSize = capsize; 
d.MarkerSize = 8; 
d = errorbar(ElevationAngle_1417(39:41,:), FitCoeff_FLMS195681_O3_Voigt_223K_HP500_PPMM_1417(39:41,:), FitCoeffError_FLMS195681_O3_Voigt_223K_HP500_PPMM_1417(39:41,:), FitCoeffError_FLMS195681_O3_Voigt_223K_HP500_PPMM_1417(39:41,:), 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);
d.MarkerFaceColor = markerfacecolorD; 
d.MarkerEdgeColor = markeredgecolor; 
d.CapSize = capsize; 
d.MarkerSize = 8; 
d = errorbar(ElevationAngle_1417(42:47,:), FitCoeff_FLMS195681_O3_Voigt_223K_HP500_PPMM_1417(42:47,:), FitCoeffError_FLMS195681_O3_Voigt_223K_HP500_PPMM_1417(42:47,:), FitCoeffError_FLMS195681_O3_Voigt_223K_HP500_PPMM_1417(42:47,:), 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);
d.MarkerFaceColor = markerfacecolorD; 
d.MarkerEdgeColor = markeredgecolor; 
d.CapSize = capsize; 
d.MarkerSize = 10; 
% xlim([-95 95]);
xlim([-78 78]);
% Insert formatting 
Fig = gca; Fig.FontSize = 12; set(gcf,'color','w'); Fig.Box = 'on';
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

subplot(2,2,3)
yyaxis left 
ylabel('mEst_R_i_n_g'); hold on 
m = errorbar(mEstLM_AW3_angle_1417(:,1), mEstLM_AW3_mEst_1417(4,1), mEstLM_AW3_stdErr_1417(4,1), mEstLM_AW3_stdErr_1417(4,1), 'vertical', 'LineStyle', 'none', 'Marker', '^', 'Color', markerfacecolor_AW3, 'linewidth', errorbarlinewidthV); hold on; 
m.CapSize = capsize; 
m.MarkerFaceColor = markerfacecolor_AW3; 
m.MarkerEdgeColor = markeredgecolor; 
m = errorbar(mEstLM_AW3_angle_1417(:,5:end), mEstLM_AW3_mEst_1417(4,5:end), mEstLM_AW3_stdErr_1417(4,5:end), mEstLM_AW3_stdErr_1417(4,5:end), 'vertical', 'LineStyle', 'none', 'Marker', '^', 'Color', markerfacecolor_AW3, 'linewidth', errorbarlinewidthV); hold on; 
m.CapSize = capsize; 
m.MarkerFaceColor = markerfacecolor_AW3; 
m.MarkerEdgeColor = markeredgecolor; 
ylim([min(mEstLM_AW3_mEst_1417(4,5:end)) max(mEstLM_AW3_mEst_1417(4,5:end))]);
yyaxis right 
ylabel('Ring fit coefficient'); hold on 
d = errorbar(ElevationAngle_1417(5:47,:), FitCoeff_Ring_1417(5:47,:), FitCoeffError_Ring_1417(5:47,:), FitCoeffError_Ring_1417(5:47,:), 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);hold on 
d.MarkerFaceColor = markerfacecolorD; 
d.MarkerEdgeColor = markeredgecolor; 
d.CapSize = capsize; 
d = errorbar(ElevationAngle_1417(32:34,:), FitCoeff_Ring_1417(32:34,:), FitCoeffError_Ring_1417(32:34,:), FitCoeffError_Ring_1417(32:34,:), 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);
d.MarkerFaceColor = markerfacecolorD; 
d.MarkerEdgeColor = markeredgecolor; 
d.CapSize = capsize; 
d.MarkerSize = 8; 
d = errorbar(ElevationAngle_1417(39:41,:), FitCoeff_Ring_1417(39:41,:), FitCoeffError_Ring_1417(39:41,:), FitCoeffError_Ring_1417(39:41,:), 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);
d.MarkerFaceColor = markerfacecolorD; 
d.MarkerEdgeColor = markeredgecolor; 
d.CapSize = capsize; 
d.MarkerSize = 8; 
d = errorbar(ElevationAngle_1417(42:47,:), FitCoeff_Ring_1417(42:47,:), FitCoeffError_Ring_1417(42:47,:), FitCoeffError_Ring_1417(42:47,:), 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);
d.MarkerFaceColor = markerfacecolorD; 
d.MarkerEdgeColor = markeredgecolor; 
d.CapSize = capsize; 
d.MarkerSize = 10; 
% xlim([-95 95]);
xlim([-78 78]);
% Insert formatting 
Fig = gca; Fig.FontSize = 12; set(gcf,'color','w'); Fig.Box = 'on';
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

subplot(2,2,4)
yyaxis left 
ylabel('mEst_B_s_h_i_f_t'); hold on 
m = errorbar(mEstLM_AW3_angle_1417(:,1), mEstLM_AW3_mEst_1417(5,1), mEstLM_AW3_stdErr_1417(5,1), mEstLM_AW3_stdErr_1417(5,1), 'vertical', 'LineStyle', 'none', 'Marker', '^', 'Color', markerfacecolor_AW3, 'linewidth', errorbarlinewidthV); hold on; 
m.CapSize = capsize; 
m.MarkerFaceColor = markerfacecolor_AW3; 
m.MarkerEdgeColor = markeredgecolor; 
m = errorbar(mEstLM_AW3_angle_1417(:,5:end), mEstLM_AW3_mEst_1417(5,5:end), mEstLM_AW3_stdErr_1417(5,5:end), mEstLM_AW3_stdErr_1417(5,5:end), 'vertical', 'LineStyle', 'none', 'Marker', '^', 'Color', markerfacecolor_AW3, 'linewidth', errorbarlinewidthV); hold on; 
m.CapSize = capsize; 
m.MarkerFaceColor = markerfacecolor_AW3; 
m.MarkerEdgeColor = markeredgecolor; 
ylim([min(mEstLM_AW3_mEst_1417(5,5:end)) max(mEstLM_AW3_mEst_1417(5,5:end))]);

yyaxis right 
ylabel('Shift (nm)'); hold on 
dummy = nan(length(Shift_FLMS195681_SO2_Bogumil_293K_HP500_PPMM_1417), 1); 
d = errorbar(ElevationAngle_1417(5:47,:), Shift_FLMS195681_SO2_Bogumil_293K_HP500_PPMM_1417(5:47,:), dummy(5:47,:), dummy(5:47,:), 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);hold on 
d.MarkerFaceColor = markerfacecolorD; 
d.MarkerEdgeColor = markeredgecolor; 
d.CapSize = capsize; 
d = errorbar(ElevationAngle_1417(32:34,:), Shift_FLMS195681_SO2_Bogumil_293K_HP500_PPMM_1417(32:34,:), dummy(32:34,:), dummy(32:34,:), 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);
d.MarkerFaceColor = markerfacecolorD; 
d.MarkerEdgeColor = markeredgecolor; 
d.CapSize = capsize; 
d.MarkerSize = 8; 
d = errorbar(ElevationAngle_1417(39:41,:), Shift_FLMS195681_SO2_Bogumil_293K_HP500_PPMM_1417(39:41,:), dummy(39:41,:), dummy(39:41,:), 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);
d.MarkerFaceColor = markerfacecolorD; 
d.MarkerEdgeColor = markeredgecolor; 
d.CapSize = capsize; 
d.MarkerSize = 8; 
d = errorbar(ElevationAngle_1417(42:47,:), Shift_FLMS195681_SO2_Bogumil_293K_HP500_PPMM_1417(42:47,:), dummy(42:47,:), dummy(42:47,:), 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);
d.MarkerFaceColor = markerfacecolorD; 
d.MarkerEdgeColor = markeredgecolor; 
d.CapSize = capsize; 
d.MarkerSize = 10; 
% xlim([-95 95]);
xlim([-78 78]);
Fig = gca; Fig.FontSize = 12; set(gcf,'color','w'); Fig.Box = 'on';
xlabel('Scan elevation angle (°)') 
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
% Squeeze_FLMS195681_SO2_Bogumil_293K_HP500_PPMM_1941


%% 1941 
% FIGURE 1 - AW2 with stdErr only 
figure
% AW3 
% 1941 
% subplot(2,2,1)
yyaxis left 
ylim([min(mEstLM_AW3_mEst_1941(1,5:end)) max(mEstLM_AW3_mEst_1941(1,5:end))]);
m = errorbar(mEstLM_AW3_angle_1941(:,1), mEstLM_AW3_mEst_1941(1,1), mEstLM_AW3_stdErr_1941(1,1), mEstLM_AW3_stdErr_1941(1,1), 'vertical', 'LineStyle', 'none', 'Marker', '^', 'Color', markerfacecolor_AW3, 'linewidth', errorbarlinewidthV); hold on; 
m.CapSize = capsize; 
m.MarkerFaceColor = markerfacecolor_AW3; 
m.MarkerEdgeColor = markeredgecolor; 
m = errorbar(mEstLM_AW3_angle_1941(:,5:end), mEstLM_AW3_mEst_1941(1,5:end), mEstLM_AW3_stdErr_1941(1,5:end), mEstLM_AW3_stdErr_1941(1,5:end), 'vertical', 'LineStyle', 'none', 'Marker', '^', 'Color', markerfacecolor_AW3, 'linewidth', errorbarlinewidthV); hold on; 
m.CapSize = capsize; 
m.MarkerFaceColor = markerfacecolor_AW3; 
m.MarkerEdgeColor = markeredgecolor; 
yyaxis right 
d = errorbar(ElevationAngle_1941(5:47,:), FitCoeff_6101_1941(5:47,:), FitCoeffError_6101_1941(5:47,:), FitCoeffError_6101_1941(5:47,:), 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV); hold on 
d.MarkerFaceColor = markerfacecolorD; 
d.MarkerEdgeColor = markeredgecolor; 
d.CapSize = capsize; 
d = errorbar(ElevationAngle_1941(32:34,:), FitCoeff_6101_1941(32:34,:), FitCoeffError_6101_1941(32:34,:), FitCoeffError_6101_1941(32:34,:), 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV); hold on 
d.MarkerFaceColor = markerfacecolorD; 
d.MarkerEdgeColor = markeredgecolor; 
d.CapSize = capsize; 
d.MarkerSize = 8; 
d = errorbar(ElevationAngle_1941(39:41,:), FitCoeff_6101_1941(39:41,:), FitCoeffError_6101_1941(39:41,:), FitCoeffError_6101_1941(39:41,:), 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);
d.MarkerFaceColor = markerfacecolorD; 
d.MarkerEdgeColor = markeredgecolor; 
d.CapSize = capsize; 
d.MarkerSize = 8; 
d = errorbar(ElevationAngle_1941(42:47,:), FitCoeff_6101_1941(42:47,:), FitCoeffError_6101_1941(42:47,:), FitCoeffError_6101_1941(42:47,:), 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);
d.MarkerFaceColor = markerfacecolorD; 
d.MarkerEdgeColor = markeredgecolor; 
d.CapSize = capsize; 
d.MarkerSize = 10; 
ylim([min(FitCoeff_6101_1941) max(FitCoeff_6101_1941)]);
% ax = gca;
% ax.YAxis(1).Color = 'k';
% ax.YAxis(2).Color = 'k';
% xlim([-95 95]);
xlim([-78 78]);
% Insert formatting 
% Fig = gca; Fig.FontSize = 12; set(gcf,'color','w'); Fig.Box = 'on';

% subplot(2,2,2)
yyaxis left 
% ylabel('mEst_O_3 (molec/cm^2)'); hold on 
m = errorbar(mEstLM_AW3_angle_1941(:,1), mEstLM_AW3_mEst_1941(3,1), mEstLM_AW3_stdErr_1941(3,1), mEstLM_AW3_stdErr_1941(3,1), 'vertical', 'LineStyle', 'none', 'Marker', '^', 'Color', markerfacecolor_AW3, 'linewidth', errorbarlinewidthV); hold on; 
m.CapSize = capsize; 
m.MarkerFaceColor = markerfacecolor_AW3; 
m.MarkerEdgeColor = markeredgecolor; 
m = errorbar(mEstLM_AW3_angle_1941(:,5:end), mEstLM_AW3_mEst_1941(3,5:end), mEstLM_AW3_stdErr_1941(3,5:end), mEstLM_AW3_stdErr_1941(3,5:end), 'vertical', 'LineStyle', 'none', 'Marker', '^', 'Color', markerfacecolor_AW3, 'linewidth', errorbarlinewidthV); hold on; 
m.CapSize = capsize; 
m.MarkerFaceColor = markerfacecolor_AW3; 
m.MarkerEdgeColor = markeredgecolor; 
ylim([min(mEstLM_AW3_mEst_1941(3,5:end)) max(mEstLM_AW3_mEst_1941(3,5:end))]);

yyaxis right 
% ylabel('O_3 fit coefficient (molec/cm^2)'); hold on 
d = errorbar(ElevationAngle_1941(5:47,:), FitCoeff_FLMS195681_O3_Voigt_223K_HP500_PPMM_1941(5:47,:), FitCoeffError_FLMS195681_O3_Voigt_223K_HP500_PPMM_1941(5:47,:), FitCoeffError_FLMS195681_O3_Voigt_223K_HP500_PPMM_1941(5:47,:), 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);hold on 
d.MarkerFaceColor = markerfacecolorD; 
d.MarkerEdgeColor = markeredgecolor; 
d.CapSize = capsize; 
d = errorbar(ElevationAngle_1941(32:34,:), FitCoeff_FLMS195681_O3_Voigt_223K_HP500_PPMM_1941(32:34,:), FitCoeffError_FLMS195681_O3_Voigt_223K_HP500_PPMM_1941(32:34,:), FitCoeffError_FLMS195681_O3_Voigt_223K_HP500_PPMM_1941(32:34,:), 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);
d.MarkerFaceColor = markerfacecolorD; 
d.MarkerEdgeColor = markeredgecolor; 
d.CapSize = capsize; 
d.MarkerSize = 8; 
d = errorbar(ElevationAngle_1941(39:41,:), FitCoeff_FLMS195681_O3_Voigt_223K_HP500_PPMM_1941(39:41,:), FitCoeffError_FLMS195681_O3_Voigt_223K_HP500_PPMM_1941(39:41,:), FitCoeffError_FLMS195681_O3_Voigt_223K_HP500_PPMM_1941(39:41,:), 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);
d.MarkerFaceColor = markerfacecolorD; 
d.MarkerEdgeColor = markeredgecolor; 
d.CapSize = capsize; 
d.MarkerSize = 8; 
d = errorbar(ElevationAngle_1941(42:47,:), FitCoeff_FLMS195681_O3_Voigt_223K_HP500_PPMM_1941(42:47,:), FitCoeffError_FLMS195681_O3_Voigt_223K_HP500_PPMM_1941(42:47,:), FitCoeffError_FLMS195681_O3_Voigt_223K_HP500_PPMM_1941(42:47,:), 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);
d.MarkerFaceColor = markerfacecolorD; 
d.MarkerEdgeColor = markeredgecolor; 
d.CapSize = capsize; 
d.MarkerSize = 10; 
% xlim([-95 95]);
xlim([-78 78]);
% Insert formatting 
% Fig = gca; Fig.FontSize = 12; set(gcf,'color','w'); Fig.Box = 'on';
% ax = gca;
% ax.YAxis(1).Color = 'k';
% ax.YAxis(2).Color = 'k';

subplot(2,2,3)
yyaxis left 
% ylabel('mEst_R_i_n_g'); hold on 
m = errorbar(mEstLM_AW3_angle_1941(:,1), mEstLM_AW3_mEst_1941(4,1), mEstLM_AW3_stdErr_1941(4,1), mEstLM_AW3_stdErr_1941(4,1), 'vertical', 'LineStyle', 'none', 'Marker', '^', 'Color', markerfacecolor_AW3, 'linewidth', errorbarlinewidthV); hold on; 
m.CapSize = capsize; 
m.MarkerFaceColor = markerfacecolor_AW3; 
m.MarkerEdgeColor = markeredgecolor; 
m = errorbar(mEstLM_AW3_angle_1941(:,5:end), mEstLM_AW3_mEst_1941(4,5:end), mEstLM_AW3_stdErr_1941(4,5:end), mEstLM_AW3_stdErr_1941(4,5:end), 'vertical', 'LineStyle', 'none', 'Marker', '^', 'Color', markerfacecolor_AW3, 'linewidth', errorbarlinewidthV); hold on; 
m.CapSize = capsize; 
m.MarkerFaceColor = markerfacecolor_AW3; 
m.MarkerEdgeColor = markeredgecolor; 
ylim([min(mEstLM_AW3_mEst_1941(4,5:end)) max(mEstLM_AW3_mEst_1941(4,5:end))]);
yyaxis right 
% ylabel('Ring fit coefficient'); hold on 
d = errorbar(ElevationAngle_1941(5:47,:), FitCoeff_Ring_1941(5:47,:), FitCoeffError_Ring_1941(5:47,:), FitCoeffError_Ring_1941(5:47,:), 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);hold on 
d.MarkerFaceColor = markerfacecolorD; 
d.MarkerEdgeColor = markeredgecolor; 
d.CapSize = capsize; 
d = errorbar(ElevationAngle_1941(32:34,:), FitCoeff_Ring_1941(32:34,:), FitCoeffError_Ring_1941(32:34,:), FitCoeffError_Ring_1941(32:34,:), 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);
d.MarkerFaceColor = markerfacecolorD; 
d.MarkerEdgeColor = markeredgecolor; 
d.CapSize = capsize; 
d.MarkerSize = 8; 
d = errorbar(ElevationAngle_1941(39:41,:), FitCoeff_Ring_1941(39:41,:), FitCoeffError_Ring_1941(39:41,:), FitCoeffError_Ring_1941(39:41,:), 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);
d.MarkerFaceColor = markerfacecolorD; 
d.MarkerEdgeColor = markeredgecolor; 
d.CapSize = capsize; 
d.MarkerSize = 8; 
d = errorbar(ElevationAngle_1941(42:47,:), FitCoeff_Ring_1941(42:47,:), FitCoeffError_Ring_1941(42:47,:), FitCoeffError_Ring_1941(42:47,:), 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);
d.MarkerFaceColor = markerfacecolorD; 
d.MarkerEdgeColor = markeredgecolor; 
d.CapSize = capsize; 
d.MarkerSize = 10; 
% xlim([-95 95]);
xlim([-78 78]);
% Insert formatting 
% Fig = gca; Fig.FontSize = 12; set(gcf,'color','w'); Fig.Box = 'on';
% ax = gca;
% ax.YAxis(1).Color = 'k';
% ax.YAxis(2).Color = 'k';

% subplot(2,2,4)
yyaxis left 
% ylabel('mEst_B_s_h_i_f_t'); hold on 
m = errorbar(mEstLM_AW3_angle_1941(:,1), mEstLM_AW3_mEst_1941(5,1), mEstLM_AW3_stdErr_1941(5,1), mEstLM_AW3_stdErr_1941(5,1), 'vertical', 'LineStyle', 'none', 'Marker', '^', 'Color', markerfacecolor_AW3, 'linewidth', errorbarlinewidthV); hold on; 
m.CapSize = capsize; 
m.MarkerFaceColor = markerfacecolor_AW3; 
m.MarkerEdgeColor = markeredgecolor; 
m = errorbar(mEstLM_AW3_angle_1941(:,5:end), mEstLM_AW3_mEst_1941(5,5:end), mEstLM_AW3_stdErr_1941(5,5:end), mEstLM_AW3_stdErr_1941(5,5:end), 'vertical', 'LineStyle', 'none', 'Marker', '^', 'Color', markerfacecolor_AW3, 'linewidth', errorbarlinewidthV); hold on; 
m.CapSize = capsize; 
m.MarkerFaceColor = markerfacecolor_AW3; 
m.MarkerEdgeColor = markeredgecolor; 
ylim([min(mEstLM_AW3_mEst_1941(5,5:end)) max(mEstLM_AW3_mEst_1941(5,5:end))]);

yyaxis right 
% ylabel('Shift (nm)'); hold on 
dummy = nan(length(Shift_FLMS195681_SO2_Bogumil_293K_HP500_PPMM_1941), 1); 
d = errorbar(ElevationAngle_1941(5:47,:), Shift_FLMS195681_SO2_Bogumil_293K_HP500_PPMM_1941(5:47,:), dummy(5:47,:), dummy(5:47,:), 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);hold on 
d.MarkerFaceColor = markerfacecolorD; 
d.MarkerEdgeColor = markeredgecolor; 
d.CapSize = capsize; 
d = errorbar(ElevationAngle_1941(32:34,:), Shift_FLMS195681_SO2_Bogumil_293K_HP500_PPMM_1941(32:34,:), dummy(32:34,:), dummy(32:34,:), 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);
d.MarkerFaceColor = markerfacecolorD; 
d.MarkerEdgeColor = markeredgecolor; 
d.CapSize = capsize; 
d.MarkerSize = 8; 
d = errorbar(ElevationAngle_1941(39:41,:), Shift_FLMS195681_SO2_Bogumil_293K_HP500_PPMM_1941(39:41,:), dummy(39:41,:), dummy(39:41,:), 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);
d.MarkerFaceColor = markerfacecolorD; 
d.MarkerEdgeColor = markeredgecolor; 
d.CapSize = capsize; 
d.MarkerSize = 8; 
d = errorbar(ElevationAngle_1941(42:47,:), Shift_FLMS195681_SO2_Bogumil_293K_HP500_PPMM_1941(42:47,:), dummy(42:47,:), dummy(42:47,:), 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);
d.MarkerFaceColor = markerfacecolorD; 
d.MarkerEdgeColor = markeredgecolor; 
d.CapSize = capsize; 
d.MarkerSize = 10; 
% xlim([-95 95]);
xlim([-78 78]);
% Fig = gca; Fig.FontSize = 12; set(gcf,'color','w'); Fig.Box = 'on';
% xlabel('Scan elevation angle (°)') 
% ax = gca;
% ax.YAxis(1).Color = 'k';
% ax.YAxis(2).Color = 'k';


























% 
% 
% 
% subplot(1,2,1)
% m = errorbar(mEstLM_AW3_angle_1417, mEstLM_AW3_mEst_1417(2,:), mEstLM_AW3_stdErr_1417(2,:), mEstLM_AW3_stdErr_1417(2,:), 'vertical', 'LineStyle', 'none', 'Marker', '^', 'Color', markerfacecolor_AW3, 'linewidth', errorbarlinewidthV); hold on; 
% m.CapSize = capsize; 
% m.MarkerFaceColor = markerfacecolor_AW3; 
% m.MarkerEdgeColor = markeredgecolor; 
% title('14:17 UTC'); hold on; 
% xlim([-95 95]);
% ylim([-9.9e16 1.5e18]);
% ylabel('SCD SO_2 (molec/cm^2)')
% set(gca,'YMinorTick','on')
% yline(linLim_AW3, ':k')
% Fig = gca; Fig.FontSize = 12; set(gcf,'color','w'); Fig.Box = 'on';
% % 1941
% subplot(1,2,2)
% m = errorbar(mEstLM_AW3_angle_1941, mEstLM_AW3_mEst_1941(2,:), mEstLM_AW3_stdErr_1941(2,:), mEstLM_AW3_stdErr_1941(2,:), 'vertical', 'LineStyle', 'none', 'Marker', '^', 'Color', markerfacecolor_AW3, 'linewidth', errorbarlinewidthV); hold on; 
% m.CapSize = capsize; 
% m.MarkerFaceColor = markerfacecolor_AW3; 
% m.MarkerEdgeColor = markeredgecolor; 
% title('19:41 UTC'); hold on; 
% xlim([-95 95]);
% ylim([-9.9e16 1.5e18]);
% yline(linLim_AW3, ':k')
% set(gca,'YMinorTick','on')
% xlabel('Scan elevation angle (°)') 
% Fig = gca; Fig.FontSize = 12; set(gcf,'color','w'); Fig.Box = 'on';
% %% Make any adjustments 
% % SAVES FIGURE
% % fnOut = 'Fig1_AW3';
% % fname = fullfile(outDirFig, fnOut);
% % saveas(gca, fname)



















% % % Plotting results of reference spectra 
% % 
% % addpath '/Users/charlotteb/Documents/Chapter 2/Reference spectra /DOAS_Eval_CK/'
% % addpath '/Users/charlotteb/Documents/Chapter 2/Reference spectra /Output/'
% % outDir = '/Users/charlotteb/Documents/Chapter 2/Reference spectra /outFiles/';
% % outDirFig = '/Users/charlotteb/Documents/Chapter 2/Reference spectra /';
% % 
% % % Loads DOAS results from CK
% % inDirD = '/Users/charlotteb/Documents/Chapter 2/Reference spectra /DOAS_Eval_CK/';
% % inFileD = 'mDOASResults_2022-02-02T11-47-08_310-340nm_VARIABLE.csv'; 
% % inD = fullfile(inDirD, inFileD); 
% % 
% % hlinesD = 1; 
% % 
% % % Loads all data in results file
% % [SpectrumID, FitCoeff_6101,	FitCoeffError_6101,	Shift_6101,	Squeeze_6101, FitCoeff_Ring, FitCoeffError_Ring,	Shift_Ring,...
% %     Squeeze_Ring, FitCoeff_FLMS195681_O3_Voigt_223K_HP500_PPMM,	FitCoeffError_FLMS195681_O3_Voigt_223K_HP500_PPMM, Shift_FLMS195681_O3_Voigt_223K_HP500_PPMM,...
% %     Squeeze_FLMS195681_O3_Voigt_223K_HP500_PPMM, FitCoeff_FLMS195681_SO2_Bogumil_293K_HP500_MOLECCM2, FitCoeff_FLMS195681_SO2_Bogumil_293K_HP500_PPMM,...
% %     FitCoeffError_FLMS195681_SO2_Bogumil_293K_HP500_PPMM, SO2_Error_PPMM, Shift_FLMS195681_SO2_Bogumil_293K_HP500_PPMM,	Squeeze_FLMS195681_SO2_Bogumil_293K_HP500_PPMM,...
% %     WavelengthRangeLow,	WavelengthRangeHigh, FitChi2, TimeStamp, Latitude, Longitude, Altitude,	Speed, Course, GPSWarnCode,	GPSQuality,	ElevationAngle,	AzimuthAngle,...
% %     ExposureTime, Exposures, MaxIntensity, MaxIntensityFitRange, FileName, SpectrometerType, SpectrometerSerialNumber,...
% %     SpectrometerChannel, Remark] = textread(inD,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%s%s%s%s%f%s', 'headerlines', hlinesD, 'delimiter', ',');  
% %  
% % % Isolates SO2 FIT COEFFICIENT and FIT ERROR 
% % SO2D = FitCoeff_FLMS195681_SO2_Bogumil_293K_HP500_MOLECCM2;  
% % SO2ErrD = FitCoeffError_FLMS195681_SO2_Bogumil_293K_HP500_PPMM; 
% % 
% % % % Extracts known SO2 concentration from 'remark'
% % % findStr = "_FLMS195681";
% % % tempStr = string(extractBefore(Remark, findStr)); 
% % % findStr = "_";
% % % elevAnD = str2double(extractBefore(tempStr, findStr)); % Elevation angle of scan 
% % % cellConD = str2double(extractAfter(tempStr, findStr)); % Gas cell concentration 
% % % 
% % % % Gas cell error is +/- 10%
% % % % Conversion from ppmm to molec/cm^2 for expected trend line
% % % conv = 2500000000000000; % copied form DOAS results table 
% % % cellConDconv = cellConD * conv; 
% % % cellConDErr = (cellConD/100) * 10; 
% % % cellConDErrconv = cellConDErr * conv;
% % % 
% % % % Loads results files form linear model 
% % % % MODEL VERSION 1 
% % % inDirLM_AW1 = '/Users/charlotteb/Documents/Chapter 2/Reference spectra /Output/Model version 1 AW1cF/ReRun'; % Later add ReRun; % MODEL VERSION 1 IS AW1 cropped frequencies
% % % inFileLM_AW1_1417 = 'ReferenceSpectra_AW1cF_mEst_manualsaveall.mat'; 
% % % inLM_AW1_1417 = fullfile(inDirLM_AW1, inFileLM_AW1_1417); 
% % % % mEst 
% % % mEstLM_AW1 = load(inLM_AW1_1417,'Results_mEst_real', 'cellCon'); 
% % % cellConLM = mEstLM_AW1.cellCon;  
% % % mEstLM_AW1 = mEstLM_AW1.Results_mEst_real; 
% % % mEstLM_v1_realimag = load(inLM_AW1_1417,'Results_mEst_realimag'); 
% % % mEstLM_v1_realimag = mEstLM_v1_realimag.Results_mEst_realimag; 
% % % % Standard error 
% % % stdErrLM_v1_real = load(inLM_AW1_1417,'Results_stdErr_real'); 
% % % stdErrLM_v1_real = stdErrLM_v1_real.Results_stdErr_real; 
% % % stdErrLM_v1_realimag = load(inLM_AW1_1417,'Results_stdErr_realimag'); 
% % % stdErrLM_v1_realimag = stdErrLM_v1_realimag.Results_stdErr_realimag; 
% % % 
% % % % MODEL VERSION 2
% % % inDirLM_v2 = '/Users/charlotteb/Documents/Chapter 2/Reference spectra /Output/Model version 4 AW2cF bonus high SO2/'; % MODEL VERSION 2 IS AW2 (bonus)
% % % inFileLM_v2 = 'ReferenceSpectra_AW2cF_mEst_manualsaveall.mat'; 
% % % inLM_v2 = fullfile(inDirLM_v2, inFileLM_v2); 
% % % % mEst 
% % % mEstLM_v2_real = load(inLM_v2,'Results_mEst_real'); 
% % % mEstLM_v2_real = mEstLM_v2_real.Results_mEst_real; 
% % % mEstLM_v2_realimag = load(inLM_v2,'Results_mEst_realimag'); 
% % % mEstLM_v2_realimag = mEstLM_v2_realimag.Results_mEst_realimag; 
% % % % Standard error 
% % % stdErrLM_v2_real = load(inLM_v2,'Results_stdErr_real'); 
% % % stdErrLM_v2_real = stdErrLM_v2_real.Results_stdErr_real; 
% % % stdErrLM_v2_realimag = load(inLM_v2,'Results_stdErr_realimag'); 
% % % stdErrLM_v2_realimag = stdErrLM_v2_realimag.Results_stdErr_realimag; 
% % % 
% % % % MODEL VERSION 3
% % % inDirLM_v3 = '/Users/charlotteb/Documents/Chapter 2/Reference spectra /Output/Model version 2 AW3cF/ReRun/'; % MODEL VERSION 3 IS AW3 
% % % inFileLM_v3 = 'ReferenceSpectra_AW3cF_mEst.mat'; 
% % % inLM_v3 = fullfile(inDirLM_v3, inFileLM_v3); 
% % % % mEst 
% % % mEstLM_v3_real = load(inLM_v3,'Results_mEst_real'); 
% % % mEstLM_v3_real = mEstLM_v3_real.Results_mEst_real; 
% % % mEstLM_v3_realimag = load(inLM_v3,'Results_mEst_realimag'); 
% % % mEstLM_v3_realimag = mEstLM_v3_realimag.Results_mEst_realimag; 
% % % % Standard error 
% % % stdErrLM_v3_real = load(inLM_v3,'Results_stdErr_real'); 
% % % stdErrLM_v3_real = stdErrLM_v3_real.Results_stdErr_real; 
% % % stdErrLM_v3_realimag = load(inLM_v3,'Results_stdErr_realimag'); 
% % % stdErrLM_v3_realimag = stdErrLM_v3_realimag.Results_stdErr_realimag; 
% % % 
% % % % MODEL VERSION 4 
% % % inDirLM_v4 = '/Users/charlotteb/Documents/Chapter 2/Reference spectra /Output/Model version 3 AW1aF/'; % MODEL VERSION 4 IS AW1 all frequencies 
% % % inFileLM_v4 = 'ReferenceSpectra_AW1aF_mEst_manualsaveall.mat'; 
% % % inLM_v4 = fullfile(inDirLM_v4, inFileLM_v4); 
% % % % mEst 
% % % mEstLM_v4_real = load(inLM_v4,'Results_mEst_real'); 
% % % mEstLM_v4_real = mEstLM_v4_real.Results_mEst_real; 
% % % mEstLM_v4_realimag = load(inLM_v4,'Results_mEst_realimag'); 
% % % mEstLM_v4_realimag = mEstLM_v4_realimag.Results_mEst_realimag; 
% % % % Standard error 
% % % stdErrLM_v4_real = load(inLM_v4,'Results_stdErr_real'); 
% % % stdErrLM_v4_real = stdErrLM_v4_real.Results_stdErr_real; 
% % % stdErrLM_v4_realimag = load(inLM_v4,'Results_stdErr_realimag'); 
% % % stdErrLM_v4_realimag = stdErrLM_v4_realimag.Results_stdErr_realimag; 
% % % 
% % % % Error calculation using linear model data 
% % % % Gas cell error is +/- 10%
% % % % Conversion from ppmm to molec/cm^2 for expected trend line
% % % cellConLMconv = cellConLM * conv; 
% % % cellConLMErr = (cellConLM/100) * 10; 
% % % cellConLMErrconv = cellConLMErr * conv;
% % % % Adds bounds to expected trend line to relating to gas cell error  
% % % maxExp = cellConLMconv + cellConLMErrconv; 
% % % minExp = cellConLMconv - cellConLMErrconv; 
% % % 
% % % % Variables for formatting plot 
% % % markeredgecolor = 'black'; 
% % % markerfacecolorD =  'white';
% % % markerfacecolor_AW1 = [0.11 0.24 0.81]; % AW1 
% % % markerfacecolor_v1cm = [0.62 0.68 0.95];
% % % markerfacecolor_v4 = [0.11 0.24 0.81]; % AW1
% % % markerfacecolor_AW2 = [0.91 0.05 0.05]; % AW2
% % % markerfacecolor_AW3 = [0.09 0.54 0.03]; % AW3
% % % markersize = 6; 
% % % capsize = 0; 
% % % errorbarcolor = [0.92,0.92,0.92]; % Horizontal 
% % % errorbarlinewidthH = 0.1; 
% % % errorbarlinewidthV = 0.5; 
% % % minX = -70; 
% % % maxX = 2330; 
% % % trendlinewidthmain = 0.1; 
% % % % trendlinewidthlimits = 0.1; 
% % % % Axis limits 
% % % minY = min(cellConD)-1e18;
% % % maxY = 8.2e18; % max(cellConDconv)+3.5e18;
% % % 
% % % % Plots SO2 fit results from DOAS analysis 
% % % figure
% % % subplot(2,2,1) % MODEL VERSION 1 
% % % % Gas cell error bars 
% % % m = errorbar(cellConD, SO2D, cellConDErr, cellConDErr, 'horizontal', 'LineStyle', 'none', 'Marker', 'none', 'Color', errorbarcolor, 'linewidth', errorbarlinewidthH); hold on; 
% % % m.CapSize = capsize; 
% % % % Expected 
% % % et = plot(cellConLM, cellConLMconv, ':k'); et.LineWidth = trendlinewidthmain; hold on; 
% % % SO2LMr = mEstLM_AW1(2,:); % SO2 is G2 (I0, SO2, O3, Ring, Bshift) 
% % % SO2LMri = mEstLM_v1_realimag(2,:); % SO2 is G2 (I0, SO2, O3, Ring, Bshift) 
% % % SO2ErrLMr = stdErrLM_v1_real(2,:); 
% % % SO2ErrLMri = stdErrLM_v1_realimag(2,:); 
% % % % Plots mEst 
% % % m = errorbar(cellConLM, SO2LMr, cellConLMErr, cellConLMErr, 'horizontal', 'LineStyle', 'none', 'Marker', 'none', 'Color', errorbarcolor, 'linewidth', errorbarlinewidthV); hold on; 
% % % m.CapSize = capsize; 
% % % m = errorbar(cellConLM, SO2LMri, cellConLMErr, cellConLMErr, 'horizontal', 'LineStyle', 'none', 'Marker', 'none', 'Color', errorbarcolor, 'linewidth', errorbarlinewidthV); hold on; 
% % % m.CapSize = capsize; 
% % % % etM = plot(cellConLM, maxExp, '-k'); etM.LineWidth = trendlinewidthlimits;
% % % % etm = plot(cellConLM, minExp, '-k'); etm.LineWidth = trendlinewidthlimits;
% % % % DOAS 
% % % d = errorbar(cellConD, SO2D, SO2ErrD, SO2ErrD, 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);
% % % d.MarkerFaceColor = markerfacecolorD; 
% % % d.MarkerEdgeColor = markeredgecolor; 
% % % d.CapSize = capsize; 
% % % xlim([minX maxX]);
% % % ylim([minY maxY]);
% % % % ylabel('(d)SCD SO_2 (molec/cm^2)')
% % % % xlabel('Gas cell concentration (ppm·m)')
% % % Fig = gca; Fig.FontSize = 14; set(gcf,'color','w'); Fig.Box = 'on'; 
% % % % mEst and standard error 
% % % % Real and complex
% % % % SO2LMr = mEstLM_v1_real(2,:); % SO2 is G2 (I0, SO2, O3, Ring, Bshift) 
% % % % SO2LMri = mEstLM_v1_realimag(2,:); % SO2 is G2 (I0, SO2, O3, Ring, Bshift) 
% % % % SO2ErrLMr = stdErrLM_v1_real(2,:); 
% % % % SO2ErrLMri = stdErrLM_v1_realimag(2,:); 
% % % % % Plots mEst 
% % % % e = errorbar(cellConLM, SO2LMr, cellConLMErr, cellConLMErr, 'horizontal', 'LineStyle', 'none', 'Marker', 'none', 'Color', errorbarcolor, 'linewidth', errorbarlinewidthV); hold on; 
% % % % e.CapSize = capsize; 
% % % % e = errorbar(cellConLM, SO2LMri, cellConLMErr, cellConLMErr, 'horizontal', 'LineStyle', 'none', 'Marker', 'none', 'Color', errorbarcolor, 'linewidth', errorbarlinewidthV); hold on; 
% % % e.CapSize = capsize; 
% % % Real
% % lr = errorbar(cellConLM, SO2LMr, SO2ErrLMr, SO2ErrLMr, 'vertical', 'LineStyle', 'none', 'Marker', '^', 'Color', markerfacecolor_AW1, 'linewidth', errorbarlinewidthV);
% % lr.MarkerFaceColor = markerfacecolor_AW1; 
% % lr.MarkerEdgeColor = markeredgecolor; 
% % lr.CapSize = capsize; 
% % % Real imaginary 
% % lri = errorbar(cellConLM, SO2LMri, SO2ErrLMri, SO2ErrLMri, 'vertical', 'LineStyle', 'none', 'Marker', 'v', 'Color', markerfacecolor_v1cm, 'linewidth', errorbarlinewidthV);
% % lri.MarkerFaceColor = markerfacecolor_v1cm; 
% % lri.MarkerEdgeColor = markeredgecolor;  % markeredgecolor; 
% % lri.CapSize = capsize;
% % % Offset 
% % nogascellInd = find(cellConLM == 0); 
% % nogascell = cellConLM(nogascellInd); 
% % nogascellmEstr = SO2LMr(nogascellInd); 
% % nogascellmEstri = SO2LMri(nogascellInd); 
% % % figure 
% % % plot(nogascell, nogascellmEstr', 'ok') 
% % disp("Version 1")
% % disp("check negative")
% % max(nogascellmEstr) 
% % max(nogascellmEstri) 
% % disp("find mean offset")
% % meanNogascellmEstr = mean(nogascellmEstr)
% % meanNogascellmEstri = mean(nogascellmEstri)
% % 
% % % MODEL VERSION 2
% % subplot(2,2,2) 
% % % Gas cell error bars 
% % m = errorbar(cellConD, SO2D, cellConDErr, cellConDErr, 'horizontal', 'LineStyle', 'none', 'Marker', 'none', 'Color', errorbarcolor, 'linewidth', errorbarlinewidthH); hold on; 
% % m.CapSize = capsize; 
% % % Expected 
% % et = plot(cellConLM, cellConLMconv, ':k'); et.LineWidth = trendlinewidthmain; hold on; 
% % % etM = plot(cellConLM, maxExp, '-k'); etM.LineWidth = trendlinewidthlimits;
% % % etm = plot(cellConLM, minExp, '-k'); etm.LineWidth = trendlinewidthlimits;
% % SO2LMr = mEstLM_v2_real(2,:); % SO2 is G2 (I0, SO2, O3, Ring, Bshift) 
% % % SO2LMri = mEstLM_v2_realimag(2,:); % SO2 is G2 (I0, SO2, O3, Ring, Bshift) 
% % SO2ErrLMr = stdErrLM_v2_real(2,:); 
% % % SO2ErrLMri = stdErrLM_v2_realimag(2,:); 
% % % SO2ErrLMr = NaN(size(SO2LMr)); % Placeholder
% % % SO2ErrLMri = NaN(size(SO2LMri)); % Placeholder
% % % Plots mEst 
% % m = errorbar(cellConLM, SO2LMr, cellConLMErr, cellConLMErr, 'horizontal', 'LineStyle', 'none', 'Marker', 'none', 'Color', errorbarcolor, 'linewidth', errorbarlinewidthV); hold on; 
% % m.CapSize = capsize; 
% % % DOAS 
% % d = errorbar(cellConD, SO2D, SO2ErrD, SO2ErrD, 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);
% % d.MarkerFaceColor = markerfacecolorD; 
% % d.MarkerEdgeColor = markeredgecolor; 
% % d.CapSize = capsize; 
% % xlim([minX maxX]);
% % ylim([minY maxY]);
% % % ylabel('(d)SCD SO_2 (molec/cm^2)')
% % % xlabel('Gas cell concentration (ppm·m)')
% % Fig = gca; Fig.FontSize = 14; set(gcf,'color','w'); Fig.Box = 'on'; 
% % % mEst and standard error 
% % % Real and complex
% % % SO2LMr = mEstLM_v2_real(2,:); % SO2 is G2 (I0, SO2, O3, Ring, Bshift) 
% % % % SO2LMri = mEstLM_v2_realimag(2,:); % SO2 is G2 (I0, SO2, O3, Ring, Bshift) 
% % % SO2ErrLMr = stdErrLM_v2_real(2,:); 
% % % % SO2ErrLMri = stdErrLM_v2_realimag(2,:); 
% % % % SO2ErrLMr = NaN(size(SO2LMr)); % Placeholder
% % % % SO2ErrLMri = NaN(size(SO2LMri)); % Placeholder
% % % % Plots mEst 
% % % e = errorbar(cellConLM, SO2LMr, cellConLMErr, cellConLMErr, 'horizontal', 'LineStyle', 'none', 'Marker', 'none', 'Color', errorbarcolor, 'linewidth', errorbarlinewidthV); hold on; 
% % % e.CapSize = capsize; 
% % % e = errorbar(cellConLM, SO2LMri, cellConLMErr, cellConLMErr, 'horizontal', 'LineStyle', 'none', 'Marker', 'none', 'Color', errorbarcolor, 'linewidth', errorbarlinewidthV); hold on; 
% % % e.CapSize = capsize; 
% % % Real
% % lr = errorbar(cellConLM, SO2LMr, SO2ErrLMr, SO2ErrLMr, 'vertical', 'LineStyle', 'none', 'Marker', '^', 'Color', markerfacecolor_AW2, 'linewidth', errorbarlinewidthV);
% % lr.MarkerFaceColor = markerfacecolor_AW2; 
% % lr.MarkerEdgeColor = markeredgecolor;  % markeredgecolor; 
% % lr.CapSize = capsize; 
% % % % Real imaginary 
% % % lri = errorbar(cellConLM, SO2LMri, SO2ErrLMri, SO2ErrLMri, 'vertical', 'LineStyle', 'none', 'Marker', 'v', 'Color', markerfacecolor_v2, 'linewidth', errorbarlinewidthV);
% % % lri.MarkerFaceColor = markerfacecolorLMAW3; 
% % % lri.MarkerEdgeColor = 'none';  % markeredgecolor; 
% % % lri.CapSize = capsize;
% % % Offset 
% % nogascellInd = find(cellConLM == 0); 
% % nogascell = cellConLM(nogascellInd); 
% % nogascellmEstr = SO2LMr(nogascellInd); 
% % % nogascellmEstri = SO2LMri(nogascellInd); 
% % % figure 
% % % plot(nogascell, nogascellmEstr', 'ok') 
% % disp("Version 2")
% % disp("check negative")
% % max(nogascellmEstr) 
% % % max(nogascellmEstri) 
% % disp("find mean offset")
% % meanNogascellmEstr = mean(nogascellmEstr)
% % % meanNogascellmEstri = mean(nogascellmEstri)
% % 
% % % MODEL VERSION 3
% % subplot(2,2,3) 
% % % Gas cell error bars 
% % m = errorbar(cellConD, SO2D, cellConDErr, cellConDErr, 'horizontal', 'LineStyle', 'none', 'Marker', 'none', 'Color', errorbarcolor, 'linewidth', errorbarlinewidthH); hold on; 
% % m.CapSize = capsize; 
% % % Expected 
% % et = plot(cellConLM, cellConLMconv, ':k'); et.LineWidth = trendlinewidthmain; hold on; 
% % % etM = plot(cellConLM, maxExp, '-k'); etM.LineWidth = trendlinewidthlimits;
% % % etm = plot(cellConLM, minExp, '-k'); etm.LineWidth = trendlinewidthlimits;
% % % mEst and standard error 
% % % Real and complex
% % SO2LMr = mEstLM_v3_real(2,:); % SO2 is G2 (I0, SO2, O3, Ring, Bshift) 
% % % SO2LMri = mEstLM_v3_realimag(2,:); % SO2 is G2 (I0, SO2, O3, Ring, Bshift) 
% % SO2ErrLMr = stdErrLM_v3_real(2,:); 
% % % SO2ErrLMri = stdErrLM_v3_realimag(2,:); 
% % % Plots mEst 
% % m = errorbar(cellConLM, SO2LMr, cellConLMErr, cellConLMErr, 'horizontal', 'LineStyle', 'none', 'Marker', 'none', 'Color', errorbarcolor, 'linewidth', errorbarlinewidthV); hold on; 
% % m.CapSize = capsize; 
% % % DOAS 
% % d = errorbar(cellConD, SO2D, SO2ErrD, SO2ErrD, 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);
% % d.MarkerFaceColor = markerfacecolorD; 
% % d.MarkerEdgeColor = markeredgecolor; 
% d.CapSize = capsize; 
% xlim([minX maxX]);
% ylim([minY maxY]);
% ylabel('(d)SCD SO_2 (molec/cm^2)')
% % xlabel('Gas cell concentration (ppm·m)')
% Fig = gca; Fig.FontSize = 14; set(gcf,'color','w'); Fig.Box = 'on';
% % % mEst and standard error 
% % % Real and complex
% % SO2LMr = mEstLM_v3_real(2,:); % SO2 is G2 (I0, SO2, O3, Ring, Bshift) 
% % % SO2LMri = mEstLM_v3_realimag(2,:); % SO2 is G2 (I0, SO2, O3, Ring, Bshift) 
% % SO2ErrLMr = stdErrLM_v3_real(2,:); 
% % % SO2ErrLMri = stdErrLM_v3_realimag(2,:); 
% % % Plots mEst 
% % e = errorbar(cellConLM, SO2LMr, cellConLMErr, cellConLMErr, 'horizontal', 'LineStyle', 'none', 'Marker', 'none', 'Color', errorbarcolor, 'linewidth', errorbarlinewidthV); hold on; 
% % e.CapSize = capsize; 
% % Real
% lr = errorbar(cellConLM, SO2LMr, SO2ErrLMr, SO2ErrLMr, 'vertical', 'LineStyle', 'none', 'Marker', '^', 'Color', markerfacecolor_AW3, 'linewidth', errorbarlinewidthV);
% lr.MarkerFaceColor = markerfacecolor_AW3; 
% lr.MarkerEdgeColor = markeredgecolor; 
% lr.CapSize = capsize; 
% % Real imaginary 
% % lri = errorbar(cellConLM, SO2LMri, SO2ErrLMri, SO2ErrLMri, 'vertical', 'LineStyle', 'none', 'Marker', 'v', 'Color', markerfacecolor_v3, 'linewidth', errorbarlinewidthV);
% % lri.MarkerFaceColor = 'none'; 
% % lri.MarkerEdgeColor = markerfacecolorLMAW1; 
% % lri.CapSize = capsize;
% % Offset 
% nogascellInd = find(cellConLM == 0); 
% nogascell = cellConLM(nogascellInd); 
% nogascellmEstr = SO2LMr(nogascellInd); 
% % nogascellmEstri = SO2LMri(nogascellInd); 
% % figure 
% % plot(nogascell, nogascellmEstr', 'ok') 
% disp("Version 3")
% disp("check negative")
% max(nogascellmEstr) 
% % max(nogascellmEstri) 
% disp("find mean offset")
% meanNogascellmEstr = mean(nogascellmEstr)
% % meanNogascellmEstri = mean(nogascellmEstri)
% 
% % MODEL VERSION 4
% subplot(2,2,4) 
% % Gas cell error bars 
% m = errorbar(cellConD, SO2D, cellConDErr, cellConDErr, 'horizontal', 'LineStyle', 'none', 'Marker', 'none', 'Color', errorbarcolor, 'linewidth', errorbarlinewidthH); hold on; 
% m.CapSize = capsize; 
% % Expected 
% et = plot(cellConLM, cellConLMconv, ':k'); et.LineWidth = trendlinewidthmain; hold on; 
% % etM = plot(cellConLM, maxExp, '-k'); etM.LineWidth = trendlinewidthlimits;
% % etm = plot(cellConLM, minExp, '-k'); etm.LineWidth = trendlinewidthlimits;
% % mEst and standard error 
% % Real and complex
% SO2LMr = mEstLM_v4_real(2,:); % SO2 is G2 (I0, SO2, O3, Ring, Bshift) 
% % SO2LMri = mEstLM_v4_realimag(2,:); % SO2 is G2 (I0, SO2, O3, Ring, Bshift) 
% SO2ErrLMr = stdErrLM_v4_real(2,:); 
% % SO2ErrLMri = stdErrLM_v4_realimag(2,:); 
% % Plots mEst 
% m = errorbar(cellConLM, SO2LMr, cellConLMErr, cellConLMErr, 'horizontal', 'LineStyle', 'none', 'Marker', 'none', 'Color', errorbarcolor, 'linewidth', errorbarlinewidthV); hold on; 
% m.CapSize = capsize; 
% m = errorbar(cellConLM, SO2LMri, cellConLMErr, cellConLMErr, 'horizontal', 'LineStyle', 'none', 'Marker', 'none', 'Color', errorbarcolor, 'linewidth', errorbarlinewidthV); hold on; 
% m.CapSize = capsize; 
% % DOAS 
% d = errorbar(cellConD, SO2D, SO2ErrD, SO2ErrD, 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);
% d.MarkerFaceColor = markerfacecolorD; 
% d.MarkerEdgeColor = markeredgecolor; 
% d.CapSize = capsize; 
% xlim([minX maxX]);
% ylim([minY maxY]);
% % ylabel('(d)SCD SO_2 (molec/cm^2)')
% xlabel('Gas cell concentration (ppm·m)')
% Fig = gca; Fig.FontSize = 14; set(gcf,'color','w'); Fig.Box = 'on';
% % % mEst and standard error 
% % % Real and complex
% % SO2LMr = mEstLM_v4_real(2,:); % SO2 is G2 (I0, SO2, O3, Ring, Bshift) 
% % % SO2LMri = mEstLM_v4_realimag(2,:); % SO2 is G2 (I0, SO2, O3, Ring, Bshift) 
% % SO2ErrLMr = stdErrLM_v4_real(2,:); 
% % % SO2ErrLMri = stdErrLM_v4_realimag(2,:); 
% % % Plots mEst 
% % e = errorbar(cellConLM, SO2LMr, cellConLMErr, cellConLMErr, 'horizontal', 'LineStyle', 'none', 'Marker', 'none', 'Color', errorbarcolor, 'linewidth', errorbarlinewidthV); hold on; 
% % e.CapSize = capsize; 
% % e = errorbar(cellConLM, SO2LMri, cellConLMErr, cellConLMErr, 'horizontal', 'LineStyle', 'none', 'Marker', 'none', 'Color', errorbarcolor, 'linewidth', errorbarlinewidthV); hold on; 
% % e.CapSize = capsize; 
% % Real
% lr = errorbar(cellConLM, SO2LMr, SO2ErrLMr, SO2ErrLMr, 'vertical', 'LineStyle', 'none', 'Marker', '^', 'Color', markerfacecolor_AW1, 'linewidth', errorbarlinewidthV);
% lr.MarkerFaceColor = markerfacecolor_AW1; 
% lr.MarkerEdgeColor = markeredgecolor; 
% lr.CapSize = capsize; 
% % Real imaginary 
% lri = errorbar(cellConLM, SO2LMri, SO2ErrLMri, SO2ErrLMri, 'vertical', 'LineStyle', 'none', 'Marker', 'v', 'Color', markerfacecolor_v1cm, 'linewidth', errorbarlinewidthV);
% lri.MarkerFaceColor = markerfacecolor_v1cm; 
% lri.MarkerEdgeColor = markeredgecolor;  % markeredgecolor; 
% lri.CapSize = capsize;
% % Offset 
% nogascellInd = find(cellConLM == 0); 
% nogascell = cellConLM(nogascellInd); 
% nogascellmEstr = SO2LMr(nogascellInd); 
% nogascellmEstri = SO2LMri(nogascellInd); 
% % figure 
% % plot(nogascell, nogascellmEstr', 'ok') 
% disp("Version 4")
% disp("check negative")
% max(nogascellmEstr) 
% max(nogascellmEstri) 
% disp("find mean offset")
% meanNogascellmEstr = mean(nogascellmEstr)
% meanNogascellmEstri = mean(nogascellmEstri)
% % Save 
% fname = fullfile(outDir, 'results');
% save(fname); % Saves all variables  
% %% Make any adjustments 
% % SAVES FIGURE
% % fnOut = 'Results';
% % fname = fullfile(outDirFig, fnOut);
% % saveas(gca, fname)
% 
% % % Check complex magnitude of AW2 
% % % MODEL VERSION 2
% % figure
% % % Gas cell error bars 
% % e = errorbar(cellConD, SO2D, cellConDErr, cellConDErr, 'horizontal', 'LineStyle', 'none', 'Marker', 'none', 'Color', errorbarcolor, 'linewidth', errorbarlinewidthH); hold on; 
% % e.CapSize = capsize; 
% % % Expected 
% % et = plot(cellConLM, cellConLMconv, ':k'); et.LineWidth = trendlinewidthmain; hold on; 
% % % etM = plot(cellConLM, maxExp, '-k'); etM.LineWidth = trendlinewidthlimits;
% % % etm = plot(cellConLM, minExp, '-k'); etm.LineWidth = trendlinewidthlimits;
% % SO2LMr = mEstLM_v2_real(2,:); % SO2 is G2 (I0, SO2, O3, Ring, Bshift) 
% % SO2LMri = mEstLM_v2_realimag(2,:); % SO2 is G2 (I0, SO2, O3, Ring, Bshift) 
% % SO2ErrLMr = stdErrLM_v2_real(2,:); 
% % SO2ErrLMri = stdErrLM_v2_realimag(2,:); 
% % % SO2ErrLMr = NaN(size(SO2LMr)); % Placeholder
% % % SO2ErrLMri = NaN(size(SO2LMri)); % Placeholder
% % % Plots mEst 
% % e = errorbar(cellConLM, SO2LMr, cellConLMErr, cellConLMErr, 'horizontal', 'LineStyle', 'none', 'Marker', 'none', 'Color', errorbarcolor, 'linewidth', errorbarlinewidthV); hold on; 
% % e.CapSize = capsize; 
% % % DOAS 
% % d = errorbar(cellConD, SO2D, SO2ErrD, SO2ErrD, 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);
% % d.MarkerFaceColor = markerfacecolorD; 
% % d.MarkerEdgeColor = markeredgecolor; 
% % d.CapSize = capsize; 
% % xlim([minX maxX]);
% % ylim([minY maxY]);
% % % ylabel('(d)SCD SO_2 (molec/cm^2)')
% % % xlabel('Gas cell concentration (ppm·m)')
% % Fig = gca; Fig.FontSize = 14; set(gcf,'color','w'); Fig.Box = 'on'; 
% % % mEst and standard error 
% % % Real and complex
% % % SO2LMr = mEstLM_v2_real(2,:); % SO2 is G2 (I0, SO2, O3, Ring, Bshift) 
% % % % SO2LMri = mEstLM_v2_realimag(2,:); % SO2 is G2 (I0, SO2, O3, Ring, Bshift) 
% % % SO2ErrLMr = stdErrLM_v2_real(2,:); 
% % % % SO2ErrLMri = stdErrLM_v2_realimag(2,:); 
% % % % SO2ErrLMr = NaN(size(SO2LMr)); % Placeholder
% % % % SO2ErrLMri = NaN(size(SO2LMri)); % Placeholder
% % % % Plots mEst 
% % % e = errorbar(cellConLM, SO2LMr, cellConLMErr, cellConLMErr, 'horizontal', 'LineStyle', 'none', 'Marker', 'none', 'Color', errorbarcolor, 'linewidth', errorbarlinewidthV); hold on; 
% % % e.CapSize = capsize; 
% % % e = errorbar(cellConLM, SO2LMri, cellConLMErr, cellConLMErr, 'horizontal', 'LineStyle', 'none', 'Marker', 'none', 'Color', errorbarcolor, 'linewidth', errorbarlinewidthV); hold on; 
% % % e.CapSize = capsize; 
% % % Real
% % % lr = errorbar(cellConLM, SO2LMr, SO2ErrLMr, SO2ErrLMr, 'vertical', 'LineStyle', 'none', 'Marker', '^', 'Color', markerfacecolor_v2, 'linewidth', errorbarlinewidthV);
% % % lr.MarkerFaceColor = markerfacecolor_v2; 
% % % lr.MarkerEdgeColor = markeredgecolor;  % markeredgecolor; 
% % % lr.CapSize = capsize; 
% % % % Real imaginary 
% % lri = errorbar(cellConLM, SO2LMri, SO2ErrLMri, SO2ErrLMri, 'vertical', 'LineStyle', 'none', 'Marker', 'v', 'Color', markerfacecolor_v2, 'linewidth', errorbarlinewidthV);
% % lri.MarkerFaceColor = markerfacecolorLMAW3; 
% % lri.MarkerEdgeColor = 'none';  % markeredgecolor; 
% % lri.CapSize = capsize;
% 
% % %% Figure 19 - mEst by elevation angle 
% % offset = 2.7e17; % Table 3 
% % % mEst and standard error 
% % % Real and complex
% % SO2LMr = mEstLM_v3_real(2,:); % SO2 is G2 (I0, SO2, O3, Ring, Bshift) % MODEL VERSION 3
% % % SO2LMri = mEstLM_v3_realimag(2,:); % SO2 is G2 (I0, SO2, O3, Ring, Bshift) 
% % SO2ErrLMr = stdErrLM_v3_real(2,:); 
% % % SO2ErrLMri = stdErrLM_v3_realimag(2,:); 
% % SO2LMrOffset = SO2LMr+offset;  
% % %% Plot by elevation angle of scan 
% % unAngle = unique(elevAnD); % Angles taken from spectra file name and refer to angle from horizon. So 90 is zenith and 30 is close to horizon. 
% % unAngle = flip(unAngle); % Reverse so zenith is first 
% % % measpAngle = 1000; 
% % % doasSO2 = nan(measpAngle, length(unAngle)); 
% % % doasSO2Err = nan(measpAngle, length(unAngle)); 
% % % mEstSO2 = nan(measpAngle, length(unAngle)); 
% % % mEstSO2Err = nan(measpAngle, length(unAngle)); 
% % % Organsise data by elevation angle 
% % figure
% % markersize = 4; 
% % markerfacecolor_v3 = [0.09 0.54 0.03]; % AW3
% % errorbarlinewidthV = 0.5; 
% % for i = 1:length(unAngle)
% %     subplot(4,2,i) 
% %     ind = find(elevAnD == unAngle(i)); 
% %     % Gas cell error bars
% %     % Turn colour into variable 
% %     % Zenith is markerfacecolor_v3 and goe slighter accoridng to increased elevation angle 
% %     % col = i/7; 
% %     % markerfacecolor_v3 = [0.09 0.54 0.03]; % AW3
% %     if i == 1 
% %         markerfacecolor_v3 = [0.09 0.32 0.07]; % AW3
% %     elseif i == 2
% %         markerfacecolor_v3 = [0.16 0.49 0.13]; 
% %          elseif i == 3
% %                 markerfacecolor_v3 = [0.18 0.52 0.15]; 
% %               elseif i == 4
% %                   markerfacecolor_v3 = [0.19 0.60 0.17]; 
% %                    elseif i == 5
% %                        markerfacecolor_v3 = [0.39 0.81 0.35]; 
% %                     elseif i == 6 
% %                             markerfacecolor_v3 = [0.61 0.83 0.59]; 
% %                        elseif i == 7 
% %                                 markerfacecolor_v3 = [0.70,0.92,0.67]; 
% %     end
% %     % For DOAS data 
% %     e = errorbar(cellConD(ind), SO2D(ind), cellConDErr(ind), cellConDErr(ind), 'horizontal', 'LineStyle', 'none', 'Marker', 'none', 'Color', errorbarcolor, 'linewidth', errorbarlinewidthH); hold on; 
% %     e.CapSize = capsize; 
% %     % Expected 
% %     et = plot(cellConLM(ind), cellConLMconv(ind), ':k'); et.LineWidth = trendlinewidthmain; hold on; 
% %     % For linear model data 
% %     e = errorbar(cellConLM(ind), SO2LMrOffset(ind), cellConLMErr(ind), cellConLMErr(ind), 'horizontal', 'LineStyle', 'none', 'Marker', 'none', 'Color', errorbarcolor, 'linewidth', errorbarlinewidthV); hold on; 
% %     e.CapSize = capsize; 
% %     % DOAS SO2 fit coefficients 
% %     d = errorbar(cellConD(ind), SO2D(ind), SO2ErrD(ind), SO2ErrD(ind), 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);
% %     d.MarkerSize = markersize; 
% %     d.MarkerFaceColor = markerfacecolorD; 
% %     d.MarkerEdgeColor = markeredgecolor; 
% %     d.CapSize = capsize; 
% %     xlim([minX maxX]);
% %     minY = -5e+17;
% %     ylim([minY maxY]);
% %     %     ylabel('(d)SCD SO_2 (molec/cm^2)')
% %     PrintEA = 90 - unAngle(i); 
% %     h = title(sprintf('%d°' , PrintEA)); h.FontSize = 12; 
% %     % xlabel('Gas cell concentration (ppm·m)')
% %     % % mEst and standard error 
% %     % % Real and complex
% %     % SO2LMr = mEstLM_v3_real(2,:); % SO2 is G2 (I0, SO2, O3, Ring, Bshift) 
% %     % % SO2LMri = mEstLM_v3_realimag(2,:); % SO2 is G2 (I0, SO2, O3, Ring, Bshift) 
% %     % SO2ErrLMr = stdErrLM_v3_real(2,:); 
% %     % % SO2ErrLMri = stdErrLM_v3_realimag(2,:); 
% %     % % Plots mEst 
% %     % e = errorbar(cellConLM, SO2LMr, cellConLMErr, cellConLMErr, 'horizontal', 'LineStyle', 'none', 'Marker', 'none', 'Color', errorbarcolor, 'linewidth', errorbarlinewidthV); hold on; 
% %     % e.CapSize = capsize; 
% %     % Real
% %     lr = errorbar(cellConLM(ind), SO2LMrOffset(ind), SO2ErrLMr(ind), SO2ErrLMr(ind), 'vertical', 'LineStyle', 'none', 'Marker', 'D', 'Color', markerfacecolor_v3, 'linewidth', errorbarlinewidthV);
% %     lr.MarkerFaceColor = markerfacecolor_v3; 
% %     lr.MarkerEdgeColor = markeredgecolor; 
% %     lr.MarkerSize = markersize; 
% %     lr.CapSize = capsize; 
% %     % Real imaginary 
% %     % lri = errorbar(cellConLM, SO2LMri, SO2ErrLMri, SO2ErrLMri, 'vertical', 'LineStyle', 'none', 'Marker', 'v', 'Color', markerfacecolor_v3, 'linewidth', errorbarlinewidthV);
% %     % lri.MarkerFaceColor = 'none'; 
% %     % lri.MarkerEdgeColor = markerfacecolorLMAW1; 
% %     if i == 3
% %         ylabel('(d)SCD SO_2 (molec/cm^2)')
% %     else
% %     end
% %     if i == 7 
% %        xlabel('Gas cell concentration (ppm·m)')
% %     end 
% %     Fig = gca; Fig.FontSize = 12; set(gcf,'color','w'); Fig.Box = 'on';
% % end 
% % %% Make any adjustments 
% % % % SAVES FIGURE
% % % fnOut = 'Results_angle';
% % % fname = fullfile(outDirFig, fnOut);
% % % saveas(gca, fname)
