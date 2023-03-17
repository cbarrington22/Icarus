% Synthetic spectra - analysis window (all frequencies)

% Must run SyntheticSpectra_AllWavelengths.m first to find WT of d and G
% This script extracts the data from that plotted in Fig. 6 
 
outDir = '/Users/charlotteb/Documents/Chapter 2/Synthetic spectra/Analysis window cropped frequnecy Fig 10 and Fig 11/'; 
outDirFig = '/Users/charlotteb/Documents/Chapter 2/Synthetic spectra/Analysis window cropped frequnecy Fig 10 and Fig 11/'; 

% Crops spatial frequency  
% Uses... WTd_2p5e15_AW1 etc. from SyntheticSpectra_analysisWindow_allF
% Selects F range 
yMin = 6e-6; 
yMax = 6e-5;
% Allows check
figure 
subplot(1,3,1)
pcolor(lambda, F, real(WT_SO2n)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
xlim([min(lambda) max(lambda)]); set(gca,'YScale','log'); % ylabel('Spatial frequency (cycles/\lambda)', 'FontSize', 10);
cb = colorbar; cb.Label.Position(1) = 3.5;
ylabel(cb,'Amplitude', 'Rotation', -90, 'FontSize', 10);
xlabel('\lambda (nm)', 'FontSize', 14); 
set(gcf,'color','w');
ylim([yMin yMax])
subplot(1,3,2)
pcolor(lambda, F, real(WT_I0)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
xlim([min(lambda) max(lambda)]); set(gca,'YScale','log'); % ylabel('Spatial frequency (cycles/\lambda)', 'FontSize', 10);
cb = colorbar; cb.Label.Position(1) = 3.5;
ylabel(cb,'Amplitude', 'Rotation', -90, 'FontSize', 10);
xlabel('\lambda (nm)', 'FontSize', 14); 
set(gcf,'color','w');
ylim([yMin yMax])
subplot(1,3,3)
pcolor(lambda, F, real(WTd_2p5e18)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
xlim([min(lambda) max(lambda)]); set(gca,'YScale','log'); % ylabel('Spatial frequency (cycles/\lambda)', 'FontSize', 10);
cb = colorbar; cb.Label.Position(1) = 3.5;
ylabel(cb,'Amplitude', 'Rotation', -90, 'FontSize', 10);
xlabel('\lambda (nm)', 'FontSize', 14); 
set(gcf,'color','w');
ylim([yMin yMax])
% Defines frequency limits accoridng to displayed axes
% Finds corresponding index
[~,iyMin] = min(abs(F-yMin));
[~,iyMax] = min(abs(F-yMax));
% Extracts window for each AW1-3 (wavelength) 
% G
% I0 
WT_I0_AW1_F = WT_I0_AW1(iyMax:iyMin, :); % AW1 
WT_I0_AW2_F = WT_I0_AW2(iyMax:iyMin, :); % AW2 
WT_I0_AW3_F = WT_I0_AW3(iyMax:iyMin, :); % AW3 
% SO2 
WT_SO2n_AW1_F = WT_SO2n_AW1(iyMax:iyMin, :); % AW1 
WT_SO2n_AW2_F = WT_SO2n_AW2(iyMax:iyMin, :); % AW2 
WT_SO2n_AW3_F = WT_SO2n_AW3(iyMax:iyMin, :); % AW3 
% d 
WTd_2p5e15_AW1_F = WTd_2p5e15_AW1(iyMax:iyMin, :); % AW1 
WTd_7p5e15_AW1_F = WTd_7p5e15_AW1(iyMax:iyMin, :);
WTd_2p5e16_AW1_F = WTd_2p5e16_AW1(iyMax:iyMin, :);
WTd_7p5e16_AW1_F = WTd_7p5e16_AW1(iyMax:iyMin, :);
WTd_2p5e17_AW1_F = WTd_2p5e17_AW1(iyMax:iyMin, :);
WTd_7p5e17_AW1_F = WTd_7p5e17_AW1(iyMax:iyMin, :);
WTd_2p5e18_AW1_F = WTd_2p5e18_AW1(iyMax:iyMin, :);
WTd_7p5e18_AW1_F = WTd_7p5e18_AW1(iyMax:iyMin, :);
WTd_2p5e19_AW1_F = WTd_2p5e19_AW1(iyMax:iyMin, :);
WTd_7p5e19_AW1_F = WTd_7p5e19_AW1(iyMax:iyMin, :);
WTd_2p5e15_AW2_F = WTd_2p5e15_AW2(iyMax:iyMin, :); % AW2 
WTd_7p5e15_AW2_F = WTd_7p5e15_AW2(iyMax:iyMin, :);
WTd_2p5e16_AW2_F = WTd_2p5e16_AW2(iyMax:iyMin, :);
WTd_7p5e16_AW2_F = WTd_7p5e16_AW2(iyMax:iyMin, :);
WTd_2p5e17_AW2_F = WTd_2p5e17_AW2(iyMax:iyMin, :);
WTd_7p5e17_AW2_F = WTd_7p5e17_AW2(iyMax:iyMin, :);
WTd_2p5e18_AW2_F = WTd_2p5e18_AW2(iyMax:iyMin, :);
WTd_7p5e18_AW2_F = WTd_7p5e18_AW2(iyMax:iyMin, :);
WTd_2p5e19_AW2_F = WTd_2p5e19_AW2(iyMax:iyMin, :);
WTd_7p5e19_AW2_F = WTd_7p5e19_AW2(iyMax:iyMin, :);
WTd_2p5e15_AW3_F = WTd_2p5e15_AW3(iyMax:iyMin, :); % AW3 
WTd_7p5e15_AW3_F = WTd_7p5e15_AW3(iyMax:iyMin, :);
WTd_2p5e16_AW3_F = WTd_2p5e16_AW3(iyMax:iyMin, :);
WTd_7p5e16_AW3_F = WTd_7p5e16_AW3(iyMax:iyMin, :);
WTd_2p5e17_AW3_F = WTd_2p5e17_AW3(iyMax:iyMin, :);
WTd_7p5e17_AW3_F = WTd_7p5e17_AW3(iyMax:iyMin, :);
WTd_2p5e18_AW3_F = WTd_2p5e18_AW3(iyMax:iyMin, :);
WTd_7p5e18_AW3_F = WTd_7p5e18_AW3(iyMax:iyMin, :);
WTd_2p5e19_AW3_F = WTd_2p5e19_AW3(iyMax:iyMin, :);
WTd_7p5e19_AW3_F = WTd_7p5e19_AW3(iyMax:iyMin, :);
% F 
Fc = F(iyMax:iyMin);
% Figure 10 - plot linear model for analysis windows now with spatial frequency cropped 
figure 
% Subplots 
% First column of figure is d (plot for 2.5e18) 
subplot(3,3,1) 
pcolor(lambda_AW1, Fc, real(WTd_2p5e18_AW1_F)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
xlim([min(lambda_AW1) max(lambda_AW1)]); set(gca,'YScale','log'); % ylabel('Spatial frequency (cycles/\lambda)', 'FontSize', 10);
% cb = colorbar; cb.Label.Position(1) = 3.5;
% ylabel(cb,'Amplitude', 'Rotation', -90, 'FontSize', 10);
set(gcf,'color','w');
title('d')
subtitle('J(\lambda) prime')
subplot(3,3,4) 
pcolor(lambda_AW2, Fc, real(WTd_2p5e18_AW2_F)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
xlim([min(lambda_AW2) max(lambda_AW2)]); set(gca,'YScale','log'); ylabel('Spatial frequency (cycles/\lambda)', 'FontSize', 14);
% cb = colorbar; cb.Label.Position(1) = 3.5;
% ylabel(cb,'Amplitude', 'Rotation', -90, 'FontSize', 10);
set(gcf,'color','w');
subplot(3,3,7) 
pcolor(lambda_AW3, Fc, real(WTd_2p5e18_AW3_F)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
xlim([min(lambda_AW3) max(lambda_AW3)]); set(gca,'YScale','log'); % ylabel('Spatial frequency (cycles/\lambda)', 'FontSize', 10);
% cb = colorbar; cb.Label.Position(1) = 3.5;
% ylabel(cb,'Amplitude', 'Rotation', -90, 'FontSize', 10);
set(gcf,'color','w');
subplot(3,3,2) 
pcolor(lambda_AW1, Fc, real(WT_I0_AW1_F)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
xlim([min(lambda_AW1) max(lambda_AW1)]); set(gca,'YScale','log'); % ylabel('Spatial frequency (cycles/\lambda)', 'FontSize', 10);
% cb = colorbar; cb.Label.Position(1) = 3.5;
% ylabel(cb,'Amplitude', 'Rotation', -90, 'FontSize', 10);
set(gcf,'color','w');
title('G') 
subtitle('J_0(\lambda) prime')
subplot(3,3,5) 
pcolor(lambda_AW2, Fc, real(WT_I0_AW2_F)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
xlim([min(lambda_AW2) max(lambda_AW2)]); set(gca,'YScale','log'); % ylabel('Spatial frequency (cycles/\lambda)', 'FontSize', 10);
% cb = colorbar; cb.Label.Position(1) = 3.5;
% ylabel(cb,'Amplitude', 'Rotation', -90, 'FontSize', 10);
set(gcf,'color','w');
subplot(3,3,8) 
pcolor(lambda_AW3, Fc, real(WT_I0_AW3_F)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
xlim([min(lambda_AW3) max(lambda_AW3)]); set(gca,'YScale','log'); % ylabel('Spatial frequency (cycles/\lambda)', 'FontSize', 10);
% cb = colorbar; cb.Label.Position(1) = 3.5;
% ylabel(cb,'Amplitude', 'Rotation', -90, 'FontSize', 10);
xlabel('\lambda (nm)', 'FontSize', 14); 
set(gcf,'color','w');
subplot(3,3,3) 
pcolor(lambda_AW1, Fc, real(WT_SO2n_AW1_F)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
xlim([min(lambda_AW1) max(lambda_AW1)]); set(gca,'YScale','log'); % ylabel('Spatial frequency (cycles/\lambda)', 'FontSize', 10);
% cb = colorbar; cb.Label.Position(1) = 3.5;
% ylabel(cb,'Amplitude', 'Rotation', -90, 'FontSize', 10);
subtitle('\sigma_S_O_2 prime')
set(gcf,'color','w');
subplot(3,3,6) 
pcolor(lambda_AW2, Fc, real(WT_SO2n_AW2_F)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
xlim([min(lambda_AW2) max(lambda_AW2)]); set(gca,'YScale','log'); % ylabel('Spatial frequency (cycles/\lambda)', 'FontSize', 10);
cb = colorbar; cb.Label.Position(1) = 3.5;
ylabel(cb,'Amplitude', 'Rotation', -90, 'FontSize', 14);
set(gcf,'color','w');
subplot(3,3,9) 
pcolor(lambda_AW3, Fc, real(WT_SO2n_AW3_F)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
xlim([min(lambda_AW3) max(lambda_AW3)]); set(gca,'YScale','log'); % ylabel('Spatial frequency (cycles/\lambda)', 'FontSize', 10); 
% cb = colorbar; cb.Label.Position(1) = 3.5;
% ylabel(cb,'Amplitude', 'Rotation', -90, 'FontSize', 10);
set(gcf,'color','w');
% Save 
fname = fullfile(outDir, 'mEst_real_imag_AW_F');
save(fname); % Saves all variables  
%% Make any adjustments 
% SAVES FIGURE
% fnOut = 'Fig10_dG_AW_F';
% fname = fullfile(outDirFig, fnOut);
% saveas(gca, fname)

% Finds mEst for AW1_F 
% Takes real and imaginary parts: I0 and SO2
% Real
WT_I0_AW1_F_real = real(WT_I0_AW1_F); 
WT_SO2n_AW1_F_real = real(WT_SO2n_AW1_F); 
% Imaginary 
WT_I0_AW1_F_imag = imag(WT_I0_AW1_F); 
WT_SO2n_AW1_F_imag = imag(WT_SO2n_AW1_F); 
% Vectorises 
% Real
WT_I0_AW1_F_real_v = WT_I0_AW1_F_real(:); 
WT_SO2n_AW1_F_real_v = WT_SO2n_AW1_F_real(:); 
% Imaginary
WT_I0_AW1_F_imag_v = WT_I0_AW1_F_imag(:); 
WT_SO2n_AW1_F_imag_v = WT_SO2n_AW1_F_imag(:); 
% Combines real and complex magnitude 
% I0
WT_I0_AW1_F_realimag = [WT_I0_AW1_F_real_v; WT_I0_AW1_F_imag_v]; 
% SO2
WT_SO2n_AW1_F_realimag = [WT_SO2n_AW1_F_real_v; WT_SO2n_AW1_F_imag_v]; 
% Creates design matrix, G 
G_AW1_F_real = [WT_I0_AW1_F_real_v, WT_SO2n_AW1_F_real_v]; 
G_AW1_F_realimag = [WT_I0_AW1_F_realimag, WT_SO2n_AW1_F_realimag]; 
% Remove NaN 
G_AW1_F_real = G_AW1_F_real(all(~isnan(G_AW1_F_real), 2),:);
G_AW1_F_realimag = G_AW1_F_realimag(all(~isnan(G_AW1_F_realimag), 2),:);
% G transponse 
Gt_AW1_F_real = G_AW1_F_real'; 
Gt_AW1_F_realimag = G_AW1_F_realimag'; 

% Takes real and imaginary parts: I0 and SO2
% Real
WTd_2p5e15_AW1_F_real = real(WTd_2p5e15_AW1_F); 
WTd_7p5e15_AW1_F_real = real(WTd_7p5e15_AW1_F); 
WTd_2p5e16_AW1_F_real = real(WTd_2p5e16_AW1_F); 
WTd_7p5e16_AW1_F_real = real(WTd_7p5e16_AW1_F); 
WTd_2p5e17_AW1_F_real = real(WTd_2p5e17_AW1_F); 
WTd_7p5e17_AW1_F_real = real(WTd_7p5e17_AW1_F); 
WTd_2p5e18_AW1_F_real = real(WTd_2p5e18_AW1_F); 
WTd_7p5e18_AW1_F_real = real(WTd_7p5e18_AW1_F); 
WTd_2p5e19_AW1_F_real = real(WTd_2p5e19_AW1_F); 
WTd_7p5e19_AW1_F_real = real(WTd_7p5e19_AW1_F); 
% Imaginary 
WTd_2p5e15_AW1_F_imag = imag(WTd_2p5e15_AW1_F); 
WTd_7p5e15_AW1_F_imag = imag(WTd_7p5e15_AW1_F); 
WTd_2p5e16_AW1_F_imag = imag(WTd_2p5e16_AW1_F); 
WTd_7p5e16_AW1_F_imag = imag(WTd_7p5e16_AW1_F); 
WTd_2p5e17_AW1_F_imag = imag(WTd_2p5e17_AW1_F); 
WTd_7p5e17_AW1_F_imag = imag(WTd_7p5e17_AW1_F); 
WTd_2p5e18_AW1_F_imag = imag(WTd_2p5e18_AW1_F); 
WTd_7p5e18_AW1_F_imag = imag(WTd_7p5e18_AW1_F); 
WTd_2p5e19_AW1_F_imag = imag(WTd_2p5e19_AW1_F); 
WTd_7p5e19_AW1_F_imag = imag(WTd_7p5e19_AW1_F); 
% Vectorises 
% Real
WTd_2p5e15_AW1_F_real_v = WTd_2p5e15_AW1_F_real(:); 
WTd_7p5e15_AW1_F_real_v = WTd_7p5e15_AW1_F_real(:); 
WTd_2p5e16_AW1_F_real_v = WTd_2p5e16_AW1_F_real(:); 
WTd_7p5e16_AW1_F_real_v = WTd_7p5e16_AW1_F_real(:); 
WTd_2p5e17_AW1_F_real_v = WTd_2p5e17_AW1_F_real(:); 
WTd_7p5e17_AW1_F_real_v = WTd_7p5e17_AW1_F_real(:); 
WTd_2p5e18_AW1_F_real_v = WTd_2p5e18_AW1_F_real(:); 
WTd_7p5e18_AW1_F_real_v = WTd_7p5e18_AW1_F_real(:); 
WTd_2p5e19_AW1_F_real_v = WTd_2p5e19_AW1_F_real(:); 
WTd_7p5e19_AW1_F_real_v = WTd_7p5e19_AW1_F_real(:); 
% Imaginary
WTd_2p5e15_AW1_F_imag_v = WTd_2p5e15_AW1_F_imag(:); 
WTd_7p5e15_AW1_F_imag_v = WTd_7p5e15_AW1_F_imag(:); 
WTd_2p5e16_AW1_F_imag_v = WTd_2p5e16_AW1_F_imag(:); 
WTd_7p5e16_AW1_F_imag_v = WTd_7p5e16_AW1_F_imag(:); 
WTd_2p5e17_AW1_F_imag_v = WTd_2p5e17_AW1_F_imag(:); 
WTd_7p5e17_AW1_F_imag_v = WTd_7p5e17_AW1_F_imag(:); 
WTd_2p5e18_AW1_F_imag_v = WTd_2p5e18_AW1_F_imag(:); 
WTd_7p5e18_AW1_F_imag_v = WTd_7p5e18_AW1_F_imag(:); 
WTd_2p5e19_AW1_F_imag_v = WTd_2p5e19_AW1_F_imag(:); 
WTd_7p5e19_AW1_F_imag_v = WTd_7p5e19_AW1_F_imag(:); 
% Combines real and complex magnitude 
WTd_2p5e15_AW1_F_realimag_v = [WTd_2p5e15_AW1_F_real_v; WTd_2p5e15_AW1_F_imag_v]; 
WTd_7p5e15_AW1_F_realimag_v = [WTd_7p5e15_AW1_F_real_v; WTd_7p5e15_AW1_F_imag_v]; 
WTd_2p5e16_AW1_F_realimag_v = [WTd_2p5e16_AW1_F_real_v; WTd_2p5e16_AW1_F_imag_v]; 
WTd_7p5e16_AW1_F_realimag_v = [WTd_7p5e16_AW1_F_real_v; WTd_7p5e16_AW1_F_imag_v]; 
WTd_2p5e17_AW1_F_realimag_v = [WTd_2p5e17_AW1_F_real_v; WTd_2p5e17_AW1_F_imag_v]; 
WTd_7p5e17_AW1_F_realimag_v = [WTd_7p5e17_AW1_F_real_v; WTd_7p5e17_AW1_F_imag_v]; 
WTd_2p5e18_AW1_F_realimag_v = [WTd_2p5e18_AW1_F_real_v; WTd_2p5e18_AW1_F_imag_v]; 
WTd_7p5e18_AW1_F_realimag_v = [WTd_7p5e18_AW1_F_real_v; WTd_7p5e18_AW1_F_imag_v]; 
WTd_2p5e19_AW1_F_realimag_v = [WTd_2p5e19_AW1_F_real_v; WTd_2p5e19_AW1_F_imag_v]; 
WTd_7p5e19_AW1_F_realimag_v = [WTd_7p5e19_AW1_F_real_v; WTd_7p5e19_AW1_F_imag_v]; 
% Remove NaN 
% Real only 
d_2p5e15_AW1_F_real_v = WTd_2p5e15_AW1_F_real_v(all(~isnan(WTd_2p5e15_AW1_F_real_v), 2),:);
d_7p5e15_AW1_F_real_v = WTd_7p5e15_AW1_F_real_v(all(~isnan(WTd_7p5e15_AW1_F_real_v), 2),:);
d_2p5e16_AW1_F_real_v = WTd_2p5e16_AW1_F_real_v(all(~isnan(WTd_2p5e16_AW1_F_real_v), 2),:);
d_7p5e16_AW1_F_real_v = WTd_7p5e16_AW1_F_real_v(all(~isnan(WTd_7p5e16_AW1_F_real_v), 2),:);
d_2p5e17_AW1_F_real_v = WTd_2p5e17_AW1_F_real_v(all(~isnan(WTd_2p5e17_AW1_F_real_v), 2),:);
d_7p5e17_AW1_F_real_v = WTd_7p5e17_AW1_F_real_v(all(~isnan(WTd_7p5e17_AW1_F_real_v), 2),:);
d_2p5e18_AW1_F_real_v = WTd_2p5e18_AW1_F_real_v(all(~isnan(WTd_2p5e18_AW1_F_real_v), 2),:);
d_7p5e18_AW1_F_real_v = WTd_7p5e18_AW1_F_real_v(all(~isnan(WTd_7p5e18_AW1_F_real_v), 2),:);
d_2p5e19_AW1_F_real_v = WTd_2p5e19_AW1_F_real_v(all(~isnan(WTd_2p5e19_AW1_F_real_v), 2),:);
d_7p5e19_AW1_F_real_v = WTd_7p5e19_AW1_F_real_v(all(~isnan(WTd_7p5e19_AW1_F_real_v), 2),:);
% Real and imag (complex magnitude) 
d_2p5e15_AW1_F_realimag_v = WTd_2p5e15_AW1_F_realimag_v(all(~isnan(WTd_2p5e15_AW1_F_realimag_v), 2),:);
d_7p5e15_AW1_F_realimag_v = WTd_7p5e15_AW1_F_realimag_v(all(~isnan(WTd_7p5e15_AW1_F_realimag_v), 2),:);
d_2p5e16_AW1_F_realimag_v = WTd_2p5e16_AW1_F_realimag_v(all(~isnan(WTd_2p5e16_AW1_F_realimag_v), 2),:);
d_7p5e16_AW1_F_realimag_v = WTd_7p5e16_AW1_F_realimag_v(all(~isnan(WTd_7p5e16_AW1_F_realimag_v), 2),:);
d_2p5e17_AW1_F_realimag_v = WTd_2p5e17_AW1_F_realimag_v(all(~isnan(WTd_2p5e17_AW1_F_realimag_v), 2),:);
d_7p5e17_AW1_F_realimag_v = WTd_7p5e17_AW1_F_realimag_v(all(~isnan(WTd_7p5e17_AW1_F_realimag_v), 2),:);
d_2p5e18_AW1_F_realimag_v = WTd_2p5e18_AW1_F_realimag_v(all(~isnan(WTd_2p5e18_AW1_F_realimag_v), 2),:);
d_7p5e18_AW1_F_realimag_v = WTd_7p5e18_AW1_F_realimag_v(all(~isnan(WTd_7p5e18_AW1_F_realimag_v), 2),:);
d_2p5e19_AW1_F_realimag_v = WTd_2p5e19_AW1_F_realimag_v(all(~isnan(WTd_2p5e19_AW1_F_realimag_v), 2),:);
d_7p5e19_AW1_F_realimag_v = WTd_7p5e19_AW1_F_realimag_v(all(~isnan(WTd_7p5e19_AW1_F_realimag_v), 2),:);
% Find mEst 
% A 
A_AW1_F_real = Gt_AW1_F_real * G_AW1_F_real;
A_AW1_F_realimag = Gt_AW1_F_realimag * G_AW1_F_realimag;
% B 
% Real only
B_AW1_F_real_d2p5e15 = Gt_AW1_F_real * d_2p5e15_AW1_F_real_v;
B_AW1_F_real_d7p5e15 = Gt_AW1_F_real * d_7p5e15_AW1_F_real_v;
B_AW1_F_real_d2p5e16 = Gt_AW1_F_real * d_2p5e16_AW1_F_real_v;
B_AW1_F_real_d7p5e16 = Gt_AW1_F_real * d_7p5e16_AW1_F_real_v;
B_AW1_F_real_d2p5e17 = Gt_AW1_F_real * d_2p5e17_AW1_F_real_v;
B_AW1_F_real_d7p5e17 = Gt_AW1_F_real * d_7p5e17_AW1_F_real_v;
B_AW1_F_real_d2p5e18 = Gt_AW1_F_real * d_2p5e18_AW1_F_real_v;
B_AW1_F_real_d7p5e18 = Gt_AW1_F_real * d_7p5e18_AW1_F_real_v;
B_AW1_F_real_d2p5e19 = Gt_AW1_F_real * d_2p5e19_AW1_F_real_v;
B_AW1_F_real_d7p5e19 = Gt_AW1_F_real * d_7p5e19_AW1_F_real_v;
% Real and imag (complex magnitude) 
B_AW1_F_realimag_d2p5e15 = Gt_AW1_F_realimag * d_2p5e15_AW1_F_realimag_v;
B_AW1_F_realimag_d7p5e15 = Gt_AW1_F_realimag * d_7p5e15_AW1_F_realimag_v;
B_AW1_F_realimag_d2p5e16 = Gt_AW1_F_realimag * d_2p5e16_AW1_F_realimag_v;
B_AW1_F_realimag_d7p5e16 = Gt_AW1_F_realimag * d_7p5e16_AW1_F_realimag_v;
B_AW1_F_realimag_d2p5e17 = Gt_AW1_F_realimag * d_2p5e17_AW1_F_realimag_v;
B_AW1_F_realimag_d7p5e17 = Gt_AW1_F_realimag * d_7p5e17_AW1_F_realimag_v;
B_AW1_F_realimag_d2p5e18 = Gt_AW1_F_realimag * d_2p5e18_AW1_F_realimag_v;
B_AW1_F_realimag_d7p5e18 = Gt_AW1_F_realimag * d_7p5e18_AW1_F_realimag_v;
B_AW1_F_realimag_d2p5e19 = Gt_AW1_F_realimag * d_2p5e19_AW1_F_realimag_v;
B_AW1_F_realimag_d7p5e19 = Gt_AW1_F_realimag * d_7p5e19_AW1_F_realimag_v;
% mEst 
% Real only 
mEst_AW1_F_real_2p5e15 = A_AW1_F_real\B_AW1_F_real_d2p5e15;
mEst_AW1_F_real_7p5e15 = A_AW1_F_real\B_AW1_F_real_d7p5e15;
mEst_AW1_F_real_2p5e16 = A_AW1_F_real\B_AW1_F_real_d2p5e16;
mEst_AW1_F_real_7p5e16 = A_AW1_F_real\B_AW1_F_real_d7p5e16;
mEst_AW1_F_real_2p5e17 = A_AW1_F_real\B_AW1_F_real_d2p5e17;
mEst_AW1_F_real_7p5e17 = A_AW1_F_real\B_AW1_F_real_d7p5e17;
mEst_AW1_F_real_2p5e18 = A_AW1_F_real\B_AW1_F_real_d2p5e18;
mEst_AW1_F_real_7p5e18 = A_AW1_F_real\B_AW1_F_real_d7p5e18;
mEst_AW1_F_real_2p5e19 = A_AW1_F_real\B_AW1_F_real_d2p5e19;
mEst_AW1_F_real_7p5e19 = A_AW1_F_real\B_AW1_F_real_d7p5e19;
% Real and imag (complex magnitude) 
mEst_AW1_F_realimag_2p5e15 = A_AW1_F_realimag\B_AW1_F_realimag_d2p5e15;
mEst_AW1_F_realimag_7p5e15 = A_AW1_F_realimag\B_AW1_F_realimag_d7p5e15;
mEst_AW1_F_realimag_2p5e16 = A_AW1_F_realimag\B_AW1_F_realimag_d2p5e16;
mEst_AW1_F_realimag_7p5e16 = A_AW1_F_realimag\B_AW1_F_realimag_d7p5e16;
mEst_AW1_F_realimag_2p5e17 = A_AW1_F_realimag\B_AW1_F_realimag_d2p5e17;
mEst_AW1_F_realimag_7p5e17 = A_AW1_F_realimag\B_AW1_F_realimag_d7p5e17;
mEst_AW1_F_realimag_2p5e18 = A_AW1_F_realimag\B_AW1_F_realimag_d2p5e18;
mEst_AW1_F_realimag_7p5e18 = A_AW1_F_realimag\B_AW1_F_realimag_d7p5e18;
mEst_AW1_F_realimag_2p5e19 = A_AW1_F_realimag\B_AW1_F_realimag_d2p5e19;
mEst_AW1_F_realimag_7p5e19 = A_AW1_F_realimag\B_AW1_F_realimag_d7p5e19;
% Results 
% concentration_SO2 = [1.5e15; 7.5e15; 2.5e16; 7.5e16; 2.5e17; 7.5e17; 2.5e18; 7.5e18; 2.5e19; 7.5e19]; 
mEst_AW1_F_real = [mEst_AW1_F_real_2p5e15, mEst_AW1_F_real_7p5e15, mEst_AW1_F_real_2p5e16, mEst_AW1_F_real_7p5e16, mEst_AW1_F_real_2p5e17, mEst_AW1_F_real_7p5e17, mEst_AW1_F_real_2p5e18, mEst_AW1_F_real_7p5e18, mEst_AW1_F_real_2p5e19, mEst_AW1_F_real_7p5e19];  
mEst_AW1_F_realimag = [mEst_AW1_F_realimag_2p5e15, mEst_AW1_F_realimag_7p5e15, mEst_AW1_F_realimag_2p5e16, mEst_AW1_F_realimag_7p5e16, mEst_AW1_F_realimag_2p5e17, mEst_AW1_F_realimag_7p5e17, mEst_AW1_F_realimag_2p5e18, mEst_AW1_F_realimag_7p5e18, mEst_AW1_F_realimag_2p5e19, mEst_AW1_F_realimag_7p5e19];  

% Finds mEst for AW2_F 
% Takes real and imaginary parts: I0 and SO2
% Real
WT_I0_AW2_F_real = real(WT_I0_AW2_F); 
WT_SO2n_AW2_F_real = real(WT_SO2n_AW2_F); 
% Imaginary 
WT_I0_AW2_F_imag = imag(WT_I0_AW2_F); 
WT_SO2n_AW2_F_imag = imag(WT_SO2n_AW2_F); 
% Vectorises 
% Real
WT_I0_AW2_F_real_v = WT_I0_AW2_F_real(:); 
WT_SO2n_AW2_F_real_v = WT_SO2n_AW2_F_real(:); 
% Imaginary
WT_I0_AW2_F_imag_v = WT_I0_AW2_F_imag(:); 
WT_SO2n_AW2_F_imag_v = WT_SO2n_AW2_F_imag(:); 
% Combines real and complex magnitude 
% I0
WT_I0_AW2_F_realimag = [WT_I0_AW2_F_real_v; WT_I0_AW2_F_imag_v]; 
% SO2
WT_SO2n_AW2_F_realimag = [WT_SO2n_AW2_F_real_v; WT_SO2n_AW2_F_imag_v]; 
% Creates design matrix, G 
G_AW2_F_real = [WT_I0_AW2_F_real_v, WT_SO2n_AW2_F_real_v]; 
G_AW2_F_realimag = [WT_I0_AW2_F_realimag, WT_SO2n_AW2_F_realimag]; 
% Remove NaN 
G_AW2_F_real = G_AW2_F_real(all(~isnan(G_AW2_F_real), 2),:);
G_AW2_F_realimag = G_AW2_F_realimag(all(~isnan(G_AW2_F_realimag), 2),:);
% G transponse 
Gt_AW2_F_real = G_AW2_F_real'; 
Gt_AW2_F_realimag = G_AW2_F_realimag'; 
% Takes real and imaginary parts: I0 and SO2
% Real
WTd_2p5e15_AW2_F_real = real(WTd_2p5e15_AW2_F); 
WTd_7p5e15_AW2_F_real = real(WTd_7p5e15_AW2_F); 
WTd_2p5e16_AW2_F_real = real(WTd_2p5e16_AW2_F); 
WTd_7p5e16_AW2_F_real = real(WTd_7p5e16_AW2_F); 
WTd_2p5e17_AW2_F_real = real(WTd_2p5e17_AW2_F); 
WTd_7p5e17_AW2_F_real = real(WTd_7p5e17_AW2_F); 
WTd_2p5e18_AW2_F_real = real(WTd_2p5e18_AW2_F); 
WTd_7p5e18_AW2_F_real = real(WTd_7p5e18_AW2_F); 
WTd_2p5e19_AW2_F_real = real(WTd_2p5e19_AW2_F); 
WTd_7p5e19_AW2_F_real = real(WTd_7p5e19_AW2_F); 
% Imaginary 
WTd_2p5e15_AW2_F_imag = imag(WTd_2p5e15_AW2_F); 
WTd_7p5e15_AW2_F_imag = imag(WTd_7p5e15_AW2_F); 
WTd_2p5e16_AW2_F_imag = imag(WTd_2p5e16_AW2_F); 
WTd_7p5e16_AW2_F_imag = imag(WTd_7p5e16_AW2_F); 
WTd_2p5e17_AW2_F_imag = imag(WTd_2p5e17_AW2_F); 
WTd_7p5e17_AW2_F_imag = imag(WTd_7p5e17_AW2_F); 
WTd_2p5e18_AW2_F_imag = imag(WTd_2p5e18_AW2_F); 
WTd_7p5e18_AW2_F_imag = imag(WTd_7p5e18_AW2_F); 
WTd_2p5e19_AW2_F_imag = imag(WTd_2p5e19_AW2_F); 
WTd_7p5e19_AW2_F_imag = imag(WTd_7p5e19_AW2_F); 
% Vectorises 
% Real
WTd_2p5e15_AW2_F_real_v = WTd_2p5e15_AW2_F_real(:); 
WTd_7p5e15_AW2_F_real_v = WTd_7p5e15_AW2_F_real(:); 
WTd_2p5e16_AW2_F_real_v = WTd_2p5e16_AW2_F_real(:); 
WTd_7p5e16_AW2_F_real_v = WTd_7p5e16_AW2_F_real(:); 
WTd_2p5e17_AW2_F_real_v = WTd_2p5e17_AW2_F_real(:); 
WTd_7p5e17_AW2_F_real_v = WTd_7p5e17_AW2_F_real(:); 
WTd_2p5e18_AW2_F_real_v = WTd_2p5e18_AW2_F_real(:); 
WTd_7p5e18_AW2_F_real_v = WTd_7p5e18_AW2_F_real(:); 
WTd_2p5e19_AW2_F_real_v = WTd_2p5e19_AW2_F_real(:); 
WTd_7p5e19_AW2_F_real_v = WTd_7p5e19_AW2_F_real(:); 
% Imaginary
WTd_2p5e15_AW2_F_imag_v = WTd_2p5e15_AW2_F_imag(:); 
WTd_7p5e15_AW2_F_imag_v = WTd_7p5e15_AW2_F_imag(:); 
WTd_2p5e16_AW2_F_imag_v = WTd_2p5e16_AW2_F_imag(:); 
WTd_7p5e16_AW2_F_imag_v = WTd_7p5e16_AW2_F_imag(:); 
WTd_2p5e17_AW2_F_imag_v = WTd_2p5e17_AW2_F_imag(:); 
WTd_7p5e17_AW2_F_imag_v = WTd_7p5e17_AW2_F_imag(:); 
WTd_2p5e18_AW2_F_imag_v = WTd_2p5e18_AW2_F_imag(:); 
WTd_7p5e18_AW2_F_imag_v = WTd_7p5e18_AW2_F_imag(:); 
WTd_2p5e19_AW2_F_imag_v = WTd_2p5e19_AW2_F_imag(:); 
WTd_7p5e19_AW2_F_imag_v = WTd_7p5e19_AW2_F_imag(:); 
% Combines real and complex magnitude 
WTd_2p5e15_AW2_F_realimag_v = [WTd_2p5e15_AW2_F_real_v; WTd_2p5e15_AW2_F_imag_v]; 
WTd_7p5e15_AW2_F_realimag_v = [WTd_7p5e15_AW2_F_real_v; WTd_7p5e15_AW2_F_imag_v]; 
WTd_2p5e16_AW2_F_realimag_v = [WTd_2p5e16_AW2_F_real_v; WTd_2p5e16_AW2_F_imag_v]; 
WTd_7p5e16_AW2_F_realimag_v = [WTd_7p5e16_AW2_F_real_v; WTd_7p5e16_AW2_F_imag_v]; 
WTd_2p5e17_AW2_F_realimag_v = [WTd_2p5e17_AW2_F_real_v; WTd_2p5e17_AW2_F_imag_v]; 
WTd_7p5e17_AW2_F_realimag_v = [WTd_7p5e17_AW2_F_real_v; WTd_7p5e17_AW2_F_imag_v]; 
WTd_2p5e18_AW2_F_realimag_v = [WTd_2p5e18_AW2_F_real_v; WTd_2p5e18_AW2_F_imag_v]; 
WTd_7p5e18_AW2_F_realimag_v = [WTd_7p5e18_AW2_F_real_v; WTd_7p5e18_AW2_F_imag_v]; 
WTd_2p5e19_AW2_F_realimag_v = [WTd_2p5e19_AW2_F_real_v; WTd_2p5e19_AW2_F_imag_v]; 
WTd_7p5e19_AW2_F_realimag_v = [WTd_7p5e19_AW2_F_real_v; WTd_7p5e19_AW2_F_imag_v]; 
% Remove NaN 
% Real only 
d_2p5e15_AW2_F_real_v = WTd_2p5e15_AW2_F_real_v(all(~isnan(WTd_2p5e15_AW2_F_real_v), 2),:);
d_7p5e15_AW2_F_real_v = WTd_7p5e15_AW2_F_real_v(all(~isnan(WTd_7p5e15_AW2_F_real_v), 2),:);
d_2p5e16_AW2_F_real_v = WTd_2p5e16_AW2_F_real_v(all(~isnan(WTd_2p5e16_AW2_F_real_v), 2),:);
d_7p5e16_AW2_F_real_v = WTd_7p5e16_AW2_F_real_v(all(~isnan(WTd_7p5e16_AW2_F_real_v), 2),:);
d_2p5e17_AW2_F_real_v = WTd_2p5e17_AW2_F_real_v(all(~isnan(WTd_2p5e17_AW2_F_real_v), 2),:);
d_7p5e17_AW2_F_real_v = WTd_7p5e17_AW2_F_real_v(all(~isnan(WTd_7p5e17_AW2_F_real_v), 2),:);
d_2p5e18_AW2_F_real_v = WTd_2p5e18_AW2_F_real_v(all(~isnan(WTd_2p5e18_AW2_F_real_v), 2),:);
d_7p5e18_AW2_F_real_v = WTd_7p5e18_AW2_F_real_v(all(~isnan(WTd_7p5e18_AW2_F_real_v), 2),:);
d_2p5e19_AW2_F_real_v = WTd_2p5e19_AW2_F_real_v(all(~isnan(WTd_2p5e19_AW2_F_real_v), 2),:);
d_7p5e19_AW2_F_real_v = WTd_7p5e19_AW2_F_real_v(all(~isnan(WTd_7p5e19_AW2_F_real_v), 2),:);
% Real and imag (complex magnitude) 
d_2p5e15_AW2_F_realimag_v = WTd_2p5e15_AW2_F_realimag_v(all(~isnan(WTd_2p5e15_AW2_F_realimag_v), 2),:);
d_7p5e15_AW2_F_realimag_v = WTd_7p5e15_AW2_F_realimag_v(all(~isnan(WTd_7p5e15_AW2_F_realimag_v), 2),:);
d_2p5e16_AW2_F_realimag_v = WTd_2p5e16_AW2_F_realimag_v(all(~isnan(WTd_2p5e16_AW2_F_realimag_v), 2),:);
d_7p5e16_AW2_F_realimag_v = WTd_7p5e16_AW2_F_realimag_v(all(~isnan(WTd_7p5e16_AW2_F_realimag_v), 2),:);
d_2p5e17_AW2_F_realimag_v = WTd_2p5e17_AW2_F_realimag_v(all(~isnan(WTd_2p5e17_AW2_F_realimag_v), 2),:);
d_7p5e17_AW2_F_realimag_v = WTd_7p5e17_AW2_F_realimag_v(all(~isnan(WTd_7p5e17_AW2_F_realimag_v), 2),:);
d_2p5e18_AW2_F_realimag_v = WTd_2p5e18_AW2_F_realimag_v(all(~isnan(WTd_2p5e18_AW2_F_realimag_v), 2),:);
d_7p5e18_AW2_F_realimag_v = WTd_7p5e18_AW2_F_realimag_v(all(~isnan(WTd_7p5e18_AW2_F_realimag_v), 2),:);
d_2p5e19_AW2_F_realimag_v = WTd_2p5e19_AW2_F_realimag_v(all(~isnan(WTd_2p5e19_AW2_F_realimag_v), 2),:);
d_7p5e19_AW2_F_realimag_v = WTd_7p5e19_AW2_F_realimag_v(all(~isnan(WTd_7p5e19_AW2_F_realimag_v), 2),:);
% Find mEst 
% A 
A_AW2_F_real = Gt_AW2_F_real * G_AW2_F_real;
A_AW2_F_realimag = Gt_AW2_F_realimag * G_AW2_F_realimag;
% B 
% Real only
B_AW2_F_real_d2p5e15 = Gt_AW2_F_real * d_2p5e15_AW2_F_real_v;
B_AW2_F_real_d7p5e15 = Gt_AW2_F_real * d_7p5e15_AW2_F_real_v;
B_AW2_F_real_d2p5e16 = Gt_AW2_F_real * d_2p5e16_AW2_F_real_v;
B_AW2_F_real_d7p5e16 = Gt_AW2_F_real * d_7p5e16_AW2_F_real_v;
B_AW2_F_real_d2p5e17 = Gt_AW2_F_real * d_2p5e17_AW2_F_real_v;
B_AW2_F_real_d7p5e17 = Gt_AW2_F_real * d_7p5e17_AW2_F_real_v;
B_AW2_F_real_d2p5e18 = Gt_AW2_F_real * d_2p5e18_AW2_F_real_v;
B_AW2_F_real_d7p5e18 = Gt_AW2_F_real * d_7p5e18_AW2_F_real_v;
B_AW2_F_real_d2p5e19 = Gt_AW2_F_real * d_2p5e19_AW2_F_real_v;
B_AW2_F_real_d7p5e19 = Gt_AW2_F_real * d_7p5e19_AW2_F_real_v;
% Real and imag (complex magnitude) 
B_AW2_F_realimag_d2p5e15 = Gt_AW2_F_realimag * d_2p5e15_AW2_F_realimag_v;
B_AW2_F_realimag_d7p5e15 = Gt_AW2_F_realimag * d_7p5e15_AW2_F_realimag_v;
B_AW2_F_realimag_d2p5e16 = Gt_AW2_F_realimag * d_2p5e16_AW2_F_realimag_v;
B_AW2_F_realimag_d7p5e16 = Gt_AW2_F_realimag * d_7p5e16_AW2_F_realimag_v;
B_AW2_F_realimag_d2p5e17 = Gt_AW2_F_realimag * d_2p5e17_AW2_F_realimag_v;
B_AW2_F_realimag_d7p5e17 = Gt_AW2_F_realimag * d_7p5e17_AW2_F_realimag_v;
B_AW2_F_realimag_d2p5e18 = Gt_AW2_F_realimag * d_2p5e18_AW2_F_realimag_v;
B_AW2_F_realimag_d7p5e18 = Gt_AW2_F_realimag * d_7p5e18_AW2_F_realimag_v;
B_AW2_F_realimag_d2p5e19 = Gt_AW2_F_realimag * d_2p5e19_AW2_F_realimag_v;
B_AW2_F_realimag_d7p5e19 = Gt_AW2_F_realimag * d_7p5e19_AW2_F_realimag_v;
% mEst 
% Real only 
mEst_AW2_F_real_2p5e15 = A_AW2_F_real\B_AW2_F_real_d2p5e15;
mEst_AW2_F_real_7p5e15 = A_AW2_F_real\B_AW2_F_real_d7p5e15;
mEst_AW2_F_real_2p5e16 = A_AW2_F_real\B_AW2_F_real_d2p5e16;
mEst_AW2_F_real_7p5e16 = A_AW2_F_real\B_AW2_F_real_d7p5e16;
mEst_AW2_F_real_2p5e17 = A_AW2_F_real\B_AW2_F_real_d2p5e17;
mEst_AW2_F_real_7p5e17 = A_AW2_F_real\B_AW2_F_real_d7p5e17;
mEst_AW2_F_real_2p5e18 = A_AW2_F_real\B_AW2_F_real_d2p5e18;
mEst_AW2_F_real_7p5e18 = A_AW2_F_real\B_AW2_F_real_d7p5e18;
mEst_AW2_F_real_2p5e19 = A_AW2_F_real\B_AW2_F_real_d2p5e19;
mEst_AW2_F_real_7p5e19 = A_AW2_F_real\B_AW2_F_real_d7p5e19;
% Real and imag (complex magnitude) 
mEst_AW2_F_realimag_2p5e15 = A_AW2_F_realimag\B_AW2_F_realimag_d2p5e15;
mEst_AW2_F_realimag_7p5e15 = A_AW2_F_realimag\B_AW2_F_realimag_d7p5e15;
mEst_AW2_F_realimag_2p5e16 = A_AW2_F_realimag\B_AW2_F_realimag_d2p5e16;
mEst_AW2_F_realimag_7p5e16 = A_AW2_F_realimag\B_AW2_F_realimag_d7p5e16;
mEst_AW2_F_realimag_2p5e17 = A_AW2_F_realimag\B_AW2_F_realimag_d2p5e17;
mEst_AW2_F_realimag_7p5e17 = A_AW2_F_realimag\B_AW2_F_realimag_d7p5e17;
mEst_AW2_F_realimag_2p5e18 = A_AW2_F_realimag\B_AW2_F_realimag_d2p5e18;
mEst_AW2_F_realimag_7p5e18 = A_AW2_F_realimag\B_AW2_F_realimag_d7p5e18;
mEst_AW2_F_realimag_2p5e19 = A_AW2_F_realimag\B_AW2_F_realimag_d2p5e19;
mEst_AW2_F_realimag_7p5e19 = A_AW2_F_realimag\B_AW2_F_realimag_d7p5e19;
% Results 
% concentration_SO2 = [1.5e15; 7.5e15; 2.5e16; 7.5e16; 2.5e17; 7.5e17; 2.5e18; 7.5e18; 2.5e19; 7.5e19]; 
mEst_AW2_F_real = [mEst_AW2_F_real_2p5e15, mEst_AW2_F_real_7p5e15, mEst_AW2_F_real_2p5e16, mEst_AW2_F_real_7p5e16, mEst_AW2_F_real_2p5e17, mEst_AW2_F_real_7p5e17, mEst_AW2_F_real_2p5e18, mEst_AW2_F_real_7p5e18, mEst_AW2_F_real_2p5e19, mEst_AW2_F_real_7p5e19];  
mEst_AW2_F_realimag = [mEst_AW2_F_realimag_2p5e15, mEst_AW2_F_realimag_7p5e15, mEst_AW2_F_realimag_2p5e16, mEst_AW2_F_realimag_7p5e16, mEst_AW2_F_realimag_2p5e17, mEst_AW2_F_realimag_7p5e17, mEst_AW2_F_realimag_2p5e18, mEst_AW2_F_realimag_7p5e18, mEst_AW2_F_realimag_2p5e19, mEst_AW2_F_realimag_7p5e19];  

% Finds mEst for AW3_F 
% Takes real and imaginary parts: I0 and SO2
% Real
WT_I0_AW3_F_real = real(WT_I0_AW3_F); 
WT_SO2n_AW3_F_real = real(WT_SO2n_AW3_F); 
% Imaginary 
WT_I0_AW3_F_imag = imag(WT_I0_AW3_F); 
WT_SO2n_AW3_F_imag = imag(WT_SO2n_AW3_F); 
% Vectorises 
% Real
WT_I0_AW3_F_real_v = WT_I0_AW3_F_real(:); 
WT_SO2n_AW3_F_real_v = WT_SO2n_AW3_F_real(:); 
% Imaginary
WT_I0_AW3_F_imag_v = WT_I0_AW3_F_imag(:); 
WT_SO2n_AW3_F_imag_v = WT_SO2n_AW3_F_imag(:); 
% Combines real and complex magnitude 
% I0
WT_I0_AW3_F_realimag = [WT_I0_AW3_F_real_v; WT_I0_AW3_F_imag_v]; 
% SO2
WT_SO2n_AW3_F_realimag = [WT_SO2n_AW3_F_real_v; WT_SO2n_AW3_F_imag_v]; 
% Creates design matrix, G 
G_AW3_F_real = [WT_I0_AW3_F_real_v, WT_SO2n_AW3_F_real_v]; 
G_AW3_F_realimag = [WT_I0_AW3_F_realimag, WT_SO2n_AW3_F_realimag]; 
% Remove NaN 
G_AW3_F_real = G_AW3_F_real(all(~isnan(G_AW3_F_real), 2),:);
G_AW3_F_realimag = G_AW3_F_realimag(all(~isnan(G_AW3_F_realimag), 2),:);
% G transponse 
Gt_AW3_F_real = G_AW3_F_real'; 
Gt_AW3_F_realimag = G_AW3_F_realimag'; 

% Takes real and imaginary parts: I0 and SO2
% Real
WTd_2p5e15_AW3_F_real = real(WTd_2p5e15_AW3_F); 
WTd_7p5e15_AW3_F_real = real(WTd_7p5e15_AW3_F); 
WTd_2p5e16_AW3_F_real = real(WTd_2p5e16_AW3_F); 
WTd_7p5e16_AW3_F_real = real(WTd_7p5e16_AW3_F); 
WTd_2p5e17_AW3_F_real = real(WTd_2p5e17_AW3_F); 
WTd_7p5e17_AW3_F_real = real(WTd_7p5e17_AW3_F); 
WTd_2p5e18_AW3_F_real = real(WTd_2p5e18_AW3_F); 
WTd_7p5e18_AW3_F_real = real(WTd_7p5e18_AW3_F); 
WTd_2p5e19_AW3_F_real = real(WTd_2p5e19_AW3_F); 
WTd_7p5e19_AW3_F_real = real(WTd_7p5e19_AW3_F); 
% Imaginary 
WTd_2p5e15_AW3_F_imag = imag(WTd_2p5e15_AW3_F); 
WTd_7p5e15_AW3_F_imag = imag(WTd_7p5e15_AW3_F); 
WTd_2p5e16_AW3_F_imag = imag(WTd_2p5e16_AW3_F); 
WTd_7p5e16_AW3_F_imag = imag(WTd_7p5e16_AW3_F); 
WTd_2p5e17_AW3_F_imag = imag(WTd_2p5e17_AW3_F); 
WTd_7p5e17_AW3_F_imag = imag(WTd_7p5e17_AW3_F); 
WTd_2p5e18_AW3_F_imag = imag(WTd_2p5e18_AW3_F); 
WTd_7p5e18_AW3_F_imag = imag(WTd_7p5e18_AW3_F); 
WTd_2p5e19_AW3_F_imag = imag(WTd_2p5e19_AW3_F); 
WTd_7p5e19_AW3_F_imag = imag(WTd_7p5e19_AW3_F); 
% Vectorises 
% Real
WTd_2p5e15_AW3_F_real_v = WTd_2p5e15_AW3_F_real(:); 
WTd_7p5e15_AW3_F_real_v = WTd_7p5e15_AW3_F_real(:); 
WTd_2p5e16_AW3_F_real_v = WTd_2p5e16_AW3_F_real(:); 
WTd_7p5e16_AW3_F_real_v = WTd_7p5e16_AW3_F_real(:); 
WTd_2p5e17_AW3_F_real_v = WTd_2p5e17_AW3_F_real(:); 
WTd_7p5e17_AW3_F_real_v = WTd_7p5e17_AW3_F_real(:); 
WTd_2p5e18_AW3_F_real_v = WTd_2p5e18_AW3_F_real(:); 
WTd_7p5e18_AW3_F_real_v = WTd_7p5e18_AW3_F_real(:); 
WTd_2p5e19_AW3_F_real_v = WTd_2p5e19_AW3_F_real(:); 
WTd_7p5e19_AW3_F_real_v = WTd_7p5e19_AW3_F_real(:); 
% Imaginary
WTd_2p5e15_AW3_F_imag_v = WTd_2p5e15_AW3_F_imag(:); 
WTd_7p5e15_AW3_F_imag_v = WTd_7p5e15_AW3_F_imag(:); 
WTd_2p5e16_AW3_F_imag_v = WTd_2p5e16_AW3_F_imag(:); 
WTd_7p5e16_AW3_F_imag_v = WTd_7p5e16_AW3_F_imag(:); 
WTd_2p5e17_AW3_F_imag_v = WTd_2p5e17_AW3_F_imag(:); 
WTd_7p5e17_AW3_F_imag_v = WTd_7p5e17_AW3_F_imag(:); 
WTd_2p5e18_AW3_F_imag_v = WTd_2p5e18_AW3_F_imag(:); 
WTd_7p5e18_AW3_F_imag_v = WTd_7p5e18_AW3_F_imag(:); 
WTd_2p5e19_AW3_F_imag_v = WTd_2p5e19_AW3_F_imag(:); 
WTd_7p5e19_AW3_F_imag_v = WTd_7p5e19_AW3_F_imag(:); 
% Combines real and complex magnitude 
WTd_2p5e15_AW3_F_realimag_v = [WTd_2p5e15_AW3_F_real_v; WTd_2p5e15_AW3_F_imag_v]; 
WTd_7p5e15_AW3_F_realimag_v = [WTd_7p5e15_AW3_F_real_v; WTd_7p5e15_AW3_F_imag_v]; 
WTd_2p5e16_AW3_F_realimag_v = [WTd_2p5e16_AW3_F_real_v; WTd_2p5e16_AW3_F_imag_v]; 
WTd_7p5e16_AW3_F_realimag_v = [WTd_7p5e16_AW3_F_real_v; WTd_7p5e16_AW3_F_imag_v]; 
WTd_2p5e17_AW3_F_realimag_v = [WTd_2p5e17_AW3_F_real_v; WTd_2p5e17_AW3_F_imag_v]; 
WTd_7p5e17_AW3_F_realimag_v = [WTd_7p5e17_AW3_F_real_v; WTd_7p5e17_AW3_F_imag_v]; 
WTd_2p5e18_AW3_F_realimag_v = [WTd_2p5e18_AW3_F_real_v; WTd_2p5e18_AW3_F_imag_v]; 
WTd_7p5e18_AW3_F_realimag_v = [WTd_7p5e18_AW3_F_real_v; WTd_7p5e18_AW3_F_imag_v]; 
WTd_2p5e19_AW3_F_realimag_v = [WTd_2p5e19_AW3_F_real_v; WTd_2p5e19_AW3_F_imag_v]; 
WTd_7p5e19_AW3_F_realimag_v = [WTd_7p5e19_AW3_F_real_v; WTd_7p5e19_AW3_F_imag_v]; 
% Remove NaN 
% Real only 
d_2p5e15_AW3_F_real_v = WTd_2p5e15_AW3_F_real_v(all(~isnan(WTd_2p5e15_AW3_F_real_v), 2),:);
d_7p5e15_AW3_F_real_v = WTd_7p5e15_AW3_F_real_v(all(~isnan(WTd_7p5e15_AW3_F_real_v), 2),:);
d_2p5e16_AW3_F_real_v = WTd_2p5e16_AW3_F_real_v(all(~isnan(WTd_2p5e16_AW3_F_real_v), 2),:);
d_7p5e16_AW3_F_real_v = WTd_7p5e16_AW3_F_real_v(all(~isnan(WTd_7p5e16_AW3_F_real_v), 2),:);
d_2p5e17_AW3_F_real_v = WTd_2p5e17_AW3_F_real_v(all(~isnan(WTd_2p5e17_AW3_F_real_v), 2),:);
d_7p5e17_AW3_F_real_v = WTd_7p5e17_AW3_F_real_v(all(~isnan(WTd_7p5e17_AW3_F_real_v), 2),:);
d_2p5e18_AW3_F_real_v = WTd_2p5e18_AW3_F_real_v(all(~isnan(WTd_2p5e18_AW3_F_real_v), 2),:);
d_7p5e18_AW3_F_real_v = WTd_7p5e18_AW3_F_real_v(all(~isnan(WTd_7p5e18_AW3_F_real_v), 2),:);
d_2p5e19_AW3_F_real_v = WTd_2p5e19_AW3_F_real_v(all(~isnan(WTd_2p5e19_AW3_F_real_v), 2),:);
d_7p5e19_AW3_F_real_v = WTd_7p5e19_AW3_F_real_v(all(~isnan(WTd_7p5e19_AW3_F_real_v), 2),:);
% Real and imag (complex magnitude) 
d_2p5e15_AW3_F_realimag_v = WTd_2p5e15_AW3_F_realimag_v(all(~isnan(WTd_2p5e15_AW3_F_realimag_v), 2),:);
d_7p5e15_AW3_F_realimag_v = WTd_7p5e15_AW3_F_realimag_v(all(~isnan(WTd_7p5e15_AW3_F_realimag_v), 2),:);
d_2p5e16_AW3_F_realimag_v = WTd_2p5e16_AW3_F_realimag_v(all(~isnan(WTd_2p5e16_AW3_F_realimag_v), 2),:);
d_7p5e16_AW3_F_realimag_v = WTd_7p5e16_AW3_F_realimag_v(all(~isnan(WTd_7p5e16_AW3_F_realimag_v), 2),:);
d_2p5e17_AW3_F_realimag_v = WTd_2p5e17_AW3_F_realimag_v(all(~isnan(WTd_2p5e17_AW3_F_realimag_v), 2),:);
d_7p5e17_AW3_F_realimag_v = WTd_7p5e17_AW3_F_realimag_v(all(~isnan(WTd_7p5e17_AW3_F_realimag_v), 2),:);
d_2p5e18_AW3_F_realimag_v = WTd_2p5e18_AW3_F_realimag_v(all(~isnan(WTd_2p5e18_AW3_F_realimag_v), 2),:);
d_7p5e18_AW3_F_realimag_v = WTd_7p5e18_AW3_F_realimag_v(all(~isnan(WTd_7p5e18_AW3_F_realimag_v), 2),:);
d_2p5e19_AW3_F_realimag_v = WTd_2p5e19_AW3_F_realimag_v(all(~isnan(WTd_2p5e19_AW3_F_realimag_v), 2),:);
d_7p5e19_AW3_F_realimag_v = WTd_7p5e19_AW3_F_realimag_v(all(~isnan(WTd_7p5e19_AW3_F_realimag_v), 2),:);
% Find mEst 
% A 
A_AW3_F_real = Gt_AW3_F_real * G_AW3_F_real;
A_AW3_F_realimag = Gt_AW3_F_realimag * G_AW3_F_realimag;
% B 
% Real only
B_AW3_F_real_d2p5e15 = Gt_AW3_F_real * d_2p5e15_AW3_F_real_v;
B_AW3_F_real_d7p5e15 = Gt_AW3_F_real * d_7p5e15_AW3_F_real_v;
B_AW3_F_real_d2p5e16 = Gt_AW3_F_real * d_2p5e16_AW3_F_real_v;
B_AW3_F_real_d7p5e16 = Gt_AW3_F_real * d_7p5e16_AW3_F_real_v;
B_AW3_F_real_d2p5e17 = Gt_AW3_F_real * d_2p5e17_AW3_F_real_v;
B_AW3_F_real_d7p5e17 = Gt_AW3_F_real * d_7p5e17_AW3_F_real_v;
B_AW3_F_real_d2p5e18 = Gt_AW3_F_real * d_2p5e18_AW3_F_real_v;
B_AW3_F_real_d7p5e18 = Gt_AW3_F_real * d_7p5e18_AW3_F_real_v;
B_AW3_F_real_d2p5e19 = Gt_AW3_F_real * d_2p5e19_AW3_F_real_v;
B_AW3_F_real_d7p5e19 = Gt_AW3_F_real * d_7p5e19_AW3_F_real_v;
% Real and imag (complex magnitude) 
B_AW3_F_realimag_d2p5e15 = Gt_AW3_F_realimag * d_2p5e15_AW3_F_realimag_v;
B_AW3_F_realimag_d7p5e15 = Gt_AW3_F_realimag * d_7p5e15_AW3_F_realimag_v;
B_AW3_F_realimag_d2p5e16 = Gt_AW3_F_realimag * d_2p5e16_AW3_F_realimag_v;
B_AW3_F_realimag_d7p5e16 = Gt_AW3_F_realimag * d_7p5e16_AW3_F_realimag_v;
B_AW3_F_realimag_d2p5e17 = Gt_AW3_F_realimag * d_2p5e17_AW3_F_realimag_v;
B_AW3_F_realimag_d7p5e17 = Gt_AW3_F_realimag * d_7p5e17_AW3_F_realimag_v;
B_AW3_F_realimag_d2p5e18 = Gt_AW3_F_realimag * d_2p5e18_AW3_F_realimag_v;
B_AW3_F_realimag_d7p5e18 = Gt_AW3_F_realimag * d_7p5e18_AW3_F_realimag_v;
B_AW3_F_realimag_d2p5e19 = Gt_AW3_F_realimag * d_2p5e19_AW3_F_realimag_v;
B_AW3_F_realimag_d7p5e19 = Gt_AW3_F_realimag * d_7p5e19_AW3_F_realimag_v;
% mEst 
% Real only 
mEst_AW3_F_real_2p5e15 = A_AW3_F_real\B_AW3_F_real_d2p5e15;
mEst_AW3_F_real_7p5e15 = A_AW3_F_real\B_AW3_F_real_d7p5e15;
mEst_AW3_F_real_2p5e16 = A_AW3_F_real\B_AW3_F_real_d2p5e16;
mEst_AW3_F_real_7p5e16 = A_AW3_F_real\B_AW3_F_real_d7p5e16;
mEst_AW3_F_real_2p5e17 = A_AW3_F_real\B_AW3_F_real_d2p5e17;
mEst_AW3_F_real_7p5e17 = A_AW3_F_real\B_AW3_F_real_d7p5e17;
mEst_AW3_F_real_2p5e18 = A_AW3_F_real\B_AW3_F_real_d2p5e18;
mEst_AW3_F_real_7p5e18 = A_AW3_F_real\B_AW3_F_real_d7p5e18;
mEst_AW3_F_real_2p5e19 = A_AW3_F_real\B_AW3_F_real_d2p5e19;
mEst_AW3_F_real_7p5e19 = A_AW3_F_real\B_AW3_F_real_d7p5e19;
% Real and imag (complex magnitude) 
mEst_AW3_F_realimag_2p5e15 = A_AW3_F_realimag\B_AW3_F_realimag_d2p5e15;
mEst_AW3_F_realimag_7p5e15 = A_AW3_F_realimag\B_AW3_F_realimag_d7p5e15;
mEst_AW3_F_realimag_2p5e16 = A_AW3_F_realimag\B_AW3_F_realimag_d2p5e16;
mEst_AW3_F_realimag_7p5e16 = A_AW3_F_realimag\B_AW3_F_realimag_d7p5e16;
mEst_AW3_F_realimag_2p5e17 = A_AW3_F_realimag\B_AW3_F_realimag_d2p5e17;
mEst_AW3_F_realimag_7p5e17 = A_AW3_F_realimag\B_AW3_F_realimag_d7p5e17;
mEst_AW3_F_realimag_2p5e18 = A_AW3_F_realimag\B_AW3_F_realimag_d2p5e18;
mEst_AW3_F_realimag_7p5e18 = A_AW3_F_realimag\B_AW3_F_realimag_d7p5e18;
mEst_AW3_F_realimag_2p5e19 = A_AW3_F_realimag\B_AW3_F_realimag_d2p5e19;
mEst_AW3_F_realimag_7p5e19 = A_AW3_F_realimag\B_AW3_F_realimag_d7p5e19;
% Results 
% concentration_SO2 = [1.5e15; 7.5e15; 2.5e16; 7.5e16; 2.5e17; 7.5e17; 2.5e18; 7.5e18; 2.5e19; 7.5e19]; 
mEst_AW3_F_real = [mEst_AW3_F_real_2p5e15, mEst_AW3_F_real_7p5e15, mEst_AW3_F_real_2p5e16, mEst_AW3_F_real_7p5e16, mEst_AW3_F_real_2p5e17, mEst_AW3_F_real_7p5e17, mEst_AW3_F_real_2p5e18, mEst_AW3_F_real_7p5e18, mEst_AW3_F_real_2p5e19, mEst_AW3_F_real_7p5e19];  
mEst_AW3_F_realimag = [mEst_AW3_F_realimag_2p5e15, mEst_AW3_F_realimag_7p5e15, mEst_AW3_F_realimag_2p5e16, mEst_AW3_F_realimag_7p5e16, mEst_AW3_F_realimag_2p5e17, mEst_AW3_F_realimag_7p5e17, mEst_AW3_F_realimag_2p5e18, mEst_AW3_F_realimag_7p5e18, mEst_AW3_F_realimag_2p5e19, mEst_AW3_F_realimag_7p5e19];  
% Plot 
figure % Figure 9  
subplot(2,1,1)
% AW1_F
p = plot(concentration_SO2, mEst_AW1_F_real(2,:), '^'); hold on; p.Color = AW1col; p.MarkerSize = 6; p.MarkerFaceColor = AW1col; 
p = plot(concentration_SO2, mEst_AW1_F_realimag(2,:), 'v'); p.Color = AW1col; p.MarkerSize = 6; p.MarkerFaceColor = AW1col; 
% AW2_F
p = plot(concentration_SO2, mEst_AW2_F_real(2,:), '^'); hold on; p.Color = AW2col; p.MarkerSize = 6; p.MarkerFaceColor = AW2col; 
p = plot(concentration_SO2, mEst_AW2_F_realimag(2,:), 'v'); p.Color = AW2col; p.MarkerSize = 6; p.MarkerFaceColor = AW2col; 
% AW3_F
p = plot(concentration_SO2, mEst_AW3_F_real(2,:), '^'); hold on; p.Color = AW3col; p.MarkerSize = 6; p.MarkerFaceColor = AW3col; 
p = plot(concentration_SO2, mEst_AW3_F_realimag(2,:), 'v'); p.Color = AW3col; p.MarkerSize = 6; p.MarkerFaceColor = AW3col; 
%
plot(concentration_SO2, concentration_SO2, ':k');
ylabel('mEst_S_O_2 (molec/cm^2)')
xlabel('SO_2 column density (molec/cm^2)')
set(gcf,'color','w');
fig = gca; fig.FontSize = 14; 
xlim([-0.1e19 8e19])
ylim([-0.5e19 8e19])
subplot(2,1,2)
% AW1_F
p = plot(concentration_SO2, mEst_AW1_F_real(2,:), '^'); hold on; p.Color = AW1col; p.MarkerSize = 6; p.MarkerFaceColor = AW1col; 
p = plot(concentration_SO2, mEst_AW1_F_realimag(2,:), 'v'); p.Color = AW1col; p.MarkerSize = 6; p.MarkerFaceColor = AW1col; 
% AW2_F
p = plot(concentration_SO2, mEst_AW2_F_real(2,:), '^'); hold on; p.Color = AW2col; p.MarkerSize = 6; p.MarkerFaceColor = AW2col; 
p = plot(concentration_SO2, mEst_AW2_F_realimag(2,:), 'v'); p.Color = AW2col; p.MarkerSize = 6; p.MarkerFaceColor = AW2col; 
% AW3_F
p = plot(concentration_SO2, mEst_AW3_F_real(2,:), '^'); hold on; p.Color = AW3col; p.MarkerSize = 6; p.MarkerFaceColor = AW3col; 
p = plot(concentration_SO2, mEst_AW3_F_realimag(2,:), 'v'); p.Color = AW3col; p.MarkerSize = 6; p.MarkerFaceColor = AW3col; 
%
plot(concentration_SO2, concentration_SO2, ':k');
ylabel('mEst_S_O_2 (molec/cm^2)')
xlabel('SO_2 column density (molec/cm^2)')
set(gcf,'color','w');
fig = gca; fig.FontSize = 14; 
xlim([-0.3e18 8e18])
ylim([-2e18 8e18])
% Save convoluted cross sections to be used later 
fname = fullfile(outDir, 'mEst_real_imag_analysis_window_F');
save(fname); % Saves AW1_F variables  
%% Make any adjustments 
% SAVES FIGURE
% fnOut = 'Fig11_mEst_real_imag_AW_F';
% fname = fullfile(outDirFig, fnOut);
% saveas(gca, fname)













