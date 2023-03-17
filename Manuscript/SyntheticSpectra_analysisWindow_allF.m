% Synthetic spectra - analysis window (all frequencies)

% Must run SyntheticSpectra_AllWavelengths.m first to find WT of d and G
% This script extracts the data from that plotted in Fig. 6 
 
outDir = '/Users/charlotteb/Documents/Chapter 2/Synthetic spectra/Analysis window all frequency Fig 8 and Fig 9/'; 
outDirFig = '/Users/charlotteb/Documents/Chapter 2/Synthetic spectra/Analysis window all frequency Fig 8 and Fig 9/'; 

% Extract data for analysis window 1 - AW 1
% Index (channel) corresponding to AW1_l1 and AW1_l2 
% With reference to lambda (defined in SyntheticSpectra_AllWavelengths.m)
% AW1_l1 (all defined in SyntheticSpectra_AllWavelengths.m)
% AW1_l2 (all defined in SyntheticSpectra_AllWavelengths.m)
indAW1_l1 = 6913; 
indAW1_l2 = 12606; 
% G 
WT_I0_AW1 = WT_I0(:, indAW1_l1:indAW1_l2); % All frequencies (rows), extract wavlengths (columns)
WT_SO2n_AW1 = WT_SO2n(:, indAW1_l1:indAW1_l2);
% d
WTd_2p5e15_AW1 = WTd_2p5e15(:, indAW1_l1:indAW1_l2);
WTd_7p5e15_AW1 = WTd_7p5e15(:, indAW1_l1:indAW1_l2); 
WTd_2p5e16_AW1 = WTd_2p5e16(:, indAW1_l1:indAW1_l2);
WTd_7p5e16_AW1 = WTd_7p5e16(:, indAW1_l1:indAW1_l2);
WTd_2p5e17_AW1 = WTd_2p5e17(:, indAW1_l1:indAW1_l2);
WTd_7p5e17_AW1 = WTd_7p5e17(:, indAW1_l1:indAW1_l2);
WTd_2p5e18_AW1 = WTd_2p5e18(:, indAW1_l1:indAW1_l2);
WTd_7p5e18_AW1 = WTd_7p5e18(:, indAW1_l1:indAW1_l2);
WTd_2p5e19_AW1 = WTd_2p5e19(:, indAW1_l1:indAW1_l2);
WTd_7p5e19_AW1 = WTd_7p5e19(:, indAW1_l1:indAW1_l2);
% 
lambda_AW1 = lambda(indAW1_l1:indAW1_l2);

% Repeat for AW 2 and 3
% Index (channel) corresponding to AW2_l1 and AW2_l2 
% With reference to lambda (defined in SyntheticSpectra_AllWavelengths.m)
% AW2
% AW2_l1 (all defined in SyntheticSpectra_AllWavelengths.m)
% AW2_l2 (all defined in SyntheticSpectra_AllWavelengths.m)
indAW2_l1 = 8929; 
indAW2_l2 = 12606; 
% G 
WT_I0_AW2 = WT_I0(:, indAW2_l1:indAW2_l2); % All frequencies (rows), extract wavlengths (columns)
WT_SO2n_AW2 = WT_SO2n(:, indAW2_l1:indAW2_l2);
% d
WTd_2p5e15_AW2 = WTd_2p5e15(:, indAW2_l1:indAW2_l2);
WTd_7p5e15_AW2 = WTd_7p5e15(:, indAW2_l1:indAW2_l2); 
WTd_2p5e16_AW2 = WTd_2p5e16(:, indAW2_l1:indAW2_l2);
WTd_7p5e16_AW2 = WTd_7p5e16(:, indAW2_l1:indAW2_l2);
WTd_2p5e17_AW2 = WTd_2p5e17(:, indAW2_l1:indAW2_l2);
WTd_7p5e17_AW2 = WTd_7p5e17(:, indAW2_l1:indAW2_l2);
WTd_2p5e18_AW2 = WTd_2p5e18(:, indAW2_l1:indAW2_l2);
WTd_7p5e18_AW2 = WTd_7p5e18(:, indAW2_l1:indAW2_l2);
WTd_2p5e19_AW2 = WTd_2p5e19(:, indAW2_l1:indAW2_l2);
WTd_7p5e19_AW2 = WTd_7p5e19(:, indAW2_l1:indAW2_l2);
% 
lambda_AW2 = lambda(indAW2_l1:indAW2_l2);
% AW3
% Index (channel) corresponding to AW3_l1 and AW3_l2 
% With reference to lambda (defined in SyntheticSpectra_AllWavelengths.m)
% AW3_l1 (all defined in SyntheticSpectra_AllWavelengths.m)
% AW3_l2 (all defined in SyntheticSpectra_AllWavelengths.m)
indAW3_l1 = 7735; 
indAW3_l2 = 10823; 
% G 
WT_I0_AW3 = WT_I0(:, indAW3_l1:indAW3_l2); % All frequencies (rows), extract wavlengths (columns)
WT_SO2n_AW3 = WT_SO2n(:, indAW3_l1:indAW3_l2);
% d
WTd_2p5e15_AW3 = WTd_2p5e15(:, indAW3_l1:indAW3_l2);
WTd_7p5e15_AW3 = WTd_7p5e15(:, indAW3_l1:indAW3_l2); 
WTd_2p5e16_AW3 = WTd_2p5e16(:, indAW3_l1:indAW3_l2);
WTd_7p5e16_AW3 = WTd_7p5e16(:, indAW3_l1:indAW3_l2);
WTd_2p5e17_AW3 = WTd_2p5e17(:, indAW3_l1:indAW3_l2);
WTd_7p5e17_AW3 = WTd_7p5e17(:, indAW3_l1:indAW3_l2);
WTd_2p5e18_AW3 = WTd_2p5e18(:, indAW3_l1:indAW3_l2);
WTd_7p5e18_AW3 = WTd_7p5e18(:, indAW3_l1:indAW3_l2);
WTd_2p5e19_AW3 = WTd_2p5e19(:, indAW3_l1:indAW3_l2);
WTd_7p5e19_AW3 = WTd_7p5e19(:, indAW3_l1:indAW3_l2);
% 
lambda_AW3 = lambda(indAW3_l1:indAW3_l2);

% Figure 8 - plot linear model for analysis windows 
figure 
% Subplots 
% First column of figure is d (plot for 2.5e18) 
subplot(3,3,1) 
pcolor(lambda_AW1, F, real(WTd_2p5e18_AW1)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
xlim([min(lambda_AW1) max(lambda_AW1)]); set(gca,'YScale','log'); % ylabel('Spatial frequency (cycles/\lambda)', 'FontSize', 10);
% cb = colorbar; cb.Label.Position(1) = 3.5;
% ylabel(cb,'Amplitude', 'Rotation', -90, 'FontSize', 10);
set(gcf,'color','w');
title('d')
subtitle('J(\lambda) prime')
subplot(3,3,4) 
pcolor(lambda_AW2, F, real(WTd_2p5e18_AW2)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
xlim([min(lambda_AW2) max(lambda_AW2)]); set(gca,'YScale','log'); ylabel('Spatial frequency (cycles/\lambda)', 'FontSize', 14);
% cb = colorbar; cb.Label.Position(1) = 3.5;
% ylabel(cb,'Amplitude', 'Rotation', -90, 'FontSize', 10);
set(gcf,'color','w');
subplot(3,3,7) 
pcolor(lambda_AW3, F, real(WTd_2p5e18_AW3)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
xlim([min(lambda_AW3) max(lambda_AW3)]); set(gca,'YScale','log'); % ylabel('Spatial frequency (cycles/\lambda)', 'FontSize', 10);
% cb = colorbar; cb.Label.Position(1) = 3.5;
% ylabel(cb,'Amplitude', 'Rotation', -90, 'FontSize', 10);
set(gcf,'color','w');
subplot(3,3,2) 
pcolor(lambda_AW1, F, real(WT_I0_AW1)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
xlim([min(lambda_AW1) max(lambda_AW1)]); set(gca,'YScale','log'); % ylabel('Spatial frequency (cycles/\lambda)', 'FontSize', 10);
% cb = colorbar; cb.Label.Position(1) = 3.5;
% ylabel(cb,'Amplitude', 'Rotation', -90, 'FontSize', 10);
set(gcf,'color','w');
title('G') 
subtitle('J_0(\lambda) prime')
subplot(3,3,5) 
pcolor(lambda_AW2, F, real(WT_I0_AW2)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
xlim([min(lambda_AW2) max(lambda_AW2)]); set(gca,'YScale','log'); % ylabel('Spatial frequency (cycles/\lambda)', 'FontSize', 10);
% cb = colorbar; cb.Label.Position(1) = 3.5;
% ylabel(cb,'Amplitude', 'Rotation', -90, 'FontSize', 10);
set(gcf,'color','w');
subplot(3,3,8) 
pcolor(lambda_AW3, F, real(WT_I0_AW3)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
xlim([min(lambda_AW3) max(lambda_AW3)]); set(gca,'YScale','log'); % ylabel('Spatial frequency (cycles/\lambda)', 'FontSize', 10);
% cb = colorbar; cb.Label.Position(1) = 3.5;
% ylabel(cb,'Amplitude', 'Rotation', -90, 'FontSize', 10);
xlabel('\lambda (nm)', 'FontSize', 14); 
set(gcf,'color','w');
subplot(3,3,3) 
pcolor(lambda_AW1, F, real(WT_SO2n_AW1)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
xlim([min(lambda_AW1) max(lambda_AW1)]); set(gca,'YScale','log'); % ylabel('Spatial frequency (cycles/\lambda)', 'FontSize', 10);
% cb = colorbar; cb.Label.Position(1) = 3.5;
% ylabel(cb,'Amplitude', 'Rotation', -90, 'FontSize', 10);
subtitle('\sigma_S_O_2 prime')
set(gcf,'color','w');
subplot(3,3,6) 
pcolor(lambda_AW2, F, real(WT_SO2n_AW2)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
xlim([min(lambda_AW2) max(lambda_AW2)]); set(gca,'YScale','log'); % ylabel('Spatial frequency (cycles/\lambda)', 'FontSize', 10);
cb = colorbar; cb.Label.Position(1) = 3.5;
ylabel(cb,'Amplitude', 'Rotation', -90, 'FontSize', 14);
set(gcf,'color','w');
subplot(3,3,9) 
pcolor(lambda_AW3, F, real(WT_SO2n_AW3)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
xlim([min(lambda_AW3) max(lambda_AW3)]); set(gca,'YScale','log'); % ylabel('Spatial frequency (cycles/\lambda)', 'FontSize', 10); 
% cb = colorbar; cb.Label.Position(1) = 3.5;
% ylabel(cb,'Amplitude', 'Rotation', -90, 'FontSize', 10);
set(gcf,'color','w');
% Save 
fname = fullfile(outDir, 'mEst_real_imag_AW');
save(fname); % Saves all variables  
%% Make any adjustments 
% SAVES FIGURE
% fnOut = 'Fig8_dG_AW';
% fname = fullfile(outDirFig, fnOut);
% saveas(gca, fname)

% Finds mEst for AW1 
% Takes real and imaginary parts: I0 and SO2
% Real
WT_I0_AW1_real = real(WT_I0_AW1); 
WT_SO2n_AW1_real = real(WT_SO2n_AW1); 
% Imaginary 
WT_I0_AW1_imag = imag(WT_I0_AW1); 
WT_SO2n_AW1_imag = imag(WT_SO2n_AW1); 
% Vectorises 
% Real
WT_I0_AW1_real_v = WT_I0_AW1_real(:); 
WT_SO2n_AW1_real_v = WT_SO2n_AW1_real(:); 
% Imaginary
WT_I0_AW1_imag_v = WT_I0_AW1_imag(:); 
WT_SO2n_AW1_imag_v = WT_SO2n_AW1_imag(:); 
% Combines real and complex magnitude 
% I0
WT_I0_AW1_realimag = [WT_I0_AW1_real_v; WT_I0_AW1_imag_v]; 
% SO2
WT_SO2n_AW1_realimag = [WT_SO2n_AW1_real_v; WT_SO2n_AW1_imag_v]; 
% Creates design matrix, G 
G_AW1_real = [WT_I0_AW1_real_v, WT_SO2n_AW1_real_v]; 
G_AW1_realimag = [WT_I0_AW1_realimag, WT_SO2n_AW1_realimag]; 
% Remove NaN 
G_AW1_real = G_AW1_real(all(~isnan(G_AW1_real), 2),:);
G_AW1_realimag = G_AW1_realimag(all(~isnan(G_AW1_realimag), 2),:);
% G transponse 
Gt_AW1_real = G_AW1_real'; 
Gt_AW1_realimag = G_AW1_realimag'; 

% Takes real and imaginary parts: I0 and SO2
% Real
WTd_2p5e15_AW1_real = real(WTd_2p5e15_AW1); 
WTd_7p5e15_AW1_real = real(WTd_7p5e15_AW1); 
WTd_2p5e16_AW1_real = real(WTd_2p5e16_AW1); 
WTd_7p5e16_AW1_real = real(WTd_7p5e16_AW1); 
WTd_2p5e17_AW1_real = real(WTd_2p5e17_AW1); 
WTd_7p5e17_AW1_real = real(WTd_7p5e17_AW1); 
WTd_2p5e18_AW1_real = real(WTd_2p5e18_AW1); 
WTd_7p5e18_AW1_real = real(WTd_7p5e18_AW1); 
WTd_2p5e19_AW1_real = real(WTd_2p5e19_AW1); 
WTd_7p5e19_AW1_real = real(WTd_7p5e19_AW1); 
% Imaginary 
WTd_2p5e15_AW1_imag = imag(WTd_2p5e15_AW1); 
WTd_7p5e15_AW1_imag = imag(WTd_7p5e15_AW1); 
WTd_2p5e16_AW1_imag = imag(WTd_2p5e16_AW1); 
WTd_7p5e16_AW1_imag = imag(WTd_7p5e16_AW1); 
WTd_2p5e17_AW1_imag = imag(WTd_2p5e17_AW1); 
WTd_7p5e17_AW1_imag = imag(WTd_7p5e17_AW1); 
WTd_2p5e18_AW1_imag = imag(WTd_2p5e18_AW1); 
WTd_7p5e18_AW1_imag = imag(WTd_7p5e18_AW1); 
WTd_2p5e19_AW1_imag = imag(WTd_2p5e19_AW1); 
WTd_7p5e19_AW1_imag = imag(WTd_7p5e19_AW1); 
% Vectorises 
% Real
WTd_2p5e15_AW1_real_v = WTd_2p5e15_AW1_real(:); 
WTd_7p5e15_AW1_real_v = WTd_7p5e15_AW1_real(:); 
WTd_2p5e16_AW1_real_v = WTd_2p5e16_AW1_real(:); 
WTd_7p5e16_AW1_real_v = WTd_7p5e16_AW1_real(:); 
WTd_2p5e17_AW1_real_v = WTd_2p5e17_AW1_real(:); 
WTd_7p5e17_AW1_real_v = WTd_7p5e17_AW1_real(:); 
WTd_2p5e18_AW1_real_v = WTd_2p5e18_AW1_real(:); 
WTd_7p5e18_AW1_real_v = WTd_7p5e18_AW1_real(:); 
WTd_2p5e19_AW1_real_v = WTd_2p5e19_AW1_real(:); 
WTd_7p5e19_AW1_real_v = WTd_7p5e19_AW1_real(:); 
% Imaginary
WTd_2p5e15_AW1_imag_v = WTd_2p5e15_AW1_imag(:); 
WTd_7p5e15_AW1_imag_v = WTd_7p5e15_AW1_imag(:); 
WTd_2p5e16_AW1_imag_v = WTd_2p5e16_AW1_imag(:); 
WTd_7p5e16_AW1_imag_v = WTd_7p5e16_AW1_imag(:); 
WTd_2p5e17_AW1_imag_v = WTd_2p5e17_AW1_imag(:); 
WTd_7p5e17_AW1_imag_v = WTd_7p5e17_AW1_imag(:); 
WTd_2p5e18_AW1_imag_v = WTd_2p5e18_AW1_imag(:); 
WTd_7p5e18_AW1_imag_v = WTd_7p5e18_AW1_imag(:); 
WTd_2p5e19_AW1_imag_v = WTd_2p5e19_AW1_imag(:); 
WTd_7p5e19_AW1_imag_v = WTd_7p5e19_AW1_imag(:); 
% Combines real and complex magnitude 
WTd_2p5e15_AW1_realimag_v = [WTd_2p5e15_AW1_real_v; WTd_2p5e15_AW1_imag_v]; 
WTd_7p5e15_AW1_realimag_v = [WTd_7p5e15_AW1_real_v; WTd_7p5e15_AW1_imag_v]; 
WTd_2p5e16_AW1_realimag_v = [WTd_2p5e16_AW1_real_v; WTd_2p5e16_AW1_imag_v]; 
WTd_7p5e16_AW1_realimag_v = [WTd_7p5e16_AW1_real_v; WTd_7p5e16_AW1_imag_v]; 
WTd_2p5e17_AW1_realimag_v = [WTd_2p5e17_AW1_real_v; WTd_2p5e17_AW1_imag_v]; 
WTd_7p5e17_AW1_realimag_v = [WTd_7p5e17_AW1_real_v; WTd_7p5e17_AW1_imag_v]; 
WTd_2p5e18_AW1_realimag_v = [WTd_2p5e18_AW1_real_v; WTd_2p5e18_AW1_imag_v]; 
WTd_7p5e18_AW1_realimag_v = [WTd_7p5e18_AW1_real_v; WTd_7p5e18_AW1_imag_v]; 
WTd_2p5e19_AW1_realimag_v = [WTd_2p5e19_AW1_real_v; WTd_2p5e19_AW1_imag_v]; 
WTd_7p5e19_AW1_realimag_v = [WTd_7p5e19_AW1_real_v; WTd_7p5e19_AW1_imag_v]; 
% Remove NaN 
% Real only 
d_2p5e15_AW1_real_v = WTd_2p5e15_AW1_real_v(all(~isnan(WTd_2p5e15_AW1_real_v), 2),:);
d_7p5e15_AW1_real_v = WTd_7p5e15_AW1_real_v(all(~isnan(WTd_7p5e15_AW1_real_v), 2),:);
d_2p5e16_AW1_real_v = WTd_2p5e16_AW1_real_v(all(~isnan(WTd_2p5e16_AW1_real_v), 2),:);
d_7p5e16_AW1_real_v = WTd_7p5e16_AW1_real_v(all(~isnan(WTd_7p5e16_AW1_real_v), 2),:);
d_2p5e17_AW1_real_v = WTd_2p5e17_AW1_real_v(all(~isnan(WTd_2p5e17_AW1_real_v), 2),:);
d_7p5e17_AW1_real_v = WTd_7p5e17_AW1_real_v(all(~isnan(WTd_7p5e17_AW1_real_v), 2),:);
d_2p5e18_AW1_real_v = WTd_2p5e18_AW1_real_v(all(~isnan(WTd_2p5e18_AW1_real_v), 2),:);
d_7p5e18_AW1_real_v = WTd_7p5e18_AW1_real_v(all(~isnan(WTd_7p5e18_AW1_real_v), 2),:);
d_2p5e19_AW1_real_v = WTd_2p5e19_AW1_real_v(all(~isnan(WTd_2p5e19_AW1_real_v), 2),:);
d_7p5e19_AW1_real_v = WTd_7p5e19_AW1_real_v(all(~isnan(WTd_7p5e19_AW1_real_v), 2),:);
% Real and imag (complex magnitude) 
d_2p5e15_AW1_realimag_v = WTd_2p5e15_AW1_realimag_v(all(~isnan(WTd_2p5e15_AW1_realimag_v), 2),:);
d_7p5e15_AW1_realimag_v = WTd_7p5e15_AW1_realimag_v(all(~isnan(WTd_7p5e15_AW1_realimag_v), 2),:);
d_2p5e16_AW1_realimag_v = WTd_2p5e16_AW1_realimag_v(all(~isnan(WTd_2p5e16_AW1_realimag_v), 2),:);
d_7p5e16_AW1_realimag_v = WTd_7p5e16_AW1_realimag_v(all(~isnan(WTd_7p5e16_AW1_realimag_v), 2),:);
d_2p5e17_AW1_realimag_v = WTd_2p5e17_AW1_realimag_v(all(~isnan(WTd_2p5e17_AW1_realimag_v), 2),:);
d_7p5e17_AW1_realimag_v = WTd_7p5e17_AW1_realimag_v(all(~isnan(WTd_7p5e17_AW1_realimag_v), 2),:);
d_2p5e18_AW1_realimag_v = WTd_2p5e18_AW1_realimag_v(all(~isnan(WTd_2p5e18_AW1_realimag_v), 2),:);
d_7p5e18_AW1_realimag_v = WTd_7p5e18_AW1_realimag_v(all(~isnan(WTd_7p5e18_AW1_realimag_v), 2),:);
d_2p5e19_AW1_realimag_v = WTd_2p5e19_AW1_realimag_v(all(~isnan(WTd_2p5e19_AW1_realimag_v), 2),:);
d_7p5e19_AW1_realimag_v = WTd_7p5e19_AW1_realimag_v(all(~isnan(WTd_7p5e19_AW1_realimag_v), 2),:);
% Find mEst 
% A 
A_AW1_real = Gt_AW1_real * G_AW1_real;
A_AW1_realimag = Gt_AW1_realimag * G_AW1_realimag;
% B 
% Real only
B_AW1_real_d2p5e15 = Gt_AW1_real * d_2p5e15_AW1_real_v;
B_AW1_real_d7p5e15 = Gt_AW1_real * d_7p5e15_AW1_real_v;
B_AW1_real_d2p5e16 = Gt_AW1_real * d_2p5e16_AW1_real_v;
B_AW1_real_d7p5e16 = Gt_AW1_real * d_7p5e16_AW1_real_v;
B_AW1_real_d2p5e17 = Gt_AW1_real * d_2p5e17_AW1_real_v;
B_AW1_real_d7p5e17 = Gt_AW1_real * d_7p5e17_AW1_real_v;
B_AW1_real_d2p5e18 = Gt_AW1_real * d_2p5e18_AW1_real_v;
B_AW1_real_d7p5e18 = Gt_AW1_real * d_7p5e18_AW1_real_v;
B_AW1_real_d2p5e19 = Gt_AW1_real * d_2p5e19_AW1_real_v;
B_AW1_real_d7p5e19 = Gt_AW1_real * d_7p5e19_AW1_real_v;
% Real and imag (complex magnitude) 
B_AW1_realimag_d2p5e15 = Gt_AW1_realimag * d_2p5e15_AW1_realimag_v;
B_AW1_realimag_d7p5e15 = Gt_AW1_realimag * d_7p5e15_AW1_realimag_v;
B_AW1_realimag_d2p5e16 = Gt_AW1_realimag * d_2p5e16_AW1_realimag_v;
B_AW1_realimag_d7p5e16 = Gt_AW1_realimag * d_7p5e16_AW1_realimag_v;
B_AW1_realimag_d2p5e17 = Gt_AW1_realimag * d_2p5e17_AW1_realimag_v;
B_AW1_realimag_d7p5e17 = Gt_AW1_realimag * d_7p5e17_AW1_realimag_v;
B_AW1_realimag_d2p5e18 = Gt_AW1_realimag * d_2p5e18_AW1_realimag_v;
B_AW1_realimag_d7p5e18 = Gt_AW1_realimag * d_7p5e18_AW1_realimag_v;
B_AW1_realimag_d2p5e19 = Gt_AW1_realimag * d_2p5e19_AW1_realimag_v;
B_AW1_realimag_d7p5e19 = Gt_AW1_realimag * d_7p5e19_AW1_realimag_v;
% mEst 
% Real only 
mEst_AW1_real_2p5e15 = A_AW1_real\B_AW1_real_d2p5e15;
mEst_AW1_real_7p5e15 = A_AW1_real\B_AW1_real_d7p5e15;
mEst_AW1_real_2p5e16 = A_AW1_real\B_AW1_real_d2p5e16;
mEst_AW1_real_7p5e16 = A_AW1_real\B_AW1_real_d7p5e16;
mEst_AW1_real_2p5e17 = A_AW1_real\B_AW1_real_d2p5e17;
mEst_AW1_real_7p5e17 = A_AW1_real\B_AW1_real_d7p5e17;
mEst_AW1_real_2p5e18 = A_AW1_real\B_AW1_real_d2p5e18;
mEst_AW1_real_7p5e18 = A_AW1_real\B_AW1_real_d7p5e18;
mEst_AW1_real_2p5e19 = A_AW1_real\B_AW1_real_d2p5e19;
mEst_AW1_real_7p5e19 = A_AW1_real\B_AW1_real_d7p5e19;
% Real and imag (complex magnitude) 
mEst_AW1_realimag_2p5e15 = A_AW1_realimag\B_AW1_realimag_d2p5e15;
mEst_AW1_realimag_7p5e15 = A_AW1_realimag\B_AW1_realimag_d7p5e15;
mEst_AW1_realimag_2p5e16 = A_AW1_realimag\B_AW1_realimag_d2p5e16;
mEst_AW1_realimag_7p5e16 = A_AW1_realimag\B_AW1_realimag_d7p5e16;
mEst_AW1_realimag_2p5e17 = A_AW1_realimag\B_AW1_realimag_d2p5e17;
mEst_AW1_realimag_7p5e17 = A_AW1_realimag\B_AW1_realimag_d7p5e17;
mEst_AW1_realimag_2p5e18 = A_AW1_realimag\B_AW1_realimag_d2p5e18;
mEst_AW1_realimag_7p5e18 = A_AW1_realimag\B_AW1_realimag_d7p5e18;
mEst_AW1_realimag_2p5e19 = A_AW1_realimag\B_AW1_realimag_d2p5e19;
mEst_AW1_realimag_7p5e19 = A_AW1_realimag\B_AW1_realimag_d7p5e19;
% Results 
% concentration_SO2 = [1.5e15; 7.5e15; 2.5e16; 7.5e16; 2.5e17; 7.5e17; 2.5e18; 7.5e18; 2.5e19; 7.5e19]; 
mEst_AW1_real = [mEst_AW1_real_2p5e15, mEst_AW1_real_7p5e15, mEst_AW1_real_2p5e16, mEst_AW1_real_7p5e16, mEst_AW1_real_2p5e17, mEst_AW1_real_7p5e17, mEst_AW1_real_2p5e18, mEst_AW1_real_7p5e18, mEst_AW1_real_2p5e19, mEst_AW1_real_7p5e19];  
mEst_AW1_realimag = [mEst_AW1_realimag_2p5e15, mEst_AW1_realimag_7p5e15, mEst_AW1_realimag_2p5e16, mEst_AW1_realimag_7p5e16, mEst_AW1_realimag_2p5e17, mEst_AW1_realimag_7p5e17, mEst_AW1_realimag_2p5e18, mEst_AW1_realimag_7p5e18, mEst_AW1_realimag_2p5e19, mEst_AW1_realimag_7p5e19];  

% Finds mEst for AW2 
% Takes real and imaginary parts: I0 and SO2
% Real
WT_I0_AW2_real = real(WT_I0_AW2); 
WT_SO2n_AW2_real = real(WT_SO2n_AW2); 
% Imaginary 
WT_I0_AW2_imag = imag(WT_I0_AW2); 
WT_SO2n_AW2_imag = imag(WT_SO2n_AW2); 
% Vectorises 
% Real
WT_I0_AW2_real_v = WT_I0_AW2_real(:); 
WT_SO2n_AW2_real_v = WT_SO2n_AW2_real(:); 
% Imaginary
WT_I0_AW2_imag_v = WT_I0_AW2_imag(:); 
WT_SO2n_AW2_imag_v = WT_SO2n_AW2_imag(:); 
% Combines real and complex magnitude 
% I0
WT_I0_AW2_realimag = [WT_I0_AW2_real_v; WT_I0_AW2_imag_v]; 
% SO2
WT_SO2n_AW2_realimag = [WT_SO2n_AW2_real_v; WT_SO2n_AW2_imag_v]; 
% Creates design matrix, G 
G_AW2_real = [WT_I0_AW2_real_v, WT_SO2n_AW2_real_v]; 
G_AW2_realimag = [WT_I0_AW2_realimag, WT_SO2n_AW2_realimag]; 
% Remove NaN 
G_AW2_real = G_AW2_real(all(~isnan(G_AW2_real), 2),:);
G_AW2_realimag = G_AW2_realimag(all(~isnan(G_AW2_realimag), 2),:);
% G transponse 
Gt_AW2_real = G_AW2_real'; 
Gt_AW2_realimag = G_AW2_realimag'; 

% Takes real and imaginary parts: I0 and SO2
% Real
WTd_2p5e15_AW2_real = real(WTd_2p5e15_AW2); 
WTd_7p5e15_AW2_real = real(WTd_7p5e15_AW2); 
WTd_2p5e16_AW2_real = real(WTd_2p5e16_AW2); 
WTd_7p5e16_AW2_real = real(WTd_7p5e16_AW2); 
WTd_2p5e17_AW2_real = real(WTd_2p5e17_AW2); 
WTd_7p5e17_AW2_real = real(WTd_7p5e17_AW2); 
WTd_2p5e18_AW2_real = real(WTd_2p5e18_AW2); 
WTd_7p5e18_AW2_real = real(WTd_7p5e18_AW2); 
WTd_2p5e19_AW2_real = real(WTd_2p5e19_AW2); 
WTd_7p5e19_AW2_real = real(WTd_7p5e19_AW2); 
% Imaginary 
WTd_2p5e15_AW2_imag = imag(WTd_2p5e15_AW2); 
WTd_7p5e15_AW2_imag = imag(WTd_7p5e15_AW2); 
WTd_2p5e16_AW2_imag = imag(WTd_2p5e16_AW2); 
WTd_7p5e16_AW2_imag = imag(WTd_7p5e16_AW2); 
WTd_2p5e17_AW2_imag = imag(WTd_2p5e17_AW2); 
WTd_7p5e17_AW2_imag = imag(WTd_7p5e17_AW2); 
WTd_2p5e18_AW2_imag = imag(WTd_2p5e18_AW2); 
WTd_7p5e18_AW2_imag = imag(WTd_7p5e18_AW2); 
WTd_2p5e19_AW2_imag = imag(WTd_2p5e19_AW2); 
WTd_7p5e19_AW2_imag = imag(WTd_7p5e19_AW2); 
% Vectorises 
% Real
WTd_2p5e15_AW2_real_v = WTd_2p5e15_AW2_real(:); 
WTd_7p5e15_AW2_real_v = WTd_7p5e15_AW2_real(:); 
WTd_2p5e16_AW2_real_v = WTd_2p5e16_AW2_real(:); 
WTd_7p5e16_AW2_real_v = WTd_7p5e16_AW2_real(:); 
WTd_2p5e17_AW2_real_v = WTd_2p5e17_AW2_real(:); 
WTd_7p5e17_AW2_real_v = WTd_7p5e17_AW2_real(:); 
WTd_2p5e18_AW2_real_v = WTd_2p5e18_AW2_real(:); 
WTd_7p5e18_AW2_real_v = WTd_7p5e18_AW2_real(:); 
WTd_2p5e19_AW2_real_v = WTd_2p5e19_AW2_real(:); 
WTd_7p5e19_AW2_real_v = WTd_7p5e19_AW2_real(:); 
% Imaginary
WTd_2p5e15_AW2_imag_v = WTd_2p5e15_AW2_imag(:); 
WTd_7p5e15_AW2_imag_v = WTd_7p5e15_AW2_imag(:); 
WTd_2p5e16_AW2_imag_v = WTd_2p5e16_AW2_imag(:); 
WTd_7p5e16_AW2_imag_v = WTd_7p5e16_AW2_imag(:); 
WTd_2p5e17_AW2_imag_v = WTd_2p5e17_AW2_imag(:); 
WTd_7p5e17_AW2_imag_v = WTd_7p5e17_AW2_imag(:); 
WTd_2p5e18_AW2_imag_v = WTd_2p5e18_AW2_imag(:); 
WTd_7p5e18_AW2_imag_v = WTd_7p5e18_AW2_imag(:); 
WTd_2p5e19_AW2_imag_v = WTd_2p5e19_AW2_imag(:); 
WTd_7p5e19_AW2_imag_v = WTd_7p5e19_AW2_imag(:); 
% Combines real and complex magnitude 
WTd_2p5e15_AW2_realimag_v = [WTd_2p5e15_AW2_real_v; WTd_2p5e15_AW2_imag_v]; 
WTd_7p5e15_AW2_realimag_v = [WTd_7p5e15_AW2_real_v; WTd_7p5e15_AW2_imag_v]; 
WTd_2p5e16_AW2_realimag_v = [WTd_2p5e16_AW2_real_v; WTd_2p5e16_AW2_imag_v]; 
WTd_7p5e16_AW2_realimag_v = [WTd_7p5e16_AW2_real_v; WTd_7p5e16_AW2_imag_v]; 
WTd_2p5e17_AW2_realimag_v = [WTd_2p5e17_AW2_real_v; WTd_2p5e17_AW2_imag_v]; 
WTd_7p5e17_AW2_realimag_v = [WTd_7p5e17_AW2_real_v; WTd_7p5e17_AW2_imag_v]; 
WTd_2p5e18_AW2_realimag_v = [WTd_2p5e18_AW2_real_v; WTd_2p5e18_AW2_imag_v]; 
WTd_7p5e18_AW2_realimag_v = [WTd_7p5e18_AW2_real_v; WTd_7p5e18_AW2_imag_v]; 
WTd_2p5e19_AW2_realimag_v = [WTd_2p5e19_AW2_real_v; WTd_2p5e19_AW2_imag_v]; 
WTd_7p5e19_AW2_realimag_v = [WTd_7p5e19_AW2_real_v; WTd_7p5e19_AW2_imag_v]; 
% Remove NaN 
% Real only 
d_2p5e15_AW2_real_v = WTd_2p5e15_AW2_real_v(all(~isnan(WTd_2p5e15_AW2_real_v), 2),:);
d_7p5e15_AW2_real_v = WTd_7p5e15_AW2_real_v(all(~isnan(WTd_7p5e15_AW2_real_v), 2),:);
d_2p5e16_AW2_real_v = WTd_2p5e16_AW2_real_v(all(~isnan(WTd_2p5e16_AW2_real_v), 2),:);
d_7p5e16_AW2_real_v = WTd_7p5e16_AW2_real_v(all(~isnan(WTd_7p5e16_AW2_real_v), 2),:);
d_2p5e17_AW2_real_v = WTd_2p5e17_AW2_real_v(all(~isnan(WTd_2p5e17_AW2_real_v), 2),:);
d_7p5e17_AW2_real_v = WTd_7p5e17_AW2_real_v(all(~isnan(WTd_7p5e17_AW2_real_v), 2),:);
d_2p5e18_AW2_real_v = WTd_2p5e18_AW2_real_v(all(~isnan(WTd_2p5e18_AW2_real_v), 2),:);
d_7p5e18_AW2_real_v = WTd_7p5e18_AW2_real_v(all(~isnan(WTd_7p5e18_AW2_real_v), 2),:);
d_2p5e19_AW2_real_v = WTd_2p5e19_AW2_real_v(all(~isnan(WTd_2p5e19_AW2_real_v), 2),:);
d_7p5e19_AW2_real_v = WTd_7p5e19_AW2_real_v(all(~isnan(WTd_7p5e19_AW2_real_v), 2),:);
% Real and imag (complex magnitude) 
d_2p5e15_AW2_realimag_v = WTd_2p5e15_AW2_realimag_v(all(~isnan(WTd_2p5e15_AW2_realimag_v), 2),:);
d_7p5e15_AW2_realimag_v = WTd_7p5e15_AW2_realimag_v(all(~isnan(WTd_7p5e15_AW2_realimag_v), 2),:);
d_2p5e16_AW2_realimag_v = WTd_2p5e16_AW2_realimag_v(all(~isnan(WTd_2p5e16_AW2_realimag_v), 2),:);
d_7p5e16_AW2_realimag_v = WTd_7p5e16_AW2_realimag_v(all(~isnan(WTd_7p5e16_AW2_realimag_v), 2),:);
d_2p5e17_AW2_realimag_v = WTd_2p5e17_AW2_realimag_v(all(~isnan(WTd_2p5e17_AW2_realimag_v), 2),:);
d_7p5e17_AW2_realimag_v = WTd_7p5e17_AW2_realimag_v(all(~isnan(WTd_7p5e17_AW2_realimag_v), 2),:);
d_2p5e18_AW2_realimag_v = WTd_2p5e18_AW2_realimag_v(all(~isnan(WTd_2p5e18_AW2_realimag_v), 2),:);
d_7p5e18_AW2_realimag_v = WTd_7p5e18_AW2_realimag_v(all(~isnan(WTd_7p5e18_AW2_realimag_v), 2),:);
d_2p5e19_AW2_realimag_v = WTd_2p5e19_AW2_realimag_v(all(~isnan(WTd_2p5e19_AW2_realimag_v), 2),:);
d_7p5e19_AW2_realimag_v = WTd_7p5e19_AW2_realimag_v(all(~isnan(WTd_7p5e19_AW2_realimag_v), 2),:);
% Find mEst 
% A 
A_AW2_real = Gt_AW2_real * G_AW2_real;
A_AW2_realimag = Gt_AW2_realimag * G_AW2_realimag;
% B 
% Real only
B_AW2_real_d2p5e15 = Gt_AW2_real * d_2p5e15_AW2_real_v;
B_AW2_real_d7p5e15 = Gt_AW2_real * d_7p5e15_AW2_real_v;
B_AW2_real_d2p5e16 = Gt_AW2_real * d_2p5e16_AW2_real_v;
B_AW2_real_d7p5e16 = Gt_AW2_real * d_7p5e16_AW2_real_v;
B_AW2_real_d2p5e17 = Gt_AW2_real * d_2p5e17_AW2_real_v;
B_AW2_real_d7p5e17 = Gt_AW2_real * d_7p5e17_AW2_real_v;
B_AW2_real_d2p5e18 = Gt_AW2_real * d_2p5e18_AW2_real_v;
B_AW2_real_d7p5e18 = Gt_AW2_real * d_7p5e18_AW2_real_v;
B_AW2_real_d2p5e19 = Gt_AW2_real * d_2p5e19_AW2_real_v;
B_AW2_real_d7p5e19 = Gt_AW2_real * d_7p5e19_AW2_real_v;
% Real and imag (complex magnitude) 
B_AW2_realimag_d2p5e15 = Gt_AW2_realimag * d_2p5e15_AW2_realimag_v;
B_AW2_realimag_d7p5e15 = Gt_AW2_realimag * d_7p5e15_AW2_realimag_v;
B_AW2_realimag_d2p5e16 = Gt_AW2_realimag * d_2p5e16_AW2_realimag_v;
B_AW2_realimag_d7p5e16 = Gt_AW2_realimag * d_7p5e16_AW2_realimag_v;
B_AW2_realimag_d2p5e17 = Gt_AW2_realimag * d_2p5e17_AW2_realimag_v;
B_AW2_realimag_d7p5e17 = Gt_AW2_realimag * d_7p5e17_AW2_realimag_v;
B_AW2_realimag_d2p5e18 = Gt_AW2_realimag * d_2p5e18_AW2_realimag_v;
B_AW2_realimag_d7p5e18 = Gt_AW2_realimag * d_7p5e18_AW2_realimag_v;
B_AW2_realimag_d2p5e19 = Gt_AW2_realimag * d_2p5e19_AW2_realimag_v;
B_AW2_realimag_d7p5e19 = Gt_AW2_realimag * d_7p5e19_AW2_realimag_v;
% mEst 
% Real only 
mEst_AW2_real_2p5e15 = A_AW2_real\B_AW2_real_d2p5e15;
mEst_AW2_real_7p5e15 = A_AW2_real\B_AW2_real_d7p5e15;
mEst_AW2_real_2p5e16 = A_AW2_real\B_AW2_real_d2p5e16;
mEst_AW2_real_7p5e16 = A_AW2_real\B_AW2_real_d7p5e16;
mEst_AW2_real_2p5e17 = A_AW2_real\B_AW2_real_d2p5e17;
mEst_AW2_real_7p5e17 = A_AW2_real\B_AW2_real_d7p5e17;
mEst_AW2_real_2p5e18 = A_AW2_real\B_AW2_real_d2p5e18;
mEst_AW2_real_7p5e18 = A_AW2_real\B_AW2_real_d7p5e18;
mEst_AW2_real_2p5e19 = A_AW2_real\B_AW2_real_d2p5e19;
mEst_AW2_real_7p5e19 = A_AW2_real\B_AW2_real_d7p5e19;
% Real and imag (complex magnitude) 
mEst_AW2_realimag_2p5e15 = A_AW2_realimag\B_AW2_realimag_d2p5e15;
mEst_AW2_realimag_7p5e15 = A_AW2_realimag\B_AW2_realimag_d7p5e15;
mEst_AW2_realimag_2p5e16 = A_AW2_realimag\B_AW2_realimag_d2p5e16;
mEst_AW2_realimag_7p5e16 = A_AW2_realimag\B_AW2_realimag_d7p5e16;
mEst_AW2_realimag_2p5e17 = A_AW2_realimag\B_AW2_realimag_d2p5e17;
mEst_AW2_realimag_7p5e17 = A_AW2_realimag\B_AW2_realimag_d7p5e17;
mEst_AW2_realimag_2p5e18 = A_AW2_realimag\B_AW2_realimag_d2p5e18;
mEst_AW2_realimag_7p5e18 = A_AW2_realimag\B_AW2_realimag_d7p5e18;
mEst_AW2_realimag_2p5e19 = A_AW2_realimag\B_AW2_realimag_d2p5e19;
mEst_AW2_realimag_7p5e19 = A_AW2_realimag\B_AW2_realimag_d7p5e19;
% Results 
% concentration_SO2 = [1.5e15; 7.5e15; 2.5e16; 7.5e16; 2.5e17; 7.5e17; 2.5e18; 7.5e18; 2.5e19; 7.5e19]; 
mEst_AW2_real = [mEst_AW2_real_2p5e15, mEst_AW2_real_7p5e15, mEst_AW2_real_2p5e16, mEst_AW2_real_7p5e16, mEst_AW2_real_2p5e17, mEst_AW2_real_7p5e17, mEst_AW2_real_2p5e18, mEst_AW2_real_7p5e18, mEst_AW2_real_2p5e19, mEst_AW2_real_7p5e19];  
mEst_AW2_realimag = [mEst_AW2_realimag_2p5e15, mEst_AW2_realimag_7p5e15, mEst_AW2_realimag_2p5e16, mEst_AW2_realimag_7p5e16, mEst_AW2_realimag_2p5e17, mEst_AW2_realimag_7p5e17, mEst_AW2_realimag_2p5e18, mEst_AW2_realimag_7p5e18, mEst_AW2_realimag_2p5e19, mEst_AW2_realimag_7p5e19];  

% Finds mEst for AW3 
% Takes real and imaginary parts: I0 and SO2
% Real
WT_I0_AW3_real = real(WT_I0_AW3); 
WT_SO2n_AW3_real = real(WT_SO2n_AW3); 
% Imaginary 
WT_I0_AW3_imag = imag(WT_I0_AW3); 
WT_SO2n_AW3_imag = imag(WT_SO2n_AW3); 
% Vectorises 
% Real
WT_I0_AW3_real_v = WT_I0_AW3_real(:); 
WT_SO2n_AW3_real_v = WT_SO2n_AW3_real(:); 
% Imaginary
WT_I0_AW3_imag_v = WT_I0_AW3_imag(:); 
WT_SO2n_AW3_imag_v = WT_SO2n_AW3_imag(:); 
% Combines real and complex magnitude 
% I0
WT_I0_AW3_realimag = [WT_I0_AW3_real_v; WT_I0_AW3_imag_v]; 
% SO2
WT_SO2n_AW3_realimag = [WT_SO2n_AW3_real_v; WT_SO2n_AW3_imag_v]; 
% Creates design matrix, G 
G_AW3_real = [WT_I0_AW3_real_v, WT_SO2n_AW3_real_v]; 
G_AW3_realimag = [WT_I0_AW3_realimag, WT_SO2n_AW3_realimag]; 
% Remove NaN 
G_AW3_real = G_AW3_real(all(~isnan(G_AW3_real), 2),:);
G_AW3_realimag = G_AW3_realimag(all(~isnan(G_AW3_realimag), 2),:);
% G transponse 
Gt_AW3_real = G_AW3_real'; 
Gt_AW3_realimag = G_AW3_realimag'; 

% Takes real and imaginary parts: I0 and SO2
% Real
WTd_2p5e15_AW3_real = real(WTd_2p5e15_AW3); 
WTd_7p5e15_AW3_real = real(WTd_7p5e15_AW3); 
WTd_2p5e16_AW3_real = real(WTd_2p5e16_AW3); 
WTd_7p5e16_AW3_real = real(WTd_7p5e16_AW3); 
WTd_2p5e17_AW3_real = real(WTd_2p5e17_AW3); 
WTd_7p5e17_AW3_real = real(WTd_7p5e17_AW3); 
WTd_2p5e18_AW3_real = real(WTd_2p5e18_AW3); 
WTd_7p5e18_AW3_real = real(WTd_7p5e18_AW3); 
WTd_2p5e19_AW3_real = real(WTd_2p5e19_AW3); 
WTd_7p5e19_AW3_real = real(WTd_7p5e19_AW3); 
% Imaginary 
WTd_2p5e15_AW3_imag = imag(WTd_2p5e15_AW3); 
WTd_7p5e15_AW3_imag = imag(WTd_7p5e15_AW3); 
WTd_2p5e16_AW3_imag = imag(WTd_2p5e16_AW3); 
WTd_7p5e16_AW3_imag = imag(WTd_7p5e16_AW3); 
WTd_2p5e17_AW3_imag = imag(WTd_2p5e17_AW3); 
WTd_7p5e17_AW3_imag = imag(WTd_7p5e17_AW3); 
WTd_2p5e18_AW3_imag = imag(WTd_2p5e18_AW3); 
WTd_7p5e18_AW3_imag = imag(WTd_7p5e18_AW3); 
WTd_2p5e19_AW3_imag = imag(WTd_2p5e19_AW3); 
WTd_7p5e19_AW3_imag = imag(WTd_7p5e19_AW3); 
% Vectorises 
% Real
WTd_2p5e15_AW3_real_v = WTd_2p5e15_AW3_real(:); 
WTd_7p5e15_AW3_real_v = WTd_7p5e15_AW3_real(:); 
WTd_2p5e16_AW3_real_v = WTd_2p5e16_AW3_real(:); 
WTd_7p5e16_AW3_real_v = WTd_7p5e16_AW3_real(:); 
WTd_2p5e17_AW3_real_v = WTd_2p5e17_AW3_real(:); 
WTd_7p5e17_AW3_real_v = WTd_7p5e17_AW3_real(:); 
WTd_2p5e18_AW3_real_v = WTd_2p5e18_AW3_real(:); 
WTd_7p5e18_AW3_real_v = WTd_7p5e18_AW3_real(:); 
WTd_2p5e19_AW3_real_v = WTd_2p5e19_AW3_real(:); 
WTd_7p5e19_AW3_real_v = WTd_7p5e19_AW3_real(:); 
% Imaginary
WTd_2p5e15_AW3_imag_v = WTd_2p5e15_AW3_imag(:); 
WTd_7p5e15_AW3_imag_v = WTd_7p5e15_AW3_imag(:); 
WTd_2p5e16_AW3_imag_v = WTd_2p5e16_AW3_imag(:); 
WTd_7p5e16_AW3_imag_v = WTd_7p5e16_AW3_imag(:); 
WTd_2p5e17_AW3_imag_v = WTd_2p5e17_AW3_imag(:); 
WTd_7p5e17_AW3_imag_v = WTd_7p5e17_AW3_imag(:); 
WTd_2p5e18_AW3_imag_v = WTd_2p5e18_AW3_imag(:); 
WTd_7p5e18_AW3_imag_v = WTd_7p5e18_AW3_imag(:); 
WTd_2p5e19_AW3_imag_v = WTd_2p5e19_AW3_imag(:); 
WTd_7p5e19_AW3_imag_v = WTd_7p5e19_AW3_imag(:); 
% Combines real and complex magnitude 
WTd_2p5e15_AW3_realimag_v = [WTd_2p5e15_AW3_real_v; WTd_2p5e15_AW3_imag_v]; 
WTd_7p5e15_AW3_realimag_v = [WTd_7p5e15_AW3_real_v; WTd_7p5e15_AW3_imag_v]; 
WTd_2p5e16_AW3_realimag_v = [WTd_2p5e16_AW3_real_v; WTd_2p5e16_AW3_imag_v]; 
WTd_7p5e16_AW3_realimag_v = [WTd_7p5e16_AW3_real_v; WTd_7p5e16_AW3_imag_v]; 
WTd_2p5e17_AW3_realimag_v = [WTd_2p5e17_AW3_real_v; WTd_2p5e17_AW3_imag_v]; 
WTd_7p5e17_AW3_realimag_v = [WTd_7p5e17_AW3_real_v; WTd_7p5e17_AW3_imag_v]; 
WTd_2p5e18_AW3_realimag_v = [WTd_2p5e18_AW3_real_v; WTd_2p5e18_AW3_imag_v]; 
WTd_7p5e18_AW3_realimag_v = [WTd_7p5e18_AW3_real_v; WTd_7p5e18_AW3_imag_v]; 
WTd_2p5e19_AW3_realimag_v = [WTd_2p5e19_AW3_real_v; WTd_2p5e19_AW3_imag_v]; 
WTd_7p5e19_AW3_realimag_v = [WTd_7p5e19_AW3_real_v; WTd_7p5e19_AW3_imag_v]; 
% Remove NaN 
% Real only 
d_2p5e15_AW3_real_v = WTd_2p5e15_AW3_real_v(all(~isnan(WTd_2p5e15_AW3_real_v), 2),:);
d_7p5e15_AW3_real_v = WTd_7p5e15_AW3_real_v(all(~isnan(WTd_7p5e15_AW3_real_v), 2),:);
d_2p5e16_AW3_real_v = WTd_2p5e16_AW3_real_v(all(~isnan(WTd_2p5e16_AW3_real_v), 2),:);
d_7p5e16_AW3_real_v = WTd_7p5e16_AW3_real_v(all(~isnan(WTd_7p5e16_AW3_real_v), 2),:);
d_2p5e17_AW3_real_v = WTd_2p5e17_AW3_real_v(all(~isnan(WTd_2p5e17_AW3_real_v), 2),:);
d_7p5e17_AW3_real_v = WTd_7p5e17_AW3_real_v(all(~isnan(WTd_7p5e17_AW3_real_v), 2),:);
d_2p5e18_AW3_real_v = WTd_2p5e18_AW3_real_v(all(~isnan(WTd_2p5e18_AW3_real_v), 2),:);
d_7p5e18_AW3_real_v = WTd_7p5e18_AW3_real_v(all(~isnan(WTd_7p5e18_AW3_real_v), 2),:);
d_2p5e19_AW3_real_v = WTd_2p5e19_AW3_real_v(all(~isnan(WTd_2p5e19_AW3_real_v), 2),:);
d_7p5e19_AW3_real_v = WTd_7p5e19_AW3_real_v(all(~isnan(WTd_7p5e19_AW3_real_v), 2),:);
% Real and imag (complex magnitude) 
d_2p5e15_AW3_realimag_v = WTd_2p5e15_AW3_realimag_v(all(~isnan(WTd_2p5e15_AW3_realimag_v), 2),:);
d_7p5e15_AW3_realimag_v = WTd_7p5e15_AW3_realimag_v(all(~isnan(WTd_7p5e15_AW3_realimag_v), 2),:);
d_2p5e16_AW3_realimag_v = WTd_2p5e16_AW3_realimag_v(all(~isnan(WTd_2p5e16_AW3_realimag_v), 2),:);
d_7p5e16_AW3_realimag_v = WTd_7p5e16_AW3_realimag_v(all(~isnan(WTd_7p5e16_AW3_realimag_v), 2),:);
d_2p5e17_AW3_realimag_v = WTd_2p5e17_AW3_realimag_v(all(~isnan(WTd_2p5e17_AW3_realimag_v), 2),:);
d_7p5e17_AW3_realimag_v = WTd_7p5e17_AW3_realimag_v(all(~isnan(WTd_7p5e17_AW3_realimag_v), 2),:);
d_2p5e18_AW3_realimag_v = WTd_2p5e18_AW3_realimag_v(all(~isnan(WTd_2p5e18_AW3_realimag_v), 2),:);
d_7p5e18_AW3_realimag_v = WTd_7p5e18_AW3_realimag_v(all(~isnan(WTd_7p5e18_AW3_realimag_v), 2),:);
d_2p5e19_AW3_realimag_v = WTd_2p5e19_AW3_realimag_v(all(~isnan(WTd_2p5e19_AW3_realimag_v), 2),:);
d_7p5e19_AW3_realimag_v = WTd_7p5e19_AW3_realimag_v(all(~isnan(WTd_7p5e19_AW3_realimag_v), 2),:);
% Find mEst 
% A 
A_AW3_real = Gt_AW3_real * G_AW3_real;
A_AW3_realimag = Gt_AW3_realimag * G_AW3_realimag;
% B 
% Real only
B_AW3_real_d2p5e15 = Gt_AW3_real * d_2p5e15_AW3_real_v;
B_AW3_real_d7p5e15 = Gt_AW3_real * d_7p5e15_AW3_real_v;
B_AW3_real_d2p5e16 = Gt_AW3_real * d_2p5e16_AW3_real_v;
B_AW3_real_d7p5e16 = Gt_AW3_real * d_7p5e16_AW3_real_v;
B_AW3_real_d2p5e17 = Gt_AW3_real * d_2p5e17_AW3_real_v;
B_AW3_real_d7p5e17 = Gt_AW3_real * d_7p5e17_AW3_real_v;
B_AW3_real_d2p5e18 = Gt_AW3_real * d_2p5e18_AW3_real_v;
B_AW3_real_d7p5e18 = Gt_AW3_real * d_7p5e18_AW3_real_v;
B_AW3_real_d2p5e19 = Gt_AW3_real * d_2p5e19_AW3_real_v;
B_AW3_real_d7p5e19 = Gt_AW3_real * d_7p5e19_AW3_real_v;
% Real and imag (complex magnitude) 
B_AW3_realimag_d2p5e15 = Gt_AW3_realimag * d_2p5e15_AW3_realimag_v;
B_AW3_realimag_d7p5e15 = Gt_AW3_realimag * d_7p5e15_AW3_realimag_v;
B_AW3_realimag_d2p5e16 = Gt_AW3_realimag * d_2p5e16_AW3_realimag_v;
B_AW3_realimag_d7p5e16 = Gt_AW3_realimag * d_7p5e16_AW3_realimag_v;
B_AW3_realimag_d2p5e17 = Gt_AW3_realimag * d_2p5e17_AW3_realimag_v;
B_AW3_realimag_d7p5e17 = Gt_AW3_realimag * d_7p5e17_AW3_realimag_v;
B_AW3_realimag_d2p5e18 = Gt_AW3_realimag * d_2p5e18_AW3_realimag_v;
B_AW3_realimag_d7p5e18 = Gt_AW3_realimag * d_7p5e18_AW3_realimag_v;
B_AW3_realimag_d2p5e19 = Gt_AW3_realimag * d_2p5e19_AW3_realimag_v;
B_AW3_realimag_d7p5e19 = Gt_AW3_realimag * d_7p5e19_AW3_realimag_v;
% mEst 
% Real only 
mEst_AW3_real_2p5e15 = A_AW3_real\B_AW3_real_d2p5e15;
mEst_AW3_real_7p5e15 = A_AW3_real\B_AW3_real_d7p5e15;
mEst_AW3_real_2p5e16 = A_AW3_real\B_AW3_real_d2p5e16;
mEst_AW3_real_7p5e16 = A_AW3_real\B_AW3_real_d7p5e16;
mEst_AW3_real_2p5e17 = A_AW3_real\B_AW3_real_d2p5e17;
mEst_AW3_real_7p5e17 = A_AW3_real\B_AW3_real_d7p5e17;
mEst_AW3_real_2p5e18 = A_AW3_real\B_AW3_real_d2p5e18;
mEst_AW3_real_7p5e18 = A_AW3_real\B_AW3_real_d7p5e18;
mEst_AW3_real_2p5e19 = A_AW3_real\B_AW3_real_d2p5e19;
mEst_AW3_real_7p5e19 = A_AW3_real\B_AW3_real_d7p5e19;
% Real and imag (complex magnitude) 
mEst_AW3_realimag_2p5e15 = A_AW3_realimag\B_AW3_realimag_d2p5e15;
mEst_AW3_realimag_7p5e15 = A_AW3_realimag\B_AW3_realimag_d7p5e15;
mEst_AW3_realimag_2p5e16 = A_AW3_realimag\B_AW3_realimag_d2p5e16;
mEst_AW3_realimag_7p5e16 = A_AW3_realimag\B_AW3_realimag_d7p5e16;
mEst_AW3_realimag_2p5e17 = A_AW3_realimag\B_AW3_realimag_d2p5e17;
mEst_AW3_realimag_7p5e17 = A_AW3_realimag\B_AW3_realimag_d7p5e17;
mEst_AW3_realimag_2p5e18 = A_AW3_realimag\B_AW3_realimag_d2p5e18;
mEst_AW3_realimag_7p5e18 = A_AW3_realimag\B_AW3_realimag_d7p5e18;
mEst_AW3_realimag_2p5e19 = A_AW3_realimag\B_AW3_realimag_d2p5e19;
mEst_AW3_realimag_7p5e19 = A_AW3_realimag\B_AW3_realimag_d7p5e19;
% Results 
% concentration_SO2 = [1.5e15; 7.5e15; 2.5e16; 7.5e16; 2.5e17; 7.5e17; 2.5e18; 7.5e18; 2.5e19; 7.5e19]; 
mEst_AW3_real = [mEst_AW3_real_2p5e15, mEst_AW3_real_7p5e15, mEst_AW3_real_2p5e16, mEst_AW3_real_7p5e16, mEst_AW3_real_2p5e17, mEst_AW3_real_7p5e17, mEst_AW3_real_2p5e18, mEst_AW3_real_7p5e18, mEst_AW3_real_2p5e19, mEst_AW3_real_7p5e19];  
mEst_AW3_realimag = [mEst_AW3_realimag_2p5e15, mEst_AW3_realimag_7p5e15, mEst_AW3_realimag_2p5e16, mEst_AW3_realimag_7p5e16, mEst_AW3_realimag_2p5e17, mEst_AW3_realimag_7p5e17, mEst_AW3_realimag_2p5e18, mEst_AW3_realimag_7p5e18, mEst_AW3_realimag_2p5e19, mEst_AW3_realimag_7p5e19];  

% Plot 
figure % Figure 9  
subplot(2,1,1)
% AW1
p = plot(concentration_SO2, mEst_AW1_real(2,:), '^'); hold on; p.Color = AW1col; p.MarkerSize = 6; 
p = plot(concentration_SO2, mEst_AW1_realimag(2,:), 'v'); p.Color = AW1col; p.MarkerSize = 6; 
% AW2
p = plot(concentration_SO2, mEst_AW2_real(2,:), '^'); hold on; p.Color = AW2col; p.MarkerSize = 6; 
p = plot(concentration_SO2, mEst_AW2_realimag(2,:), 'v'); p.Color = AW2col; p.MarkerSize = 6; 
% AW3
p = plot(concentration_SO2, mEst_AW3_real(2,:), '^'); hold on; p.Color = AW3col; p.MarkerSize = 6; 
p = plot(concentration_SO2, mEst_AW3_realimag(2,:), 'v'); p.Color = AW3col; p.MarkerSize = 6; 
%
plot(concentration_SO2, concentration_SO2, ':k');
ylabel('mEst_S_O_2 (molec/cm^2)')
xlabel('Concentration (molec/cm^2)')
set(gcf,'color','w');
fig = gca; fig.FontSize = 14; 
xlim([-0.1e19 8e19])
ylim([-0.5e19 8e19])
subplot(2,1,2)
% AW1
p = plot(concentration_SO2, mEst_AW1_real(2,:), '^'); hold on; p.Color = AW1col; p.MarkerSize = 6; 
p = plot(concentration_SO2, mEst_AW1_realimag(2,:), 'v'); p.Color = AW1col; p.MarkerSize = 6; 
% AW2
p = plot(concentration_SO2, mEst_AW2_real(2,:), '^'); hold on; p.Color = AW2col; p.MarkerSize = 6;  
p = plot(concentration_SO2, mEst_AW2_realimag(2,:), 'v'); p.Color = AW2col; p.MarkerSize = 6;  
% AW3
p = plot(concentration_SO2, mEst_AW3_real(2,:), '^'); hold on; p.Color = AW3col; p.MarkerSize = 6; 
p = plot(concentration_SO2, mEst_AW3_realimag(2,:), 'v'); p.Color = AW3col; p.MarkerSize = 6;  
%
plot(concentration_SO2, concentration_SO2, ':k');
ylabel('mEst_S_O_2 (molec/cm^2)')
xlabel('Concentration (molec/cm^2)')
set(gcf,'color','w');
fig = gca; fig.FontSize = 14; 
xlim([-0.3e18 8e18])
ylim([-1.5e18 8e18])
% Save convoluted cross sections to be used later 
fname = fullfile(outDir, 'mEst_real_imag_analysis_window');
save(fname); % Saves AW1 variables  
%% Make any adjustments 
% SAVES FIGURE
% fnOut = 'Fig9_mEst_real_imag_AW';
% fname = fullfile(outDirFig, fnOut);
% saveas(gca, fname)
























