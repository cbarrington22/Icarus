 
% Synthetic spectra - all wavelength ranges and all frequencies

% Must run SyntheticSpectra_AllWavelengths.m first to find WT of d and G
% This script extracts the data from that plotted in Fig. 6 
% Must run SyntheticSpectra_analysisWindow_croppedF.m to obtain indicies for spatial frequency cropping  
 
outDir = '/Users/charlotteb/Documents/Chapter 2/Synthetic spectra/All wavelengths cropped frequency Fig 12 and Fig 13/'; 
outDirFig = '/Users/charlotteb/Documents/Chapter 2/Synthetic spectra/All wavelengths cropped frequency Fig 12 and Fig 13/'; 

% Crops spatial frequency  
% Uses... WTd_2p5e15_all etc. from SyntheticSpectra_AllWavelengths.m
% Uses iyMin and iyMax from SyntheticSpectra_analysisWindow_allF.m
% Crops spatial frequency of WT of d 
WTd_2p5e15_all_F = WTd_2p5e15_all(iyMax:iyMin, :); 
WTd_7p5e15_all_F = WTd_7p5e15_all(iyMax:iyMin, :); 
WTd_2p5e16_all_F = WTd_2p5e16_all(iyMax:iyMin, :); 
WTd_7p5e16_all_F = WTd_7p5e16_all(iyMax:iyMin, :); 
WTd_2p5e17_all_F = WTd_2p5e17_all(iyMax:iyMin, :); 
WTd_7p5e17_all_F = WTd_7p5e17_all(iyMax:iyMin, :); 
WTd_2p5e18_all_F = WTd_2p5e18_all(iyMax:iyMin, :); 
WTd_7p5e18_all_F = WTd_7p5e18_all(iyMax:iyMin, :); 
WTd_2p5e19_all_F = WTd_2p5e19_all(iyMax:iyMin, :); 
WTd_7p5e19_all_F = WTd_7p5e19_all(iyMax:iyMin, :);  
% G 
WT_I0_all_F = WT_I0_all(iyMax:iyMin, :);  
WT_SO2n_all_F = WT_SO2n_all(iyMax:iyMin, :);   
% Figure 12 - plot linear model for all wavelengths now with spatial frequency cropped 
figure 
% Subplots 
% First column of figure is d (plot for 2.5e18) 
subplot(3,1,1) 
pcolor(lambda, Fc, real(WTd_2p5e18_all_F)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
xlim([min(lambda) max(lambda)]); set(gca,'YScale','log'); ylabel('Spatial frequency (cycles/\lambda)', 'FontSize', 10);
% cb = colorbar; cb.Label.Position(1) = 3.5;
% ylabel(cb,'Amplitude', 'Rotation', -90, 'FontSize', 10);
set(gcf,'color','w');
title('d')
subtitle('J(\lambda) prime')
subplot(3,1,2) 
pcolor(lambda, Fc, real(WT_I0_all_F)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
xlim([min(lambda) max(lambda)]); set(gca,'YScale','log'); % ylabel('Spatial frequency (cycles/\lambda)', 'FontSize', 10);
% cb = colorbar; cb.Label.Position(1) = 3.5;
% ylabel(cb,'Amplitude', 'Rotation', -90, 'FontSize', 10);
xlabel('\lambda (nm)', 'FontSize', 14); 
set(gcf,'color','w');
title('G') 
subtitle('J_0(\lambda) prime')
subplot(3,1,3) 
pcolor(lambda, Fc, real(WT_SO2n_all_F)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
xlim([min(lambda) max(lambda)]); set(gca,'YScale','log'); % ylabel('Spatial frequency (cycles/\lambda)', 'FontSize', 10);
cb = colorbar; cb.Label.Position(1) = 3.5;
ylabel(cb,'Amplitude', 'Rotation', -90, 'FontSize', 10);
subtitle('\sigma_S_O_2 prime')
set(gcf,'color','w');
% Save 
fname = fullfile(outDir, 'mEst_real_imag_ALL_F');
save(fname); % Saves all variables  
%% Make any adjustments 
% SAVES FIGURE
% fnOut = 'Fig12_dG_ALL_F';
% fname = fullfile(outDirFig, fnOut);
% saveas(gca, fname)
% Takes real and imaginary parts: I0 and SO2
% Real
WT_I0_all_F_real = real(WT_I0_all_F); 
WT_SO2n_all_F_real = real(WT_SO2n_all_F); 
% Imaginary 
WT_I0_all_F_imag = imag(WT_I0_all_F); 
WT_SO2n_all_F_imag = imag(WT_SO2n_all_F); 
% Vectorises 
% Real
WT_I0_all_F_real_v = WT_I0_all_F_real(:); 
WT_SO2n_all_F_real_v = WT_SO2n_all_F_real(:); 
% Imaginary
WT_I0_all_F_imag_v = WT_I0_all_F_imag(:); 
WT_SO2n_all_F_imag_v = WT_SO2n_all_F_imag(:); 
% Combines real and complex magnitude 
% I0
WT_I0_all_F_realimag = [WT_I0_all_F_real_v; WT_I0_all_F_imag_v]; 
% SO2
WT_SO2n_all_F_realimag = [WT_SO2n_all_F_real_v; WT_SO2n_all_F_imag_v]; 
% Creates design matrix, G 
G_all_F_real = [WT_I0_all_F_real_v, WT_SO2n_all_F_real_v]; 
G_all_F_realimag = [WT_I0_all_F_realimag, WT_SO2n_all_F_realimag]; 
% Remove NaN 
G_all_F_real = G_all_F_real(all(~isnan(G_all_F_real), 2),:);
G_all_F_realimag = G_all_F_realimag(all(~isnan(G_all_F_realimag), 2),:);
% G transponse 
Gt_all_F_real = G_all_F_real'; 
Gt_all_F_realimag = G_all_F_realimag'; 
% d 
% Takes real and imaginary parts: I0 and SO2
% Real
WTd_2p5e15_all_F_real = real(WTd_2p5e15_all_F); 
WTd_7p5e15_all_F_real = real(WTd_7p5e15_all_F); 
WTd_2p5e16_all_F_real = real(WTd_2p5e16_all_F); 
WTd_7p5e16_all_F_real = real(WTd_7p5e16_all_F); 
WTd_2p5e17_all_F_real = real(WTd_2p5e17_all_F); 
WTd_7p5e17_all_F_real = real(WTd_7p5e17_all_F); 
WTd_2p5e18_all_F_real = real(WTd_2p5e18_all_F); 
WTd_7p5e18_all_F_real = real(WTd_7p5e18_all_F); 
WTd_2p5e19_all_F_real = real(WTd_2p5e19_all_F); 
WTd_7p5e19_all_F_real = real(WTd_7p5e19_all_F); 
% Imaginary 
WTd_2p5e15_all_F_imag = imag(WTd_2p5e15_all_F); 
WTd_7p5e15_all_F_imag = imag(WTd_7p5e15_all_F); 
WTd_2p5e16_all_F_imag = imag(WTd_2p5e16_all_F); 
WTd_7p5e16_all_F_imag = imag(WTd_7p5e16_all_F); 
WTd_2p5e17_all_F_imag = imag(WTd_2p5e17_all_F); 
WTd_7p5e17_all_F_imag = imag(WTd_7p5e17_all_F); 
WTd_2p5e18_all_F_imag = imag(WTd_2p5e18_all_F); 
WTd_7p5e18_all_F_imag = imag(WTd_7p5e18_all_F); 
WTd_2p5e19_all_F_imag = imag(WTd_2p5e19_all_F); 
WTd_7p5e19_all_F_imag = imag(WTd_7p5e19_all_F); 
% Vectorises 
% Real
WTd_2p5e15_all_F_real_v = WTd_2p5e15_all_F_real(:); 
WTd_7p5e15_all_F_real_v = WTd_7p5e15_all_F_real(:); 
WTd_2p5e16_all_F_real_v = WTd_2p5e16_all_F_real(:); 
WTd_7p5e16_all_F_real_v = WTd_7p5e16_all_F_real(:); 
WTd_2p5e17_all_F_real_v = WTd_2p5e17_all_F_real(:); 
WTd_7p5e17_all_F_real_v = WTd_7p5e17_all_F_real(:); 
WTd_2p5e18_all_F_real_v = WTd_2p5e18_all_F_real(:); 
WTd_7p5e18_all_F_real_v = WTd_7p5e18_all_F_real(:); 
WTd_2p5e19_all_F_real_v = WTd_2p5e19_all_F_real(:); 
WTd_7p5e19_all_F_real_v = WTd_7p5e19_all_F_real(:); 
% Imaginary
WTd_2p5e15_all_F_imag_v = WTd_2p5e15_all_F_imag(:); 
WTd_7p5e15_all_F_imag_v = WTd_7p5e15_all_F_imag(:); 
WTd_2p5e16_all_F_imag_v = WTd_2p5e16_all_F_imag(:); 
WTd_7p5e16_all_F_imag_v = WTd_7p5e16_all_F_imag(:); 
WTd_2p5e17_all_F_imag_v = WTd_2p5e17_all_F_imag(:); 
WTd_7p5e17_all_F_imag_v = WTd_7p5e17_all_F_imag(:); 
WTd_2p5e18_all_F_imag_v = WTd_2p5e18_all_F_imag(:); 
WTd_7p5e18_all_F_imag_v = WTd_7p5e18_all_F_imag(:); 
WTd_2p5e19_all_F_imag_v = WTd_2p5e19_all_F_imag(:); 
WTd_7p5e19_all_F_imag_v = WTd_7p5e19_all_F_imag(:); 
% Combines real and complex magnitude 
WTd_2p5e15_all_F_realimag_v = [WTd_2p5e15_all_F_real_v; WTd_2p5e15_all_F_imag_v]; 
WTd_7p5e15_all_F_realimag_v = [WTd_7p5e15_all_F_real_v; WTd_7p5e15_all_F_imag_v]; 
WTd_2p5e16_all_F_realimag_v = [WTd_2p5e16_all_F_real_v; WTd_2p5e16_all_F_imag_v]; 
WTd_7p5e16_all_F_realimag_v = [WTd_7p5e16_all_F_real_v; WTd_7p5e16_all_F_imag_v]; 
WTd_2p5e17_all_F_realimag_v = [WTd_2p5e17_all_F_real_v; WTd_2p5e17_all_F_imag_v]; 
WTd_7p5e17_all_F_realimag_v = [WTd_7p5e17_all_F_real_v; WTd_7p5e17_all_F_imag_v]; 
WTd_2p5e18_all_F_realimag_v = [WTd_2p5e18_all_F_real_v; WTd_2p5e18_all_F_imag_v]; 
WTd_7p5e18_all_F_realimag_v = [WTd_7p5e18_all_F_real_v; WTd_7p5e18_all_F_imag_v]; 
WTd_2p5e19_all_F_realimag_v = [WTd_2p5e19_all_F_real_v; WTd_2p5e19_all_F_imag_v]; 
WTd_7p5e19_all_F_realimag_v = [WTd_7p5e19_all_F_real_v; WTd_7p5e19_all_F_imag_v]; 
% Remove NaN 
% Real only 
d_2p5e15_all_F_real_v = WTd_2p5e15_all_F_real_v(all(~isnan(WTd_2p5e15_all_F_real_v), 2),:);
d_7p5e15_all_F_real_v = WTd_7p5e15_all_F_real_v(all(~isnan(WTd_7p5e15_all_F_real_v), 2),:);
d_2p5e16_all_F_real_v = WTd_2p5e16_all_F_real_v(all(~isnan(WTd_2p5e16_all_F_real_v), 2),:);
d_7p5e16_all_F_real_v = WTd_7p5e16_all_F_real_v(all(~isnan(WTd_7p5e16_all_F_real_v), 2),:);
d_2p5e17_all_F_real_v = WTd_2p5e17_all_F_real_v(all(~isnan(WTd_2p5e17_all_F_real_v), 2),:);
d_7p5e17_all_F_real_v = WTd_7p5e17_all_F_real_v(all(~isnan(WTd_7p5e17_all_F_real_v), 2),:);
d_2p5e18_all_F_real_v = WTd_2p5e18_all_F_real_v(all(~isnan(WTd_2p5e18_all_F_real_v), 2),:);
d_7p5e18_all_F_real_v = WTd_7p5e18_all_F_real_v(all(~isnan(WTd_7p5e18_all_F_real_v), 2),:);
d_2p5e19_all_F_real_v = WTd_2p5e19_all_F_real_v(all(~isnan(WTd_2p5e19_all_F_real_v), 2),:);
d_7p5e19_all_F_real_v = WTd_7p5e19_all_F_real_v(all(~isnan(WTd_7p5e19_all_F_real_v), 2),:);
% Real and imag (complex magnitude) 
d_2p5e15_all_F_realimag_v = WTd_2p5e15_all_F_realimag_v(all(~isnan(WTd_2p5e15_all_F_realimag_v), 2),:);
d_7p5e15_all_F_realimag_v = WTd_7p5e15_all_F_realimag_v(all(~isnan(WTd_7p5e15_all_F_realimag_v), 2),:);
d_2p5e16_all_F_realimag_v = WTd_2p5e16_all_F_realimag_v(all(~isnan(WTd_2p5e16_all_F_realimag_v), 2),:);
d_7p5e16_all_F_realimag_v = WTd_7p5e16_all_F_realimag_v(all(~isnan(WTd_7p5e16_all_F_realimag_v), 2),:);
d_2p5e17_all_F_realimag_v = WTd_2p5e17_all_F_realimag_v(all(~isnan(WTd_2p5e17_all_F_realimag_v), 2),:);
d_7p5e17_all_F_realimag_v = WTd_7p5e17_all_F_realimag_v(all(~isnan(WTd_7p5e17_all_F_realimag_v), 2),:);
d_2p5e18_all_F_realimag_v = WTd_2p5e18_all_F_realimag_v(all(~isnan(WTd_2p5e18_all_F_realimag_v), 2),:);
d_7p5e18_all_F_realimag_v = WTd_7p5e18_all_F_realimag_v(all(~isnan(WTd_7p5e18_all_F_realimag_v), 2),:);
d_2p5e19_all_F_realimag_v = WTd_2p5e19_all_F_realimag_v(all(~isnan(WTd_2p5e19_all_F_realimag_v), 2),:);
d_7p5e19_all_F_realimag_v = WTd_7p5e19_all_F_realimag_v(all(~isnan(WTd_7p5e19_all_F_realimag_v), 2),:);
% Find mEst 
% A 
A_all_F_real = Gt_all_F_real * G_all_F_real;
A_all_F_realimag = Gt_all_F_realimag * G_all_F_realimag;
% B 
% Real only
B_all_F_real_d2p5e15 = Gt_all_F_real * d_2p5e15_all_F_real_v;
B_all_F_real_d7p5e15 = Gt_all_F_real * d_7p5e15_all_F_real_v;
B_all_F_real_d2p5e16 = Gt_all_F_real * d_2p5e16_all_F_real_v;
B_all_F_real_d7p5e16 = Gt_all_F_real * d_7p5e16_all_F_real_v;
B_all_F_real_d2p5e17 = Gt_all_F_real * d_2p5e17_all_F_real_v;
B_all_F_real_d7p5e17 = Gt_all_F_real * d_7p5e17_all_F_real_v;
B_all_F_real_d2p5e18 = Gt_all_F_real * d_2p5e18_all_F_real_v;
B_all_F_real_d7p5e18 = Gt_all_F_real * d_7p5e18_all_F_real_v;
B_all_F_real_d2p5e19 = Gt_all_F_real * d_2p5e19_all_F_real_v;
B_all_F_real_d7p5e19 = Gt_all_F_real * d_7p5e19_all_F_real_v;
% Real and imag (complex magnitude) 
B_all_F_realimag_d2p5e15 = Gt_all_F_realimag * d_2p5e15_all_F_realimag_v;
B_all_F_realimag_d7p5e15 = Gt_all_F_realimag * d_7p5e15_all_F_realimag_v;
B_all_F_realimag_d2p5e16 = Gt_all_F_realimag * d_2p5e16_all_F_realimag_v;
B_all_F_realimag_d7p5e16 = Gt_all_F_realimag * d_7p5e16_all_F_realimag_v;
B_all_F_realimag_d2p5e17 = Gt_all_F_realimag * d_2p5e17_all_F_realimag_v;
B_all_F_realimag_d7p5e17 = Gt_all_F_realimag * d_7p5e17_all_F_realimag_v;
B_all_F_realimag_d2p5e18 = Gt_all_F_realimag * d_2p5e18_all_F_realimag_v;
B_all_F_realimag_d7p5e18 = Gt_all_F_realimag * d_7p5e18_all_F_realimag_v;
B_all_F_realimag_d2p5e19 = Gt_all_F_realimag * d_2p5e19_all_F_realimag_v;
B_all_F_realimag_d7p5e19 = Gt_all_F_realimag * d_7p5e19_all_F_realimag_v;
% mEst 
% Real only 
mEst_all_F_real_2p5e15 = A_all_F_real\B_all_F_real_d2p5e15;
mEst_all_F_real_7p5e15 = A_all_F_real\B_all_F_real_d7p5e15;
mEst_all_F_real_2p5e16 = A_all_F_real\B_all_F_real_d2p5e16;
mEst_all_F_real_7p5e16 = A_all_F_real\B_all_F_real_d7p5e16;
mEst_all_F_real_2p5e17 = A_all_F_real\B_all_F_real_d2p5e17;
mEst_all_F_real_7p5e17 = A_all_F_real\B_all_F_real_d7p5e17;
mEst_all_F_real_2p5e18 = A_all_F_real\B_all_F_real_d2p5e18;
mEst_all_F_real_7p5e18 = A_all_F_real\B_all_F_real_d7p5e18;
mEst_all_F_real_2p5e19 = A_all_F_real\B_all_F_real_d2p5e19;
mEst_all_F_real_7p5e19 = A_all_F_real\B_all_F_real_d7p5e19;
% Real and imag (complex magnitude) 
mEst_all_F_realimag_2p5e15 = A_all_F_realimag\B_all_F_realimag_d2p5e15;
mEst_all_F_realimag_7p5e15 = A_all_F_realimag\B_all_F_realimag_d7p5e15;
mEst_all_F_realimag_2p5e16 = A_all_F_realimag\B_all_F_realimag_d2p5e16;
mEst_all_F_realimag_7p5e16 = A_all_F_realimag\B_all_F_realimag_d7p5e16;
mEst_all_F_realimag_2p5e17 = A_all_F_realimag\B_all_F_realimag_d2p5e17;
mEst_all_F_realimag_7p5e17 = A_all_F_realimag\B_all_F_realimag_d7p5e17;
mEst_all_F_realimag_2p5e18 = A_all_F_realimag\B_all_F_realimag_d2p5e18;
mEst_all_F_realimag_7p5e18 = A_all_F_realimag\B_all_F_realimag_d7p5e18;
mEst_all_F_realimag_2p5e19 = A_all_F_realimag\B_all_F_realimag_d2p5e19;
mEst_all_F_realimag_7p5e19 = A_all_F_realimag\B_all_F_realimag_d7p5e19;
% Results 
concentration_SO2 = [2.5e15; 7.5e15; 2.5e16; 7.5e16; 2.5e17; 7.5e17; 2.5e18; 7.5e18; 2.5e19; 7.5e19]; 
mEst_all_F_real = [mEst_all_F_real_2p5e15, mEst_all_F_real_7p5e15, mEst_all_F_real_2p5e16, mEst_all_F_real_7p5e16, mEst_all_F_real_2p5e17, mEst_all_F_real_7p5e17, mEst_all_F_real_2p5e18, mEst_all_F_real_7p5e18, mEst_all_F_real_2p5e19, mEst_all_F_real_7p5e19];  
mEst_all_F_realimag = [mEst_all_F_realimag_2p5e15, mEst_all_F_realimag_7p5e15, mEst_all_F_realimag_2p5e16, mEst_all_F_realimag_7p5e16, mEst_all_F_realimag_2p5e17, mEst_all_F_realimag_7p5e17, mEst_all_F_realimag_2p5e18, mEst_all_F_realimag_7p5e18, mEst_all_F_realimag_2p5e19, mEst_all_F_realimag_7p5e19];  
% Plot 
figure % Figure 13 
subplot(2,1,1)
p = plot(concentration_SO2, mEst_all_F_real(2,:), '^k'); hold on; p.MarkerSize = 8; p.MarkerFaceColor = 'k';  
p = plot(concentration_SO2, mEst_all_F_realimag(2,:), 'vk'); p.MarkerSize = 8; p.MarkerFaceColor = 'k'; 
plot(concentration_SO2, concentration_SO2, ':k');
ylabel('mEst_S_O_2 (molec/cm^2)')
xlabel('SO_2 column density (molec/cm^2)')
set(gcf,'color','w');
fig = gca; fig.FontSize = 14; 
xlim([-0.1e19 8e19])
ylim([-0.1e19 8e19])
subplot(2,1,2)
p = plot(concentration_SO2, mEst_all_F_real(2,:), '^k'); hold on; p.MarkerSize = 8; p.MarkerFaceColor = 'k'; 
p = plot(concentration_SO2, mEst_all_F_realimag(2,:), 'vk'); p.MarkerSize = 8; p.MarkerFaceColor = 'k'; 
plot(concentration_SO2, concentration_SO2, ':k');
ylabel('mEst_S_O_2 (molec/cm^2)')
xlabel('SO_2 column density (molec/cm^2)')
set(gcf,'color','w');
fig = gca; fig.FontSize = 14; 
xlim([-0.3e18 8e18])
ylim([-0.3e18 8e18])
% Save convoluted cross sections to be used later 
fname = fullfile(outDir, 'mEst_real_imag_ALL_F');
save(fname); % Saves all variables  
%% Make any adjustments 
% SAVES FIGURE
% fnOut = 'Fig13_mEst_real_imag_ALL_F';
% fname = fullfile(outDirFig, fnOut);
% saveas(gca, fname)

























