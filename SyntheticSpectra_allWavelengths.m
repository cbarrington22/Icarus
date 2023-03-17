% Synthetic spectra 

% Directory to synthetic spectra 
inDir = '/Users/charlotteb/Documents/Chapter 2/Synthetic spectra/inFiles/'; 
outDir = '/Users/charlotteb/Documents/Chapter 2/Synthetic spectra/All wavelengths Fig 5 Fig 6 and Fig 7/'; 
outDirFig = '/Users/charlotteb/Documents/Chapter 2/Synthetic spectra/All wavelengths Fig 5 Fig 6 and Fig 7/'; 
% Loads synthetic spectra for d
% 2.5e15
inI = fullfile(inDir, 'SO2_2.5e15_conv_0.5nm.txt'); 
[wI, iI, x] = textread(inI,'%s%s%s'); 
wI = str2double(wI);
iI = str2double(iI);            
d_2p5e15 = [wI, iI]; 
% 7.5e15
inI = fullfile(inDir, 'SO2_7.5e15_conv_0.5nm.txt'); 
[wI, iI] = textread(inI,'%s%s'); 
wI = str2double(wI);
iI = str2double(iI);            
d_7p5e15 = [wI, iI]; 
% 2.5e16
inI = fullfile(inDir, 'SO2_2.5e16_conv_0.5nm.txt'); 
[wI, iI] = textread(inI,'%s%s'); 
wI = str2double(wI);
iI = str2double(iI);            
d_2p5e16 = [wI, iI]; 
% 7.5e16
inI = fullfile(inDir, 'SO2_7.5e16_conv_0.5nm.txt'); 
[wI, iI] = textread(inI,'%s%s'); 
wI = str2double(wI);
iI = str2double(iI);            
d_7p5e16 = [wI, iI]; 
% 2.5e17
inI = fullfile(inDir, 'SO2_2.5e17_conv_0.5nm.txt'); 
[wI, iI] = textread(inI,'%s%s'); 
wI = str2double(wI);
iI = str2double(iI);            
d_2p5e17 = [wI, iI]; 
% 7.5e17
inI = fullfile(inDir, 'SO2_7.5e17_conv_0.5nm.txt'); 
[wI, iI] = textread(inI,'%s%s'); 
wI = str2double(wI);
iI = str2double(iI);            
d_7p5e17 = [wI, iI]; 
% 2.5e18
inI = fullfile(inDir, 'SO2_2.5e18_conv_0.5nm.txt'); 
[wI, iI] = textread(inI,'%s%s'); 
wI = str2double(wI);
iI = str2double(iI);            
d_2p5e18 = [wI, iI]; 
% 7.5e18
inI = fullfile(inDir, 'SO2_7.5e18_conv_0.5nm.txt'); 
[wI, iI] = textread(inI,'%s%s'); 
wI = str2double(wI);
iI = str2double(iI);            
d_7p5e18 = [wI, iI]; 
% 2.5e19
inI = fullfile(inDir, 'SO2_2.5e19_conv_0.5nm.txt'); 
[wI, iI] = textread(inI,'%s%s'); 
wI = str2double(wI);
iI = str2double(iI);            
d_2p5e19 = [wI, iI]; 
% 7.5e19
inI = fullfile(inDir, 'SO2_7.5e19_conv_0.5e19.txt'); 
[wI, iI] = textread(inI,'%s%s'); 
wI = str2double(wI);
iI = str2double(iI);            
d_7p5e19 = [wI, iI]; 
% Figure 5 - Synthetic spectra (a) magnitude (b) frequency 
% Plots d in magnitude domain
figure
% Spectra colour 
col = [0 0 0]; 
% Analysis windows 
% AW1 - x
AW1col = [0.11 0.24 0.81]; 
AW1_l1 = 310; 
AW1_l2 = 340; 
% AW2 - red 
AW2col = [0.91 0.05 0.05]; 
AW2_l1 = 320; 
AW2_l2 = 340; 
% AW3 - green 
AW3col = [0.09 0.54 0.03]; 
AW3_l1 = 314; 
AW3_l2 = 330; 
% (a) magnitude 
% Subplots 
subplot(5,2,1)
s = plot(d_2p5e15(:,1), d_2p5e15(:,2)); hold on; s.Color = col; 
xlim([280 max(wI)])
set(gcf,'color','w');
ylabel('Intensity')
xlabel('\lambda (nm)')
subtitle('2.5 · 10^1^5')
x = xline(AW1_l1, '-'); x.Color = AW1col; x.LineWidth = 0.7; 
x = xline(AW1_l2, '-'); x.Color = AW1col; x.LineWidth = 0.7; 
x = xline(AW2_l1, ':'); x.Color = AW2col; x.LineWidth = 1; 
x = xline(AW2_l2, ':'); x.Color = AW2col; x.LineWidth = 1; 
x = xline(AW3_l1, '--'); x.Color = AW3col; x.LineWidth = 0.7; 
x = xline(AW3_l2, '--'); x.Color = AW3col; x.LineWidth = 0.7; 
subplot(5,2,2)
s = plot(d_7p5e15(:,1), d_7p5e15(:,2)); hold on; s.Color = col; 
xlim([280 max(wI)])
set(gcf,'color','w');
ylabel('Intensity')
xlabel('\lambda (nm)')
subtitle('7.5 · 10^1^5')
x = xline(AW1_l1, '-'); x.Color = AW1col; x.LineWidth = 0.7; 
x = xline(AW1_l2, '-'); x.Color = AW1col; x.LineWidth = 0.7; 
x = xline(AW2_l1, ':'); x.Color = AW2col; x.LineWidth = 1; 
x = xline(AW2_l2, ':'); x.Color = AW2col; x.LineWidth = 1; 
x = xline(AW3_l1, '--'); x.Color = AW3col; x.LineWidth = 0.7; 
x = xline(AW3_l2, '--'); x.Color = AW3col; x.LineWidth = 0.7; 
subplot(5,2,3)
s = plot(d_2p5e16(:,1), d_2p5e16(:,2)); hold on; s.Color = col; 
xlim([280 max(wI)])
set(gcf,'color','w');
ylabel('Intensity')
xlabel('\lambda (nm)')
subtitle('2.5 · 10^1^6')
x = xline(AW1_l1, '-'); x.Color = AW1col; x.LineWidth = 0.7; 
x = xline(AW1_l2, '-'); x.Color = AW1col; x.LineWidth = 0.7; 
x = xline(AW2_l1, ':'); x.Color = AW2col; x.LineWidth = 1; 
x = xline(AW2_l2, ':'); x.Color = AW2col; x.LineWidth = 1; 
x = xline(AW3_l1, '--'); x.Color = AW3col; x.LineWidth = 0.7; 
x = xline(AW3_l2, '--'); x.Color = AW3col; x.LineWidth = 0.7; 
subplot(5,2,4)
s = plot(d_7p5e16(:,1), d_7p5e16(:,2)); hold on; s.Color = col; 
xlim([280 max(wI)])
set(gcf,'color','w');
ylabel('Intensity')
xlabel('\lambda (nm)')
subtitle('7.5 · 10^1^6')
x = xline(AW1_l1, '-'); x.Color = AW1col; x.LineWidth = 0.7; 
x = xline(AW1_l2, '-'); x.Color = AW1col; x.LineWidth = 0.7; 
x = xline(AW2_l1, ':'); x.Color = AW2col; x.LineWidth = 1; 
x = xline(AW2_l2, ':'); x.Color = AW2col; x.LineWidth = 1; 
x = xline(AW3_l1, '--'); x.Color = AW3col; x.LineWidth = 0.7; 
x = xline(AW3_l2, '--'); x.Color = AW3col; x.LineWidth = 0.7; 
subplot(5,2,5)
s = plot(d_2p5e17(:,1), d_2p5e17(:,2)); hold on; s.Color = col; 
xlim([280 max(wI)])
set(gcf,'color','w');
ylabel('Intensity')
xlabel('\lambda (nm)')
subtitle('2.5 · 10^1^7')
x = xline(AW1_l1, '-'); x.Color = AW1col; x.LineWidth = 0.7; 
x = xline(AW1_l2, '-'); x.Color = AW1col; x.LineWidth = 0.7; 
x = xline(AW2_l1, ':'); x.Color = AW2col; x.LineWidth = 1; 
x = xline(AW2_l2, ':'); x.Color = AW2col; x.LineWidth = 1; 
x = xline(AW3_l1, '--'); x.Color = AW3col; x.LineWidth = 0.7; 
x = xline(AW3_l2, '--'); x.Color = AW3col; x.LineWidth = 0.7; 
subplot(5,2,6)
s = plot(d_7p5e17(:,1), d_7p5e17(:,2)); hold on; s.Color = col; 
xlim([280 max(wI)])
set(gcf,'color','w');
ylabel('Intensity')
xlabel('\lambda (nm)')
subtitle('7.5 · 10^1^7')
x = xline(AW1_l1, '-'); x.Color = AW1col; x.LineWidth = 0.7; 
x = xline(AW1_l2, '-'); x.Color = AW1col; x.LineWidth = 0.7; 
x = xline(AW2_l1, ':'); x.Color = AW2col; x.LineWidth = 1; 
x = xline(AW2_l2, ':'); x.Color = AW2col; x.LineWidth = 1; 
x = xline(AW3_l1, '--'); x.Color = AW3col; x.LineWidth = 0.7; 
x = xline(AW3_l2, '--'); x.Color = AW3col; x.LineWidth = 0.7; 
subplot(5,2,7)
s = plot(d_2p5e18(:,1), d_7p5e18(:,2)); hold on; s.Color = col; 
xlim([280 max(wI)])
set(gcf,'color','w');
ylabel('Intensity')
xlabel('\lambda (nm)')
subtitle('2.5 · 10^1^8')
x = xline(AW1_l1, '-'); x.Color = AW1col; x.LineWidth = 0.7; 
x = xline(AW1_l2, '-'); x.Color = AW1col; x.LineWidth = 0.7; 
x = xline(AW2_l1, ':'); x.Color = AW2col; x.LineWidth = 1; 
x = xline(AW2_l2, ':'); x.Color = AW2col; x.LineWidth = 1; 
x = xline(AW3_l1, '--'); x.Color = AW3col; x.LineWidth = 0.7; 
x = xline(AW3_l2, '--'); x.Color = AW3col; x.LineWidth = 0.7; 
subplot(5,2,8)
s = plot(d_7p5e18(:,1), d_7p5e18(:,2)); hold on; s.Color = col; 
xlim([280 max(wI)])
set(gcf,'color','w');
ylabel('Intensity')
xlabel('\lambda (nm)')
subtitle('7.5 · 10^1^8')
x = xline(AW1_l1, '-'); x.Color = AW1col; x.LineWidth = 0.7; 
x = xline(AW1_l2, '-'); x.Color = AW1col; x.LineWidth = 0.7; 
x = xline(AW2_l1, ':'); x.Color = AW2col; x.LineWidth = 1; 
x = xline(AW2_l2, ':'); x.Color = AW2col; x.LineWidth = 1; 
x = xline(AW3_l1, '--'); x.Color = AW3col; x.LineWidth = 0.7; 
x = xline(AW3_l2, '--'); x.Color = AW3col; x.LineWidth = 0.7; 
subplot(5,2,9)
s = plot(d_2p5e19(:,1), d_2p5e19(:,2)); hold on; s.Color = col; 
xlim([280 max(wI)])
set(gcf,'color','w');
ylabel('Intensity')
xlabel('\lambda (nm)')
subtitle('2.5 · 10^1^9')
x = xline(AW1_l1, '-'); x.Color = AW1col; x.LineWidth = 0.7; 
x = xline(AW1_l2, '-'); x.Color = AW1col; x.LineWidth = 0.7; 
x = xline(AW2_l1, ':'); x.Color = AW2col; x.LineWidth = 1; 
x = xline(AW2_l2, ':'); x.Color = AW2col; x.LineWidth = 1; 
x = xline(AW3_l1, '--'); x.Color = AW3col; x.LineWidth = 0.7; 
x = xline(AW3_l2, '--'); x.Color = AW3col; x.LineWidth = 0.7; 
subplot(5,2,10)
s = plot(d_7p5e19(:,1), d_7p5e19(:,2)); hold on; s.Color = col; 
xlim([280 max(wI)])
set(gcf,'color','w');
ylabel('Intensity')
xlabel('\lambda (nm)')
subtitle('7.5 · 10^1^9')
x = xline(AW1_l1, '-'); x.Color = AW1col; x.LineWidth = 0.7; 
x = xline(AW1_l2, '-'); x.Color = AW1col; x.LineWidth = 0.7; 
x = xline(AW2_l1, ':'); x.Color = AW2col; x.LineWidth = 1; 
x = xline(AW2_l2, ':'); x.Color = AW2col; x.LineWidth = 1; 
x = xline(AW3_l1, '--'); x.Color = AW3col; x.LineWidth = 0.7; 
x = xline(AW3_l2, '--'); x.Color = AW3col; x.LineWidth = 0.7; 
% Save convoluted cross sections to be used later 
fname = fullfile(outDir, 'SyntheticSpectra_magnitude');
save(fname); % Saves all variables  
%% Make any adjustments 
% SAVES FIGURE
% fnOut = 'Fig5_SyntheticMagnitude';
% fname = fullfile(outDirFig, fnOut);
% saveas(gca, fname)
% (b) frequency  
lambda = d_2p5e19(:,1); 
chnls = length(lambda); % Number of spectrometer channels  
minlambda = min(lambda); % Starting wavelength 
maxlambda = max(lambda); % Final wavelength 
rangeW = maxlambda - minlambda; % Determines wavelength range of spectrometer 
Fs = rangeW/chnls; % Fs is the sampling rate in wavelength, each data point (or channel) is equivalent to Fs wavelengths 
% Extracts data betwwen 280 and 420 nm (max)
% Index (channel) corresponding to 280 and 420 
ind280 = 16572; 
d_2p5e15 = d_2p5e15(ind280:end,:);
d_7p5e15 = d_7p5e15(ind280:end,:);
d_2p5e16 = d_2p5e16(ind280:end,:);
d_7p5e16 = d_7p5e16(ind280:end,:);
d_2p5e17 = d_2p5e17(ind280:end,:);
d_7p5e17 = d_7p5e17(ind280:end,:);
d_2p5e18 = d_2p5e18(ind280:end,:);
d_7p5e18 = d_7p5e18(ind280:end,:);
d_2p5e19 = d_2p5e19(ind280:end,:);
d_7p5e19 = d_7p5e19(ind280:end,:);
lambda = d_7p5e19(:,1);
% ln(d) 
d_ln_2p5e15 = log(d_2p5e15(:,2)); 
d_ln_7p5e15 = log(d_7p5e15(:,2)); 
d_ln_2p5e16 = log(d_2p5e16(:,2)); 
d_ln_7p5e16 = log(d_7p5e16(:,2)); 
d_ln_2p5e17 = log(d_2p5e17(:,2)); 
d_ln_7p5e17 = log(d_7p5e17(:,2)); 
d_ln_2p5e18 = log(d_2p5e18(:,2));
d_ln_7p5e18 = log(d_7p5e18(:,2));
d_ln_2p5e19 = log(d_2p5e19(:,2));  
d_ln_7p5e19 = log(d_7p5e19(:,2));  
% CWT 
[WTd_2p5e15, ~, ~] = cwt(d_ln_2p5e15, Fs); 
[WTd_7p5e15, ~, ~] = cwt(d_ln_7p5e15, Fs); 
[WTd_2p5e16, ~, ~] = cwt(d_ln_2p5e16, Fs);  
[WTd_7p5e16, ~, ~] = cwt(d_ln_7p5e16, Fs); 
[WTd_2p5e17,  ~, ~] = cwt(d_ln_2p5e17, Fs); 
[WTd_7p5e17,  ~, ~] = cwt(d_ln_7p5e17, Fs);  
[WTd_2p5e18,  ~, ~] = cwt(d_ln_2p5e18, Fs); 
[WTd_7p5e18,  ~, ~] = cwt(d_ln_7p5e18, Fs); 
[WTd_2p5e19, ~, ~] = cwt(d_ln_2p5e19, Fs); 
[WTd_7p5e19, F, COI] = cwt(d_ln_7p5e19, Fs); 
% Removes COI
szCoi = length(COI); szF = length(F); matF = repelem(F, 1, szCoi); colCoi = COI'; matCoi = repelem(colCoi, szF, 1); idxCoi = matF <= matCoi; % Removes COI
WTd_2p5e15(idxCoi) = NaN; 
WTd_7p5e15(idxCoi) = NaN; 
WTd_2p5e16(idxCoi) = NaN; 
WTd_7p5e16(idxCoi) = NaN; 
WTd_2p5e17(idxCoi) = NaN; 
WTd_7p5e17(idxCoi) = NaN; 
WTd_2p5e18(idxCoi) = NaN; 
WTd_7p5e18(idxCoi) = NaN; 
WTd_2p5e19(idxCoi) = NaN; 
WTd_7p5e19(idxCoi) = NaN; 
figure
% (b) frequency 
% Subplots 
subplot(5,2,1)
pcolor(lambda, F, real(WTd_2p5e15)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
xlim([min(lambda) max(lambda)]); ylabel('Spatial frequency (cycles/\lambda)', 'FontSize', 10);set(gca,'YScale','log'); 
cb = colorbar; cb.Label.Position(1) = 3.5;
ylabel(cb,'Amplitude', 'Rotation', -90, 'FontSize', 10);
set(gcf,'color','w');
subtitle('2.5 · 10^1^5')
subplot(5,2,2)
pcolor(lambda, F, real(WTd_7p5e15)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
xlim([min(lambda) max(lambda)]); ylabel('Spatial frequency (cycles/\lambda)', 'FontSize', 10);set(gca,'YScale','log'); 
cb = colorbar; cb.Label.Position(1) = 3.5;
ylabel(cb,'Amplitude', 'Rotation', -90, 'FontSize', 10);
set(gcf,'color','w');
subtitle('7.5 · 10^1^5')
subplot(5,2,3)
pcolor(lambda, F, real(WTd_2p5e16)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
xlim([min(lambda) max(lambda)]); ylabel('Spatial frequency (cycles/\lambda)', 'FontSize', 10);set(gca,'YScale','log'); 
cb = colorbar; cb.Label.Position(1) = 3.5;
ylabel(cb,'Amplitude', 'Rotation', -90, 'FontSize', 10);
set(gcf,'color','w');
subtitle('2.5 · 10^1^6')
subplot(5,2,4)
pcolor(lambda, F, real(WTd_7p5e16)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
xlim([min(lambda) max(lambda)]); ylabel('Spatial frequency (cycles/\lambda)', 'FontSize', 10);set(gca,'YScale','log'); 
cb = colorbar; cb.Label.Position(1) = 3.5;
ylabel(cb,'Amplitude', 'Rotation', -90, 'FontSize', 10);
set(gcf,'color','w');
subtitle('7.5 · 10^1^6')
subplot(5,2,5)
pcolor(lambda, F, real(WTd_2p5e17)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
xlim([min(lambda) max(lambda)]); ylabel('Spatial frequency (cycles/\lambda)', 'FontSize', 10);set(gca,'YScale','log'); 
cb = colorbar; cb.Label.Position(1) = 3.5;
ylabel(cb,'Amplitude', 'Rotation', -90, 'FontSize', 10);
set(gcf,'color','w');
subtitle('2.5 · 10^1^7')
subplot(5,2,6)
pcolor(lambda, F, real(WTd_7p5e17)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
xlim([min(lambda) max(lambda)]); ylabel('Spatial frequency (cycles/\lambda)', 'FontSize', 10);set(gca,'YScale','log'); 
cb = colorbar; cb.Label.Position(1) = 3.5;
ylabel(cb,'Amplitude', 'Rotation', -90, 'FontSize', 10);
set(gcf,'color','w');
subtitle('7.5 · 10^1^7')
subplot(5,2,7)
pcolor(lambda, F, real(WTd_2p5e18)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
xlim([min(lambda) max(lambda)]); ylabel('Spatial frequency (cycles/\lambda)', 'FontSize', 10);set(gca,'YScale','log'); 
cb = colorbar; cb.Label.Position(1) = 3.5;
ylabel(cb,'Amplitude', 'Rotation', -90, 'FontSize', 10);
set(gcf,'color','w');
subtitle('2.5 · 10^1^8')
subplot(5,2,8)
pcolor(lambda, F, real(WTd_7p5e18)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
xlim([min(lambda) max(lambda)]); ylabel('Spatial frequency (cycles/\lambda)', 'FontSize', 10);set(gca,'YScale','log'); 
cb = colorbar; cb.Label.Position(1) = 3.5;
ylabel(cb,'Amplitude', 'Rotation', -90, 'FontSize', 10);
set(gcf,'color','w');
subtitle('7.5 · 10^1^8')
subplot(5,2,9)
pcolor(lambda, F, real(WTd_2p5e19)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
xlim([min(lambda) max(lambda)]); ylabel('Spatial frequency (cycles/\lambda)', 'FontSize', 10);set(gca,'YScale','log'); 
cb = colorbar; cb.Label.Position(1) = 3.5;
ylabel(cb,'Amplitude', 'Rotation', -90, 'FontSize', 10);
set(gcf,'color','w');
subtitle('2.5 · 10^1^9')
subplot(5,2,10)
pcolor(lambda, F, real(WTd_7p5e19)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
xlim([min(lambda) max(lambda)]); ylabel('Spatial frequency (cycles/\lambda)', 'FontSize', 10);set(gca,'YScale','log'); 
cb = colorbar; cb.Label.Position(1) = 3.5;
ylabel(cb,'Amplitude', 'Rotation', -90, 'FontSize', 10);
set(gcf,'color','w');
subtitle('7.5 · 10^1^9')
% Save convoluted cross sections to be used later 
fname = fullfile(outDir, 'SyntheticSpectra_frequency');
save(fname); % Saves all variables  
%% Make any adjustments 
% SAVES FIGURE
% fnOut = 'Fig5_SyntheticFrequency';
% fname = fullfile(outDirFig, fnOut);
% saveas(gca, fname)

% Figure 6 - deisgn matrix G 
% G 
% loads I0 
inI = fullfile(inDir, 'solar_ChanceKurucz_2010_unk_conv.txt'); 
[wI, iI] = textread(inI,'%s%s'); 
wI = str2double(wI);
iI = str2double(iI);            
I0 = [wI, iI]; 
% Extracts data betwwen 280 and 420 nm (max)
I0 = I0(ind280:end,2);
% 
I0 = log(I0);
% CWT 
[WT_I0, ~, ~] = cwt(I0, Fs);  
% Removes COI 
WT_I0(idxCoi) = NaN; 
figure 
subplot(2,1,1) % a 
pcolor(lambda, F, real(WT_I0)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
xlim([min(lambda) max(lambda)]); ylabel('Spatial frequency (cycles/\lambda)', 'FontSize', 14);set(gca,'YScale','log'); 
xlabel('\lambda (nm)', 'FontSize', 14); 
cb = colorbar; cb.Label.Position(1) = 3.5;
ylabel(cb,'Amplitude', 'Rotation', -90, 'FontSize', 14);
set(gcf,'color','w'); fig = gca; fig.FontSize = 14; 
subtitle('J_0(\lambda) prime')
% loads cross section SO2 
inI = fullfile(inDir, 'SO2_Vandaele_2009_298K_unk_conv.txt'); 
[wI, iI] = textread(inI,'%s%s'); 
wI = str2double(wI);
iI = str2double(iI);            
SO2 = [wI, iI]; 
% Extracts data betwwen 280 and 420 nm (max)
SO2 = SO2(ind280:end,2);
% Negative cross section
SO2n = -(SO2); % Negative 
% CWT 
[WT_SO2n, ~, ~] = cwt(SO2n); 
% Removes COI
WT_SO2n(idxCoi) = NaN; 
subplot(2,1,2) % b 
pcolor(lambda, F, real(WT_SO2n)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
xlim([min(lambda) max(lambda)]); ylabel('Spatial frequency (cycles/\lambda)', 'FontSize', 14);set(gca,'YScale','log'); 
xlabel('\lambda (nm)', 'FontSize', 14); 
cb = colorbar; cb.Label.Position(1) = 3.5;
ylabel(cb,'Amplitude', 'Rotation', -90, 'FontSize', 14);
set(gcf,'color','w'); fig = gca; fig.FontSize = 14; 
subtitle('\sigma_S_O_2 prime')
% Save convoluted cross sections to be used later 
fname = fullfile(outDir, 'G_frequency');
save(fname); % Saves all variables  
%% Make any adjustments 
% SAVES FIGURE
% fnOut = 'Fig6_G';
% fname = fullfile(outDirFig, fnOut);
% saveas(gca, fname)

% Figure 7 - mEstSO2 
% Prepare linear model 
% Makes copy of WT_G
WT_I0_all = WT_I0;
WT_SO2n_all = WT_SO2n; 
% Takes real and imaginary parts: I0 and SO2
% Real
WT_I0_all_real = real(WT_I0_all); 
WT_SO2n_all_real = real(WT_SO2n_all); 
% Imaginary 
WT_I0_all_imag = imag(WT_I0_all); 
WT_SO2n_all_imag = imag(WT_SO2n_all); 
% Vectorises 
% Real
WT_I0_all_real_v = WT_I0_all_real(:); 
WT_SO2n_all_real_v = WT_SO2n_all_real(:); 
% Imaginary
WT_I0_all_imag_v = WT_I0_all_imag(:); 
WT_SO2n_all_imag_v = WT_SO2n_all_imag(:); 
% Combines real and complex magnitude 
% I0
WT_I0_all_realimag = [WT_I0_all_real_v; WT_I0_all_imag_v]; 
% SO2
WT_SO2n_all_realimag = [WT_SO2n_all_real_v; WT_SO2n_all_imag_v]; 
% Creates design matrix, G 
G_all_real = [WT_I0_all_real_v, WT_SO2n_all_real_v]; 
G_all_realimag = [WT_I0_all_realimag, WT_SO2n_all_realimag]; 
% Remove NaN 
G_all_real = G_all_real(all(~isnan(G_all_real), 2),:);
G_all_realimag = G_all_realimag(all(~isnan(G_all_realimag), 2),:);
% G transponse 
Gt_all_real = G_all_real'; 
Gt_all_realimag = G_all_realimag'; 
% d 
% Makes copy of WT_d
WTd_2p5e15_all = WTd_2p5e15; 
WTd_7p5e15_all = WTd_7p5e15; 
WTd_2p5e16_all = WTd_2p5e16; 
WTd_7p5e16_all = WTd_7p5e16; 
WTd_2p5e17_all = WTd_2p5e17; 
WTd_7p5e17_all = WTd_7p5e17; 
WTd_2p5e18_all = WTd_2p5e18; 
WTd_7p5e18_all = WTd_7p5e18; 
WTd_2p5e19_all = WTd_2p5e19; 
WTd_7p5e19_all = WTd_7p5e19; 
% Takes real and imaginary parts: I0 and SO2
% Real
WTd_2p5e15_all_real = real(WTd_2p5e15_all); 
WTd_7p5e15_all_real = real(WTd_7p5e15_all); 
WTd_2p5e16_all_real = real(WTd_2p5e16_all); 
WTd_7p5e16_all_real = real(WTd_7p5e16_all); 
WTd_2p5e17_all_real = real(WTd_2p5e17_all); 
WTd_7p5e17_all_real = real(WTd_7p5e17_all); 
WTd_2p5e18_all_real = real(WTd_2p5e18_all); 
WTd_7p5e18_all_real = real(WTd_7p5e18_all); 
WTd_2p5e19_all_real = real(WTd_2p5e19_all); 
WTd_7p5e19_all_real = real(WTd_7p5e19_all); 
% Imaginary 
WTd_2p5e15_all_imag = imag(WTd_2p5e15_all); 
WTd_7p5e15_all_imag = imag(WTd_7p5e15_all); 
WTd_2p5e16_all_imag = imag(WTd_2p5e16_all); 
WTd_7p5e16_all_imag = imag(WTd_7p5e16_all); 
WTd_2p5e17_all_imag = imag(WTd_2p5e17_all); 
WTd_7p5e17_all_imag = imag(WTd_7p5e17_all); 
WTd_2p5e18_all_imag = imag(WTd_2p5e18_all); 
WTd_7p5e18_all_imag = imag(WTd_7p5e18_all); 
WTd_2p5e19_all_imag = imag(WTd_2p5e19_all); 
WTd_7p5e19_all_imag = imag(WTd_7p5e19_all); 
% Vectorises 
% Real
WTd_2p5e15_all_real_v = WTd_2p5e15_all_real(:); 
WTd_7p5e15_all_real_v = WTd_7p5e15_all_real(:); 
WTd_2p5e16_all_real_v = WTd_2p5e16_all_real(:); 
WTd_7p5e16_all_real_v = WTd_7p5e16_all_real(:); 
WTd_2p5e17_all_real_v = WTd_2p5e17_all_real(:); 
WTd_7p5e17_all_real_v = WTd_7p5e17_all_real(:); 
WTd_2p5e18_all_real_v = WTd_2p5e18_all_real(:); 
WTd_7p5e18_all_real_v = WTd_7p5e18_all_real(:); 
WTd_2p5e19_all_real_v = WTd_2p5e19_all_real(:); 
WTd_7p5e19_all_real_v = WTd_7p5e19_all_real(:); 
% Imaginary
WTd_2p5e15_all_imag_v = WTd_2p5e15_all_imag(:); 
WTd_7p5e15_all_imag_v = WTd_7p5e15_all_imag(:); 
WTd_2p5e16_all_imag_v = WTd_2p5e16_all_imag(:); 
WTd_7p5e16_all_imag_v = WTd_7p5e16_all_imag(:); 
WTd_2p5e17_all_imag_v = WTd_2p5e17_all_imag(:); 
WTd_7p5e17_all_imag_v = WTd_7p5e17_all_imag(:); 
WTd_2p5e18_all_imag_v = WTd_2p5e18_all_imag(:); 
WTd_7p5e18_all_imag_v = WTd_7p5e18_all_imag(:); 
WTd_2p5e19_all_imag_v = WTd_2p5e19_all_imag(:); 
WTd_7p5e19_all_imag_v = WTd_7p5e19_all_imag(:); 
% Combines real and complex magnitude 
WTd_2p5e15_all_realimag_v = [WTd_2p5e15_all_real_v; WTd_2p5e15_all_imag_v]; 
WTd_7p5e15_all_realimag_v = [WTd_7p5e15_all_real_v; WTd_7p5e15_all_imag_v]; 
WTd_2p5e16_all_realimag_v = [WTd_2p5e16_all_real_v; WTd_2p5e16_all_imag_v]; 
WTd_7p5e16_all_realimag_v = [WTd_7p5e16_all_real_v; WTd_7p5e16_all_imag_v]; 
WTd_2p5e17_all_realimag_v = [WTd_2p5e17_all_real_v; WTd_2p5e17_all_imag_v]; 
WTd_7p5e17_all_realimag_v = [WTd_7p5e17_all_real_v; WTd_7p5e17_all_imag_v]; 
WTd_2p5e18_all_realimag_v = [WTd_2p5e18_all_real_v; WTd_2p5e18_all_imag_v]; 
WTd_7p5e18_all_realimag_v = [WTd_7p5e18_all_real_v; WTd_7p5e18_all_imag_v]; 
WTd_2p5e19_all_realimag_v = [WTd_2p5e19_all_real_v; WTd_2p5e19_all_imag_v]; 
WTd_7p5e19_all_realimag_v = [WTd_7p5e19_all_real_v; WTd_7p5e19_all_imag_v]; 
% Remove NaN 
% Real only 
d_2p5e15_all_real_v = WTd_2p5e15_all_real_v(all(~isnan(WTd_2p5e15_all_real_v), 2),:);
d_7p5e15_all_real_v = WTd_7p5e15_all_real_v(all(~isnan(WTd_7p5e15_all_real_v), 2),:);
d_2p5e16_all_real_v = WTd_2p5e16_all_real_v(all(~isnan(WTd_2p5e16_all_real_v), 2),:);
d_7p5e16_all_real_v = WTd_7p5e16_all_real_v(all(~isnan(WTd_7p5e16_all_real_v), 2),:);
d_2p5e17_all_real_v = WTd_2p5e17_all_real_v(all(~isnan(WTd_2p5e17_all_real_v), 2),:);
d_7p5e17_all_real_v = WTd_7p5e17_all_real_v(all(~isnan(WTd_7p5e17_all_real_v), 2),:);
d_2p5e18_all_real_v = WTd_2p5e18_all_real_v(all(~isnan(WTd_2p5e18_all_real_v), 2),:);
d_7p5e18_all_real_v = WTd_7p5e18_all_real_v(all(~isnan(WTd_7p5e18_all_real_v), 2),:);
d_2p5e19_all_real_v = WTd_2p5e19_all_real_v(all(~isnan(WTd_2p5e19_all_real_v), 2),:);
d_7p5e19_all_real_v = WTd_7p5e19_all_real_v(all(~isnan(WTd_7p5e19_all_real_v), 2),:);
% Real and imag (complex magnitude) 
d_2p5e15_all_realimag_v = WTd_2p5e15_all_realimag_v(all(~isnan(WTd_2p5e15_all_realimag_v), 2),:);
d_7p5e15_all_realimag_v = WTd_7p5e15_all_realimag_v(all(~isnan(WTd_7p5e15_all_realimag_v), 2),:);
d_2p5e16_all_realimag_v = WTd_2p5e16_all_realimag_v(all(~isnan(WTd_2p5e16_all_realimag_v), 2),:);
d_7p5e16_all_realimag_v = WTd_7p5e16_all_realimag_v(all(~isnan(WTd_7p5e16_all_realimag_v), 2),:);
d_2p5e17_all_realimag_v = WTd_2p5e17_all_realimag_v(all(~isnan(WTd_2p5e17_all_realimag_v), 2),:);
d_7p5e17_all_realimag_v = WTd_7p5e17_all_realimag_v(all(~isnan(WTd_7p5e17_all_realimag_v), 2),:);
d_2p5e18_all_realimag_v = WTd_2p5e18_all_realimag_v(all(~isnan(WTd_2p5e18_all_realimag_v), 2),:);
d_7p5e18_all_realimag_v = WTd_7p5e18_all_realimag_v(all(~isnan(WTd_7p5e18_all_realimag_v), 2),:);
d_2p5e19_all_realimag_v = WTd_2p5e19_all_realimag_v(all(~isnan(WTd_2p5e19_all_realimag_v), 2),:);
d_7p5e19_all_realimag_v = WTd_7p5e19_all_realimag_v(all(~isnan(WTd_7p5e19_all_realimag_v), 2),:);
% Find mEst 
% A 
A_all_real = Gt_all_real * G_all_real;
A_all_realimag = Gt_all_realimag * G_all_realimag;
% B 
% Real only
B_all_real_d2p5e15 = Gt_all_real * d_2p5e15_all_real_v;
B_all_real_d7p5e15 = Gt_all_real * d_7p5e15_all_real_v;
B_all_real_d2p5e16 = Gt_all_real * d_2p5e16_all_real_v;
B_all_real_d7p5e16 = Gt_all_real * d_7p5e16_all_real_v;
B_all_real_d2p5e17 = Gt_all_real * d_2p5e17_all_real_v;
B_all_real_d7p5e17 = Gt_all_real * d_7p5e17_all_real_v;
B_all_real_d2p5e18 = Gt_all_real * d_2p5e18_all_real_v;
B_all_real_d7p5e18 = Gt_all_real * d_7p5e18_all_real_v;
B_all_real_d2p5e19 = Gt_all_real * d_2p5e19_all_real_v;
B_all_real_d7p5e19 = Gt_all_real * d_7p5e19_all_real_v;
% Real and imag (complex magnitude) 
B_all_realimag_d2p5e15 = Gt_all_realimag * d_2p5e15_all_realimag_v;
B_all_realimag_d7p5e15 = Gt_all_realimag * d_7p5e15_all_realimag_v;
B_all_realimag_d2p5e16 = Gt_all_realimag * d_2p5e16_all_realimag_v;
B_all_realimag_d7p5e16 = Gt_all_realimag * d_7p5e16_all_realimag_v;
B_all_realimag_d2p5e17 = Gt_all_realimag * d_2p5e17_all_realimag_v;
B_all_realimag_d7p5e17 = Gt_all_realimag * d_7p5e17_all_realimag_v;
B_all_realimag_d2p5e18 = Gt_all_realimag * d_2p5e18_all_realimag_v;
B_all_realimag_d7p5e18 = Gt_all_realimag * d_7p5e18_all_realimag_v;
B_all_realimag_d2p5e19 = Gt_all_realimag * d_2p5e19_all_realimag_v;
B_all_realimag_d7p5e19 = Gt_all_realimag * d_7p5e19_all_realimag_v;
% mEst 
% Real only 
mEst_all_real_2p5e15 = A_all_real\B_all_real_d2p5e15;
mEst_all_real_7p5e15 = A_all_real\B_all_real_d7p5e15;
mEst_all_real_2p5e16 = A_all_real\B_all_real_d2p5e16;
mEst_all_real_7p5e16 = A_all_real\B_all_real_d7p5e16;
mEst_all_real_2p5e17 = A_all_real\B_all_real_d2p5e17;
mEst_all_real_7p5e17 = A_all_real\B_all_real_d7p5e17;
mEst_all_real_2p5e18 = A_all_real\B_all_real_d2p5e18;
mEst_all_real_7p5e18 = A_all_real\B_all_real_d7p5e18;
mEst_all_real_2p5e19 = A_all_real\B_all_real_d2p5e19;
mEst_all_real_7p5e19 = A_all_real\B_all_real_d7p5e19;
% Real and imag (complex magnitude) 
mEst_all_realimag_2p5e15 = A_all_realimag\B_all_realimag_d2p5e15;
mEst_all_realimag_7p5e15 = A_all_realimag\B_all_realimag_d7p5e15;
mEst_all_realimag_2p5e16 = A_all_realimag\B_all_realimag_d2p5e16;
mEst_all_realimag_7p5e16 = A_all_realimag\B_all_realimag_d7p5e16;
mEst_all_realimag_2p5e17 = A_all_realimag\B_all_realimag_d2p5e17;
mEst_all_realimag_7p5e17 = A_all_realimag\B_all_realimag_d7p5e17;
mEst_all_realimag_2p5e18 = A_all_realimag\B_all_realimag_d2p5e18;
mEst_all_realimag_7p5e18 = A_all_realimag\B_all_realimag_d7p5e18;
mEst_all_realimag_2p5e19 = A_all_realimag\B_all_realimag_d2p5e19;
mEst_all_realimag_7p5e19 = A_all_realimag\B_all_realimag_d7p5e19;
% Results 
concentration_SO2 = [2.5e15; 7.5e15; 2.5e16; 7.5e16; 2.5e17; 7.5e17; 2.5e18; 7.5e18; 2.5e19; 7.5e19]; 
mEst_all_real = [mEst_all_real_2p5e15, mEst_all_real_7p5e15, mEst_all_real_2p5e16, mEst_all_real_7p5e16, mEst_all_real_2p5e17, mEst_all_real_7p5e17, mEst_all_real_2p5e18, mEst_all_real_7p5e18, mEst_all_real_2p5e19, mEst_all_real_7p5e19];  
mEst_all_realimag = [mEst_all_realimag_2p5e15, mEst_all_realimag_7p5e15, mEst_all_realimag_2p5e16, mEst_all_realimag_7p5e16, mEst_all_realimag_2p5e17, mEst_all_realimag_7p5e17, mEst_all_realimag_2p5e18, mEst_all_realimag_7p5e18, mEst_all_realimag_2p5e19, mEst_all_realimag_7p5e19];  
% Plot 
figure % Figure 7 
subplot(2,1,1)
p = plot(concentration_SO2, mEst_all_real(2,:), '^k'); hold on; p.MarkerSize = 8; 
p = plot(concentration_SO2, mEst_all_realimag(2,:), 'vk'); p.MarkerSize = 8; 
plot(concentration_SO2, concentration_SO2, ':k');
ylabel('mEst_S_O_2 (molec/cm^2)')
xlabel('SO_2 column density (molec/cm^2)')
set(gcf,'color','w');
fig = gca; fig.FontSize = 14; 
xlim([-0.1e19 8e19])
ylim([-0.1e19 8e19])
subplot(2,1,2)
p = plot(concentration_SO2, mEst_all_real(2,:), '^k'); hold on; p.MarkerSize = 8; 
p = plot(concentration_SO2, mEst_all_realimag(2,:), 'vk'); p.MarkerSize = 8; 
plot(concentration_SO2, concentration_SO2, ':k');
ylabel('mEst_S_O_2 (molec/cm^2)')
xlabel('SO_2 column density (molec/cm^2)')
set(gcf,'color','w');
fig = gca; fig.FontSize = 14; 
xlim([-0.3e18 8e18])
ylim([-0.3e18 8e18])
% Save convoluted cross sections to be used later 
fname = fullfile(outDir, 'mEst_real_imag');
save(fname); % Saves all variables  
%% Make any adjustments 
% SAVES FIGURE
% fnOut = 'Fig7_mEst_real_imag';
% fname = fullfile(outDirFig, fnOut);
% saveas(gca, fname)
























