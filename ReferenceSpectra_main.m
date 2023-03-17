% Reference spectra 
outDir = '/Users/charlotteb/Documents/Chapter 2/Reference spectra /outFiles/';
outDirFig = '/Users/charlotteb/Documents/Chapter 2/Reference spectra /Figures/';

% LOADS EXAMPLE d TO CHECK: 
% Extracted wavelength range (no negative intensity counts) 
% Spatial frequency range (determined based on spatial frequencies of SO2 cross setion and high frequency noise in d) 
inDird = '/Users/charlotteb/Documents/Chapter 2/Reference spectra /inFiles/Spectra/'; % Directory to recorded spectra
fnInd1 = '30_0_FLMS195681__0__11-26-46-904.txt'; % Examples to load 
fnInd2 = '30_2080_FLMS195681__99__11-30-10-698.txt'; 
fnInd3 = '90_0_FLMS195681__0__10-13-49-165.txt'; 
fnInd4 = '90_2080_FLMS195681__99__10-16-47-535.txt'; % Target 

ind1 = fullfile(inDird, fnInd1); 
ind2 = fullfile(inDird, fnInd2); 
ind3 = fullfile(inDird, fnInd3); 
ind4 = fullfile(inDird, fnInd4); 

[wId1, iId1,] = textread(ind1,'%s%s'); 
[wId2, iId2,] = textread(ind2,'%s%s');
[wId3, iId3,] = textread(ind3,'%s%s');
[wId4, iId4,] = textread(ind4,'%s%s');

wId1 = str2double(wId1); iId1 = str2double(iId1); d_1 = [wId1, iId1]; 
wId2 = str2double(wId2); iId2 = str2double(iId2); d_2 = [wId2, iId2]; 
wId3 = str2double(wId3); iId3 = str2double(iId3); d_3 = [wId3, iId3]; 
wId4 = str2double(wId4); iId4 = str2double(iId4); d_4 = [wId4, iId4]; 

figure
plot(d_1(:,1), d_1(:, 2)); hold on 
plot(d_2(:,1), d_2(:, 2));
plot(d_3(:,1), d_3(:, 2));
plot(d_4(:,1), d_4(:, 2));
ylim([0 7e4]); 

%% Visual check: all values > 0 after 279.0 nm 
% We use 280 to 420 to extract information 
l = 280;
h = 420; 
dLambda = d_1(:,1); 
[rId, ~] = find(dLambda > l & dLambda < h); 
dLambda = dLambda(rId,:); 
d_1 = d_1(rId,:); 
d_2 = d_2(rId,:); 
d_3 = d_3(rId,:); 
d_4 = d_4(rId,:); 

figure
plot(d_1(:,1), d_1(:, 2)); hold on 
plot(d_2(:,1), d_2(:, 2));
plot(d_3(:,1), d_3(:, 2));
plot(d_4(:,1), d_4(:, 2));
ylim([0 7e4]); 

% Takes intensity information 
d_1i = d_1(:, 2);
d_2i = d_2(:, 2);
d_3i = d_3(:, 2);
d_4i = d_4(:, 2);

% Check CWT computes 
% Take ln  
d_1ln = log(d_1i);
d_2ln = log(d_2i);
d_3ln = log(d_3i);
d_4ln = log(d_4i);

figure
plot(d_1(:, 1), d_1ln); hold on 
plot(d_2(:, 1), d_2ln);
plot(d_3(:, 1), d_3ln);
plot(d_4(:, 1), d_4ln);

% CWT 
%% Find Fs below 0.0700
[WT_d_1ln, Fd, ~] = cwt(d_1ln, Fs); 
[WT_d_2ln, ~, ~] = cwt(d_2ln, Fs); 
[WT_d_3ln, ~, ~] = cwt(d_3ln, Fs); 
[WT_d_4ln, ~, ~] = cwt(d_4ln, Fs); 

% Additional figure 
% Checks spatial frequency limits - will be Fig. 17 
figure 
subplot(3,2,1)
title('No gas cell'); hold on; 
plot(d_3(:,1), d_3(:, 2));
xlim([min(dLambda) max(dLambda)])
box on
subplot(3,2,2)
title('2080 ppm·m'); hold on; 
plot(d_4(:,1), d_4(:, 2));
xlim([min(dLambda) max(dLambda)])
box on
% Defines limits (used in next figure) to plot box 
AW_Ld = 300;
AW_Hd = 330;
AW_Ftd = 0.01; 
AW_Fbd = 0.0015;
pos = [AW_Ld AW_Fbd AW_Hd-AW_Ld AW_Ftd-AW_Fbd]; 
subplot(3,2,3)
title('No gas cell'); hold on; 
pcolor(dLambda, Fd, real(WT_d_3ln)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
xlim([min(dLambda) max(dLambda)]); set(gca,'YScale','log'); ylabel('Spatial frequency (cycles/\lambda)', 'FontSize', 14);
set(gcf,'color','w'); ylim([min(Fd) max(Fd)]);
rectangle('Position',pos)
box on
subplot(3,2,4)
title('2080 ppm·m'); hold on; 
pcolor(dLambda, Fd, real(WT_d_4ln)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
xlim([min(dLambda) max(dLambda)]); set(gca,'YScale','log');
set(gcf,'color','w'); ylim([min(Fd) max(Fd)]);
rectangle('Position',pos)
cb = colorbar; cb.Label.Position(1) = 3.5;
ylabel(cb,'Amplitude', 'Rotation', -90, 'FontSize', 10);
box on
% Zoom in 
subplot(3,2,5)
title('No gas cell'); hold on; 
pcolor(dLambda, Fd, real(WT_d_3ln)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
xlim([AW_Ld AW_Hd]); set(gca,'YScale','log'); 
set(gcf,'color','w'); ylim([AW_Fbd AW_Ftd]);
xlabel('\lambda (nm)', 'FontSize', 14); 
box on
subplot(3,2,6)
title('2080 ppm·m'); hold on; 
pcolor(dLambda, Fd, real(WT_d_4ln)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
xlim([AW_Ld AW_Hd]); set(gca,'YScale','log'); 
set(gcf,'color','w'); ylim([AW_Fbd AW_Ftd]);
box on
%% Make any adjustments 
% SAVES FIGURE
% fnOut = 'Fig17_dExample_referenceSpectra';
% fname = fullfile(outDirFig, fnOut);
% saveas(gca, fname)
% SAVES COMPLETE CWT 
fname = fullfile(outDir, 'dExample');
save(fname) 

% LOADS COMPONENTS OF G 
% G 
% Run Introduction_Figures.m and FigD1.m first 
inDirG = '/Users/charlotteb/Documents/Chapter 2/Reference spectra /Bshift/outFiles/'; 
fnInG = 'Fig1_Fig3_FigD1_conv_I0_SO2_O3_Ring_Bshift.mat'; 
inFile = fullfile(inDirG, fnInG); 
load(inFile)

% PREPARES EACH ELEMENT OF G (APPENDIX B)  
chnls = length(lambda); % Number of spectrometer channels  
minlambda = min(lambda); % Starting wavelength 
maxlambda = max(lambda); % Final wavelength 
rangeW = maxlambda - minlambda; % Determines wavelength range of spectrometer 
Fs = rangeW/chnls; % Fs is the sampling rate in wav

% EXTRACTS WAVELENGTHS 
% min(lambda) max(lambda) 
% l = min(lambda);
% h = max(lambda); 
[rI, ~] = find(lambda > l & lambda < h); 
lambda = lambda(rI,:); 
I0 = I0(rI,:); 
SO2 = SO2(rI,:); 
O3 = O3(rI,:); 
Ring = Ring(rI,:); 
Bshift = Bshift(rI,:); 

% APPLIES LOG OR NEGATIVE 
I0 = log(I0); % ln 
SO2 = -(SO2); % Negative 
O3 = -(O3); % Negative 

% CWT 
[WT_I0, F, COI] = cwt(I0, Fs); 
[WT_SO2, ~, ~] = cwt(SO2, Fs); 
[WT_O3, ~, ~] = cwt(O3, Fs); 
[WT_Ring, ~, ~] = cwt(Ring, Fs); 
[WT_Bshift, ~, ~] = cwt(Bshift, Fs); 
% REMOVES COI
szCoi = length(COI); szF = length(F); matF = repelem(F, 1, szCoi); colCoi = COI'; matCoi = repelem(colCoi, szF, 1); idxCoi = matF <= matCoi; % Removes COI
WT_I0(idxCoi) = NaN; 
WT_SO2(idxCoi) = NaN; 
WT_O3(idxCoi) = NaN; 
WT_Ring(idxCoi) = NaN; 
WT_Bshift(idxCoi) = NaN; 

% SAVES COMPLETE CWT 
fname = fullfile(outDir, 'WT_G');
save(fname, 'fnInd4', 'WT_I0', 'WT_SO2', 'WT_O3', 'WT_Ring', 'WT_Bshift', 'lambda', 'F', 'COI' , 'Fs', 'rI', 'l', 'h', 'inFile', 'idxCoi'); 

% % DEFINES ANALYSIS WINDOW 
% % WAVELENGTH 
% AW_L = 310;
% AW_H = 340;
% [~, idxL]= min(abs(lambda-AW_L));
% wavL = lambda(idxL); % Low wavelength limit of AW
% [~, idxH]= min(abs(lambda-AW_H));
% wavH = lambda(idxH);  % High wavelength limit of AW
% % SPATIAL FREQUENCY 
% AW_Ft = 0.01; % max(F);
% AW_Fb = 0.0015; % min(F);
% [~, idxFt]= min(abs(F-AW_Ft));
% Ft = F(idxFt); % Low wavelength limit of AW
% [~, idxFb]= min(abs(F-AW_Fb));
% Fb = F(idxFb);  % High wavelength limit of AW
% 
% % EXTRACTS ANALYSIS WINDOW 
% F = F(idxFt:idxFb);
% lambda = lambda(idxL:idxH);
% WT_I0 = WT_I0(idxFt:idxFb, idxL:idxH);
% WT_SO2 = WT_SO2(idxFt:idxFb, idxL:idxH);
% WT_O3 = WT_O3(idxFt:idxFb, idxL:idxH);
% WT_Ring = WT_Ring(idxFt:idxFb, idxL:idxH);
% WT_Bshift = WT_Bshift(idxFt:idxFb, idxL:idxH);
% 
% figure % figure 16 
% subplot(3,2,1) 
% subtitle('J_0(\lambda) prime'); hold on; 
% pcolor(lambda, F, real(WT_I0)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
% xlim([min(lambda) max(lambda)]); set(gca,'YScale','log'); ylim([min(F) max(F)]);
% set(gcf,'color','w'); 
% 
% subplot(3,2,2) 
% subtitle('\sigma_S_O_2 prime'); hold on; 
% pcolor(lambda, F, real(WT_SO2)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
% xlim([min(lambda) max(lambda)]); set(gca,'YScale','log'); ylim([min(F) max(F)]);
% set(gcf,'color','w'); 
% 
% subplot(3,2,3) 
% subtitle('\sigma_O_3 prime'); hold on; 
% pcolor(lambda, F, real(WT_O3)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
% xlim([min(lambda) max(lambda)]); set(gca,'YScale','log'); ylabel('Spatial frequency (cycles/\lambda)', 'FontSize', 14);
% set(gcf,'color','w'); ylim([min(F) max(F)]);
% 
% subplot(3,2,4) 
% subtitle('Ring prime'); hold on; 
% pcolor(lambda, F, real(WT_Ring)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
% xlim([min(lambda) max(lambda)]); set(gca,'YScale','log'); ylim([min(F) max(F)]);
% cb = colorbar; cb.Label.Position(1) = 3.5;
% ylabel(cb,'Amplitude', 'Rotation', -90, 'FontSize', 10);
% set(gcf,'color','w'); 
% 
% subplot(3,2,5) 
% subtitle('B_s_h_i_f_t prime'); hold on; 
% pcolor(lambda, F, real(WT_Bshift)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
% xlim([min(lambda) max(lambda)]); set(gca,'YScale','log'); ylim([min(F) max(F)]);
% set(gcf,'color','w'); xlabel('\lambda (nm)', 'FontSize', 14); 
% 
% %% Make any adjustments 
% % SAVES FIGURE
% % fnOut = 'Fig16_G_referenceSpectra';
% % fname = fullfile(outDirFig, fnOut);
% % saveas(gca, fname)