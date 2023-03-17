% Stray light offset   

outDir = '/Users/charlotteb/Documents/Chapter 2/NOVAC spectra/outFiles/';
outDirFig = '/Users/charlotteb/Documents/Chapter 2/NOVAC spectra/outFiles/Figures/';

% LOADS EXAMPLE d TO CHECK: - this is done for one scan in Introduction.m 
% Extracted wavelength range (no negative intensity counts) 
% Spatial frequency range (determined based on spatial frequencies of SO2 cross setion and high frequency noise in d) 
% inDird = '/Users/charlotteb/Documents/Chapter 2/Reference spectra /inFiles/Spectra/'; % Directory to recorded spectra
% Run Introduction_Figures.m first to I0 
inDir = '/Users/charlotteb/Documents/Chapter 2/NOVAC spectra/outFiles/'; 
fnIn = 'Fig1_Fig3_conv_I0_SO2_O3_Ring_NOVAC_straylight.mat'; 
inFile = fullfile(inDir, fnIn); 
load(inFile)
d_1 = iMasaya1; 
dLambda = lambda; % Here we don't have additional wavelength data in spectrum 
% fnInd1 = '30_0_FLMS195681__0__11-26-46-904.txt'; % Examples to load 
% fnInd2 = '30_2080_FLMS195681__99__11-30-10-698.txt'; 
% fnInd3 = '90_0_FLMS195681__0__10-13-49-165.txt'; 
% fnInd4 = '90_2080_FLMS195681__99__10-16-47-535.txt'; % Target 
% 
% ind1 = fullfile(inDird, fnInd1); 
% ind2 = fullfile(inDird, fnInd2); 
% ind3 = fullfile(inDird, fnInd3); 
% ind4 = fullfile(inDird, fnInd4); 
% 
% [wId1, iId1,] = textread(ind1,'%s%s'); 
% [wId2, iId2,] = textread(ind2,'%s%s');
% [wId3, iId3,] = textread(ind3,'%s%s');
% [wId4, iId4,] = textread(ind4,'%s%s');
% 
% wId1 = str2double(wId1); iId1 = str2double(iId1); d_1 = [wId1, iId1]; 
% wId2 = str2double(wId2); iId2 = str2double(iId2); d_2 = [wId2, iId2]; 
% wId3 = str2double(wId3); iId3 = str2double(iId3); d_3 = [wId3, iId3]; 
% wId4 = str2double(wId4); iId4 = str2double(iId4); d_4 = [wId4, iId4]; 

% figure
% plot(d_1(:,1), d_1(:, 2)); hold on 
% plot(d_2(:,1), d_2(:, 2));
% plot(d_3(:,1), d_3(:, 2));
% plot(d_4(:,1), d_4(:, 2));
% ylim([0 7e4]); 

% Crop value already determined in Introduction_NOVAC.m 
% %% Visual check: all values > 0 after 279.0 nm 
% % We use 280 to 420 to extract information - now using all data but leave here incase encounter issue later (all NOVAC data have offset)
% l = min(lambda);
% h = max(lambda); 
% [rId, ~] = find(dLambda >= l & dLambda <= h); % Use all wavelengths 
dLambda = dLambda(ind:end); 
d_1 = d_1(ind:end); 

figure
subplot(4,2,1)
title('No offset'); hold on; 
plot(dLambda, d_1); hold on 
offsetInd = (30:149); 
% plot(d_2(:,1), d_2(:, 2));
% plot(d_3(:,1), d_3(:, 2));
% plot(d_4(:,1), d_4(:, 2));
ylim([0 5e4]); 
xlim([min(dLambda) max(dLambda)]); 
ylabel('Intensity (counts)')
set(gcf,'color','w'); 
box on 

subplot(4,2,2)
offset = mean(d_1(offsetInd));
d_1Offset = d_1-offset;
title('Offset removed'); hold on; 
plot(dLambda, d_1Offset); hold on 
ylim([-1e4 5e4]); 
xlim([min(dLambda) max(dLambda)]); 
box on 
% ylabel('Intensity (counts)')
% Takes intensity information 
% d_1i = d_1(:, 2);

% Check CWT computes 
% Take ln  
d_1ln = log(d_1);
d_1Offsetln = log(d_1Offset);

subplot(4,2,3)
plot(dLambda, d_1ln); hold on 
xlim([min(dLambda) max(dLambda)]); 
ylabel('Intensity (ln(counts))')
ylim([0 12]); 

subplot(4,2,4)
plot(dLambda, d_1Offsetln); hold on 
% ylim([-1e4 5e4]); 
xlim([min(dLambda) max(dLambda)]); 
% ylabel('Intensity (ln(counts))')
ylim([0 12]); 

% figure
% plot(dLambda, d_1ln); hold on 
d_1ln(isnan(d_1ln)) = 0; 
d_1Offsetln(isnan(d_1Offsetln)) = 0; 

% figure 

% CWT 
%% Find Fs below 0.0700
[WT_d_1ln, Fd, ~] = cwt(d_1ln, Fs); 
[WT_d_1Offsetln, ~, ~] = cwt(real(d_1Offsetln), Fs); 

% AW_Ld = dLambda(1,1);
% AW_Hd = dLambda(end,1);
AW_Ftd = 0.01;
AW_Fbd = 0.0015;
pos = [300 AW_Fbd 40 AW_Ftd-AW_Fbd]; 

subplot(4,2,5)
pcolor(dLambda, Fd, real(WT_d_1ln)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
xlim([min(dLambda) max(dLambda)]); set(gca,'YScale','log'); % ylabel('Spatial frequency (cycles/\lambda)', 'FontSize', 14);
rectangle('Position',pos)
ylabel('Spatial frequency (cycles/\lambda)', 'FontSize', 14);

subplot(4,2,6)
pcolor(dLambda, Fd, real(WT_d_1Offsetln)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
xlim([min(dLambda) max(dLambda)]); set(gca,'YScale','log'); % ylabel('Spatial frequency (cycles/\lambda)', 'FontSize', 14);
rectangle('Position',pos)
cb = colorbar; cb.Label.Position(1) = 3.5;
ylabel(cb,'Amplitude', 'Rotation', -90); 

lambda = dLambda;
AW_L = 310;
AW_H = 340;

% DEFINES ANALYSIS WINDOW 
% WAVELENGTH 
[~, idxL]= min(abs(lambda-AW_L));
wavL = lambda(idxL); % Low wavelength limit of AW
[~, idxH]= min(abs(lambda-AW_H));
wavH = lambda(idxH);  % High wavelength limit of AW
% % SPATIAL FREQUENCY 
AW_Ft = 0.01; % max(F);
AW_Fb = 0.0015; % min(F);
[~, idxFt]= min(abs(Fd-AW_Ft));
Ft = Fd(idxFt); % Low wavelength limit of AW
[~, idxFb]= min(abs(Fd-AW_Fb));
Fb = Fd(idxFb);  % High wavelength limit of AW
% 
% % EXTRACTS ANALYSIS WINDOW 
F = Fd(idxFt:idxFb);
lambda = lambda(idxL:idxH);
WT_d_1ln = WT_d_1ln(idxFt:idxFb, idxL:idxH);
WT_d_1Offsetln = WT_d_1Offsetln(idxFt:idxFb, idxL:idxH);

subplot(4,2,7)
pcolor(lambda, F, real(WT_d_1ln)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
xlim([min(lambda) max(lambda)]); set(gca,'YScale','log'); % ylabel('Spatial frequency (cycles/\lambda)', 'FontSize', 14);
subplot(4,2,8)
pcolor(lambda, F, real(WT_d_1Offsetln)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
xlim([min(lambda) max(lambda)]); set(gca,'YScale','log'); % ylabel('Spatial frequency (cycles/\lambda)', 'FontSize', 14);
xlabel('\lambda (nm)', 'FontSize', 14);




