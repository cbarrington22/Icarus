% Figures for Introduction - modified for NOVAC spectra

% Fig 1 (a and b) - now convolution of SO2 and O3 for NOVAC spectra 
% a 
% Plots convoluted trace gas cross sections 
inDirCS = '/Users/charlotteb/Documents/Chapter 2/Cross sections/General/Used in figure 1/';
inDirHg = '/Users/charlotteb/Documents/Chapter 2/Cross sections/Hg spectra/'; 
inDirFS = '/Users/charlotteb/Documents/Chapter 2/NOVAC spectra/inFiles/Masaya/D2J2375_140319_1941_0/';
inDirSR = '/Users/charlotteb/Documents/Chapter 2/Cross sections/Vac (for use with reference and NOVAC spectra)/'; 
outDir = '/Users/charlotteb/Documents/Chapter 2/NOVAC spectra/outFiles/'; 
outDirFig = '/Users/charlotteb/Documents/Chapter 2/NOVAC spectra/outFiles/Figures/'; 
% IMPORTANT: The Hg-spectrum MUST be calibrated so that the wavelength information associated with each channel of the spectrometer is correct. All files should be in '.txt' format and contain 
% wavelength information in the first column and absorption OR intensity information the second. For the trace gas cross sections (i and ii) and the solar reference (iii) it is assumed that the
% .txt files DO NOT contain header or footer lines but only the wavelength and absorption OR intensity data. Since the Hg-spectrum should be measured, it is assumed that the .txt file will contain header 
% lines. Please determine the number of header lines (hl) below:
hl = 0; % Change to the row number before start of data. If there are no header lines enter '0'. 
% Reference spectra 
fnHg = 'NOVAC_HG_D2J2375_0.txt'; % Hg-spectrum 
fnSO2 = 'SO2_Bogumil_2003_293K_vac.txt'; 
fnO3 = 'O3_Gorshelev_2014_243K_vac.txt'; 
% fnBRO = 'BrO_Fleischmann_2004_298K_vac.txt'; 
% fnO4 = 'O4_ThalmanVolkamer_2013_293K_air.txt'; 
% fnNO2 = 'NO2_Burrows_1998_293K_vac.txt';
% fnOCLO = 'OClO_Kromminga_2003_293K_vac.txt'; 
% fnCH2O = 'CH2O_MellerMoortgat_2000_298K_air.txt'; 
% DEFINES WAVELENGTH AND ABSOPRTION (OR INTENSITY) 
% Defines file directory, converts to character array and loads wavelength and absorption OR intensity information into the workspace 
dHg = fullfile(inDirHg, fnHg); inHg = char(dHg); % Hg-spectrum 
if hl >= 1
   [wHg, iHg] = textread(inHg,'%f %f', 'headerlines', hl); % Skips header lines if exist 
else
   [wHg, iHg] = textread(inHg,'%f %f'); 
end  
dSO2 = fullfile(inDirCS, fnSO2); inSO2 = char(dSO2); [wSO2, aSO2] = textread(inSO2,'%f %f');
dO3 = fullfile(inDirCS, fnO3); inO3 = char(dO3); [wO3, aO3] = textread(inO3,'%f %f'); 
% dBRO = fullfile(inDirCS, fnBRO); inBRO = char(dBRO); [wBRO, aBRO] = textread(inBRO,'%f %f'); 
% dO4 = fullfile(inDirCS, fnO4); inO4 = char(dO4); [wO4, aO4] = textread(inO4,'%f %f'); 
% dNO2 = fullfile(inDirCS, fnNO2); inNO2 = char(dNO2); [wNO2, aNO2] = textread(inNO2,'%f %f'); 
% dOCLO = fullfile(inDirCS, fnOCLO); inOCLO = char(dOCLO); [wOCLO, aOCLO] = textread(inOCLO,'%f %f');
% dCH2O = fullfile(inDirCS, fnCH2O); inCH2O = char(dCH2O); [wCH2O, aCH2O] = textread(inCH2O,'%f %f');  
% % CONVOLUTION
% Convolutes resampled trace gas cross sections and solar reference with data kernel to match ILS of spectrometer 
figure('Renderer', 'painters', 'Position', [900 900 900 600]) % Plots Hg-spectrum 
dp = plot(wHg, iHg); hold on; dp.LineWidth = 1;
xlim([min(wHg) max(wHg)]); xtickn = (290:10:420); xticks(xtickn);
set(gca,'XMinorTick','on'); 
xlabel('\lambda (nm)'); ylabel('Intensity (counts)');
legend({'Recorded Hg-spectrum','Wavelength range of analysis window'},'FontSize',14);  
title('Selecting Hg-peak to use for data kernel k(x)'); 
Fig = gca; Fig.FontSize = 14; set(gcf,'color','w');
% Defines data kernel k(x)
pL = sprintf('Select a single peak of the Hg-spectrum to be used for the data kernel.\n Select a peak which is not saturated but in close proximity to the analysis window.\n Input the x-axis value corresponding to the start of the selected Hg-peak (l):\n');
l = 288.50; 
pH = sprintf('Input the x-axis value corresponding to the end of the selected Hg-peak (h):\n');
h = 290.00; 
dp = xline(l,'-r'); dp.LineWidth = 1;
dp = xline(h, '-r'); dp.LineWidth = 1;
% Extracts data s(x) from k(x)
[rI, ~] = find(wHg > l&wHg < h); % Returns index of wavelengths between defined limits of kx (iL and iH)
sxW = wHg(rI, :); % Extracts wavelength data for kernel 
sxI = iHg(rI, :); % Extracts intensity data for kernel  
% Plots data kernel s(x) with Smin
figure('Renderer','painters','Position',[900 900 900 600]) 
subplot(1,3,1) % Plots data kernel s(x) with Smin as red cross
dp = plot(sxW, sxI); hold on; dp.LineWidth=1;
[smin, sminInd] = min(sxI);
plot(sxW(sminInd, :), smin, 'xr');
xlim([min(sxW) max(sxW)]); xlabel('\lambda (nm)'); ylabel('Intensity (counts)'); legend({'s(x) from l to h'},'FontSize', 14, 'Location', 'NorthOutside');    
Fig = gca; Fig.FontSize = 15; set(gcf,'color','w');
subplot(1,3,2) % Plots data kernel s(x) - Smin 
sxISmin = sxI - smin;
dp = plot(sxW, sxISmin); hold on; dp.LineWidth=1;
xlim([min(sxW) max(sxW)]);
xlabel('\lambda (nm)'); ylabel('Intensity (counts)'); legend({'s(x) - S_m_i_n from l to h'},'FontSize', 14, 'Location', 'NorthOutside');    
Fig = gca; Fig.FontSize = 15; set(gcf,'color','w');
subplot(1,3,3) % Plots normalised data kernel k(x) =  s(x) - Smin / ∫l^h s(x) - Smin
% 'trapz' computes the approximate integral of Y via the trapezoidal method with unit spacing. If Y is a vector, then trapz(Y) is the approximate integral of Y
kx = sxISmin/(trapz(sxISmin)); 
dp = plot(sxW,kx); hold on; dp.LineWidth = 1;
xlim([min(sxW) max(sxW)]); xlabel('\lambda (nm)'); ylabel('Normalised intensity');
legend({'s(x) - S_m_i_n / ∫_l^h s(x) - S_m_i_n'},'FontSize', 14,'Location', 'NorthOutside');   
t = sgtitle('Defining data kernel k(x) for convolution'); t.FontWeight = 'bold'; 
Fig = gca; Fig.FontSize = 15; set(gcf,'color','w');
% Convolutes resampled trace gas cross sections
lambda = wHg; 
cSO2 = conv(aSO2, kx, 'same'); 
cO3 = conv(aO3, kx, 'same'); 
% cBRO = conv(aBRO, kx, 'same');
% cO4 = conv(aO4, kx, 'same');
% cNO2 = conv(aNO2, kx, 'same');
% cOCLO = conv(aOCLO, kx, 'same');
% cCH2O = conv(aCH2O, kx, 'same');
% % Resamples
SO2 = interp1(wSO2, cSO2, lambda,'linear'); 
O3 = interp1(wO3, cO3, lambda,'linear'); 
% BRO = interp1(wBRO, cBRO, lambda,'linear'); 
% O4 = interp1(wO4, cO4, lambda,'linear'); 
% NO2 = interp1(wNO2, cNO2, lambda,'linear'); 
% OCLO = interp1(wOCLO, cOCLO, lambda,'linear'); 
% CH2O = interp1(wCH2O, cCH2O, lambda,'linear'); 
% To convert to vacuum wavelengths we blueshift (to shorter wavelengths) by 0.1 nm
% wHgCH2O = lambda+0.1; 
% wHgO4 = lambda+0.1; 
% PLOTS RESAMPLED AND CONVOLUTED ABSORPTION CROSS SECTIONS AND SOLAR REFERENCE 
% warning('off','MATLAB:Axes:NegativeDataInLogAxis') % Turns off warning ignoring negative data in log axis 
colSO2 = [1 0.15 0]; % Defines colours for each trace gas
colO3 = [0 0 0];
% colBRO = [0.05 0.19 0.82];
% colO4 = [0.98 0.49 0.03];
% colNO2 = [0.38 0.66 0.01];
% colOCLO = [0 0.72 1];
% colCH2O = [0.70 0.02 0.70];
figure
subplot(2,1,1)
% yyaxis left % Plots convoluted cross sections 
% dp = plot(lambda, BRO, '-'); hold on; dp.LineWidth = 1.5; dp.Color = colBRO;
% dp = plot(lambda, OCLO, '-'); dp.LineWidth = 1.5; dp.Color = colOCLO;
dp = plot(lambda, SO2, '-'); dp.LineWidth = 1.5; dp.Color = colSO2; hold on 
% dp = plot(wHgCH2O, CH2O, '-'); dp.LineWidth = 1.5; dp.Color = colCH2O;
% dp = plot(lambda, NO2, '-'); dp.LineWidth = 1.5; dp.Color = colNO2;
dp = plot(lambda, O3, '-'); dp.LineWidth = 1.5; dp.Color = colO3;
ylim([0 3e-17]);
set(gca, 'YScale', 'log'); 
ylabel('Mean absoprtion (cm^2/molec)', 'FontSize', 12); 
% yyaxis right % Plots O4
% dp = plot(wHgO4, O4, '-'); hold on; dp.LineWidth = 1.5; dp.Color = colO4;
% ylim([0 4.5e-46]); 
% ylabel('Mean absoprtion (cm^5/molec^2)'); set(get(gca,'ylabel'),'rotation', -90);
xlabel('\lambda (nm)'); hold on 
xlim([280 420]);
Fig = gca; Fig.FontSize = 13; set(gcf,'color','w');
Fig.Box = 'on'; 
ax = gca;
ax.YAxis(1).Color = 'k';
% ax.YAxis(2).Color = 'k';
ax.XColor = 'k';

% b - Now plots two spectra from Masaya 
% Plots Fraunhofer spectra 
sky = '00000_0.STD';
dark = '00001_0.STD';
% DEFINES WAVELENGTH AND ABSOPRTION (OR INTENSITY) 
% Defines file directory, converts to character array and loads wavelength and absorption OR intensity information into the workspace 
dSolar30 = fullfile(inDirFS, sky); inSolar30 = char(dSolar30); 
dSolar90 = fullfile(inDirFS, dark); inSolar90 = char(dSolar90);  
dHg = fullfile(inDirHg, fnHg); inHg = char(dHg); % Hg-spectrum 
hlMasaya = 3; 
flMasaya = 2048;
[iMasaya1] = textread(inSolar30,'%s', 'headerlines', hlMasaya); % Skips header lines if exist 
[iMasaya2] = textread(inSolar90,'%s', 'headerlines', hlMasaya); % Skips header lines if exist 
iMasaya1 = iMasaya1(1:flMasaya, :); 
iMasaya2 = iMasaya2(1:flMasaya, :); 
iMasaya1 = str2double(iMasaya1);
iMasaya2 = str2double(iMasaya2);
subplot(2,1,2)
dp = plot(lambda, iMasaya1, '-k'); hold on; dp.LineWidth = 1; % Sky 
dp = plot(lambda, iMasaya2, '-k'); hold on; dp.LineWidth = 1; % Dark
% ylim([0 6.2e04]); 
ylabel('Intensity (counts)', 'FontSize', 12); set(get(gca,'ylabel'),'position', [271.6620125960061,29354.20369232715,-1]);
% x-axis 
xlabel('\lambda (nm)'); hold on 
xlim([280 420]); 
Fig = gca; Fig.FontSize = 13; set(gcf,'color','w');
Fig.Box = 'on'; 
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XColor = 'k';
%% Make any adjustments 
% SAVES FIGURE
fnOut = 'Fig1_CrossSections_NOVAC';
fname = fullfile(outDirFig, fnOut);
saveas(gca, fname)

% % Fig 2 (a and b) 
% lambda = wHg; Use Hg spectrum
chnls = length(lambda); % Number of spectrometer channels  
minlambda = min(lambda); % Starting wavelength 
maxlambda = max(lambda); % Final wavelength 
rangeW = maxlambda - minlambda; % Determines wavelength range of spectrometer 
Fs = rangeW/chnls; % Fs is the sampling rate in wavelength, each data point (or channel) is equivalent to Fs wavelengths 
% CWT 
% Test CWT for d 
iMasaya1(isnan(iMasaya1)) = 0; % Don't extract a wavelength range since there is an offset 
iMasaya1_ln = log(iMasaya1); % Produces error in CWT since log(0) is inf... 
% Therefore must crop by one data point to remove first zero entry 
ind = 2; % Saves in mat file to use in main_NOVAC.m 
% [WTiMasaya1, ~, ~] = cwt(iMasaya1_ln, Fs);  
% SO2 and O3 for deisgn matrix - now for Masaya 
SO2(isnan(SO2)) = 0; O3(isnan(O3)) = 0; % Converts any NAN values to zeros
[WTcSO2, ~, ~] = cwt(SO2, Fs);  
[WTcO3, F, COI] = cwt(O3, Fs);  
% COI
szCoi = length(COI); szF = length(F); matF = repelem(F, 1, szCoi); colCoi = COI'; matCoi = repelem(colCoi, szF, 1); idxCoi = matF <= matCoi; % Removes COI
WTcSO2(idxCoi) = NaN; 
WTcO3(idxCoi) = NaN; 
figure
subplot(2,1,1)
pcolor(lambda, F, real(WTcSO2)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w');
xlim([min(lambda) max(lambda)]); xlabel('\lambda (nm)'); ylabel('Spatial frequency (cycles/\lambda)'); set(gca,'YScale','log');
cb = colorbar; cb.Label.Position(1) = 3.5;
ylabel(cb,'Amplitude', 'Rotation', -90); 
Fig = gca; Fig.FontSize = 13; set(gcf,'color','w');
subplot(2,1,2)
pcolor(lambda, F, real(WTcO3)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w');
xlim([min(lambda) max(lambda)]); xlabel('\lambda (nm)'); ylabel('Spatial frequency (cycles/\lambda)'); set(gca,'YScale','log');
cb = colorbar; cb.Label.Position(1) = 3.5;
ylabel(cb,'Amplitude', 'Rotation', -90); 
Fig = gca; Fig.FontSize = 13; set(gcf,'color','w');
% Save convoluted cross sections to be used later 
%% Make any adjustments 
% SAVES FIGURE
fnOut = 'Fig2_SpatialFrequency_NOVAC';
fname = fullfile(outDirFig, fnOut);
saveas(gca, fname)

% Fig 3 (a and b) 
% Loads solar reference 
fnSolar = 'solar_ChanceKurucz_2010_vac.txt';
inSolar = fullfile(inDirSR, fnSolar); inSolar = char(inSolar); 
[wSolar, aSolar] = textread(inSolar,'%f %f');
solar = interp1(wSolar, aSolar, wHg,'linear'); % Reversed as high resolution 
csolar = conv(solar, kx, 'same');
fnRing = 'Ring_solar_ChanceKurucz_2010_unk.txt';
inRing = fullfile(inDirSR, fnRing); inRing = char(inRing); 
[wRing, aRing] = textread(inRing,'%f %f');
Ring = interp1(wRing, aRing, wHg,'linear'); % Reversed as high resolution 
cRing = conv(Ring, kx, 'same');
figure 
subplot(2,1,1)
dp = plot(wHg, csolar, '-'); hold on; dp.LineWidth = 1; dp.Color = 'k'; 
% ylim([0 3e-17]);
ylabel('Intensity', 'FontSize', 12); 
xlabel('\lambda (nm)'); hold on 
xlim([280 420]);
Fig = gca; Fig.FontSize = 13; set(gcf,'color','w');
Fig.Box = 'on'; 
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XColor = 'k';
subplot(2,1,2)
lncsolar = log(csolar); 
csolar(isnan(lncsolar)) = 0; % Converts any NAN values to zeros
[WTlncsolar, F, COI] = cwt(lncsolar, Fs);  
% COI
szCoi = length(COI); szF = length(F); matF = repelem(F, 1, szCoi); colCoi = COI'; matCoi = repelem(colCoi, szF, 1); idxCoi = matF <= matCoi; % Removes COI
WTlncsolar(idxCoi) = NaN; 
pcolor(lambda, F, real(WTlncsolar)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w');
xlim([min(lambda) max(lambda)]); xlabel('\lambda (nm)'); ylabel('Spatial frequency (cycles/\lambda)'); set(gca,'YScale','log');
cb = colorbar; cb.Label.Position(1) = 3.5;
ylabel(cb,'Amplitude', 'Rotation', -90); 
Fig = gca; Fig.FontSize = 13; set(gcf,'color','w');
% Save convoluted cross sections to be used later 
fname = fullfile(outDir, 'Introduction_Figures_NOVAC');
save(fname); % Saves all variables  
% Changes name
I0 = csolar; 
Ring = cRing; 
fname = fullfile(outDir, 'Fig1_Fig3_conv_I0_SO2_O3_Ring_NOVAC');
save(fname, 'fnHg', 'hlMasaya', 'flMasaya','fnSolar', 'I0', 'fnSO2', 'SO2', 'fnO3', 'O3', 'fnRing' , 'Ring', 'lambda', 'kx', 'iMasaya1', 'ind', 'Fs'); 
%% Make any adjustments 
% SAVES FIGURE
fnOut = 'Fig3_FraunhoferLines_NOVAC';
fname = fullfile(outDirFig, fnOut);
saveas(gca, fname)
 
