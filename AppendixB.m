% Icarus: exploIting the frequenCy signature of trace gAs absoRption in Uv-Spectra
 
% Appendix B: Converting prepared solar reference (I_0), trace gas absorption cross sections (σ_so2, σ_o3) and Ring spectrum (Ring) to spatial frequency (I_0', σ_so2', σ_o3' and Ring') for use in the linear model 
 
% This code is provided as part of the following paper: 'Barrington et al., Exploiting the frequency signature of trace gas absorption: A new model to quantify SO2 from UV spectra of volcanic plumes'
 
% Written by C. Barrington August 2022
 
% DESCRIPTION of code: 
% This code determines the CWT of the prepared trace gas cross sections and solar reference produced in Appendix A and creates the deisgn matrix (G) of the linear model. 
 
% INSTRUCTIONS for user: 
% NOTE: users who have run AppendixA.m may skip ahead to step (3).
% (1) Create an 'Icarus' directory at: /Users/<username>/ e.g., /Users/<username>/Icarus/ and modify the path names below: 
% addpath /Users/charlotteb/Icarus/
addpath '/Users/charlotteb/Desktop/Frontiers SIRS manuscript submission/'
% dir = '/Users/charlotteb/Icarus/';
dir = '/Users/charlotteb/Desktop/Frontiers SIRS manuscript submission/Code (Github)/Icarus/';
% inDir = '/Users/charlotteb/Icarus/inFiles/';
inDir = '/Users/charlotteb/Desktop/Frontiers SIRS manuscript submission/Code (Github)/Icarus/inFiles/';
 
% (2) If you have not used AppendixA.m to prepare trace gas cross sections and solar reference, please add a .mat file named 'AppendixA' to the path '/Users/charlotteb/Icarus/inFiles/' which contains the following variables: 
    % (a) 'solar' - solar reference which matches the wavelength range and Instrument Line Shape (ILS) of the spectrometer used to record the measurement spectra to be analysed (I_0)
    % (b) 'SO2' - absorption cross section of SO2 which matches the wavelength range and Instrument Line Shape (ILS) of the spectrometer used to record the measurement spectra to be analysed (σ_so2)
    % (c) 'O3' - absorption cross section of O3 which matches the wavelength range and Instrument Line Shape (ILS) of the spectrometer used to record the measurement spectra to be analysed (σ_o3)
    % (d) 'ring' - Rings spectrum which matches the wavelength range and Instrument Line Shape (ILS) of the spectrometer used to record the measurement spectra to be analysed (Ring)
    % (e) 'wHg' - wavelength information of the spectrometer used to record the measurement spectra to be analysed
    % (f) 'lAw' - lower wavelength limit of analysis window 
    % (g) 'uAw' - upper wavelength limit of analysis window
% NOTE: Variables a to e should be a single vector whose length is equivalent to the number of spectrometer channels. Wavelength information in e should be correct for a to d. And should correspond to channels of the spectrometer 
% which record positive intensities (e.g., not less than 280 nm). Variables f and g are single values representing the wavelength range of analysis window. These are set to 315.6 and 326.8 respectively in AppendixA.m. 
 
% (3) Save this script (AppendixB.m) in the main directory: /Users/<username>/Icarus/  
 
% (4) Run code
 
% (6) Find the components of and the created deisgn matrix (G): I_0', σ_so2', σ_o3' and Ring' wavelength and frequency information, sampling frequency (Fs) limits of
% analysis window used and indices (both wavelength and frequency) together with G, the transpose of G (Gt) as .mat file in:'/Users/<user>/Icarus/outFiles/' and figures in /Users/<user>/Icarus/outFiles/Figures/
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOT INTENDED TO BE MODIFIED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
fprintf('Creating output directory..\n'); % Displays message to user 
 
% Creates output directories
fnOut1 = 'outFiles'; fname = fullfile(dir, fnOut1); mkdir(fname)
str1 = 'outFiles/'; outDir = [dir str1]; 
fnOut2 = 'Figures'; fname = fullfile(outDir, fnOut2); mkdir(fname); 
str2 = 'Figures/'; outDirFig = [outDir str2]; 
 
fprintf('Loading solar reference (I_0), trace gas absorption cross sections (σ_so2, σ_o3) and Ring spectrum (Ring) from AppendixA.mat..\n'); % Displays message to user 
 
% LOADS MODIFIED CROSS SECTIONS AND SOLAR REFERENCE WITH WAVELENGTH INFORMATION 
fnIn = 'AppendixA.mat'; fname = fullfile(outDir, fnIn); load(fname); 
solar(isnan(solar)) = 0; SO2(isnan(SO2)) = 0; O3(isnan(O3)) = 0; ring(isnan(ring)) = 0; % Converts any NAN values to zeros
 
fprintf('Determining sampling frequency..\n'); % Displays message to user 
 
% DETERMINES SAMPLING FREQUENCY (Fs)
% [...,F] = cwt(...,Fs) specifies the sampling frequency, Fs, in hertz as a positive scalar and returns the scale-to-frequency conversions in hertz, F. 
% If the sampling frequency is not specified, cwt returns F in cycles/sample (i.e. channel). 
% If the sampling frequency is specified (in wavelength), cwt returns F in cycles/wavelength. 
chnls = length(wHg); % Number of spectrometer channels  
minW = min(wHg); % Starting wavelength 
maxW = max(wHg); % End wavelength 
rangeW = maxW - minW; % Determines wavelength range of spectrometer 
Fs = rangeW/chnls; % Fs is the sampling rate in wavelength, each data point (or channel) is equivalent to Fs wavelengths 
 
fprintf('Computing CWT..\n'); % Displays message to user 
 
% COMPUTES CWT
figure
subplot(2,2,1) 
[wtsolar, F, COI] = cwt(solar, Fs); % Solar reference
szCoi = length(COI); szF = length(F); matF = repelem(F, 1, szCoi); colCoi = COI'; matCoi = repelem(colCoi, szF, 1); idxCoi = matF <= matCoi; % Removes COI
wtsolar(idxCoi) = NaN; 
pcolor(wHg, F, abs(wtsolar)); shading interp; colorbar; colormap(parula(100)); set(gcf,'color','w');
xlim([min(wHg) max(wHg)]); xlabel('\lambda (nm)'); ylabel('Spatial frequency (cycles/\lambda)'); set(gca,'YScale','log');
cb = colorbar; 
ylabel(cb,'Solar intensity', 'FontSize', 16.5, 'Rotation', -90); 
cb.Label.Position(1) = 3.5;
Fig = gca; Fig.FontSize = 16.5; set(gcf,'color','w');
title('Solar reference')
 
subplot(2,2,2) 
[wtSO2, F, COI] = cwt(SO2, Fs); % SO2
wtSO2(idxCoi) = NaN; 
pcolor(wHg, F, abs(wtSO2)); shading interp; colorbar; colormap(parula(100)); set(gcf,'color','w');
xlim([min(wHg) max(wHg)]); xlabel('\lambda (nm)'); ylabel('Spatial frequency (cycles/\lambda)'); set(gca,'YScale','log');
cb = colorbar; 
ylabel(cb,'Magnitude of absorption (cm^2/molec)', 'FontSize', 16.5, 'Rotation', -90); 
cb.Label.Position(1) = 4.5;
Fig = gca; Fig.FontSize = 16.5; set(gcf,'color','w');
title('SO_2')
 
subplot(2,2,3) 
[wtO3, ~ , ~] = cwt(O3, Fs); % O3
wtO3(idxCoi) = NaN; 
pcolor(wHg, F, abs(wtO3)); shading interp; colorbar; colormap(parula(100)); set(gcf,'color','w');
xlim([min(wHg) max(wHg)]); xlabel('\lambda (nm)'); ylabel('Spatial frequency (cycles/\lambda)'); set(gca,'YScale','log');
cb = colorbar; 
ylabel(cb,'Magnitude of absorption (cm^2/molec)', 'FontSize', 16.5, 'Rotation', -90); 
cb.Label.Position(1) = 4.5;
Fig = gca; Fig.FontSize = 16.5; set(gcf,'color','w');
title('O_3')
 
subplot(2,2,4) 
[wtring, ~, ~] = cwt(ring, Fs); % Solar reference
wtring(idxCoi) = NaN; 
pcolor(wHg, F, abs(wtring)); shading interp; colorbar; colormap(parula(100)); set(gcf,'color','w');
xlim([min(wHg) max(wHg)]); xlabel('\lambda (nm)'); ylabel('Spatial frequency (cycles/\lambda)'); set(gca,'YScale','log');
cb = colorbar; 
ylabel(cb,'Solar intensity', 'FontSize', 16.5, 'Rotation', -90); 
cb.Label.Position(1) = 3.5;
Fig = gca; Fig.FontSize = 16.5; set(gcf,'color','w');
title('Ring')
% Saves figure
fnOut = 'CWT.png';
fname = fullfile(outDirFig, fnOut);
saveas(gcf, fname)
 
fprintf('Extracting analysis window..\n'); % Displays message to user 
 
% DEFINES ANALYSIS WINDOW
% WAVELENGTH 
% Finds channel (wavelength row index) which is closest to upper and lower wavelength limits of the analysis window 
[~, idxWlAw] = min(abs(wHg - lAw)); % Lower 
WlAw = wHg(idxWlAw); 
[~, idxWuAw] = min(abs(wHg - uAw)); % Upper
WuAw = wHg(idxWuAw); 
% FREQUENCY
% Frequency levels in F are descending, so low row indices represent higher frequencies and higher row indices represent low frequencies 
Aw1 = 0.01; % First limit on frequency (highest frequency)
Aw2 = 0.001; % Second limit on frequency (lowest frequency) 
[~, idxF1Aw] = min(abs(F - Aw1)); % Highest 
F1Aw = F(idxF1Aw); 
[~, idxF2Aw] = min(abs(F - Aw2)); % Lowest 
F2Aw = F(idxF2Aw); 
% EXTRACTS ANALYSIS WINDOW
wtsolar = abs(wtsolar(idxF1Aw:idxF2Aw,idxWlAw:idxWuAw)); 
wtSO2 = abs(wtSO2(idxF1Aw:idxF2Aw,idxWlAw:idxWuAw)); 
wtO3 = abs(wtO3(idxF1Aw:idxF2Aw,idxWlAw:idxWuAw)); 
wtring = abs(wtring(idxF1Aw:idxF2Aw,idxWlAw:idxWuAw)); 
wHgOrig = wHg; % Saves original wHg for use in preparing d
wHg = wHg(idxWlAw:idxWuAw); % Corresponding wavelength information
Forig = F; % Saves original F for use in preparing d
F = F(idxF1Aw:idxF2Aw); % Corresponding frequency information
 
% PLOTS I_0', σ_so2', σ_o3' and Ring'
figure
subplot(2,2,1) 
pcolor(wHg, F, wtsolar); shading interp; colorbar; colormap(parula(100)); set(gcf,'color','w');
xlim([min(wHg) max(wHg)]); xlabel('\lambda (nm)'); ylabel('Spatial frequency (cycles/\lambda)'); set(gca,'YScale','log');
cb = colorbar; 
ylabel(cb,'Solar intensity', 'FontSize', 16.5, 'Rotation', -90); 
cb.Label.Position(1) = 3.5;
Fig = gca; Fig.FontSize = 16.5; set(gcf,'color','w');
title('I_0 prime')
 
subplot(2,2,2) 
pcolor(wHg, F, wtSO2); shading interp; colorbar; colormap(parula(100)); set(gcf,'color','w');
xlim([min(wHg) max(wHg)]); xlabel('\lambda (nm)'); ylabel('Spatial frequency (cycles/\lambda)'); set(gca,'YScale','log');
cb = colorbar; 
ylabel(cb,'Magnitude of absorption (cm^2/molec)', 'FontSize', 16.5, 'Rotation', -90); 
cb.Label.Position(1) = 4.5;
Fig = gca; Fig.FontSize = 16.5; set(gcf,'color','w');
title('SO_2 prime')
 
subplot(2,2,3) 
pcolor(wHg, F, wtO3); shading interp; colorbar; colormap(parula(100)); set(gcf,'color','w');
xlim([min(wHg) max(wHg)]); xlabel('\lambda (nm)'); ylabel('Spatial frequency (cycles/\lambda)'); set(gca,'YScale','log');
cb = colorbar; 
ylabel(cb,'Magnitude of absorption (cm^2/molec)', 'FontSize', 16.5, 'Rotation', -90); 
cb.Label.Position(1) = 4.5;
Fig = gca; Fig.FontSize = 16.5; set(gcf,'color','w');
title('O_3 prime')
 
subplot(2,2,4) 
pcolor(wHg, F, wtring); shading interp; colorbar; colormap(parula(100)); set(gcf,'color','w');
xlim([min(wHg) max(wHg)]); xlabel('\lambda (nm)'); ylabel('Spatial frequency (cycles/\lambda)'); set(gca,'YScale','log');
cb = colorbar; 
ylabel(cb,'Solar intensity', 'FontSize', 16.5, 'Rotation', -90); 
cb.Label.Position(1) = 3.5;
Fig = gca; Fig.FontSize = 16.5; set(gcf,'color','w');
title('Ring prime')
% Saves figure
fnOut = 'G.png';
fname = fullfile(outDirFig, fnOut);
saveas(gcf, fname)
 
fprintf('Analysis window extracted. Vectorising CWT data ready for incorporation in design matrix (G)\n'); % Displays message to user 
 
% VECTORISES
wtsolar = wtsolar(:); wtSO2 = wtSO2(:); wtO3 = wtO3(:); wtring = wtring(:);
 
% REMOVES ROWS CORRESPONDING TO COI (NaN) 
wtsolar = wtsolar(all(~isnan(wtsolar), 2),:); wtSO2 = wtSO2(all(~isnan(wtSO2), 2),:); wtO3 = wtO3(all(~isnan(wtO3), 2),:); wtring = wtring(all(~isnan(wtring), 2),:);
 
fprintf('Creating design matrix (G)..\n'); % Displays message to user 
 
% DESIGN MATRIX (G)
G = [wtsolar, wtSO2, wtO3, wtring];
szG = size(G); % Size of G
fprintf('G is matrix %d by %d\n', szG(1,1), szG(1,2)); % Displays message to user 
Gt = G'; % Transpose of G
 
% SAVES MODIFIED ABSORPTION CROSS SECTIONS AND SOLAR REFERENCE 
fnOut = 'AppendixB';
fname = fullfile(outDir, fnOut);
save(fname,'wtsolar', 'wtSO2', 'wtO3', 'wtring', 'wHgOrig', 'wHg', 'Forig', 'F', 'Fs', 'lAw', 'uAw', 'F1Aw', 'F2Aw', 'idxCoi', 'idxF1Aw', 'idxF2Aw', 'idxWlAw', 'idxWuAw', 'G', 'Gt'); % Saves I_0', σ_so2', σ_o3' and Ring', wavelength 
% and frequency information (both original and analysis window), sampling frequency (Fs), range of the analysis window limits and indices (both wavelength and frequency) and the COI index, together with G, the transpose of G (Gt) as .mat file in: 
% '/Users/<user>/Icarus/outFiles/' and figures in /Users/<user>/Icarus/outFiles/Figures/
 
fprintf('Design matrix (G) and associated variables saved in %s\n', outDir); % Displays message to user 
 
clearvars
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%