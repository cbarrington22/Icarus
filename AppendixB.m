% Icarus: exploIting the frequenCy signature of trace gAs absoRption in Uv-Spectra

% Appendix B: Preparation of trace gas cross sections and solar reference for design matrix G

% This code forms part of the supplementary material in the following paper: 'Barrington et al., Exploiting the frequency signature of trace gas absorption: A new model to quantify SO2 from UV spectra of volcanic plumes'.

% Written by C. Barrington August 2022
 
% DESCRIPTION of code: 
% This code determines the CWT of the prepared trace gas cross sections and solar reference produced in Appendix A, ready to be included in the deisgn matrix (G) of the linear model. 

% INSTRUCTIONS for user: 
% NOTE: users who have run AppendixA.m may skip ahead to step (3).
% (1) Create an 'Icarus' directory at: /Users/<username>/ e.g., /Users/<username>/Icarus/ and modify the path names below: 
addpath /Users/<user>/Icarus/
outDir = '/Users/<user>/Icarus/outFiles/';

% (2) If you have not used AppendixA.m to prepare trace gas cross sections and solar reference, please add a .mat file with the name 'AppendixA' to outDir which contains the following variables: 
% (a) 'cSO2'absorption cross section of SO2 which matches the wavelength range and Instrument Line Shape (ILS) of the spectrometer used to record the measurement spectra to be analysed
% (b) 'cO3' absorption cross section of O3 which matches the wavelength range and Instrument Line Shape (ILS) of the spectrometer used to record the measurement spectra to be analysed
% (c) 'cSolar' solar reference which matches the wavelength range and Instrument Line Shape (ILS) of the spectrometer used to record the measurement spectra to be analysed
% (d) 'wHg' wavelength information 
% (e) 'lAw' lower wavelength limit of analysis window 
% (f) 'uAw' upper wavelength limit of analysis window
% NOTE: Variables a to d should be a single vector whos length is equivalent to the number of spectrometer channels. Variables e and f are single values representing the wavelength range of analysis window. These are set to 
% 315.6 and 326.8 respectively in AppendixA.m. 

% (3) Save this script (AppendixB.m) in the main directory: /Users/<username>/Icarus/  

% (4) If any initial channels of the spectrometer record negative values, indiate the channel number at which only posiive values follow.
% For example if the first 5 channels of the spectrometer produce 0 -3 5 -5 -2 5 followed by only positive values, 6 should be assigned to cropVal.
% NOTE: This step is to allow all spectra, cross sections and the solar reference to be cropped, in order for the CWT to be computed using positive values only. This value may be easily found by looking at the intensity data 
% from the Hg-spectrum.  
cropVal = 0; 

% (5) Run code

% (6) Find the components of deisgn matrix (G): sigma_SO2, sigma_O3 and I_solar saved as .mat file in outDir. 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOT INTENDED TO BE MODIFIED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Loading modified trace gas cross sections and solar reference..\n'); % Displays message to user 

% LOADS MODIFIED CROSS SECTIONS AND SOLAR REFERENCE WITH WAVELENGTH INFORMATION 
fnIn = 'AppendixA.mat';
fname = fullfile(outDir, fnIn);
load(fname); 
SO2(isnan(cSO2)) = 0; O3(isnan(cO3)) = 0; solar(isnan(cSolar)) = 0; % Converts any NAN values to zeros

fprintf('Determining sampling frequency..\n'); % Displays message to user 

% DETERMINES SAMPLING FREQUENCY (Fs)
% [...,F] = cwt(...,Fs) specifies the sampling frequency, Fs, in hertz as a positive scalar and returns the scale-to-frequency conversions in hertz, F. 
% If the sampling frequency is not specified, cwt returns F in cycles/sample (i.e. channel). 
% If the sampling frequnecy is specified (in wavelength), cwt returns F in cycles/wavelength. 
chnls = length(wHg); % Number of spectrometer channels  
minW = min(wHg); % Starting wavelength 
maxW = max(wHg); % End wavelength 
rangeW = maxW - minW; % Determines wavelength range of spectrometer 
Fs = rangeW/chnls; % Fs is the sampling rate in wavelength, each data point (or channel) is equivalent to Fs wavelengths 

fprintf('Cropping data to exclude any negative values..\n'); % Displays message to user 

% CROPS DATA TO EXCLUDE NEGATIVE VALUES 
if cropVal >0    
    SO2 = SO2(cropVal:end,:); % SO2 cross section
    O3 = O3(cropVal:end,:); % O3 cross section
    solar = solar(cropVal:end,:); % Solar reference
    wHg = wHg(cropVal:end,:); % Wavelength information  
else  
end 

fprintf('Computing CWT..\n'); % Displays message to user 

% COMPUTES CWT
figure
subplot(1,3,1)
[wtSO2, F, COI] = cwt(SO2, Fs); % SO2
szCoi = length(COI); szF = length(F); matF = repelem(F, 1, szCoi); colCoi = COI'; matCoi = repelem(colCoi, szF, 1); idxCoi = matF <= matCoi; % Removes COI
wtSO2(idxCoi) = NaN; 
pcolor(wHg, F, abs(wtSO2)); shading interp; colorbar; colormap(parula(100)); set(gcf,'color','w');
xlim([min(wHg) max(wHg)]); xlabel('\lambda (nm)'); ylabel('Spatial frequency (cycles/\lambda)'); set(gca,'YScale','log');
cb = colorbar; 
ylabel(cb,'Magnitude of absorption (cm^2/molec)', 'FontSize', 16.5, 'Rotation', -90); 
cb.Label.Position(1) = 4.5;
Fig = gca; Fig.FontSize = 16.5; set(gcf,'color','w');
title('SO_2')
subplot(1,3,2)
[wtO3, ~ , ~] = cwt(O3, Fs); % O3
wtO3(idxCoi) = NaN; 
pcolor(wHg, F, abs(wtO3)); shading interp; colorbar; colormap(parula(100)); set(gcf,'color','w');
xlim([min(wHg) max(wHg)]); xlabel('\lambda (nm)'); ylabel('Spatial frequency (cycles/\lambda)'); set(gca,'YScale','log');
cb = colorbar; 
ylabel(cb,'Magnitude of absorption (cm^2/molec)', 'FontSize', 16.5, 'Rotation', -90); 
cb.Label.Position(1) = 4.5;
Fig = gca; Fig.FontSize = 16.5; set(gcf,'color','w');
title('O_3')
subplot(1,3,3) 
[wtSolar, ~, ~] = cwt(solar, Fs); % Solar reference
wtSolar(idxCoi) = NaN; 
pcolor(wHg, F, abs(wtSolar)); shading interp; colorbar; colormap(parula(100)); set(gcf,'color','w');
xlim([min(wHg) max(wHg)]); xlabel('\lambda (nm)'); ylabel('Spatial frequency (cycles/\lambda)'); set(gca,'YScale','log');
cb = colorbar; 
ylabel(cb,'Solar intensity', 'FontSize', 16.5, 'Rotation', -90); 
cb.Label.Position(1) = 3.5;
Fig = gca; Fig.FontSize = 16.5; set(gcf,'color','w');
title('Solar reference')

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
wtSO2 = wtSO2(idxF1Aw:idxF2Aw,idxWlAw:idxWuAw); 
wtO3 = wtO3(idxF1Aw:idxF2Aw,idxWlAw:idxWuAw); 
wtSolar = wtSolar(idxF1Aw:idxF2Aw,idxWlAw:idxWuAw); 
wHg = wHg(idxWlAw:idxWuAw); % Corresponding wavelength information

fprintf('Analysis window extracted. Vectorising CWT data ready for incorporation in design matrix (G)\n'); % Displays message to user 

% VECTORISES
wtSO2 = wtSO2(:);
wtO3 = wtO3(:);
wtSolar = wtSolar(:);

% REMOVES ROWS CORRESPONDING TO COI (NaN) 
wtSO2 = wtSO2(all(~isnan(wtSO2), 2),:);
wtO3 = wtO3(all(~isnan(wtO3), 2),:);
wtSolar = wtSolar(all(~isnan(wtSolar), 2),:);

% SAVES MODIFIED ABSORPTION CROSS SECTIONS AND SOLAR REFERENCE 
fnOut = 'AppendixB';
fname = fullfile(outDir, fnOut);
save(fname,'wtSO2', 'wtO3', 'wtSolar', 'wHg', 'lAw', 'uAw', 'F1Aw', 'F2Aw'); % Saves CWT of (i) SO2, (ii) O3, (iii) solar reference and (iv) wavelength information for (i to iii) as well as the upper and lower wavelength ranges 
% of the analysis window (both wavelength and frequency)

fprintf('Prepared CWT of absorption cross sections and solar reference saved in %s\n', outDir); % Displays message to user 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
