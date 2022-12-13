% Icarus: exploIting the frequenCy signature of trace gAs absoRption in Uv-Spectra
 
% Appendix A: Convolution of reference spectra and trace gas cross sections 

% This code is provided as part of Chapter 2 of the thesis 'Using Volcanic Gases to Understand Open-vent Volcanoes' submitted to the Nanyang Technological University
% in partial fulfillment of the requirements for the degree of Doctor of Philosophy. 

% Written by C. Barrington December 2022
 
% DESCRIPTION of code: 

% This code uses a Hg spectrum which contains correct wavelength information to resample and convolute (i) a high-resolution solar reference and absorption cross sections for (iI) SO2 and (iii) O3 taken from the literature as 
% well as (iv) a Ring spectrum. 
 
% INSTRUCTIONS for user: 
% (1) Create an 'Icarus' directory at: /Users/<username>/ e.g., /Users/<username>/Icarus/ and modify the path names below: 
% addpath /Users/<username>/Icarus/  
% dir = '/Users/<username>/Icarus/';
% inDir = '/Users/<username>/Icarus/inFiles/';
 
% (2) Save this script (AppendixA.m) in the main directory: /Users/<username>/Icarus/  
 
% (3) Create a folder within the main directory, called 'inFile' e.g., /Users/<username>/Icarus/inFile and /Users/<username>/Icarus/inFile 
 
% (4) Place the following files inside the 'inFile' which was created in step 4: 
% High resolution solar reference 
% SO2 absorption cross section 
% O3 absorption cross section 
% Ring spectrum
% Hg spectrum 

% IMPORTANT: 
% As default, convoluted trace gas cross sections and reference spectra are resampled to the wavelength range of the Hg-spectrum. 
% If the user wishes to resample to an alternative wavelength range, define lambda below 
% All files should be in '.txt' format and contain wavelength information in the 
% first column (in nm) and intensity (for the solar reference, Ring and Hg spectra) OR absorption (for cross sections for SO2 and O3) in the second. It is assumed that the .txt files do not contain header or footer lines except 
% for the Hg spectrum. Since this should be measured, it is expected that the .txt file will contain header lines. Please change 'hlines' according to reflect the number of header lines contained in the Hg spectrum .txt file. 
hlines = 14; % If there are no header lines enter '0'
% lambda = <wavelength>; 
 
% (5) Change the file names below to reflect the file names placed inside 'inFiles': 
fnSolar = 'solar_ChanceKurucz_2010_vac.txt'; % Solar reference
fnSO2 = 'SO2_Bogumil_2003_293K_vac.txt'; % SO2 absorption cross section
fnO3 = 'O3_Gorshelev_2014_243K_vac.txt'; % O3 absorption cross section
fnRing = 'Ring_solar_ChanceKurucz_2010_unk.txt'; % Ring spectrum
fnHg = 'FLMS195681__0__12-14-16-726.txt'; % Hg spectrum 
 
% (6) Run code 
 
% (7) Find the solar reference (I_0), SO2 and O3 trace gas absorption cross sections (σ_so2, σ_o3) and Ring spectrum (Ring) saved as .mat file in: '/Users/<user>/Icarus/outFiles/' and figures in /Users/<user>/Icarus/outFiles/Figures/
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOT INTENDED TO BE MODIFIED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Turns off warnings 
warning('off','MATLAB:MKDIR:DirectoryExists') 
 
fprintf('Creating output directory..\n'); % Displays message to user 
 
% Creates output directories
fnOut1 = 'outFiles'; fname = fullfile(dir, fnOut1); mkdir(fname)
str1 = 'outFiles/'; outDir = [dir str1]; 
fnOut2 = 'Figures'; fname = fullfile(outDir, fnOut2); mkdir(fname); 
str2 = 'Figures/'; outDirFig = [outDir str2]; 
 
fprintf('Loading files from %s..\n', inDir); % Displays message to user 
 
% DEFINES WAVELENGTH AND INTENSITY/ABSORPTION 
% Defines file directory, converts to character array and loads both wavelength and intensity/absorption data into the workspace 
dSolar = fullfile(inDir, fnSolar); inSolar = char(dSolar); [wSolar, aSolar] = textread(inSolar,'%f %f'); % Solar reference
dSO2 = fullfile(inDir, fnSO2); inSO2 = char(dSO2); [wSO2, aSO2] = textread(inSO2,'%f %f'); % SO2 absorption cross section
dO3 = fullfile(inDir, fnO3); inO3 = char(dO3); [wO3, aO3] = textread(inO3,'%f %f'); % O3 absorption cross section
dRing = fullfile(inDir, fnRing); inRing = char(dRing); [wRing, aRing] = textread(inRing,'%f %f'); % Ring spectrum
dHg = fullfile(inDir, fnHg); inHg = char(dHg); % Hg spectrum 
if hlines >= 1
   [wHg, iHg] = textread(inHg,'%f %f', 'headerlines', hlines); % Skips header lines if exist 
else
   [wHg, iHg] = textread(inHg,'%f %f'); 
end  
 
fprintf('Files loaded. Define data kernel..\n'); % Displays message to user 
 
% DEFINES DATA KERNEL k(x) 
figure('Renderer', 'painters', 'Position', [900 900 900 600]) % Plots Hg-spectrum 
title('Hg-spectrum'); 
dp = plot(wHg, iHg); hold on; dp.LineWidth = 1;
xlim([min(wHg) max(wHg)]); xtickn = (290:20:420); xticks(xtickn);
set(gca,'XMinorTick','on'); 
xlabel('\lambda (nm)'); ylabel('Intensity (counts)');
Fig = gca; Fig.FontSize = 14; set(gcf,'color','w');
xlim([280 420]);
pL = sprintf('Select a single peak of the Hg-spectrum to be used for the data kernel.\n Select a peak which is not saturated but in close proximity to the analysis window.\n Input the x-axis value corresponding to the start of the selected Hg-peak (l):\n');
l = input(pL); % e.g., 296.17
pH = sprintf('Input the x-axis value corresponding to the end of the selected Hg-peak (h):\n');
h = input(pH); % e.g., 297.73; 
dp = xline(l,'-r'); dp.LineWidth = 0.5;
dp = xline(h, '-r'); dp.LineWidth = 0.5;
% Extracts data s(x) 
[rI, ~] = find(wHg > l&wHg < h); % Returns index of wavelengths between defined limits l and h
sxW = wHg(rI, :); % Extracts wavelength data for kernel 
sxI = iHg(rI, :); % Extracts intensity data for kernel 
% Plots data kernel s(x)
figure('Renderer','painters','Position',[900 900 900 600]) 
title('Convolution data kernel k(x)'); 
[smin, ~] = min(sxI);
sxISmin = sxI - smin;
kx = sxISmin/(trapz(sxISmin)); 
dp = plot(sxW, kx); hold on; dp.LineWidth = 1;
xtickn = (290:0.5:420); xticks(xtickn);
xlim([min(sxW) max(sxW)])
xlabel('\lambda (nm)'); ylabel('Normalised s(x)');
Fig = gca; Fig.FontSize = 14; set(gcf,'color','w');

fprintf('Convoluting trace gas absorption cross sections..\n'); % Displays message to user 

% RESAMPLING 
% Checks if lambda is defined in the workspace 
if exist('lambda','var') == 1
else 
    lambda = wHg; % If lambda is not defined by the user, trace gas cross sections and reference spectra are resampled to the wavelength range of the Hg-spectrum
end 
% CONVOLUTION OF TRACE GAS CROSS SECTIONS
cSO2 = conv(aSO2, kx, 'same'); 
cO3 = conv(aO3, kx, 'same'); 
cBRO = conv(aBRO, kx, 'same');
cO4 = conv(aO4, kx, 'same');
cNO2 = conv(aNO2, kx, 'same');
cOCLO = conv(aOCLO, kx, 'same');
cCH2O = conv(aCH2O, kx, 'same');
fprintf('Resampling convoluted cross sections ..\n'); % Displays message to user 
% RESAMPLING OF TRACE GAS CROSS SECTIONS 
SO2 = interp1(wSO2, cSO2, lambda,'linear'); 
O3 = interp1(wO3, cO3, lambda,'linear'); 
BRO = interp1(wBRO, cBRO, lambda,'linear'); 
O4 = interp1(wO4, cO4, lambda,'linear'); 
NO2 = interp1(wNO2, cNO2, lambda,'linear'); 
OCLO = interp1(wOCLO, cOCLO, lambda,'linear'); 
CH2O = interp1(wCH2O, cCH2O, lambda,'linear'); 
fprintf('Convoluting and resampling high resolution solar reference and Ring spectrum..\n'); % Displays message to user 
% RESAMPLING OF HIGH RESOLUTION SOLAR REFERENCE (AND RING)
solar = interp1(wSolar, aSolar, lambda,'linear'); 
csolar = conv(solar, kx, 'same');
Ring = interp1(wRing, aRing, lambda,'linear'); 
cRing = conv(Ring, kx, 'same');
% CHANGES VARIABLE NAME 
I0 = csolar; 
Ring = cRing; 

% PLOTS I0, SO2, O3 AND RING 
figure
subplot(2,2,1)
plot(lambda, I0, '-r'); hold on; 
ylabel('Intensity')
xlabel('\lambda (nm)'); 
xlim([280 420]);
Fig = gca; Fig.Box = 'on'; Fig.FontSize = 14; set(gcf,'color','w');
subplot(2,2,2)
plot(lambda, SO2, '-r'); hold on; 
ylabel('Mean absoprtion (cm^2/molec)'); 
xlabel('\lambda (nm)'); 
xlim([280 420]);
Fig = gca; Fig.Box = 'on'; Fig.FontSize = 14; set(gcf,'color','w');
subplot(2,2,3)
plot(lambda, O3, '-r'); hold on; 
ylabel('Mean absoprtion (cm^2/molec)'); 
xlabel('\lambda (nm)'); 
xlim([280 420]);
Fig = gca; Fig.Box = 'on'; Fig.FontSize = 14; set(gcf,'color','w');
subplot(2,2,4)
plot(lambda, Ring, '-r'); hold on; 
ylabel('Intensity')
xlabel('\lambda (nm)'); 
xlim([280 420]);
Fig = gca; Fig.Box = 'on'; Fig.FontSize = 14; set(gcf,'color','w');
% Saves figure 
fnOut = 'Output.png';
fname = fullfile(outDirFig, fnOut);
saveas(gcf, fname)

% SAVES MODIFIED ABSORPTION CROSS SECTIONS AND SOLAR REFERENCE 
fnOut = 'AppendixA';
fname = fullfile(outDir, fnOut);
save(fname, 'fnHg', 'fnSolar', 'I0', 'fnSO2', 'SO2', 'fnO3', 'O3', 'fnRing' , 'Ring', 'lambda', 'kx'); 

fprintf('Convoluted trace gas cross sections and reference spectra saved in %s\n', outDir); % Displays message to user 
 
% clearvars 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
