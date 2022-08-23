% Icarus: exploIting the frequenCy signature of trace gAs absoRption in Uv-Spectra
 
% Appendix A: Preparation of the solar reference (I_0), SO2 and O3 trace gas absorption cross sections (σ_so2, σ_o3) and Ring spectrum (Ring)
 
% This code is provided as part of the following paper: 'Barrington et al., Exploiting the frequency signature of trace gas absorption: A new model to quantify SO2 from UV spectra of volcanic plumes'
 
% Written by C. Barrington August 2022
 
% DESCRIPTION of code: 
% This code uses a Hg spectrum which contains correct wavelength information to resample and convolute (i) a high-resolution solar reference and absorption cross sections for (iI) SO2 and (iii) O3 taken from the literature as 
% well as (iv) a Ring spectrum. 
 
% INSTRUCTIONS for user: 
% (1) Create an 'Icarus' directory at: /Users/<username>/ e.g., /Users/<username>/Icarus/ and modify the path names below: 
addpath /Users/<user>/Icarus/
dir = '/Users/<user>/Icarus/';
inDir = '/Users/<user>/Icarus/inFiles/';
 
% (2) Save this script (AppendixA.m) in the main directory: /Users/<username>/Icarus/  
 
% (3) Create a folder within the main directory, called 'inFile' e.g., /Users/<username>/Icarus/inFile and /Users/<username>/Icarus/inFile 
 
% (4) Place the following files inside the 'inFile' which was created in step 4: 
% High resolution solar reference 
% SO2 absorption cross section 
% O3 absorption cross section 
% Ring spectrum
% Hg spectrum 
 
% IMPORTANT: The Hg spectrum MUST be calibrated so that the wavelength information associated with each channel of the spectrometer is correct. All files should be in '.txt' format and contain wavelength information in the 
% first column (in nm) and intensity (for the solar reference, Ring and Hg spectra) OR absorption (for cross sections for SO2 and O3) in the second. It is assumed that the .txt files do not contain header or footer lines except 
% for the Hg spectrum. Since this should be measured, it is expected that the .txt file will contain header lines. Please change 'hlines' according to reflect the number of header lines contained in the Hg spectrum .txt file. 
hlines = 14; % If there are no header lines enter '0'
 
% (5) Change the file names below to reflect the file names placed inside 'inFiles': 
fnSolar = 'SAO2010_Chance_Kurucz.txt'; % Solar reference
fnSO2 = 'SO2_VandaeleHermansFally(2009)_298K_227.275-416.658nm.txt'; % SO2 absorption cross section
fnO3 = 'O3_Burrows(1999)_293K_230-794nm(air).txt'; % O3 absorption cross section
fnRing = 'Ring_SAO2010_Chance_Kurucz.txt'; % Ring spectrum
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
 
fprintf('Files loaded. Now resampling data according to wavelength information in Hg spectrum..\n'); % Displays message to user 
 
% RESAMPLES DATA AT WAVELENGTHS DEFINED BY CALIBRATED Hg spectrum
% Crops wHg to limit wavelengths to those 280 nm or above
[ind, ~] = find(wHg >= 280.00);
wHg = wHg(ind, 1); 
iHg = iHg(ind, 1); 
 
% Uses interp(x, v, xq) to query trace gas cross sections and solar reference at wavelengths defined by calibrated Hg spectrum   
% 'interp' returns interpolated values of a 1-D function at specific query points (xq) using linear interpolation 
% Vector x contains the sample points - original wavelength information 
% Vector v contains the corresponding values, v(x) - original absorption (OR intensity) information 
solar = interp1(wSolar, aSolar, wHg,'linear'); SO2 = interp1(wSO2, aSO2, wHg,'linear'); O3 = interp1(wO3, aO3, wHg,'linear'); ring = interp1(wRing, aRing, wHg,'linear'); 
 
fprintf('Resampling complete.\n'); % Displays message to user 
 
% CONVOLUTION
% Convolutes resampled trace gas cross sections and solar reference with data kernel to match ILS of spectrometer 
% Defines wavelength of analysis window 
lAw = 315.6; % Lower wavelength limit 
uAw = 326.8; % Upper wavelength limit 
figure('Renderer', 'painters', 'Position', [900 900 900 600]) % Plots Hg spectrum 
dp = plot(wHg, iHg); hold on; dp.LineWidth = 1;
xlim([min(wHg) max(wHg)]); xtickn = (280:10:420); xticks(xtickn);
set(gca,'XMinorTick','on'); 
dp = xline(lAw, '-k'); dp.LineWidth = 1; % Places vertical line at lower wavelength range of analysis window 
dp = xline(uAw, '-k'); dp.LineWidth = 1; % Places vertical line at upper wavelength range of analysis window 
xlabel('\lambda (nm)'); ylabel('Intensity (counts)');
legend({'Recorded Hg spectrum','Wavelength range of analysis window'},'FontSize',14);  
title('Selecting Hg-peak to use for data kernel k(x)'); 
Fig = gca; Fig.FontSize = 14; set(gcf,'color','w');
% Defines data kernel k(x)
pL = sprintf('Select a single peak of the Hg spectrum to be used for the data kernel.\n Select a peak which is not saturated but in close proximity to the analysis window.\n Input the x-axis value corresponding to the start of the selected Hg-peak (l):\n');
l = input(pL);
pH = sprintf('Input the x-axis value corresponding to the end of the selected Hg-peak (h):\n');
h = input(pH);
dp = xline(lAw,'-k'); dp.LineWidth = 1;
dp = xline(l,'-r'); dp.LineWidth = 1;
dp = xline(uAw, '-k'); dp.LineWidth = 1;
dp = xline(h, '-r'); dp.LineWidth = 1;
 
fprintf('Preparing data kernel..\n'); % Displays message to user 
 
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
subplot(1,3,3) % Plots normalised data kernel k(x) =  s(x) - Smin / ‚à´l^h s(x) - Smin
% 'trapz' computes the approximate integral of Y via the trapezoidal method with unit spacing. If Y is a vector, then trapz(Y) is the approximate integral of Y
kx = sxISmin/(trapz(sxISmin)); 
dp = plot(sxW,kx); hold on; dp.LineWidth = 1;
xlim([min(sxW) max(sxW)]); xlabel('\lambda (nm)'); ylabel('Normalised intensity');
legend({'s(x) - S_m_i_n / ‚à´_l^h s(x) - S_m_i_n'},'FontSize', 14,'Location', 'NorthOutside');   
t = sgtitle('Defining data kernel k(x) for convolution'); t.FontWeight = 'bold'; 
Fig = gca; Fig.FontSize = 15; set(gcf,'color','w');
% Saves figure
fnOut = 'K(x).png';
fname = fullfile(outDirFig, fnOut);
saveas(gcf, fname)
 
fprintf('Convoluting spectra..\n'); % Displays message to user 
 
% Convolutes resampled trace gas cross sections
% w = conv(u, v, shape) returns a subsection of the convolution, as specified by shape. For example, conv(u, v, 'same') returns only the central part of the convolution, the same size as u
cSolar = conv(solar, kx, 'same'); cSO2 = conv(SO2, kx, 'same'); cO3 = conv(O3, kx, 'same'); cRing = conv(ring, kx, 'same');
 
% PLOTS RESAMPLED AND CONVOLUTED ABSORPTION CROSS SECTIONS AND SOLAR REFERENCE 
warning('off','MATLAB:Axes:NegativeDataInLogAxis') % Turns off warning ignoring negative data in log axis 
colSO2 = [0.8 0 0]; % Defines colours for each trace gas
colO3 = [0.2 0.1 1];
colsolar = [0 0 0];
figure 
subplot(2,2,1)
dp = plot(wHg, cSolar, '-'); hold on; dp.LineWidth = 1; dp.Color = colsolar;
ylim([0 max(cSolar)]);
ylabel('Intensity (counts)'); xlim([min(wHg) max(wHg)]); xlabel('\lambda (nm)');
dp = xline(lAw, ':k'); dp.LineWidth = 1; % Places vertical line at lower wavelength range of analysis window 
dp = xline(uAw, ':k'); dp.LineWidth = 1; % Places vertical line at upper wavelength range of analysis window 
title('I_0'); 
Fig = gca; Fig.FontSize = 14; set(gcf,'color','w');
 
subplot(2,2,2)
dp = plot(wHg, cSO2, '-'); hold on; dp.LineWidth = 1; dp.Color = colSO2; % SO2
ylim([0 max(cSO2)]);
ylabel('Mean absorption (molec/cm^2)'); xlim([min(wHg) max(wHg)]); xlabel('\lambda (nm)');
dp = xline(lAw, ':k'); dp.LineWidth = 1; % Places vertical line at lower wavelength range of analysis window 
dp = xline(uAw, ':k'); dp.LineWidth = 1; % Places vertical line at upper wavelength range of analysis window 
title('\sigma_S_O_2'); 
Fig = gca; Fig.FontSize = 14; set(gcf,'color','w');
 
subplot(2,2,3)
dp = plot(wHg, cO3, '-'); hold on; dp.LineWidth = 1; dp.Color = colO3; % O3
ylim([0 max(cO3)]);
ylabel('Mean absorption (molec/cm^2)'); xlim([min(wHg) max(wHg)]); xlabel('\lambda (nm)');
dp = xline(lAw, ':k'); dp.LineWidth = 1; % Places vertical line at lower wavelength range of analysis window 
dp = xline(uAw, ':k'); dp.LineWidth = 1; % Places vertical line at upper wavelength range of analysis window 
title('\sigma_O_3');
Fig = gca; Fig.FontSize = 14; set(gcf,'color','w');
 
subplot(2,2,4)
dp = plot(wHg, cRing, '-'); hold on; dp.LineWidth = 1; dp.Color = colsolar;
ylim([0 max(cRing)]);
ylabel('Intensity (counts)'); xlim([min(wHg) max(wHg)]); xlabel('\lambda (nm)');
dp = xline(lAw, ':k'); dp.LineWidth = 1; % Places vertical line at lower wavelength range of analysis window 
dp = xline(uAw, ':k'); dp.LineWidth = 1; % Places vertical line at upper wavelength range of analysis window 
title('Ring'); 
Fig = gca; Fig.FontSize = 14; set(gcf,'color','w');
% Saves figure 
fnOut = 'Output.png';
fname = fullfile(outDirFig, fnOut);
saveas(gcf, fname)
 
% SAVES MODIFIED ABSORPTION CROSS SECTIONS AND SOLAR REFERENCE 
fnOut = 'AppendixA';
fname = fullfile(outDir, fnOut);
 
solar = cSolar; SO2 = cSO2; O3 = cO3; ring = cRing; % Re-defines variable names for saving 
 
save(fname,'SO2', 'O3', 'solar', 'ring', 'wHg', 'lAw', 'uAw', 'hlines', 'ind', 'l', 'h', 'kx'); % Saves I_0, σ_so2, σ_o3 and Ring as well as wavelength information, lower and upper wavelength range of the analysis window, number of 
% header lines, index for wavelengths used, wavelengths used to produce kx and kx itself
 
fprintf('Modified absorption cross sections and solar reference saved in %s\n', outDir); % Displays message to user 
 
clearvars 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
