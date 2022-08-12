
% Icarus: exploIting the frequenCy signature of trace gAs absoRption in Uv-Spectra

% Appendix A: Preparation of trace gas absorption cross sections and solar reference

% This code forms part of the supplementary material in the following paper: 'Barrington et al., Exploiting the frequency signature of trace gas absorption: A new model to quantify SO2 from UV spectra of volcanic plumes'.

% Written by C. Barrington August 2022
 
% DESCRIPTION of code: 
% This code uses a Hg-spectrum which contains correct wavelength information to resample and convolute absorption cross sections for (i) SO2 and (ii) O3 as well as a (iii) high resolution solar reference taken from the 
% literature.

% INSTRUCTIONS for user: 
% (1) Create an 'Icarus' directory at: /Users/<username>/ e.g., /Users/<username>/Icarus/ and modify the path names below: 
addpath /Users/<user>/Icarus/
inDir = '/Users/<user>/Icarus/inFiles/';
outDir ='/Users/<user>/Icarus/outFiles/';

% (2) Save this script (AppendixA.m) in the main directory: /Users/<username>/Icarus/  

% (3) Create two folders within the main directory, an 'inFile' and 'outFile' e.g., /Users/<username>/Icarus/inFile and /Users/<username>/Icarus/outFile 

% (4) Place the (i) SO2 and (ii) O3 trace gas cross sections, (iii) solar reference taken from the literature together with the (iv) Hg-spectrum inside the 'inFile' folder.
% IMPORTANT: The Hg-spectrum MUST be calibrated so that the wavelength information associated with each channel of the spectrometer is correct. All files should be in '.txt' format and contain 
% wavelength information in the first column and absorption OR intensity information the second. For the trace gas cross sections (i and ii) and the solar reference (iii) it is assumed that the
% .txt files DO NOT contain header or footer lines but only the wavelength and absorption OR intensity data. Since the Hg-spectrum should be measured, it is assumed that the .txt file will contain header 
% lines. Please determine the number of header lines (hl) below:
hl = 14; % Change to the row number before start of data. If there are no header lines enter '0'. 

% (5) Change the file names below to reflect the file names inside 'inFiles': 
fnSO2 = 'SO2_VandaeleHermansFally(2009)_298K_227.275-416.658nm.txt'; % SO2 absorption cross section
fnO3 = 'O3_Burrows(1999)_293K_230-794nm(air).txt'; % O3 absorption cross section
fnSolar = 'Chance-Kurucz-solar2010-JQSRT.txt'; % Solar reference
fnHg = 'FLMS195681__0__12-14-16-726.txt'; % Hg-spectrum 

% (6) Run code 

% (7) Find the prepared absorption cross sections for (i) SO2 and (ii) O3, together with  (iii) solar reference and (iv) wavelength information for (i to iii) saved as .mat file in outDir. 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOT INTENDED TO BE MODIFIED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Loading files..\n'); % Displays message to user 

% DEFINES WAVELENGTH AND ABSOPRTION (OR INTENSITY) 
% Defines file directory, converts to character array and loads wavelength and absorption OR intensity information into the workspace 
dSO2 = fullfile(inDir, fnSO2); inSO2 = char(dSO2); [wSO2, aSO2] = textread(inSO2,'%f %f'); % SO2 absorption cross section
dO3 = fullfile(inDir, fnO3); inO3 = char(dO3); [wO3, aO3] = textread(inO3,'%f %f'); % O3 absorption cross section
dSolar = fullfile(inDir, fnSolar); inSolar = char(dSolar); [wSolar, aSolar] = textread(inSolar,'%f %f'); % Solar reference
dHg = fullfile(inDir, fnHg); inHg = char(dHg); % Hg-spectrum 
if hl >= 1
   [wHg, iHg] = textread(inHg,'%f %f', 'headerlines', hl); % Skips header lines if exist 
else
   [wHg, iHg] = textread(inHg,'%f %f'); 
end  

fprintf('Files loaded. Resampling data..\n'); % Displays message to user 

% RESAMPLES DATA AT WAVELENGTHS DEFINED BY CALIBRATED HG-SPECTRUM
% Uses interp(x, v, xq) to query trace gas cross sections and solar reference at wavelengths defined by calibrated Hg-spectrum   
% 'interp' returns interpolated values of a 1-D function at specific query points (xq) using linear interpolation 
% Vector x contains the sample points - original wavelength information 
% Vector v contains the correpsonding values, v(x) - original absorption (OR intensity) information 
SO2 = interp1(wSO2, aSO2, wHg,'linear'); O3 = interp1(wO3, aO3, wHg,'linear'); solar = interp1(wSolar, aSolar, wHg,'linear'); 

fprintf('Resampling complete.\n'); % Displays message to user 

% CONVOLUTION
% Convolutes resampled trace gas cross sections and solar reference with data kernel to match ILS of spectrometer 
% Defines wavelength of analysis window 
lAw = 315.6; % Lower wavelength limit 
uAw = 326.8; % Upper wavelength limit 
figure('Renderer', 'painters', 'Position', [900 900 900 600]) % Plots Hg-spectrum 
dp = plot(wHg, iHg); hold on; dp.LineWidth = 1;
xlim([min(wHg) max(wHg)]); xtickn = (280:10:420); xticks(xtickn);
set(gca,'XMinorTick','on'); 
dp = xline(lAw, '-k'); dp.LineWidth = 1; % Places vertical line at lower wavelength range of analysis window 
dp = xline(uAw, '-k'); dp.LineWidth = 1; % Places vertical line at upper wavelength range of analysis window 
xlabel('\lambda (nm)'); ylabel('Intensity (counts)');
legend({'Recorded Hg-spectrum','Wavelength range of analysis window'},'FontSize',14);  
title('Selecting Hg-peak to use for data kernel k(x)'); 
Fig = gca; Fig.FontSize = 14; set(gcf,'color','w');
% Defines data kernel k(x)
pL = sprintf('Select a single peak of the Hg-spectrum to be used for the data kernel.\n Select a peak which is not saturated but in close proximity to the analysis window.\n Input the x-axis value corresponding to the start of the selected Hg-peak (l):\n');
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
subplot(1,3,3) % Plots normalised data kernel k(x) =  s(x) - Smin / ∫l^h s(x) - Smin
% 'trapz' computes the approximate integral of Y via the trapezoidal method with unit spacing. If Y is a vector, then trapz(Y) is the approximate integral of Y
kx = sxISmin/(trapz(sxISmin)); 
dp = plot(sxW,kx); hold on; dp.LineWidth = 1;
xlim([min(sxW) max(sxW)]); xlabel('\lambda (nm)'); ylabel('Normalised intensity');
legend({'s(x) - S_m_i_n / ∫_l^h s(x) - S_m_i_n'},'FontSize', 14,'Location', 'NorthOutside');   
t = sgtitle('Defining data kernel k(x) for convolution'); t.FontWeight = 'bold'; 
Fig = gca; Fig.FontSize = 15; set(gcf,'color','w');

fprintf('Convoluting spectra..\n'); % Displays message to user 

% Convolutes resampled trace gas cross sections
% w = conv(u, v, shape) returns a subsection of the convolution, as specified by shape. For example, conv(u, v, 'same') returns only the central part of the convolution, the same size as u
cSO2 = conv(SO2, kx, 'same');
cO3 = conv(O3, kx, 'same');
cSolar = conv(solar, kx, 'same');

% PLOTS RESAMPLED AND CONVOLUTED ABSORPTION CROSS SECTIONS AND SOLAR REFERENCE 
warning('off','MATLAB:Axes:NegativeDataInLogAxis') % Turns off warning ignoring negative data in log axis 
colSO2 = [0.8 0 0]; % Defines colours for each trace gas
colO3 = [0.2 0.1 1];
colsolar = [0 0 0];
figure 
yyaxis left % Plots convoluted cross sections 
dp = plot(wHg, cSO2, '-'); hold on; dp.LineWidth = 1; dp.Color = colSO2;
dp = plot(wHg, cO3, '-'); dp.LineWidth = 1; dp.Color = colO3;
ylim([0 max(cO3)]);
set(gca, 'YScale', 'log'); ylabel('Mean absoprtion (molec/cm^2)'); xlim([min(wHg) max(wHg)]); xlabel('\lambda (nm)');
yyaxis right % Plots convoluted solar reference 
dp = plot(wHg, cSolar, '-'); hold on; dp.LineWidth = 1; dp.Color = colsolar;
ylim([0 max(cSolar)]); ylabel('Intensity (counts)');
dp = xline(lAw, ':k'); dp.LineWidth = 1; % Places vertical line at lower wavelength range of analysis window 
dp = xline(uAw, ':k'); dp.LineWidth = 1; % Places vertical line at upper wavelength range of analysis window 
title('Modified cross sections and solar reference'); 
legend({'SO_2','O_3', 'Solar reference', 'Wavelength range of analysis window'},'FontSize', 14, 'Location', 'eastoutside');Fig = gca; Fig.FontSize = 14; set(gcf,'color','w');

% SAVES MODIFIED ABSORPTION CROSS SECTIONS AND SOLAR REFERENCE 
fnOut = 'AppendixA';
fname = fullfile(outDir, fnOut);
save(fname,'cSO2', 'cO3', 'cSolar', 'wHg'); % Saves (i) SO2, (ii) O3, (iii) solar reference and (iv) wavelength information for (i to iii) 

fprintf('Modified absoprtion cross sections and solar reference saved in %s\n', outDir); % Displays message to user 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
