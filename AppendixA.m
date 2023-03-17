% Icarus: exploIting the frequenCy signature of trace gAs absoRption in Uv-Spectra
 
% (Appendix A) Appendix B: Steps for assembling the design matrix G

% This code is provided as part of Chapter 2 of the thesis 'Using Volcanic Gases to Understand Open-vent Volcanoes' submitted to the Nanyang Technological University
% in partial fulfillment of the requirements for the degree of Doctor of Philosophy. 

% Written by C. Barrington December 2022
 
% DESCRIPTION of code: 

% This code assembles design matrix G from I0, SO2, O3, Ring and Bshift (followng AppendixA.m and shiftSpectrum.m)
 
% INSTRUCTIONS for user: 
% (1) This script follows on from AppendixA.m and will load the variables saved in AppendixA.mat
inDir = '/Users/<username>/Icarus/outFiles/';

% (2) In addition, the uthe Bshift spectrum returned from running shiftSpectrum.m should be saved within AppendixA.mat

% (3) Define the wavelength range over which the CWT will be computed 
l = min(lambda); % lower wavelength limit  
h = max(lambda); % upper wavelength limit 
% IMPORTANT: If the measurement spectrum (which will form d) has negative intensity counts in this wavelength range, and error will occur when computing the CWT

% (4) Define analysis window to be used for the lienar model:
% Wavelength range (AW1 as default): 
AW_L = 310; % Low 
AW_H = 340; % High 
% Spatial frequency range (default for Fs = 0.0700):
AW_Ft = 0.01; % Highest 
AW_Fb = 0.0015; % Lowest 
% IMPORTANT: Spatial frequency will depend on the resolution of the measured spectra and may require the user to plot d and SO2 prime first in order to select the appropriate spatial frequency range 
 
% (5) Run code 
 
% (6) Find deisgn matrix (G) for both the real and complex magnitude saved as .mat file in: '/Users/<user>/Icarus/outFiles/designMatrix/'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOT INTENDED TO BE MODIFIED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Creates output directories
fnOut1 = 'designMatrix'; fname = fullfile(inDir, fnOut1); mkdir(fname)
str1 = 'designMatrix/'; outDir = [dir str1]; 

% EXTRACTS WAVELENGTHS 
[rI, ~] = find(lambda > l & lambda < h); 
lambda = lambda(rI,:); 
I0 = I0(rI,:); 
SO2 = SO2(rI,:); 
O3 = O3(rI,:); 
Ring = Ring(rI,:); 
Bshift = Bshift(rI,:); 
% Defines Fs
chnls = length(lambda); % Number of spectrometer channels  
minlambda = min(lambda); % Starting wavelength 
maxlambda = max(lambda); % Final wavelength 
rangeW = maxlambda - minlambda; % Determines wavelength range of spectrometer 
Fs = rangeW/chnls; % Fs is the sampling rate in wav

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

% DEFINES ANALYSIS WINDOW 
% WAVELENGTH 
[~, idxL]= min(abs(lambda-AW_L));
wavL = lambda(idxL); % Low wavelength limit of AW
[~, idxH]= min(abs(lambda-AW_H));
wavH = lambda(idxH);  % High wavelength limit of AW
% SPATIAL FREQUENCY 
[~, idxFt]= min(abs(F-AW_Ft));
Ft = F(idxFt); % Low wavelength limit of AW
[~, idxFb]= min(abs(F-AW_Fb));
Fb = F(idxFb);  % High wavelength limit of AW

% EXTRACTS ANALYSIS WINDOW 
F = F(idxFt:idxFb);
lambda = lambda(idxL:idxH);
WT_I0 = WT_I0(idxFt:idxFb, idxL:idxH);
WT_SO2 = WT_SO2(idxFt:idxFb, idxL:idxH);
WT_O3 = WT_O3(idxFt:idxFb, idxL:idxH);
WT_Ring = WT_Ring(idxFt:idxFb, idxL:idxH);
WT_Bshift = WT_Bshift(idxFt:idxFb, idxL:idxH);

% DEFINES REAL AND IMAGINARY PARTS OF CWT 
WT_I0_real  = real(WT_I0); 
WT_I0_imag  = imag(WT_I0); 
WT_SO2_real  = real(WT_SO2); 
WT_SO2_imag  = imag(WT_SO2); 
WT_O3_real  = real(WT_O3); 
WT_O3_imag  = imag(WT_O3); 
WT_Ring_real  = real(WT_Ring); 
WT_Ring_imag  = imag(WT_Ring); 
WT_Bshift_real  = real(WT_Bshift); 
WT_Bshift_imag  = imag(WT_Bshift); 

% VECTORISES 
% Real
WT_I0_real_v = WT_I0_real(:); 
WT_SO2_real_v = WT_SO2_real(:); 
WT_O3_real_v = WT_O3_real(:); 
WT_Ring_real_v = WT_Ring_real(:); 
WT_Bshift_real_v = WT_Bshift_real(:); 
% imag
WT_I0_imag_v = WT_I0_imag(:); 
WT_SO2_imag_v = WT_SO2_imag(:); 
WT_O3_imag_v = WT_O3_imag(:); 
WT_Ring_imag_v = WT_Ring_imag(:); 
WT_Bshift_imag_v = WT_Bshift_imag(:); 
% Combines real and imaginary part for complex magnitude 
WT_I0_realimag = [WT_I0_real_v; WT_I0_imag_v]; 
WT_SO2_realimag = [WT_SO2_real_v; WT_SO2_imag_v]; 
WT_O3_realimag = [WT_O3_real_v; WT_O3_imag_v]; 
WT_Ring_realimag = [WT_Ring_real_v; WT_Ring_imag_v]; 
WT_Bshift_realimag = [WT_Bshift_real_v; WT_Bshift_imag_v]; 

% DESIGN MATRIX, G 
G_real = [WT_I0_real_v, WT_SO2_real_v, WT_O3_real_v, WT_Ring_real_v, WT_Bshift_real_v]; 
G_realimag = [WT_I0_realimag, WT_SO2_realimag, WT_O3_realimag, WT_Ring_realimag, WT_Bshift_realimag]; 
% Remove NaN 
G_real = G_real(all(~isnan(G_real), 2),:);
G_realimag = G_realimag(all(~isnan(G_realimag), 2),:);
% G transponse 
Gt_real = G_real'; 
Gt_realimag = G_realimag'; 

% SAVES variables 
fname = fullfile(outDir, 'AppendixB');
save(fname, 'G_real', 'G_realimag', 'Gt_real', 'Gt_realimag', 'rI', 'lambda', 'WT_I0', 'WT_SO2', 'WT_O3', 'WT_Ring', 'WT_Bshift', 'l', 'h' , 'AW_L', 'AW_H', 'AW_Ft', 'AW_Fb', 'idxCoi', 'wavL', 'wavH', 'Ft', 'Fb','idxFt', 'idxFb', 'idxL', 'idxH'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







