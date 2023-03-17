
% Final figure 

% Plotting results for NOVAC analysis 
outDir = '/Users/charlotteb/Documents/Chapter 2/NOVAC spectra/Results/';
outDirFig = '/Users/charlotteb/Documents/Chapter 2/NOVAC spectra/Results/Figures/';

% Loads results files form linear model 
% AW3
inDirLM_AW3 = '/Users/charlotteb/Documents/Chapter 2/NOVAC spectra/AW3'; % Later add ReRun; % MODEL VERSION 1 IS AW3 cropped frequencies
% 1417 
inFileLM_AW3_1417 = 'ReferenceSpectra_AW3cF_mEst_NOVAC_Masaya_1417.mat'; 
inLM_AW3_1417 = fullfile(inDirLM_AW3, inFileLM_AW3_1417); 
% mEst 
mEstLM_AW3_1417 = load(inLM_AW3_1417,'Results_mEst_real', 'Results_stdErr_real', 'elevAn', 'Results_res_real'); 
mEstLM_AW3_mEst_1417 = mEstLM_AW3_1417.Results_mEst_real;  
mEstLM_AW3_stdErr_1417 = mEstLM_AW3_1417.Results_stdErr_real;  
mEstLM_AW3_angle_1417 = mEstLM_AW3_1417.elevAn; 
Res_AW3_1417 = mEstLM_AW3_1417.Results_res_real;
% 1941
inFileLM_AW3_1941 = 'ReferenceSpectra_AW3cF_mEst_NOVAC_Masaya_1941.mat'; 
inLM_AW3_1941 = fullfile(inDirLM_AW3, inFileLM_AW3_1941); 
% mEst 
mEstLM_AW3_1941 = load(inLM_AW3_1941,'Results_mEst_real', 'Results_stdErr_real', 'elevAn', 'Results_res_real'); 
mEstLM_AW3_mEst_1941 = mEstLM_AW3_1941.Results_mEst_real;  
mEstLM_AW3_stdErr_1941 = mEstLM_AW3_1941.Results_stdErr_real;  
mEstLM_AW3_angle_1941 = mEstLM_AW3_1941.elevAn;  
Res_AW3_1941 = mEstLM_AW3_1941.Results_res_real;

% AW2
inDirLM_AW2 = '/Users/charlotteb/Documents/Chapter 2/NOVAC spectra/AW2'; % Later add ReRun; % MODEL VERSION 1 IS AW2 cropped frequencies
% 1417 
inFileLM_AW2_1417 = 'ReferenceSpectra_AW2cF_mEst_NOVAC_Masaya_1417.mat'; 
inLM_AW2_1417 = fullfile(inDirLM_AW2, inFileLM_AW2_1417); 
% mEst 
mEstLM_AW2_1417 = load(inLM_AW2_1417,'Results_mEst_real', 'Results_stdErr_real', 'elevAn'); 
mEstLM_AW2_mEst_1417 = mEstLM_AW2_1417.Results_mEst_real;  
mEstLM_AW2_stdErr_1417 = mEstLM_AW2_1417.Results_stdErr_real;  
mEstLM_AW2_angle_1417 = mEstLM_AW2_1417.elevAn;  
% 1941
inFileLM_AW2_1941 = 'ReferenceSpectra_AW2cF_mEst_NOVAC_Masaya_1941.mat'; 
inLM_AW2_1941 = fullfile(inDirLM_AW2, inFileLM_AW2_1941); 
% mEst 
mEstLM_AW2_1941 = load(inLM_AW2_1941,'Results_mEst_real', 'Results_stdErr_real', 'elevAn'); 
mEstLM_AW2_mEst_1941 = mEstLM_AW2_1941.Results_mEst_real;  
mEstLM_AW2_stdErr_1941 = mEstLM_AW2_1941.Results_stdErr_real;  
elevAn = mEstLM_AW2_1941.elevAn;  

% Variables for plotting 
markeredgecolor = 'black'; 
markerfacecolor_AW2 = [0.91 0.05 0.05]; % AW2
markerfacecolor_AW3 = [0.09 0.54 0.03]; % AW3
markersize = 4; 
capsize = 0; 
errorbarlinewidthV = 0.5; 
% Include only scans where we have equivalent DOAS result 
elevAn = elevAn(3:end);
mEstLM_AW2_mEst_1417 = mEstLM_AW2_mEst_1417(2,3:end);
mEstLM_AW2_stdErr_1417 = mEstLM_AW2_stdErr_1417(2,3:end);
mEstLM_AW2_mEst_1941 = mEstLM_AW2_mEst_1941(2,3:end);
mEstLM_AW2_stdErr_1941 = mEstLM_AW2_stdErr_1941(2,3:end);
mEstLM_AW3_mEst_1417 = mEstLM_AW3_mEst_1417(2,3:end);
mEstLM_AW3_stdErr_1417 = mEstLM_AW3_stdErr_1417(2,3:end);
mEstLM_AW3_mEst_1941 = mEstLM_AW3_mEst_1941(2,3:end);
mEstLM_AW3_stdErr_1941 = mEstLM_AW3_stdErr_1941(2,3:end);

% Loads DOAS data 
addpath '/Users/charlotteb/Documents/Chapter 2/NOVAC spectra/'
inDirD = '/Users/charlotteb/Documents/Chapter 2/NOVAC spectra/'; 
inFileD_1417 = 'mDOASResults_2022-09-27T11-37-44_D2J2375_140319_1417_0.csv'; 
inFileD_1941 = 'mDOASResults_2022-09-27T11-40-34_D2J2375_140319_1941_0.csv'; 
hlinesD = 1; 
% 1417 
inD = fullfile(inDirD, inFileD_1417); 
% Loads all data in results file
[SpectrumID_1417, FitCoeff_6101_1417,	FitCoeffError_6101_1417,	Shift_6101_1417,	Squeeze_6101_1417, FitCoeff_Ring_1417, FitCoeffError_Ring_1417,	Shift_Ring_1417,...
    Squeeze_Ring_1417, FitCoeff_FLMS195681_O3_Voigt_223K_HP500_PPMM_1417,	FitCoeffError_FLMS195681_O3_Voigt_223K_HP500_PPMM_1417, Shift_FLMS195681_O3_Voigt_223K_HP500_PPMM_1417,...
    Squeeze_FLMS195681_O3_Voigt_223K_HP500_PPMM_1417, FitCoeff_FLMS195681_SO2_Bogumil_293K_HP500_MOLECCM2_1417,...
    FitCoeffError_FLMS195681_SO2_Bogumil_293K_HP500_PPMM_1417, nSO2D_1417, Shift_FLMS195681_SO2_Bogumil_293K_HP500_PPMM_1417,	Squeeze_FLMS195681_SO2_Bogumil_293K_HP500_PPMM_1417,...
    WavelengthRangeLow_1417,	WavelengthRangeHigh_1417, FitChi2_1417, TimeStamp_1417, Latitude_1417, Longitude_1417, Altitude_1417,	Speed_1417, Course_1417, GPSWarnCode_1417,	GPSQuality_1417,	ElevationAngle_1417,	AzimuthAngle_1417,...
    ExposureTime_1417, Exposures_1417, MaxIntensity_1417, MaxIntensityFitRange_1417, FileName_1417, SpectrometerType_1417, SpectrometerSerialNumber_1417,...
    SpectrometerChannel_1417, Remark_1417] = textread(inD,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%s%s%s%s%f%s', 'headerlines', hlinesD, 'delimiter', ',');  
% 1941 
inD = fullfile(inDirD, inFileD_1941); 
[SpectrumID_1941, FitCoeff_6101_1941,	FitCoeffError_6101_1941,	Shift_6101_1941,	Squeeze_6101_1941, FitCoeff_Ring_1941, FitCoeffError_Ring_1941,	Shift_Ring_1941,...
    Squeeze_Ring_1941, FitCoeff_FLMS195681_O3_Voigt_223K_HP500_PPMM_1941,	FitCoeffError_FLMS195681_O3_Voigt_223K_HP500_PPMM_1941, Shift_FLMS195681_O3_Voigt_223K_HP500_PPMM_1941,...
    Squeeze_FLMS195681_O3_Voigt_223K_HP500_PPMM_1941, FitCoeff_FLMS195681_SO2_Bogumil_293K_HP500_MOLECCM2_1941,...
    FitCoeffError_FLMS195681_SO2_Bogumil_293K_HP500_PPMM_1941, nSO2D_1941, Shift_FLMS195681_SO2_Bogumil_293K_HP500_PPMM_1941,	Squeeze_FLMS195681_SO2_Bogumil_293K_HP500_PPMM_1941,...
    WavelengthRangeLow_1941,	WavelengthRangeHigh_1941, FitChi2_1941, TimeStamp_1941, Latitude_1941, Longitude_1941, Altitude_1941,	Speed_1941, Course_1941, GPSWarnCode_1941,	GPSQuality_1941,	ElevationAngle_1941,	AzimuthAngle_1941,...
    ExposureTime_1941, Exposures_1941, MaxIntensity_1941, MaxIntensityFitRange_1941, FileName_1941, SpectrometerType_1941, SpectrometerSerialNumber_1941,...
    SpectrometerChannel_1941, Remark_1941] = textread(inD,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%s%s%s%s%f%s', 'headerlines', hlinesD, 'delimiter', ',');  
 
SO2D_1417 = FitCoeff_FLMS195681_SO2_Bogumil_293K_HP500_MOLECCM2_1417;  
SO2ErrD_1417 = FitCoeffError_FLMS195681_SO2_Bogumil_293K_HP500_PPMM_1417; 
% nSO2D_1417
% ElevationAngle_1417
SO2D_1941 = FitCoeff_FLMS195681_SO2_Bogumil_293K_HP500_MOLECCM2_1941;  
SO2ErrD_1941 = FitCoeffError_FLMS195681_SO2_Bogumil_293K_HP500_PPMM_1941; 
% nSO2D_1941
% ElevationAngle_1941

% Exclude data not within 75 degrees of zenith 
% 5 to 47 
minInd = 5; 
maxInd = 47; 
% Linear model 
elevAn = elevAn(minInd:maxInd);
mEstLM_AW2_mEst_1417 = mEstLM_AW2_mEst_1417(minInd:maxInd);
mEstLM_AW2_mEst_1941 = mEstLM_AW2_mEst_1941(minInd:maxInd);
mEstLM_AW3_mEst_1417 = mEstLM_AW3_mEst_1417(minInd:maxInd);
mEstLM_AW3_mEst_1941 = mEstLM_AW3_mEst_1941(minInd:maxInd);
% error 
mEstLM_AW2_stdErr_1417 = mEstLM_AW2_stdErr_1417(minInd:maxInd);
mEstLM_AW2_stdErr_1941 = mEstLM_AW2_stdErr_1941(minInd:maxInd);
mEstLM_AW3_stdErr_1417 = mEstLM_AW3_stdErr_1417(minInd:maxInd);
mEstLM_AW3_stdErr_1941 = mEstLM_AW3_stdErr_1941(minInd:maxInd);

% DOAS 
ElevationAngle_1417 = ElevationAngle_1417(minInd:maxInd);
SO2D_1417 = SO2D_1417(minInd:maxInd);
SO2ErrD_1417 = SO2ErrD_1417(minInd:maxInd);
nSO2D_1417 = nSO2D_1417(minInd:maxInd);
ElevationAngle_1941 = ElevationAngle_1941(minInd:maxInd);
SO2D_1941 = SO2D_1941(minInd:maxInd);
SO2ErrD_1941 = SO2ErrD_1941(minInd:maxInd);
nSO2D_1941 = nSO2D_1941(minInd:maxInd);

figure
% Combine all raw data 
mEstLM_AW3 = [mEstLM_AW3_mEst_1417'; mEstLM_AW3_mEst_1941']; % Using AW3 
mEstErr_AW3 = [mEstLM_AW3_stdErr_1417'; mEstLM_AW3_stdErr_1941']; % Using AW3 
mEstErr_AW2 = [mEstLM_AW2_stdErr_1417'; mEstLM_AW2_stdErr_1941']; % Using AW2
SO2D = [SO2D_1417; SO2D_1941];
SO2DErr = [SO2ErrD_1417; SO2ErrD_1941];
nSO2D = [nSO2D_1417; nSO2D_1941];
subplot(2,1,1)
% First plot all second scan at standard wavelength 
in1 = 44;
in2 = 86; 
errorbar(SO2D(in1:in2), mEstLM_AW3(in1:in2), mEstErr_AW3(in1:in2), mEstErr_AW3(in1:in2), SO2DErr(in1:in2), SO2DErr(in1:in2), 'LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV, 'CapSize', 0, 'MarkerFaceColor', markerfacecolor_AW3, 'MarkerSize', 6); hold on; 
% Different wavelengths for first scan 
in1 = 1;
in2 = 28; 
errorbar(SO2D(in1:in2), mEstLM_AW3(in1:in2), mEstErr_AW3(in1:in2), mEstErr_AW3(in1:in2), SO2DErr(in1:in2), SO2DErr(in1:in2), 'LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV, 'CapSize', 0, 'MarkerFaceColor', markerfacecolor_AW3, 'MarkerSize', 6); hold on; 
in1 = 30;
in2 = 35; 
errorbar(SO2D(in1:in2), mEstLM_AW3(in1:in2), mEstErr_AW3(in1:in2), mEstErr_AW3(in1:in2), SO2DErr(in1:in2), SO2DErr(in1:in2), 'LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV, 'CapSize', 0, 'MarkerFaceColor', markerfacecolor_AW3, 'MarkerSize', 6); hold on; 
% Plot those at mid-wavelengths 
in1 = 28;
in2 = 30; 
errorbar(SO2D(in1:in2), mEstLM_AW3(in1:in2), mEstErr_AW3(in1:in2), mEstErr_AW3(in1:in2), SO2DErr(in1:in2), SO2DErr(in1:in2), 'LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV, 'CapSize', 0, 'MarkerFaceColor', markerfacecolor_AW3, 'MarkerSize', 8); hold on; 
in1 = 35;
in2 = 37; 
errorbar(SO2D(in1:in2), mEstLM_AW3(in1:in2), mEstErr_AW3(in1:in2), mEstErr_AW3(in1:in2), SO2DErr(in1:in2), SO2DErr(in1:in2), 'LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV, 'CapSize', 0, 'MarkerFaceColor', markerfacecolor_AW3, 'MarkerSize', 8); hold on; 
% Plot those at high wavelengths 
in1 = 38;
in2 = 43;
errorbar(SO2D(in1:in2), mEstLM_AW3(in1:in2), mEstErr_AW3(in1:in2), mEstErr_AW3(in1:in2), SO2DErr(in1:in2), SO2DErr(in1:in2), 'LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV, 'CapSize', 0, 'MarkerFaceColor', markerfacecolor_AW3, 'MarkerSize', 10); hold on; 
%% 
ylabel('mEst_S_O_2 (molec/cm^2)') 
xlabel('DOAS SO_2 fit coefficient (molec/cm^2)') 
xlim([-8e17 2.7e18]);
ylim([-8e17 2.7e18]);
Fig = gca; Fig.FontSize = 12; set(gcf,'color','w'); Fig.Box = 'on';
%% 
subplot(2,1,2)
% Add offset to AW3 
offset_AW3_1417 = 5.56e14;
offset_AW3_1941 = 2.69e16;
offset_AW2_1417 = 4.05e18;
offset_AW2_1941 = 4.08e18;
% 
OFFSETmEstLM_AW3_mEst_1417 = mEstLM_AW3_mEst_1417+offset_AW3_1417;
OFFSETmEstLM_AW3_mEst_1941 = mEstLM_AW3_mEst_1941+offset_AW3_1941;
OFFSETmEstLM_AW2_mEst_1417 = mEstLM_AW2_mEst_1417+offset_AW2_1417;
OFFSETmEstLM_AW2_mEst_1941 = mEstLM_AW2_mEst_1941+offset_AW2_1941;
% 
OFFSETmEstLM_AW3 = [OFFSETmEstLM_AW3_mEst_1417'; OFFSETmEstLM_AW3_mEst_1941']; % Using AW3 
OFFSETmEstLM_AW2 = [OFFSETmEstLM_AW2_mEst_1417'; OFFSETmEstLM_AW2_mEst_1941']; % Using AW3 
% Plot normalised DOAS 
% First plot all at standard wavelength 
% First plot all second scan at standard wavelength 
in1 = 44;
in2 = 86; 
errorbar(nSO2D(in1:in2), mEstLM_AW3(in1:in2), mEstErr_AW3(in1:in2), mEstErr_AW3(in1:in2), SO2DErr(in1:in2), SO2DErr(in1:in2), 'LineStyle', 'none', 'Marker', 's', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV, 'CapSize', 0, 'MarkerFaceColor', markerfacecolor_AW3, 'MarkerSize', 6); hold on; 
% Different wavelengths for first scan 
in1 = 1;
in2 = 28; 
errorbar(nSO2D(in1:in2), OFFSETmEstLM_AW3(in1:in2), mEstErr_AW3(in1:in2), mEstErr_AW3(in1:in2), SO2DErr(in1:in2), SO2DErr(in1:in2), 'LineStyle', 'none', 'Marker', 's', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV, 'CapSize', 0, 'MarkerFaceColor', markerfacecolor_AW3, 'MarkerSize', 6); hold on; 
in1 = 30;
in2 = 35; 
errorbar(nSO2D(in1:in2), OFFSETmEstLM_AW3(in1:in2), mEstErr_AW3(in1:in2), mEstErr_AW3(in1:in2), SO2DErr(in1:in2), SO2DErr(in1:in2), 'LineStyle', 'none', 'Marker', 's', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV, 'CapSize', 0, 'MarkerFaceColor', markerfacecolor_AW3, 'MarkerSize', 6); hold on; 
% Plot those at mid-wavelengths 
in1 = 28;
in2 = 30; 
errorbar(nSO2D(in1:in2), OFFSETmEstLM_AW3(in1:in2), mEstErr_AW3(in1:in2), mEstErr_AW3(in1:in2), SO2DErr(in1:in2), SO2DErr(in1:in2), 'LineStyle', 'none', 'Marker', 's', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV, 'CapSize', 0, 'MarkerFaceColor', markerfacecolor_AW3, 'MarkerSize', 8); hold on; 
in1 = 35;
in2 = 37; 
errorbar(nSO2D(in1:in2), OFFSETmEstLM_AW3(in1:in2), mEstErr_AW3(in1:in2), mEstErr_AW3(in1:in2), SO2DErr(in1:in2), SO2DErr(in1:in2), 'LineStyle', 'none', 'Marker', 's', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV, 'CapSize', 0, 'MarkerFaceColor', markerfacecolor_AW3, 'MarkerSize', 8); hold on; 
% Then we have two scan at 57 and 61 degrees which is higher DOAS wavelengths but still AW3
in1 = 38;
in2 = 39;
errorbar(nSO2D(in1:in2), OFFSETmEstLM_AW3(in1:in2), mEstErr_AW3(in1:in2), mEstErr_AW3(in1:in2), SO2DErr(in1:in2), SO2DErr(in1:in2), 'LineStyle', 'none', 'Marker', 's', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV, 'CapSize', 0, 'MarkerFaceColor', markerfacecolor_AW3, 'MarkerSize', 10); hold on; 
% Rest are higher wavelengths and AW3 
in1 = 40;
in2 = 43;
errorbar(nSO2D(in1:in2), OFFSETmEstLM_AW2(in1:in2), mEstErr_AW2(in1:in2), mEstErr_AW2(in1:in2), SO2DErr(in1:in2), SO2DErr(in1:in2), 'LineStyle', 'none', 'Marker', 's', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV, 'CapSize', 0, 'MarkerFaceColor', markerfacecolor_AW2, 'MarkerSize', 10); hold on; 
%% 
ylabel('Offset mEst_S_O_2 (molec/cm^2)') 
xlabel('Normalised DOAS SO_2 fit coefficient (molec/cm^2)') 
xlim([0 2.7e18]);
ylim([0 2.7e18]);
Fig = gca; Fig.FontSize = 12; set(gcf,'color','w'); Fig.Box = 'on';
%% Make any adjustments 
% SAVES FIGURE
fnOut = 'Fig_final_27';
fname = fullfile(outDirFig, fnOut);

