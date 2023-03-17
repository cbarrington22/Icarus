
% Plotting results of reference spectra 

addpath '/Users/charlotteb/Documents/Chapter 2/Reference spectra /DOAS_Eval_CK/'
addpath '/Users/charlotteb/Documents/Chapter 2/Reference spectra /Output/'
outDir = '/Users/charlotteb/Documents/Chapter 2/Reference spectra /outFiles/';
outDirFig = '/Users/charlotteb/Documents/Chapter 2/Reference spectra /';

% Loads DOAS results from CK
inDirD = '/Users/charlotteb/Documents/Chapter 2/Reference spectra /DOAS_Eval_CK/';
inFileD = 'mDOASResults_2022-02-02T11-47-08_310-340nm_VARIABLE.csv'; 
inD = fullfile(inDirD, inFileD); 

hlinesD = 1; 

% Loads all data in results file
[SpectrumID, FitCoeff_6101,	FitCoeffError_6101,	Shift_6101,	Squeeze_6101, FitCoeff_Ring, FitCoeffError_Ring,	Shift_Ring,...
    Squeeze_Ring, FitCoeff_FLMS195681_O3_Voigt_223K_HP500_PPMM,	FitCoeffError_FLMS195681_O3_Voigt_223K_HP500_PPMM, Shift_FLMS195681_O3_Voigt_223K_HP500_PPMM,...
    Squeeze_FLMS195681_O3_Voigt_223K_HP500_PPMM, FitCoeff_FLMS195681_SO2_Bogumil_293K_HP500_MOLECCM2, FitCoeff_FLMS195681_SO2_Bogumil_293K_HP500_PPMM,...
    FitCoeffError_FLMS195681_SO2_Bogumil_293K_HP500_PPMM, SO2_Error_PPMM, Shift_FLMS195681_SO2_Bogumil_293K_HP500_PPMM,	Squeeze_FLMS195681_SO2_Bogumil_293K_HP500_PPMM,...
    WavelengthRangeLow,	WavelengthRangeHigh, FitChi2, TimeStamp, Latitude, Longitude, Altitude,	Speed, Course, GPSWarnCode,	GPSQuality,	ElevationAngle,	AzimuthAngle,...
    ExposureTime, Exposures, MaxIntensity, MaxIntensityFitRange, FileName, SpectrometerType, SpectrometerSerialNumber,...
    SpectrometerChannel, Remark] = textread(inD,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%s%s%s%s%f%s', 'headerlines', hlinesD, 'delimiter', ',');  
 
% Isolates SO2 FIT COEFFICIENT and FIT ERROR 
SO2D = FitCoeff_FLMS195681_SO2_Bogumil_293K_HP500_MOLECCM2;  
SO2ErrD = FitCoeffError_FLMS195681_SO2_Bogumil_293K_HP500_PPMM; 

% Extracts known SO2 concentration from 'remark'
findStr = "_FLMS195681";
tempStr = string(extractBefore(Remark, findStr)); 
findStr = "_";
elevAnD = str2double(extractBefore(tempStr, findStr)); % Elevation angle of scan 
cellConD = str2double(extractAfter(tempStr, findStr)); % Gas cell concentration 

% Gas cell error is +/- 10%
% Conversion from ppmm to molec/cm^2 for expected trend line
conv = 2500000000000000; % copied form DOAS results table 
cellConDconv = cellConD * conv; 
cellConDErr = (cellConD/100) * 10; 
cellConDErrconv = cellConDErr * conv;

% Loads results files form linear model 
% MODEL VERSION 1 
inDirLM_v1 = '/Users/charlotteb/Documents/Chapter 2/Reference spectra /Output/Model version 1 AW1cF/ReRun'; % Later add ReRun; % MODEL VERSION 1 IS AW1 cropped frequencies
inFileLM_v1 = 'ReferenceSpectra_AW1cF_mEst_manualsaveall.mat'; 
inLM_v1 = fullfile(inDirLM_v1, inFileLM_v1); 
% mEst 
mEstLM_v1_real = load(inLM_v1,'Results_mEst_real', 'cellCon'); 
cellConLM = mEstLM_v1_real.cellCon;  
mEstLM_v1_real = mEstLM_v1_real.Results_mEst_real; 
mEstLM_v1_realimag = load(inLM_v1,'Results_mEst_realimag'); 
mEstLM_v1_realimag = mEstLM_v1_realimag.Results_mEst_realimag; 
% Standard error 
stdErrLM_v1_real = load(inLM_v1,'Results_stdErr_real'); 
stdErrLM_v1_real = stdErrLM_v1_real.Results_stdErr_real; 
stdErrLM_v1_realimag = load(inLM_v1,'Results_stdErr_realimag'); 
stdErrLM_v1_realimag = stdErrLM_v1_realimag.Results_stdErr_realimag; 

% MODEL VERSION 2
inDirLM_v2 = '/Users/charlotteb/Documents/Chapter 2/Reference spectra /Output/Model version 4 AW2cF bonus high SO2/'; % MODEL VERSION 2 IS AW2 (bonus)
inFileLM_v2 = 'ReferenceSpectra_AW2cF_mEst_manualsaveall.mat'; 
inLM_v2 = fullfile(inDirLM_v2, inFileLM_v2); 
% mEst 
mEstLM_v2_real = load(inLM_v2,'Results_mEst_real'); 
mEstLM_v2_real = mEstLM_v2_real.Results_mEst_real; 
mEstLM_v2_realimag = load(inLM_v2,'Results_mEst_realimag'); 
mEstLM_v2_realimag = mEstLM_v2_realimag.Results_mEst_realimag; 
% Standard error 
stdErrLM_v2_real = load(inLM_v2,'Results_stdErr_real'); 
stdErrLM_v2_real = stdErrLM_v2_real.Results_stdErr_real; 
stdErrLM_v2_realimag = load(inLM_v2,'Results_stdErr_realimag'); 
stdErrLM_v2_realimag = stdErrLM_v2_realimag.Results_stdErr_realimag; 

% MODEL VERSION 3
inDirLM_v3 = '/Users/charlotteb/Documents/Chapter 2/Reference spectra /Output/Model version 2 AW3cF/ReRun/'; % MODEL VERSION 3 IS AW3 
inFileLM_v3 = 'ReferenceSpectra_AW3cF_mEst.mat'; 
inLM_v3 = fullfile(inDirLM_v3, inFileLM_v3); 
% mEst 
mEstLM_v3_real = load(inLM_v3,'Results_mEst_real'); 
mEstLM_v3_real = mEstLM_v3_real.Results_mEst_real; 
mEstLM_v3_realimag = load(inLM_v3,'Results_mEst_realimag'); 
mEstLM_v3_realimag = mEstLM_v3_realimag.Results_mEst_realimag; 
% Standard error 
stdErrLM_v3_real = load(inLM_v3,'Results_stdErr_real'); 
stdErrLM_v3_real = stdErrLM_v3_real.Results_stdErr_real; 
stdErrLM_v3_realimag = load(inLM_v3,'Results_stdErr_realimag'); 
stdErrLM_v3_realimag = stdErrLM_v3_realimag.Results_stdErr_realimag; 

% MODEL VERSION 4 
inDirLM_v4 = '/Users/charlotteb/Documents/Chapter 2/Reference spectra /Output/Model version 3 AW1aF/'; % MODEL VERSION 4 IS AW1 all frequencies 
inFileLM_v4 = 'ReferenceSpectra_AW1aF_mEst_manualsaveall.mat'; 
inLM_v4 = fullfile(inDirLM_v4, inFileLM_v4); 
% mEst 
mEstLM_v4_real = load(inLM_v4,'Results_mEst_real'); 
mEstLM_v4_real = mEstLM_v4_real.Results_mEst_real; 
mEstLM_v4_realimag = load(inLM_v4,'Results_mEst_realimag'); 
mEstLM_v4_realimag = mEstLM_v4_realimag.Results_mEst_realimag; 
% Standard error 
stdErrLM_v4_real = load(inLM_v4,'Results_stdErr_real'); 
stdErrLM_v4_real = stdErrLM_v4_real.Results_stdErr_real; 
stdErrLM_v4_realimag = load(inLM_v4,'Results_stdErr_realimag'); 
stdErrLM_v4_realimag = stdErrLM_v4_realimag.Results_stdErr_realimag; 

% Error calculation using linear model data 
% Gas cell error is +/- 10%
% Conversion from ppmm to molec/cm^2 for expected trend line
cellConLMconv = cellConLM * conv; 
cellConLMErr = (cellConLM/100) * 10; 
cellConLMErrconv = cellConLMErr * conv;
% Adds bounds to expected trend line to relating to gas cell error  
maxExp = cellConLMconv + cellConLMErrconv; 
minExp = cellConLMconv - cellConLMErrconv; 

% Variables for formatting plot 
markeredgecolor = 'black'; 
markerfacecolorD =  'white';
markerfacecolor_v1 = [0.11 0.24 0.81]; % AW1 
markerfacecolor_v1cm = [0.62 0.68 0.95];
markerfacecolor_v4 = [0.11 0.24 0.81]; % AW1
markerfacecolor_v2 = [0.91 0.05 0.05]; % AW2
markerfacecolor_v3 = [0.09 0.54 0.03]; % AW3
markersize = 6; 
capsize = 0; 
errorbarcolor = [0.92,0.92,0.92]; % Horizontal 
errorbarlinewidthH = 0.1; 
errorbarlinewidthV = 0.5; 
minX = -70; 
maxX = 2330; 
trendlinewidthmain = 0.1; 
% trendlinewidthlimits = 0.1; 
% Axis limits 
minY = min(cellConD)-1e18;
maxY = 8.2e18; % max(cellConDconv)+3.5e18;

% Plots SO2 fit results from DOAS analysis 
figure
subplot(2,2,1) % MODEL VERSION 1 
% Gas cell error bars 
e = errorbar(cellConD, SO2D, cellConDErr, cellConDErr, 'horizontal', 'LineStyle', 'none', 'Marker', 'none', 'Color', errorbarcolor, 'linewidth', errorbarlinewidthH); hold on; 
e.CapSize = capsize; 
% Expected 
et = plot(cellConLM, cellConLMconv, ':k'); et.LineWidth = trendlinewidthmain; hold on; 
SO2LMr = mEstLM_v1_real(2,:); % SO2 is G2 (I0, SO2, O3, Ring, Bshift) 
SO2LMri = mEstLM_v1_realimag(2,:); % SO2 is G2 (I0, SO2, O3, Ring, Bshift) 
SO2ErrLMr = stdErrLM_v1_real(2,:); 
SO2ErrLMri = stdErrLM_v1_realimag(2,:); 
% Plots mEst 
e = errorbar(cellConLM, SO2LMr, cellConLMErr, cellConLMErr, 'horizontal', 'LineStyle', 'none', 'Marker', 'none', 'Color', errorbarcolor, 'linewidth', errorbarlinewidthV); hold on; 
e.CapSize = capsize; 
e = errorbar(cellConLM, SO2LMri, cellConLMErr, cellConLMErr, 'horizontal', 'LineStyle', 'none', 'Marker', 'none', 'Color', errorbarcolor, 'linewidth', errorbarlinewidthV); hold on; 
e.CapSize = capsize; 
% etM = plot(cellConLM, maxExp, '-k'); etM.LineWidth = trendlinewidthlimits;
% etm = plot(cellConLM, minExp, '-k'); etm.LineWidth = trendlinewidthlimits;
% DOAS 
d = errorbar(cellConD, SO2D, SO2ErrD, SO2ErrD, 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);
d.MarkerFaceColor = markerfacecolorD; 
d.MarkerEdgeColor = markeredgecolor; 
d.CapSize = capsize; 
xlim([minX maxX]);
ylim([minY maxY]);
% ylabel('(d)SCD SO_2 (molec/cm^2)')
% xlabel('Gas cell concentration (ppm·m)')
Fig = gca; Fig.FontSize = 14; set(gcf,'color','w'); Fig.Box = 'on'; 
% mEst and standard error 
% Real and complex
% SO2LMr = mEstLM_v1_real(2,:); % SO2 is G2 (I0, SO2, O3, Ring, Bshift) 
% SO2LMri = mEstLM_v1_realimag(2,:); % SO2 is G2 (I0, SO2, O3, Ring, Bshift) 
% SO2ErrLMr = stdErrLM_v1_real(2,:); 
% SO2ErrLMri = stdErrLM_v1_realimag(2,:); 
% % Plots mEst 
% e = errorbar(cellConLM, SO2LMr, cellConLMErr, cellConLMErr, 'horizontal', 'LineStyle', 'none', 'Marker', 'none', 'Color', errorbarcolor, 'linewidth', errorbarlinewidthV); hold on; 
% e.CapSize = capsize; 
% e = errorbar(cellConLM, SO2LMri, cellConLMErr, cellConLMErr, 'horizontal', 'LineStyle', 'none', 'Marker', 'none', 'Color', errorbarcolor, 'linewidth', errorbarlinewidthV); hold on; 
% e.CapSize = capsize; 
% Real
lr = errorbar(cellConLM, SO2LMr, SO2ErrLMr, SO2ErrLMr, 'vertical', 'LineStyle', 'none', 'Marker', '^', 'Color', markerfacecolor_v1, 'linewidth', errorbarlinewidthV);
lr.MarkerFaceColor = markerfacecolor_v1; 
lr.MarkerEdgeColor = markeredgecolor; 
lr.CapSize = capsize; 
% Real imaginary 
lri = errorbar(cellConLM, SO2LMri, SO2ErrLMri, SO2ErrLMri, 'vertical', 'LineStyle', 'none', 'Marker', 'v', 'Color', markerfacecolor_v1cm, 'linewidth', errorbarlinewidthV);
lri.MarkerFaceColor = markerfacecolor_v1cm; 
lri.MarkerEdgeColor = markeredgecolor;  % markeredgecolor; 
lri.CapSize = capsize;
% Offset 
nogascellInd = find(cellConLM == 0); 
nogascell = cellConLM(nogascellInd); 
nogascellmEstr = SO2LMr(nogascellInd); 
nogascellmEstri = SO2LMri(nogascellInd); 
% figure 
% plot(nogascell, nogascellmEstr', 'ok') 
disp("Version 1")
disp("check negative")
max(nogascellmEstr) 
max(nogascellmEstri) 
disp("find mean offset")
meanNogascellmEstr = mean(nogascellmEstr)
meanNogascellmEstri = mean(nogascellmEstri)

% MODEL VERSION 2
subplot(2,2,2) 
% Gas cell error bars 
e = errorbar(cellConD, SO2D, cellConDErr, cellConDErr, 'horizontal', 'LineStyle', 'none', 'Marker', 'none', 'Color', errorbarcolor, 'linewidth', errorbarlinewidthH); hold on; 
e.CapSize = capsize; 
% Expected 
et = plot(cellConLM, cellConLMconv, ':k'); et.LineWidth = trendlinewidthmain; hold on; 
% etM = plot(cellConLM, maxExp, '-k'); etM.LineWidth = trendlinewidthlimits;
% etm = plot(cellConLM, minExp, '-k'); etm.LineWidth = trendlinewidthlimits;
SO2LMr = mEstLM_v2_real(2,:); % SO2 is G2 (I0, SO2, O3, Ring, Bshift) 
% SO2LMri = mEstLM_v2_realimag(2,:); % SO2 is G2 (I0, SO2, O3, Ring, Bshift) 
SO2ErrLMr = stdErrLM_v2_real(2,:); 
% SO2ErrLMri = stdErrLM_v2_realimag(2,:); 
% SO2ErrLMr = NaN(size(SO2LMr)); % Placeholder
% SO2ErrLMri = NaN(size(SO2LMri)); % Placeholder
% Plots mEst 
e = errorbar(cellConLM, SO2LMr, cellConLMErr, cellConLMErr, 'horizontal', 'LineStyle', 'none', 'Marker', 'none', 'Color', errorbarcolor, 'linewidth', errorbarlinewidthV); hold on; 
e.CapSize = capsize; 
% DOAS 
d = errorbar(cellConD, SO2D, SO2ErrD, SO2ErrD, 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);
d.MarkerFaceColor = markerfacecolorD; 
d.MarkerEdgeColor = markeredgecolor; 
d.CapSize = capsize; 
xlim([minX maxX]);
ylim([minY maxY]);
% ylabel('(d)SCD SO_2 (molec/cm^2)')
% xlabel('Gas cell concentration (ppm·m)')
Fig = gca; Fig.FontSize = 14; set(gcf,'color','w'); Fig.Box = 'on'; 
% mEst and standard error 
% Real and complex
% SO2LMr = mEstLM_v2_real(2,:); % SO2 is G2 (I0, SO2, O3, Ring, Bshift) 
% % SO2LMri = mEstLM_v2_realimag(2,:); % SO2 is G2 (I0, SO2, O3, Ring, Bshift) 
% SO2ErrLMr = stdErrLM_v2_real(2,:); 
% % SO2ErrLMri = stdErrLM_v2_realimag(2,:); 
% % SO2ErrLMr = NaN(size(SO2LMr)); % Placeholder
% % SO2ErrLMri = NaN(size(SO2LMri)); % Placeholder
% % Plots mEst 
% e = errorbar(cellConLM, SO2LMr, cellConLMErr, cellConLMErr, 'horizontal', 'LineStyle', 'none', 'Marker', 'none', 'Color', errorbarcolor, 'linewidth', errorbarlinewidthV); hold on; 
% e.CapSize = capsize; 
% e = errorbar(cellConLM, SO2LMri, cellConLMErr, cellConLMErr, 'horizontal', 'LineStyle', 'none', 'Marker', 'none', 'Color', errorbarcolor, 'linewidth', errorbarlinewidthV); hold on; 
% e.CapSize = capsize; 
% Real
lr = errorbar(cellConLM, SO2LMr, SO2ErrLMr, SO2ErrLMr, 'vertical', 'LineStyle', 'none', 'Marker', '^', 'Color', markerfacecolor_v2, 'linewidth', errorbarlinewidthV);
lr.MarkerFaceColor = markerfacecolor_v2; 
lr.MarkerEdgeColor = markeredgecolor;  % markeredgecolor; 
lr.CapSize = capsize; 
% % Real imaginary 
% lri = errorbar(cellConLM, SO2LMri, SO2ErrLMri, SO2ErrLMri, 'vertical', 'LineStyle', 'none', 'Marker', 'v', 'Color', markerfacecolor_v2, 'linewidth', errorbarlinewidthV);
% lri.MarkerFaceColor = markerfacecolorLMAW3; 
% lri.MarkerEdgeColor = 'none';  % markeredgecolor; 
% lri.CapSize = capsize;
% Offset 
nogascellInd = find(cellConLM == 0); 
nogascell = cellConLM(nogascellInd); 
nogascellmEstr = SO2LMr(nogascellInd); 
% nogascellmEstri = SO2LMri(nogascellInd); 
% figure 
% plot(nogascell, nogascellmEstr', 'ok') 
disp("Version 2")
disp("check negative")
max(nogascellmEstr) 
% max(nogascellmEstri) 
disp("find mean offset")
meanNogascellmEstr = mean(nogascellmEstr)
% meanNogascellmEstri = mean(nogascellmEstri)

% MODEL VERSION 3
subplot(2,2,3) 
% Gas cell error bars 
e = errorbar(cellConD, SO2D, cellConDErr, cellConDErr, 'horizontal', 'LineStyle', 'none', 'Marker', 'none', 'Color', errorbarcolor, 'linewidth', errorbarlinewidthH); hold on; 
e.CapSize = capsize; 
% Expected 
et = plot(cellConLM, cellConLMconv, ':k'); et.LineWidth = trendlinewidthmain; hold on; 
% etM = plot(cellConLM, maxExp, '-k'); etM.LineWidth = trendlinewidthlimits;
% etm = plot(cellConLM, minExp, '-k'); etm.LineWidth = trendlinewidthlimits;
% mEst and standard error 
% Real and complex
SO2LMr = mEstLM_v3_real(2,:); % SO2 is G2 (I0, SO2, O3, Ring, Bshift) 
% SO2LMri = mEstLM_v3_realimag(2,:); % SO2 is G2 (I0, SO2, O3, Ring, Bshift) 
SO2ErrLMr = stdErrLM_v3_real(2,:); 
% SO2ErrLMri = stdErrLM_v3_realimag(2,:); 
% Plots mEst 
e = errorbar(cellConLM, SO2LMr, cellConLMErr, cellConLMErr, 'horizontal', 'LineStyle', 'none', 'Marker', 'none', 'Color', errorbarcolor, 'linewidth', errorbarlinewidthV); hold on; 
e.CapSize = capsize; 
% DOAS 
d = errorbar(cellConD, SO2D, SO2ErrD, SO2ErrD, 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);
d.MarkerFaceColor = markerfacecolorD; 
d.MarkerEdgeColor = markeredgecolor; 
d.CapSize = capsize; 
xlim([minX maxX]);
ylim([minY maxY]);
ylabel('(d)SCD SO_2 (molec/cm^2)')
% xlabel('Gas cell concentration (ppm·m)')
Fig = gca; Fig.FontSize = 14; set(gcf,'color','w'); Fig.Box = 'on';
% % mEst and standard error 
% % Real and complex
% SO2LMr = mEstLM_v3_real(2,:); % SO2 is G2 (I0, SO2, O3, Ring, Bshift) 
% % SO2LMri = mEstLM_v3_realimag(2,:); % SO2 is G2 (I0, SO2, O3, Ring, Bshift) 
% SO2ErrLMr = stdErrLM_v3_real(2,:); 
% % SO2ErrLMri = stdErrLM_v3_realimag(2,:); 
% % Plots mEst 
% e = errorbar(cellConLM, SO2LMr, cellConLMErr, cellConLMErr, 'horizontal', 'LineStyle', 'none', 'Marker', 'none', 'Color', errorbarcolor, 'linewidth', errorbarlinewidthV); hold on; 
% e.CapSize = capsize; 
% Real
lr = errorbar(cellConLM, SO2LMr, SO2ErrLMr, SO2ErrLMr, 'vertical', 'LineStyle', 'none', 'Marker', '^', 'Color', markerfacecolor_v3, 'linewidth', errorbarlinewidthV);
lr.MarkerFaceColor = markerfacecolor_v3; 
lr.MarkerEdgeColor = markeredgecolor; 
lr.CapSize = capsize; 
% Real imaginary 
% lri = errorbar(cellConLM, SO2LMri, SO2ErrLMri, SO2ErrLMri, 'vertical', 'LineStyle', 'none', 'Marker', 'v', 'Color', markerfacecolor_v3, 'linewidth', errorbarlinewidthV);
% lri.MarkerFaceColor = 'none'; 
% lri.MarkerEdgeColor = markerfacecolorLMAW1; 
% lri.CapSize = capsize;
% Offset 
nogascellInd = find(cellConLM == 0); 
nogascell = cellConLM(nogascellInd); 
nogascellmEstr = SO2LMr(nogascellInd); 
% nogascellmEstri = SO2LMri(nogascellInd); 
% figure 
% plot(nogascell, nogascellmEstr', 'ok') 
disp("Version 3")
disp("check negative")
max(nogascellmEstr) 
% max(nogascellmEstri) 
disp("find mean offset")
meanNogascellmEstr = mean(nogascellmEstr)
% meanNogascellmEstri = mean(nogascellmEstri)

% MODEL VERSION 4
subplot(2,2,4) 
% Gas cell error bars 
e = errorbar(cellConD, SO2D, cellConDErr, cellConDErr, 'horizontal', 'LineStyle', 'none', 'Marker', 'none', 'Color', errorbarcolor, 'linewidth', errorbarlinewidthH); hold on; 
e.CapSize = capsize; 
% Expected 
et = plot(cellConLM, cellConLMconv, ':k'); et.LineWidth = trendlinewidthmain; hold on; 
% etM = plot(cellConLM, maxExp, '-k'); etM.LineWidth = trendlinewidthlimits;
% etm = plot(cellConLM, minExp, '-k'); etm.LineWidth = trendlinewidthlimits;
% mEst and standard error 
% Real and complex
SO2LMr = mEstLM_v4_real(2,:); % SO2 is G2 (I0, SO2, O3, Ring, Bshift) 
% SO2LMri = mEstLM_v4_realimag(2,:); % SO2 is G2 (I0, SO2, O3, Ring, Bshift) 
SO2ErrLMr = stdErrLM_v4_real(2,:); 
% SO2ErrLMri = stdErrLM_v4_realimag(2,:); 
% Plots mEst 
e = errorbar(cellConLM, SO2LMr, cellConLMErr, cellConLMErr, 'horizontal', 'LineStyle', 'none', 'Marker', 'none', 'Color', errorbarcolor, 'linewidth', errorbarlinewidthV); hold on; 
e.CapSize = capsize; 
e = errorbar(cellConLM, SO2LMri, cellConLMErr, cellConLMErr, 'horizontal', 'LineStyle', 'none', 'Marker', 'none', 'Color', errorbarcolor, 'linewidth', errorbarlinewidthV); hold on; 
e.CapSize = capsize; 
% DOAS 
d = errorbar(cellConD, SO2D, SO2ErrD, SO2ErrD, 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);
d.MarkerFaceColor = markerfacecolorD; 
d.MarkerEdgeColor = markeredgecolor; 
d.CapSize = capsize; 
xlim([minX maxX]);
ylim([minY maxY]);
% ylabel('(d)SCD SO_2 (molec/cm^2)')
xlabel('Gas cell concentration (ppm·m)')
Fig = gca; Fig.FontSize = 14; set(gcf,'color','w'); Fig.Box = 'on';
% % mEst and standard error 
% % Real and complex
% SO2LMr = mEstLM_v4_real(2,:); % SO2 is G2 (I0, SO2, O3, Ring, Bshift) 
% % SO2LMri = mEstLM_v4_realimag(2,:); % SO2 is G2 (I0, SO2, O3, Ring, Bshift) 
% SO2ErrLMr = stdErrLM_v4_real(2,:); 
% % SO2ErrLMri = stdErrLM_v4_realimag(2,:); 
% % Plots mEst 
% e = errorbar(cellConLM, SO2LMr, cellConLMErr, cellConLMErr, 'horizontal', 'LineStyle', 'none', 'Marker', 'none', 'Color', errorbarcolor, 'linewidth', errorbarlinewidthV); hold on; 
% e.CapSize = capsize; 
% e = errorbar(cellConLM, SO2LMri, cellConLMErr, cellConLMErr, 'horizontal', 'LineStyle', 'none', 'Marker', 'none', 'Color', errorbarcolor, 'linewidth', errorbarlinewidthV); hold on; 
% e.CapSize = capsize; 
% Real
lr = errorbar(cellConLM, SO2LMr, SO2ErrLMr, SO2ErrLMr, 'vertical', 'LineStyle', 'none', 'Marker', '^', 'Color', markerfacecolor_v1, 'linewidth', errorbarlinewidthV);
lr.MarkerFaceColor = markerfacecolor_v1; 
lr.MarkerEdgeColor = markeredgecolor; 
lr.CapSize = capsize; 
% Real imaginary 
lri = errorbar(cellConLM, SO2LMri, SO2ErrLMri, SO2ErrLMri, 'vertical', 'LineStyle', 'none', 'Marker', 'v', 'Color', markerfacecolor_v1cm, 'linewidth', errorbarlinewidthV);
lri.MarkerFaceColor = markerfacecolor_v1cm; 
lri.MarkerEdgeColor = markeredgecolor;  % markeredgecolor; 
lri.CapSize = capsize;
% Offset 
nogascellInd = find(cellConLM == 0); 
nogascell = cellConLM(nogascellInd); 
nogascellmEstr = SO2LMr(nogascellInd); 
nogascellmEstri = SO2LMri(nogascellInd); 
% figure 
% plot(nogascell, nogascellmEstr', 'ok') 
disp("Version 4")
disp("check negative")
max(nogascellmEstr) 
max(nogascellmEstri) 
disp("find mean offset")
meanNogascellmEstr = mean(nogascellmEstr)
meanNogascellmEstri = mean(nogascellmEstri)
% Save 
fname = fullfile(outDir, 'results');
save(fname); % Saves all variables  
%% Make any adjustments 
% SAVES FIGURE
% fnOut = 'Results';
% fname = fullfile(outDirFig, fnOut);
% saveas(gca, fname)

% % Check complex magnitude of AW2 
% % MODEL VERSION 2
% figure
% % Gas cell error bars 
% e = errorbar(cellConD, SO2D, cellConDErr, cellConDErr, 'horizontal', 'LineStyle', 'none', 'Marker', 'none', 'Color', errorbarcolor, 'linewidth', errorbarlinewidthH); hold on; 
% e.CapSize = capsize; 
% % Expected 
% et = plot(cellConLM, cellConLMconv, ':k'); et.LineWidth = trendlinewidthmain; hold on; 
% % etM = plot(cellConLM, maxExp, '-k'); etM.LineWidth = trendlinewidthlimits;
% % etm = plot(cellConLM, minExp, '-k'); etm.LineWidth = trendlinewidthlimits;
% SO2LMr = mEstLM_v2_real(2,:); % SO2 is G2 (I0, SO2, O3, Ring, Bshift) 
% SO2LMri = mEstLM_v2_realimag(2,:); % SO2 is G2 (I0, SO2, O3, Ring, Bshift) 
% SO2ErrLMr = stdErrLM_v2_real(2,:); 
% SO2ErrLMri = stdErrLM_v2_realimag(2,:); 
% % SO2ErrLMr = NaN(size(SO2LMr)); % Placeholder
% % SO2ErrLMri = NaN(size(SO2LMri)); % Placeholder
% % Plots mEst 
% e = errorbar(cellConLM, SO2LMr, cellConLMErr, cellConLMErr, 'horizontal', 'LineStyle', 'none', 'Marker', 'none', 'Color', errorbarcolor, 'linewidth', errorbarlinewidthV); hold on; 
% e.CapSize = capsize; 
% % DOAS 
% d = errorbar(cellConD, SO2D, SO2ErrD, SO2ErrD, 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);
% d.MarkerFaceColor = markerfacecolorD; 
% d.MarkerEdgeColor = markeredgecolor; 
% d.CapSize = capsize; 
% xlim([minX maxX]);
% ylim([minY maxY]);
% % ylabel('(d)SCD SO_2 (molec/cm^2)')
% % xlabel('Gas cell concentration (ppm·m)')
% Fig = gca; Fig.FontSize = 14; set(gcf,'color','w'); Fig.Box = 'on'; 
% % mEst and standard error 
% % Real and complex
% % SO2LMr = mEstLM_v2_real(2,:); % SO2 is G2 (I0, SO2, O3, Ring, Bshift) 
% % % SO2LMri = mEstLM_v2_realimag(2,:); % SO2 is G2 (I0, SO2, O3, Ring, Bshift) 
% % SO2ErrLMr = stdErrLM_v2_real(2,:); 
% % % SO2ErrLMri = stdErrLM_v2_realimag(2,:); 
% % % SO2ErrLMr = NaN(size(SO2LMr)); % Placeholder
% % % SO2ErrLMri = NaN(size(SO2LMri)); % Placeholder
% % % Plots mEst 
% % e = errorbar(cellConLM, SO2LMr, cellConLMErr, cellConLMErr, 'horizontal', 'LineStyle', 'none', 'Marker', 'none', 'Color', errorbarcolor, 'linewidth', errorbarlinewidthV); hold on; 
% % e.CapSize = capsize; 
% % e = errorbar(cellConLM, SO2LMri, cellConLMErr, cellConLMErr, 'horizontal', 'LineStyle', 'none', 'Marker', 'none', 'Color', errorbarcolor, 'linewidth', errorbarlinewidthV); hold on; 
% % e.CapSize = capsize; 
% % Real
% % lr = errorbar(cellConLM, SO2LMr, SO2ErrLMr, SO2ErrLMr, 'vertical', 'LineStyle', 'none', 'Marker', '^', 'Color', markerfacecolor_v2, 'linewidth', errorbarlinewidthV);
% % lr.MarkerFaceColor = markerfacecolor_v2; 
% % lr.MarkerEdgeColor = markeredgecolor;  % markeredgecolor; 
% % lr.CapSize = capsize; 
% % % Real imaginary 
% lri = errorbar(cellConLM, SO2LMri, SO2ErrLMri, SO2ErrLMri, 'vertical', 'LineStyle', 'none', 'Marker', 'v', 'Color', markerfacecolor_v2, 'linewidth', errorbarlinewidthV);
% lri.MarkerFaceColor = markerfacecolorLMAW3; 
% lri.MarkerEdgeColor = 'none';  % markeredgecolor; 
% lri.CapSize = capsize;

% %% Figure 19 - mEst by elevation angle 
% offset = 2.7e17; % Table 3 
% % mEst and standard error 
% % Real and complex
% SO2LMr = mEstLM_v3_real(2,:); % SO2 is G2 (I0, SO2, O3, Ring, Bshift) % MODEL VERSION 3
% % SO2LMri = mEstLM_v3_realimag(2,:); % SO2 is G2 (I0, SO2, O3, Ring, Bshift) 
% SO2ErrLMr = stdErrLM_v3_real(2,:); 
% % SO2ErrLMri = stdErrLM_v3_realimag(2,:); 
% SO2LMrOffset = SO2LMr+offset;  
% %% Plot by elevation angle of scan 
% unAngle = unique(elevAnD); % Angles taken from spectra file name and refer to angle from horizon. So 90 is zenith and 30 is close to horizon. 
% unAngle = flip(unAngle); % Reverse so zenith is first 
% % measpAngle = 1000; 
% % doasSO2 = nan(measpAngle, length(unAngle)); 
% % doasSO2Err = nan(measpAngle, length(unAngle)); 
% % mEstSO2 = nan(measpAngle, length(unAngle)); 
% % mEstSO2Err = nan(measpAngle, length(unAngle)); 
% % Organsise data by elevation angle 
% figure
% markersize = 4; 
% markerfacecolor_v3 = [0.09 0.54 0.03]; % AW3
% errorbarlinewidthV = 0.5; 
% for i = 1:length(unAngle)
%     subplot(4,2,i) 
%     ind = find(elevAnD == unAngle(i)); 
%     % Gas cell error bars
%     % Turn colour into variable 
%     % Zenith is markerfacecolor_v3 and goe slighter accoridng to increased elevation angle 
%     % col = i/7; 
%     % markerfacecolor_v3 = [0.09 0.54 0.03]; % AW3
%     if i == 1 
%         markerfacecolor_v3 = [0.09 0.32 0.07]; % AW3
%     elseif i == 2
%         markerfacecolor_v3 = [0.16 0.49 0.13]; 
%          elseif i == 3
%                 markerfacecolor_v3 = [0.18 0.52 0.15]; 
%               elseif i == 4
%                   markerfacecolor_v3 = [0.19 0.60 0.17]; 
%                    elseif i == 5
%                        markerfacecolor_v3 = [0.39 0.81 0.35]; 
%                     elseif i == 6 
%                             markerfacecolor_v3 = [0.61 0.83 0.59]; 
%                        elseif i == 7 
%                                 markerfacecolor_v3 = [0.70,0.92,0.67]; 
%     end
%     % For DOAS data 
%     e = errorbar(cellConD(ind), SO2D(ind), cellConDErr(ind), cellConDErr(ind), 'horizontal', 'LineStyle', 'none', 'Marker', 'none', 'Color', errorbarcolor, 'linewidth', errorbarlinewidthH); hold on; 
%     e.CapSize = capsize; 
%     % Expected 
%     et = plot(cellConLM(ind), cellConLMconv(ind), ':k'); et.LineWidth = trendlinewidthmain; hold on; 
%     % For linear model data 
%     e = errorbar(cellConLM(ind), SO2LMrOffset(ind), cellConLMErr(ind), cellConLMErr(ind), 'horizontal', 'LineStyle', 'none', 'Marker', 'none', 'Color', errorbarcolor, 'linewidth', errorbarlinewidthV); hold on; 
%     e.CapSize = capsize; 
%     % DOAS SO2 fit coefficients 
%     d = errorbar(cellConD(ind), SO2D(ind), SO2ErrD(ind), SO2ErrD(ind), 'vertical','LineStyle', 'none', 'Marker', 'o', 'Color', markeredgecolor, 'linewidth', errorbarlinewidthV);
%     d.MarkerSize = markersize; 
%     d.MarkerFaceColor = markerfacecolorD; 
%     d.MarkerEdgeColor = markeredgecolor; 
%     d.CapSize = capsize; 
%     xlim([minX maxX]);
%     minY = -5e+17;
%     ylim([minY maxY]);
%     %     ylabel('(d)SCD SO_2 (molec/cm^2)')
%     PrintEA = 90 - unAngle(i); 
%     h = title(sprintf('%d°' , PrintEA)); h.FontSize = 12; 
%     % xlabel('Gas cell concentration (ppm·m)')
%     % % mEst and standard error 
%     % % Real and complex
%     % SO2LMr = mEstLM_v3_real(2,:); % SO2 is G2 (I0, SO2, O3, Ring, Bshift) 
%     % % SO2LMri = mEstLM_v3_realimag(2,:); % SO2 is G2 (I0, SO2, O3, Ring, Bshift) 
%     % SO2ErrLMr = stdErrLM_v3_real(2,:); 
%     % % SO2ErrLMri = stdErrLM_v3_realimag(2,:); 
%     % % Plots mEst 
%     % e = errorbar(cellConLM, SO2LMr, cellConLMErr, cellConLMErr, 'horizontal', 'LineStyle', 'none', 'Marker', 'none', 'Color', errorbarcolor, 'linewidth', errorbarlinewidthV); hold on; 
%     % e.CapSize = capsize; 
%     % Real
%     lr = errorbar(cellConLM(ind), SO2LMrOffset(ind), SO2ErrLMr(ind), SO2ErrLMr(ind), 'vertical', 'LineStyle', 'none', 'Marker', 'D', 'Color', markerfacecolor_v3, 'linewidth', errorbarlinewidthV);
%     lr.MarkerFaceColor = markerfacecolor_v3; 
%     lr.MarkerEdgeColor = markeredgecolor; 
%     lr.MarkerSize = markersize; 
%     lr.CapSize = capsize; 
%     % Real imaginary 
%     % lri = errorbar(cellConLM, SO2LMri, SO2ErrLMri, SO2ErrLMri, 'vertical', 'LineStyle', 'none', 'Marker', 'v', 'Color', markerfacecolor_v3, 'linewidth', errorbarlinewidthV);
%     % lri.MarkerFaceColor = 'none'; 
%     % lri.MarkerEdgeColor = markerfacecolorLMAW1; 
%     if i == 3
%         ylabel('(d)SCD SO_2 (molec/cm^2)')
%     else
%     end
%     if i == 7 
%        xlabel('Gas cell concentration (ppm·m)')
%     end 
%     Fig = gca; Fig.FontSize = 12; set(gcf,'color','w'); Fig.Box = 'on';
% end 
% %% Make any adjustments 
% % % SAVES FIGURE
% % fnOut = 'Results_angle';
% % fname = fullfile(outDirFig, fnOut);
% % saveas(gca, fname)
