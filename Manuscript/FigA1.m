% % Figures for Introduction 

% Fig 1 (a and b) 
% a 
% Plots convoluted trace gas cross sections 
inDirCS = '/Users/charlotteb/Documents/Chapter 2/Cross sections/General/Used in figure 1/';

inDirCSx = '/Users/charlotteb/Documents/Chapter 2/Synthetic spectra/inFiles/';
% 
inDirCSDOASIS = '/Users/charlotteb/Documents/Chapter 2/Introduction/DOASIS/Output/'; 
inDirHg = '/Users/charlotteb/Documents/Chapter 2/Cross sections/Hg spectra/'; 
outDir = '/Users/charlotteb/Documents/Chapter 2/Introduction/outFiles/'; 
outDirFig = '/Users/charlotteb/Documents/Chapter 2/Introduction/outFiles/Figures/'; 
% IMPORTANT: The Hg-spectrum MUST be calibrated so that the wavelength information associated with each channel of the spectrometer is correct. All files should be in '.txt' format and contain 
% wavelength information in the first column and absorption OR intensity information the second. For the trace gas cross sections (i and ii) and the solar reference (iii) it is assumed that the
% .txt files DO NOT contain header or footer lines but only the wavelength and absorption OR intensity data. Since the Hg-spectrum should be measured, it is assumed that the .txt file will contain header 
% lines. Please determine the number of header lines (hl) below:
hl = 14; % Change to the row number before start of data. If there are no header lines enter '0'. 
% Reference spectra 
fnHg = 'Reference_FLMS195681__0__12-14-16-726.txt'; % Hg-spectrum 
fnSO2 = 'SO2_Bogumil_2003_293K_vac.txt'; 
fnO3 = 'O3_Gorshelev_2014_243K_vac.txt'; 
fnSO2DOASIS = 'SO2_Bogumil_2003_293K_vac_DOASISconv.txt'; 
fnO3DOASIS = 'O3_Gorshelev_2014_243K_vac_DOASISconv.txt'; 
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
dSO2DOASIS = fullfile(inDirCSDOASIS, fnSO2DOASIS); inSO2DOASIS = char(dSO2DOASIS); [wSO2DOASIS, aSO2DOASIS] = textread(inSO2DOASIS,'%f %f');
dO3DOASIS = fullfile(inDirCSDOASIS, fnO3DOASIS); inO3DOASIS = char(dO3DOASIS); [wO3DOASIS, aO3DOASIS] = textread(inO3DOASIS,'%f %f'); 
% % Defines wavelength range of interest 
% lAw = 279; % Set slightly lower so x-axis on figure may be set to 280 nm 
% uAw = 421; 
% % RESAMPLES DATA AT WAVELENGTHS DEFINED BY CALIBRATED Hg spectrum
% % Crops wHg to limit wavelengths to those 280 nm or above
% [ind, ~] = find(wHg >= lAw);
% wHg = wHg(ind, 1); 
% iHg = iHg(ind, 1); 
% Crops wHg to limit wavelengths to those 420 nm or below
% [ind, ~] = find(wHg <= uAw);
% wHg = wHg(ind, 1); 
% iHg = iHg(ind, 1); 
% 
% [ind, ~] = find(wO3 <= uAw);
% wO3 = wO3(ind, 1); 
% aO3 = aO3(ind, 1); 
% 
% [ind, ~] = find(wSO2 <= uAw);
% wSO2 = wSO2(ind, 1); 
% aSO2 = aSO2(ind, 1); 

% RESAMPLES DATA AT WAVELENGTHS DEFINED BY CALIBRATED HG-SPECTRUM
% Uses interp(x, v, xq) to query trace gas cross sections and solar reference at wavelengths defined by calibrated Hg-spectrum   
% 'interp' returns interpolated values of a 1-D function at specific query points (xq) using linear interpolation 
% Vector x contains the sample points - original wavelength information 
% Vector v contains the correpsonding values, v(x) - original absorption (OR intensity) information 
% SO2 = interp1(wSO2, aSO2, wHg,'cubic'); 
% O3 = interp1(wO3, aO3, wHg,'cubic'); 
% % CONVOLUTION
% % Convolutes resampled trace gas cross sections and solar reference with data kernel to match ILS of spectrometer 
% figure('Renderer', 'painters', 'Position', [900 900 900 600]) % Plots Hg-spectrum 
% dp = plot(wHg, iHg); hold on; dp.LineWidth = 1;
% xlim([min(wHg) max(wHg)]); xtickn = (290:10:420); xticks(xtickn);
% set(gca,'XMinorTick','on'); 
% xlabel('\lambda (nm)'); ylabel('Intensity (counts)');
% legend({'Recorded Hg-spectrum','Wavelength range of analysis window'},'FontSize',14);  
% title('Selecting Hg-peak to use for data kernel k(x)'); 
% Fig = gca; Fig.FontSize = 14; set(gcf,'color','w');
% % Defines data kernel k(x)
% pL = sprintf('Select a single peak of the Hg-spectrum to be used for the data kernel.\n Select a peak which is not saturated but in close proximity to the analysis window.\n Input the x-axis value corresponding to the start of the selected Hg-peak (l):\n');
% l = input(pL); % use 311
% pH = sprintf('Input the x-axis value corresponding to the end of the selected Hg-peak (h):\n');
% h = input(pH); % use 315 
% dp = xline(l,'-r'); dp.LineWidth = 1;
% dp = xline(h, '-r'); dp.LineWidth = 1;
% % Extracts data s(x) from k(x)
% [rI, ~] = find(wHg > l&wHg < h); % Returns index of wavelengths between defined limits of kx (iL and iH)
% sxW = wHg(rI, :); % Extracts wavelength data for kernel 
% sxI = iHg(rI, :); % Extracts intensity data for kernel  
% % Plots data kernel s(x) with Smin
% figure('Renderer','painters','Position',[900 900 900 600]) 
% subplot(1,3,1) % Plots data kernel s(x) with Smin as red cross
% dp = plot(sxW, sxI); hold on; dp.LineWidth=1;
% [smin, sminInd] = min(sxI);
% plot(sxW(sminInd, :), smin, 'xr');
% xlim([min(sxW) max(sxW)]); xlabel('\lambda (nm)'); ylabel('Intensity (counts)'); legend({'s(x) from l to h'},'FontSize', 14, 'Location', 'NorthOutside');    
% Fig = gca; Fig.FontSize = 15; set(gcf,'color','w');
% subplot(1,3,2) % Plots data kernel s(x) - Smin 
% sxISmin = sxI - smin;
% dp = plot(sxW, sxISmin); hold on; dp.LineWidth=1;
% xlim([min(sxW) max(sxW)]);
% xlabel('\lambda (nm)'); ylabel('Intensity (counts)'); legend({'s(x) - S_m_i_n from l to h'},'FontSize', 14, 'Location', 'NorthOutside');    
% Fig = gca; Fig.FontSize = 15; set(gcf,'color','w');
% subplot(1,3,3) % Plots normalised data kernel k(x) =  s(x) - Smin / ∫l^h s(x) - Smin
% % 'trapz' computes the approximate integral of Y via the trapezoidal method with unit spacing. If Y is a vector, then trapz(Y) is the approximate integral of Y
% kx = sxISmin/(trapz(sxISmin)); 
% Nr = normalize(sxISmin,'range');
% 
% dp = plot(sxW,kx); hold on; dp.LineWidth = 1;
% xlim([min(sxW) max(sxW)]); xlabel('\lambda (nm)'); ylabel('Normalised intensity');
% legend({'s(x) - S_m_i_n / ∫_l^h s(x) - S_m_i_n'},'FontSize', 14,'Location', 'NorthOutside');   
% t = sgtitle('Defining data kernel k(x) for convolution'); t.FontWeight = 'bold'; 
% Fig = gca; Fig.FontSize = 15; set(gcf,'color','w');
%  
% % figure 
% % dp = plot(sxW, Nr); hold on; dp.LineWidth=1;
% % 
% 
% % % % Convolutes resampled trace gas cross sections
% % % % w = conv(u, v, shape) returns a subsection of the convolution, as specified by shape. For example, conv(u, v, 'same') returns only the central part of the convolution, the same size as u
% cSO2 = conv(aSO2, kx); % , 'same'); 
% cO3 = conv(aO3, kx); % 'same'); 
% 
% % 
% % cSO2 = conv(aSO2, Nr, 'same'); 
% % cO3 = conv(aO3, Nr, 'same'); 
% 
% 
% 
% figure 
% dp = plot(cSO2, '-'); dp.LineWidth = 1.5; dp.Color = colSO2
% 
% 
% % 
% %  
% % figure
% % subplot(2,1,1)
% % dp = plot(cSO2, '-'); dp.LineWidth = 1.5; dp.Color = colSO2; hold on 
% % % dp = plot(wSO2DOASIS, aSO2DOASIS, '-K'); dp.LineWidth = 1.5; 
% % ylabel('Mean absoprtion (cm^2/molec)', 'FontSize', 12); 
% % xlabel('\lambda (nm)'); hold on 
% % % ylim([0 3e-17]);
% % xlim([280 420]);
% % Fig = gca; Fig.FontSize = 13; set(gcf,'color','w');
% % Fig.Box = 'on'; 
% % subplot(2,1,2)
% % dp = plot(cO3, '-'); dp.LineWidth = 1.5; dp.Color = colO3; hold on 
% % % dp = plot(wO3DOASIS, aO3DOASIS, '-K'); dp.LineWidth = 1.5; 
% % ylabel('Mean absoprtion (cm^2/molec)', 'FontSize', 12); 
% % xlabel('\lambda (nm)'); hold on 
% % % ylim([0 3e-17]);
% % xlim([280 420]);
% % Fig = gca; Fig.FontSize = 13; set(gcf,'color','w');
% % Fig.Box = 'on'; 
% % % Save convoluted cross sections to be used later 
% % fname = fullfile(outDir, 'FigA1_Convolute_conv');
% % save(fname); % Saves all variables  
% % % 
% % % 
% % % 
% % % % 
% % % % % cSO2 = conv(aSO2, kx, 'same'); 
% % % % % cO3 = conv(aO3, kx, 'same'); 
% % % % % 
% % % % % SO2 = interp1(wSO2, cSO2, wHg,'linear'); 
% % % % % O3 = interp1(wO3, cO3, wHg,'linear');  
% % % % 
% % % % % PLOTS RESAMPLED AND CONVOLUTED ABSORPTION CROSS SECTIONS AND SOLAR REFERENCE 
% % % % warning('off','MATLAB:Axes:NegativeDataInLogAxis') % Turns off warning ignoring negative data in log axis 
% % % % colSO2 = 'r'; % Defines colours for each trace gas
% % % % colO3 = 'r';
% % % % figure
% % % % subplot(2,1,1)
% % % % dp = plot(wHg, cSO2, '-'); dp.LineWidth = 1.5; dp.Color = colSO2; hold on 
% % % % dp = plot(wSO2DOASIS, aSO2DOASIS, '-K'); dp.LineWidth = 1.5; 
% % % % ylabel('Mean absoprtion (cm^2/molec)', 'FontSize', 12); 
% % % % xlabel('\lambda (nm)'); hold on 
% % % % % ylim([0 3e-17]);
% % % % xlim([280 420]);
% % % % Fig = gca; Fig.FontSize = 13; set(gcf,'color','w');
% % % % Fig.Box = 'on'; 
% % % % subplot(2,1,2)
% % % % dp = plot(wHg, cO3, '-'); dp.LineWidth = 1.5; dp.Color = colO3; hold on 
% % % % dp = plot(wO3DOASIS, aO3DOASIS, '-K'); dp.LineWidth = 1.5; 
% % % % ylabel('Mean absoprtion (cm^2/molec)', 'FontSize', 12); 
% % % % xlabel('\lambda (nm)'); hold on 
% % % % % ylim([0 3e-17]);
% % % % xlim([280 420]);
% % % % Fig = gca; Fig.FontSize = 13; set(gcf,'color','w');
% % % % Fig.Box = 'on'; 
% % % % % Save convoluted cross sections to be used later 
% % % % fname = fullfile(outDir, 'FigA1_Convolute_conv');
% % % % save(fname); % Saves all variables  
% % % % %% Make any adjustments 
% % % % % SAVES FIGURE
% % % % fnOut = 'FigA1_Convolute';
% % % % fname = fullfile(outDirFig, fnOut);
% % % % saveas(gca, fname)
% % % 
% % % 
