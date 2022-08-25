% convolute.m 
% This function convolutes high resolution trace gas absorption cross sections to match the lower optical resolution of the spectrometer used to record a Hg-spectrum
% 
% C. Barrington August 2022
% 
% HOW TO USE this function: 
% Run this function from the command line with 'convolute(<dir>, <hg>, <hlines>, <varargin>)' where:
% <dir> is replaced by the path to where the input files <hg> and <varargin> can be found e.g., '/Users/username/Documents/inFiles/'
% <hg> is replaced by the file name (excluding the file extension) of the recorded Hg-spectrum e.g., 'hg’
%   The Hg-spectrum file must be .txt format and consist of two columns where the first column contains wavelength information and the second column the measured intensity 
%   The Hg-spectrum must contain correct wavelength information (calibrate this spectrum before using this function should the wavelength information be incorrect) 
% <hlines> is the number of header lines contained in the <hg> before the start of the spectral data e.g., 14 (<hlines> should be 0 if there are no header lines)
% <varargin> should be replaced by any number of file names (excluding the file extension) corresponding to higher resolution trace gas absorption cross sections saved in <dir> separated by comma e.g., 'SO2', 'O3', 'BrO' 
%    All trace gas absorption cross sections files must be .txt format and consist of two columns where the first column contains wavelength information and the second column the measured intensity 
% The output variable [cxs] is a n-by-m matrix where n is equal to the number of channels of the spectrometer used to record the Hg-spectrum and m the number of input variables used for <varargin> + 1 where the first column is wavelength information and the 
%   Subsequent columns the convoluted trace gas absorption cross sections for each trace gas in the order they were input 
%   E.g., if <varargin> is 'SO2', 'O3', 'BrO', [cxs] will have four columns [lamba, SO2, O3, BrO] and the number of rows will equal the number of rows in hg (excluding hlines)
% 
% Detailed DESCRIPTION of function: 
% This function loads higher resolution trace gas cross sections according to <dir> and <varargin> and queries them according to the wavelength information in <hg> using linear interpolation of the values at neighbouring grid points
% The function plots the recorded Hg-spectrum and prompts the to define the lower (l) and upper (h) wavelength limits used to define the convolution data kernel k(x)
% The function defines k(x) according to the same equation used in DOAS Intelligent System (DOASIS, Kraus, 2006): 
%   k(x) = (s(x)-S_min)/(∫_l^h s(x)- S_min dx) where S_min = min s(x)|l ≤ x ≤ h
% The high-resolution trace gas absorption cross sections are then convoluted with k(x) to produce trace gas cross sections which match the lower optical resolution the spectrometer used to record the Hg-spectrum


% START OF FUNCTION
function [cxs] = convolute(dir, hg, hlines, varargin)

% LOADS HG-SPECTRUM 
suf = '.txt'; % Defines file extension

inHg = fullfile(dir, [hg suf]); inHg = char(inHg); % Defines Hg-spectrum filename according to <dir> and <hg>

    if hlines >= 1 % If headerlines exist in the Hg-spectrum 
       [lambda, iHg] = textread(inHg,'%f %f', 'headerlines', hlines); % Skip header lines and load file
    else
       [lambda, iHg] = textread(inHg,'%f %f'); % If there are no header lines load file 
    end  

% PLOTS HG-SPECTRUM    
figure('Renderer', 'painters', 'Position', [900 900 900 600]) 
p = plot(lambda, iHg, '-b'); hold on; p.LineWidth = 1; % Plots hg-spectrum as blue line 
xlim([min(lambda) max(lambda)]); set(gca,'XMinorTick','on'); xlabel('\lambda (nm)'); % X-axis 
ylabel('Intensity (counts)'); % Y-axis 
legend({'Recorded Hg-spectrum'},'FontSize', 14); title('Select Hg-peak to use for data kernel k(x)'); Fig = gca; Fig.FontSize = 14; set(gcf,'color','w');

% DEFINES DATA KERNAL k(x) 
% Promts user to input lower wavelength limit  
pl = sprintf('Select a single peak of the Hg spectrum to be used for the data kernel.\n Select a peak which is not saturated but in close proximity to the analysis window.\n Input the x-axis value corresponding to the start of the selected Hg-peak (l):\n');
l = input(pl);

% Promts user to input higher wavelength limit  
ph = sprintf('Input the x-axis value corresponding to the end of the selected Hg-peak (h):\n');
h = input(ph);

% Plots lower and higher wavelength limits
p = xline(l,'-r'); p.LineWidth = 1; p = xline(h, '-r'); p.LineWidth = 1;

% Extracts data s(x) from k(x)
[r, ~] = find(lambda > l&lambda < h); % Returns index of wavelengths between defined limits of kx (iL and iH)
sxW = lambda(r, :); % Extracts wavelength data for kernel 
sxI = iHg(r, :); % Extracts intensity data for kernel  

% Plots data kernel s(x) with Smin
figure('Renderer','painters','Position',[900 900 900 600]) 

subplot(1,3,1) % Plots data kernel s(x) with Smin as red cross
p = plot(sxW, sxI); hold on; p.LineWidth=1;
[smin, sminInd] = min(sxI);
plot(sxW(sminInd, :), smin, 'xr');
xlim([min(sxW) max(sxW)]); xlabel('\lambda (nm)');  % X-axis 
ylabel('Intensity (counts)'); % Y-axis 
legend({'s(x) from l to h'},'FontSize', 14, 'Location', 'NorthOutside'); Fig = gca; Fig.FontSize = 14; set(gcf,'color','w');

subplot(1,3,2) % Plots data kernel s(x) - Smin 
sxISmin = sxI - smin;
p = plot(sxW, sxISmin); hold on; p.LineWidth=1;
xlim([min(sxW) max(sxW)]); xlabel('\lambda (nm)');% X-axis 
ylabel('Intensity (counts)'); % Y-axis
legend({'s(x) - S_m_i_n from l to h'},'FontSize', 14, 'Location', 'NorthOutside');    
Fig = gca; Fig.FontSize = 15; set(gcf,'color','w');

subplot(1,3,3) % Plots normalised data kernel k(x) = (s(x)-S_min)/(∫_l^h s(x)- S_min)
% 'trapz' computes the approximate integral of Y via the trapezoidal method with unit spacing. If Y is a vector, then trapz(Y) is the approximate integral of Y
kx = sxISmin/(trapz(sxISmin)); 
p = plot(sxW,kx); hold on; p.LineWidth = 1;
xlim([min(sxW) max(sxW)]); xlabel('\lambda (nm)'); % X-axis 
ylabel('Normalised intensity'); % Y-axis 
legend({'s(x) - S_m_i_n / ‚à´_l^h s(x) - S_m_i_n'},'FontSize', 14,'Location', 'NorthOutside');   
t = sgtitle('Defining data kernel k(x) for convolution'); t.FontWeight = 'bold'; Fig = gca; Fig.FontSize = 14; set(gcf,'color','w');

% CREATES OUTPUT VARIABLE (cxs) 
cxs = zeros(length(lambda), length(varargin)+1); % Pre-allocates space for output variable based on wavelength data and number of trace gas cross sections
cxs(:, 1) = lambda; % Populates output variable with wavelength information

    for i = 1:length(varargin) % For each trace gas cross section 
        tgs = char(varargin{i}); % Selects trace gas filename according to <varargin>
        dSolar = fullfile(dir, [tgs suf]); inxs = char(dSolar); [wXs, aXs] = textread(inxs,'%f %f'); % Defines filename according to <dir> and <varargin>
        xs = interp1(wXs, aXs, lambda, 'linear'); % Queries absorption data according to the wavelength information in <hg> using linear interpolation of the values at neighbouring grid points
        cxsi = conv(xs, kx, 'same'); % Absorption data is convoluted with k(x) to produce trace gas cross section which matches the lower optical resolution the spectrometer used to record the Hg-spectrum 
        cxs(:, i+1) = cxsi; % Populates output variable with convoluted trace gas absorption cross section data 
    end 
end
% END OF FUNCTION 








