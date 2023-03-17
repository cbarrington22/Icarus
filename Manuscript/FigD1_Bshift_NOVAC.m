% Bshift 
% Run Introduction_Figures.m first to I0 
inDir = '/Users/charlotteb/Documents/Chapter 2/NOVAC spectra/outFiles/'; 
fnIn = 'Fig1_Fig3_conv_I0_SO2_O3_Ring_NOVAC.mat'; 
outDir = '/Users/charlotteb/Documents/Chapter 2/NOVAC spectra/outFiles/';
outDirFig = '/Users/charlotteb/Documents/Chapter 2/NOVAC spectra/outFiles/Bshift/';
inFile = fullfile(inDir, fnIn); 
load(inFile)
% Calculates Bshift spectrum from convoluted solar reference (Introduction_Figures.m)
% B3
outyP = zeros(length(I0), 1); 
h = 1; % Channels 
y = I0; 
for yInd = 3:length(y)
    yPrime = (1/(2*h))*(y(yInd-2)-8*y(yInd-1)+8*y(yInd)-(yInd+1));
    outyP(yInd, :) = yPrime;
end
figure 
plot(lambda, outyP); % Plots y' according to B3
xlabel('Channel')
ylabel('dy')
Fig = gca; Fig.FontSize = 14; set(gcf,'color','w');
% Bshift 
yPrime = outyP;
Bshift = yPrime./y; 
figure 
dp = plot(lambda, Bshift, '-k'); hold on; dp.LineWidth = 0.5; % Plots Bshift 
ylim([0 max(Bshift)]);
ylabel('Intensity'); xlim([min(lambda) max(lambda)]); xlabel('\lambda (nm)');
Fig = gca; Fig.FontSize = 14; set(gcf,'color','w');
% Save convoluted cross sections to be used later 
fname = fullfile(outDir, 'Bshift');
save(fname); % Saves all variables  
% Changes name
fname = fullfile(outDir, 'Fig1_Fig3_FigD1_conv_I0_SO2_O3_Ring_Bshift_NOVAC');
save(fname, 'fnHg', 'fnSolar', 'I0', 'fnSO2', 'SO2', 'fnO3', 'O3', 'fnRing' , 'Ring', 'lambda', 'kx', 'Bshift', 'hlMasaya', 'flMasaya','iMasaya1', 'ind', 'Fs'); 
%% Make any adjustments 
% SAVES FIGURE
fnOut = 'FigD1_Bshift_NOVAC';
fname = fullfile(outDirFig, fnOut);
saveas(gca, fname)
