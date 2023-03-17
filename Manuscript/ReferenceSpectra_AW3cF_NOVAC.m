
% Reference spectra - now modified for NOVAC scans (Masaya) 
% Model version 
% (1) AW3 cropped spatial frequencies 
% - Real 
% - Complex magnitude - now only using real part 

addpath '/Users/charlotteb/Documents/Chapter 2/NOVAC spectra/BootstrapFunction/'

% LOADS ELEMENTS OF DESGN MATRIX G from ReferenceSpectra_main.m 
inDirG = '/Users/charlotteb/Documents/Chapter 2/NOVAC spectra/outFiles/';
outDir = '/Users/charlotteb/Documents/Chapter 2/NOVAC spectra/AW3/';
outDirFig = '/Users/charlotteb/Documents/Chapter 2/NOVAC spectra/AW3/Figures/';
inFile = 'WT_G_NOVAC.mat'; 
inFile = fullfile(inDirG, inFile); 
load(inFile)

% SAVES LAMBDA - WAVLENEGTH RANGE OVER WHICH CWT IS COMPUTED
lambdaS = lambda; 
FS = F; 

% Analysis window
% DEFINES ANALYSIS WINDOW 
% WAVELENGTH 
AW_L = 314;
AW_H = 330;
[~, idxL]= min(abs(lambda-AW_L));
wavL = lambda(idxL); % Low wavelength limit of AW
[~, idxH]= min(abs(lambda-AW_H));
wavH = lambda(idxH);  % High wavelength limit of AW
% SPATIAL FREQUENCY 
AW_Ft = 0.01; % max(F);
AW_Fb = 0.0015; % min(F);
[~, idxFt]= min(abs(F-AW_Ft));
Ft = F(idxFt); % Low wavelength limit of AW
[~, idxFb]= min(abs(F-AW_Fb));
Fb = F(idxFb);  % High wavelength limit of AW

% G 
% EXTRACTS ANALYSIS WINDOW 
F = F(idxFt:idxFb);
lambda = lambda(idxL:idxH);
WT_I0 = WT_I0(idxFt:idxFb, idxL:idxH);
WT_SO2 = WT_SO2(idxFt:idxFb, idxL:idxH);
WT_O3 = WT_O3(idxFt:idxFb, idxL:idxH);
WT_Ring = WT_Ring(idxFt:idxFb, idxL:idxH);
WT_Bshift = WT_Bshift(idxFt:idxFb, idxL:idxH);

% figure % Plots design matrix (for checking only)
% subplot(3,2,1) 
% subtitle('J_0(\lambda) prime'); hold on; 
% pcolor(lambda, F, real(WT_I0)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
% xlim([min(lambda) max(lambda)]); set(gca,'YScale','log'); ylim([min(F) max(F)]);
% set(gcf,'color','w'); 
% subplot(3,2,2) 
% subtitle('\sigma_S_O_2 prime'); hold on; 
% pcolor(lambda, F, real(WT_SO2)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
% xlim([min(lambda) max(lambda)]); set(gca,'YScale','log'); ylim([min(F) max(F)]);
% set(gcf,'color','w'); 
% subplot(3,2,3) 
% subtitle('\sigma_O_3 prime'); hold on; 
% pcolor(lambda, F, real(WT_O3)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
% xlim([min(lambda) max(lambda)]); set(gca,'YScale','log'); ylabel('Spatial frequency (cycles/\lambda)', 'FontSize', 14);
% set(gcf,'color','w'); ylim([min(F) max(F)]);
% subplot(3,2,4) 
% subtitle('Ring prime'); hold on; 
% pcolor(lambda, F, real(WT_Ring)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
% xlim([min(lambda) max(lambda)]); set(gca,'YScale','log'); ylim([min(F) max(F)]);
% cb = colorbar; cb.Label.Position(1) = 3.5;
% ylabel(cb,'Amplitude', 'Rotation', -90, 'FontSize', 10);
% set(gcf,'color','w'); 
% subplot(3,2,5) 
% subtitle('B_s_h_i_f_t prime'); hold on; 
% pcolor(lambda, F, real(WT_Bshift)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
% xlim([min(lambda) max(lambda)]); set(gca,'YScale','log'); ylim([min(F) max(F)]);
% set(gcf,'color','w'); xlabel('\lambda (nm)', 'FontSize', 14); 

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
fname = fullfile(outDir, 'ReferenceSpectra_AW3cF_G_NOVAC');
save(fname, 'G_real', 'G_realimag', 'Gt_real', 'Gt_realimag', 'rI', 'lambda', 'WT_I0', 'WT_SO2', 'WT_O3', 'WT_Ring', 'WT_Bshift', 'l', 'h' , 'AW_L', 'AW_H', 'AW_Ft', 'AW_Fb', 'idxCoi', 'wavL', 'wavH', 'Ft', 'Fb','idxFt', 'idxFb', 'idxL', 'idxH'); 

% d
% TWO FULL INSTRUMENT SCANS FROM MASAYA (RUN TWO LOOPS)
% SCAN 1417 
% Loads example d 
inDird = '/Users/charlotteb/Documents/Chapter 2/NOVAC spectra/inFiles/Masaya/D2J2375_140319_1417_0/'; % Directory to recorded spectra
clear dir 
subFolders=dir(inDird); % Returns data on subfolders to be analysed
subFoldersData=subFolders(~ismember({subFolders(:).name},{'.','..','.DS_Store'})); % Removes '.' and '..' and '.DS_Store' contained in subFolder varaible
spectraList=extractfield(subFoldersData,'name');spectraList=spectraList'; % Extracts file names to list all spectra to be analysed
spectraList=spectraList(1:end); % Removes 'subset' folder 
% EXTRACTS KNOWN CELL SO2 CONCENTRATION FROM FILENAME - infomation now inside spectra file (move later) 
% findStr = "_FLMS195681";
% tempStr = string(extractBefore(spectraList, findStr)); 
% findStr = "_";
% elevAn = str2double(extractBefore(tempStr, findStr)); % ELEVATION ANGLE OF SCAN 
% cellCon = str2double(extractAfter(tempStr, findStr)); % GAS CELL SO2 CONCENTRATION
% hlines = 14; 
% ALLOCATES SPACE 
Results_mEst_real = NaN(width(G_real), length(spectraList)); 
% Results_mEst_realimag = NaN(width(G_realimag), length(spectraList)); 
Results_res_real = NaN(1, length(spectraList)); 
% Results_res_realimag = zeros(1, length(spectraList)); 
Results_stdErr_real = NaN(width(G_real), length(spectraList)); 
elevAn = NaN(1, length(spectraList)); 
% Results_stdErr_realimag = NaN(width(G_realimag), length(spectraList)); 
% ANALYSIS LOOP 
% % Turns off warnings 
% x1 = [];
% x2 = [];
% warning('off','MATLAB:nearlySingularMatrix') 
% tic
% figure 
% for s = 1:length(spectraList)
%     if s == 2 % don't analyse dark scan with is "00001_0.STD", second in list 
%         else
%         inI = fullfile(inDird, spectraList(s)); inI = char(inI);    
%         if hlMasaya >= 1
%            [iMasaya] = textread(inI,'%s', 'headerlines', hlMasaya); % Skips header lines if exist 
%         else
%            [iMasaya] = textread(inI,'%s');
%         end
%         elevAn(:, s) = str2double(iMasaya(2090,:));  
%         iMasaya = iMasaya(1:flMasaya, :); iMasaya = str2double(iMasaya);
%         % EXTRACTS WAVELENGTHS used to build G in Reference_Spectra_main.m  
%         I = iMasaya(ind:end); % Should be equal to lambdaS
%         tf = isequal(length(lambdaS), length(I)); 
%         if tf == 1
%         else 
%            disp('Warning wavelength information of measurement spectrum differs to that in G')  
%         end 
%         % APPLIES LOG 
%         J = log(I);
%         % CWT 
%         [WT_J, FI, ~] = cwt(J, Fs); 
%         tf = isequal(FS, FI); 
%         if tf == 1
%         else 
%            disp('Frequency scales for measurement spectrum differ to that in G')  
%         end 
%         % REMOVES COI
%         WT_J(idxCoi) = NaN; 
%         % EXTRACTS ANALYSIS WINDOW 
%         WT_J = WT_J(idxFt:idxFb, idxL:idxH);   
% %         if s == 28 % Index for target d (2080 ppmm at zenith) - same used for figures - now index for zenith scan "00027_0.STD"
% %            figure
% %            title('0°'); hold on; 
% %            pcolor(lambda, F, real(WT_J)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
% %            xlim([min(lambda) max(lambda)]); set(gca,'YScale','log'); ylabel('Spatial frequency (cycles/\lambda)', 'FontSize', 14);
% %            set(gcf,'color','w'); ylim([min(F) max(F)]); xlabel('\lambda (nm)', 'FontSize', 14);
% %            cb = colorbar; cb.Label.Position(1) = 3.5;
% %            ylabel(cb,'Amplitude', 'Rotation', -90, 'FontSize', 10);
% %            box on
% %            % SAVES FIGURE
% %            fnOut = 'ReferenceSpectra_AW3cF_dExample_NOVAC_zenith';
% %            fname = fullfile(outDirFig, fnOut);
% %            saveas(gca, fname)
% %         else 
% %         end
%         % DEFINES REAL AND IMAGINARY PARTS OF CWT 
%         WT_J_real = real(WT_J); 
%         % WT_J_imag = imag(WT_J); 
%         % CALCULATION OF ERROR WITH BOOTSTRAP FUNCTION 
%         % REAL 
%         dr = WT_J_real'; 
%         G1r = WT_I0_real'; 
%         G2r = WT_SO2_real';
%         G3r = WT_O3_real'; 
%         G4r = WT_Ring_real';
%         G5r = WT_Bshift_real'; 
%         BOOTSTAT = bootstrp(1000, @standardErrorReal_NOVAC, dr, G1r, G2r, G3r, G4r, G5r);
%         G1rErr = std(BOOTSTAT(:, 1)); 
%         G2rErr = std(BOOTSTAT(:, 2)); 
%         G3rErr = std(BOOTSTAT(:, 3)); 
%         G4rErr = std(BOOTSTAT(:, 4));
%         G5rErr = std(BOOTSTAT(:, 5)); 
%         Results_stdErr_real(:,s) = [G1rErr; G2rErr; G3rErr; G4rErr; G5rErr]; 
%         % REALIMAG
%         % dr, G1r, G2r, G3r, G4r, G5r already defined 
%         % di = WT_J_imag'; 
%         % G1i = WT_I0_imag'; 
%         % G2i = WT_SO2_imag';
%         % G3i = WT_O3_imag'; 
%         % G4i = WT_Ring_imag';
%         % G5i = WT_Bshift_imag'; 
%         % BOOTSTAT = bootstrp(1000, @standardErrorRealimag, dr, G1r, G2r, G3r, G4r, G5r, di, G1i, G2i, G3i, G4i, G5i);
%         % G1riErr = std(BOOTSTAT(:, 1)); 
%         % G2riErr = std(BOOTSTAT(:, 2)); 
%         % G3riErr = std(BOOTSTAT(:, 3)); 
%         % G4riErr = std(BOOTSTAT(:, 4));
%         % G5riErr = std(BOOTSTAT(:, 5)); 
%         % Results_stdErr_realimag(:,s) = [G1riErr; G2riErr; G3riErr; G4riErr; G5riErr];          
%         % VECTORISES 
%         % Real
%         WT_J_real_v = WT_J_real(:); 
%         % imag
%         % WT_J_imag_v = WT_J_imag(:); 
%         % Combines real and imaginary part for complex magnitude 
%         % WT_J_realimag = [WT_J_real_v; WT_J_imag_v]; 
%         % d 
%         d_real = WT_J_real_v;
%         % d_realimag = WT_J_realimag; 
%         % Remove NaN 
%         d_real = d_real(all(~isnan(d_real), 2),:);
%         % d_realimag = d_realimag(all(~isnan(d_realimag), 2),:);
%         % FINDS MEST 
%         % A 
%         A_real = Gt_real * G_real;
%         % A_realimag = Gt_realimag * G_realimag;
%         % B 
%         % Real only
%         B_real = Gt_real * d_real;
%         % Real and imag (complex magnitude) 
%         % B_realimag = Gt_realimag * d_realimag;
%         % SAVES to allocated space 
%         % mEst 
%         % Real only 
%         mEst_real = A_real\B_real;
%         % Real and imag (complex magnitude) 
%         % mEst_realimag = A_realimag\B_realimag;
%         % Residuals (e and E)
%         % Real only 
%         e_real = d_real - (G_real * mEst_real);
%         x1 = [x1, e_real]; 
%         res_real = e_real' * e_real;   
% %         % PLOT RESIDUALS FOR CERTAIN SCANS - this code only works for real part!
% %         if s == 14 || s == 28 || s == 35 || s == 49 % Plots residuals for -39, zenith, 39 and 79 degrees 
% %             % CALCULATE RESIDUALS 
% %             idxCoiT = idxCoi(idxFt:idxFb, idxL:idxH); % Extracts analysis window from COI template 
% %             idxCoiT = idxCoiT(:); 
% %             idxCoiT = double(idxCoiT); % Converts to double
% %             rIndR = 1; % Index for moving through rows of R
% %             for rIndRv = 1:length(idxCoiT) % For each row of idxCoiT
% %               if idxCoiT(rIndRv,:) == 1 % If cell is inside COI
% %                  idxCoiT(rIndRv,:) = NaN; % Convert to 'NaN' 
% %               elseif idxCoiT(rIndRv,:) == 0 % But if cell is outside COI
% %                  y = rIndR; % Assigns value 1
% %                  rIndR = rIndR + 1; % Adds on to index to move on one place
% %                  idxCoiT(rIndRv, :) = e_real(y, :); % Populates template with the corresponding residual value
% %               else
% %               end
% %             end
% %            residuals = reshape(idxCoiT, length(F), length(lambda)); % Reshapes compiled residuals back into a matrix for plotting
% %            % PLOTS RESIDUALS 
% %            if s == 14
% %                subplot(2,2,1); hold on 
% %            elseif s == 28
% %                subplot(2,2,2); hold on 
% %            elseif s == 35
% %                subplot(2,2,3); hold on 
% %            elseif s == 49
% %                subplot(2,2,4); hold on  
% %            end 
% %            title(sprintf('%d°', elevAn(:, s))); hold on;     
% %            pcolor(lambda, F, residuals); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
% %            xlim([min(lambda) max(lambda)]); set(gca,'YScale','log'); ylabel('Spatial frequency (cycles/\lambda)', 'FontSize', 14);
% %            set(gcf,'color','w'); ylim([min(F) max(F)]); xlabel('\lambda (nm)', 'FontSize', 14);
% %            cb = colorbar; cb.Label.Position(1) = 3.5;
% %            ylabel(cb,'Amplitude', 'Rotation', -90, 'FontSize', 10);
% %            box on
% %            % SAVES FIGURE
% %            fnOut = 'ReferenceSpectra_AW3cF_Residual_NOVAC';
% %            fname = fullfile(outDirFig, fnOut);
% %            saveas(gca, fname)
% %         else
% %         end
%         % Real and imag (complex magnitude) 
%         % e_realimag = d_realimag - (G_realimag * mEst_realimag);
%         % res_realimag = e_realimag' * e_realimag; 
%         % SAVES to preallocated space 
%         Results_mEst_real(:,s) = mEst_real; 
%         % Results_mEst_realimag(:,s) = mEst_realimag; 
%         Results_res_real(:,s) = res_real; 
%         % Results_res_realimag(:,s) = res_realimag; 
%     end 
% end 
% % SAVES variables 
% fname = fullfile(outDir, 'ReferenceSpectra_AW3cF_mEst_NOVAC_Masaya_1417');
% save(fname, 'G_real', 'G_realimag', 'elevAn', 'lambda', 'F', 'Results_mEst_real', 'Results_res_real', 'Results_stdErr_real', 'wavL', 'wavH', 'Ft', 'Fb'); 
% toc 
% 
% % Plot results to check look OK 
% figure
% title('Masaya 1417'); hold on;     
% plot(elevAn, Results_mEst_real(2,:), 'o'); 
% xlabel('Scan elevation angle (° from zenith)'); 
% ylabel('mEst_S_O_2 (molec/cm^2)');
% xlim([min(elevAn) max(elevAn)])

% SCAN 1941 
% Loads example d 
inDird = '/Users/charlotteb/Documents/Chapter 2/NOVAC spectra/inFiles/Masaya/D2J2375_140319_1941_0/'; % Directory to recorded spectra
clear dir 
subFolders=dir(inDird); % Returns data on subfolders to be analysed
subFoldersData=subFolders(~ismember({subFolders(:).name},{'.','..','.DS_Store'})); % Removes '.' and '..' and '.DS_Store' contained in subFolder varaible
spectraList=extractfield(subFoldersData,'name');spectraList=spectraList'; % Extracts file names to list all spectra to be analysed
spectraList=spectraList(1:end); % Removes 'subset' folder 
% EXTRACTS KNOWN CELL SO2 CONCENTRATION FROM FILENAME - infomation now inside spectra file (move later) 
% findStr = "_FLMS195681";
% tempStr = string(extractBefore(spectraList, findStr)); 
% findStr = "_";
% elevAn = str2double(extractBefore(tempStr, findStr)); % ELEVATION ANGLE OF SCAN 
% cellCon = str2double(extractAfter(tempStr, findStr)); % GAS CELL SO2 CONCENTRATION
% hlines = 14; 
% ALLOCATES SPACE 
Results_mEst_real = NaN(width(G_real), length(spectraList)); 
% Results_mEst_realimag = NaN(width(G_realimag), length(spectraList)); 
Results_res_real = NaN(1, length(spectraList)); 
% Results_res_realimag = zeros(1, length(spectraList)); 
Results_stdErr_real = NaN(width(G_real), length(spectraList)); 
elevAn = NaN(1, length(spectraList)); 
% Results_stdErr_realimag = NaN(width(G_realimag), length(spectraList)); 
% ANALYSIS LOOP 
% Turns off warnings 
warning('off','MATLAB:nearlySingularMatrix') 
tic
figure
for s = 1:length(spectraList)
    if s == 2 % don't analyse dark scan with is "00001_0.STD", second in list 
        else
        inI = fullfile(inDird, spectraList(s)); inI = char(inI);    
        if hlMasaya >= 1
           [iMasaya] = textread(inI,'%s', 'headerlines', hlMasaya); % Skips header lines if exist 
        else
           [iMasaya] = textread(inI,'%s');
        end
        elevAn(:, s) = str2double(iMasaya(2090,:));  
        iMasaya = iMasaya(1:flMasaya, :); iMasaya = str2double(iMasaya);
        % EXTRACTS WAVELENGTHS used to build G in Reference_Spectra_main.m  
        I = iMasaya(ind:end); % Should be equal to lambdaS
        tf = isequal(length(lambdaS), length(I)); 
        if tf == 1
        else 
           disp('Warning wavelength information of measurement spectrum differs to that in G')  
        end 
        % APPLIES LOG 
        J = log(I);
        % CWT 
        [WT_J, FI, ~] = cwt(J, Fs); 
        tf = isequal(FS, FI); 
        if tf == 1
        else 
           disp('Frequency scales for measurement spectrum differ to that in G')  
        end 
        % REMOVES COI
        WT_J(idxCoi) = NaN; 
        % EXTRACTS ANALYSIS WINDOW 
        WT_J = WT_J(idxFt:idxFb, idxL:idxH);   
%         if s == 28 % Index for target d (2080 ppmm at zenith) - same used for figures - now index for zenith scan "00027_0.STD"
%            figure
%            title('0°'); hold on; 
%            pcolor(lambda, F, real(WT_J)); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
%            xlim([min(lambda) max(lambda)]); set(gca,'YScale','log'); ylabel('Spatial frequency (cycles/\lambda)', 'FontSize', 14);
%            set(gcf,'color','w'); ylim([min(F) max(F)]); xlabel('\lambda (nm)', 'FontSize', 14);
%            cb = colorbar; cb.Label.Position(1) = 3.5;
%            ylabel(cb,'Amplitude', 'Rotation', -90, 'FontSize', 10);
%            box on
%            % SAVES FIGURE
%            fnOut = 'ReferenceSpectra_AW3cF_dExample_NOVAC_zenith';
%            fname = fullfile(outDirFig, fnOut);
%            saveas(gca, fname)
%         else 
%         end
        % DEFINES REAL AND IMAGINARY PARTS OF CWT 
        WT_J_real = real(WT_J); 
        % WT_J_imag = imag(WT_J); 
        % CALCULATION OF ERROR WITH BOOTSTRAP FUNCTION 
        % REAL 
        dr = WT_J_real'; 
        G1r = WT_I0_real'; 
        G2r = WT_SO2_real';
        G3r = WT_O3_real'; 
        G4r = WT_Ring_real';
        G5r = WT_Bshift_real'; 
        BOOTSTAT = bootstrp(1000, @standardErrorReal_NOVAC, dr, G1r, G2r, G3r, G4r, G5r);
        G1rErr = std(BOOTSTAT(:, 1)); 
        G2rErr = std(BOOTSTAT(:, 2)); 
        G3rErr = std(BOOTSTAT(:, 3)); 
        G4rErr = std(BOOTSTAT(:, 4));
        G5rErr = std(BOOTSTAT(:, 5)); 
        Results_stdErr_real(:,s) = [G1rErr; G2rErr; G3rErr; G4rErr; G5rErr]; 
        % REALIMAG
        % dr, G1r, G2r, G3r, G4r, G5r already defined 
        % di = WT_J_imag'; 
        % G1i = WT_I0_imag'; 
        % G2i = WT_SO2_imag';
        % G3i = WT_O3_imag'; 
        % G4i = WT_Ring_imag';
        % G5i = WT_Bshift_imag'; 
        % BOOTSTAT = bootstrp(1000, @standardErrorRealimag, dr, G1r, G2r, G3r, G4r, G5r, di, G1i, G2i, G3i, G4i, G5i);
        % G1riErr = std(BOOTSTAT(:, 1)); 
        % G2riErr = std(BOOTSTAT(:, 2)); 
        % G3riErr = std(BOOTSTAT(:, 3)); 
        % G4riErr = std(BOOTSTAT(:, 4));
        % G5riErr = std(BOOTSTAT(:, 5)); 
        % Results_stdErr_realimag(:,s) = [G1riErr; G2riErr; G3riErr; G4riErr; G5riErr];          
        % VECTORISES 
        % Real
        WT_J_real_v = WT_J_real(:); 
        % imag
        % WT_J_imag_v = WT_J_imag(:); 
        % Combines real and imaginary part for complex magnitude 
        % WT_J_realimag = [WT_J_real_v; WT_J_imag_v]; 
        % d 
        d_real = WT_J_real_v;
        % d_realimag = WT_J_realimag; 
        % Remove NaN 
        d_real = d_real(all(~isnan(d_real), 2),:);
        % d_realimag = d_realimag(all(~isnan(d_realimag), 2),:);
        % FINDS MEST 
        % A 
        A_real = Gt_real * G_real;
        % A_realimag = Gt_realimag * G_realimag;
        % B 
        % Real only
        B_real = Gt_real * d_real;
        % Real and imag (complex magnitude) 
        % B_realimag = Gt_realimag * d_realimag;
        % SAVES to allocated space 
        % mEst 
        % Real only 
        mEst_real = A_real\B_real;
        % Real and imag (complex magnitude) 
        % mEst_realimag = A_realimag\B_realimag;
        % Residuals (e and E)
        % Real only 
        e_real = d_real - (G_real * mEst_real);
        % x2 = [x2, e_real]; 
        % res_real = e_real' * e_real;    
        res_real = e_real' * e_real;   
        % PLOT RESIDUALS FOR CERTAIN SCANS - this code only works for real part!
%        if s == 14 || s == 28 || s == 35 || s == 49 % Plots residuals for -39, zenith, 39 and 79 degrees 
%             % CALCULATE RESIDUALS 
%             idxCoiT = idxCoi(idxFt:idxFb, idxL:idxH); % Extracts analysis window from COI template 
%             idxCoiT = idxCoiT(:); 
%             idxCoiT = double(idxCoiT); % Converts to double
%             rIndR = 1; % Index for moving through rows of R
%             for rIndRv = 1:length(idxCoiT) % For each row of idxCoiT
%               if idxCoiT(rIndRv,:) == 1 % If cell is inside COI
%                  idxCoiT(rIndRv,:) = NaN; % Convert to 'NaN' 
%               elseif idxCoiT(rIndRv,:) == 0 % But if cell is outside COI
%                  y = rIndR; % Assigns value 1
%                  rIndR = rIndR + 1; % Adds on to index to move on one place
%                  idxCoiT(rIndRv, :) = e_real(y, :); % Populates template with the corresponding residual value
%               else
%               end
%             end
%            residuals = reshape(idxCoiT, length(F), length(lambda)); % Reshapes compiled residuals back into a matrix for plotting
%            % PLOTS RESIDUALS 
%            if s == 14
%                subplot(2,2,1); hold on 
%            elseif s == 28
%                subplot(2,2,2)
%            elseif s == 35
%                subplot(2,2,3)
%            else 
%                subplot(2,2,4) 
%            end 
%            title(sprintf('%d°', elevAn(:, s))); hold on;     
%            pcolor(lambda, F, residuals); shading interp; colorbar; colormap(jet(1000)); set(gcf,'color','w'); 
%            xlim([min(lambda) max(lambda)]); set(gca,'YScale','log'); ylabel('Spatial frequency (cycles/\lambda)', 'FontSize', 14);
%            set(gcf,'color','w'); ylim([min(F) max(F)]); xlabel('\lambda (nm)', 'FontSize', 14);
%            cb = colorbar; cb.Label.Position(1) = 3.5;
%            ylabel(cb,'Amplitude', 'Rotation', -90, 'FontSize', 10);
%            box on
%            % SAVES FIGURE
%            fnOut = 'ReferenceSpectra_AW3cF_Residual_NOVAC';
%            fname = fullfile(outDirFig, fnOut);
%            saveas(gca, fname)
%         else
%         end
        % Real and imag (complex magnitude) 
        % e_realimag = d_realimag - (G_realimag * mEst_realimag);
        % res_realimag = e_realimag' * e_realimag; 
        % SAVES to preallocated space 
        Results_mEst_real(:,s) = mEst_real; 
        % Results_mEst_realimag(:,s) = mEst_realimag; 
        Results_res_real(:,s) = res_real; 
        % Results_res_realimag(:,s) = res_realimag; 
    end 
end 
% SAVES variables 
fname = fullfile(outDir, 'ReferenceSpectra_AW3cF_mEst_NOVAC_Masaya_1941');
save(fname, 'G_real', 'G_realimag', 'elevAn', 'lambda', 'F', 'Results_mEst_real', 'Results_res_real', 'Results_stdErr_real', 'wavL', 'wavH', 'Ft', 'Fb'); 
toc 

% Plot results to check look OK 
figure
title('Masaya 1941'); hold on;     
plot(elevAn, Results_mEst_real(2,:), 'o'); 
xlabel('Scan elevation angle (° from zenith)'); 
ylabel('mEst_S_O_2 (molec/cm^2)');
xlim([min(elevAn) max(elevAn)])













































