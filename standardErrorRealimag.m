% Error calculation 
% realimag 

function [mEst_realimag] = standardErrorRealimag(dr, G1r, G2r, G3r, G4r, G5r, di, G1i, G2i, G3i, G4i, G5i)

% dr - WT_J_real 

% G1r - WT_I0_real

% G2r - WT_SO2_real

% G3r - WT_O3_real

% G4r - WT_Ring_real

% G5r - WT_Bshift_real

% di - WT_J_imag 

% G1i - WT_I0_imag

% G2i - WT_SO2_imag

% G3i - WT_O3_imag

% G4i - WT_Ring_imag

% G5i - WT_Bshift_imag

% Once wavelengths have been randomly sampled, calculate TRANSPOSE to return d and G 
drt = dr'; 
G1rt = G1r'; 
G2rt = G2r'; 
G3rt = G3r'; 
G4rt = G4r'; 
G5rt = G5r';
dit = di'; 
G1it = G1i'; 
G2it = G2i'; 
G3it = G3i'; 
G4it = G4i'; 
G5it = G5i';

% Determines, d  
% Real
drt_v = drt(:); 
% Imag
dit_v = dit(:); 
% Combines real and imaginary part for complex magnitude 
d_realimag = [drt_v; dit_v];
d_realimag = d_realimag(all(~isnan(d_realimag), 2),:);

% Re-builds G
% Real
G1rt_v = G1rt(:); 
G2rt_v = G2rt(:); 
G3rt_v = G3rt(:); 
G4rt_v = G4rt(:); 
G5rt_v = G5rt(:); 
% Imag
G1it_v = G1it(:); 
G2it_v = G2it(:); 
G3it_v = G3it(:); 
G4it_v = G4it(:); 
G5it_v = G5it(:); 
% Combines real and imaginary part for complex magnitude 
G1rit_v = [G1rt_v; G1it_v];
G2rit_v = [G2rt_v; G2it_v];
G3rit_v = [G3rt_v; G3it_v];
G4rit_v = [G4rt_v; G4it_v];
G5rit_v = [G5rt_v; G5it_v];
% G
G_realimag = [G1rit_v, G2rit_v, G3rit_v, G4rit_v, G5rit_v]; 
G_realimag = G_realimag(all(~isnan(G_realimag), 2),:);
Gt_realimag = G_realimag'; 

% Finds mEst 
A_realimag = Gt_realimag * G_realimag;
B_realimag = Gt_realimag * d_realimag;
mEst_realimag = A_realimag\B_realimag;
end 
