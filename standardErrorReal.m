% Error calculation 
% real 

function [mEst_real] = standardErrorReal(dr, G1r, G2r, G3r, G4r, G5r)

% dr - WT_J_real 

% G1r - WT_I0_real

% G2r - WT_SO2_real

% G3r - WT_O3_real

% G4r - WT_Ring_real

% G5r - WT_Bshift_real

% Once wavelengths have been randomly sampled, calculate TRANSPOSE to return d and G 
drt = dr'; 
G1rt = G1r'; 
G2rt = G2r'; 
G3rt = G3r'; 
G4rt = G4r'; 
G5rt = G5r';
 
% Determines, d  
drt_v = drt(:); 
d_real = drt_v;
d_real = d_real(all(~isnan(d_real), 2),:);

% Re-builds G
G1rt_v = G1rt(:); 
G2rt_v = G2rt(:); 
G3rt_v = G3rt(:); 
G4rt_v = G4rt(:); 
G5rt_v = G5rt(:); 
G_real = [G1rt_v, G2rt_v, G3rt_v, G4rt_v, G5rt_v]; 
G_real = G_real(all(~isnan(G_real), 2),:);
Gt_real = G_real'; 

% Finds mEst 
A_real = Gt_real * G_real;
B_real = Gt_real * d_real;
mEst_real = A_real\B_real;
end 
