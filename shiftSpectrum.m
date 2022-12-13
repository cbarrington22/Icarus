% Icarus: exploIting the frequenCy signature of trace gAs absoRption in Uv-Spectra
%  
% Appendix D: Calculation of Bshift 
%
% This function is provided as part of Chapter 2 of the thesis 'Using Volcanic Gases to Understand Open-vent Volcanoes' submitted to the Nanyang Technological University
% in partial fulfillment of the requirements for the degree of Doctor of Philosophy. 
%
% Written by C. Barrington December 2022
%
% DESCRIPTION of function: 
% 
% This function calculates Bshift spectrum from convoluted I0
% 
% INSTRUCTIONS for user: 
% 
% (1) Use function with the following INPUTS: 
% I0 - is intenisty information of the (resampled) and convoluted solar reference. It must be vector where the number of rows is equal to the number of channels of the spectrometer.
%
% (2) The function OUTPUTS: 
% Bshift - is the Bshift spectrum 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOT INTENDED TO BE MODIFIED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Bshift] = shiftSpectrum(I0)
        outyP = zeros(length(I0), 1); 
        h = 1; % Channels 
        y = I0; 
        % y prime
        for yInd = 3:length(y)
            yPrime = (1/(2*h))*(y(yInd-2)-8*y(yInd-1)+8*y(yInd)-(yInd+1));
            outyP(yInd, :) = yPrime;
        end
        % Bshift 
        yPrime = outyP;
        Bshift = yPrime./y; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%