function complexCoherence = calcComplexCoherence(slc1,slc2,filter2D)
% CALCCOMPLEXCOHERENCE calculates the complex coherence between two SLC
%   data sets. I.e., it can be used to calculate the complex interfero-
%   metric coherence between two SLCs as well as the complex polarimetric
%   coherence between polarimetric channels. 
%   As a third input parameter a 2-d filter window is required to apply
%   multilooking by means of spatial averaging.
%
%     Usage:
%      complexCoherence = calcComplexCoherence(slc1,slc2,filter2D);
%
%     Input:     
%       slc1             :   reference SLC
%       slc2             :   co-registered SLC
%       filter2D         :   2-D filter window
%
%           Examples for creating filter windows:
%             filter2D = 1/36*ones(6);
%
%            or:
%             fsize = 25
%             sigma = round(fsize/2);
%             filter2D = fspecial('gaussian', fsize, sigma);
%
%     Output:
%       complexCoherence :   2-D array of complex coherences
%
%
%     Last modified: 08.Sep.2017
%     Copyright: 2017 Gamma Remote Sensing AG
%                Othmar Frey, frey@gamma-rs.ch
%
  
complexCoherence = slc1.*conj(slc2);
complexCoherence =  (filter2(filter2D,complexCoherence))./(sqrt( filter2(filter2D,(abs(slc1)).^2) .* filter2(filter2D,(abs(slc2)).^2)));

end
