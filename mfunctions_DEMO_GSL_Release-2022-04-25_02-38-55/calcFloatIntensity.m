function floatIntensity = calcFloatIntensity(slc,filter2D)
% calcFloatIntensity calculates the intenisity of a complex SLC
%   As a second input parameter a 2-d filter window is required to apply
%   multilooking by means of spatial averaging.
%
%     Usage:
%      floatIntensity = calcFloatIntensity(slc1,filter2D);
%
%     Input:     
%       slc              :   SLC image
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
%       floatIntensity :   2-D array of float intensity values
%
%
%     Last modified: 20-Sep-2019
%     Copyright: 2017 Gamma Remote Sensing AG
%                Othmar Frey, frey@gamma-rs.ch 
%                Charles Werner, cw@gamma-rs.ch
%

floatIntensity = real(slc).*real(slc) + imag(slc).*imag(slc);
floatIntensity = filter2(filter2D,floatIntensity);

end
