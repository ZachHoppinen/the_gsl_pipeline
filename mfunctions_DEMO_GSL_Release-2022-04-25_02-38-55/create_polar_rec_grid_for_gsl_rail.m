function strucPolGrid = create_polar_rec_grid_for_gsl_rail(inRawParFilename,az_3dB_bw,outPolarRecGridFilename,az_ovr_samp)
% strucPolGrid = create_polar_rec_grid_for_gsl_rail creates a polar reconstruction grid
%   for rail-based GAMMA L-band SAR data focusing onto a polar reconstruction grid.
%   The polar reconstruction grid can then be used with the TDBP processor of the 
%   GAMMA software.
%
%   USAGE:
%       strucPolGrid = create_polar_rec_grid_for_gsl_rail(inRawParFilename,az_3dB_bw,outPolarRecGridFilename,az_ovr_samp)
% 
%
%   SEE ALSO:
%       create_slc_par_for_gsl_rail.m, create_pos_vel_for_rail.m 
%
%
%   Created:         2018-07-17 by Othmar Frey <frey@gamma-rs.ch>
%   Modified:        2019-09-26 by Othmar Frey <frey@gamma-rs.ch> 
%                           - Angular spread of the polar grid increased from 3dB beamwidth 
%                            to 2 times the 3dB beamwidth 
%                           - added output parameter strucPolGrid
%                             (used for SLC_par information)
%
%   Modifed:        2019-10-11 by 
%                   Changed:
%                       rps = (rawParStruct.ADC_sample_rate/rawParStruct.CGEN_num_samp*c/2.)/rawParStruct.RF_chirp_rate; % range pixel spacing
%                   to
%                       rps = (rawParStruct.ADC_sample_rate/rawParStruct.CGEN_num_samp*c/2.)/abs(rawParStruct.RF_chirp_rate); % range pixel spacing
%                   to mitigate the issue caused by the changed of sign of chirp rate of the GAMMA L-band SAR (raw_par file).
%                   Also changed:
%                    To convert to a GPRI compatible clockwise (left-to-right) image 
%                       the sign of the angles is changed (new: -polar_angles ) for 
%                       x = slr*cos(-polar_angles);
%                       y = slr*sin(-polar_angles);
%                   
%   Copyright:  2019 Gamma Remote Sensing AG
%               Othmar Frey <frey@gamma-rs.ch>
%


c = 299792458.0;

verbose = 0;

% Read GS-L raw parameter file
rawParStruct = readGammaParFile(inRawParFilename);

start_pos = rawParStruct.RAIL_start_pos;
end_pos = rawParStruct.RAIL_end_pos; 
vel = rawParStruct.RAIL_velocity; 
acc = rawParStruct.RAIL_acceleration; 
%    else if start_pos == None or end_pos == None or vel == None or acc == None:
%    raise SystemExit('ERROR: RAIL_motion_params: missing position, velocity, or acceleration parameters!')

SA_mm = end_pos - start_pos;
SA = SA_mm/1000.0;

lambda = c/rawParStruct.RF_center_freq;

ang_res = lambda/(2*SA);

ang_samp = ang_res*az_ovr_samp;

% The total angular spread of the polar 
% reconstruction grid is now
% now 2 times the 3dB beamwidth
az_agular_spread_rad = 2.0*az_3dB_bw*pi/180;

polar_angles_rad_left = -ang_samp:-ang_samp:-(az_agular_spread_rad/2);
polar_angles_rad_right = 0:ang_samp:(az_agular_spread_rad/2);

polar_angles = [fliplr(polar_angles_rad_left) polar_angles_rad_right];

if(verbose)
    figure, plot(polar_angles)
end

RANGE_OFFSET = 0;
% list of slant range pixel numbers
pn1 = 0:(floor(rawParStruct.CGEN_num_samp/2)-1); 
% range pixel spacing
rps = (rawParStruct.ADC_sample_rate/rawParStruct.CGEN_num_samp*c/2.)/abs(rawParStruct.RF_chirp_rate); % range pixel spacing
%slant range for each sample
slr = (pn1.' *rps) + RANGE_OFFSET; % array of slant range distances

% To convert to a GPRI compatible clockwise (left-to-right) image 
% the sign of the angles is changed here
x = slr*cos(-polar_angles);
y = slr*sin(-polar_angles);
z = zeros(size(x));

if(verbose)
    figure, plot(x(1:100:end,:),y(1:100:end,:),'LineStyle','none','Marker','o')
    axis equal
end

writeMatrixNoHeader(outPolarRecGridFilename,[x(:); y(:); z(:)],'double');

strucPolGrid.polar_angles = polar_angles*180.0/pi;
strucPolGrid.ang_samp = ang_samp*180.0/pi;
strucPolGrid.az_agular_spread = az_agular_spread_rad*180.0/pi;
strucPolGrid.slr = slr;
strucPolGrid.rps = rps;
strucPolGrid.pn1 = pn1;
