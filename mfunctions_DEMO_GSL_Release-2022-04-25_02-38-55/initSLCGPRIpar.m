function paramStruct = initSLCGPRIpar
% INITSLCPAR initializes a matlab structure variable
%   containing the GAMMA-style SLC par keywords as fields
%   including the additional keywords used for the
%   terrestrial radars (GPRI, rail-based SAR)
%
% Usage:
%   paramStruct = initSLCGPRIpar
%
%
%  Created: otfrey@ethz.ch, 26. Sep 2019
%

paramStruct.title = ' ';
paramStruct.sensor = ' ';
paramStruct.date   = [0000 00 00];
paramStruct.start_time =                0.0;  % s
paramStruct.center_time =               0.0;  % s
paramStruct.end_time  =                 0.0;  % s
paramStruct.azimuth_line_time =         0.0;  % s
paramStruct.line_header_size =          0;
paramStruct.range_samples =             0;
paramStruct.azimuth_lines =             0;
paramStruct.range_looks =               0;
paramStruct.azimuth_looks =             0;
paramStruct.image_format =            ' ';
paramStruct.image_geometry =          ' ';
paramStruct.range_scale_factor =        0.0000000e+00;
paramStruct.azimuth_scale_factor =      0.0000000e+00;
paramStruct.center_latitude =           0.0;   %degrees
paramStruct.center_longitude =          0.0;   %degrees
paramStruct.heading =                   0.0;   %degrees
paramStruct.range_pixel_spacing   =     0.0;   %m
paramStruct.azimuth_pixel_spacing =     0.0;   %m
paramStruct.near_range_slc =            0.0;   % m
paramStruct.center_range_slc =          0.0;   % m
paramStruct.far_range_slc =             0.0;   % m
paramStruct.first_slant_range_polynomial   = [0.00000 0.00000  0.00000e+00  0.00000e+00  0.00000e+00  0.00000e+00];
paramStruct.center_slant_range_polynomial  = [0.00000 0.00000  0.00000e+00  0.00000e+00  0.00000e+00  0.00000e+00];
paramStruct.last_slant_range_polynomial    = [0.00000 0.00000  0.00000e+00  0.00000e+00  0.00000e+00  0.00000e+00];
paramStruct.incidence_angle =           0.0;  % degrees
paramStruct.azimuth_deskew =  'ON';
paramStruct.azimuth_angle  =            0.0;  %   degrees
paramStruct.radar_frequency =           0.0;  %   Hz
paramStruct.adc_sampling_rate =         0.0;  %   Hz
paramStruct.chirp_bandwidth =           0.0;  %   Hz
paramStruct.prf =                       0.0;  %   Hz
paramStruct.azimuth_proc_bandwidth =    0.0;  %   Hz

paramStruct.doppler_polynomial = [0.00000  0.00000e+00  0.00000e+00  0.00000e+00];
paramStruct.doppler_poly_dot =   [0.00000e+00  0.00000e+00  0.00000e+00  0.00000e+00];
paramStruct.doppler_poly_ddot =  [0.00000e+00  0.00000e+00  0.00000e+00  0.00000e+00];
paramStruct.receiver_gain = 0.0000;
paramStruct.calibration_gain = 0.0000;
paramStruct.sar_to_earth_center = 0.0;
paramStruct.earth_radius_below_sensor = 0.0;
paramStruct.earth_semi_major_axis = 0.0;
paramStruct.earth_semi_minor_axis = 0.0;
paramStruct.number_of_state_vectors =  0;
%paramStruct.time_of_first_state_vector = 0.0;
%paramStruct.state_vector_interval = 0.0;
%MAX_STATE = 64;
%for i = 1:MAX_STATE
%  paramStruct.state_vector_position(i,:) = [0.0   0.0  0.0];
%  paramStruct.state_vector_velocity(i,:) = [0.0   0.0  0.0];
%end
% GPRI data structure 20080825 clw
paramStruct.GPRI_TX_mode = ' ';         % transmitter mode  H, V, HV, CTI, or None */
paramStruct.GPRI_TX_antenna = ' ';      % transmitter antenna used for the current image (H, V, Upper, Lower) */
paramStruct.GPRI_az_start_angle = 0.0;  % external info /* starting rotation angle (degrees, + is clockwise looking down the rotation axis) */
paramStruct.GPRI_az_angle_step = 0.0;   % external from pre_proc m file negative sign for rail  /* angular step degrees (+angle is clockwise rotation looking down the rotation axis) */
paramStruct.GPRI_ant_elev_angle = 0.0;  % antenna elevation angle, + angles are up,*/
paramStruct.GPRI_ref_north = 0.0;       % reference point northing or latitude in the projection and datum of the DEM */
paramStruct.GPRI_ref_east = 0.0;        % from raw_par /* reference point easting or longitude in the projection and datum of the DEM */
paramStruct.GPRI_ref_alt = 0.0;         % from raw_par /* reference point altitude in the projection and datum of the DEM */
paramStruct.GPRI_geoid = 0.0;           % height of the geoid relative to the WGS84 ellipsoid */
paramStruct.GPRI_scan_heading = 0.0;    % heading of central sweep scan line relative to north, + is clockwise from north
                                   % looking down the tower rotation axis. The tower rotation axis is aligned with vertical down
                                   % towards the earth center. Down is one component of the North, East, Down (NED) coordinate system.
paramStruct.GPRI_tx_coord = [0.0 0.0 0.0];       % VEC:  transmit antenna phase center coordinates (XYZ meters) in the local tower coordinate system */
paramStruct.GPRI_rx1_coord = [0.0 0.0 0.0];      % VEC 0 /* receive antenna 1 phase center coordinates (XYZ meters) in the local tower coordinate system */
paramStruct.GPRI_rx2_coord = [0.0 0.0 0.0];      % VEC /* receive antenna 2 phase center coordinates (XYZ meters) in the local tower coordinate system */
paramStruct.GPRI_tower_roll = 0.0;               % tower roll angle rotation about the local Y axis */
paramStruct.GPRI_tower_pitch = 0.0;              % tower pitch angle rotation about the local X axis */
paramStruct.GPRI_phase_offset = 0.0;             % interferogram phase offset s.t. ph = -4pi/lam* (r2 - r1) + phase_offset */
% GPRI_PAR;




