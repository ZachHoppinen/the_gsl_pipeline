function paramStruct = initPROCpar
% INITPROCPAR intitializes the GAMMA software-style PROC_par parameter structure 
% and returns the structure in the output paramter: paramStruct
%
% Usage:
%   paramStruct = initPROCpar
%
% Created: Othmar Frey <frey@gamma-rs.ch>, 03. July 2012
%

paramStruct.title                           = 'SAR';
paramStruct.date                            = [0000 00 00];
paramStruct.raw_data_start_time             = [00 00 0.0000];
paramStruct.channel_mode                    = 'VV';

paramStruct.earth_semi_major_axis           = 6378137.0000;
paramStruct.earth_semi_minor_axis           = 6356752.3141;
paramStruct.scene_center_latitude           = 0.0;
paramStruct.scene_center_longitude          = 0.0;
paramStruct.track_angle                     = 0.0;
paramStruct.platform_altitude               = 0.0;
paramStruct.terrain_height                  = 0.0;
paramStruct.sensor_position_vector(1)       = 0.0;
paramStruct.sensor_position_vector(2)       = 0.0;
paramStruct.sensor_position_vector(3)       = 0.0;
paramStruct.sensor_velocity_vector(1)       = 0.0;
paramStruct.sensor_velocity_vector(2)       = 0.0;
paramStruct.sensor_velocity_vector(3)       = 0.0;
paramStruct.sensor_acceleration_vector(1)   = 0.0;
paramStruct.sensor_acceleration_vector(2)   = 0.0;
paramStruct.sensor_acceleration_vector(3)   = 0.0;
paramStruct.pulse_repetition_frequency      = 0.0;
paramStruct.I_bias                          = 0.0;
paramStruct.Q_bias                          = 0.0;
paramStruct.I_sigma                         = 0.0;
paramStruct.Q_sigma                         = 0.0;
paramStruct.IQ_corr                         = 0.0;

paramStruct.SNR_range_spectrum              = 0.0;
paramStruct.DAR_doppler                     = 0.0;   % unambiguous Doppler estimate
paramStruct.DAR_snr                         = 0.0;   % unambiguous Doppler estimate SNR

paramStruct.doppler_polynomial(1)           = 0.0;
paramStruct.doppler_polynomial(2)           = 0.0;
paramStruct.doppler_polynomial(3)           = 0.0;
paramStruct.doppler_polynomial(4)           = 0.0;

paramStruct.doppler_poly_dot(1)             = 0.0;
paramStruct.doppler_poly_dot(2)             = 0.0;
paramStruct.doppler_poly_dot(3)             = 0.0;
paramStruct.doppler_poly_dot(4)             = 0.0;

paramStruct.doppler_poly_ddot(1)            = 0.0;
paramStruct.doppler_poly_ddot(2)            = 0.0;
paramStruct.doppler_poly_ddot(3)            = 0.0;
paramStruct.doppler_poly_ddot(4)            = 0.0;

paramStruct.echo_time_delay                 = 0.0;

paramStruct.receiver_gain                   = 0.0;
paramStruct.calibration_gain                = 0.0;

paramStruct.near_range_raw                  = 0.0;
paramStruct.center_range_raw                = 0.0;
paramStruct.far_range_raw                   = 0.0;

paramStruct.near_range_slc                  = 0.0;
paramStruct.center_range_slc                = 0.0;
paramStruct.far_range_slc                   = 0.0;
 
paramStruct.range_pixel_spacing             = 0.0;
paramStruct.range_resolution                = 0.0;

paramStruct.sec_range_migration             = 'OFF'; % secondary range migration	                
paramStruct.azimuth_deskew                  = 'ON';  % azimuth deskew
paramStruct.autofocus_snr                   = 0.0;   % autofocus SNR 
paramStruct.azimuth_bandwidth_fraction      = 0.8;   % default fraction of the PRF bandwidth to process
paramStruct.azimuth_presum_factor           = 1;     % default no prefilter
paramStruct.total_raw_echoes                = 0;     % total number of echoes in the raw data file
paramStruct.offset_to_first_echo_to_process = 0;     % start at as close to the beginning of the raw data as possible
paramStruct.echoes_to_process               = 0;
paramStruct.range_offset                    = 0;
paramStruct.raw_range_samples               = 0;
paramStruct.near_range_extension            = 0;
paramStruct.far_range_extension             = 0;
paramStruct.pre_azimuth_extension           = 0;
paramStruct.post_azimuth_extension          = 0;
paramStruct.range_looks                     = 1;
paramStruct.azimuth_looks                   = 1;
paramStruct.azimuth_offset                  = 0;
paramStruct.azimuth_pixel_spacing           = 0.0;
paramStruct.azimuth_resolution              = 0.0;
paramStruct.range_pixels                    = 0;
paramStruct.azimuth_pixels                  = 0;
paramStruct.image_format                    = 'SCOMPLEX';
paramStruct.sensor_latitude                 = 0.0;
paramStruct.sensor_longitude                = 0.0;
paramStruct.sensor_track_angle              = 0.0;

paramStruct.map_coordinate_1(1)            = 0.0;
paramStruct.map_coordinate_1(2)            = 0.0;
paramStruct.map_coordinate_1(3)            = 0.0;

paramStruct.map_coordinate_2(1)            = 0.0;
paramStruct.map_coordinate_2(2)            = 0.0;
paramStruct.map_coordinate_2(3)            = 0.0;

paramStruct.map_coordinate_3(1)            = 0.0;
paramStruct.map_coordinate_3(2)            = 0.0;
paramStruct.map_coordinate_3(3)            = 0.0;

paramStruct.map_coordinate_4(1)            = 0.0;
paramStruct.map_coordinate_4(2)            = 0.0;
paramStruct.map_coordinate_4(3)            = 0.0;

paramStruct.map_coordinate_5(1)            = 0.0;
paramStruct.map_coordinate_5(2)            = 0.0;
paramStruct.map_coordinate_5(3)            = 0.0;

paramStruct.number_of_state_vectors         = 1;
  
paramStruct.time_of_first_state_vector      = 0.0;
paramStruct.state_vector_interval           = 0.0;

paramStruct.state_vector_position(1,1)      = 0.0;
paramStruct.state_vector_position(1,2)      = 0.0;
paramStruct.state_vector_position(1,3)      = 0.0;

paramStruct.state_vector_velocity(1,1)      = 0.0;
paramStruct.state_vector_velocity(1,2)      = 0.0;
paramStruct.state_vector_velocity(1,3)      = 0.0;

