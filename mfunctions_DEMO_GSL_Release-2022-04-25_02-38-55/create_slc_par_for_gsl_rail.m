function create_slc_par_for_gsl_rail(strucPolGrid,strucLocation,RAWparFilename,PROCparFilename,SARparFilename,outSLCparFilename)
% create_slc_par_for_gsl_rail writes a GAMMA-style SLC parameter file
%  (*.slc.par) for GAMMA L-band SAR data acquired on the GAMMA linear rail.
%   
%  It requires  
%   
%
% in addition to the GAMMA 
%
%   Usage:
%       create_slc_par_for_gsl_rail(strucPolGrid,strucLocation,RAWparFilename,PROCparFilename,SARparFilename,outSLCparFilename)
%
%       where:
%           strucPolGrid    : Structure variable (output variable of create_polar_rec_grid_for_gsl_rail)
%                               with fields
%                                   polar_angles        array with polar angles of polar
%                                                        reconstruction grid [deg]
%                                   ang_samp            angular sampling of polar grid [deg]
%                                   az_angular_spread   angular total spread of polar grid [deg]
%                                   slr                 slant range distances (array) [m]
%                                   rps                 range pixel spacing [m]
%                                   pn1                 range pixel numbers (starting with index 0)
%                            
%
%
%   SEE ALSO:
%       create_polar_rec_grid_for_gsl_rail.m, create_pos_vel_for_rail.m 
%
%
%   Created: by Othmar Frey <frey@gamma-rs.ch>, 30. Sep. 2019
%
%   Copyright: 2019 Gamma Remote Sensing AG
%              Othmar Frey <frey@gamma-rs.ch>
%
%
%

% Read raw_par, PROC_par and SAR_par file and initialize SLC_par structure
rawPar = readGammaParFile(RAWparFilename);
procPar = readGammaParFile(PROCparFilename);
sarPar = readGammaParFile(SARparFilename);
slcPar = initSLCGPRIpar;

%keyboard
% Fill in the values into slcPar structure:
slcPar.title = [procPar.title ' ' procPar.channel_mode];
slcPar.sensor = sarPar.sensor_name;
slcPar.date = procPar.date;

slcPar.range_samples = procPar.range_pixels;
slcPar.azimuth_lines = length(strucPolGrid.polar_angles);

%slcPar.start_time = [num2str(procPar.raw_data_start_time(1)) ' ' num2str(procPar.raw_data_start_time(2)) ' ' num2str(procPar.raw_data_start_time(3))]
slcPar.start_time = 0.0;
slcPar.end_time =(procPar.total_raw_echoes-1)/procPar.pulse_repetition_frequency;
slcPar.center_time = (slcPar.start_time + slcPar.end_time)/2.0;

slcPar.azimuth_line_time = (slcPar.end_time-slcPar.start_time)/(slcPar.azimuth_lines-1);
%                  line_header_size: 0
slcPar.range_looks = 1;
slcPar.azimuth_looks = 1;
slcPar.image_format = procPar.image_format;
slcPar.image_geometry = 'SLANT_RANGE';
slcPar.range_scale_factor = 1.0;
slcPar.azimuth_scale_factor = 1.0;
slcPar.range_pixel_spacing = procPar.range_pixel_spacing; 
    %slcPar.azimuth_pixel_spacing = 0.0; 
slcPar.near_range_slc = procPar.near_range_slc;
slcPar.center_range_slc = procPar.center_range_slc;
slcPar.far_range_slc = procPar.far_range_slc;

slcPar.azimuth_deskew = 'ON'; % ON: radar data is in zero Doppler geometry; OFF: not in z.D. geometry
slcPar.azimuth_angle = 90.0000; % decides whether radar is right-looking (90) of left-looking (-90)    degrees

slcPar.radar_frequency = sarPar.SAR_center_frequency;
slcPar.adc_sampling_rate     = sarPar.ADC_sampling_frequency;
slcPar.chirp_bandwidth       = sarPar.chirp_bandwidth;
slcPar.prf                   = procPar.pulse_repetition_frequency;
    %receiver_gain               = 0.0
slcPar.earth_semi_major_axis = procPar.earth_semi_major_axis;
slcPar.earth_semi_minor_axis = procPar.earth_semi_minor_axis;
slcPar.number_of_state_vectors = 0;
    %GPRI_TX_mode: None
slcPar.GPRI_TX_antenna = procPar.channel_mode;
slcPar.GPRI_az_start_angle = strucPolGrid.polar_angles(1); %  degrees
slcPar.GPRI_az_angle_step = strucPolGrid.ang_samp;
slcPar.GPRI_ant_elev_angle = 0.0;       %  degrees

%if external information on the positioning is availble 
% use that
% otherwise use the coordinates taken from the raw_par file.
if(isfield(strucLocation,'GPRI_ref_north'))
    slcPar.GPRI_ref_north = strucLocation.GPRI_ref_north;
else
    slcPar.GPRI_ref_north = rawPar.geographic_coordinates(1);
end
if(isfield(strucLocation,'GPRI_ref_east'))
    slcPar.GPRI_ref_east = strucLocation.GPRI_ref_east;
else
    slcPar.GPRI_ref_east = rawPar.geographic_coordinates(2);
end
if(isfield(strucLocation,'GPRI_ref_alt'))
    slcPar.GPRI_ref_alt = strucLocation.GPRI_ref_alt;
else
    slcPar.GPRI_ref_alt = rawPar.geographic_coordinates(3);
end
if(isfield(strucLocation,'GPRI_geoid'))
    slcPar.GPRI_geoid = strucLocation.GPRI_geoid;
else
    slcPar.GPRI_geoid = rawPar.geographic_coordinates(4);% m
end
% The heading (pointing direction in degree,
% clockwise with 0 pointing towards north)
slcPar.GPRI_scan_heading = strucLocation.GPRI_scan_heading; % deg
    %GPRI_tx_coord:     0.00000    0.00000   0.00000  m m m
    %GPRI_rx1_coord:    0.00000    0.00000   0.00000  m m m
    %GPRI_rx2_coord:    0.00000    0.00000   0.00000  m m m
    %GPRI_tower_roll:   0.00000  degrees
    %GPRI_tower_pitch:  0.00000  degrees
    %GPRI_phase_offset: 0.00000  radians
% Write SLC parameter file to disk
writeSLCpar(outSLCparFilename,slcPar);
end