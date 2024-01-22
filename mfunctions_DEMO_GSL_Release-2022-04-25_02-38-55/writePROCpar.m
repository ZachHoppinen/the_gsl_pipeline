function writePROCpar(procParFileName,paramStruct)
% WRITEPROCPAR writes the the GAMMA software-style PROC_par parameter
% structure stored in the structure variable paramStruct to the
% processing parameter file with the name given in procParFileName.
%
% Usage:
%   writePROCpar(procParFileName,paramStruct)
%
% SEE ALSO:
%   initPROCpar.m, and further related: initSARpar.m, writeSARpar.m, initSLCpar.m, writeSLCpar.m
%
% Created:  Othmar Frey <frey@gamma-rs.ch>, 03. July 2012
% Modified: Othmar Frey <frey@gamma-rs.ch>, 17. Sep. 2019
%               if (paramStruct.number_of_state_vectors == 0)
%               --> fprintf(fid,'\n');
%               --> fclose(fid);
%                   return;
%               end
%

fid=fopen(procParFileName,'wt');
message = ferror(fid);
if(~isempty(message))
	error(message);
	return;
end

fprintf(fid, 'GAMMA Modular SAR Processor (MSP) - Processing Parameter File v2.5\n\n');
fprintf(fid, 'title:                  %s\n', paramStruct.title);
fprintf(fid, 'date:                   %d %d %d\n', paramStruct.date(1),paramStruct.date(2),paramStruct.date(3));
fprintf(fid, 'raw_data_start_time:    %.0f %.0f %.12f\n', paramStruct.raw_data_start_time(1),paramStruct.raw_data_start_time(2),paramStruct.raw_data_start_time(3));
fprintf(fid, 'channel/mode:           %s\n', paramStruct.channel_mode);

fprintf(fid, 'earth_semi_major_axis:        %12.4f   m\n', paramStruct.earth_semi_major_axis);
fprintf(fid, 'earth_semi_minor_axis:        %12.4f   m\n', paramStruct.earth_semi_minor_axis);
fprintf(fid, 'scene_center_latitude:        %12.6f   decimal degrees\n', paramStruct.scene_center_latitude);
fprintf(fid, 'scene_center_longitude:       %12.6f   decimal degrees\n', paramStruct.scene_center_longitude);
fprintf(fid, 'track_angle:                  %12.6f   degrees\n', paramStruct.track_angle);
fprintf(fid, 'platform_altitude:            %12.4f   m\n', paramStruct.platform_altitude);
fprintf(fid, 'terrain_height:               %12.4f   m\n', paramStruct.terrain_height);
fprintf(fid, 'sensor_position_vector:     %14.6f %14.6f %14.6f   m     m      m\n',    paramStruct.sensor_position_vector(1),    paramStruct.sensor_position_vector(2),    paramStruct.sensor_position_vector(3));
fprintf(fid, 'sensor_velocity_vector:     %14.6f %14.6f %14.6f   m/s   m/s   m/s  \n', paramStruct.sensor_velocity_vector(1),    paramStruct.sensor_velocity_vector(2),    paramStruct.sensor_velocity_vector(3));
fprintf(fid, 'sensor_acceleration_vector: %14.6f %14.6f %14.6f   m/s^2 m/s^2 m/s^2\n', paramStruct.sensor_acceleration_vector(1),paramStruct.sensor_acceleration_vector(2),paramStruct.sensor_acceleration_vector(3));
fprintf(fid, 'pulse_repetition_frequency: %14.6f   Hz\n', paramStruct.pulse_repetition_frequency);
fprintf(fid, 'I_bias:                     %14.6f\n', paramStruct.I_bias);
fprintf(fid, 'Q_bias:                     %14.6f\n', paramStruct.Q_bias);
fprintf(fid, 'I_sigma:                    %14.6f\n', paramStruct.I_sigma);
fprintf(fid, 'Q_sigma:                    %14.6f\n', paramStruct.Q_sigma);
fprintf(fid, 'IQ_corr:                    %14.6f\n', paramStruct.IQ_corr);

fprintf(fid, 'SNR_range_spectrum:         %14.3f\n', paramStruct.SNR_range_spectrum);
fprintf(fid, 'DAR_doppler:                %14.3f   Hz\n', paramStruct.DAR_doppler);
fprintf(fid, 'DAR_snr:                    %14.3f\n', paramStruct.DAR_snr);

fprintf(fid, 'doppler_polynomial:  %12.5e %12.5e %12.5e %12.5e   Hz     Hz/m     Hz/m^2     Hz/m^3\n',    paramStruct.doppler_polynomial(1),paramStruct.doppler_polynomial(2),paramStruct.doppler_polynomial(3),paramStruct.doppler_polynomial(4));
fprintf(fid, 'doppler_poly_dot:    %12.5e %12.5e %12.5e %12.5e   Hz/s   Hz/s/m   Hz/s/m^2   Hz/s/m^3\n',  paramStruct.doppler_poly_dot(1),  paramStruct.doppler_poly_dot(2),  paramStruct.doppler_poly_dot(3),  paramStruct.doppler_poly_dot(4));
fprintf(fid, 'doppler_poly_ddot:   %12.5e %12.5e %12.5e %12.5e   Hz/s^2 Hz/s^2/m Hz/s^2/m^2 Hz/s^2/m^3\n',paramStruct.doppler_poly_ddot(1), paramStruct.doppler_poly_ddot(2), paramStruct.doppler_poly_ddot(3), paramStruct.doppler_poly_ddot(4));

fprintf(fid, 'echo_time_delay:            %14.6e   s\n', paramStruct.echo_time_delay);

fprintf(fid, 'receiver_gain:              %14.4f   dB\n', paramStruct.receiver_gain);
fprintf(fid, 'calibration_gain:           %14.4f   dB\n', paramStruct.calibration_gain);

fprintf(fid, 'near_range_raw:             %.9f   m\n', paramStruct.near_range_raw);
fprintf(fid, 'center_range_raw:           %.9f   m\n', paramStruct.center_range_raw);
fprintf(fid, 'far_range_raw:              %.9f   m\n', paramStruct.far_range_raw);

fprintf(fid, 'near_range_slc:             %.9f   m\n', paramStruct.near_range_slc);
fprintf(fid, 'center_range_slc:           %.9f   m\n', paramStruct.center_range_slc);
fprintf(fid, 'far_range_slc:              %.9f   m\n', paramStruct.far_range_slc);

fprintf(fid, 'range_pixel_spacing:          %15.12f   m\n', paramStruct.range_pixel_spacing);
fprintf(fid, 'range_resolution:             %15.12f   m\n', paramStruct.range_resolution);

fprintf(fid, 'sec_range_migration:                 %s\n', paramStruct.sec_range_migration);
fprintf(fid, 'azimuth_deskew:                      %s\n', paramStruct.azimuth_deskew);
fprintf(fid, 'autofocus_snr:                    %8.4f\n', paramStruct.autofocus_snr);
fprintf(fid, 'azimuth_bandwidth_fraction:       %8.4f\n', paramStruct.azimuth_bandwidth_fraction);
fprintf(fid, 'azimuth_presum_factor:              %6d\n', paramStruct.azimuth_presum_factor);
fprintf(fid, 'total_raw_echoes:                   %6d   echoes\n', paramStruct.total_raw_echoes);
fprintf(fid, 'offset_to_first_echo_to_process:    %6d   echoes\n', paramStruct.offset_to_first_echo_to_process);
fprintf(fid, 'echoes_to_process:                  %6d   echoes\n', paramStruct.echoes_to_process);
fprintf(fid, 'range_offset:                       %6d   samples\n', paramStruct.range_offset);
fprintf(fid, 'raw_range_samples:                  %6d   samples\n', paramStruct.raw_range_samples);
fprintf(fid, 'near_range_extension:               %6d   samples\n', paramStruct.near_range_extension);
fprintf(fid, 'far_range_extension:                %6d   samples\n', paramStruct.far_range_extension);
fprintf(fid, 'pre_azimuth_extension:              %6d   echoes\n', paramStruct.pre_azimuth_extension);
fprintf(fid, 'post_azimuth_extension:             %6d   echoes\n', paramStruct.post_azimuth_extension);
fprintf(fid, 'range_looks:                        %6d   looks\n', paramStruct.range_looks);
fprintf(fid, 'azimuth_looks:                      %6d   looks\n', paramStruct.azimuth_looks);
fprintf(fid, 'azimuth_offset:               %12.6f   s\n', paramStruct.azimuth_offset);
fprintf(fid, 'azimuth_pixel_spacing:        %12.6f   m\n', paramStruct.azimuth_pixel_spacing);
fprintf(fid, 'azimuth_resolution:           %12.3f   m\n', paramStruct.azimuth_resolution);
fprintf(fid, 'range_pixels:                       %6d   image output samples\n', paramStruct.range_pixels);
fprintf(fid, 'azimuth_pixels:                     %6d   image output lines\n', paramStruct.azimuth_pixels);
fprintf(fid, 'image_format:                     %s\n', paramStruct.image_format);
fprintf(fid, 'sensor_latitude:              %12.6f   decimal degrees\n', paramStruct.sensor_latitude);
fprintf(fid, 'sensor_longitude:             %12.6f   decimal degrees\n', paramStruct.sensor_longitude);
fprintf(fid, 'sensor_track_angle:           %12.6f   decimal degrees\n', paramStruct.sensor_track_angle);
fprintf(fid, 'map_coordinate_1:             %12.7f   %12.7f   %12.4f   deg.  deg.  m\n', paramStruct.map_coordinate_1(1), paramStruct.map_coordinate_1(2), paramStruct.map_coordinate_1(3));
fprintf(fid, 'map_coordinate_2:             %12.7f   %12.7f   %12.4f   deg.  deg.  m\n', paramStruct.map_coordinate_2(1), paramStruct.map_coordinate_2(2), paramStruct.map_coordinate_2(3));
fprintf(fid, 'map_coordinate_3:             %12.7f   %12.7f   %12.4f   deg.  deg.  m\n', paramStruct.map_coordinate_3(1), paramStruct.map_coordinate_3(2), paramStruct.map_coordinate_3(3));
fprintf(fid, 'map_coordinate_4:             %12.7f   %12.7f   %12.4f   deg.  deg.  m\n', paramStruct.map_coordinate_4(1), paramStruct.map_coordinate_4(2), paramStruct.map_coordinate_4(3));
fprintf(fid, 'map_coordinate_5:             %12.7f   %12.7f   %12.4f   deg.  deg.  m\n', paramStruct.map_coordinate_5(1), paramStruct.map_coordinate_5(2), paramStruct.map_coordinate_5(3));

fprintf(fid, 'number_of_state_vectors:            %6d\n', paramStruct.number_of_state_vectors);

if (paramStruct.number_of_state_vectors == 0)
    fprintf(fid,'\n');
    fclose(fid);
    return;
end
  
fprintf(fid, 'time_of_first_state_vector:   %12.5f   s\n', paramStruct.time_of_first_state_vector);
fprintf(fid, 'state_vector_interval:        %12.5f   s\n', paramStruct.state_vector_interval);

[firstDim, secDim] = size(paramStruct.state_vector_position);
for i=1:firstDim
    fprintf(fid,'%s: %14.4f  %14.4f  %14.4f   m   m   m\n',['state_vector_position_' num2str(i)],paramStruct.state_vector_position(i,1),paramStruct.state_vector_position(i,2),paramStruct.state_vector_position(i,3));
    fprintf(fid,'%s: %14.5f  %14.5f  %14.5f   m/s m/s m/s\n',['state_vector_velocity_' num2str(i)],paramStruct.state_vector_velocity(i,1),paramStruct.state_vector_velocity(i,2),paramStruct.state_vector_velocity(i,3));
end
fprintf(fid, '\n');

fclose(fid);
