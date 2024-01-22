function writeSLCpar(slcParFileName,paramStruct)
%WRITESLCPAR writes the the GAMMA software-style SLC_par parameter
% structure stored in the structure variable paramStruct to the
% processing parameter file with the name given in procParFileName.
%
% Usage:
%   writeSLCpar(slcParFileName,paramStruct)
%
% SEE ALSO:
%  initSLCpar.m, and further related: initSARpar.m, writeSARpar.m, initPROCpar.m,writePROCpar.m
%
% Created:  Othmar Frey <frey@gamma-rs.ch> 10. July 2012
% Modified:  Othmar Frey <frey@gamma-rs.ch>  27. Nov. 2012
% Modified : Othmar Frey <frey@gamma-rs.ch> 26. Sep 2019
%                   Added GPRI parameter set to writeSLCpar
%
%   Copyright: 2019 Gamma Remote Sensing AG
%              Othmar Frey <frey@gamma-rs.ch>
%

fid=fopen(slcParFileName,'wt');
message = ferror(fid);
if(~isempty(message))
	error(message);
	return;
end

fprintf(fid,'Interferometric SAR Processor (ISP) - Image Parameter File\n\n');

fprintf(fid,'title:     %s\n',paramStruct.title);
fprintf(fid,'sensor:    %s\n',paramStruct.sensor);
fprintf(fid,'date:                %d %d %d\n', paramStruct.date(1),paramStruct.date(2),paramStruct.date(3));
fprintf(fid,'start_time:          %15.6f   s\n',paramStruct.start_time);
fprintf(fid,'center_time:         %15.6f   s\n',paramStruct.center_time);
fprintf(fid,'end_time:            %15.6f   s\n',paramStruct.end_time);
fprintf(fid,'azimuth_line_time:   %15.7e  s\n',paramStruct.azimuth_line_time);
fprintf(fid,'line_header_size:       %12d\n',paramStruct.line_header_size);
fprintf(fid,'range_samples:          %12d\n',paramStruct.range_samples);
fprintf(fid,'azimuth_lines:          %12d\n',paramStruct.azimuth_lines);
fprintf(fid,'range_looks:            %12d\n',paramStruct.range_looks);
fprintf(fid,'azimuth_looks:          %12d\n',paramStruct.azimuth_looks);
fprintf(fid,'image_format:               %s\n',paramStruct.image_format);
fprintf(fid,'image_geometry:             %s\n',paramStruct.image_geometry);
fprintf(fid,'range_scale_factor:   %15.7e\n',paramStruct.range_scale_factor);
fprintf(fid,'azimuth_scale_factor: %15.7e\n',paramStruct.azimuth_scale_factor);
fprintf(fid,'center_latitude:      %14.7f   degrees\n',paramStruct.center_latitude);
fprintf(fid,'center_longitude:     %14.7f   degrees\n',paramStruct.center_longitude);
fprintf(fid,'heading:              %14.7f   degrees\n',paramStruct.heading);
fprintf(fid,'range_pixel_spacing:    %12.6f   m\n',paramStruct.range_pixel_spacing);
fprintf(fid,'azimuth_pixel_spacing:  %12.6f   m\n',paramStruct.azimuth_pixel_spacing);
fprintf(fid,'near_range_slc:       %15.4f   m\n',paramStruct.near_range_slc);
fprintf(fid,'center_range_slc:     %15.4f   m\n',paramStruct.center_range_slc);
fprintf(fid,'far_range_slc:        %15.4f   m\n',paramStruct.far_range_slc);
fprintf(fid,'first_slant_range_polynomial:   %12.5f %12.5f %12.5e %12.5e %12.5e %12.5e  s m 1 m^-1 m^-2 m^-3\n',paramStruct.first_slant_range_polynomial(1),paramStruct.first_slant_range_polynomial(2),paramStruct.first_slant_range_polynomial(3),paramStruct.first_slant_range_polynomial(4),paramStruct.first_slant_range_polynomial(5),paramStruct.first_slant_range_polynomial(6));
fprintf(fid,'center_slant_range_polynomial:  %12.5f %12.5f %12.5e %12.5e %12.5e %12.5e  s m 1 m^-1 m^-2 m^-3\n',paramStruct.center_slant_range_polynomial(1),paramStruct.center_slant_range_polynomial(2),paramStruct.center_slant_range_polynomial(3),paramStruct.center_slant_range_polynomial(4),paramStruct.center_slant_range_polynomial(5),paramStruct.center_slant_range_polynomial(6));
fprintf(fid,'last_slant_range_polynomial:    %12.5f %12.5f %12.5e %12.5e %12.5e %12.5e  s m 1 m^-1 m^-2 m^-3\n',paramStruct.last_slant_range_polynomial(1),paramStruct.last_slant_range_polynomial(2),paramStruct.last_slant_range_polynomial(3),paramStruct.last_slant_range_polynomial(4),paramStruct.last_slant_range_polynomial(5),paramStruct.last_slant_range_polynomial(6));
fprintf(fid,'incidence_angle:        %12.4f   degrees\n',paramStruct.incidence_angle);
fprintf(fid,'azimuth_deskew:          %s\n',paramStruct.azimuth_deskew);
fprintf(fid,'azimuth_angle:          %12.4f   degrees\n',paramStruct.azimuth_angle);
fprintf(fid,'radar_frequency:      %15.7e   Hz\n',paramStruct.radar_frequency);
fprintf(fid,'adc_sampling_rate:    %15.7e   Hz\n',paramStruct.adc_sampling_rate);
fprintf(fid,'chirp_bandwidth:      %15.7e   Hz\n',paramStruct.chirp_bandwidth);
fprintf(fid,'prf:                    %12.6f   Hz\n',paramStruct.prf);
fprintf(fid,'azimuth_proc_bandwidth: %12.5f   Hz\n',paramStruct.azimuth_proc_bandwidth);
fprintf(fid,'doppler_polynomial:     %12.5f %12.5e %12.5e %12.5e  Hz     Hz/m     Hz/m^2     Hz/m^3\n',paramStruct.doppler_polynomial(1),paramStruct.doppler_polynomial(2),paramStruct.doppler_polynomial(3),paramStruct.doppler_polynomial(4));
fprintf(fid,'doppler_poly_dot:       %12.5e %12.5e %12.5e %12.5e  Hz/s   Hz/s/m   Hz/s/m^2   Hz/s/m^3\n',paramStruct.doppler_poly_dot(1),paramStruct.doppler_poly_dot(2),paramStruct.doppler_poly_dot(3),paramStruct.doppler_poly_dot(4));
fprintf(fid,'doppler_poly_ddot:      %12.5e %12.5e %12.5e %12.5e  Hz/s^2 Hz/s^2/m Hz/s^2/m^2 Hz/s^2/m^3\n',paramStruct.doppler_poly_ddot(1),paramStruct.doppler_poly_ddot(2),paramStruct.doppler_poly_ddot(3),paramStruct.doppler_poly_ddot(4));
fprintf(fid,'receiver_gain:                %15.4f   dB\n',paramStruct.receiver_gain);
fprintf(fid,'calibration_gain:             %15.4f   dB\n',paramStruct.calibration_gain);
fprintf(fid,'sar_to_earth_center:          %15.4f   m\n',paramStruct.sar_to_earth_center);
fprintf(fid,'earth_radius_below_sensor:    %15.4f   m\n',paramStruct.earth_radius_below_sensor);
fprintf(fid,'earth_semi_major_axis:        %15.4f   m\n',paramStruct.earth_semi_major_axis);
fprintf(fid,'earth_semi_minor_axis:        %15.4f   m\n',paramStruct.earth_semi_minor_axis);
fprintf(fid,'number_of_state_vectors:      %15d\n',paramStruct.number_of_state_vectors);
if (paramStruct.number_of_state_vectors>0)
    fprintf(fid,'time_of_first_state_vector:   %15.6f   s\n',paramStruct.time_of_first_state_vector);
    fprintf(fid,'state_vector_interval:        %15.6f   s\n',paramStruct.state_vector_interval);
    [firstDim, secDim] = size(paramStruct.state_vector_position)
    for i=1:firstDim
        fprintf(fid,'%s: %14.4f  %14.4f  %14.4f   m   m   m\n',['state_vector_position_' num2str(i)],paramStruct.state_vector_position(i,1),paramStruct.state_vector_position(i,2),paramStruct.state_vector_position(i,3));
        fprintf(fid,'%s: %14.5f  %14.5f  %14.5f   m/s m/s m/s\n',['state_vector_velocity_' num2str(i)],paramStruct.state_vector_velocity(i,1),paramStruct.state_vector_velocity(i,2),paramStruct.state_vector_velocity(i,3));
    end
end
if(isfield(paramStruct,'GPRI_TX_mode'))
    fprintf(fid,'GPRI_TX_mode:          %s\n',paramStruct.GPRI_TX_mode);
    fprintf(fid,'GPRI_TX_antenna:       %s\n',paramStruct.GPRI_TX_antenna);
    fprintf(fid,'GPRI_az_start_angle:   %15.6f  degrees\n',paramStruct.GPRI_az_start_angle);
    fprintf(fid,'GPRI_az_angle_step:    %15.6f  degrees\n',paramStruct.GPRI_az_angle_step);
    fprintf(fid,'GPRI_ant_elev_angle:   %15.6f  degrees\n',paramStruct.GPRI_ant_elev_angle);
    fprintf(fid,'GPRI_ref_north:        %15.12f\n',paramStruct.GPRI_ref_north);
    fprintf(fid,'GPRI_ref_east:         %15.12f\n',paramStruct.GPRI_ref_east);
    fprintf(fid,'GPRI_ref_alt:          %15.6f m\n',paramStruct.GPRI_ref_alt);
    fprintf(fid,'GPRI_geoid:            %15.6f m\n',paramStruct.GPRI_geoid);
    fprintf(fid,'GPRI_scan_heading:     %15.6f degrees\n',paramStruct.GPRI_scan_heading);
    fprintf(fid,'GPRI_tx_coord:     %15.6f    %15.6f   %15.6f  m m m\n',paramStruct.GPRI_tx_coord);
    fprintf(fid,'GPRI_rx1_coord:    %15.6f    %15.6f   %15.6f  m m m\n',paramStruct.GPRI_rx1_coord);
    fprintf(fid,'GPRI_rx2_coord:    %15.6f    %15.6f   %15.6f  m m m\n',paramStruct.GPRI_rx2_coord);
    fprintf(fid,'GPRI_tower_roll:   %15.6f  degrees\n',paramStruct.GPRI_tower_roll);
    fprintf(fid,'GPRI_tower_pitch:  %15.6f  degrees\n',paramStruct.GPRI_tower_pitch);
    fprintf(fid,'GPRI_phase_offset: %15.6f  radians\n',paramStruct.GPRI_phase_offset);
end

fclose(fid);