function writeSARpar(sarParFileName,paramStruct)
% WRITESARPAR writes the the GAMMA software-style SAR parameter
% structure stored in the structure variable paramStruct to the
% SAR sensor parameter file (SAR_par) with the name given in sarParFileName.
%
% Usage:
%   writeSARpar(sarParFileName,paramStruct)
%
% SEE ALSO:
%   initPROCpar.m, and further related: initSARpar.m, writeSARpar.m, initSLCpar.m, writeSLCpar.m
%
% Created:  Othmar Frey <frey@gamma-rs.ch>, 10. Aug 2012
% Modified: Othmar Frey <frey@gamma-rs.ch>, 27. Nov. 2012
%


fid=fopen(sarParFileName,'wt');
message = ferror(fid);
if(~isempty(message))
	error(message);
	return;
end

fprintf(fid, 'GAMMA Modular SAR Processor (MSP) - SAR Sensor Parameter file v2.1\n');
fprintf(fid, 'title:                         %s\n', paramStruct.title);
fprintf(fid, 'sensor_name:                   %s\n', paramStruct.sensor_name);
fprintf(fid, 'chirp_direction:               %s\n', paramStruct.chirp_direction);
fprintf(fid, 'receiver_adc_mode:             %s\n', paramStruct.receiver_adc_mode);
fprintf(fid, 'sample_type:                   %s\n', paramStruct.sample_type);
fprintf(fid, 'receiver_spectrum_type:        %s\n', paramStruct.receiver_spectrum_type);
%fprintf(fid, 'SAR_center_frequency:          %13.6e  Hz\n', paramStruct.SAR_center_frequency);
%fprintf(fid, 'chirp_bandwidth:               %13.6e  Hz\n', paramStruct.chirp_bandwidth);
%fprintf(fid, 'chirp_duration:                %13.6e   s\n', paramStruct.chirp_duration);
%fprintf(fid, 'ADC_sampling_frequency:        %13.6e  Hz\n', paramStruct.ADC_sampling_frequency);
fprintf(fid, 'SAR_center_frequency:          %20.15e  Hz\n', paramStruct.SAR_center_frequency);
fprintf(fid, 'chirp_bandwidth:               %20.15e  Hz\n', paramStruct.chirp_bandwidth);
fprintf(fid, 'chirp_duration:                %20.15e   s\n', paramStruct.chirp_duration);
fprintf(fid, 'ADC_sampling_frequency:        %20.15e  Hz\n', paramStruct.ADC_sampling_frequency);
fprintf(fid, 'file_header_size:              %6d         bytes\n', paramStruct.file_header_size);
fprintf(fid, 'record_length:                 %6d         bytes\n', paramStruct.record_length);
fprintf(fid, 'record_header_size:            %6d         bytes\n', paramStruct.record_header_size);
fprintf(fid, 'samples_per_record:            %6d\n', paramStruct.samples_per_record);
fprintf(fid, 'antenna_azimuth_3dB_beamwidth: %10.4f     degrees\n', paramStruct.antenna_azimuth_3dB_beamwidth);
fprintf(fid, 'antenna_range_3dB_beamwidth:   %10.4f     degrees\n', paramStruct.antenna_range_3dB_beamwidth);
fprintf(fid, 'nominal_antenna_azimuth_angle: %10.4f     degrees\n', paramStruct.nominal_antenna_azimuth_angle);
fprintf(fid, 'nominal_antenna_look_angle:    %10.4f     degrees\n', paramStruct.nominal_antenna_look_angle);
fprintf(fid, 'nominal_platform_pitch_angle:  %10.4f     degrees\n', paramStruct.nominal_platform_pitch_angle);
fprintf(fid, 'antenna_pattern_filename:      %s\n', paramStruct.antenna_pattern_filename);
fprintf(fid, '\n');

fclose(fid);