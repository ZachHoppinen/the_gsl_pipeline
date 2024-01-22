function paramStruct = initSARpar
% INITSARPAR intitializes the GAMMA software-style SAR_par parameter structure 
% and returns the structure in the output paramter: paramStruct
%
% Usage:
%   paramStruct = initSARpar
%
% Created: Othmar Frey <frey@gamma-rs.ch>, 03. July 2012
%

paramStruct.title                          = ' ';             % ascii text description of parameter file
paramStruct.sensor_name                    = ' ';             % sensor name (RADARSAT, ERS1, ERS2, JERS-1, PALSAR,...)
paramStruct.chirp_direction                = 'UP_CHIRP';      % chirp direction flag values: (UP_CHIRP, DOWN_CHIRP) 
paramStruct.receiver_adc_mode              = 'NORMAL';        % Receiver ADC mode: (REAL, FMCW, IQ), REAL and FMCW denote offset-video sampling, while
                                                              % IQ is in-phase/quadrature ADC sampling, 2 samples/point
paramStruct.sample_type                    = 'FLOAT';         % sample type (FLOAT, BYTE), floats are 4 bytes/value and bytes
                                                              % are of type unsigned char, 1 byte/value
paramStruct.receiver_spectrum_type         = 'NORMAL';        % SAR receiver spectrum: (NORMAL, INVERT), inverted spectrum caused by the
                                                              % recevier L.O. frequency above the chirp signal spectrum
paramStruct.SAR_center_frequency           = 0.0;             % SAR center frequency (Hz)
paramStruct.chirp_bandwidth                = 0.0;             % chirp bandwidth (Hz)
paramStruct.chirp_duration                 = 0.0;             % chirp duration (sec)
paramStruct.ADC_sampling_frequency         = 0.0;             % range ADC sampling frequency (Hz)
paramStruct.file_header_size               = 0.0;             % SAR raw data file header size (bytes)
paramStruct.record_length                  = 0.0;             % length of each record in bytes
paramStruct.record_header_size             = 0.0;             % record header size in bytes
paramStruct.samples_per_record             = 0.0;             % number of samples/record (IQ pair counts as one sample)
paramStruct.antenna_azimuth_3dB_beamwidth  = 0.0;             % azimuth antenna 3 dB (half-power) beamwidth (decimal deg.)
paramStruct.antenna_range_3dB_beamwidth    = 0.0;             % range antenna 3 dB (half-power) beamwidth (decimal deg.)
                                                              %
                                                              %
                                                              % COORDINATE SYSTEM FOR AZIMUTH, LOOK AND PITCH ANGLES:
                                                              %
                                                              %   N = DOWN
                                                              %   C = N x V /|N x V|
                                                              %   T = C x N
                                                              %
                                                              %   ROLL is measured about T axis CW + (right wing down)
                                                              %   PITCH is measured CW about the temp C axis (nose up +)
                                                              %   YAW is measured CW about +N axis (looking down).
                                                              %
paramStruct.nominal_antenna_azimuth_angle  = 0.0;             % nominal azimuth antenna angle (including ave. squint for simulation)
                                                              % (decimal degrees CW about N, right looking SAR 90.0, left looking: -90.0)*/
paramStruct.nominal_antenna_look_angle     = 0.0;             % nominal antenna look angle (decimal degrees CW (includes ave. roll for simulation)
                                                              % about the final T axis, rotation dependent on right or left looking
paramStruct.nominal_platform_pitch_angle   = 0.0;             % nominal platform pitch angle, CW rot. about the temp. C axis, nose up +
paramStruct.antenna_pattern_filename       = ' ';             % antenna pattern filename (relative to peak, one-way gain) (not in dB!)
