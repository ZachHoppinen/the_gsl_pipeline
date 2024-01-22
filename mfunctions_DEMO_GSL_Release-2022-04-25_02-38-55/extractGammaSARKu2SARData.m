function extractGammaSARKu2SARData(inRawFilename,inRawParFilename,outRawCpxFilename,outPROCparFilename,outSARparFilename,rg_3dB_bw,az_3dB_bw, polarization,self_zero,rg_bw_frac,rg_kaiser_par,az_bw_frac,az_kaiser_par)
% extractGammaSARKu2SARData extracts the given channnel of the real-valued 
% Gamma GPRI-2 raw data (new version of raw data with 2 separate raw data files
% for each of the 2 simultaneous receive channels and then
% 1. performs a Hilbert transform and range spectral windowing and then
%    writes the complex-valued beat signal data as float-complex data to
%    disk,
% 2. it also writes the corresponding SAR_par and PROC_par files.
%
% Usage:
%    extractGammaSARKu2SARData(inRawFilename,inRawParFilename,outRawCpxFilename,outPROCparFilename,outSARparFilename,rg_3dB_bw,az_3dB_bw, polarization,self_zero,rg_bw_frac,rg_kaiser_par,az_bw_frac,az_kaiser_par)
%
%
%
%    where:
%       inRawFilename       : GAMMA SAR (GPRI-II) raw, real-valued beat-signal (*.raw)
%       inRawParFilename    : corresponding raw parameter file (*.raw_par)
%       outRawCpxFilename   : complex-valued beat signal (*.raw_cplx)
%       outPROCparFilename  : GAMMA-style processing parameter file (PROC_par)
%       outSARparFilename   : GAMMA-style SAR system parameter file (SAR_par)
%       rg_3dB_bw           : 3dB antenna beamwidth in elevation direction
%       az_3dB_bw           : 3dB antenna beamwidth in azimuth direction
%       polarization        : polarization
%       self_zero           : number of samples to set to zero at the
%                             beginning
%       rg_bw_frac          : centered fraction of the total range spectrum
%                             (0.x - 1.0) used  for windowing the range spectrum
%       rg_kaiser_par       : kaiser window parameter for range window
%       az_bw_frac          : centered fraction of the total azimuth spectrum (0.x - 1.0)
%                             (0.x - 1.0) used  for windowing the azimuth spectrum
%       az_kaiser_par       : kaiser window parameter for azimuth window
%
%
% Created:         2021-09-20 by Othmar Frey, adapted from extractGPRI2SARDataDerampFFT 1;
%                   copy-paste-modified from extractGammaSARL2SARData.m
%                   
%
%   Copyright: 2021 Gamma Remote Sensing AG
%              Othmar Frey <frey@gamma-rs.ch>
%

orig_blocksize = 2000;

% set to 1 for debugging and data visualization
verbose = 0;

% Default values for optional parameters
if(~exist('rg_bw_frac'))
    rg_bw_frac = 0.8;
    disp('Setting default value: rg_bw_frac = 0.8')
end
if(~exist('rg_kaiser_par'))
    rg_kaiser_par = 4.0;
    disp('Setting default value: rg_kaiser_par = 4.0')
end
if(~exist('self_zero'))
    self_zero = 300;
    disp('Setting default value: self_zero = 300')
end

cj = sqrt(-1);
c = 299792458.0;

ra = 6378137.0000;    %WGS-84 semi-major axis
rb = 6356752.3142;    %WGS-84 semi-minor axis

EARTH.A_ELLIPSOID = 6378137.0;
EARTH.F_ELLIPSOID = 1.0/298.257223563;
EARTH.MU          = 3.986004418E14;
%EARTH.OMEGA       = [0.0 0.0 7.292115E-5];

% Read GS-L raw parameter file
rawParStruct  = readGammaParFile(inRawParFilename);

rangeDim = rawParStruct.CHP_num_samp;
block_length = rawParStruct.CHP_num_samp+1;

SamplingRate = rawParStruct.ADC_sample_rate/2;   % sampling rate of complex-valued beat signal (= deramped signal) 
NrOfRgPix    =  floor(rawParStruct.CHP_num_samp/2);     % number of samples of complex-valued beat signal (= deramped signal) 

rawParStruct.TX_pol_mode = 'VV'
%if (strcmp(rawParStruct.TX_pol_mode,'HV'))	%check if alternating TX antennas
%    tcycle = 2.0*block_length/rawParStruct.ADC_sample_rate;    % time per transmit cycle equals 2* IPP in alt. transmit mode
%else
tcycle = block_length/rawParStruct.ADC_sample_rate;
%end 
fprintf('inter-pulse period (s): %e\n',(block_length/rawParStruct.ADC_sample_rate));
fprintf('transmit cycle time (s): %e\n', tcycle);
fprintf('ADC capture time (s): %.5f\n',rawParStruct.ADC_capture_time);

valtype = 'int16';
outvaltype = 'float32';
isComplex = 0;
firstPix = 1;
%sizeof_data = 8;        %number of bytes per complex sample (for float complex: 8 bytes)
sizeof_data = 4;        %number of bytes per complex sample (for float complex: 8 bytes)
%verbose = 0;


fileinfo = dir(inRawFilename);
filesize = fileinfo.bytes;

sizeof_data = 2;     % number of bytes per channel sample
num_channels = 1;    % number of data channels in one data set
                     % #number of receiver ADC channels,
                     %  GPRI had 2 channels in one data set
                     %  GS-L SAR is different. Has only 1 channel per data set.
                     % corrected on 2018-07-22

bytes_per_record = num_channels * sizeof_data * block_length;  %number of bytes per echo
nl_tot = floor(filesize/bytes_per_record);
NrOfEchoes = nl_tot;
aziDim = NrOfEchoes;

RANGE_OFFSET = 0;

% list of slant range pixel numbers
pn1 = 0:(floor(rangeDim/2)-1); 
% range pixel spacing (corrected version talking abs() of negative chirp rate for downchirp)
rps = (rawParStruct.ADC_sample_rate/rawParStruct.CHP_num_samp*c/2.)/abs(rawParStruct.RF_chirp_rate); % range pixel spacing
%slant range for each sample
slr = (pn1 *rps) + RANGE_OFFSET;
%scale = ((abs(slr)./slr(round(rawParStruct.CHP_num_samp/8))).^1.5)';
%scale = ones(size(scale));
%keyboard

% Write PROC_par parameter file
sub_rout_writePROCpar(rawParStruct,outPROCparFilename,aziDim,tcycle,ra,rb,slr,rps,NrOfRgPix,polarization);
% Write SAR_par parameter file
sub_rout_writeSARpar(rawParStruct,outSARparFilename,c,rps,bytes_per_record,rg_3dB_bw,az_3dB_bw);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract and process raw GS-L data azimuth-blockwise

% Blockwise calculate overall DC bias 
DC_BIAS = 0.0;
blocksize = orig_blocksize;
numOfBlocks = floor(aziDim/blocksize);
sizeRemainder = rem(aziDim,blocksize);
if(sizeRemainder>0)
    numOfBlocks = numOfBlocks+1;
end
for i=1:numOfBlocks
    firstEcho = (i-1)*blocksize+1;
    if(sizeRemainder>0 && i==numOfBlocks) % if there are remaining echoes then change the AziBlockSize in the last round of the for loop.
        blocksize = sizeRemainder;
    end
    numEchoes = blocksize;
    rawData = readMatrixNoHeaderLittleEndian(inRawFilename,block_length, aziDim, valtype,isComplex, firstPix,block_length,firstEcho,numEchoes);
    DC_BIAS = mean(mean(rawData((self_zero+1):end,:)))*blocksize + DC_BIAS;
end
DC_BIAS = DC_BIAS/aziDim;

% Blockwise perform Hilbert transform and range spectral windowing
%outRawFilename = regexprep(inRawFilename,'.raw','.raw_cplx_noSpecShiftNew')
fidout_raw_cpx_noshift = fopen(outRawCpxFilename,'wb','ieee-be');
blocksize = orig_blocksize;
numOfBlocks = floor(aziDim/blocksize);
sizeRemainder = rem(aziDim,blocksize);
if(sizeRemainder>0)
    numOfBlocks = numOfBlocks+1;
end
for i=1:numOfBlocks
    firstEcho = (i-1)*blocksize+1;
    if(sizeRemainder>0 && i==numOfBlocks) % if there are remaining echoes then change the AziBlockSize in the last round of the for loop.
        blocksize = sizeRemainder;
    end
    numEchoes = blocksize;
    rawData = readMatrixNoHeaderLittleEndian(inRawFilename,block_length, aziDim, valtype,isComplex, firstPix,block_length,firstEcho,numEchoes);

    % Range-compressed = fft of Real-valued raw data plots
    rawData_first300toZero = rawData;
    rawData_first300toZero(1:self_zero,:) = 0.0;
    
    self_win = kaiserWindow(block_length,rg_bw_frac,0,rg_kaiser_par);

    % Remove DC bias
    rawData_first300toZero = rawData_first300toZero-DC_BIAS;

    % Prepare a windowed but NON-SPECTRUM-SHIFTED raw_cplx data set
    range_compressed = fft(rawData_first300toZero.*(self_win*ones(1,numEchoes)),[],1);
    raw_cplx = ifft(range_compressed(1:NrOfRgPix,:));
    
    if(exist('az_bw_frac') && exist('az_kaiser_par'))
        % New azimuth-windowing in the azimuth frequency domain to remove artifacts.     
        self_win_az_spec = kaiserWindow(numEchoes,az_bw_frac,0,az_kaiser_par).';
        
        % Azimuth FFT and windowing in azimuth of range-windowed beat signal
        az_filtered_data = fftshift(fft(raw_cplx,[],2),2) .* (ones(NrOfRgPix,1)*self_win_az_spec);
        %figure, imagesc(20*log10(abs(az_filtered_data)))
    
        % Inverse azimuth FFT to obtain azimuth filtered (and range window filtered) raw_cplx data  
        raw_cplx_az_filt = ifft(ifftshift(az_filtered_data,2),[],2);
        %figure, imagesc(20*log10(abs(raw_cplx_az_filt))), colormap jet

        successful = writeComplexDataBlockwise(fidout_raw_cpx_noshift, raw_cplx_az_filt, NrOfRgPix, numEchoes, outvaltype);
        % End of azimuth windowing of az. spectrum
    else
        successful = writeComplexDataBlockwise(fidout_raw_cpx_noshift, raw_cplx, NrOfRgPix, numEchoes, outvaltype);
    end
end

fclose(fidout_raw_cpx_noshift);



% Local subroutines to write the SAR_par and the PROC_par files
% Write SAR_par and PROC_par files (after gpri2_rc_of.py or gpri2_rc.py respectively) 

function sub_rout_writePROCpar(rawParStruct,outPROCparFilename,aziDim,tcycle,ra,rb,slr,rps,NrOfRgPix,polarization)

strucPROCPar = initPROCpar;

ts = rawParStruct.time_start;
strucPROCPar.title = ts;
strucPROCPar.sensor = 'Gamma SAR GS-L';

fprintf('\nraw data start time: %s\n',rawParStruct.time_start);

cell_time_start = strsplit(rawParStruct.time_start, ' ');   % split the line on whitespace
ymd = strsplit(cell_time_start{1}, '-');                    % then parse using '-' for date
hms_tz = strsplit(cell_time_start{2}, '+');                 % split into HMS and time zone information
if(length(hms_tz)==1) 
    hms_tz{2} = '00';
end
fprintf('hms_tz: %s %s\n',hms_tz{1},hms_tz{2});
hms = strsplit(hms_tz{1},':');                              % split HMS string using :
if(length(hms)==1) 
    hms{2} = '00'; 
    hms{3} = '00.00';
end
sod = str2double(hms{1})*3600.0 + str2double(hms{2})*60.0 + str2double(hms{3});     % raw data starting time, seconds of day
        
% Calculate the new raw_data_start_time of the selected subset of echoes 
% in seconds of the day 

strucPROCPar.total_raw_echoes = aziDim;
strucPROCPar.offset_to_first_echo_to_process = 0;
azs = strucPROCPar.offset_to_first_echo_to_process;
strucPROCPar.echoes_to_process = aziDim;

sod_sel = sod + azs*tcycle;
        
% Convert the new raw_data_start_time to the hh:mm:ss.ssss format of the PROC_par file 
ss = rem(sod_sel,60);
rem_mm = floor(sod_sel/60);
mm = rem(rem_mm,60);
hh = floor(rem_mm/60);

fprintf('\nGPS seconds of day (sod) of raw data set: %12.10f',sod);
fprintf('\nGPS seconds of day (sod) of selected start echo: %12.10f',sod_sel)
fprintf('\nExtracting and range-compressing of raw data:');
fprintf('\nAzimuth index of start echo (starting from 0): %d',azs);
fprintf('\nAzimuth index of last echo  (starting from 0): %d', aziDim-1);

strucPROCPar.date(1) = str2double(ymd{1});
if(~isnan(str2double(ymd{2})))
    strucPROCPar.date(2) = str2double(ymd{2});
else
    strucPROCPar.date(2) = monthname2monthnumber(upper(ymd{2}))
end
strucPROCPar.date(3) = str2double(ymd{3});

strucPROCPar.raw_data_start_time(1) = str2double(hms{1});
strucPROCPar.raw_data_start_time(2) = str2double(hms{2});
strucPROCPar.raw_data_start_time(3) = str2double(hms{3});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This needs to be adjusted accordingly for other polarization modes
strucPROCPar.channel_mode = polarization;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

strucPROCPar.earth_semi_major_axis = ra;
strucPROCPar.earth_semi_minor_axis = rb;
%       strucPROCPar('scene_center_latitude:   0.0000   decimal degrees')
%       strucPROCPar('scene_center_longitude:  0.0000   decimal degrees')
%       strucPROCPar('track_angle:             0.0000   degrees')
%       strucPROCPar('platform_altitude:       0.0000   m')
%       strucPROCPar('terrain_height:          0.0000   m')
%       strucPROCPar('sensor_position_vector:     0.000000       0.000000       0.000000   m     m      m')
%       strucPROCPar('sensor_velocity_vector:  %f   0.000000       0.000000   m/s   m/s   m/s '%self.ave_vel)
%       strucPROCPar('sensor_acceleration_vector: 0.000000       0.000000       0.000000   m/s^2 m/s^2 m/s^2')
       
prf = abs(1.0/tcycle);
strucPROCPar.pulse_repetition_frequency = prf;

%        strucPROCPar('I_bias:     0.0')
%        strucPROCPar('Q_bias:     0.0')
%        strucPROCPar('I_sigma:    0.0')
%        strucPROCPar('Q_sigma:    0.0')
%        strucPROCPar('IQ_corr:    0.0')
%        strucPROCPar('SNR_range_spectrum:   0.0')
%        strucPROCPar('DAR_doppler:          0.0   Hz')
%        strucPROCPar('DAR_snr:              0.0')
%        strucPROCPar('doppler_polynomial:   0.00000e+01  0.00000e+00  0.00000e+00  0.00000e+00   Hz     Hz/m     Hz/m^2     Hz/m^3')
%        strucPROCPar('doppler_poly_dot:     0.00000e+00  0.00000e+00  0.00000e+00  0.00000e+00   Hz/s   Hz/s/m   Hz/s/m^2   Hz/s/m^3')
%        strucPROCPar('doppler_poly_ddot:    0.00000e+00  0.00000e+00  0.00000e+00  0.00000e+00   Hz/s^2 Hz/s^2/m Hz/s^2/m^2 Hz/s^2/m^3')
%        strucPROCPar('echo_time_delay:      0.000000   s')
%        strucPROCPar('receiver_gain:        %f  dB'%(60 - self.grp.atten_dB))	#gain is 60 - attenuatio


%strucPROCPar.receiver_gain = 60.0 - ???;  % Still to be set: waiting for Charles' input
strucPROCPar.calibration_gain = 0.0; %  dB
strucPROCPar.near_range_raw    = slr(1);
strucPROCPar.center_range_raw  = (slr(end) + slr(1))/2.0;
strucPROCPar.far_range_raw     = slr(end);
strucPROCPar.near_range_slc    = slr(1);
strucPROCPar.center_range_slc  = (slr(end) + slr(1))/2.0;
strucPROCPar.far_range_slc     = slr(end);
strucPROCPar.range_pixel_spacing = rps;
%        strucPROCPar('range_resolution:     %f   m'%(1.2*self.rps))
%        strucPROCPar('sec_range_migration:  OFF')
%        strucPROCPar('azimuth_deskew:       OFF')
%        strucPROCPar('autofocus_snr:        0.0000')
%        strucPROCPar.azimuth_bandwidth_fraction: 0.8')
strucPROCPar.azimuth_presum_factor = 1;
strucPROCPar.total_raw_echoes = aziDim;
%strucPROCPar.offset_to_first_echo_to_process = 0;
strucPROCPar.echoes_to_process = aziDim;
strucPROCPar.range_offset = 0;
strucPROCPar.raw_range_samples = NrOfRgPix;

strucPROCPar.range_looks    = 1;
strucPROCPar.azimuth_looks  = 1;
%        strucPROCPar('azimuth_offset:          0.0   s')
%        strucPROCPar('azimuth_pixel_spacing:   0.0   m')
%       strucPROCPar('azimuth_resolution:      0.0   m')
strucPROCPar.range_pixels   = NrOfRgPix;            %d  image output samples'%self.ns_out)
strucPROCPar.azimuth_pixels = aziDim;               %d  image output lines'%self.naz)       
strucPROCPar.image_format   = 'FCOMPLEX';
%        strucPROCPar('sensor_latitude:         0.0   decimal degrees')
%        strucPROCPar('sensor_longitude:        0.0   decimal degrees')
%        strucPROCPar('sensor_track_angle:      0.0  decimal degrees')
strucPROCPar.number_of_state_vectors = 0;
strucPROCPar.time_of_first_state_vector = 0.0;
strucPROCPar.state_vector_interval = 0.0;
strucPROCPar.GPRI_TX_mode = rawParStruct.TX_pol_mode;

%TX_antenna = 'VV' % this has to be given as an option at the command line. check again with full pol variant:

% Add some GPRI specific information to the PROC_par file
%strucPROCPar.GPRI_TX_antenna = TX_antenna;
%strucPROCPar.GPRI_ant_elev_angle = rawParStruct.antenna_elevation;  % antenna elevation angle in degrees
% #    for str in p[0:5]:		#write the rest of the MSP processing parameter file lines
% #       self.proc_par1.write('%s\n'%str)
% #       self.proc_par2.write('%s\n'%str)
% #
% #       self.proc_par1.write('%s %s\n'%(p[5],'CH1 lower'))
% #       self.proc_par2.write('%s %s\n'%(p[5],'CH2 upper'))
% #
% #    for str in p[6:]:		#write the rest of the MSP processing parameter file lines
% #       self.proc_par1.write('%s\n'%str)
% #       self.proc_par2.write('%s\n'%str)
% #
% #       self.proc_par1.close()
% #       self.proc_par2.close()
% #

if (exist(outPROCparFilename, 'file')==2)
  error('File aleady exists! Abort writing PROC_par files.');
else
    writePROCpar(outPROCparFilename,strucPROCPar);
end

function sub_rout_writeSARpar(rawParStruct,outSARparFilename,c,rps,bytes_per_record,rg_3dB_bw,az_3dB_bw);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write SAR_par files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
strucSARPar = initSARpar;

fadc = c/(2.0*rps);                                 % effective ADC sampling rate for a pulsed radar rather than an FM-CW system
cbw = rawParStruct.RF_freq_max - rawParStruct.RF_freq_min;
ts = rawParStruct.time_start;

chirp_duration = rawParStruct.CHP_num_samp/rawParStruct.ADC_sample_rate;

strucSARPar.title = ts;
strucSARPar.sensor_name = 'Gamma SAR Ku-band (GPRI-II)';
strucSARPar.chirp_direction = 'UP_CHIRP';
strucSARPar.receiver_adc_mode = 'FMCW';
strucSARPar.sample_type = 'SHORT';
strucSARPar.receiver_spectrum_type = 'NORMAL';
%strucSARPar.SAR_center_frequency = rawParStruct.RF_center_freq; % Hz
strucSARPar.SAR_center_frequency = (rawParStruct.RF_freq_max + rawParStruct.RF_freq_min)/2.0;
strucSARPar.chirp_bandwidth = cbw;                  % Hz
strucSARPar.chirp_duration = chirp_duration;        %  sec 
strucSARPar.ADC_sampling_frequency = fadc;          % Hz'%fadc)
strucSARPar.file_header_size = 0;                   % bytes
strucSARPar.record_length = bytes_per_record;       % bytes
strucSARPar.record_header_size = 0;                 % bytes
strucSARPar.samples_per_record = rawParStruct.CHP_num_samp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The beamwidth will have to be adjusted accordingly depending on the
% antennas used.
strucSARPar.antenna_azimuth_3dB_beamwidth = az_3dB_bw;    % degrees
strucSARPar.antenna_range_3dB_beamwidth   = rg_3dB_bw;    % degrees
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

strucSARPar.nominal_antenna_azimuth_angle = 90.0; % degrees
strucSARPar.nominal_antenna_look_angle    = 90.0; % degrees
strucSARPar.nominal_platform_pitch_angle =  0.0;  % degrees
strucSARPar.antenna_pattern_filename = 'constant_antenna.gain';
      
if(exist(outSARparFilename, 'file')==2)
  error('File aleady exists! Abort writing SAR_par files.');
else
    writeSARpar(outSARparFilename,strucSARPar);
end
