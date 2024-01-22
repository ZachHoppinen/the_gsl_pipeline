function prepareGAMMASARHoneyWellHGuideMoco(procParFile,honeywellNaviTXTFile, leverArmToAPC, outputMocoFile, outputPosXYZ_ECEF_File, outputVelXYZ_ECEF_File,verbose,geoidUndulation)
% prepareGAMMASARHoneyWellHGuideMoco creates a *.moco (ASCII file)
% and ECEF XYZ positions and ECEF V_x,V_y, V_z velocity files navigation data
% at each azimuth echo and optionally plots the data
%
%   Usage:
%       prepareGAMMASARHoneyWellHGuideMoco(procParFile,honeywellNaviTXTFile, leverArmToAPC, outputMocoFile, outputPosXYZ_ECEF_File, outputVelXYZ_ECEF_File,verbose,geoidUndulation)
%
%           procParFile:            GAMMA-style processing parameter file
%           honeywellNaviTXTFile:   txt file with 
%                                   GPS Time, Lat, lon, Height, Roll, Pitch, Yaw
%           leverArmToAPC:          3x1 array [xbody; ybody;zbody] : lever arm from INS origin to antenna phase center of radar antenna
%                                   given in body frame coordinates (as indicated on INS: x = forward, y = left, z = up)
%           outputMocoFile:         *.moco file (ASCII) with data entries at each azimuth receive time
%           outputPosXYZ_ECEF_File: ASCII file containing the antenna phase center position (X_ECEF,Y_ECEF,Z_ECEF) of the radar
%           outputVelXYZ_ECEF_File: ASCII file containing the antenna phase center velocity (VX_ECEF,VY_ECEF,VZ_ECEF) of the radar
%           verbose:                verbose = 1 to plot all the data
%           geoidUndulation:        set the geoid undulation in meter for kml file heights compatible with 
%                                   Google Earth.
%                                   If not given at input, the geoidUndulation is set to 0.0 by default.
%                                   Google Earth uses something close to heights above the
%                                   geoid ("sea level heights"), which may be different from the
%                                   heights above the global WGS84 ellipsoid by +/-several dozens
%                                   of meters depending on the lat/lon on earth.   
%
%   The input navigation (text) file must have the follwing columns:
%   1. GPS Time [s]
%   2. Lat [rad]
%   3. Lon  [rad]
%   4. Height [m]
%   5. Roll [rad]
%   6. Pitch [rad]
%   7. Yaw [rad]
%
% Created: 2016-08-22  frey@gamma-rs.ch
% Modified: 2016-09-20 frey@gamma-rs.ch
%        Updated function (inherited from prepare HPSTT input)
%        1. Now accounts for leap seconds in UTC time when converting to
%           GPS time. Currently valied through 2017.
%           After that, update the lookup table leapsecs.dat in 
%           yourMatlabPathToSARProcessing/GPSTimeUtilitiesODTBX/leapsecs.dat
%           accordingly.
%           See also function utc2gps.m in
%           yourMatlabPathToSARProcessing/GPSTimeUtilitiesODTBX/utc2gps.m
%        2. Also changed:
%           t2_gps_sec_of_week = t0_gps_sec_of_week + procPar.total_raw_echoes*sarPar.chirp_duration;
%           to:
%           t2_gps_sec_of_week = t0_gps_sec_of_week + procPar.total_raw_echoes/procPar.pulse_repetition_frequency;
%        3. As a consequence of change no. 2 the *.sar_par parameter file is no longer
%           required. Changed usage/parameters accordingly.
% Modified: 2016-09-23 frey@gamma-rs.ch  (inherited from prepare HPSTT input)
%       Added check for azimuth times of radar acquisition versus times
%       available in the HPSTT scan log file before the interpolation is
%       taking place.
% Modified: 2018-08-23 frey@gamma-rs.ch (inherited from prepare HPSTT input)
%       kml file generation of tracks are now only with every 1000th
%       position to avoid display bug in Google Earth
%
% Modified: 2019-03-08 frey@gamma-rs.ch: azimuth receive times in GPS time 
%    dlmwrite([outputPosXYZ_ECEF_File '.azi_times'], AziReceiveTimes, 'delimiter', ',', 'precision', '%.10f');    
%
% Modified: 2020-01-09 frey@gamma-rs.ch
%       Added additional checks and more robust handling of year month day
%       or day month year datum format in coming from the date field. 
%      and
%       replaced kmlwriteline by our custom function InputOutput/mykmlwriteline
%       for compatibility with octave.
%
% Modified: 2020-02-24 frey@gamma-rs.ch
%       - Cleaned up help text and code.
%       - additional input parameter geoidUndulation
%         if not given at input, the geoidUndulation is set to 0.0 by default.
%
% Modified: 2021-10-19 frey@gamma-rs.ch
%       - replaced
%           nav_data = load(honeywellNaviTXTFile);
%         by 
%           nav_data = dlmread(honeywellNaviTXTFile);
%         to also support space-delimited in addition to the comma-separated input ascii file.
% 

    % Lever arm of radar antenna phase center (APC)
    % with respect to INS coordinate system
    xbody = leverArmToAPC(1);
    ybody = leverArmToAPC(2);
    zbody = leverArmToAPC(3);
    
    if(~exist('verbose'))
        verbose = 0;
    end
    
    if(~exist('geoidUndulation'))
        geoidUndulation = 0.0;
    end

    % Read the Honeywell navigation output data by AeroScout
    nav_data = dlmread(honeywellNaviTXTFile);
    
    % Read synthetic aperture radar parameter files
    procPar = readGammaParFile(procParFile);

    % Convert data and time to GPS second of the week
    if(procPar.date(1)>1900 && procPar.date(2)>=1 && procPar.date(2)<=12 && procPar.date(3)>=1 && procPar.date(3)<=31) 
        year = procPar.date(1);
        month = procPar.date(2);
        day = procPar.date(3);
    elseif(procPar.date(3)>1900 && procPar.date(2)>=1 && procPar.date(2)<=12 && procPar.date(1)>=1 && procPar.date(1)<=31)
        year = procPar.date(3);
        month = procPar.date(2);
        day = procPar.date(1);
    else
        error('Wrong year month day format in file: %s',procParFile);
    end
        
    hrs = procPar.raw_data_start_time(1);
    minutes = procPar.raw_data_start_time(2);
    sec = procPar.raw_data_start_time(3);
    
    %julianday = cal2jd(year,month,day);
    %[gpsweek,gps_sec_of_week_at_0h,rollover]=jd2gps(julianday);
    
    UTC_time = [year month day hrs minutes sec];
    [GPS_week, GPS_sec_of_week, GPS_day] = utc2gps(UTC_time); % Accounts for leap seconds in UTC time: currently 17sec (Sep. 2016)

    %t0_gps_sec_of_week = gps_sec_of_week_at_0h + hrs * 3600. + minutes * 60. + sec + procPar.azimuth_offset;
    t0_gps_sec_of_week = GPS_sec_of_week + procPar.azimuth_offset;
    %t2_gps_sec_of_week = t0_gps_sec_of_week + procPar.total_raw_echoes*sarPar.chirp_duration;
    t2_gps_sec_of_week = t0_gps_sec_of_week + (procPar.total_raw_echoes-1)/procPar.pulse_repetition_frequency;
    t1_gps_sec_of_week = (t2_gps_sec_of_week+t0_gps_sec_of_week)/2.;

    AziReceiveTimes = (t0_gps_sec_of_week:(1.0/procPar.pulse_repetition_frequency):t2_gps_sec_of_week).';
    
    disp(' ');
    disp(['Radar raw data start time (UTC): YYYY-MM-DD : ' num2str(year) '-' num2str(month) '-' num2str(day) ...
    ', HH-MM-SS.SSSS :'  num2str(hrs) '-' num2str(minutes) '-' num2str(sec)]);
    disp(['Radar raw data start and end time (GPS sec of week):       ' num2str(t0_gps_sec_of_week) ' [sec] '  num2str(t2_gps_sec_of_week) ' [sec]']);
    disp(['HPSTT Scan log file start and end times (GPS sec of week): ' num2str(nav_data(1,1)) ' [sec] ' num2str(nav_data(end,1)) ' [sec]']);
    disp(' ');
    
%    =============================================
%  	Interpolate all position and attitude data for each echo receive time by
%  	spline interpolation of the IMU/GPS data:
	disp('interpolating position, velocity and attitude at echo times');	
	
    % Check whether the echo times are not exceeding (before and/or after)
    % the times at which the INS/GNSS data of the Scan log file atre
    if(AziReceiveTimes(1)<nav_data(1,1))
        error('The GPS time of the first radar echo:         %20.10f [sec] \nis earlier than the \nfirst GPS time of the INS/GNSS Scan log file: %20.10f [sec].\nMake sure you are using the corresponding PROC_par and HPSTT_Scan_log files.\nNow you are using:\n%s and \n%s\n',AziReceiveTimes(1),nav_data(1,1),procParFile,honeywellNaviTXTFile); 
        return;
    end
    if(AziReceiveTimes(end)>nav_data(end,1))
        error('The GPS time of the last radar echo:          %20.10f [sec] \nis later than the \nlast GPS time of the INS/GNSS Scan log file:    %20.10f [sec].\nMake sure you are using the corresponding PROC_par and HPSTT_Scan_log files.\nNow you are using:\n%s and \n%s\n',AziReceiveTimes(end),nav_data(end,1),procParFile,honeywellNaviTXTFile); 
        return;
    end
    
	interptype = 'pchip';
	IntpGPSLat		    = interp1(nav_data(:,1),nav_data(:,2),AziReceiveTimes,interptype);
	IntpGPSLon		    = interp1(nav_data(:,1),nav_data(:,3),AziReceiveTimes,interptype);
	IntpGPSHeight 		= interp1(nav_data(:,1),nav_data(:,4),AziReceiveTimes,interptype);
	IntpRoll	     	= interp1(nav_data(:,1),nav_data(:,5),AziReceiveTimes,interptype);
	IntpPitch           = interp1(nav_data(:,1),nav_data(:,6),AziReceiveTimes,interptype);
	IntpYaw             = interp1(nav_data(:,1),nav_data(:,7),AziReceiveTimes,interptype);
	  
    % Convert geodetic coordinates to Cartesian coordinates
    EARTH.A_ELLIPSOID = 6378137.0;
	EARTH.F_ELLIPSOID = 1.0/298.257223560;
	EARTH.MU          = 3.986004418e14;
    EARTH.OMEGA       = [0.0 0.0 0.0];
    % =============================================	
    % transform GPS/INS positions from geodetic to Cartesian coords.
	disp('transform GPS/INS positions from geodetic to Cartesian coords.');
    [X_INS_ECEF,Y_INS_ECEF,Z_INS_ECEF] = geod2cartesianForArrays(IntpGPSLat,IntpGPSLon,IntpGPSHeight,EARTH);

    %[VX_ECEF,VY_ECEF,VZ_ECEF] = topocentric2geocentric(IntpGPSVelNorthing,IntpGPSVelEasting,-IntpGPSVelDown,IntpGPSLat,IntpGPSLon);
    VX_ECEF = zeros(size(X_INS_ECEF));
    VY_ECEF = zeros(size(Y_INS_ECEF));
    VZ_ECEF = zeros(size(Z_INS_ECEF));
        
    VX_ECEF(1:(end-1)) = diff(X_INS_ECEF);
    VY_ECEF(1:(end-1)) = diff(Y_INS_ECEF);
    VZ_ECEF(1:(end-1)) = diff(Z_INS_ECEF);
    
    VX_ECEF(end) = VX_ECEF(end-1);
    VY_ECEF(end) = VY_ECEF(end-1);
    VZ_ECEF(end) = VZ_ECEF(end-1);
    
    dt = AziReceiveTimes(2,1)-AziReceiveTimes(1,1);
    
    VX_ECEF = VX_ECEF./dt;
    VY_ECEF = VY_ECEF./dt;
    VZ_ECEF = VZ_ECEF./dt;
    
    % Lever arm between INS coord origin and radar antenna phase center
    
    % =============================================		
    % transform lever arm coords given in the B-frame to the NED-frame
	disp('transform lever arm coords given in the B-frame to the NED-frame');
	
	[dNorth, dEast, dDown] = BFrame2NEDFrame(xbody,ybody,zbody,IntpRoll,IntpPitch,IntpYaw);
    
      
    % transform lever arm coords from topocentric to the geocentric system
	disp('transform lever arm coords from topocentric to the geocentric system');
	
	[dX,dY,dZ] = topocentric2geocentric(dNorth,dEast,-dDown,IntpGPSLat,IntpGPSLon);
	

    % =============================================	
    % add difference vector of lever arm to the GPS/INS positions to get the
    % position of the antenna phase center (APC)
	disp('Add difference vector of lever arm to the GPS/INS positions'); 
	disp('to get the position of the phase center of the antenna'); 
	X_APC_ECEF = X_INS_ECEF+dX;
	Y_APC_ECEF = Y_INS_ECEF+dY;
	Z_APC_ECEF = Z_INS_ECEF+dZ;

    % Write the WGS84 Positioning data for radar antenna phase center
    dlmwrite(outputPosXYZ_ECEF_File, [X_APC_ECEF, Y_APC_ECEF, Z_APC_ECEF], 'delimiter', ',', 'precision', '%.10f');
    % Write the WGS84 Velocity data
    dlmwrite(outputVelXYZ_ECEF_File, [VX_ECEF, VY_ECEF, VZ_ECEF], 'delimiter', ',', 'precision', '%.10f');

    % NEW 2019-03-08: azimuth receive times in GPS time 
    dlmwrite([outputPosXYZ_ECEF_File '.azi_times'], AziReceiveTimes, 'delimiter', ',', 'precision', '%.10f');    
    
    APC_lat = zeros(size(X_APC_ECEF));
    APC_lon = zeros(size(X_APC_ECEF));
    APC_H   = zeros(size(X_APC_ECEF));
     % Write the WGS84 Positioning data for radar antenna phase center in
     % geodetic coordinates
    for i=1:length(X_APC_ECEF)
       [APC_lat(i), APC_lon(i), APC_H(i)] = cartesian2geod([X_APC_ECEF(i), Y_APC_ECEF(i), Z_APC_ECEF(i)],EARTH);
    end
    
    dlmwrite([outputPosXYZ_ECEF_File '.geod'], [APC_lat.*180/pi, APC_lon.*180/pi, APC_H],'delimiter', ',', 'precision', '%.12f');
     
    mykmlwriteline([outputPosXYZ_ECEF_File 'geod.kml'], APC_lat(1:100:end).*180/pi, APC_lon(1:100:end).*180/pi,APC_H(1:100:end),geoidUndulation,3,'FF0000FF');
   
    if(verbose)
        % Plot to check data
        figure, plot(AziReceiveTimes,'LineStyle','none','Marker','.');
        xlabel('Data Sample Number');
        ylabel('GPS sec of week [sec] (col 1)');
        grid on;
            
        figure
        line(IntpGPSLon,IntpGPSLat,'Color','red');
        line(IntpGPSLon,IntpGPSLat,'LineStyle','none','Marker','.');
        title('INS/GNSS Positions (Lon/Lat) after HPSTT (Scan file)');
        xlabel('Lon [deg]');
        ylabel('Lat [deg]');
        grid on;
        axis equal;
    
        figure, plot(AziReceiveTimes,IntpGPSHeight,'LineStyle','none','Marker','.');
        xlabel('GPS sec of week [sec] (col 1)');
        ylabel('Altitude [m] (col 4)');
        grid on;
        
        figure, plot(AziReceiveTimes,[IntpYaw, IntpPitch,  IntpRoll],'LineStyle','none','Marker','.');
        xlabel('GPS sec of week [sec] (col 1)');
        ylabel('Attitude euler angles [deg]  ');
        legend('yaw','pitch','roll');
        grid on;
       
        figure, plot3(dNorth, dEast, dDown,'LineStyle','none','Marker','.');
        xlabel('dNorth [m]');
        ylabel('dEast [m]');
        zlabel('dDown [m]');
        grid on;
        axis equal;
        
        figure,
        line(AziReceiveTimes, dNorth,'LineStyle','none','Marker','.');
        xlabel('GPS sec of week [sec] (col 1)');
        ylabel('dNorth [m]');
        
        figure,
        line(AziReceiveTimes, dEast,'LineStyle','none','Marker','.');
        xlabel('GPS sec of week [sec] (col 1)');
        ylabel('dEast [m]');
        
        figure,
        line(AziReceiveTimes, dDown,'LineStyle','none','Marker','.');
        xlabel('GPS sec of week [sec] (col 1)');
        ylabel('dDown [m]');
        
    end
end
