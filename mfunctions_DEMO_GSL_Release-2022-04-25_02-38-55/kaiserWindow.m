function kaiserWin = kaiserWindow(arraySize,BandWidthFraction, DisplacementFromCenterOfArray,betaValue)
%KAISERWINDOW creates a Kaiser window within
%	the 'BandWidthFraction' of the 'arraySize'.
%	The window is centered around:
%	center index + DisplacementFromCenterOfArray .
%	A positive/negative 'DisplacementFromCenterOfArray'
%	is equivalent with a circular shift to the right/left from
% 	the center of the array.
%
%	The window is centered around:
%	center index + DisplacementFromCenterOfArray .
%	A positive/negative 'DisplacementFromCenterOfArray'
%	is equivalent with a circular shift to the right/left from
% 	the center of the array.
%
%   Usage:
%      kaiserWin = kaiserWindow(arraySize,BandWidthFraction, DisplacementFromCenterOfArray,betaValue)
%
%	BETA = 2.120; // -30 dB stopband ripple in Kaiser window filter
%	see also ref. [1] pp.544-548
%
%	Kaiser Window Parameter:
%
%	   betaValue    passband ripple   stopband ripple            PSLR [dB]
%                                         (alpha (positive sign
%                                         for eq. below))
%	***********************************************************************
%	   1.000      .86                -20 dB                      -13.3 dB
%	   2.120      .27                -30 dB                      -19.0 dB
%	   2.500                                                     -21.0 dB
%      3.0                                                       -23.5 dB                        
%	   3.384      .0864              -40 dB                      -26   dB
%      4.0                                                       -30   dB
%	   4.538      .0274              -50 dB                      -33.5 dB
%	   5.658      .00868             -60 dB                      -41   dB
%	   6.764      .00275             -70 dB                      -49   dB
%	   7.865      .000868            -80 dB                      -57.5 dB
%	   8.960      .000275            -90 dB                      -66   dB
%	***********************************************************************
%
%
%	To obtain a Kaiser window that designs a FIR filter with
%	of a given stop band ripple alpha [dB], calculate betaValue as follows:
%
%	       
%       0.1102*(alpha-8.7),                                       for alpha > 50 dB
%	betaValue = 0.5842*(alpha-21)^0.4 + 0.07886*(alpha-21),   for 50 >= alpha >= 21 dB
%       0,                                                        for alpha < 21 dB
%	
%	See ref. [2].
%        
%	Increasing betaValue widens the main lobe and decreases the
%	amplitude of the sidelobes (increases the attenuation).
%
%
%	References:
%		[1] Oppenheim, A.V., and R.W. Schafer,
%	        Zeitdiskrete Signalverarbeitung, 3. durchges. Auflage,
%	        Muenchen, Wien: Oldenbourg, 1999
%
%	    [2] http://www.mathworks.com/access/helpdesk/help/toolbox/signal/kaiser.html
%	        
%
%	Created:    2004/2005   Othmar Frey <ofrey@geo.unizh.ch>	
%	Modified:   24.11.2005  Othmar Frey <ofrey@geo.unizh.ch>
%
%   Copyright: 2017 Gamma Remote Sensing AG
%              Othmar Frey <frey@gamma-rs.ch>
%


% Parse the inputs
error(nargchk(3,4,nargin));
if(nargin==3)
	% if no betaValue is given, it is set to the default value 2.120
	betaValue = 2.120;
end
if(arraySize <= 0 || rem(arraySize,round(arraySize)))
	error(sprintf('arraySize must be a positive integer value. (Given: arraySize = %f)',arraySize));
end
if (BandWidthFraction <= 0.0 || BandWidthFraction > 1.0)
	error(sprintf('Given BandWidthFraction = %f. Required: BandWidthFraction in [0.0,1.0]',BandWidthFraction));
end
if (betaValue < 0.0)
	error(sprintf('Given betaValue coeffient = %f. Required: betaValue >= 0.0',betaValue));
end

n = arraySize;
n_frac = round(BandWidthFraction*n);
if(~rem(n,2)) % if n is even
	if(rem(n_frac,2))	% if n_frac is odd
    	n_frac = n_frac-1;
    end
end
if(rem(n,2)) % if n is odd    
	if(~rem(n_frac,2))	% if n_frac is even
    	n_frac = n_frac-1;
    end
end

fracwindow = kaiser(n_frac,betaValue);
kaiserWin = [zeros((n-n_frac)/2,1);fracwindow;zeros((n-n_frac)/2,1)];

kaiserWin = circshift(kaiserWin,DisplacementFromCenterOfArray);



%D = (attenuation - 7.95)/(2*pi*2.285)
%df = abs(freq2 - freq1);
%L = D/df + 1;
%%		freq1: passband cutoff freq (NORMALIZED)
%% 		freq2: stopband cutoff freq (NORMALIZED)
%%		dev1: passband ripple (DESIRED)
%%		dev2: stopband attenuation (not in dB)
%
%L % = filter Length (# of samples)
%
%D = (L-1)*df 
%
%attenuation = D*(2*pi*2.285)+7.95;
%
%
%
%
%[W1,f] = freqz(fracwindow/sum(fracwindow),1,512,2);
%plot(f,20*log10(abs(W1))); grid;
%legend(sprintf('beta = %f',betaValue),3)
