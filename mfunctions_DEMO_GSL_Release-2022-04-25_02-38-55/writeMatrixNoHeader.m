function writeMatrixNoHeader(filename,datamatrix,valtype)
% WRITEMATRIX writes complex and non-complex raster data given
%   in the matrix 'datamatrix' WITHOUT any data header!!!
%   The data type is defined by 'valtype', which must contain
%   a Matlab-style data type (platform-independent version);
%   e.g. 'uint16', 'float32' etc.
%   Whether the data are complex or not is automatically
%   determined from the input matrix 'datamatrix'
%
%   Usage: 
%       writeMatrixNoHeader(filename,datamatrix,valtype)
%
%               'datamatrix has to be of the form:'
%               datamatrix(:,k) is the range echo at the k-th azimuth position
%               datamatrix(j,:)  is the azimuth line at the j-th range bin 
%
%   Created:    05. Dec 2005 by Othmar Frey <ofrey@geo.unizh.ch>
%
%   Copyright: 2017 Gamma Remote Sensing AG
%              Othmar Frey <frey@gamma-rs.ch>
%


% Checking number of input arguments
	error(nargchk(3,3,nargin));
	[NrOfRgPix NrOfEchoes] = size(datamatrix);

% Check whether input matrix datamatrix is complex or real. 
	if(isreal(datamatrix))
		isComplex = 0;
	else    % datamatrix is complex
		isComplex = 1;
	end

% Open raster file and determine the dimension and type from 32 byte header
	fid=fopen(filename,'wb','ieee-be');
	message = ferror(fid);
	if(~isempty(message))
		error(message);
		return;
    end

    if(isComplex)
        fprintf('     Writing complex data matrix of dimensions (%d x %d)\n',NrOfRgPix,NrOfEchoes);
    else
        fprintf('     Writing data matrix of dimensions (%d x %d)\n',NrOfRgPix,NrOfEchoes);
    end
    fprintf('     to the file %s .\n',filename);
	disp('     Byte order: "ieee-be" (MSBFirst)');
	fprintf('     Precision: %s .\n',valtype);
    
    if(isComplex)
        for k=1:NrOfEchoes
            curEcho = datamatrix(:,k);
            real_index = 1:2:(NrOfRgPix*2 -1);
            imag_index = 2:2:NrOfRgPix*2;
            writeEcho(real_index) = real(curEcho);
            writeEcho(imag_index) = imag(curEcho);
            write_count = fwrite(fid,writeEcho,valtype);
            if(write_count ~= NrOfRgPix*2)
                error(sprintf('Error: %d numbers written instead of %d !\n',write_count, NrOfRgPix*2));
                return;
            end
            fprintf('     %d echoes of %d written...\r',k,NrOfEchoes);
        end
        fprintf('     All echoes (%d) written to file %s\n',NrOfEchoes,filename)
        fclose(fid);
        return
    else
        write_count = fwrite(fid,datamatrix,valtype);
        if(write_count ~= NrOfRgPix*NrOfEchoes)
            error(sprintf('Error: %d numbers written instead of %d !\n',write_count, NrOfRgPix*NrOfEchoes));
            return;
        end
        fprintf('     All echoes (%d) written to file %s\n',NrOfEchoes,filename)
        fclose(fid);
        return
    end

