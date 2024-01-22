function partmatrix=readMatrixNoHeader(filename,rangeDim, aziDim, valtype,isComplex, firstPix,numPix,firstEcho,numEchoes)
% READMATRIXNOHEADER extracts a complex raster data set from a
%	binary data set.  NO header in front of the data!!!
%	readMatrixNoHeader extracts a part
%	firstPix:(firstPix+numPix)),(firstEcho:(firstEcho+numEchoes)
%	of a 2-D binary data file of dimensions:
%	rangeDim x aziDim.	
%
%	Usage:
%       partmatrix=readMatrixNoHeader(filename,rangeDim, aziDim, valtype,isComplex, firstPix,numPix,firstEcho,numEchoes)
%	
%	Created:	17. June 2008 by Othmar Frey <ofrey@geo.unizh.ch>
%   Modified:   04. Oct. 2011 by Othmar Frey <othmar.frey@gmail.com>
%
%
%  Modified:    25. Apr 2022 by Othmar Frey  <frey@gamma-rs.ch>
%               added
%                   elseif(strcmp(valtype,'int32'))
%	                    bytesPerPixel = 4;       
%                   elseif(strcmp(valtype,'uint32'))
%           	        bytesPerPixel = 4;
%
%   Copyright: 2022 Gamma Remote Sensing AG
%              Othmar Frey <frey@gamma-rs.ch>
%


% Checking number of input arguments
error(nargchk(9,9,nargin));

% Open raster file and determine the dimension and type from 32 byte header
fid=fopen(filename,'rb','ieee-be');
message = ferror(fid);
if(~isempty(message))
	error(message);
    return;
end

if(strcmp(valtype,'float32'))
	bytesPerPixel = 4;
elseif(strcmp(valtype,'double'))
	bytesPerPixel = 8;
elseif(strcmp(valtype,'uchar'))
	bytesPerPixel = 1;
elseif(strcmp(valtype,'char'))
	bytesPerPixel = 1;
elseif(strcmp(valtype,'int8'))
	bytesPerPixel = 1;
elseif(strcmp(valtype,'int16'))
	bytesPerPixel = 2;       
elseif(strcmp(valtype,'int32'))
	bytesPerPixel = 4;       
elseif(strcmp(valtype,'uint32'))
	bytesPerPixel = 4;
elseif(strcmp(valtype,'uint8'))
	bytesPerPixel = 1;  
elseif(strcmp(valtype,'short'))
	bytesPerPixel = 2;
elseif(strcmp(valtype,'ushort'))
	bytesPerPixel = 2;
else
	error('Unsupported data type')
end

% Allocate matrix for data
if(isComplex)
    partmatrix=complex(zeros(numPix,numEchoes));
    fact = 2;
else
    partmatrix=zeros(numPix,numEchoes);
    fact = 1;
end

% Skip all echos before firstEcho
status = fseek(fid, bytesPerPixel*fact*rangeDim*(firstEcho-1), 'bof');

if(isComplex)
    fprintf('     Extracting a complex data matrix of dimensions (%d x %d)\n',numPix,numEchoes);
else
    fprintf('     Extracting a data matrix of dimensions (%d x %d)\n',numPix,numEchoes);
end
fprintf('     from the file %s .\n',filename);
fprintf('     Dimension of complete file: (1:%d,1:%d)\n',rangeDim,aziDim);
fprintf('     Global indices of extracted part: (%d:%d,%d:%d)\n',firstPix,firstPix+numPix-1,firstEcho,firstEcho+numEchoes-1);
disp('     Byte order: "ieee-be" (MSBFirst)');
fprintf('     Precision: %s .\n',valtype);

if(isComplex)    
    for k=1:numEchoes
        [curEcho,read_count]=fread(fid,[1,rangeDim*2],valtype);
        if read_count<rangeDim*2
            error(['read only ',num2str(read_count),' numbers instead of ',...
                num2str(rangeDim*2)])
        end
        im_range=2*firstPix:2:2*(firstPix-1+numPix);
        re_range=im_range-1;
        partmatrix(:,k)=(curEcho(re_range)+curEcho(im_range).*j);
        %fprintf('     %d echoes of %d read...\r',k,numEchoes);
    end
    fprintf('     All echoes (%d) read.              \n',numEchoes);
else
    if(firstPix==1 && numPix==rangeDim)    
        [partmatrix,read_count] = fread(fid,[numPix,numEchoes],valtype);
        if read_count<(numPix*numEchoes)
            error(['read only ',num2str(read_count),' numbers instead of ',...
                num2str(numPix*numEchoes)])
        end
        fprintf('     All echoes (%d) read.              \n',numEchoes);
    else
        for k=1:numEchoes
            [curEcho,read_count]=fread(fid,[1,rangeDim],valtype);
            if read_count<rangeDim
                error(['read only ',num2str(read_count),' numbers instead of ',...
                    num2str(rangeDim)])
            end
		partmatrix(:,k) = curEcho(firstPix:(firstPix-1+numPix));
            %fprintf('     %d echoes of %d read...\r',k,numEchoes);
        end
        fprintf('     All echoes (%d) read.              \n',numEchoes);
    end
end

fclose(fid);
return