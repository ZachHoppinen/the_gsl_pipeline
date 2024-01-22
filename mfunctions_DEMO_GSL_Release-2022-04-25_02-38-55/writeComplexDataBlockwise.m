function successful = writeComplexDataBlockwise(outfid, rc_data, NrOfRgPix, AziBlockSize, valtype)
%  	
%  	IMPORTANT !!! 	rc_data must be a complex matrix of dimension [NrOfRgPix,AziBlockSize],
%  			so that the first echo is stored in rc_data(:,1), the second echo
%  			in rc_data(:,2) etc.

	successful = 1;
	
	writeBlock = zeros(NrOfRgPix*2,AziBlockSize);
	
%  	Writing block of range compressed data to file
	real_index = 1:2:(NrOfRgPix*2 -1);
	imag_index = 2:2:NrOfRgPix*2;
	writeBlock(real_index,:) = real(rc_data);
	writeBlock(imag_index,:) = imag(rc_data);
	
	write_count = fwrite(outfid,writeBlock,valtype);
	if(write_count ~= NrOfRgPix*2*AziBlockSize)
		successful = 0; % i.e. NOT successful
		message = sprintf('Error: %d numbers written instead of %d !\n',write_count, NrOfRgPix*2*AziBlockSize);
			return;
	end
	
	
	
