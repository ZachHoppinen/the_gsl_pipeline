function  gcf = plotPhaseDifferenceHSVScheme(complexCoherence,refSLC,pixCoordRefDim1,pixCoordRefDim2,extDim1,extDim2,axLabelDim1,axLabelDim2,sizeOfFont,dBThreshold)
% PLOTPHASEDIFFERENCEHSVScheme plots the 2-D map of phase differences
% contained in the 2-D array of complex coherences relative to a phase
% difference at reference location specified by the pixel coordinates
% pixCoordRefDim1, pixCoordRefDim2.
% The complexCoherence has to be spatially filtered adequately beforehand
% when it is calculated.
%
% The function can be used to plot 2-D maps of interferometric or
% polarimetric phase difference depending on the input data.
% 
% The phase of the complex coherence is combined with intensity in dB as
% a HSV color scheme    
%
%     Usage:
%       gcf = plotPhaseDifferenceHSVScheme(complexCoherence,refSLC,pixCoordRefDim1,pixCoordRefDim2,extDim1,extDim2,axLabelDim1,axLabelDim2,sizeOfFont,dBThreshold)
%
%     Input:
%       complexCoherence :   2-D array of complex coherences
%       refSLC           :   2-D array with SLC or MLI for intensity parameter
%       pixCoordRefDim1  :   pixel coord. of renference point
%       pixCoordRefDim2  :   " "
%       extDim1          :   1-D array with x-axis coordinates
%       extDim2          :   1-D array with y-axis coordinates
%       axLabelDim1      :   string with axis label for x-axis
%       axLabelDim2      :   string with axis label for x-axis
%       sizeOfFont       :
%       dBThreshold      :
%                 
%     Output:
%       gcf :  figure plot handle
%     
%     Last modified: 17.Sep.2019
%     Copyright: 2019 Gamma Remote Sensing AG
%                Othmar Frey, frey@gamma-rs.ch
%

    % Allocate memory for HSV (3-D) array
    [size1, size2] = size(complexCoherence);
    HSV = zeros(size1,size2,3);

    % Phase difference
    HSV(:,:,1) = angle(complexCoherence.*conj(complexCoherence(pixCoordRefDim1,pixCoordRefDim2)))./(2*pi) + 0.5;
    % Hue : color varies from red through yellow, green, cyan, blue, and magenta, and returns to red 

    HSV(:,:,2) = ones(size(HSV(:,:,1)));
    % Saturation: (no variation)
    % when HSV(:,2) is 0, the colors are unsaturated (i.e., shades of gray).
    % when HSV(:,2) is 1, the colors are fully saturated (i.e., they contain no white component).

    intensity_dB = 20*log10(abs(refSLC));
    intensity_dB(intensity_dB<dBThreshold) =dBThreshold;
    HSV(:,:,3) = (intensity_dB - min(intensity_dB(:)))./(max(intensity_dB(:)) - min(intensity_dB(:)));
    %HSV(:,:,3) = sqrt(abs(refSLC)./max(abs(refSLC(:))));
    %Value (V): As HSV(:,:,3) varies from 0 to 1, the brightness increases.

    % Convert from HSV to RGB scheme for plotting    
    rgb_image = hsv2rgb(HSV);
    
    % Plot the 2-D map of phase differences combined with intensity in dB as a HSV color scheme
    gcf = figure, image(extDim2,flipud(extDim1),rgb_image);
    colormap(hsv);
    cb_handle = colorbar;
    %keyboard
    set(cb_handle,'YTick',(0:9)*45/360);
    set(cb_handle,'YTickLabel',{'-180^o','-135^o','-90^o','-45^o','0^o','45^o','90^o','135^o','180^o'})
    axis xy;
    grid on;
    xlabel(axLabelDim2);
    ylabel(axLabelDim1);
    set(gca,'FontSize',sizeOfFont);
    set(get(gca,'XLabel'),'FontSize',sizeOfFont); 
    set(get(gca,'YLabel'),'FontSize',sizeOfFont);
    axis tight;
    hold on;
    plot(extDim2(pixCoordRefDim2),extDim1(pixCoordRefDim1),'r','MarkerSize',15,'Marker','+'); 
    hold off;
    colormap(hsv);
    drawnow;
end
