function[I] = XRT_2DAnlys_AddScalebar(I,ImagePixelSize,Opt)   
%XRT_2DAnlys_AddScalebar Add scalebar to image.
%   [I] = XRT_2DAnlys_AddScalebar(I,ImagePixelSize,Opt) function to burn 
%   a scalebar in the 2D image (I) with a minimum size of 200 x 200 px. 
%   Automatically tries to optimise scalebar size. ImagePixelSize in um/px.
%
%   OPTIONS:
%      'Print_Log_On'      Enable fprintf outputs to monitor/record
%                          progress (default: false)
%
%   Author
%   --------
%   Frederik Doerr, Aug 2020 (MATLAB R2019b)
%   frederik.doerr(at)strath.ac.uk | CMAC (http://www.cmac.ac.uk/)
%   Github: https://github.com/frederik-d/XRT_3DAnlys


if ~isfield(Opt,'Print_Log_On')
    Opt.Print_Log_On = false;
end

[numRow,numCols,~] = size(I);
    
if size(I,1) > 200 || size(I,2) > 200
    
        %% Dimension Scale Bar
        row1 = floor(0.9*numRow);
        row2 = row1 + floor(numRow/150);
        column1 = floor(0.8*numCols); 
        column2 = column1 + floor(numCols/10);

        numPixel = length(column1:column2);
        ScaleSize = numPixel*ImagePixelSize; % scale in um

        % Round to next best 100rd
        ScaleSize_2 = round(ScaleSize,-2);
        numPixel_2 = round(ScaleSize / ImagePixelSize);

        column2 = column1 + numPixel_2;

        %% Write Scale Bar in Image
        if  islogical(I)
            %Background
            I(row1-10:row2+60, column1-20:column2+20) = 0;
            I(row1:row2, column1:column2) = 1;       
        elseif numel(size(I)) > 2
            I(row1-10:row2+60, column1-20:column2+20) = 255;
            I(row1:row2,column1:column2,:) = 255;
        else
            I(row1:row2, column1:column2) = 255; % Write white bar into image.
        end

    try
        I = im2uint8(I);
        text_str = sprintf('%.0f µm',ScaleSize_2);
        position = [column1,(row2+20)];
        FontSize_scaled = round(numCols/60);
        I = insertText(I,position,text_str,'FontSize',FontSize_scaled,'BoxColor','black','TextColor','white');
    catch 
        warning('%s: im2uint8 / insert scalebar failed \n',mfilename())
    end
    
    try
        FontSize_scaled = round(numCols/60);
        text_str_conv = sprintf('%.2f µm/px',round(ImagePixelSize,2));
        position = [floor(0.1*numCols),floor(0.05*numRow)];
        I = insertText(I,position,text_str_conv,'FontSize',FontSize_scaled,'BoxColor','black','TextColor','white');
    catch 
        warning('%s: Imagepixelsize String failed \n',mfilename())
    end
else
    if Opt.Print_Log_On
        warning('%s failed: Image smaller than specs (200px)\n',mfilename())
    end
end

    