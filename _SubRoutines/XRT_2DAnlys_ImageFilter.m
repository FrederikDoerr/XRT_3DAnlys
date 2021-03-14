function [I,R] = XRT_2DAnlys_ImageFilter(I,Method_imFilter,Opt)
%XRT_2DAnlys_ImageFilter Apply 2D image filter.
%   [I,R] = XRT_2DAnlys_ImageFilter(I,Method_imFilter,Opt) function to ease 
%   the routine application of image filters. Select specific metod using
%   Method_imFilter:
%       'wiener2'
%       'imgaussfilt'
%       'medfilt2'
%       'kuwahara'
%       'imdiffusefilt'
%       'anisodiff2d'
%       'localcontrast'
%       'dncnn'
%
%   OPTIONS:
%      'ImgPrc_imFilter_Size'           Filter local neighbourhood size
%      'ImgPrc_imFilter_Parameter_1'    Flexible image filter parameter.
%                                       See method.
%      'Print_Log_On'       Enable fprintf outputs to monitor/record
%                           progress (default: false)
%      'Print_on'           Print to paramter file.
%      'ExpShorthand'       Sample ID (default: '')
%      'AppShorthand'       Application ID (default: mfilename())
% 
%
%   Author
%   --------
%   Frederik Doerr, Aug 2020 (MATLAB R2019b)
%   frederik.doerr(at)strath.ac.uk | CMAC (http://www.cmac.ac.uk/)
%   Github: https://github.com/frederik-d/XRT_3DAnlys

if ~isfield(Opt,'ExpShorthand')
    Opt.ExpShorthand = '';
    Opt.Print_Str = Opt.ExpShorthand;
else
    Opt.Print_Str = sprintf('%s - ',Opt.ExpShorthand);
end

if ~isfield(Opt,'AppShorthand')
    Opt.AppShorthand = '';
else
    Opt.Print_Str = sprintf('%s%s - ',Opt.Print_Str,Opt.AppShorthand);
end

if ~isfield(Opt,'Print_Log_On')
    Opt.Print_Log_On = false;
end
Opt.Print_Str = sprintf('%s%s:',Opt.Print_Str,mfilename());


if ~isfield(Opt, 'Print_on') || isempty(Opt.Print_on)
    Opt.Print_on = false;
end

if ~isfield(Opt, 'ImgPrc_imFilter_Parameter_1') ...
        || isempty(Opt.ImgPrc_imFilter_Parameter_1) ...
        || isnan(Opt.ImgPrc_imFilter_Parameter_1)
    Opt.ImgPrc_imFilter_Parameter_on = false;
else
    Opt.ImgPrc_imFilter_Parameter_on = true;
end

    
    
numImg = size(I,3);

for k = 1:numImg
    
    switch lower(Method_imFilter)
        case {'wiener2'}
            if isfield(Opt,'ImgPrc_imFilter_Size') && ~isnan(Opt.ImgPrc_imFilter_Size)
                R.dim_wiener1 = Opt.ImgPrc_imFilter_Size;
                R.dim_wiener1 = Opt.ImgPrc_imFilter_Size;
            else
                R.dim_wiener1 = 10;
                R.dim_wiener1 = 10;
            end
            I(:,:,k) = wiener2(I(:,:,k),[R.dim_wiener1 R.dim_wiener2]);
            
            if Opt.Print_Log_On
                fprintf('%s wiener2 Filter, Size >> %2.0f <<\n',Opt.Print_Str,R.dim_wiener2)
            end
        case {'imgaussfilt'}
            if Opt.ImgPrc_imFilter_Parameter_on
                R.SIGMA = Opt.ImgPrc_imFilter_Parameter_1;
            else
                R.SIGMA = 1.2;
            end

            I(:,:,k) = imgaussfilt(I(:,:,k),R.SIGMA);

            if Opt.Print_Log_On
                fprintf('%s imgaussfilt Filter\n',Opt.Print_Str)
            end
        case {'medfilt2'}
            if isfield(Opt,'ImgPrc_imFilter_Size') && ~isnan(Opt.ImgPrc_imFilter_Size)
                R.dim_medfilt1 = Opt.ImgPrc_imFilter_Size;
                R.dim_medfilt2 = Opt.ImgPrc_imFilter_Size;
            else
                R.dim_medfilt1 = 10;
                R.dim_medfilt2 = 10;
            end

            I(:,:,k) = medfilt2(I(:,:,k),[R.dim_medfilt1 R.dim_medfilt2]);

            if Opt.Print_Log_On
                fprintf('%s medfilt2 Filter, Size >> %2.0f <<\n',Opt.Print_Str,R.dim_medfilt2)
            end
        case {'kuwahara'}
            addpath(fullfile(Opt.filesPath_FileExchange,'Kuwahara'))
            if isfield(Opt,'ImgPrc_imFilter_Size') && ~isnan(Opt.ImgPrc_imFilter_Size)
                R.k_Kuwahara = Opt.ImgPrc_imFilter_Size;
            else
                R.k_Kuwahara = 4;
            end
            if R.k_Kuwahara < 2
                R.k_Kuwahara = 2;
            end
            I_filt = im2double(I(:,:,k)); % convert to double precision
            I_filt = Kuwahara(I_filt,(4*R.k_Kuwahara+1));
            I(:,:,k) = im2uint8(I_filt);
            
            
            if Opt.Print_Log_On
                fprintf('%s Kuwahara Filter, Size >> %2.0f << \n',Opt.Print_Str,R.k_Kuwahara)
            end
        case {'imdiffusefilt'}    
            % https://www.mathworks.com/help/images/ref/imdiffusefilt.html
            I(:,:,k) = imdiffusefilt(I(:,:,k));
            
        case {'anisodiff2d'}
            addpath(fullfile(Opt.filesPath_FileExchange,'\anisodiff_Perona-Malik\anisodiff_Perona-Malik'))
            %       ARGUMENT DESCRIPTION:
            %               VOL      - gray scale volume data (MxNxP).
            %               NUM_ITER - number of iterations. 
            %               DELTA_T  - integration constant (0 <= delta_t <= 3/44).
            %                          Usually, due to numerical stability this 
            %                          parameter is set to its maximum value.
            %               KAPPA    - gradient modulus threshold that controls the conduction.
            %               OPTION   - conduction coefficient functions proposed by Perona & Malik:
            %                          1 - c(x,y,z,t) = exp(-(nablaI/kappa).^2),
            %                              privileges high-contrast edges over low-contrast ones. 
            %                          2 - c(x,y,z,t) = 1./(1 + (nablaI/kappa).^2),
            %                              privileges wide regions over smaller ones.
            %          VOXEL_SPACING - 3x1 vector column with the x, y and z dimensions of
            %                          the voxel (milimeters). In particular, only cubic and 
            %                          anisotropic voxels in the z-direction are considered. 
            %                          When dealing with DICOM images, the voxel spacing 
            %                          dimensions can be extracted using MATLAB's dicominfo(.).

            I_filt = im2double(I(:,:,k)); % convert to double precision
            R.num_iter = 12;
            R.delta_t = 1/7;
            R.kappa = 30; 
            R.option = 1;
            I_filt = anisodiff2D(I_filt,R.num_iter, ....
                R.delta_t,R.kappa,R.option);                         
            I(:,:,k) = im2uint8(I_filt);
                        
            if Opt.Print_Log_On
                fprintf('%s anisodiff2D Filter, num_iter >> %2.0f << \n',Opt.Print_Str,R.num_iter)
            end

        case {'localcontrast'}
            if Opt.ImgPrc_imFilter_Parameter_on
                R.edgeThreshold = Opt.ImgPrc_imFilter_Parameter_1;
            else
                R.edgeThreshold = 0.1;
            end
            R.amount = -1;
            I(:,:,k) = localcontrast(I(:,:,k), R.edgeThreshold, R.amount);
            
            if Opt.Print_Log_On
                fprintf('%s Local Contrast Filter (edgeThreshold %.2f, amount %.2f)\n', ...
                    Opt.Print_Str,R.edgeThreshold,R.amount)
            end
        case {'dncnn'}
            % https://uk.mathworks.com/help/images/deep-learning.html
            % https://uk.mathworks.com/help/images/ref/denoiseimage.html
            R.net = denoisingNetwork('DnCNN');
            R.layers = dnCNNLayers;
            I(:,:,k) = denoiseImage(I(:,:,k), net);
            
            if Opt.Print_Log_On
                fprintf('%s DnCNN Filter\n',Opt.Print_Str)
            end
        otherwise
            if Opt.Print_Log_On
                warning('%s No Image Filter Matched\n',Opt.Print_Str)
            end
    end
end

if ~exist('R','var')
    R = struct;
end