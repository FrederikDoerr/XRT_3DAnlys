function [XRT_3DAnlys_Opts,XRT_ImgPrc_Opts] = XRT_3DAnlys_ParameterLoader(XRT_Opts,Batch3DAnlys_Iter,XRT_3DAnlys_Opts,XRT_ImgPrc_Opts)
%XRT_3DAnlys_ParameterLoader Load image processing parameter.
%   [XRT_3DAnlys_Opts,XRT_ImgPrc_Opts] = ...
% XRT_3DAnlys_ParameterLoader(XRT_Opts,Batch3DAnlys_Iter,XRT_3DAnlys_Opts,XRT_ImgPrc_Opts)
%   
%   XRT_Opts.ImgPrc_3DAnlys: table with image processing details
%   Batch3DAnlys_Iter: current analysis (table row)
%   XRT_3DAnlys_Opts.ImgPrc_Param: table with image processing parameters
%   > assigned to XRT_ImgPrc_Opts
% 
%   Used for XRT_3DAnlys package.
% 
%   Author
%   --------
%   Frederik Doerr, Aug 2020 (MATLAB R2019b)
%   frederik.doerr(at)strath.ac.uk | CMAC (http://www.cmac.ac.uk/)
%   Github: https://github.com/frederik-d/XRT_3DAnlys


XRT_3DAnlys_Opts.formatSpec_s_d_s = '%s,%d,%s\n';
XRT_3DAnlys_Opts.formatSpec_s_d = '%s,%d\n';
XRT_3DAnlys_Opts.formatSpec_s_s = '%s,%s\n';
XRT_3DAnlys_Opts.formatSpec_s_d_d = '%s,%d,%d\n';

XRT_3DAnlys_Opts.ExpShorthand = XRT_Opts.ImgPrc_3DAnlys.ExpShorthand{Batch3DAnlys_Iter};
XRT_3DAnlys_Opts.Operator = XRT_Opts.ImgPrc_3DAnlys.Operator{Batch3DAnlys_Iter};
XRT_3DAnlys_Opts.ImgFormat = XRT_Opts.ImgPrc_3DAnlys.ImgFormat{Batch3DAnlys_Iter};
XRT_3DAnlys_Opts.DataResizeFactor = XRT_Opts.ImgPrc_3DAnlys.DataResizeFactor(Batch3DAnlys_Iter);
XRT_3DAnlys_Opts.ImagePixelSize = XRT_Opts.ImgPrc_3DAnlys.ImagePixelSize(Batch3DAnlys_Iter); %(um/px)
XRT_3DAnlys_Opts.ImagePixelSize_rS = XRT_3DAnlys_Opts.ImagePixelSize*XRT_3DAnlys_Opts.DataResizeFactor; %(um/px)

XRT_3DAnlys_Opts.VolumeAnlys_On = XRT_Opts.ImgPrc_3DAnlys.VolumeAnlys_On(Batch3DAnlys_Iter);
XRT_3DAnlys_Opts.PorosAnlys3D_On = XRT_Opts.ImgPrc_3DAnlys.PorosAnlys3D_On(Batch3DAnlys_Iter);
XRT_3DAnlys_Opts.BinSegment3D_On = XRT_Opts.ImgPrc_3DAnlys.BinSegment3D_On(Batch3DAnlys_Iter);
XRT_3DAnlys_Opts.VolumeAnlys_Concave_Anlys = XRT_Opts.ImgPrc_3DAnlys.VolumeAnlys_Concave_Anlys_On(Batch3DAnlys_Iter);
XRT_3DAnlys_Opts.VolumeAnlys_3DDescriptors = XRT_Opts.ImgPrc_3DAnlys.VolumeAnlys_3DDescriptors(Batch3DAnlys_Iter);

XRT_3DAnlys_Opts.VolExport_SurfRend_Check = XRT_Opts.ImgPrc_3DAnlys.VolExport_SurfRend_Check(Batch3DAnlys_Iter);
XRT_3DAnlys_Opts.VolumeBuild_on = XRT_Opts.ImgPrc_3DAnlys.VolumeBuild_on(Batch3DAnlys_Iter);
XRT_3DAnlys_Opts.ROI_Reduction_Check = XRT_Opts.ImgPrc_3DAnlys.ROI_Reduction_Check(Batch3DAnlys_Iter);

%% Assign Image File Paths

XRT_3DAnlys_Opts.Input_FolderPath = XRT_Opts.ImgPrc_3DAnlys.Input_FolderPath{Batch3DAnlys_Iter};
XRT_3DAnlys_Opts.Input_FolderPath_ROI = XRT_Opts.ImgPrc_3DAnlys.Input_FolderPath_ROI{Batch3DAnlys_Iter};
XRT_3DAnlys_Opts.Output_FolderPath = XRT_Opts.ImgPrc_3DAnlys.Output_FolderPath{Batch3DAnlys_Iter};

% filepath defined? exists? Relative filepath?
if ~exist(XRT_3DAnlys_Opts.Input_FolderPath,'dir') ...
    && ~isempty(XRT_3DAnlys_Opts.Input_FolderPath)
    if exist(fullfile(XRT_Opts.filesPath_MAIN,XRT_3DAnlys_Opts.Input_FolderPath),'dir')
        XRT_3DAnlys_Opts.Input_FolderPath = fullfile(XRT_Opts.filesPath_MAIN,XRT_3DAnlys_Opts.Input_FolderPath);
        warning('Relative Input_FolderPath detected > absolute: %s',XRT_3DAnlys_Opts.Input_FolderPath)
    end
end

if ~exist(XRT_3DAnlys_Opts.Input_FolderPath_ROI,'dir') ...
        && ~isempty(XRT_3DAnlys_Opts.Input_FolderPath_ROI)
    if exist(fullfile(XRT_Opts.filesPath_MAIN,XRT_3DAnlys_Opts.Input_FolderPath_ROI),'dir')
        XRT_3DAnlys_Opts.Input_FolderPath_ROI = fullfile(XRT_Opts.filesPath_MAIN,XRT_3DAnlys_Opts.Input_FolderPath_ROI);
        warning('Relative Input_FolderPath_ROI detected > absolute: %s',XRT_3DAnlys_Opts.Input_FolderPath_ROI)
    end
end

if ~exist(XRT_3DAnlys_Opts.Output_FolderPath,'dir') ...
        && ~isempty(XRT_3DAnlys_Opts.Output_FolderPath)
    if exist(fullfile(XRT_Opts.filesPath_MAIN,XRT_3DAnlys_Opts.Output_FolderPath),'dir')       
        XRT_3DAnlys_Opts.Output_FolderPath = fullfile(XRT_Opts.filesPath_MAIN,XRT_3DAnlys_Opts.Output_FolderPath);
        warning('Relative Output_FolderPath detected > absolute: %s',XRT_3DAnlys_Opts.Output_FolderPath)
    end
end

%% Move to Output_FolderPath
XRT_3DAnlys_Opts.Output_FolderPath_name = sprintf('%s_XRT_3DAnlys',XRT_3DAnlys_Opts.ExpShorthand);
XRT_3DAnlys_Opts.Output_FolderPath = fullfile(XRT_3DAnlys_Opts.Output_FolderPath,XRT_3DAnlys_Opts.Output_FolderPath_name);

if ~exist(XRT_3DAnlys_Opts.Output_FolderPath,'dir')
    mkdir(XRT_3DAnlys_Opts.Output_FolderPath)
end
cd(XRT_3DAnlys_Opts.Output_FolderPath)



%% Load Image Processing Parameter
XRT_3DAnlys_Opts.ImgPrc_Param_REF = XRT_Opts.ImgPrc_3DAnlys.ImgPrc_Param_REF{Batch3DAnlys_Iter};
XRT_3DAnlys_Opts.ImgPrc_Param = readtable(fullfile(XRT_Opts.filesPath_ImgPrc_Parameter,XRT_Opts.ImgPrc_Parameter_name),'Sheet','ImgPrc_Param');
idx_ImgPrc_Param_REF = find(strcmp(XRT_3DAnlys_Opts.ImgPrc_Param_REF,XRT_3DAnlys_Opts.ImgPrc_Param.ImgPrc_Param_REF)==1,1,'last');

if isempty(idx_ImgPrc_Param_REF)
    warning('%s - %s - %s: No Match for ImgPrc_Param_REF in Sheet >> ImgPrc_Param<< \n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename())

    prompt = {'Enter Image Processing Reference:'};
    prompttitle = 'ImgPrc_Param_REF';
    num_lines = 1; 
    XRT_3DAnlys_Opts.ImgPrc_Param_REF  = inputdlg(prompt,prompttitle,num_lines);
    options.Interpreter='tex';
    XRT_3DAnlys_Opts.ImgPrc_Param_REF  = char(XRT_3DAnlys_Opts.ImgPrc_Param_REF);
    idx_ImgPrc_Param_REF = find(strcmp(XRT_3DAnlys_Opts.ImgPrc_Param_REF,XRT_Opts.ImgPrc_Param.ExpShorthand)==1,1,'last');
end


% 2D Image filtering / thresholding parameters
XRT_ImgPrc_Opts.Img_Polarity = XRT_3DAnlys_Opts.ImgPrc_Param.Img_Polarity(idx_ImgPrc_Param_REF);
XRT_ImgPrc_Opts.ImgPrc_bw_thresVal_Manual = XRT_3DAnlys_Opts.ImgPrc_Param.ImgPrc_bw_thresVal_Manual(idx_ImgPrc_Param_REF);
XRT_ImgPrc_Opts.ImgPrc_bw_thresMethod = XRT_3DAnlys_Opts.ImgPrc_Param.ImgPrc_bw_thresMethod{idx_ImgPrc_Param_REF};
XRT_ImgPrc_Opts.ImgPrc_imFilter = XRT_3DAnlys_Opts.ImgPrc_Param.ImgPrc_imFilter{idx_ImgPrc_Param_REF};
if isempty(XRT_ImgPrc_Opts.ImgPrc_imFilter) || ...
        any(isnan(XRT_ImgPrc_Opts.ImgPrc_imFilter)) || ...
        strcmpi(XRT_ImgPrc_Opts.ImgPrc_imFilter,'NaN')
    XRT_ImgPrc_Opts.ImgPrc_imFilter_On = false;
else
    XRT_ImgPrc_Opts.ImgPrc_imFilter_On = true;
end
XRT_ImgPrc_Opts.ImgPrc_imFilter_Size = XRT_3DAnlys_Opts.ImgPrc_Param.ImgPrc_imFilter_Size(idx_ImgPrc_Param_REF);
XRT_ImgPrc_Opts.ImgPrc_imFilter_Parameter_1 = XRT_3DAnlys_Opts.ImgPrc_Param.ImgPrc_imFilter_Parameter_1(idx_ImgPrc_Param_REF);
if exist('ImgPrc_imFilter_Size','var') && isempty(ImgPrc_imFilter_Size)
    XRT_ImgPrc_Opts.ImgPrc_imFilter_Size = 12;
    warning('%s - %s - %s: ImgPrc_imFilter_Size not defined (default: 12) \n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename())
end

% 2D Noise reduction parameters
XRT_ImgPrc_Opts.NoiseReduction2D.ClearBorder = XRT_3DAnlys_Opts.ImgPrc_Param.NoiseReduction2D_ClearBorder(idx_ImgPrc_Param_REF);
XRT_ImgPrc_Opts.NoiseReduction2D.WhiteNoise_Size = XRT_3DAnlys_Opts.ImgPrc_Param.NoiseReduction2D_WhiteNoise_Size(idx_ImgPrc_Param_REF);
XRT_ImgPrc_Opts.NoiseReduction2D.Close_Size = XRT_3DAnlys_Opts.ImgPrc_Param.NoiseReduction2D_Close_Size(idx_ImgPrc_Param_REF);
XRT_ImgPrc_Opts.NoiseReduction2D.Open_Size = XRT_3DAnlys_Opts.ImgPrc_Param.NoiseReduction2D_Open_Size(idx_ImgPrc_Param_REF);
XRT_ImgPrc_Opts.NoiseReduction2D.BlackNoise_Size = XRT_3DAnlys_Opts.ImgPrc_Param.NoiseReduction2D_BlackNoise_Size(idx_ImgPrc_Param_REF);

% 3D Noise reduction parameters
XRT_ImgPrc_Opts.NoiseReduction3D.Method = XRT_3DAnlys_Opts.ImgPrc_Param.NoiseReduction3D_Method{Batch3DAnlys_Iter};
XRT_ImgPrc_Opts.NoiseReduction3D.ClearBorder = XRT_3DAnlys_Opts.ImgPrc_Param.NoiseReduction3D_ClearBorder(idx_ImgPrc_Param_REF);
XRT_ImgPrc_Opts.NoiseReduction3D.Close_Size_V = XRT_3DAnlys_Opts.ImgPrc_Param.NoiseReduction3D_Close_Size_V(idx_ImgPrc_Param_REF);
XRT_ImgPrc_Opts.NoiseReduction3D.Open_Size_V = XRT_3DAnlys_Opts.ImgPrc_Param.NoiseReduction3D_Open_Size_V(idx_ImgPrc_Param_REF);   
XRT_ImgPrc_Opts.NoiseReduction3D.Close_Size_V_ROI = XRT_3DAnlys_Opts.ImgPrc_Param.NoiseReduction3D_Close_Size_V_ROI(idx_ImgPrc_Param_REF);
XRT_ImgPrc_Opts.NoiseReduction3D.Open_Size_V_ROI = XRT_3DAnlys_Opts.ImgPrc_Param.NoiseReduction3D_Open_Size_V_ROI(idx_ImgPrc_Param_REF);

if ~isempty(XRT_ImgPrc_Opts.NoiseReduction3D.Method) || ...
        ~isempty(XRT_ImgPrc_Opts.NoiseReduction3D.ClearBorder) || ...
        ~isempty(XRT_ImgPrc_Opts.NoiseReduction3D.Close_Size_V) || ...
        ~isempty(XRT_ImgPrc_Opts.NoiseReduction3D.Open_Size_V) || ...
        ~isempty(XRT_ImgPrc_Opts.NoiseReduction3D.Close_Size_V_ROI) || ...
        ~isempty(XRT_ImgPrc_Opts.NoiseReduction3D.Open_Size_V_ROI)  
    XRT_3DAnlys_Opts.NoiseReduction3D = true;
else
    XRT_3DAnlys_Opts.NoiseReduction3D = false;
end


% Additional NanoCT_3DAnlys_Options
XRT_ImgPrc_Opts.Watershed.imextendedmin_Size = XRT_3DAnlys_Opts.ImgPrc_Param.Watershed_imextendedmin_Size(idx_ImgPrc_Param_REF);
XRT_ImgPrc_Opts.Watershed.bwdist_method = XRT_3DAnlys_Opts.ImgPrc_Param.Watershed_bwdist_method{idx_ImgPrc_Param_REF};
XRT_ImgPrc_Opts.Watershed.Shrinkwrap_Size  = XRT_3DAnlys_Opts.ImgPrc_Param.Watershed_Shrinkwrap_Size(idx_ImgPrc_Param_REF);
XRT_ImgPrc_Opts.Watershed.Shrinkwrap_Method = XRT_3DAnlys_Opts.ImgPrc_Param.Watershed_Shrinkwrap_Method{idx_ImgPrc_Param_REF};



%% Get Image Data from XRT_3DAnlys_Opts.Input_FolderPath

if ~isfield(XRT_3DAnlys_Opts,'filelist_imgXRT') || isempty(XRT_3DAnlys_Opts.filelist_imgXRT)
    if isfield(XRT_3DAnlys_Opts,'ImgFormat') && ~isempty(XRT_3DAnlys_Opts.ImgFormat)
        XRT_3DAnlys_Opts.filelist_imgXRT = dir(fullfile(XRT_3DAnlys_Opts.Input_FolderPath,sprintf('*.%s',XRT_3DAnlys_Opts.ImgFormat)));
    else

        % If image 
        XRT_3DAnlys_Opts.filelist_imgXRT_tif = dir(fullfile(XRT_3DAnlys_Opts.Input_FolderPath,'*.tif'));
        XRT_3DAnlys_Opts.filelist_imgXRT_bmp = dir(fullfile(XRT_3DAnlys_Opts.Input_FolderPath,'*.bmp'));

        if isempty(XRT_3DAnlys_Opts.filelist_imgXRT_tif) && ~isempty(XRT_3DAnlys_Opts.filelist_imgXRT_bmp)
            XRT_3DAnlys_Opts.filelist_imgXRT = XRT_3DAnlys_Opts.filelist_imgXRT_bmp;
            fprintf('%s - %s: ImgFormat not defined (autodetect: bmp)!\n',XRT_3DAnlys_Opts.ExpShorthand,mfilename())
        elseif isempty(XRT_3DAnlys_Opts.filelist_imgXRT_bmp) && ~isempty(XRT_3DAnlys_Opts.filelist_imgXRT_tif)
            XRT_3DAnlys_Opts.filelist_imgXRT = XRT_3DAnlys_Opts.filelist_imgXRT_tif;
            fprintf('%s - %s: ImgFormat not defined (autodetect: tif)!\n',XRT_3DAnlys_Opts.ExpShorthand,mfilename())
        else
            if isempty(XRT_3DAnlys_Opts.filelist_imgXRT_tif) && isempty(XRT_3DAnlys_Opts.filelist_imgXRT_bmp)
                prompt = {'XRT_3DAnlys_Opts.Input_FolderPath','ImgFormat'};
                prompttitle_Nuc = 'IMG Input';
                num_lines = 2; 
                c = inputdlg(prompt,prompttitle_Nuc,num_lines);
                options.Interpreter='tex';
                XRT_3DAnlys_Opts.Input_FolderPath = c{1};
                XRT_3DAnlys_Opts.ImgFormat = c{2};

                XRT_3DAnlys_Opts.filelist_imgXRT_bmp = dir(fullfile(XRT_3DAnlys_Opts.Input_FolderPath,sprintf('*.%s',XRT_3DAnlys_Opts.ImgFormat)));

            elseif length(XRT_3DAnlys_Opts.filelist_imgXRT_tif) > length(XRT_3DAnlys_Opts.filelist_imgXRT_bmp)
                XRT_3DAnlys_Opts.filelist_imgXRT = XRT_3DAnlys_Opts.filelist_imgXRT_tif;
                fprintf('%s - %s: ImgFormat not defined (autodetect: tif)!\n',XRT_3DAnlys_Opts.ExpShorthand,mfilename())
            else
                XRT_3DAnlys_Opts.filelist_imgXRT = XRT_3DAnlys_Opts.filelist_imgXRT_bmp;
                fprintf('%s - %s: ImgFormat not defined (autodetect: bmp)!\n',XRT_3DAnlys_Opts.ExpShorthand,mfilename())
            end
        end
    end
end


if ~isempty(XRT_3DAnlys_Opts.Input_FolderPath_ROI)
    XRT_3DAnlys_Opts.Input_FolderPath_ROI_Check = true;
    XRT_3DAnlys_Opts.Input_FolderPath_ROI = char(XRT_3DAnlys_Opts.Input_FolderPath_ROI);
    XRT_ImgPrc_Opts.ROI_Method = '';
    XRT_ImgPrc_Opts.ROI_Method_ShrinkWrap_Size = nan;
    XRT_ImgPrc_Opts.ROI_Method_ShrinkWrap_strel = nan;
    XRT_3DAnlys_Opts.convexhull_2d_check = false;
else
    XRT_3DAnlys_Opts.Input_FolderPath_ROI_Check = false;
    XRT_ImgPrc_Opts.ROI_Method = char(XRT_3DAnlys_Opts.ImgPrc_Param.ROI_Method(Batch3DAnlys_Iter));
    XRT_ImgPrc_Opts.ROI_Method_ShrinkWrap_Size = XRT_3DAnlys_Opts.ImgPrc_Param.ROI_Method_ShrinkWrap_Size(Batch3DAnlys_Iter);
    XRT_ImgPrc_Opts.ROI_Method_ShrinkWrap_strel = char(XRT_3DAnlys_Opts.ImgPrc_Param.ROI_Method_ShrinkWrap_strel(Batch3DAnlys_Iter));
    if strcmpi(XRT_ImgPrc_Opts.ROI_Method,'convexhull_2d')
        XRT_ImgPrc_Opts.convexhull_2d_check = true;
    else
        XRT_ImgPrc_Opts.convexhull_2d_check = false;
    end
end

if XRT_3DAnlys_Opts.Input_FolderPath_ROI_Check
    % Get ROI Images from XRT_3DAnlys_Opts.Input_FolderPath
    XRT_3DAnlys_Opts.filelist_imgXRT_ROI = dir(fullfile(XRT_3DAnlys_Opts.Input_FolderPath_ROI,'*.tif'));
    if isempty(XRT_3DAnlys_Opts.filelist_imgXRT_ROI)
        XRT_3DAnlys_Opts.filelist_imgXRT_ROI = dir(fullfile(XRT_3DAnlys_Opts.Input_FolderPath_ROI,'*.bmp'));
    end
end

