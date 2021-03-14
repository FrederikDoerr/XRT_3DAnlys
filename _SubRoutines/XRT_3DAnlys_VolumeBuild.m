%------------------------------------------------------------------------------------------------
% Code written by Frederik Doerr, Aug 2020 (MATLAB R2019b)
% Application: XRT Image and Data Analysis Framework (SubRoutine
% XRT_3DAnlys_VolumeBuild)
% 
% https://github.com/frederik-d
% Contact: frederik.doerr(at)strath.ac.uk | CMAC (http://www.cmac.ac.uk/)


%% Check if ImageStack pre-exists in Output_FolderPath
XRT_3DAnlys_Opts.Image_FolderPath_logical = fullfile(XRT_3DAnlys_Opts.Output_FolderPath,'XRT_3DAnlys_ImgStack_V_logical'); 
XRT_3DAnlys_Opts.Image_FolderPath_uint8 = fullfile(XRT_3DAnlys_Opts.Output_FolderPath,'XRT_3DAnlys_ImgStack_V_uint8'); 

if ~isfield(XRT_3DAnlys_Opts,'ROI_Reduction_Check')
    XRT_3DAnlys_Opts.ROI_Reduction_Check = false;
end

% Check if ROI_Method is defined
if isfield(XRT_3DAnlys_Opts,'ROI_Method') && ( ...
        strcmpi(XRT_3DAnlys_Opts.ROI_Method,'convexhull_2d') || ...
        strcmpi(XRT_3DAnlys_Opts.ROI_Method,'shrinkwrap_op') || ...
        strcmpi(XRT_3DAnlys_Opts.ROI_Method,'shrinkwrap'))
    XRT_3DAnlys_Opts.ROI_Method_Check = true;
elseif XRT_3DAnlys_Opts.Input_FolderPath_ROI_Check    
    XRT_3DAnlys_Opts.ROI_Method_Check = false;
    fprintf('%s - %s - %s: No ROI_Method >> Load from Directory << \n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename())
    fprintf('%s - %s - %s: ROI Folder >> %s \n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename(),XRT_3DAnlys_Opts.Input_FolderPath_ROI)
else
    XRT_3DAnlys_Opts.ROI_Method_Check = false;
    fprintf('%s - %s - %s: No Match for ROI_Method >> bw_ROI << \n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename())
end

clearvars V

if isfield(XRT_3DAnlys_Opts,'VolumeBuild_on') && XRT_3DAnlys_Opts.VolumeBuild_on % Forced VolumeBuild
    VolumeBuild_on = true;
    VolumeLoad_on = false;
elseif (~exist('V','var') && ~exist(XRT_3DAnlys_Opts.Image_FolderPath_logical,'dir') == 7) && V_binarize_Check == true
    VolumeBuild_on = true;
    VolumeLoad_on = false;
elseif (~exist('V','var') && ~exist(XRT_3DAnlys_Opts.Image_FolderPath_logical,'dir') == 7) && V_binarize_Check == true
    VolumeBuild_on = true;
    VolumeLoad_on = false;  
elseif (~exist('V','var') && exist(XRT_3DAnlys_Opts.Image_FolderPath_logical,'dir') == 7) && V_binarize_Check == true
    VolumeBuild_on = false;
    VolumeLoad_on = true;
    ImgStack_Load_Options.Image_FolderPath = fullfile(XRT_3DAnlys_Opts.Output_FolderPath,'XRT_3DAnlys_ImgStack_V_logical');
elseif (~exist('V','var') && exist(XRT_3DAnlys_Opts.Image_FolderPath_uint8,'dir') == 7) && V_binarize_Check == false
    VolumeBuild_on = false;
    VolumeLoad_on = true;
    ImgStack_Load_Options.Image_FolderPath = fullfile(XRT_3DAnlys_Opts.Output_FolderPath,'XRT_3DAnlys_ImgStack_V_uint8');
else 
    VolumeBuild_on = true;
    VolumeLoad_on = false;
end

%% Run XRT_3DAnlys_VolumeBuild

if VolumeBuild_on
    %% Create empty ImageV_binary 
    filelist_imgXRT = XRT_3DAnlys_Opts.filelist_imgXRT;
    XRT_3DAnlys_Opts.idx_Img_c = round(length(filelist_imgXRT)/2);    
    I_init_path = fullfile(filelist_imgXRT(XRT_3DAnlys_Opts.idx_Img_c).folder,filelist_imgXRT(XRT_3DAnlys_Opts.idx_Img_c).name);    
    I_init_info = imfinfo(I_init_path);
    I_init_info.BitDepth;

    if I_init_info.BitDepth > 8
        % Images need to be converted to 8bit prior to image procssing
        clear ImgConverter_Options
        ImgConverter_Options.ExpShorthand = {XRT_3DAnlys_Opts.ExpShorthand};
        ImgConverter_Options.mainFolder_FullPath = XRT_3DAnlys_Opts.Input_FolderPath(Batch3DAnlys_Iter);
        ImgConverter_Options.mainFolder_Output_FullPath = {fullfile(char(XRT_3DAnlys_Opts.Output_FolderPath(Batch3DAnlys_Iter)),sprintf('%s_XRT_3DAnlys',XRT_3DAnlys_Opts.ExpShorthand))};
        ImgConverter_Options.Prc_Method = 'imagestack';
        ImgConverter_Results = ImgConverter(ImgConverter_Options);
        
        % Update XRT_3DAnlys_Opts_3DAnlys_Input_FolderPath
        ImgConverter_Results.Image_Format = 'bmp';
        XRT_3DAnlys_Opts.Input_FolderPath(Batch3DAnlys_Iter) = ImgConverter_Results.mainFolder_FullPath_Output_List;
        filelist_imgXRT = dir(fullfile(char(XRT_3DAnlys_Opts.Input_FolderPath(Batch3DAnlys_Iter)),sprintf('*.%s',ImgConverter_Results.Image_Format)));
    end

    numImg = length(filelist_imgXRT);
    filelist_imgXRT_Resize = filelist_imgXRT(1:XRT_3DAnlys_Opts.DataResizeFactor:numImg);
    numImg_Resize = length(filelist_imgXRT_Resize);
    
    XRT_3DAnlys_Opts.idx_Img_c = round(numImg_Resize/2);
    I_init = imread(fullfile(filelist_imgXRT_Resize(XRT_3DAnlys_Opts.idx_Img_c).folder,filelist_imgXRT_Resize(XRT_3DAnlys_Opts.idx_Img_c).name));
    [mImage,nImage] = size(I_init);
    mImage = round(mImage/XRT_3DAnlys_Opts.DataResizeFactor);
    nImage = round(nImage/XRT_3DAnlys_Opts.DataResizeFactor);
    
    % Get a subset of frames to create image with representive grayscale 
    if ~islogical(I_init)
        % RC_Frac: number of frames used for comparison
        RC_Frac = 50;
        
        % RC_Frac equally spaced within filelist_imgXRT_Resize
        if 5*RC_Frac > length(filelist_imgXRT_Resize)
            RC_Frac = floor(length(filelist_imgXRT_Resize)/5);
        end
        numImg_RC_max = round(length(filelist_imgXRT_Resize)/RC_Frac);
        idx_RC = 2:1:numImg_RC_max-1;
        idx_RC = idx_RC.*RC_Frac;

        [rows_I_RC,col_I_RC] = size(imresize(I_init,[mImage nImage],'cubic'));
        
        % Copy RC_Frac images to empty image stack I_RC
        I_RC = zeros(rows_I_RC,col_I_RC,length(idx_RC));
        % Calculate image mean (foreground)
        I_RC_mean_list = zeros(1,length(idx_RC));
        for i = 1:length(idx_RC)
            I = imresize(imread(fullfile(filelist_imgXRT_Resize(idx_RC(i)).folder,filelist_imgXRT_Resize(idx_RC(i)).name)),[mImage nImage],'cubic');
            I_RC(:,:,i) = I;
            I_RC_mean_list(1,i) = mean(I(I~=0));
        end
        
        % Calculate image stack mean (foreground)
        I_RC_mean = mean(I_RC(I_RC~=0));
        if ~isnan(I_RC_mean)
            % Find image mean intensity with closest match to image stack mean 
            XRT_3DAnlys_Opts.idx_Img_c = idx_RC(min(abs(I_RC_mean_list - I_RC_mean))==abs(I_RC_mean_list - I_RC_mean));
            I_init = imresize(imread(fullfile(filelist_imgXRT(XRT_3DAnlys_Opts.idx_Img_c).folder,filelist_imgXRT(XRT_3DAnlys_Opts.idx_Img_c).name)),[mImage nImage],'cubic');
            fprintf('%s - %s - %s: Update I_init from Central Position to Grayscale Match (RC_Frac: %.0f) \n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename(),RC_Frac)
        end
    end
    
    I_init_info = imfinfo(I_init_path);
    I_init_info.BitDepth;
    I_init_resize = imresize(I_init,[mImage nImage],'cubic');
    
    % Load ROI if filelist_imgXRT_ROI specified
    if XRT_3DAnlys_Opts.Input_FolderPath_ROI_Check
        filelist_imgXRT_Resize_ROI = XRT_3DAnlys_Opts.filelist_imgXRT_ROI(1:XRT_3DAnlys_Opts.DataResizeFactor:numImg);
    else
        filelist_imgXRT_Resize_ROI = zeros(numImg_Resize);
    end
  
    
    
    %% Run XRT_2DPreProcessing
    XRT_3DAnlys_Opts.filesPath_2DImagePrc_Name = 'XRT_2DImgPrc';
    XRT_3DAnlys_Opts.filesPath_2DImagePrc = fullfile(XRT_3DAnlys_Opts.Output_FolderPath,XRT_3DAnlys_Opts.filesPath_2DImagePrc_Name);
    if ~exist(XRT_3DAnlys_Opts.filesPath_2DImagePrc,'dir')
        mkdir(XRT_3DAnlys_Opts.filesPath_2DImagePrc)
    end

    fileID = fopen(fullfile(XRT_3DAnlys_Opts.filesPath_2DImagePrc,sprintf('%s_ImgPrc_Parameter.txt',XRT_3DAnlys_Opts.ExpShorthand)),'wt');
    ImgProcessing_Parameter_header = { ...
                'ExpShorthand',(XRT_3DAnlys_Opts.ExpShorthand); ...
                'I_init',filelist_imgXRT_Resize(XRT_3DAnlys_Opts.idx_Img_c).name; ...
                };
    for i = 1:length(ImgProcessing_Parameter_header)
        fprintf(fileID,XRT_3DAnlys_Opts.formatSpec_s_s,ImgProcessing_Parameter_header{i,:});
    end
    fclose(fileID);
    
    
    if ~XRT_3DAnlys_Opts.Input_FolderPath_ROI_Check || V_binarize_Check
        if ~islogical(I_init)
           
            fprintf('%s - %s - %s: XRT_2DImgPrc\n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename())

            
            ImgName = sprintf('%s_0_I_init_resize.bmp',XRT_3DAnlys_Opts.ExpShorthand);
            imwrite(I_init_resize, fullfile(XRT_3DAnlys_Opts.filesPath_2DImagePrc,ImgName),'bmp');
            
            [Image_scaled] = XRT_2DAnlys_AddScalebar(I_init_resize,XRT_3DAnlys_Opts.ImagePixelSize_rS,XRT_3DAnlys_Opts);
            imwrite(Image_scaled, fullfile(XRT_3DAnlys_Opts.filesPath_2DImagePrc,sprintf('%s_0_I_init_resize_scaled.bmp',XRT_3DAnlys_Opts.ExpShorthand)),'bmp')
            
            if exist('Img_Polarity','var') && Img_Polarity == 0
                I_init_resize = imcomplement(I_init_resize);
                fileID = fopen(fullfile(XRT_3DAnlys_Opts.filesPath_2DImagePrc,sprintf('%s_ImgPrc_Parameter.txt',XRT_3DAnlys_Opts.ExpShorthand)),'at');
                fprintf(fileID,XRT_3DAnlys_Opts.formatSpec_s_s,'Img_Polarity', '0 - reversed');
                fclose(fileID);
            else
                fileID = fopen(fullfile(XRT_3DAnlys_Opts.filesPath_2DImagePrc,sprintf('%s_ImgPrc_Parameter.txt',XRT_3DAnlys_Opts.ExpShorthand)),'at');
                fprintf(fileID,XRT_3DAnlys_Opts.formatSpec_s_s,'Img_Polarity', '1 - raw');
                fclose(fileID);
            end
            
            fileID = fopen(fullfile(XRT_3DAnlys_Opts.filesPath_2DImagePrc,sprintf('%s_ImgPrc_Parameter.txt',XRT_3DAnlys_Opts.ExpShorthand)),'at');
            fprintf(fileID,'\n[imFilter2D]\r\n');
            if XRT_ImgPrc_Opts.ImgPrc_imFilter_On
                Image_Filter_Options = XRT_ImgPrc_Opts;
                Image_Filter_Options.ExpShorthand = XRT_3DAnlys_Opts.ExpShorthand;
                Image_Filter_Options.AppShorthand = XRT_3DAnlys_Opts.AppShorthand;
                Image_Filter_Options.filesPath_2DImagePrc = XRT_3DAnlys_Opts.filesPath_2DImagePrc;
                Image_Filter_Options.Control_on = true;
                Image_Filter_Options.Print_on = true;
                Image_Filter_Options.filesPath_FileExchange = XRT_Opts.filesPath_FileExchange;
                Image_Filter_Options.filesPath_2DImagePrc_FileName = sprintf('%s_ImgPrc_Parameter.txt',XRT_3DAnlys_Opts.ExpShorthand);
                
                [I_init_filt,XRT_2DImgPrc_Param] = XRT_2DAnlys_ImageFilter(I_init_resize,XRT_ImgPrc_Opts.ImgPrc_imFilter,Image_Filter_Options);
                Image_Filter_Options.Print_on = false;
                Image_Filter_Options.Control_on = false;
                
                fprintf(fileID,XRT_3DAnlys_Opts.formatSpec_s_s,'ImgPrc_imFilterMethod', XRT_ImgPrc_Opts.ImgPrc_imFilter);
                fprintf(fileID,XRT_3DAnlys_Opts.formatSpec_s_d,'ImgPrc_imFilter_Size', XRT_ImgPrc_Opts.ImgPrc_imFilter_Size);
                fprintf(fileID,XRT_3DAnlys_Opts.formatSpec_s_d,'ImgPrc_imFilter_Parameter_1', XRT_ImgPrc_Opts.ImgPrc_imFilter_Parameter_1);
            else
                I_init_filt = I_init_resize;
                XRT_2DImgPrc_Param = struct;
                fprintf(fileID,XRT_3DAnlys_Opts.formatSpec_s_s,'ImgPrc_imFilterMethod', 'None');
            end
            fclose(fileID);
            
            ImgName = sprintf('%s_1_I_init_filter.bmp',XRT_3DAnlys_Opts.ExpShorthand);
            imwrite(I_init_filt, fullfile(XRT_3DAnlys_Opts.filesPath_2DImagePrc,ImgName),'bmp');
            [Image_scaled] = XRT_2DAnlys_AddScalebar(I_init_filt,XRT_3DAnlys_Opts.ImagePixelSize_rS,XRT_3DAnlys_Opts);
            imwrite(Image_scaled, fullfile(XRT_3DAnlys_Opts.filesPath_2DImagePrc,sprintf('%s_1_I_init_filter_scaled.bmp',XRT_3DAnlys_Opts.ExpShorthand)),'bmp')
            
            fileID = fopen(fullfile(XRT_3DAnlys_Opts.filesPath_2DImagePrc,sprintf('%s_ImgPrc_Parameter.txt',XRT_3DAnlys_Opts.ExpShorthand)),'at');
            fprintf(fileID,'\n[Binarization2D]\n');
            switch lower(XRT_ImgPrc_Opts.ImgPrc_bw_thresMethod)
                case {'otsu'}
                    try
                        % Get a Subset of Frames to calculate Threshold
                        RC_Frac = length(filelist_imgXRT_Resize)*0.1;
                        numImg_RC_max = round(length(filelist_imgXRT_Resize)/RC_Frac);
                        idx_RC = 2:1:numImg_RC_max-1;
                        idx_RC = fix(idx_RC.*RC_Frac);
                        
                        [rows_I_RC,col_I_RC] = size(I_init_filt);
                        
                        I_RC = uint8(zeros(rows_I_RC,col_I_RC,length(idx_RC)));
                        for i = 1:length(idx_RC)
                            I_RC(:,:,i) = imresize(imread(fullfile(filelist_imgXRT_Resize(idx_RC(i)).folder,filelist_imgXRT_Resize(idx_RC(i)).name)),[mImage nImage],'cubic');
                        end
                        
                        Image_Filter_Options.Control_on = false;
                        Image_Filter_Options.Print_on = false;
                        Image_Filter_Options.filesPath_FileExchange = XRT_Opts.filesPath_FileExchange;
                        [I_RC,~] = XRT_2DAnlys_ImageFilter(I_RC,XRT_ImgPrc_Opts.ImgPrc_imFilter,Image_Filter_Options);
                    catch
                        warning('%s: RC SubStack FAILED (Multi Cross-Section Approach for Thresholding)\n',mfilename())
                        I_RC = I_init_filt;
                    end
                    [level,~] = graythresh(I_RC);
                    [~, grayLevels] = imhist(I_RC);
                    clearvars I_RC
                    XRT_2DImgPrc_Param.ImgPrc_bw_thresVal = round(grayLevels(end)*level);
                    XRT_2DImgPrc_Param.ImgPrc_bw_thresMethod = XRT_ImgPrc_Opts.ImgPrc_bw_thresMethod;
                    bw = I_init_filt > XRT_2DImgPrc_Param.ImgPrc_bw_thresVal;

                    
                    fprintf(fileID,XRT_3DAnlys_Opts.formatSpec_s_s,'ImgPrc_bw_thresMethod', XRT_2DImgPrc_Param.ImgPrc_bw_thresMethod);
                    fprintf(fileID,XRT_3DAnlys_Opts.formatSpec_s_d,'ImgPrc_bw_thresVal', XRT_2DImgPrc_Param.ImgPrc_bw_thresVal);
                    fclose(fileID);

                    fprintf('%s - %s - %s: Otsu Threshold >> %3.0f << \n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename(),XRT_2DImgPrc_Param.ImgPrc_bw_thresVal)
                case {'activecontour'}
                    I_mask = imbinarize(I_init_filt);
                    XRT_2DImgPrc_Param.bw_thresVal_Method_Iter = 300;
                    bw = activecontour(I_init_filt,I_mask,XRT_2DImgPrc_Param.bw_thresVal_Method_Iter);
                    
                    XRT_2DImgPrc_Param.ImgPrc_bw_thresMethod = XRT_ImgPrc_Opts.ImgPrc_bw_thresMethod;
                    
                    fprintf(fileID,XRT_3DAnlys_Opts.formatSpec_s_s,'ImgPrc_bw_thresMethod', XRT_2DImgPrc_Param.ImgPrc_bw_thresMethod);
                    fprintf(fileID,XRT_3DAnlys_Opts.formatSpec_s_d,'bw_thresVal_Method_Iter', XRT_2DImgPrc_Param.bw_thresVal_Method_Iter);
                    fclose(fileID);

                    fprintf('%s - %s - %s: Active Contour \n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename())
                case {'ridlercalvard'}
                    addpath(fullfile(XRT_Opts.filesPath_FileExchange,'RidlerCalvard_isodata'))
                    try
                        RC_Frac = length(filelist_imgXRT_Resize)*0.1;
                        numImg_RC_max = round(length(filelist_imgXRT_Resize)/RC_Frac);
                        idx_RC = 2:1:numImg_RC_max-1;
                        idx_RC = fix(idx_RC.*RC_Frac);
                        
                        [rows_I_RC,col_I_RC] = size(I_init_filt);
                        
                        I_RC = uint8(zeros(rows_I_RC,col_I_RC,length(idx_RC)));
                        for i = 1:length(idx_RC)
                            I_RC(:,:,i) = imresize(imread(fullfile(filelist_imgXRT_Resize(idx_RC(i)).folder,filelist_imgXRT_Resize(idx_RC(i)).name)),[mImage nImage],'cubic');
                        end
                        Image_Filter_Options.Control_on = false;
                        Image_Filter_Options.Print_on = false;
                        Image_Filter_Options.filesPath_FileExchange = XRT_Opts.filesPath_FileExchange;
                        [I_RC,~] = XRT_2DAnlys_ImageFilter(I_RC,XRT_ImgPrc_Opts.ImgPrc_imFilter,Image_Filter_Options);
                    catch
                        warning('%s: RC SubStack FAILED (Multi Cross-Section Approach for Thresholding)\n',mfilename())
                        I_RC = I_init_filt;
                    end

                    level = RidlerCalvard_isodata(I_RC);
                    [~, grayLevels] = imhist(I_RC);
                    clearvars I_RC
                    XRT_2DImgPrc_Param.ImgPrc_bw_thresVal = round(grayLevels(end)*level);
                    bw = I_init_filt > XRT_2DImgPrc_Param.ImgPrc_bw_thresVal;
                    XRT_2DImgPrc_Param.ImgPrc_bw_thresMethod = XRT_ImgPrc_Opts.ImgPrc_bw_thresMethod;
                    
                    fprintf(fileID,XRT_3DAnlys_Opts.formatSpec_s_s','ImgPrc_bw_thresMethod', XRT_2DImgPrc_Param.ImgPrc_bw_thresMethod);
                    fprintf(fileID,XRT_3DAnlys_Opts.formatSpec_s_s,'FileExchange', 'RidlerCalvard_isodata');
                    fprintf(fileID,XRT_3DAnlys_Opts.formatSpec_s_d,'ImgPrc_bw_thresVal', XRT_2DImgPrc_Param.ImgPrc_bw_thresVal);
                    fclose(fileID);

                    fprintf('%s - %s - %s: Ridler Calvard Threshold >> %3.0f << \n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename(),XRT_2DImgPrc_Param.ImgPrc_bw_thresVal)
                case {'adaptive'}
                    XRT_2DImgPrc_Param.adaptive_dim = 10;
                    XRT_2DImgPrc_Param.h = fspecial('average',[XRT_2DImgPrc_Param.adaptive_dim XRT_2DImgPrc_Param.adaptive_dim]);
                    bw = imfilter(I_init_filt, XRT_2DImgPrc_Param.h);
                    
                    XRT_2DImgPrc_Param.ImgPrc_bw_thresMethod = XRT_ImgPrc_Opts.ImgPrc_bw_thresMethod;

                    fprintf(fileID,XRT_3DAnlys_Opts.formatSpec_s_s,'ImgPrc_bw_thresMethod', XRT_2DImgPrc_Param.ImgPrc_bw_thresMethod);
                    fprintf(fileID,XRT_3DAnlys_Opts.formatSpec_s_d,'adaptive_dim', XRT_2DImgPrc_Param.adaptive_dim);
                    fclose(fileID);

                    fprintf('%s - %s - %s: Adaptive Threshold, Size >> %2.0f <<\n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename(),XRT_2DImgPrc_Param.adaptive_dim)
                otherwise
                    XRT_2DImgPrc_Param.ImgPrc_bw_thresMethod = 'Manual';
                    XRT_2DImgPrc_Param.ImgPrc_bw_thresVal = XRT_ImgPrc_Opts.ImgPrc_bw_thresVal_Manual;
                    XRT_2DImgPrc_Param.ImgPrc_bw_thresVal_Manual = true;
                    bw = I_init_filt > XRT_2DImgPrc_Param.ImgPrc_bw_thresVal;
                    
                    
                    fprintf(fileID,XRT_3DAnlys_Opts.formatSpec_s_d,'XRT_2DImgPrc_Param.ImgPrc_bw_thresVal_Manual', XRT_2DImgPrc_Param.ImgPrc_bw_thresVal_Manual);
                    fprintf(fileID,XRT_3DAnlys_Opts.formatSpec_s_d,'XRT_2DImgPrc_Param.ImgPrc_bw_thresVal', XRT_2DImgPrc_Param.ImgPrc_bw_thresVal);
                    fclose(fileID);

                    fprintf('%s - %s - %s: Manual Threshold >> %3.0f << \n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename(),XRT_2DImgPrc_Param.ImgPrc_bw_thresVal)
            end
            

            ImgName = sprintf('%s_2_I_init_bw.bmp',XRT_3DAnlys_Opts.ExpShorthand);
            imwrite(bw, fullfile(XRT_3DAnlys_Opts.filesPath_2DImagePrc,ImgName),'bmp');
            [Image_scaled] = XRT_2DAnlys_AddScalebar(bw,XRT_3DAnlys_Opts.ImagePixelSize_rS,XRT_3DAnlys_Opts);
            imwrite(Image_scaled, fullfile(XRT_3DAnlys_Opts.filesPath_2DImagePrc,sprintf('%s_2_I_init_bw_scaled.bmp',XRT_3DAnlys_Opts.ExpShorthand)),'bmp')
                                                
            binaryImage_denoised = XRT_2DAnlys_NoiseReduction(bw,XRT_ImgPrc_Opts.NoiseReduction2D);
            
            fprintf('%s - %s - %s: NoiseReduction2D Reference >> %s <<\n', ...
                XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename(),XRT_3DAnlys_Opts.ImgPrc_Param_REF)
            fprintf('%s - %s - %s: NoiseReduction2D_ClearBorder >> %.0f <<\n', ...
                XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename(),XRT_ImgPrc_Opts.NoiseReduction2D.ClearBorder)
            fprintf('%s - %s - %s: NoiseReduction2D_WhiteNoise_Size >> %.0f <<\n', ...
                XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename(),XRT_ImgPrc_Opts.NoiseReduction2D.WhiteNoise_Size)
            fprintf('%s - %s - %s: NoiseReduction2D_Close_Size >> %.0f <<\n', ...
                XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename(),XRT_ImgPrc_Opts.NoiseReduction2D.Close_Size)
            fprintf('%s - %s - %s: NoiseReduction2D_Open_Size >> %.0f <<\n', ...
                XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename(),XRT_ImgPrc_Opts.NoiseReduction2D.Open_Size)
            fprintf('%s - %s - %s: NoiseReduction2D_BlackNoise_Size >> %.0f <<\n', ...
                XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename(),XRT_ImgPrc_Opts.NoiseReduction2D.BlackNoise_Size)
            
            ImgName = sprintf('%s_3_I_init_bw_denoised.bmp',XRT_3DAnlys_Opts.ExpShorthand);
            imwrite(binaryImage_denoised, fullfile(XRT_3DAnlys_Opts.filesPath_2DImagePrc,ImgName),'bmp');                
            [Image_scaled] = XRT_2DAnlys_AddScalebar(binaryImage_denoised,XRT_3DAnlys_Opts.ImagePixelSize_rS,XRT_3DAnlys_Opts);
            imwrite(Image_scaled, fullfile(XRT_3DAnlys_Opts.filesPath_2DImagePrc,sprintf('%s_3_I_init_bw_denoised_scaled.bmp',XRT_3DAnlys_Opts.ExpShorthand)),'bmp')
            
            % Evaluate Signal-to-Noise Ratio
            I_init = imresize(imread(fullfile(filelist_imgXRT(XRT_3DAnlys_Opts.idx_Img_c).folder,filelist_imgXRT(XRT_3DAnlys_Opts.idx_Img_c).name)),[mImage nImage],'cubic'); 
            I_init_resize = imresize(I_init,[mImage nImage],'cubic');
            
            I_init_resize_mean = mean(I_init_resize(binaryImage_denoised));
            I_init_resize_std = std(double(I_init_resize(binaryImage_denoised)));
            
            I_init_filt_mean = mean(I_init_filt(binaryImage_denoised));
            I_init_filt_std = std(double(I_init_filt(binaryImage_denoised)));
        end
        
    else
        % All Image Data Supplied: No Pre-processing Required        
        Img_Polarity = XRT_Opts.Img_Polarity(idxExp);
    end
    
    
    fileID = fopen(fullfile(XRT_3DAnlys_Opts.filesPath_2DImagePrc,sprintf('%s_ImgPrc_Parameter.txt',XRT_3DAnlys_Opts.ExpShorthand)),'at');
    fprintf(fileID,'\n[NoiseReduction2D]\n');
        ImgProcessing_Parameter_dat = { ...
                'ImgPrc_Param_REF',XRT_3DAnlys_Opts.ImgPrc_Param_REF; ...
                'NoiseReduction2D_ClearBorder',XRT_ImgPrc_Opts.NoiseReduction2D.ClearBorder; ...
                'NoiseReduction2D_WhiteNoise_Size',XRT_ImgPrc_Opts.NoiseReduction2D.WhiteNoise_Size; ...
                'NoiseReduction2D_Close_Size',XRT_ImgPrc_Opts.NoiseReduction2D.Close_Size; ...
                'NoiseReduction2D_Open_Size',XRT_ImgPrc_Opts.NoiseReduction2D.Open_Size; ...
                'NoiseReduction2D_BlackNoise_Size',XRT_ImgPrc_Opts.NoiseReduction2D.BlackNoise_Size; ...
                };
    for i = 1:length(ImgProcessing_Parameter_dat)
        fprintf(fileID,XRT_3DAnlys_Opts.formatSpec_s_d,ImgProcessing_Parameter_dat{i,:});
    end  
    fclose(fileID);
    
    if exist('I_init_resize_mean','var')
        fileID = fopen(fullfile(XRT_3DAnlys_Opts.filesPath_2DImagePrc,sprintf('%s_ImgPrc_Parameter.txt',XRT_3DAnlys_Opts.ExpShorthand)),'at');
        fprintf(fileID,'\n[Noise Analysis (I_init)]\n');
        fprintf(fileID,XRT_3DAnlys_Opts.formatSpec_s_d_d,'I_init_resize (Int SNR +/- Std)',I_init_resize_mean,I_init_resize_std);
        fprintf(fileID,XRT_3DAnlys_Opts.formatSpec_s_d_d,'I_init_resize (Int SNR +/- Std)',I_init_resize_mean/I_init_resize_std);
        fprintf(fileID,XRT_3DAnlys_Opts.formatSpec_s_d_d,'I_init_filt (Int SNR +/- Std)',I_init_filt_mean,I_init_filt_std);
        fprintf(fileID,XRT_3DAnlys_Opts.formatSpec_s_d_d,'I_init_resize (Int SNR +/- Std)',I_init_resize_mean/I_init_resize_std);
        fclose(fileID);
    end
    
    
    %% Build Volume % % % % % % % % % % % % % % % % % % %     
    % Create empty volume for parallel 3DVB: V uint8 or logical
    if ~V_binarize_Check
        V = zeros(mImage,nImage,numImg_Resize,'uint8');
        V_ROI = zeros(mImage,nImage,numImg_Resize,'logical');
    else
        V = zeros(mImage,nImage,numImg_Resize,'logical');
        V_ROI = zeros(mImage,nImage,numImg_Resize,'logical');
    end

    % Determine size of image stack
    I_info = imfinfo(I_init_path);
    image_bytes = I_info.Width * I_info.Height * I_info.BitDepth;
    ImageStack_Size = image_bytes * numImg_Resize * 10^(-9);
    
    % Determine number of images per VolChunk at fixed VolChunk size
    num_VolChunks = floor(ImageStack_Size / 5); 
    if num_VolChunks < 1
        num_VolChunks = 1;
    end
    num_VolChunks_nImg = floor(numImg_Resize/num_VolChunks);
    VolChunk_ImgID_Array = 1:num_VolChunks_nImg-1:num_VolChunks*num_VolChunks_nImg;
    if (numImg_Resize - num_VolChunks * num_VolChunks_nImg) > 0
        VolChunk_ImgID_Array = [VolChunk_ImgID_Array,numImg_Resize];
    end
    
    
    % Start of parallel 3DVB
    for iter = 1:num_VolChunks

        VolChunk_Start = VolChunk_ImgID_Array(iter);
        VolChunk_End = VolChunk_ImgID_Array(iter+1);
        numImg_VolChunks = VolChunk_End - VolChunk_Start+1;

        if ~V_binarize_Check
            ImageV_binary = zeros(mImage,nImage,numImg_VolChunks,'uint8');
            ImageV_binary_ROI = zeros(mImage,nImage,numImg_VolChunks,'logical');
        else
            ImageV_binary = zeros(mImage,nImage,numImg_VolChunks,'logical');
            ImageV_binary_ROI = zeros(mImage,nImage,numImg_VolChunks,'logical');
        end

        filelist_imgXRT_VolChunk = filelist_imgXRT_Resize(VolChunk_Start:VolChunk_End);
        numImg_Resize_VolChunk = length(filelist_imgXRT_VolChunk);
        
        if XRT_3DAnlys_Opts.Input_FolderPath_ROI_Check
           filelist_imgXRT_ROI_VolChunk = filelist_imgXRT_Resize_ROI(VolChunk_Start:VolChunk_End);
        else
            filelist_imgXRT_ROI_VolChunk = filelist_imgXRT_VolChunk;
        end

        parfor k = 1:length(filelist_imgXRT_VolChunk)
            I_tmp = imread(fullfile(filelist_imgXRT_VolChunk(k).folder,filelist_imgXRT_VolChunk(k).name));

            if XRT_3DAnlys_Opts.DataResizeFactor > 1
                % Resize image
                I_tmp = imresize(I_tmp,[mImage nImage],'cubic');
            end
            
            if ~XRT_ImgPrc_Opts.Img_Polarity
                I_tmp = imcomplement(I_tmp);
            end
            
            if numel(size(I_tmp)) > 2
                I_tmp = rgb2gray(I_tmp);
            end

            % Apply Image Fiter and Binarisation
            if ~XRT_3DAnlys_Opts.Input_FolderPath_ROI_Check || V_binarize_Check
                if ~islogical(I_init)
                    if XRT_ImgPrc_Opts.ImgPrc_imFilter_On
                        [I_tmp,~] = XRT_2DAnlys_ImageFilter(I_tmp, ...
                            XRT_ImgPrc_Opts.ImgPrc_imFilter, ...
                            Image_Filter_Options);
                    end
                    
                   % Convert to binary image 
                    switch lower(XRT_2DImgPrc_Param.ImgPrc_bw_thresMethod)
                        case {'otsu'}
                            bw_tmp = I_tmp > XRT_2DImgPrc_Param.ImgPrc_bw_thresVal; %#ok<*PFBNS>
                        case {'activecontour'}
                            I_mask = imbinarize(I_tmp);
                            bw_tmp = activecontour(I_tmp,I_mask,XRT_2DImgPrc_Param.bw_thresVal_Method_Iter);
                        case {'ridlercalvard'}
                            bw_tmp = I_tmp > XRT_2DImgPrc_Param.ImgPrc_bw_thresVal;
                        case {'adaptive'}
                            bw_tmp = imfilter(I_tmp, XRT_2DImgPrc_Param.h);
                      otherwise
                            bw_tmp = I_tmp > XRT_2DImgPrc_Param.ImgPrc_bw_thresVal;
                    end
                    bw_tmp = XRT_2DAnlys_NoiseReduction(bw_tmp,XRT_ImgPrc_Opts.NoiseReduction2D);
                else
                    bw_tmp = I_tmp;
                end
            else
                % Assign zero array to bw_tmp to avoid runtime error
                bw_tmp = zero(size(I_tmp),'logical');
            end
            
            % Load ROI images in case Input_FolderPath_ROI was defined
            if XRT_3DAnlys_Opts.Input_FolderPath_ROI_Check
                bw_ROI = imread(fullfile(filelist_imgXRT_ROI_VolChunk(k).folder,filelist_imgXRT_ROI_VolChunk(k).name));
                if XRT_3DAnlys_Opts.DataResizeFactor > 1
                    % Resize image
                    bw_ROI = imresize(bw_ROI,[mImage nImage],'nearest');
                end  
            else
                bw_ROI = bw_tmp;
            end

            if ~V_binarize_Check % Add image to uint8 VolChunk
                ImageV_binary(:,:,k) = I_tmp;
            else % Add image to logical VolChunk
                ImageV_binary(:,:,k) = bw_tmp;
            end
            
            % Add ROI image to logical VolChunk
            ImageV_binary_ROI(:,:,k) = bw_ROI;
        end
        
        % Add VolChunk to volume / volume ROI
        V(:,:,VolChunk_Start:VolChunk_End) = ImageV_binary;
        V_ROI(:,:,VolChunk_Start:VolChunk_End) = ImageV_binary_ROI;
        
        fprintf('%s - %s - %s: VolChunk %.0f/%.0f \n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename(),iter,num_VolChunks)

    end
    fprintf('%s - %s - %s: >>> Complete <<< \n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename())
 
    clearvars ImageV_binary ImageV_binary_ROI
    
    %% NoiseReduction3D
    
    if isfield(XRT_3DAnlys_Opts,'NoiseReduction3D') && XRT_3DAnlys_Opts.NoiseReduction3D
        fileID = fopen(fullfile(XRT_3DAnlys_Opts.filesPath_2DImagePrc,sprintf('%s_ImgPrc_Parameter.txt',XRT_3DAnlys_Opts.ExpShorthand)),'at');
        fprintf(fileID,'\n[NoiseReduction3D]\n');
        fprintf(fileID,XRT_3DAnlys_Opts.formatSpec_s_d,'NoiseReduction3D_ClearBorder', XRT_ImgPrc_Opts.NoiseReduction3D.ClearBorder);
        fprintf(fileID,XRT_3DAnlys_Opts.formatSpec_s_d,'NoiseReduction3D_Close_Size_V', XRT_ImgPrc_Opts.NoiseReduction3D.Close_Size_V);
        fprintf(fileID,XRT_3DAnlys_Opts.formatSpec_s_d,'NoiseReduction3D_Open_Size_V', XRT_ImgPrc_Opts.NoiseReduction3D.Open_Size_V);
        fprintf(fileID,XRT_3DAnlys_Opts.formatSpec_s_d,'NoiseReduction3D_Close_Size_V_ROI', XRT_ImgPrc_Opts.NoiseReduction3D.Close_Size_V_ROI);
        fprintf(fileID,XRT_3DAnlys_Opts.formatSpec_s_d,'NoiseReduction3D_Open_Size_V_ROI', XRT_ImgPrc_Opts.NoiseReduction3D.Open_Size_V_ROI);
        
        if isfield(XRT_ImgPrc_Opts.NoiseReduction3D, 'method') 
            fprintf(fileID,XRT_3DAnlys_Opts.formatSpec_s_s,'NoiseReduction3D_method', XRT_ImgPrc_Opts.NoiseReduction3D.method);
        end
        fclose(fileID);
        
        if XRT_ImgPrc_Opts.NoiseReduction3D.ClearBorder || ... 
                XRT_ImgPrc_Opts.NoiseReduction3D.Close_Size_V > 0 || ...
                XRT_ImgPrc_Opts.NoiseReduction3D.Open_Size_V > 0
            
            XRT_3DAnlys_NoiseReduction_Opts.filesPath_FileExchange = XRT_Opts.filesPath_FileExchange;    
            XRT_3DAnlys_NoiseReduction_Opts.ClearBorder = XRT_ImgPrc_Opts.NoiseReduction3D.ClearBorder;
            if XRT_3DAnlys_NoiseReduction_Opts.ClearBorder
                ClearBorder_str = 'true';
                XRT_3DAnlys_NoiseReduction_Opts.Method = 'clearborder';
                XRT_3DAnlys_NoiseReduction_Opts.filesPath_FileExchange = XRT_Opts.filesPath_FileExchange;
                [V, ~] = XRT_3DAnlys_NoiseReduction(V,XRT_3DAnlys_NoiseReduction_Opts);
            else
                ClearBorder_str = 'false';
            end
            XRT_3DAnlys_NoiseReduction_Opts.Processing_mode = true;
            
            if XRT_ImgPrc_Opts.NoiseReduction3D.Close_Size_V > 0
                XRT_3DAnlys_NoiseReduction_Opts.Method = 'imclose';
                XRT_3DAnlys_NoiseReduction_Opts.Close_Size = XRT_ImgPrc_Opts.NoiseReduction3D.Close_Size_V;
            end
            
            if XRT_ImgPrc_Opts.NoiseReduction3D.Open_Size_V > 0
                XRT_3DAnlys_NoiseReduction_Opts.Method = 'imopen';
                XRT_3DAnlys_NoiseReduction_Opts.Open_Size = XRT_ImgPrc_Opts.NoiseReduction3D.Open_Size_V;
            end
            
            if islogical(V) && ~XRT_3DAnlys_Opts.Input_FolderPath_ROI_Check
                V_ROI = V;
            end
            XRT_3DAnlys_NoiseReduction_Opts.filesPath_FileExchange = XRT_Opts.filesPath_FileExchange;
            [V_ROI, NoiseReduction_Results] = XRT_3DAnlys_NoiseReduction(V_ROI,XRT_3DAnlys_NoiseReduction_Opts);
            if ~islogical(V)
                V = immultiply(V_ROI,V);
            else
                V  = V_ROI;
            end
            fprintf('%s - %s - %s: NoiseReduction3D >>> morph. Operation <<< 1/2 COMPLETE (ClearBorder: %s, CloseSize: %.0f, OpenSize: %.0f)\n', ...
                XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename(),ClearBorder_str, ...
                XRT_ImgPrc_Opts.NoiseReduction3D.Close_Size_V,XRT_ImgPrc_Opts.NoiseReduction3D.Open_Size_V)
        else
            fprintf('%s - %s - %s: NoiseReduction3D >>> morph. Operation <<< 1/2 PASSED \n', ...
                XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename())
        end
        
        if isfield(XRT_ImgPrc_Opts.NoiseReduction3D, 'method')
            XRT_3DAnlys_NoiseReduction_Opts.Method = XRT_ImgPrc_Opts.NoiseReduction3D.Method;
            if islogical(V) && ~XRT_3DAnlys_Opts.Input_FolderPath_ROI_Check
                V_ROI = V;
            end
            XRT_3DAnlys_NoiseReduction_Opts.filesPath_FileExchange = XRT_Opts.filesPath_FileExchange;
            [V_ROI, ~] = XRT_3DAnlys_NoiseReduction(V_ROI,XRT_3DAnlys_NoiseReduction_Opts);
            V = immultiply(V_ROI,V);
            fprintf('%s - %s - %s: NoiseReduction3D >>> %s <<< 2/2 COMPLETE \n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename(),XRT_ImgPrc_Opts.NoiseReduction3D.Method)                    
        else
            fprintf('%s - %s - %s: NoiseReduction3D >>> Main Method <<< 2/2 PASSED \n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename())                    
        end
        
        ImgName = sprintf('%s_4_I_init_bw_NoiseReduction3D.bmp',XRT_3DAnlys_Opts.ExpShorthand);
        imwrite(V(:,:,XRT_3DAnlys_Opts.idx_Img_c), fullfile(XRT_3DAnlys_Opts.filesPath_2DImagePrc,ImgName),'bmp');                
        [Image_scaled] = XRT_2DAnlys_AddScalebar(V(:,:,XRT_3DAnlys_Opts.idx_Img_c),XRT_3DAnlys_Opts.ImagePixelSize_rS,XRT_3DAnlys_Opts);
        imwrite(Image_scaled, fullfile(XRT_3DAnlys_Opts.filesPath_2DImagePrc,sprintf('%s_4_I_init_bw_NoiseReduction3D_scaled.bmp',XRT_3DAnlys_Opts.ExpShorthand)),'bmp')
    end
      
      
    
    %% Apply ROI Method
    if isfield(XRT_3DAnlys_Opts,'ROI_Method') && islogical(V)
        if ~isfield(XRT_3DAnlys_Opts,'ROI_Method_ShrinkWrap_strel') || isempty(XRT_3DAnlys_Opts.ROI_Method_ShrinkWrap_strel)
            if XRT_3DAnlys_Opts.ROI_Method_ShrinkWrap_Size > 30
                XRT_3DAnlys_Opts.ROI_Method_ShrinkWrap_strel = 'disk'; %'sphere';
            else
                XRT_3DAnlys_Opts.ROI_Method_ShrinkWrap_strel = 'sphere'; 
            end
            fprintf('%s - %s - %s: ROI_Method_ShrinkWrap_strel DEFAULT (%s)\n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename(),XRT_3DAnlys_Opts.ROI_Method_ShrinkWrap_strel)
        end
         switch lower(XRT_3DAnlys_Opts.ROI_Method)
            case {'convexhull'}
                fprintf('%s - %s - %s: Enter ROI Method >>> %s <<< \n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename(),XRT_3DAnlys_Opts.ROI_Method)
                XRT_3DAnlys_ConvexHull3_Opts.filesPath_FileExchange = XRT_Opts.filesPath_FileExchange;
                V_ROI = XRT_3DAnlys_ConvexHull3(V,XRT_3DAnlys_ConvexHull3_Opts);
                
                fprintf('%s - %s - %s: ROI >>> %s <<< COMPLETE\n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename(),XRT_3DAnlys_Opts.ROI_Method)
             case {'convexhull_2d'}
                [numRow,numCol,numImg] = size(V);
                fprintf('%s - %s - %s: Enter ROI Method >>> %s <<< \n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename(),XRT_3DAnlys_Opts.ROI_Method)
                 parfor i = 1:numImg
                     V_ROI(:,:,i) = bwconvhull(V(:,:,i));
                 end
                 fprintf('%s - %s - %s: ROI >>> %s <<< COMPLETE\n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename(),XRT_3DAnlys_Opts.ROI_Method)
             
             case {'imfill'}
                V_ROI = imfill(V,'holes');
                 fprintf('%s - %s - %s: ROI >>> %s <<< COMPLETE\n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename(),XRT_3DAnlys_Opts.ROI_Method)
           
             case {'imfill_2d'}
                [numRow,numCol,numImg] = size(V_ROI);
                fprintf('%s - %s - %s: Enter ROI Method >>> %s <<< \n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename(),XRT_3DAnlys_Opts.ROI_Method)
                 parfor i = 1:numImg
                     V_ROI(:,:,i) = imfill(V(:,:,i),'holes');
                 end
                 fprintf('%s - %s - %s: ROI >>> %s <<< COMPLETE\n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename(),XRT_3DAnlys_Opts.ROI_Method)
           
             case {'shrinkwrap_op'}
                V_ROI = imfill(V,'holes');
                ROI_Method_Options.ExpShorthand = XRT_3DAnlys_Opts.ExpShorthand;
                ROI_Method_Options.AppShorthand = sprintf('%s_ROI_ShrinkWrap_oP',mfilename());
                ROI_Method_Options.subFolderName_path = fullfile(XRT_3DAnlys_Opts.Output_FolderPath,XRT_3DAnlys_Opts.AppShorthand,'ROI_ShrinkWrap_oP');
                ROI_Method_Options.Method = 'shrinkwrap_op';
                ROI_Method_Options.imclose_strel_type = XRT_3DAnlys_Opts.ROI_Method_ShrinkWrap_strel;
                ROI_Method_Options.imclose_strel_size = XRT_3DAnlys_Opts.ROI_Method_ShrinkWrap_Size;
                ROI_Method_Options.ImagePixelSize = XRT_3DAnlys_Opts.ImagePixelSize;
                ROI_Method_Options.DataResizeFactor = XRT_3DAnlys_Opts.DataResizeFactor;
                ROI_Method_Options.Control_on = false;
                ROI_Method_Options.Control_Print_on = true;
                
                fprintf('%s - %s - %s: Enter ROI Method >>> %s <<< %s, %.f px \n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename(),XRT_3DAnlys_Opts.ROI_Method,ROI_Method_Options.imclose_strel_type,XRT_3DAnlys_Opts.ROI_Method_ShrinkWrap_Size)
                V_ROI = XRT_3DAnlys_Shrinkwrap(V_ROI,ROI_Method_Options);
                fprintf('%s - %s - %s: ROI >>> %s <<< COMPLETE\n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename(),XRT_3DAnlys_Opts.ROI_Method)
                
             case {'shrinkwrap'}
                ROI_Method_Options.ExpShorthand = XRT_3DAnlys_Opts.ExpShorthand;
                ROI_Method_Options.AppShorthand = sprintf('%s_ROI_ShrinkWrap',mfilename());
                ROI_Method_Options.subFolderName_path = fullfile(XRT_3DAnlys_Opts.Output_FolderPath,XRT_3DAnlys_Opts.AppShorthand,'ROI_ShrinkWrap');
                ROI_Method_Options.Method = 'shrinkwrap';
                ROI_Method_Options.imclose_strel_type = XRT_3DAnlys_Opts.ROI_Method_ShrinkWrap_strel;
                ROI_Method_Options.imclose_strel_size = XRT_3DAnlys_Opts.ROI_Method_ShrinkWrap_Size;
                ROI_Method_Options.ImagePixelSize = XRT_3DAnlys_Opts.ImagePixelSize;
                ROI_Method_Options.DataResizeFactor = XRT_3DAnlys_Opts.DataResizeFactor;
                ROI_Method_Options.Control_on = false;
                ROI_Method_Options.Control_Print_on = true;
                
                fprintf('%s - %s - %s: Enter ROI Method >>> %s <<< %s, %.f px \n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename(),XRT_3DAnlys_Opts.ROI_Method,ROI_Method_Options.imclose_strel_type,XRT_3DAnlys_Opts.ROI_Method_ShrinkWrap_Size)
                V_ROI = XRT_3DAnlys_Shrinkwrap(V,ROI_Method_Options);
                fprintf('%s - %s - %s: ROI >>> %s <<< COMPLETE\n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename(),XRT_3DAnlys_Opts.ROI_Method)
            otherwise
                fprintf('%s - %s - %s: ROI Method PASSED\n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename())
            
         end
         

    elseif  isfield(XRT_3DAnlys_Opts,'ROI_Method') && ~islogical(V)
        % ROI_Method_ShrinkWrap_strel
        if ~isfield(XRT_3DAnlys_Opts,'ROI_Method_ShrinkWrap_strel') && isempty(XRT_3DAnlys_Opts.ROI_Method_ShrinkWrap_strel)
            if XRT_3DAnlys_Opts.ROI_Method_ShrinkWrap_Size > 30
                XRT_3DAnlys_Opts.ROI_Method_ShrinkWrap_strel = 'disk'; %'sphere';
            else
                XRT_3DAnlys_Opts.ROI_Method_ShrinkWrap_strel = 'sphere'; 
            end
            fprintf('%s - %s - %s: ROI_Method_ShrinkWrap_strel DEFAULT (%s)\n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename(),XRT_3DAnlys_Opts.ROI_Method_ShrinkWrap_strel)
        end

        switch lower(XRT_3DAnlys_Opts.ROI_Method)
            case {'convexhull'}
                XRT_3DAnlys_ConvexHull3_Opts.filesPath_FileExchange = XRT_Opts.filesPath_FileExchange;
                V_ROI = XRT_3DAnlys_ConvexHull3(V_ROI,XRT_3DAnlys_ConvexHull3_Opts);
                
                fprintf('%s - %s - %s: ROI >>> %s <<< COMPLETE\n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename(),XRT_3DAnlys_Opts.ROI_Method)
            case {'convexhull_2d'}
                [numRow,numCol,numImg] = size(V_ROI);
                fprintf('%s - %s - %s: Enter ROI Method >>> %s <<< \n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename(),XRT_3DAnlys_Opts.ROI_Method)
                 parfor i = 1:numImg
                     V_ROI(:,:,i) = bwconvhull(V_ROI(:,:,i));
                 end
                 fprintf('%s - %s - %s: ROI >>> %s <<< COMPLETE\n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename(),XRT_3DAnlys_Opts.ROI_Method)
             case {'imfill'}
                V_ROI = imfill(V_ROI,'holes');
                 fprintf('%s - %s - %s: ROI >>> %s <<< COMPLETE\n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename(),XRT_3DAnlys_Opts.ROI_Method)
            case {'imfill_2d'}
                [numRow,numCol,numImg] = size(V_ROI);
                fprintf('%s - %s - %s: Enter ROI Method >>> %s <<< \n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename(),XRT_3DAnlys_Opts.ROI_Method)
                 parfor i = 1:numImg
                     V_ROI(:,:,i) = imfill(V_ROI(:,:,i),'holes');
                 end
                 fprintf('%s - %s - %s: ROI >>> %s <<< COMPLETE\n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename(),XRT_3DAnlys_Opts.ROI_Method)
           
            case {'shrinkwrap_op'}
                V_ROI = imfill(V_ROI,'holes');
                ROI_Method_Options.ExpShorthand = XRT_3DAnlys_Opts.ExpShorthand;
                ROI_Method_Options.AppShorthand = sprintf('%s_ROI_ShrinkWrap_oP',mfilename());
                ROI_Method_Options.subFolderName_path = fullfile(XRT_3DAnlys_Opts.Output_FolderPath,XRT_3DAnlys_Opts.AppShorthand,'ROI_ShrinkWrap_oP');
                ROI_Method_Options.Method = 'shrinkwrap_op';
                ROI_Method_Options.imclose_strel_type = XRT_3DAnlys_Opts.ROI_Method_ShrinkWrap_strel;
                ROI_Method_Options.imclose_strel_size = XRT_3DAnlys_Opts.ROI_Method_ShrinkWrap_Size;
                ROI_Method_Options.ImagePixelSize = XRT_3DAnlys_Opts.ImagePixelSize;
                ROI_Method_Options.DataResizeFactor = XRT_3DAnlys_Opts.DataResizeFactor;
                ROI_Method_Options.Control_on = false;
                ROI_Method_Options.Control_Print_on = true;
                
                fprintf('%s - %s - %s: Enter ROI Method >>> %s <<< %s, %.f px \n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename(),XRT_3DAnlys_Opts.ROI_Method,ROI_Method_Options.imclose_strel_type,XRT_3DAnlys_Opts.ROI_Method_ShrinkWrap_Size)
                V_ROI = XRT_3DAnlys_Shrinkwrap(V_ROI,ROI_Method_Options);
                fprintf('%s - %s - %s: ROI >>> %s <<< COMPLETE\n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,XRT_3DAnlys_Opts.ROI_Method)
                
             case {'shrinkwrap'}
                V_ROI = imfill(V_ROI,'holes');
                ROI_Method_Options.ExpShorthand = XRT_3DAnlys_Opts.ExpShorthand;
                ROI_Method_Options.AppShorthand = sprintf('%s_ROI_ShrinkWrap',mfilename());
                ROI_Method_Options.subFolderName_path = fullfile(XRT_3DAnlys_Opts.Output_FolderPath,XRT_3DAnlys_Opts.AppShorthand,'ROI_ShrinkWrap');
                ROI_Method_Options.Method = 'shrinkwrap';
                ROI_Method_Options.imclose_strel_type = XRT_3DAnlys_Opts.ROI_Method_ShrinkWrap_strel;
                ROI_Method_Options.imclose_strel_size = XRT_3DAnlys_Opts.ROI_Method_ShrinkWrap_Size;
                ROI_Method_Options.ImagePixelSize = XRT_3DAnlys_Opts.ImagePixelSize;
                ROI_Method_Options.DataResizeFactor = XRT_3DAnlys_Opts.DataResizeFactor;
                ROI_Method_Options.Control_on = false;
                ROI_Method_Options.Control_Print_on = true;
                
                fprintf('%s - %s - %s: Enter ROI Method >>> %s <<< %s, %.f px \n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename(),XRT_3DAnlys_Opts.ROI_Method,ROI_Method_Options.imclose_strel_type,XRT_3DAnlys_Opts.ROI_Method_ShrinkWrap_Size)
                V_ROI = XRT_3DAnlys_Shrinkwrap(V_ROI,ROI_Method_Options);
                fprintf('%s - %s - %s: ROI >>> %s <<< COMPLETE\n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename(),XRT_3DAnlys_Opts.ROI_Method)
             case {'shrinkwrap_2d'}
                V_ROI = imfill(V_ROI,'holes');
                ROI_Method_Options.ExpShorthand = XRT_3DAnlys_Opts.ExpShorthand;
                ROI_Method_Options.AppShorthand = sprintf('%s_ROI_ShrinkWrap',mfilename());
                ROI_Method_Options.subFolderName_path = fullfile(XRT_3DAnlys_Opts.Output_FolderPath,XRT_3DAnlys_Opts.AppShorthand,'ROI_ShrinkWrap');
                ROI_Method_Options.Method = 'shrinkwrap';
                ROI_Method_Options.imclose_strel_type = XRT_3DAnlys_Opts.ROI_Method_ShrinkWrap_strel;
                ROI_Method_Options.imclose_strel_size = XRT_3DAnlys_Opts.ROI_Method_ShrinkWrap_Size;
                ROI_Method_Options.ImagePixelSize = XRT_3DAnlys_Opts.ImagePixelSize;
                ROI_Method_Options.DataResizeFactor = XRT_3DAnlys_Opts.DataResizeFactor;
                ROI_Method_Options.Control_on = false;
                ROI_Method_Options.Control_Print_on = false;
                
                fprintf('%s - %s - %s: Enter ROI Method >>> %s <<< %s, %.f px \n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename(),XRT_3DAnlys_Opts.ROI_Method,ROI_Method_Options.imclose_strel_type,XRT_3DAnlys_Opts.ROI_Method_ShrinkWrap_Size)
                [numRow,numCol,numImg] = size(V_ROI);
                parfor i = 1:numImg
                    V_ROI(:,:,i) = XRT_3DAnlys_Shrinkwrap(V_ROI(:,:,i),ROI_Method_Options);
                 end
                fprintf('%s - %s - %s: ROI >>> %s <<< COMPLETE\n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename(),XRT_3DAnlys_Opts.ROI_Method)
            otherwise
                fprintf('%s - %s - %s: ROI Method PASSED\n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename())
            
        end
         
        V = immultiply(V,V_ROI);
    elseif XRT_3DAnlys_Opts.Input_FolderPath_ROI_Check && ~islogical(V)
        fprintf('%s - %s - %s: Enter ROI Update - Apply on V \n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename())
        V = immultiply(V,V_ROI);
        fprintf('%s - %s - %s: ROI Update - Apply on V \n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename())
    end
    
        ImgName = sprintf('%s_5_I_init_ROI.bmp',XRT_3DAnlys_Opts.ExpShorthand);
        imwrite(V(:,:,XRT_3DAnlys_Opts.idx_Img_c), fullfile(XRT_3DAnlys_Opts.filesPath_2DImagePrc,ImgName),'bmp');                
        [Image_scaled] = XRT_2DAnlys_AddScalebar(V(:,:,XRT_3DAnlys_Opts.idx_Img_c),XRT_3DAnlys_Opts.ImagePixelSize_rS,XRT_3DAnlys_Opts);
        imwrite(Image_scaled, fullfile(XRT_3DAnlys_Opts.filesPath_2DImagePrc,sprintf('%s_5_I_init_ROI_scaled.bmp',XRT_3DAnlys_Opts.ExpShorthand)),'bmp')
        
        ImgName = sprintf('%s_5_ROI.bmp',XRT_3DAnlys_Opts.ExpShorthand);
        imwrite(V_ROI(:,:,XRT_3DAnlys_Opts.idx_Img_c), fullfile(XRT_3DAnlys_Opts.filesPath_2DImagePrc,ImgName),'bmp');                
        [Image_scaled] = XRT_2DAnlys_AddScalebar(V_ROI(:,:,XRT_3DAnlys_Opts.idx_Img_c),XRT_3DAnlys_Opts.ImagePixelSize_rS,XRT_3DAnlys_Opts);
        imwrite(Image_scaled, fullfile(XRT_3DAnlys_Opts.filesPath_2DImagePrc,sprintf('%s_5_ROI_scaled.bmp',XRT_3DAnlys_Opts.ExpShorthand)),'bmp')
        
%% NoiseReduction3D V_ROI
    if isfield(XRT_3DAnlys_Opts,'NoiseReduction3D') && XRT_3DAnlys_Opts.NoiseReduction3D        
        
        if XRT_ImgPrc_Opts.NoiseReduction3D.Close_Size_V_ROI > 0 || ...
                XRT_ImgPrc_Opts.NoiseReduction3D.Open_Size_V_ROI > 0
                
            XRT_3DAnlys_NoiseReduction_Opts.ClearBorder = false;
            if XRT_3DAnlys_NoiseReduction_Opts.ClearBorder
                ClearBorder_str = 'true';
                XRT_3DAnlys_NoiseReduction_Opts.Method = 'clearborder';
                XRT_3DAnlys_NoiseReduction_Opts.filesPath_FileExchange = XRT_Opts.filesPath_FileExchange;
                [V_ROI, ~] = XRT_3DAnlys_NoiseReduction(V_ROI,XRT_3DAnlys_NoiseReduction_Opts);
            else
                ClearBorder_str = 'false';
            end
            XRT_3DAnlys_NoiseReduction_Opts.Processing_mode = true;
            if XRT_ImgPrc_Opts.NoiseReduction3D.Close_Size_V_ROI > 0
                XRT_3DAnlys_NoiseReduction_Opts.Method = 'imclose';
                XRT_3DAnlys_NoiseReduction_Opts.Close_Size = XRT_ImgPrc_Opts.NoiseReduction3D.Close_Size_V_ROI;
                XRT_3DAnlys_NoiseReduction_Opts.filesPath_FileExchange = XRT_Opts.filesPath_FileExchange;
                [V_ROI, ~] = XRT_3DAnlys_NoiseReduction(V_ROI,XRT_3DAnlys_NoiseReduction_Opts);
            end
            
            if XRT_ImgPrc_Opts.NoiseReduction3D.Open_Size_V_ROI > 0
                XRT_3DAnlys_NoiseReduction_Opts.Method = 'imopen';
                XRT_3DAnlys_NoiseReduction_Opts.Open_Size = XRT_ImgPrc_Opts.NoiseReduction3D.Open_Size_V_ROI;
                XRT_3DAnlys_NoiseReduction_Opts.filesPath_FileExchange = XRT_Opts.filesPath_FileExchange;
                [V_ROI, ~] = XRT_3DAnlys_NoiseReduction(V_ROI,XRT_3DAnlys_NoiseReduction_Opts);
            end
            
            V = immultiply(V_ROI,V);
            fprintf('%s - %s - %s: NoiseReduction3D_ROI >>> morph. Operation <<< COMPLETE (CloseSize: %.0f, OpenSize: %.0f)\n', ...
                XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename(),XRT_ImgPrc_Opts.NoiseReduction3D.Close_Size_V_ROI,XRT_ImgPrc_Opts.NoiseReduction3D.Open_Size_V_ROI)
        else
            fprintf('%s - %s - %s: NoiseReduction3D_ROI >>> morph. Operation <<< PASSED \n', ...
                XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename())
        end
        
        ImgName = sprintf('%s_4_I_init_bw_NoiseReduction3D_ROI.bmp',XRT_3DAnlys_Opts.ExpShorthand);
        imwrite(V(:,:,XRT_3DAnlys_Opts.idx_Img_c), fullfile(XRT_3DAnlys_Opts.filesPath_2DImagePrc,ImgName),'bmp');                
        [Image_scaled] = XRT_2DAnlys_AddScalebar(V(:,:,XRT_3DAnlys_Opts.idx_Img_c),XRT_3DAnlys_Opts.ImagePixelSize_rS,XRT_3DAnlys_Opts);
        imwrite(Image_scaled, fullfile(XRT_3DAnlys_Opts.filesPath_2DImagePrc,sprintf('%s_4_I_init_bw_NoiseReduction3D_ROI_scaled.bmp',XRT_3DAnlys_Opts.ExpShorthand)),'bmp')
    end
    
    
    %% Reduce Image Stack size from V_ROI
    if XRT_3DAnlys_Opts.ROI_Reduction_Check
        fprintf('%s - %s - %s: VOI Reduction \n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename())
        
        clearvars XRT_3DAnlys_ROI_Reduction_Opts
        XRT_3DAnlys_ROI_Reduction_Opts.numSpacer = 1;
        [V_ROI,XRT_3DAnlys_ROI_Reduction_Results] = XRT_3DAnlys_ROI_Reduction(V_ROI,XRT_3DAnlys_ROI_Reduction_Opts);
        
        XRT_3DAnlys_ROI_Reduction_Opts.numImg_Start = XRT_3DAnlys_ROI_Reduction_Results.numImg_Start;
        XRT_3DAnlys_ROI_Reduction_Opts.num_Img_End = XRT_3DAnlys_ROI_Reduction_Results.num_Img_End;
        XRT_3DAnlys_ROI_Reduction_Opts.row_Min = XRT_3DAnlys_ROI_Reduction_Results.row_Min;
        XRT_3DAnlys_ROI_Reduction_Opts.row_Max = XRT_3DAnlys_ROI_Reduction_Results.row_Max;
        XRT_3DAnlys_ROI_Reduction_Opts.column_Min = XRT_3DAnlys_ROI_Reduction_Results.column_Min;
        XRT_3DAnlys_ROI_Reduction_Opts.column_Max = XRT_3DAnlys_ROI_Reduction_Results.column_Max;
        XRT_3DAnlys_ROI_Reduction_Opts.numSpacer = XRT_3DAnlys_ROI_Reduction_Results.numSpacer;
        
        if XRT_3DAnlys_ROI_Reduction_Results.V_any
            [V,XRT_3DAnlys_ROI_Reduction_Results] = XRT_3DAnlys_ROI_Reduction(V,XRT_3DAnlys_ROI_Reduction_Opts);
        end
        
        fileID = fopen(fullfile(XRT_3DAnlys_Opts.filesPath_2DImagePrc,sprintf('%s_ImgPrc_Parameter.txt',XRT_3DAnlys_Opts.ExpShorthand)),'at');
        fprintf(fileID,'\n[VOI_Reduction]\r\n');
        fprintf(fileID,XRT_3DAnlys_Opts.formatSpec_s_s,'VOI_Reduction', 'TRUE');
        fprintf(fileID,XRT_3DAnlys_Opts.formatSpec_s_d,'VOI_Reduction_numImg_Start', XRT_3DAnlys_ROI_Reduction_Results.numImg_Start);
        fprintf(fileID,XRT_3DAnlys_Opts.formatSpec_s_d,'VOI_Reduction_numImg_Start', XRT_3DAnlys_ROI_Reduction_Results.num_Img_End);
        fprintf(fileID,XRT_3DAnlys_Opts.formatSpec_s_d,'VOI_Reduction_row_Min', XRT_3DAnlys_ROI_Reduction_Results.row_Min);
        fprintf(fileID,XRT_3DAnlys_Opts.formatSpec_s_d,'VOI_Reduction_row_Max', XRT_3DAnlys_ROI_Reduction_Results.row_Max);
        fprintf(fileID,XRT_3DAnlys_Opts.formatSpec_s_d,'VOI_Reduction_column_Min', XRT_3DAnlys_ROI_Reduction_Results.column_Min);
        fprintf(fileID,XRT_3DAnlys_Opts.formatSpec_s_d,'VOI_Reduction_column_Max', XRT_3DAnlys_ROI_Reduction_Results.column_Max);
        fprintf(fileID,XRT_3DAnlys_Opts.formatSpec_s_d,'VOI_Reduction_numSpacer', XRT_3DAnlys_ROI_Reduction_Opts.numSpacer);
        fclose(fileID);
        
        fprintf('%s - %s - %s: VOI Reduction COMPLETE \n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename())
    end
    


    %% XRT_3DAnlys_VolExport
    if ~isfield(XRT_3DAnlys_Opts,'VolExport_SurfRend_Check')
        XRT_3DAnlys_Opts.VolExport_SurfRend_Check = false;
    end
    
    if XRT_3DAnlys_Opts.VolExport_SurfRend_Check && V_binarize_Check
        fprintf('%s - %s - %s: Export Volume (Surface Rendered) \n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename())
        XRT_3DAnlys_VolExport_Options.ExpShorthand = XRT_3DAnlys_Opts.ExpShorthand;
        XRT_3DAnlys_VolExport_Options.AppShorthand = XRT_3DAnlys_Opts.AppShorthand;
        XRT_3DAnlys_VolExport_Options.path_FileExchange = XRT_Opts.filesPath_FileExchange;
        XRT_3DAnlys_VolExport_Options.path_ImgFolder = XRT_3DAnlys_Opts.Output_FolderPath;
        XRT_3DAnlys_VolExport_Options.ImgName = sprintf('%s_VolExport',XRT_3DAnlys_Opts.ExpShorthand);
        XRT_3DAnlys_VolExport_Options.AxisOff_On = false;
        XRT_3DAnlys_VolExport_Options.ParforComp = false;
        
        XRT_3DAnlys_VolExport(V,XRT_3DAnlys_VolExport_Options)
        
        XRT_3DAnlys_VolExport(V_ROI,XRT_3DAnlys_VolExport_Options)
    elseif XRT_3DAnlys_Opts.VolExport_SurfRend_Check
        fprintf('%s - %s: Export Volume (Surface Rendered) \n',XRT_3DAnlys_Opts.ExpShorthand,mfilename())
        XRT_3DAnlys_VolExport_Options.ExpShorthand = XRT_3DAnlys_Opts.ExpShorthand;
        XRT_3DAnlys_VolExport_Options.AppShorthand = XRT_3DAnlys_Opts.AppShorthand;
        XRT_3DAnlys_VolExport_Options.path_FileExchange = XRT_Opts.filesPath_FileExchange;
        XRT_3DAnlys_VolExport_Options.subFolderName = fullfile('XRT_3DAnlys_VolExport');
        XRT_3DAnlys_VolExport_Options.AxisOff_On = false;
        XRT_3DAnlys_VolExport_Options.ParforComp = false;
        XRT_3DAnlys_VolExport(V_ROI,XRT_3DAnlys_VolExport_Options)
    end

    %% Save V, V_ROI to ImgStack
    fprintf('%s - %s - %s: Save Image Stack \n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename())
    if V_binarize_Check
        ImgStack_Create_Opts.path_ImgFolder = fullfile(XRT_3DAnlys_Opts.Output_FolderPath,'XRT_3DAnlys_ImgStack_V_bin');
        ImgStack_Create_Opts.ImgName = 'V';
        ImgStack_Create_Opts.ImgFormat = 'bmp';
        ImgStack_Create_Opts.pool_mode = true;
        ImgStack_Create_Opts.Print_Log_On = false;
        ImgStack_Create_Opts.path_FileExchange = XRT_Opts.filesPath_FileExchange;
        ImgStack_Create_Opts.ExpShorthand = XRT_3DAnlys_Opts.ExpShorthand;
        ImgStack_Create_Opts.AppShorthand = XRT_3DAnlys_Opts.AppShorthand;
        ImgStack_Create_Opts.Save_IndImages_On = true;
        ImgStack_Create_Opts.ImageFormat = 'logical';
        XRT_3DAnlys_ImgStack_Create(V,ImgStack_Create_Opts);
    else
        ImgStack_Create_Opts.path_ImgFolder = fullfile(XRT_3DAnlys_Opts.Output_FolderPath,'XRT_3DAnlys_ImgStack_V_uint8');
        ImgStack_Create_Opts.ImgName = 'V';
        ImgStack_Create_Opts.ImgFormat = 'bmp';
        ImgStack_Create_Opts.pool_mode = true;
        ImgStack_Create_Opts.Print_Log_On = false;
        ImgStack_Create_Opts.path_FileExchange = XRT_Opts.filesPath_FileExchange;
        ImgStack_Create_Opts.ExpShorthand = XRT_3DAnlys_Opts.ExpShorthand;
        ImgStack_Create_Opts.AppShorthand = XRT_3DAnlys_Opts.AppShorthand;
        ImgStack_Create_Opts.Save_IndImages_On = true;
        ImgStack_Create_Opts.ImageFormat = 'unit8';
        XRT_3DAnlys_ImgStack_Create(V,ImgStack_Create_Opts);
    end
   
    ImgStack_Create_Opts.path_ImgFolder = fullfile(XRT_3DAnlys_Opts.Output_FolderPath,'XRT_3DAnlys_ImgStack_V_ROI');
    ImgStack_Create_Opts.ImgName = 'V_ROI';
    ImgStack_Create_Opts.ImgFormat = 'bmp';
    ImgStack_Create_Opts.pool_mode = true;
    ImgStack_Create_Opts.Print_Log_On = false;
    ImgStack_Create_Opts.path_FileExchange = XRT_Opts.filesPath_FileExchange;
    ImgStack_Create_Opts.ExpShorthand = XRT_3DAnlys_Opts.ExpShorthand;
    ImgStack_Create_Opts.AppShorthand = XRT_3DAnlys_Opts.AppShorthand;
    ImgStack_Create_Opts.Save_IndImages_On = true;
    ImgStack_Create_Opts.ImageFormat = 'logical';
    XRT_3DAnlys_ImgStack_Create(V_ROI,ImgStack_Create_Opts);
    
    XRT_3DAnlys_Opts.VolumeBuild_on = false;
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %   

elseif VolumeLoad_on
        ImgStack_Load_Options.ExpShorthand = XRT_3DAnlys_Opts.ExpShorthand;
        ImgStack_Load_Options.AppShorthand = XRT_3DAnlys_Opts.AppShorthand;
        ImgStack_Load_Options.pool_mode = true;
        ImgStack_Load_Options.ImgFormat = XRT_3DAnlys_Opts.ImgFormat;
        ImgStack_Create_Opts.Print_Log_On = true;
        
        if V_binarize_Check
            ImgStack_Load_Options.path_ImgFolder = fullfile(XRT_3DAnlys_Opts.Output_FolderPath,'XRT_3DAnlys_ImgStack_V_logical');
            V = XRT_3DAnlys_ImgStack_Load(ImgStack_Load_Options);
        else
            ImgStack_Load_Options.path_ImgFolder = fullfile(XRT_3DAnlys_Opts.Output_FolderPath,'XRT_3DAnlys_ImgStack_V_uint8'); 
            V = XRT_3DAnlys_ImgStack_Load(ImgStack_Load_Options);
        end
        ImgStack_Load_Options.path_ImgFolder = fullfile(XRT_3DAnlys_Opts.Output_FolderPath,'XRT_3DAnlys_ImgStack_V_ROI'); 
        V_ROI = XRT_3DAnlys_ImgStack_Load(ImgStack_Load_Options);
end


clearvars V_binarize_Check

