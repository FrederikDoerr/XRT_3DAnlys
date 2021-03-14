%------------------------------------------------------------------------------------------------
% Code written by Frederik Doerr, Aug 2020 (MATLAB R2019b)
% Application: XRT Image and Data Analysis Framework (SubRoutine
% XRT_3DAnlys_VolumeAnlys)
% 
% https://github.com/frederik-d
% Contact: frederik.doerr(at)strath.ac.uk | CMAC (http://www.cmac.ac.uk/)

clearvars -except XRT_3DAnlys_Opts XRT_3DAnlys_MasterDB Batch3DAnlys_Iter numExp pool XRT_Opts
XRT_ImgPrc_Opts = struct();
[XRT_3DAnlys_Opts,XRT_ImgPrc_Opts] = XRT_3DAnlys_ParameterLoader(XRT_Opts,Batch3DAnlys_Iter,XRT_3DAnlys_Opts,XRT_ImgPrc_Opts);
XRT_3DAnlys_Opts.AppShorthand = mfilename();

if XRT_3DAnlys_Opts.VolumeAnlys_On
    % log-file
    fprintf('- - - - - - - - - - - -\n')
    fprintf('%s - Enter %s (Elapsed time: %.0f sec)\n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,toc(XRT_3DAnlys_Opts.tic))
    fprintf('- - - - - - - - - - - -\n')
    pause(2)
    
    XRT_3DAnlys_Opts.PrcTime_Logger{end+1,1} = mfilename();
    XRT_3DAnlys_Opts.PrcTime_Logger{end,2} = 'Start';
    XRT_3DAnlys_Opts.PrcTime_Logger{end,3} = toc(XRT_3DAnlys_Opts.tic);
    
    % % % % % % % % % % % 0) XRT_3DAnlys_ReconImage_VolumeBuild % % % % % % % % % % %
    V_binarize_Check = true;
    run XRT_3DAnlys_VolumeBuild

    % Export 2D Images
    XRT_3DAnlys_ImgExport_Opts.ExpShorthand = XRT_3DAnlys_Opts.ExpShorthand;
    XRT_3DAnlys_ImgExport_Opts.path_ImgFolder = fullfile(XRT_3DAnlys_Opts.Output_FolderPath,XRT_3DAnlys_Opts.AppShorthand);
    XRT_3DAnlys_ImgExport_Opts.ImagePixelSize_rS = XRT_3DAnlys_Opts.ImagePixelSize_rS;
    XRT_3DAnlys_ImgExport_Opts.ImgName = sprintf('%s_V.bmp',XRT_3DAnlys_Opts.ExpShorthand);
    XRT_3DAnlys_ImgExport(V,XRT_3DAnlys_ImgExport_Opts)
    
    XRT_3DAnlys_ImgExport_Opts.ImgName = sprintf('%s_V_ROI.bmp',XRT_3DAnlys_Opts.ExpShorthand);
    XRT_3DAnlys_ImgExport(V_ROI,XRT_3DAnlys_ImgExport_Opts)

    % % % % % % % % % % % 
    % log-file
    fprintf('%s - %s - XRT_3DAnlys_VolumeBuild (Elapsed time: %.0f sec)\n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,toc(XRT_3DAnlys_Opts.tic))
    XRT_3DAnlys_Opts.PrcTime_Logger{end+1,1} = mfilename();
    XRT_3DAnlys_Opts.PrcTime_Logger{end,2} = 'XRT_3DAnlys_VolumeBuild';
    XRT_3DAnlys_Opts.PrcTime_Logger{end,3} = toc(XRT_3DAnlys_Opts.tic);

    %% % % % % % % % % % % 1) XRT_3DAnlys_MophDescriptors: Extract Global Size/Shape Descriptors % % % % % % % % % % %

    if isfield(XRT_3DAnlys_Opts,'VolumeAnlys_3DDescriptors') && ~isnan(XRT_3DAnlys_Opts.VolumeAnlys_3DDescriptors) && XRT_3DAnlys_Opts.VolumeAnlys_3DDescriptors
        fprintf('%s - %s - XRT_3DAnlys_Descriptors INITIATE (Elapsed time: %.0f sec)\n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,toc(XRT_3DAnlys_Opts.tic))
        
        XRT_3DAnlys_Descriptors_Opts.ExpShorthand = XRT_3DAnlys_Opts.ExpShorthand;
        XRT_3DAnlys_Descriptors_Opts.AppShorthand = XRT_3DAnlys_Opts.AppShorthand;
        XRT_3DAnlys_Descriptors_Opts.subFolderName = XRT_3DAnlys_Opts.AppShorthand;
        XRT_3DAnlys_Descriptors_Opts.ImagePixelSize = XRT_3DAnlys_Opts.ImagePixelSize;
        XRT_3DAnlys_Descriptors_Opts.DataResizeFactor = XRT_3DAnlys_Opts.DataResizeFactor;
        XRT_3DAnlys_Descriptors_Opts.filesPath_FileExchange = XRT_Opts.filesPath_FileExchange;

        XRT_3DAnlys_Descriptors_Opts.numImg_min = 0;
        XRT_3DAnlys_Descriptors_Opts.numRow_min = 0;
        XRT_3DAnlys_Descriptors_Opts.numCol_min = 0;

        [V_stats,V_stats3,V_MophDescriptors_Results] = XRT_3DAnlys_Descriptors(V,XRT_3DAnlys_Descriptors_Opts);
        [V_ROI_stats,V_ROI_stats3,V_ROI_MophDescriptors_Results] = XRT_3DAnlys_Descriptors(V_ROI,XRT_3DAnlys_Descriptors_Opts);

        save(sprintf('%s_%s',XRT_3DAnlys_Opts.AppShorthand,'XRT_3DAnlys_Descriptors'),'V_stats','V_stats3','V_ROI_stats','V_ROI_stats3')

        [numRow,numCol,numImg] = size(V);
        size_ImageBox_X = numCol*XRT_3DAnlys_Opts.ImagePixelSize_rS;
        size_ImageBox_Y = numRow*XRT_3DAnlys_Opts.ImagePixelSize_rS;
        size_ImageBox_Z = numImg*XRT_3DAnlys_Opts.ImagePixelSize_rS;
        V_stats3.Solidity_3_Max = V_stats3.Volume_3_Max/V_ROI_stats3.Volume_3_Max;
        
        fileID = fopen(fullfile(XRT_3DAnlys_Opts.AppShorthand,sprintf('%s_VolumeAnlys_stats.csv',XRT_3DAnlys_Opts.ExpShorthand)),'wt');
        Data_stats = { ...
            'size_ImageBox_X',size_ImageBox_X,'um'; ...
            'size_ImageBox_Y',size_ImageBox_Y,'um'; ...
            'size_ImageBox_Z',size_ImageBox_Z,'um'; ...
            'V_BoundingBox_MaxLength',V_stats.BoundingBox(1),'um'; ...
            'V_BoundingBox_MinLength',V_stats.BoundingBox(2),'um'; ...
            'V_ROI_BoundingBox_MaxLength',V_ROI_stats.BoundingBox(1),'um'; ...
            'V_ROI_BoundingBox_MinLength',V_ROI_stats.BoundingBox(2),'um'; ...
            'NumEl_Volume',numel(V_stats.Area),''; ...
            'V_Volume_um',sum(V_stats.Area),'um^3'; ...
            'V_ROI_Volume_um',sum(V_ROI_stats.Area),'um^3'; ...
            'V_Surface_um',V_stats.Surface,'um^2'; ...
            'V_ROI_Surface_um',V_ROI_stats.Surface,'um^2'; ...
            'V_radius_spherefit_um',V_stats.radius_spherefit,'um'; ...
            'V_ROI_radius_spherefit_um',V_ROI_stats.radius_spherefit,'um'; ...
            'V_radii_elpsfit_1_um',V_stats.radii_elpsfit(1),'um'; ...
            'V_radii_elpsfit_2_um',V_stats.radii_elpsfit(2),'um'; ...
            'V_radii_elpsfit_3_um',V_stats.radii_elpsfit(3),'um'; ...
            'V_chi2_elpsfit_um',V_stats.chi2_elpsfit,'um'; ...
            'V_ROI_radii_elpsfit_1_um',V_ROI_stats.radii_elpsfit(1),'um'; ...
            'V_ROI_radii_elpsfit_2_um',V_ROI_stats.radii_elpsfit(2),'um'; ...
            'V_ROI_radii_elpsfit_3_um',V_ROI_stats.radii_elpsfit(3),'um'; ...
            'V_ROI_chi2_elpsfit_um',V_ROI_stats.chi2_elpsfit,'um'; ...
            'V_Solidity',sum(V_stats.Area)/sum(V_ROI_stats.Area),'';
            'V_AspectRatio',V_stats.AspectRatio,''; ...
            'V_ROI_AspectRatio',V_ROI_stats.AspectRatio,''; ...
            };
        for i = 1:length(Data_stats)
            fprintf(fileID,XRT_3DAnlys_Opts.formatSpec_s_d_s,Data_stats{i,:});
        end  
        fclose(fileID);


        fileID = fopen(fullfile(XRT_3DAnlys_Opts.AppShorthand,sprintf('%s_VolumeAnlys_stats3.csv',XRT_3DAnlys_Opts.ExpShorthand)),'wt');
        Data_stats3 = { ...
            'size_ImageBox_X',numCol*XRT_3DAnlys_Opts.ImagePixelSize_rS,'um'; ...
            'size_ImageBox_Y',numRow*XRT_3DAnlys_Opts.ImagePixelSize_rS,'um'; ...
            'size_ImageBox_Z',numImg*XRT_3DAnlys_Opts.ImagePixelSize_rS,'um'; ...
            'V_Volume_um',sum(V_stats3.Volume_3_Max),'um^3'; ...
            'V_ROI_Volume_um',sum(V_ROI_stats3.Volume_3_Max),'um^3'; ...
            'V_ConvexVolume',V_stats3.ConvexVolume_3_Max,'um^3';
            'V_ROI_ConvexVolume',V_ROI_stats3.ConvexVolume_3_Max,'um^3';
            'V_Surface_um',V_stats3.SurfaceArea_3_Max,'um^2'; ...
            'V_ROI_Surface_um',V_stats3.SurfaceArea_3_Max,'um^2'; ...
            'V_PrincipalAxisLength_1_um',V_stats3.PrincipalAxisLength_3_Max(1),'um'; ...
            'V_PrincipalAxisLength_2_um',V_stats3.PrincipalAxisLength_3_Max(2),'um'; ...
            'V_PrincipalAxisLength_3_um',V_stats3.PrincipalAxisLength_3_Max(3),'um'; ...
            'V_ROI_PrincipalAxisLength_1_um',V_ROI_stats3.PrincipalAxisLength_3_Max(1),'um'; ...
            'V_ROI_PrincipalAxisLength_2_um',V_ROI_stats3.PrincipalAxisLength_3_Max(2),'um'; ...
            'V_ROI_PrincipalAxisLength_3_um',V_ROI_stats3.PrincipalAxisLength_3_Max(3),'um'; ...
            'V_Solidity',V_stats3.Solidity_3_Max,'';
            'V_Solidity_CH',V_stats3.Solidity_CH_3_Max,'';
            'V_ROI_Solidity_CH',V_ROI_stats3.Solidity_CH_3_Max,'';
            'V_Extent',V_stats3.Extent_3_Max,'';
            'V_ROI_Extent',V_ROI_stats3.Extent_3_Max,'';
            'V_AspectRatio',V_stats3.PrincipalAxisLength_3_Max(1)/V_stats3.PrincipalAxisLength_3_Max(3),''; ...
            'V_ROI_AspectRatio',V_ROI_stats3.PrincipalAxisLength_3_Max(1)/V_ROI_stats3.PrincipalAxisLength_3_Max(3),''; ...
            'V_EquivDiameter',V_stats3.EquivDiameter_3_Max,'um'; ...
            'V_ROI_EquivDiameter',V_ROI_stats3.EquivDiameter_3_Max,'um'; ...
            'V_Orientation_Phi',V_stats3.Orientation_3_Max(1),''; ...
            'V_Orientation_Theta',V_stats3.Orientation_3_Max(2),''; ...
            'V_Orientation_Psi',V_stats3.Orientation_3_Max(3),''; ...
            'V_ROI_Orientation_Phi',V_ROI_stats3.Orientation_3_Max(1),''; ...
            'V_ROI_Orientation_Theta',V_ROI_stats3.Orientation_3_Max(2),''; ...
            'V_ROI_Orientation_Psi',V_ROI_stats3.Orientation_3_Max(3),''; ...
            'V_EigenValues_1',V_stats3.EigenValues_3_Max(1),''; ...
            'V_EigenValues_2',V_stats3.EigenValues_3_Max(2),''; ...
            'V_EigenValues_3',V_stats3.EigenValues_3_Max(3),''; ...
            'V_ROI_EigenValues_1',V_ROI_stats3.EigenValues_3_Max(1),''; ...
            'V_ROI_EigenValues_2',V_ROI_stats3.EigenValues_3_Max(2),''; ...
            'V_ROI_EigenValues_3',V_ROI_stats3.EigenValues_3_Max(3),''; ...
            };

        for i = 1:length(Data_stats3)
            fprintf(fileID,XRT_3DAnlys_Opts.formatSpec_s_d_s,Data_stats3{i,:});
        end  
        fclose(fileID);

        
        %% Max Feret Diameter
        XRT_3DAnlys_FeretDiameter3_Opts.ExpShorthand = XRT_3DAnlys_Opts.ExpShorthand;
        XRT_3DAnlys_FeretDiameter3_Opts.AppShorthand = XRT_3DAnlys_Opts.AppShorthand;
        [max_FeretDiameter,~] = XRT_3DAnlys_FeretDiameter3(V_ROI,XRT_3DAnlys_FeretDiameter3_Opts);
        V_max_FeretDiameter = max_FeretDiameter*(XRT_3DAnlys_Opts.ImagePixelSize*XRT_3DAnlys_Opts.DataResizeFactor);

        fileID = fopen(fullfile(XRT_3DAnlys_Opts.AppShorthand,sprintf('%s_VolumeAnlys_stats3.txt',XRT_3DAnlys_Opts.ExpShorthand)),'at');
        fprintf(fileID,XRT_3DAnlys_Opts.formatSpec_s_d_s,'V_max_FeretDiameter',V_max_FeretDiameter,'um');
        fclose(fileID);

        % % % % % % % % % % % 
        % log-file
        fprintf('%s - %s - XRT_3DAnlys_FeretDiameter3: V_max_FeretDiameter = %.2f um (Elapsed time: %.0f sec)\n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,V_max_FeretDiameter,toc(XRT_3DAnlys_Opts.tic))
        XRT_3DAnlys_Opts.PrcTime_Logger{end+1,1} = mfilename();
        XRT_3DAnlys_Opts.PrcTime_Logger{end,2} = 'XRT_3DAnlys_FeretDiameter3';
        XRT_3DAnlys_Opts.PrcTime_Logger{end,3} = toc(XRT_3DAnlys_Opts.tic);

        
        %% Wadell sphericity
        % Wadell, Hakon (1935). "Volume, Shape and Roundness of Quartz Particles". Journal of Geology. 43 (3): 250–280. doi:10.1086/624298.
        if ~exist('SurfaceArea_3_Max','var') || isnan(V_ROI_stats3.SurfaceArea_3_Max)
             stats3 = regionprops3(V_ROI,'Volume','SurfaceArea');
             V_ROI_Surface = sum(stats3.SurfaceArea) *(XRT_3DAnlys_Opts.ImagePixelSize*XRT_3DAnlys_Opts.DataResizeFactor)^2;
             V_ROI_Volume = sum(stats3.Volume) *(XRT_3DAnlys_Opts.ImagePixelSize*XRT_3DAnlys_Opts.DataResizeFactor)^3;
        else
            V_ROI_Surface = V_ROI_stats3.SurfaceArea_3_Max;
            V_ROI_Volume = V_ROI_stats3.Volume_3_Max;
        end

        V_ROI_Sphericity = (pi^(1/3) * (6*V_ROI_Volume).^(2/3))/V_ROI_Surface;

        fileID = fopen(fullfile(XRT_3DAnlys_Opts.AppShorthand,sprintf('%s_VolumeAnlys_stats3.txt',XRT_3DAnlys_Opts.ExpShorthand)),'at');
        fprintf(fileID,XRT_3DAnlys_Opts.formatSpec_s_d_s,'V_ROI_Sphericity',V_ROI_Sphericity,'');
        fclose(fileID);

        % % % % % % % % % % % 
        % log-file
        fprintf('%s - %s - V_ROI_Sphericity %.2f (Elapsed time: %.0f sec)\n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,V_ROI_Sphericity,toc(XRT_3DAnlys_Opts.tic))
        XRT_3DAnlys_Opts.PrcTime_Logger{end+1,1} = mfilename();
        XRT_3DAnlys_Opts.PrcTime_Logger{end,2} = 'V_ROI_Sphericity';
        XRT_3DAnlys_Opts.PrcTime_Logger{end,3} = toc(XRT_3DAnlys_Opts.tic);

        % % % % % % % % % % % 
        % log-file
        fprintf('%s - %s - XRT_3DAnlys_Descriptors COMPLETE (Elapsed time: %.0f sec)\n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,toc(XRT_3DAnlys_Opts.tic))
        XRT_3DAnlys_Opts.PrcTime_Logger{end+1,1} = mfilename();
        XRT_3DAnlys_Opts.PrcTime_Logger{end,2} = 'XRT_3DAnlys_Descriptors';
        XRT_3DAnlys_Opts.PrcTime_Logger{end,3} = toc(XRT_3DAnlys_Opts.tic);

    else
        V_stats.Area = nan;
        V_ROI_stats.Area = nan;
        V_stats.Surface = nan;
        V_ROI_stats.Surface = nan;
        V_stats.BoundingBox = [nan,nan];
        V_ROI_stats.BoundingBox = [nan,nan];
        V_stats.radius_spherefit = nan;
        V_ROI_stats.radius_spherefit = nan;
        V_stats.radii_elpsfit = [nan,nan,nan];
        V_ROI_stats.radii_elpsfit = [nan,nan,nan];
        V_stats.chi2_elpsfit = nan;
        V_ROI_stats.chi2_elpsfit = nan;
        V_stats.AspectRatio = nan;
        V_ROI_stats.AspectRatio = nan;

        V_stats3.Volume_3_Max = nan;
        V_stats3.BoundingBox_3_Max = [nan,nan];
        V_stats3.SurfaceArea_3_Max = nan;
        V_stats3.Solidity_3_Max = nan;
        V_stats3.Solidity_CH_3_Max = nan;
        V_stats3.Extent_3_Max = nan;
        V_stats3.Centroid_3_Max = nan;
        V_stats3.ConvexVolume_3_Max = nan;
        V_stats3.PrincipalAxisLength_3_Max = [nan,nan,nan];
        V_stats3.EquivDiameter_3_Max	= nan;

        V_stats3.Orientation_3_Max(1) = nan;
        V_stats3.Orientation_3_Max(2) = nan;
        V_stats3.Orientation_3_Max(3) = nan;

        V_stats3.EigenValues_3_Max(1) = nan;
        V_stats3.EigenValues_3_Max(2) = nan;
        V_stats3.EigenValues_3_Max(3) = nan;

        V_ROI_stats3.Volume_3_Max = nan;
        V_ROI_stats3.BoundingBox_3_Max = [nan,nan];
        V_ROI_stats3.SurfaceArea_3_Max = nan;
        V_ROI_stats3.Solidity_CH_3_Max = nan;
        V_ROI_stats3.Extent_3_Max = nan;
        V_ROI_stats3.Centroid_3_Max = nan;
        V_ROI_stats3.ConvexVolume_3_Max = nan;
        V_ROI_stats3.PrincipalAxisLength_3_Max = nan;
        V_ROI_stats3.EquivDiameter_3_Max	= nan;

        V_ROI_stats3.Orientation_3_Max(1) = nan;
        V_ROI_stats3.Orientation_3_Max(2) = nan;
        V_ROI_stats3.Orientation_3_Max(3) = nan;

        V_ROI_stats3.EigenValues_3_Max(1) = nan;
        V_ROI_stats3.EigenValues_3_Max(2) = nan;
        V_ROI_stats3.EigenValues_3_Max(3) = nan;

        size_ImageBox_X = nan;
        size_ImageBox_Y = nan;
        size_ImageBox_Z	= nan;

        V_max_FeretDiameter = nan;
        V_ROI_max_FeretDiameter = nan;
        V_imMinkowski_breadth = nan;
        V_ROI_imMinkowski_breadth = nan;
        V_Sphericity = nan;
        V_ROI_Sphericity = nan;
    end


    if isfield(XRT_3DAnlys_Opts,'VolumeAnlys_Concave_Anlys') && ~isnan(XRT_3DAnlys_Opts.VolumeAnlys_Concave_Anlys) && XRT_3DAnlys_Opts.VolumeAnlys_Concave_Anlys
        fprintf('%s - %s - XRT_3DAnlys_ConvexHull3 INITIATE (Elapsed time: %.0f sec)\n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,toc(XRT_3DAnlys_Opts.tic))
        
        XRT_3DAnlys_ConvexHull3_Opts.filesPath_FileExchange = XRT_Opts.filesPath_FileExchange;
        [V_ROI_convhull] = XRT_3DAnlys_ConvexHull3(V_ROI,XRT_3DAnlys_ConvexHull3_Opts);
        V_Concave = immultiply(V_ROI_convhull,~V_ROI);

        % Save Stack V_ROI_convhull
        ImgStack_Create_Opts.ExpShorthand = XRT_3DAnlys_Opts.ExpShorthand;
        ImgStack_Create_Opts.AppShorthand = XRT_3DAnlys_Opts.AppShorthand;
        ImgStack_Create_Opts.Save_IndImages_On = true;
        ImgStack_Create_Opts.Save_IndImages_NRecon_On = false;
        ImgStack_Create_Opts.path_ImgFolder = fullfile(XRT_3DAnlys_Opts.Output_FolderPath,'XRT_3DAnlys_ImgStack_V_ROI_Convhull_logical');
        ImgStack_Create_Opts.TiffStack_Name_string = 'V_ROI_Convhull';
        ImgStack_Create_Opts.ImgFormat = 'bmp';
        ImgStack_Create_Opts.pool_mode = true;
        ImgStack_Create_Opts.progressPrompt_On = false;
        ImgStack_Create_Opts.ImageFormat = 'logical';
        addpath(genpath(fullfile(XRT_Opts.filesPath_Mat_SR)));
        XRT_3DAnlys_ImgStack_Create(V_ROI_convhull,ImgStack_Create_Opts);

        % Save Stack V_Concave
        ImgStack_Create_Opts.ExpShorthand = XRT_3DAnlys_Opts.ExpShorthand;
        ImgStack_Create_Opts.AppShorthand = XRT_3DAnlys_Opts.AppShorthand;
        ImgStack_Create_Opts.Save_IndImages_On = true;
        ImgStack_Create_Opts.Save_IndImages_NRecon_On = false;
        ImgStack_Create_Opts.path_ImgFolder = fullfile(XRT_3DAnlys_Opts.Output_FolderPath,'XRT_3DAnlys_ImgStack_V_Concave_logical');
        ImgStack_Create_Opts.ImgName = 'V_Concave';
        ImgStack_Create_Opts.ImgFormat = 'bmp';
        ImgStack_Create_Opts.pool_mode = true;
        ImgStack_Create_Opts.progressPrompt_On = false;
        ImgStack_Create_Opts.ImageFormat = 'logical';
        addpath(genpath(fullfile(XRT_Opts.filesPath_Mat_SR)));
        XRT_3DAnlys_ImgStack_Create(V_Concave,ImgStack_Create_Opts);

        V_Concave_sumVolume_um = sum(V_Concave(:))*XRT_3DAnlys_Opts.ImagePixelSize_rS^3;
        V_ROI_convhull_um = sum(V_ROI_convhull(:))*XRT_3DAnlys_Opts.ImagePixelSize_rS^3;

        fileID = fopen(fullfile(XRT_3DAnlys_Opts.AppShorthand,sprintf('%s_VolumeAnlys_stats3.csv',XRT_3DAnlys_Opts.ExpShorthand)),'at');
        fprintf(fileID,XRT_3DAnlys_Opts.formatSpec_s_d_s,'V_Concave_sumVolume',V_Concave_sumVolume_um,'um^3');
        fprintf(fileID,XRT_3DAnlys_Opts.formatSpec_s_d_s,'V_ROI_convhull_um',V_ROI_convhull_um,'um^3');
        fclose(fileID);

        % % % % % % % % % % % 
        % log-file
        fprintf('%s - %s - XRT_3DAnlys_ConvexHull3: %.2d um^3 (Elapsed time: %.0f sec)\n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,V_ROI_convhull_um,toc(XRT_3DAnlys_Opts.tic))
        XRT_3DAnlys_Opts.PrcTime_Logger{end+1,1} = mfilename();
        XRT_3DAnlys_Opts.PrcTime_Logger{end,2} = 'XRT_3DAnlys_ConvexHull3';
        XRT_3DAnlys_Opts.PrcTime_Logger{end,3} = toc(XRT_3DAnlys_Opts.tic);
    else
        V_ROI_convhull_um = nan;
        fprintf('%s - %s - XRT_3DAnlys_ConvexHull3: >>> PASSED <<< (Elapsed time: %.0f sec)\n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,toc(XRT_3DAnlys_Opts.tic))
    end

        
    %% Update XRT_3DAnlys_MasterDB        
    idx_Col = strcmp('size_ImageBox_X',XRT_3DAnlys_MasterDB(1,:));
    XRT_3DAnlys_MasterDB{XRT_3DAnlys_Opts.idx_newRow, idx_Col} = size_ImageBox_X;
    idx_Col = strcmp('size_ImageBox_Y',XRT_3DAnlys_MasterDB(1,:));
    XRT_3DAnlys_MasterDB{XRT_3DAnlys_Opts.idx_newRow, idx_Col} = size_ImageBox_Y;
    idx_Col = strcmp('size_ImageBox_Z',XRT_3DAnlys_MasterDB(1,:));
    XRT_3DAnlys_MasterDB{XRT_3DAnlys_Opts.idx_newRow, idx_Col} = size_ImageBox_Z;
    idx_Col = strcmp('V_BoundingBox_MaxLength',XRT_3DAnlys_MasterDB(1,:));
    XRT_3DAnlys_MasterDB{XRT_3DAnlys_Opts.idx_newRow, idx_Col} = V_stats.BoundingBox(1);
    idx_Col = strcmp('V_BoundingBox_MinLength',XRT_3DAnlys_MasterDB(1,:));
    XRT_3DAnlys_MasterDB{XRT_3DAnlys_Opts.idx_newRow, idx_Col} = V_stats.BoundingBox(2);
    idx_Col = strcmp('V_ROI_BoundingBox_MaxLength',XRT_3DAnlys_MasterDB(1,:));
    XRT_3DAnlys_MasterDB{XRT_3DAnlys_Opts.idx_newRow, idx_Col} = V_ROI_stats.BoundingBox(1);
    idx_Col = strcmp('V_ROI_BoundingBox_MinLength',XRT_3DAnlys_MasterDB(1,:));
    XRT_3DAnlys_MasterDB{XRT_3DAnlys_Opts.idx_newRow, idx_Col} = V_ROI_stats.BoundingBox(2);
    idx_Col = strcmp('NumEl_Volume',XRT_3DAnlys_MasterDB(1,:));
    XRT_3DAnlys_MasterDB{XRT_3DAnlys_Opts.idx_newRow, idx_Col} = numel(V_stats.Area);
    idx_Col = strcmp('V_Solidity',XRT_3DAnlys_MasterDB(1,:));
    XRT_3DAnlys_MasterDB{XRT_3DAnlys_Opts.idx_newRow, idx_Col} = V_stats3.Solidity_3_Max;
    idx_Col = strcmp('V_Solidity_CH',XRT_3DAnlys_MasterDB(1,:));
    XRT_3DAnlys_MasterDB{XRT_3DAnlys_Opts.idx_newRow, idx_Col} = V_stats3.Solidity_CH_3_Max;
    idx_Col = strcmp('V_ROI_Solidity_CH',XRT_3DAnlys_MasterDB(1,:));
    XRT_3DAnlys_MasterDB{XRT_3DAnlys_Opts.idx_newRow, idx_Col} = V_ROI_stats3.Solidity_CH_3_Max;
    idx_Col = strcmp('V_Extent',XRT_3DAnlys_MasterDB(1,:));
    XRT_3DAnlys_MasterDB{XRT_3DAnlys_Opts.idx_newRow, idx_Col} = V_stats3.Extent_3_Max;
    idx_Col = strcmp('V_ROI_Extent',XRT_3DAnlys_MasterDB(1,:));
    XRT_3DAnlys_MasterDB{XRT_3DAnlys_Opts.idx_newRow, idx_Col} = V_ROI_stats3.Extent_3_Max;
    idx_Col = strcmp('V_Volume_um',XRT_3DAnlys_MasterDB(1,:));
    XRT_3DAnlys_MasterDB{XRT_3DAnlys_Opts.idx_newRow, idx_Col} = sum(V_stats.Area);
    idx_Col = strcmp('V_ROI_Volume_um',XRT_3DAnlys_MasterDB(1,:));
    XRT_3DAnlys_MasterDB{XRT_3DAnlys_Opts.idx_newRow, idx_Col} = sum(V_ROI_stats.Area);
    idx_Col = strcmp('V_ROI_CH_Volume_um',XRT_3DAnlys_MasterDB(1,:));
    XRT_3DAnlys_MasterDB{XRT_3DAnlys_Opts.idx_newRow, idx_Col} = V_ROI_convhull_um;
    idx_Col = strcmp('V_Surface_um',XRT_3DAnlys_MasterDB(1,:));
    XRT_3DAnlys_MasterDB{XRT_3DAnlys_Opts.idx_newRow, idx_Col} = V_stats.Surface;
    idx_Col = strcmp('V_ROI_Surface_um',XRT_3DAnlys_MasterDB(1,:));
    XRT_3DAnlys_MasterDB{XRT_3DAnlys_Opts.idx_newRow, idx_Col} = V_ROI_stats.Surface;
    idx_Col = strcmp('V_radius_spherefit_um',XRT_3DAnlys_MasterDB(1,:));
    XRT_3DAnlys_MasterDB{XRT_3DAnlys_Opts.idx_newRow, idx_Col} = V_stats.radius_spherefit;
    idx_Col = strcmp('V_ROI_radius_spherefit_um',XRT_3DAnlys_MasterDB(1,:));
    XRT_3DAnlys_MasterDB{XRT_3DAnlys_Opts.idx_newRow, idx_Col} = V_ROI_stats.radius_spherefit;
    idx_Col = strcmp('V_radii_elpsfit_1_um',XRT_3DAnlys_MasterDB(1,:));
    XRT_3DAnlys_MasterDB{XRT_3DAnlys_Opts.idx_newRow, idx_Col} = V_stats.radii_elpsfit(1);
    idx_Col = strcmp('V_radii_elpsfit_2_um',XRT_3DAnlys_MasterDB(1,:));
    XRT_3DAnlys_MasterDB{XRT_3DAnlys_Opts.idx_newRow, idx_Col} = V_stats.radii_elpsfit(2);
    idx_Col = strcmp('V_radii_elpsfit_3_um',XRT_3DAnlys_MasterDB(1,:));
    XRT_3DAnlys_MasterDB{XRT_3DAnlys_Opts.idx_newRow, idx_Col} = V_stats.radii_elpsfit(3);
    idx_Col = strcmp('V_chi2_elpsfit_um',XRT_3DAnlys_MasterDB(1,:));
    XRT_3DAnlys_MasterDB{XRT_3DAnlys_Opts.idx_newRow, idx_Col} = V_stats.chi2_elpsfit;
    idx_Col = strcmp('V_ROI_radii_elpsfit_1_um',XRT_3DAnlys_MasterDB(1,:));
    XRT_3DAnlys_MasterDB{XRT_3DAnlys_Opts.idx_newRow, idx_Col} = V_ROI_stats.radii_elpsfit(1);
    idx_Col = strcmp('V_ROI_radii_elpsfit_2_um',XRT_3DAnlys_MasterDB(1,:));
    XRT_3DAnlys_MasterDB{XRT_3DAnlys_Opts.idx_newRow, idx_Col} = V_ROI_stats.radii_elpsfit(2);
    idx_Col = strcmp('V_ROI_radii_elpsfit_3_um',XRT_3DAnlys_MasterDB(1,:));
    XRT_3DAnlys_MasterDB{XRT_3DAnlys_Opts.idx_newRow, idx_Col} = V_ROI_stats.radii_elpsfit(3);
    idx_Col = strcmp('V_ROI_chi2_elpsfit_um',XRT_3DAnlys_MasterDB(1,:));
    XRT_3DAnlys_MasterDB{XRT_3DAnlys_Opts.idx_newRow, idx_Col} = V_ROI_stats.chi2_elpsfit;
    idx_Col = strcmp('V_AspectRatio',XRT_3DAnlys_MasterDB(1,:));
    XRT_3DAnlys_MasterDB{XRT_3DAnlys_Opts.idx_newRow, idx_Col} = V_stats.AspectRatio;
    idx_Col = strcmp('V_ROI_AspectRatio',XRT_3DAnlys_MasterDB(1,:));
    XRT_3DAnlys_MasterDB{XRT_3DAnlys_Opts.idx_newRow, idx_Col} = V_ROI_stats.AspectRatio;
    idx_Col = strcmp('V_ROI_Sphericity',XRT_3DAnlys_MasterDB(1,:));
    XRT_3DAnlys_MasterDB{XRT_3DAnlys_Opts.idx_newRow, idx_Col} = V_ROI_Sphericity;
    idx_Col = strcmp('V_max_FeretDiameter',XRT_3DAnlys_MasterDB(1,:));
    XRT_3DAnlys_MasterDB{XRT_3DAnlys_Opts.idx_newRow, idx_Col} = V_max_FeretDiameter;

    idx_Col = strcmp('V_EquivDiameter',XRT_3DAnlys_MasterDB(1,:));
    XRT_3DAnlys_MasterDB{XRT_3DAnlys_Opts.idx_newRow, idx_Col} = V_stats3.EquivDiameter_3_Max;
    idx_Col = strcmp('V_ROI_EquivDiameter',XRT_3DAnlys_MasterDB(1,:));
    XRT_3DAnlys_MasterDB{XRT_3DAnlys_Opts.idx_newRow, idx_Col} = V_ROI_stats3.EquivDiameter_3_Max;

    idx_Col = strcmp('V_Orientation_Phi',XRT_3DAnlys_MasterDB(1,:));
    XRT_3DAnlys_MasterDB{XRT_3DAnlys_Opts.idx_newRow, idx_Col} = V_stats3.Orientation_3_Max(1);
    idx_Col = strcmp('V_Orientation_Theta',XRT_3DAnlys_MasterDB(1,:));
    XRT_3DAnlys_MasterDB{XRT_3DAnlys_Opts.idx_newRow, idx_Col} = V_stats3.Orientation_3_Max(2);
    idx_Col = strcmp('V_Orientation_Psi',XRT_3DAnlys_MasterDB(1,:));
    XRT_3DAnlys_MasterDB{XRT_3DAnlys_Opts.idx_newRow, idx_Col} = V_stats3.Orientation_3_Max(3);
    idx_Col = strcmp('V_ROI_Orientation_Phi',XRT_3DAnlys_MasterDB(1,:));
    XRT_3DAnlys_MasterDB{XRT_3DAnlys_Opts.idx_newRow, idx_Col} = V_ROI_stats3.Orientation_3_Max(1);
    idx_Col = strcmp('V_ROI_Orientation_Theta',XRT_3DAnlys_MasterDB(1,:));
    XRT_3DAnlys_MasterDB{XRT_3DAnlys_Opts.idx_newRow, idx_Col} = V_ROI_stats3.Orientation_3_Max(2);
    idx_Col = strcmp('V_ROI_Orientation_Psi',XRT_3DAnlys_MasterDB(1,:));
    XRT_3DAnlys_MasterDB{XRT_3DAnlys_Opts.idx_newRow, idx_Col} = V_ROI_stats3.Orientation_3_Max(3);

        
    %% End Sequence
    fprintf('%s - %s - %s COMPLETE (Elapsed time: %.0f sec)\n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename(),toc(XRT_3DAnlys_Opts.tic))
    XRT_3DAnlys_Opts.PrcTime_Logger{end+1,1} = mfilename();
    XRT_3DAnlys_Opts.PrcTime_Logger{end,2} = sprintf('%s COMPLETE',mfilename());
    XRT_3DAnlys_Opts.PrcTime_Logger{end,3} = toc(XRT_3DAnlys_Opts.tic);
end

