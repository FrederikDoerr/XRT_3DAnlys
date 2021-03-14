%------------------------------------------------------------------------------------------------
% Code written by Frederik Doerr, Aug 2020 (MATLAB R2019b)
% Application: XRT Image and Data Analysis Framework (SubRoutine
% XRT_3DAnlys_BinSegment)
% 
% https://github.com/frederik-d
% Contact: frederik.doerr(at)strath.ac.uk | CMAC (http://www.cmac.ac.uk/)

%% XRT_3DAnlys_BinSegment
clearvars -except XRT_3DAnlys_Opts XRT_3DAnlys_MasterDB Batch3DAnlys_Iter numExp pool XRT_Opts
XRT_ImgPrc_Opts = struct();
[XRT_3DAnlys_Opts,XRT_ImgPrc_Opts] = XRT_3DAnlys_ParameterLoader(XRT_Opts,Batch3DAnlys_Iter,XRT_3DAnlys_Opts,XRT_ImgPrc_Opts);
XRT_3DAnlys_Opts.AppShorthand = mfilename();
XRT_3DAnlys_Opts.ROI_Reduction_Check = false;

if XRT_3DAnlys_Opts.BinSegment3D_On
    % log-file
    fprintf('- - - - - - - - - - - -\n')
    fprintf('%s - Enter %s (Elapsed time: %.0f sec)\n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,toc(XRT_3DAnlys_Opts.tic))
    fprintf('- - - - - - - - - - - -\n')
    pause(2)
    
    % % % % % % % % % % % 0) XRT_3DAnlys_ReconImage_VolumeBuild % % % % % % % % % % %
    V_binarize_Check = true;
    run XRT_3DAnlys_VolumeBuild
    % % % % % % % % % % % 
        
    % log-file
    fprintf('%s - %s - XRT_3DAnlys_VolumeBuild COMPLETE (Elapsed time: %.0f sec)\n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,toc(XRT_3DAnlys_Opts.tic))
    
    XRT_3DAnlys_Opts.filesPath_BinSegment = fullfile(XRT_3DAnlys_Opts.Output_FolderPath,XRT_3DAnlys_Opts.AppShorthand);
    if ~exist(XRT_3DAnlys_Opts.filesPath_BinSegment,'dir')
        mkdir(XRT_3DAnlys_Opts.filesPath_BinSegment)
    end
    
    %% % % % % % % % % % % 1) XRT_3DAnlys_FillclosedPoros (Fill small, closed inner sctructures) % % % % % % % % % % %
    V = imfill(V,'holes');

    % Export 2D Images
    ImgExport_Opts.ExpShorthand = XRT_3DAnlys_Opts.ExpShorthand;
    ImgExport_Opts.ImagePixelSize = XRT_3DAnlys_Opts.ImagePixelSize;
    ImgExport_Opts.DataResizeFactor = XRT_3DAnlys_Opts.DataResizeFactor;
    ImgExport_Opts.path_ImgFolder = XRT_3DAnlys_Opts.filesPath_BinSegment;
    ImgExport_Opts.ImgName = sprintf('%s_V_FillcP',XRT_3DAnlys_Opts.ExpShorthand);
    XRT_3DAnlys_ImgExport(V,ImgExport_Opts)
    
    
    % Save image stack
    ImgStack_Create_Opts.TiffStack_Name_string = 'V_FillcP';
    ImgStack_Create_Opts.pool_mode = true;
    ImgStack_Create_Opts.ExpShorthand = XRT_3DAnlys_Opts.ExpShorthand;
    ImgStack_Create_Opts.AppShorthand = XRT_3DAnlys_Opts.AppShorthand;
    ImgStack_Create_Opts.Save_IndImages_On = true;
    ImgStack_Create_Opts.Save_IndImages_NRecon_On = false;
    ImgStack_Create_Opts.progressPrompt_On = false;
    ImgStack_Create_Opts.ImagePixelSize = XRT_3DAnlys_Opts.ImagePixelSize;
    ImgStack_Create_Opts.DataResizeFactor = XRT_3DAnlys_Opts.DataResizeFactor;
    ImgStack_Create_Opts.ImgName = 'V_FillcP';
    ImgStack_Create_Opts.path_ImgFolder = fullfile(XRT_3DAnlys_Opts.Output_FolderPath,'XRT_3DAnlys_ImgStack_V_FillcP');
    XRT_3DAnlys_ImgStack_Create(V,ImgStack_Create_Opts);

    fprintf('%s - %s - XRT_3DAnlys_FillclosedPoros (imfill) COMPLETE (Elapsed time: %.0f sec)\n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,toc(XRT_3DAnlys_Opts.tic))

    clearvars -except V XRT_3DAnlys_Opts XRT_ImgPrc_Opts XRT_3DAnlys_MasterDB XRT_Opts Batch3DAnlys_Iter numExp pool

    %% % % % % % % % % % % 2) XRT_3DAnlys_Shrinkwrap (Fill small, open inner sctructures) % % % % % % % % % % %
    XRT_3DAnlys_Shrinkwrap_Opts.AppShorthand = XRT_3DAnlys_Opts.ExpShorthand;
    XRT_3DAnlys_Shrinkwrap_Opts.AppShorthand = XRT_3DAnlys_Opts.AppShorthand;
    XRT_3DAnlys_Shrinkwrap_Opts.subFolderName_path = fullfile(pwd,XRT_3DAnlys_Opts.AppShorthand);
    XRT_3DAnlys_Shrinkwrap_Opts.Method = XRT_ImgPrc_Opts.Watershed.Shrinkwrap_Method;
    XRT_3DAnlys_Shrinkwrap_Opts.imclose_strel_type = 'sphere';
    XRT_3DAnlys_Shrinkwrap_Opts.imclose_strel_size = XRT_ImgPrc_Opts.Watershed.Shrinkwrap_Size;
    XRT_3DAnlys_Shrinkwrap_Opts.Control_on = false;
    
    V = XRT_3DAnlys_Shrinkwrap(V,XRT_3DAnlys_Shrinkwrap_Opts);

    % Export 2D Images
    ImgExport_Opts.ExpShorthand = XRT_3DAnlys_Opts.ExpShorthand;
    ImgExport_Opts.ImagePixelSize = XRT_3DAnlys_Opts.ImagePixelSize;
    ImgExport_Opts.DataResizeFactor = XRT_3DAnlys_Opts.DataResizeFactor;
    ImgExport_Opts.path_ImgFolder = XRT_3DAnlys_Opts.filesPath_BinSegment;
    ImgExport_Opts.ImgName = sprintf('%s_V_FilloP',XRT_3DAnlys_Opts.ExpShorthand);
    XRT_3DAnlys_ImgExport(V,ImgExport_Opts)
    
    % Save image stack
    ImgStack_Create_Opts.pool_mode = true;
    ImgStack_Create_Opts.ExpShorthand = XRT_3DAnlys_Opts.ExpShorthand;
    ImgStack_Create_Opts.AppShorthand = XRT_3DAnlys_Opts.AppShorthand;
    ImgStack_Create_Opts.Save_IndImages_On = true;
    ImgStack_Create_Opts.Save_IndImages_NRecon_On = false;
    ImgStack_Create_Opts.progressPrompt_On = false;
    ImgStack_Create_Opts.ImgFormat = 'bmp';
    ImgStack_Create_Opts.path_ImgFolder = fullfile(XRT_3DAnlys_Opts.Output_FolderPath,'XRT_3DAnlys_ImgStack_V_FilloP');
    ImgStack_Create_Opts.ImgName = 'V_FilloP';
    XRT_3DAnlys_ImgStack_Create(V,ImgStack_Create_Opts);

    fprintf('%s - %s - XRT_3DAnlys_Shrinkwrap COMPLETE (Elapsed time: %.0f sec)\n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,toc(XRT_3DAnlys_Opts.tic))

    clearvars -except V XRT_3DAnlys_Opts XRT_ImgPrc_Opts XRT_3DAnlys_MasterDB XRT_Opts Batch3DAnlys_Iter numExp pool

    %% % % % % % % % % % % XRT_3DAnlys_Watershed % % % % % % % % % % %
   
    XRT_3DAnlys_Watershed_Opts.ExpShorthand = XRT_3DAnlys_Opts.ExpShorthand;
    XRT_3DAnlys_Watershed_Opts.AppShorthand = XRT_3DAnlys_Opts.AppShorthand;
    XRT_3DAnlys_Watershed_Opts.ImagePixelSize = XRT_3DAnlys_Opts.ImagePixelSize;
    XRT_3DAnlys_Watershed_Opts.DataResizeFactor = XRT_3DAnlys_Opts.DataResizeFactor;
    
%     XRT_3DAnlys_Watershed_Opts.subFolderName = fullfile(XRT_3DAnlys_Opts.AppShorthand,sprintf('XRT_3DAnlys_SegmentWatershed_ImgExpt'));
    XRT_3DAnlys_Watershed_Opts.imextendedmin_Size = XRT_ImgPrc_Opts.Watershed.imextendedmin_Size;
    XRT_3DAnlys_Watershed_Opts.bwdist_method = XRT_ImgPrc_Opts.Watershed.bwdist_method;
    XRT_3DAnlys_Watershed_Opts.Control_on = true;
    XRT_3DAnlys_Watershed_Opts.R_Ctrl_On = true;
    [V,R] = XRT_3DAnlys_Watershed(V,XRT_3DAnlys_Watershed_Opts);
    % % % % % % % % % % % 
    
    % Export 2D Images
    ImgExport_Opts.ExpShorthand = XRT_3DAnlys_Opts.ExpShorthand;
    ImgExport_Opts.ImagePixelSize = XRT_3DAnlys_Opts.ImagePixelSize;
    ImgExport_Opts.DataResizeFactor = XRT_3DAnlys_Opts.DataResizeFactor;
    ImgExport_Opts.path_ImgFolder = XRT_3DAnlys_Opts.filesPath_BinSegment;
    ImgExport_Opts.ImgName = sprintf('%s_V_Segment',XRT_3DAnlys_Opts.ExpShorthand);
    XRT_3DAnlys_ImgExport(V,ImgExport_Opts)
    
    
    % Save image stack
    ImgStack_Create_Opts.pool_mode = true;
    ImgStack_Create_Opts.progressPrompt_On = false;
    ImgStack_Create_Opts.ExpShorthand = XRT_3DAnlys_Opts.ExpShorthand;
    ImgStack_Create_Opts.AppShorthand = XRT_3DAnlys_Opts.AppShorthand;
    ImgStack_Create_Opts.Save_IndImages_On = true;
    ImgStack_Create_Opts.Save_IndImages_NRecon_On = false;
    ImgStack_Create_Opts.ImgFormat = 'bmp';
    
    ImgStack_Create_Opts.path_ImgFolder = fullfile(XRT_3DAnlys_Opts.Output_FolderPath,'XRT_3DAnlys_ImgStack_D');
    ImgStack_Create_Opts.ImgName = 'D';
    XRT_3DAnlys_ImgStack_Create(R.D_Img,ImgStack_Create_Opts);
    
    ImgStack_Create_Opts.path_ImgFolder = fullfile(XRT_3DAnlys_Opts.Output_FolderPath,'XRT_3DAnlys_ImgStack_marker');
    ImgStack_Create_Opts.ImgName = 'marker';
    XRT_3DAnlys_ImgStack_Create(R.marker_Img,ImgStack_Create_Opts);
    
    ImgStack_Create_Opts.path_ImgFolder = fullfile(XRT_3DAnlys_Opts.Output_FolderPath,'XRT_3DAnlys_ImgStack_Ld');
    ImgStack_Create_Opts.ImgName = 'Ld';
    XRT_3DAnlys_ImgStack_Create(R.Ld_Img,ImgStack_Create_Opts);
    
    ImgStack_Create_Opts.path_ImgFolder = fullfile(XRT_3DAnlys_Opts.Output_FolderPath,'XRT_3DAnlys_ImgStack_V_Segmented');
    ImgStack_Create_Opts.ImgName = 'V_Segmented';
    XRT_3DAnlys_ImgStack_Create(V,ImgStack_Create_Opts);
    
    % log-file
    fprintf('%s - %s - XRT_3DAnlys_Watershed COMPLETE (Elapsed time: %.0f sec)\n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,toc(XRT_3DAnlys_Opts.tic))

end