%------------------------------------------------------------------------------------------------
% Code written by Frederik Doerr, Aug 2020 (MATLAB R2019b)
% Application: XRT Image and Data Analysis Framework (SubRoutine
% XRT_3DAnlys_PorosAnlys)
% 
% https://github.com/frederik-d
% Contact: frederik.doerr(at)strath.ac.uk | CMAC (http://www.cmac.ac.uk/)


clearvars -except XRT_3DAnlys_Opts XRT_3DAnlys_MasterDB Batch3DAnlys_Iter numExp pool XRT_Opts
XRT_ImgPrc_Opts = struct();
[XRT_3DAnlys_Opts,XRT_ImgPrc_Opts] = XRT_3DAnlys_ParameterLoader(XRT_Opts,Batch3DAnlys_Iter,XRT_3DAnlys_Opts,XRT_ImgPrc_Opts);
XRT_3DAnlys_Opts.AppShorthand = mfilename();

if XRT_3DAnlys_Opts.PorosAnlys3D_On
    % log-file
    fprintf('- - - - - - - - - - - -\n')
    fprintf('%s - Enter %s (Elapsed time: %.0f sec)\n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,toc(XRT_3DAnlys_Opts.tic))
    fprintf('- - - - - - - - - - - -\n')
    pause(2)
    
    XRT_3DAnlys_Opts.PrcTime_Logger{end+1,1} = XRT_3DAnlys_Opts.AppShorthand;
    XRT_3DAnlys_Opts.PrcTime_Logger{end,2} = 'INITIATE';
    XRT_3DAnlys_Opts.PrcTime_Logger{end,3} = toc(XRT_3DAnlys_Opts.tic);
    
    % % % % % % % % % % % 0) XRT_3DAnlys_ReconImage_VolumeBuild % % % % % % % % % % %
    V_binarize_Check = true;
    run XRT_3DAnlys_VolumeBuild
    % % % % % % % % % % % 
    % log-file
    fprintf('%s - %s - XRT_3DAnlys_VolumeBuild Completed (Elapsed time: %.0f sec)\n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,toc(XRT_3DAnlys_Opts.tic))
    XRT_3DAnlys_Opts.PrcTime_Logger{end+1,1} = XRT_3DAnlys_Opts.AppShorthand;
    XRT_3DAnlys_Opts.PrcTime_Logger{end,2} = 'XRT_3DAnlys_VolumeBuild';
    XRT_3DAnlys_Opts.PrcTime_Logger{end,3} = toc(XRT_3DAnlys_Opts.tic);
    

    XRT_3DAnlys_Opts.filesPath_PorosAnlys = fullfile(XRT_3DAnlys_Opts.Output_FolderPath,XRT_3DAnlys_Opts.AppShorthand);
    if ~exist(XRT_3DAnlys_Opts.filesPath_PorosAnlys,'dir')
        mkdir(XRT_3DAnlys_Opts.filesPath_PorosAnlys)
    end
    
    %% (1) Calculate open - / closed - porosity
    V_Poros = immultiply(V_ROI,~V);
    V_fill = imfill(V,'holes');

    V_Poros_closed = V_fill - V;
    V_Poros_open = immultiply(V_ROI,~V_fill);   
   
    %% (2) Calculate Global 3DPorosity Descriptors

    % % % % % % % Global 3DPorosity Descriptors % % % % % % %
    V_Poros_Volume_stats = regionprops(V_Poros,'PixelIdxList','Area');
    V_Poros_Volume_list = [V_Poros_Volume_stats.Area];
    V_Poros_Volume_list = sort(V_Poros_Volume_list, 'descend');
    V_Poros_Volume = sum([V_Poros_Volume_stats.Area])*XRT_3DAnlys_Opts.ImagePixelSize_rS^3;

    V_Poros_closed_Volume_stats = regionprops(V_Poros_closed,'PixelIdxList','Area');
    V_Poros_closed_Volume_list = [V_Poros_closed_Volume_stats.Area];
    V_Poros_closed_Volume_list = sort(V_Poros_closed_Volume_list, 'descend');
    V_Poros_closed_Volume = sum([V_Poros_closed_Volume_stats.Area])*XRT_3DAnlys_Opts.ImagePixelSize_rS^3;

    V_Poros_open_Volume_stats = regionprops(V_Poros_open,'PixelIdxList','Area');
    V_Poros_open_Volume_list = [V_Poros_open_Volume_stats.Area];
    V_Poros_open_Volume_list = sort(V_Poros_open_Volume_list, 'descend');
    V_Poros_open_Volume = sum([V_Poros_open_Volume_stats.Area])*XRT_3DAnlys_Opts.ImagePixelSize_rS^3; 

    V_ROI_Volume_stats = regionprops(V_ROI,'PixelIdxList','Area');
    V_ROI_Volume_list = [V_ROI_Volume_stats.Area];
    V_ROI_Volume_list = sort(V_ROI_Volume_list, 'descend');
    V_ROI_Volume = sum([V_ROI_Volume_stats.Area])*XRT_3DAnlys_Opts.ImagePixelSize_rS^3;

    V_Volume_stats = regionprops(V,'PixelIdxList','Area');
    V_Volume_list = [V_Volume_stats.Area];
    V_Volume_list = sort(V_Volume_list, 'descend');
    V_Volume = sum([V_Volume_stats.Area])*XRT_3DAnlys_Opts.ImagePixelSize_rS^3;

    % Volume Fractions:
    V_Poros_VolumeFrac = V_Poros_Volume/V_ROI_Volume;
    V_Poros_closed_VolumeFrac = V_Poros_closed_Volume/V_ROI_Volume;
    V_Poros_open_VolumeFrac = V_Poros_open_Volume/V_ROI_Volume;
    V_VolumeFrac = V_Volume/V_ROI_Volume;
    
    
    % Export 2D Images
    XRT_3DAnlys_ImgExport_Opts.ExpShorthand = XRT_3DAnlys_Opts.ExpShorthand;
    XRT_3DAnlys_ImgExport_Opts.AppShorthand = XRT_3DAnlys_Opts.AppShorthand;
    XRT_3DAnlys_ImgExport_Opts.ImagePixelSize = XRT_3DAnlys_Opts.ImagePixelSize;
    XRT_3DAnlys_ImgExport_Opts.DataResizeFactor = XRT_3DAnlys_Opts.DataResizeFactor;
    XRT_3DAnlys_ImgExport_Opts.path_ImgFolder = fullfile(XRT_3DAnlys_Opts.Output_FolderPath,XRT_3DAnlys_Opts.AppShorthand);
    XRT_3DAnlys_ImgExport_Opts.ImgName = sprintf('%s_V_Poros',XRT_3DAnlys_Opts.ExpShorthand);
    XRT_3DAnlys_ImgExport(V_Poros,XRT_3DAnlys_ImgExport_Opts)
    XRT_3DAnlys_ImgExport_Opts.ImgName = sprintf('%s_V_Poros_closed',XRT_3DAnlys_Opts.ExpShorthand);
    XRT_3DAnlys_ImgExport(V_Poros_closed,XRT_3DAnlys_ImgExport_Opts)
    XRT_3DAnlys_ImgExport_Opts.ImgName = sprintf('%s_V_Poros_open',XRT_3DAnlys_Opts.ExpShorthand);
    XRT_3DAnlys_ImgExport(V_Poros_open,XRT_3DAnlys_ImgExport_Opts)
    
    
    %% Export data
    fileID = fopen(fullfile(XRT_3DAnlys_Opts.filesPath_PorosAnlys,sprintf('%s_%s_Descriptors.csv',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand)),'wt');
    Data_XRTAnlys_3DAnlys_dat = { ...
        'V_Volume_um',V_Volume,'um^3'; ...
        'V_ROI_Volume_um',V_ROI_Volume,'um^3'; ...
        'V_Poros_Volume_um',V_Poros_Volume,'um^3'; ...
        'V_Poros_closed_Volume_um',V_Poros_closed_Volume,'um^3'; ...
        'V_Poros_open_Volume_um',V_Poros_open_Volume,'um^3'; ...
        'V_VolumeFrac',V_VolumeFrac,'-'; ...
        'V_Poros_VolumeFrac',V_Poros_VolumeFrac,'-'; ...
        'V_Poros_closed_VolumeFrac',V_Poros_closed_VolumeFrac,'-'; ...
        'V_Poros_open_VolumeFrac',V_Poros_open_VolumeFrac,'-'};
    for i = 1:length(Data_XRTAnlys_3DAnlys_dat)
        fprintf(fileID,XRT_3DAnlys_Opts.formatSpec_s_d_s,Data_XRTAnlys_3DAnlys_dat{i,:});
    end  
    fclose(fileID);
    
    
    DataDim = max([length(V_Volume_list),length(V_ROI_Volume_list), ...
        length(V_Poros_Volume_list),length(V_Poros_open_Volume_list),length(V_Poros_closed_Volume_list)]);
    Data_3DPorosAnlys = zeros(DataDim,5);
    % Export Volume Lists
    Header_3DPorosAnlys = {'V_Volume_list','V_ROI_Volume_list','V_Poros_Volume_list','V_Poros_open_Volume_list','V_Poros_closed_Volume_list'}; 
    Data_3DPorosAnlys(1:length(V_Volume_list),1) = V_Volume_list;
    Data_3DPorosAnlys(1:length(V_ROI_Volume_list),2) = V_ROI_Volume_list;
    Data_3DPorosAnlys(1:length(V_Poros_Volume_list),3) = V_Poros_Volume_list;
    Data_3DPorosAnlys(1:length(V_Poros_open_Volume_list),4) = V_Poros_open_Volume_list;
    Data_3DPorosAnlys(1:length(V_Poros_closed_Volume_list),5) = V_Poros_closed_Volume_list;
    fid = fopen(fullfile(XRT_3DAnlys_Opts.filesPath_PorosAnlys,sprintf('%s_Data.csv',XRT_3DAnlys_Opts.AppShorthand)), 'w') ;
    fprintf(fid, '%s,', Header_3DPorosAnlys{1,1:end-1}) ;
    fprintf(fid, '%s\n', Header_3DPorosAnlys{1,end}) ;
    fclose(fid) ;
    dlmwrite(fullfile(XRT_3DAnlys_Opts.filesPath_PorosAnlys,sprintf('%s_Data.csv',XRT_3DAnlys_Opts.AppShorthand)), Data_3DPorosAnlys(1:end,:), '-append');
   

    %% Save ImgStacks
    ImgStack_Create_Options.pool_mode = true;
    ImgStack_Create_Options.Print_Log_On = false;
    ImgStack_Create_Options.Save_IndImages_On = true;
    ImgStack_Create_Options.ExpShorthand = XRT_3DAnlys_Opts.ExpShorthand;
    ImgStack_Create_Options.ImgFormat = 'bmp';

    ImgStack_Create_Options.ImgName = 'V_Poros';
    ImgStack_Create_Options.AppShorthand = sprintf('%s_%s',XRT_3DAnlys_Opts.AppShorthand,ImgStack_Create_Options.ImgName);
    ImgStack_Create_Options.path_ImgFolder = fullfile(XRT_3DAnlys_Opts.Output_FolderPath,sprintf('XRT_3DAnlys_ImgStack_%s',ImgStack_Create_Options.ImgName));
    addpath(fullfile(XRT_Opts.filesPath_Mat_SR));
    XRT_3DAnlys_ImgStack_Create(V_Poros,ImgStack_Create_Options);

    ImgStack_Create_Options.ImgName = 'V_Poros_closed';
    ImgStack_Create_Options.AppShorthand = sprintf('%s_%s',XRT_3DAnlys_Opts.AppShorthand,ImgStack_Create_Options.ImgName);
    ImgStack_Create_Options.path_ImgFolder = fullfile(XRT_3DAnlys_Opts.Output_FolderPath,sprintf('XRT_3DAnlys_ImgStack_%s',ImgStack_Create_Options.ImgName));
    addpath(fullfile(XRT_Opts.filesPath_Mat_SR));
    XRT_3DAnlys_ImgStack_Create(V_Poros_closed,ImgStack_Create_Options);

    ImgStack_Create_Options.ImgName = 'V_Poros_open';
    ImgStack_Create_Options.AppShorthand = sprintf('%s_%s',XRT_3DAnlys_Opts.AppShorthand,ImgStack_Create_Options.ImgName);
    ImgStack_Create_Options.path_ImgFolder = fullfile(XRT_3DAnlys_Opts.Output_FolderPath,sprintf('XRT_3DAnlys_ImgStack_%s',ImgStack_Create_Options.ImgName));
    addpath(fullfile(XRT_Opts.filesPath_Mat_SR));
    XRT_3DAnlys_ImgStack_Create(V_Poros_open,ImgStack_Create_Options);


    %% Update XRT_3DAnlys_MasterDB        
    idx_Col = strcmp('V_Poros_Volume_um',XRT_3DAnlys_MasterDB(1,:));
    XRT_3DAnlys_MasterDB{XRT_3DAnlys_Opts.idx_newRow, idx_Col} = V_Poros_Volume;
    idx_Col = strcmp('V_Poros_closed_Volume_um',XRT_3DAnlys_MasterDB(1,:));
    XRT_3DAnlys_MasterDB{XRT_3DAnlys_Opts.idx_newRow, idx_Col} = V_Poros_closed_Volume;
    idx_Col = strcmp('V_Poros_open_Volume_um',XRT_3DAnlys_MasterDB(1,:));
    XRT_3DAnlys_MasterDB{XRT_3DAnlys_Opts.idx_newRow, idx_Col} = V_Poros_open_Volume;


    %% End Sequence PorosAnlys3D
    fprintf('%s - %s - %s COMPLETE (Elapsed time: %.0f sec)\n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.AppShorthand,mfilename(),toc(XRT_3DAnlys_Opts.tic))
    XRT_3DAnlys_Opts.PrcTime_Logger{end+1,1} = mfilename();
    XRT_3DAnlys_Opts.PrcTime_Logger{end,2} = sprintf('%s COMPLETE',mfilename());
    XRT_3DAnlys_Opts.PrcTime_Logger{end,3} = toc(XRT_3DAnlys_Opts.tic); 
end
