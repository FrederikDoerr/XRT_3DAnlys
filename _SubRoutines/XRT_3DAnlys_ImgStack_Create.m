function XRT_3DAnlys_ImgStack_Create(V,Opt)
%XRT_3DAnlys_ImgStack_Create Write image data to repository.
%   [V] = XRT_3DAnlys_ImgStack_Create(V,Opt) function to ease the routine 
%   exporting of image data to folders.
%   Potential to use parallel processing (Opt.pool_mode = true, see parfor)
%
%   OPTIONS:
%      'path_ImgFolder'     Image folder path.
%      'ImgFormat'          Image file format (default: 'bmp').
%      'ImgName'            Image file name.
%      'pool_mode'          Use parallel processing (parfor, default:
%                           false)
%      'size_VolChunks'     Define size for indvidual VolChunks during
%                           parallel processing. Needs to be optimised
%                           according to system specs (default: 5 GB).
%      'Save_IndImages_On'  Save indvidual images instead of tiff stack.
%      'Print_Log_On'      Enable fprintf outputs to monitor/record
%                           progress (default: false)
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
ExpShorthand = Opt.ExpShorthand;

if ~isfield(Opt,'AppShorthand')
    Opt.AppShorthand = '';
else
    Opt.Print_Str = sprintf('%s%s - ',Opt.Print_Str,Opt.AppShorthand);
end

if ~isfield(Opt,'Print_Log_On') || ~Opt.Print_Log_On
    Print_Log_On = false;
elseif Opt.Print_Log_On
    Print_Log_On = true;
end
Opt.Print_Str = sprintf('%s%s:',Opt.Print_Str,mfilename());


if isfield(Opt,'pool_mode')
    pool_mode = Opt.pool_mode;
    if pool_mode
        pool = gcp;
        n_worker_parfor = pool.NumWorkers;
    else
        n_worker_parfor = 1;
    end
else
    pool_mode = false;
    n_worker_parfor = 1;
    warning('%s pool_mode NOT DEFINED (default: false)\n',Opt.Print_Str)
end


if ~isfield(Opt,'Save_IndImages_On')
    Opt.Save_IndImages_On = true;
end

if isfield(Opt,'ImgFormat')
    ImgFormat = Opt.ImgFormat;
else
    warning('%s ImgFormat not defined (default: bmp)',mfilename())
    ImgFormat = 'bmp';
end

if isfield(Opt,'ImgName')
    ImgName = Opt.ImgName;
else
    ImgName = Opt.ExpShorthand;
end


if numel(size(V)) > 3
    Save_IndImages_RGB_On = true;
    numImg = size(V,4);
else
    Save_IndImages_RGB_On = false;
    numImg = size(V,3);
end

%% Setup folder path
path_ImgFolder = Opt.path_ImgFolder;
if Print_Log_On
    fprintf('%s Subfolderpath (%s)\n',Opt.Print_Str,path_ImgFolder)
end

% See if files already exist in location
if exist(path_ImgFolder,'dir') == 7
    filePattern = fullfile(path_ImgFolder,sprintf('*.%s',Opt.ImgFormat));
    files = dir(filePattern);
    if ~isempty(files)
        for k = 1 : length(files)
          baseFileName = files(k).name;
          fullFileName = fullfile(path_ImgFolder, baseFileName);
          delete(fullFileName);
        end
    end
end

if ~exist(path_ImgFolder,'dir')
    mkdir(path_ImgFolder)
end


%% Start imwrite export
if pool_mode
    Parfor_Prog_string = sprintf('%s (NumWorkers %.0f)',Opt.Print_Str(1:end-1),n_worker_parfor);
    Opt_parfor_progress.path_fileCom = fileparts(path_ImgFolder);
    if Print_Log_On
        fprintf('%s: NumWorkers %.0f \n',Opt.Print_Str,n_worker_parfor)
        % Progress bar for parfor
        addpath(Opt.path_FileExchange)
        parfor_progress_FrDo(ExpShorthand,Parfor_Prog_string,numImg,Opt_parfor_progress);
    end

    if Opt.Save_IndImages_On
        % Save images as individual tiffs
        parfor k = 1:numImg
            if numImg < 9999
                ImgName_Stack = sprintf('%s_%04d.%s',ImgName,k,ImgFormat);
            else
                ImgName_Stack = sprintf('%s_%08d.%s',ImgName,k,ImgFormat);
            end

            if Save_IndImages_RGB_On
                imwrite(V(:,:,:,k), fullfile(path_ImgFolder,ImgName_Stack),ImgFormat);
            else
                imwrite(V(:,:,k), fullfile(path_ImgFolder,ImgName_Stack),ImgFormat);
            end

            if Print_Log_On
                parfor_progress_FrDo(ExpShorthand,Parfor_Prog_string,-1,Opt_parfor_progress);
            end
        end        
    else
        % Initiate Stack Writing
        imwrite(V(:,:,1), fullfile(path_ImgFolder,sprintf('%s.tiff',ImgName)))
        parfor k = 2:numImg
            imwrite(V(:,:,k), fullfile(path_ImgFolder,sprintf('%s.tiff',ImgName)), 'writemode', 'append');
            if Print_Log_On
                parfor_progress_FrDo(ExpShorthand,Parfor_Prog_string,-1,Opt_parfor_progress);
            end
        end
    end
    % Close waitbar
    if Print_Log_On
        parfor_progress_FrDo(ExpShorthand,Parfor_Prog_string,0,Opt_parfor_progress);
    end
else
    if Print_Log_On
        fprintf('%s: NumWorkers 1\n',Opt.Print_Str)
    end
    
    if Opt.Save_IndImages_On
        % Save images as individual tiffs
        for k = 1:numImg
            if numImg < 9999
                ImgName_Stack = sprintf('%s_%04d.%s',ImgName,k,ImgFormat);
            else
                ImgName_Stack = sprintf('%s_%08d.%s',ImgName,k,ImgFormat);
            end

            if Save_IndImages_RGB_On
                imwrite(V(:,:,:,k), fullfile(path_ImgFolder,ImgName_Stack),ImgFormat);
            else
                imwrite(V(:,:,k), fullfile(path_ImgFolder,ImgName_Stack),ImgFormat);
            end

            if Print_Log_On
                fprintf('%s: Image %.0f/%.0f \n',Opt.Print_Str,k,numImg)
            end
        end

    else
        % Initiate Stack Writing
        imwrite(V(:,:,1), fullfile(path_ImgFolder,sprintf('%s.tiff',ImgName)))
        for k = 2:numImg
            imwrite(V(:,:,k), fullfile(path_ImgFolder,sprintf('%s.tiff',ImgName)), 'writemode', 'append');
            if Print_Log_On
                fprintf('%s: Image %.0f/%.0f \n',Opt.Print_Str,k,numImg)
            end
        end

    end
end
