function [V] = XRT_3DAnlys_ImgStack_Load(Opt)
%Load image data from repository.
%   [V] = XRT_3DAnlys_ImgStack_Load(Opt) function to ease the routine 
%   loading of image data from folders.
%   Potential to use parallel processing (Opt.pool_mode = true, see parfor)
%
%   OPTIONS:
%      'path_ImgFolder'     Image folder path.
%      'ImgFormat'          Image file format.
%      'pool_mode'          Use parallel processing (parfor, default:
%                           false)
%      'size_VolChunks'     Define size for indvidual VolChunks during
%                           parallel processing. Needs to be optimised
%                           according to system specs (default: 5 GB).
%      'Print_Ctrl_On'      Enable fprintf outputs to monitor/record
%                           progress (default: false)
%      'ExpShorthand'       Sample ID (default: '')
%      'AppShorthand'       Application ID (default: mfilename())
%
%   Example
%   -------------
%   % see Demo_XRT_3DAnlys_Watershed.m
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

if ~isfield(Opt,'Print_Ctrl_On')
    Opt.Print_Ctrl_On = false;
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

if ~isfield(Opt,'size_VolChunks')
    Opt.size_VolChunks = 5;
end


filelist_img = dir(fullfile(Opt.path_ImgFolder,sprintf('*.%s',Opt.ImgFormat)));

% Create empty ImageV_binary
idx_c = round(length(filelist_img)/2);
I_init = imread(fullfile(filelist_img(idx_c).folder,filelist_img(idx_c).name));

[mImage,nImage] = size(I_init);
numImg = length(filelist_img);


if ~islogical(I_init)
    V = zeros(mImage,nImage,numImg,'uint8');
else
	V = zeros(mImage,nImage,numImg,'logical');
end


% Determine the number of VolChunks
ImageStack_Size = filelist_img(idx_c).bytes * numImg * 1e-9;
num_VolChunks = floor(ImageStack_Size / Opt.size_VolChunks);
if num_VolChunks < 1
    num_VolChunks = 1;
end
num_VolChunks_nImg = floor(numImg/num_VolChunks);
VolChunk_ImgID_Array = 1:num_VolChunks_nImg-1:num_VolChunks*num_VolChunks_nImg;
if (numImg - num_VolChunks * num_VolChunks_nImg) > 0
    VolChunk_ImgID_Array = [VolChunk_ImgID_Array,numImg];
end



if num_VolChunks > 1
    fprintf('%s num_VolChunks %.0f\n',Opt.Print_Str,num_VolChunks)
end

for iter = 1:num_VolChunks
    
    VolChunk_Start = VolChunk_ImgID_Array(iter);
    VolChunk_End = VolChunk_ImgID_Array(iter+1);
    numImg_VolChunks = VolChunk_End - VolChunk_Start;

    if ~islogical(I_init)
        ImageV_binary = zeros(mImage,nImage,numImg_VolChunks,'uint8');
    else
        ImageV_binary = zeros(mImage,nImage,numImg_VolChunks,'logical');
    end
    
    filelist_img_VolChunk = filelist_img(VolChunk_Start:VolChunk_End);
    
    if pool_mode
        numImg_VolChunk = length(filelist_img_VolChunk);
        parfor k = 1:numImg_VolChunk
            I = imread(fullfile(filelist_img_VolChunk(k).folder,filelist_img_VolChunk(k).name));
            ImageV_binary(:,:,k) = I;
        end
    else
        numImg_VolChunk = length(filelist_img_VolChunk);
        for k = 1:numImg_VolChunk
            I = imread(fullfile(filelist_img_VolChunk(k).folder,filelist_img_VolChunk(k).name));
            ImageV_binary(:,:,k) = I;
            if Opt.Print_Ctrl_On
                fprintf('%s: Loaded %.0f / %.0f\n',Opt.Print_Str,k,numImg_VolChunk)
            end
        end
    end

    V(:,:,VolChunk_Start:VolChunk_End) = ImageV_binary;

    if Opt.Print_Ctrl_On
        fprintf('%s %.0f / %.0f (NumWorkers %.0f)\n',Opt.Print_Str,iter,num_VolChunks,n_worker_parfor)
    end
end


clearvars ImageV_binary ImageV_binary_ROI
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 