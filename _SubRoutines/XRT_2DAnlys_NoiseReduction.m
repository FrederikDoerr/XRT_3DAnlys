function[bw] = XRT_2DAnlys_NoiseReduction(bw,Opt)
%2D image noise reduction.
%   [bw] = XRT_2DAnlys_NoiseReduction(bw,Opt) function to ease the routine
%   application of 2D noise reduction methods during image processing.
%
%   OPTIONS:
%      'ClearBorder'        Exceute imclearborder
%      'WhiteNoise_Size'    Despeckle white (foreground) noise
%      'Close_Size'         morphological operation: imclose
%      'Open_Size'          morphological operation: imopen 
%      'BlackNoise_Size'    Despeckle black (background) noise
% 
%
%   Author
%   --------
%   Frederik Doerr, Aug 2020 (MATLAB R2019b)
%   frederik.doerr(at)strath.ac.uk | CMAC (http://www.cmac.ac.uk/)
%   Github: https://github.com/frederik-d/XRT_3DAnlys


% Remove any blobs on the border of the image
if isfield(Opt,'ClearBorder')
    bw = imclearborder(bw);
end

% Despeckle all foreground objects containing fewer than [WhiteNoise_Size] pixels, Despeckle white noise
if isfield(Opt,'WhiteNoise_Size')
    bw = bwareaopen(bw,Opt.WhiteNoise_Size);
end

% imclose with radius = ThresSize_close
if isfield(Opt,'Close_Size')
    se = strel('disk',Opt.Close_Size);
    bw = imclose(bw,se);    
end

% imopen with radius = ThresSize_open
if isfield(Opt,'Open_Size')
    se = strel('disk',Opt.Open_Size);
    bw = imopen(bw,se);
end

% Despeckle all background objects containing fewer than [ThresSize_b] pixels, Despeckle black noise
if isfield(Opt,'BlackNoise_Size')
    bw = imcomplement(bwareaopen(imcomplement(bw),Opt.BlackNoise_Size));
end

if isfield(Opt,'SizeSweep') && Opt.SizeSweep
    if isfield(Opt,'SizeSweep_numObj')
        % Keep n largest areas
        bw = bwareafilt(bw,Opt.SizeSweep_numObj);
    else
        % Keep largest areas
        bw = bwareafilt(bw,1);
    end
end

