function XRT_3DAnlys_ImgExport(V,Opt)
%XRT_3DAnlys_ImgExport Export 2D cross-section images.
%   XRT_3DAnlys_ImgExport(V,Opt) export 2D image cross-sections normal to 
%   each axis of the current image space in V.
%
%   Dependencies on FileExchange packages/functions:
%       > XRT_2DAnlys_AddScalebar
% 
%   OPTIONS:
%      'path_ImgFolder'     Image folder path.
%      'ImgFormat'          Image file format (default: 'bmp').
%      'ImgName'            Image file name.
%      'ImagePixelSize'     Original imagepixelsize (um/px)
%      'ImagePixelSize_rS' 	Scaled imagepixelsize (um/px)	
%      'DataResizeFactor'   User-defined image data resize factor
% 
%   Author
%   --------
%   Frederik Doerr, Aug 2020 (MATLAB R2019b)
%   frederik.doerr(at)strath.ac.uk | CMAC (http://www.cmac.ac.uk/)
%   Github: https://github.com/frederik-d/XRT_3DAnlys

if ~isfield(Opt,'ImagePixelSize_rS')
    Opt.ImagePixelSize_rS = Opt.ImagePixelSize*Opt.DataResizeFactor;
end


if ~exist(Opt.path_ImgFolder,'dir')
    mkdir(Opt.path_ImgFolder);
end


[numRow,numCol,numImg] = size(V);
I = V(:,:,round(numImg/2));

try
    I = XRT_2DAnlys_AddScalebar(I,Opt.ImagePixelSize_rS,Opt);
end
imwrite(I, fullfile(Opt.path_ImgFolder,sprintf('%s_XY.bmp',Opt.ImgName)),'bmp')

I = V(round(numRow/2),:,:);I = permute(I,[3,2,1]);

try
    I = XRT_2DAnlys_AddScalebar(I,Opt.ImagePixelSize_rS,Opt);
end
imwrite(I, fullfile(Opt.path_ImgFolder,sprintf('%s_XZ.bmp',Opt.ImgName)),'bmp')

I = V(:,round(numCol/2),:);I = permute(I,[3,1,2]);

try
    I = XRT_2DAnlys_AddScalebar(I,Opt.ImagePixelSize_rS,Opt);
end
imwrite(I, fullfile(Opt.path_ImgFolder,sprintf('%s_YZ.bmp',Opt.ImgName)),'bmp')
