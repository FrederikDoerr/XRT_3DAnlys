function[V_ConvexHull3] = XRT_3DAnlys_ConvexHull3(V,Opt)
%XRT_3DAnlys_ConvexHull3 Convex Hull of 3-D volumetric images.
%   V_ConvexHull3 = XRT_3DAnlys_ConvexHull3(V,XRT_3DAnlys_ConvexHull3_Opts) 
%   creates a Convex Hull of V.V must be a binary 2-D or 3-D volumetric image.
%   Opt is a structure with properties (fields):
%       'filesPath_FileExchange' path to FileExchange packages.
% 
%   Example 1
%   ---------
%   % Create a 3D Convex Hull of objects in a 3-D volumetric image
%
%         % Create a binary image with two spheres
%         [x, y, z] = meshgrid(1:50, 1:50, 1:50);
%         bw1 = sqrt((x-10).^2 + (y-15).^2 + (z-35).^2) < 5;
%         bw2 = sqrt((x-20).^2 + (y-30).^2 + (z-15).^2) < 10;
%         V = bw1 | bw2;
%         [V_ConvexHull3] = XRT_3DAnlys_ConvexHull3(V, ...
%                               XRT_3DAnlys_ConvexHull3_Opts);
%         
%   See also REGIONPROPS3(BW,'ConvexImage').     
%         
%   % % References (open access):
%   Doerr, F. J. S., Florence, A. J. (2020)
%   A micro-XRT image analysis and machine learning methodology for the characterisation of multi-particulate capsule formulations. 
%   International Journal of Pharmaceutics: X. 
%   https://doi.org/10.1016/j.ijpx.2020.100041
% 
%   Doerr, F. J. S., Oswald, I. D. H., & Florence, A. J. (2018)
%   Quantitative investigation of particle formation of a model pharmaceutical formulation using single droplet evaporation experiments and X-ray tomography. 
%   Advanced Powder Technology, 29(12), 2996–3006. 
%   https://doi.org/10.1016/j.apt.2018.09.027
%
%   % % Related Links
%   https://uk.mathworks.com/help/matlab/ref/convhulln.html
%   https://uk.mathworks.com/matlabcentral/fileexchange/10226-inhull



% Find XYZ LIst for all Voxels of V
stats3 = regionprops3(V,'VoxelList');
numObj = size(stats3,1);
if numObj > 1
    warning('%s: More than one object detected (size(stats3,1) > 1)!',mfilename())
    V_XYZ = [];
    for k = 1:numObj
        V_XYZ = [V_XYZ;stats3.VoxelList{k}];
    end
else
    V_XYZ = stats3.VoxelList{1};
end

k = convhulln(V_XYZ);

% figure; trisurf(k,V_XYZ(:,1),V_XYZ(:,2),V_XYZ(:,3),'FaceColor','cyan')

% Check each voxel in Image Space inside K-3DConvexHull
addpath(fullfile(Opt.filesPath_FileExchange,'inhull'))
V_ConvexHull3 = ones(size(V),'logical');
V_ConvexHull3(V) = 0;
stats3_3 = regionprops3(V_ConvexHull3,'VoxelList');
V_ConvexHull3_XYZ = stats3_3.VoxelList{1};
in = inhull(V_ConvexHull3_XYZ,V_XYZ,k);  

% Define XYZ ConvexImage
V_ConvexHull3_XYZ(~in,:)=[];

% Find original Dimensions and add as corner to XYZ
maxSize_x = size(V,2);
maxSize_y = size(V,1);
maxSize_z = size(V,3);
V_ConvexHull3_XYZ = [V_ConvexHull3_XYZ;maxSize_x,maxSize_y,maxSize_z];

% Transfer XYZ to image space
V_ConvexHull3 = XRT_3DAnlys_XYZ2Image(V_ConvexHull3_XYZ);
V_ConvexHull3(maxSize_y,maxSize_x,maxSize_z) = 0;
         
V_ConvexHull3(V) = 1;


