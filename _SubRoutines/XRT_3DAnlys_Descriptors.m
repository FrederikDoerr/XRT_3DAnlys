function[R_stats,R_stats3,R] = XRT_3DAnlys_Descriptors(V,Opt)
%XRT_3DAnlys_Descriptors to analyse 3D binary image data.
%   [R_stats,R_stats3,R] = XRT_3DAnlys_Descriptors(V,Opt) function to ease 
%   the routine extraction of object descriptors/features. 
%
%   OPTIONS:
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

if ~isfield(Opt,'AppShorthand')
    Opt.AppShorthand = '';
else
    Opt.Print_Str = sprintf('%s%s - ',Opt.Print_Str,Opt.AppShorthand);
end

if ~isfield(Opt,'Print_Log_On')
    Opt.Print_Log_On = false;
end
Opt.Print_Str = sprintf('%s%s:',Opt.Print_Str,mfilename());


DataResizeFactor = Opt.DataResizeFactor;
ImagePixelSize_rS = Opt.ImagePixelSize*DataResizeFactor;

if ~isfield(Opt,'numImg_min') || isnan(Opt.numImg_min)
    Opt.numImg_min = 0;
end

if ~isfield(Opt,'numRow_min') || isnan(Opt.numRow_min)
    Opt.numRow_min = 0;
end

if ~isfield(Opt,'numCol_min') || isnan(Opt.numCol_min)
    Opt.numCol_min = 0;
end


R = struct();

[numRows,numCols,numImg] = size(V);
numDim_minSize = Opt.numImg_min*Opt.numRow_min*Opt.numCol_min;

if numImg < Opt.numImg_min || ...
        numRows < Opt.numRow_min || ...
        numCols < Opt.numCol_min

    warning('%s numDim_minSize <= %.0f, skipped\n',Opt.Print_Str,numDim_minSize)
    radius_spherefit = nan;
    radii_elpsfit = [nan,nan,nan];
    chi2_elpsfit = nan;
    Area = nan;
    Surface = nan;
    AspectRatio = nan;
    BoundingBox = [nan,nan];

    Volume_3_Max = nan;
    BoundingBox_3_Max = nan(2,1);
    SurfaceArea_3_Max  = nan;
    Solidity_CH_3_Max  = nan;
    Extent_3_Max = nan;
    Centroid_3_Max = nan;
    ConvexVolume_3_Max = nan;
    PrincipalAxisLength_3_Max = nan(3,1);

    EquivDiameter_3_Max  = nan;
    Orientation_3_Max  = nan(1,3);
    EigenVectors_3_Max  = nan(3,3);
    EigenValues_3_Max  = nan(3,1);
else

    % Get k_v_List (x,y,z) data from boundary region of image stack
    k_v_List = [];
    Surface_Count = 0;
    V_S = zeros(size(V),'logical');
    for v_iter = 1:numImg
        bw = V(:,:,v_iter);
        bw_S = bwboundaries(bw);
        if ~isempty(bw_S)
            [row_k_v,~] = size(bw_S{1});
            k_v = [bw_S{1},linspace(v_iter,v_iter,row_k_v)'];
            k_v_List = [k_v_List;k_v];
            Surface_Count = Surface_Count + length(k_v);
        end
%             V_S(:,:,v_iter) = bw_S;
        Surface = Surface_Count*ImagePixelSize_rS^2;
    end

    addpath(fullfile(Opt.filesPath_FileExchange,'ShapeFitting_MinBoundSuite'))
%         [center_spherefit,radius_spherefit] = minboundsphere(k_v_List);
    try 
        [~,radius_spherefit] = minboundsphere(k_v_List);
        radius_spherefit = radius_spherefit*ImagePixelSize_rS;
    catch ME %ME is an MException struct   
        radius_spherefit = nan;
        warning('%s minboundsphere ERROR, %s (%s)\n',Opt.Print_Str,ME.message,ME.identifier)
    end


    addpath(fullfile(Opt.filesPath_FileExchange,'ellipsoid_fit\ellipsoid_fit'))
    try 
        [~, radii_elpsfit, ~, ~, chi2_elpsfit] = ellipsoid_fit(k_v_List);
        chi2_elpsfit = chi2_elpsfit*ImagePixelSize_rS;
        radii_elpsfit = radii_elpsfit.*ImagePixelSize_rS;
        radii_elpsfit_Max = max([radii_elpsfit(1),radii_elpsfit(2),radii_elpsfit(3)]);
        radii_elpsfit_Min = min([radii_elpsfit(1),radii_elpsfit(2),radii_elpsfit(3)]); 
        AspectRatio = radii_elpsfit_Max / radii_elpsfit_Min;
    catch ME %ME is an MException struct   
        chi2_elpsfit = nan;
        radii_elpsfit = [nan,nan,nan];
        AspectRatio = nan;

        warning('%s ellipsoid_fit ERROR, %s (%s)\n',Opt.Print_Str,ME.message,ME.identifier)
    end

    % Output:
    % * center    -  ellispoid or other conic center coordinates [xc; yc; zc]
    % * radii     -  ellipsoid or other conic radii [a; b; c]
    % * evecs     -  the radii directions as columns of the 3x3 matrix
    % * v         -  the 10 parameters describing the ellipsoid / conic algebraically: 
    %                Ax^2 + By^2 + Cz^2 + 2Dxy + 2Exz + 2Fyz + 2Gx + 2Hy + 2Iz + J = 0
    % * chi2      -  residual sum of squared errors (chi^2), this chi2 is in the 
    %                coordinate frame in which the ellipsoid is a unit sphere.

    stats = regionprops(V,'Area','BoundingBox','Centroid');

    Area = [stats.Area]*ImagePixelSize_rS^3;
    Centroid = [stats.Centroid];

    [~,idx] = max([stats.Area].');
    bb = stats(idx).BoundingBox; %[ul_corner width], [x y z x_width y_width z_width]
    x_width = bb(4);
    y_width = bb(5);
    z_width = bb(6);
    bb_MaxLength = max([x_width,y_width,z_width]);
    bb_MinLength = min([x_width,y_width,z_width]);

%         rectangle('Position',[bb(1),bb(2),bb(3),bb(4)],...
% 'EdgeColor','r','LineWidth',2 ) %[x y w h]

    BoundingBox = [bb_MaxLength,bb_MinLength].*ImagePixelSize_rS;


%         % regionprops3 (Introduced in R2017b, https://uk.mathworks.com/help/images/ref/regionprops3.html)
%         'BoundingBox': 	Smallest cuboid containing the region, returned as a 1-by-6 vector of the form [ulf_x ulf_y ulf_z width_x width_y width_z]. ulf_x, ulf_y, and ulf_z specify the upper-left front corner of the cuboid. width_x, width_y, and width_z specify the width of the cuboid along each dimension.
%         'Centroid': Center of mass of the region, returned as a 1-by-3 vector of the form [centroid_x centroid_y and centroid_z]. The first element, centroid_x, is the horizontal coordinate (or x-coordinate) of the center of mass. The second element, centroid_y, is the vertical coordinate (or y-coordinate). The third element, centroid_z, is the planar coordinate (or z-coordinate).
%         'ConvexVolume':	Number of voxels in 'ConvexImage', returned as a scalar.
%         'PrincipalAxesLength':	Length (in voxels) of the major axes of the ellipsoid that have the same normalized second central moments as the region, returned as 1-by-3 vector. regionprops3 sorts the values from highest to lowest.
%         'Solidity':	Proportion of the voxels in the convex hull that are also in the region, returned as a scalar. Computed as Volume/ConvexVolume.
% 
%       'Volume', 'Centroid', 'BoundingBox', 'SubarrayIdx', 'Image', 'EquivDiameter', 'Extent', 'VoxelIdxList', 'VoxelList', 'PrincipalAxisLength',
%       'Orientation', 'EigenVectors', 'EigenValues', 'ConvexHull', 'ConvexImage', 'ConvexVolume', 'Solidity', 'SurfaceArea', 'VoxelValues',
%       'WeightedCentroid', 'MeanIntensity', 'MinIntensity', 'MaxIntensity'
    try
        stats3 = regionprops3(V,'Volume','BoundingBox','SurfaceArea','Centroid','PrincipalAxisLength','EquivDiameter','Solidity','Extent','Orientation','EigenVectors','EigenValues');
        [~,idx3] = max([stats3.Volume].');

        Volume_3_Max = stats3.Volume(idx3,:) * ImagePixelSize_rS^3;
        BoundingBox_3_Max  = stats3.BoundingBox(idx3,:).* ImagePixelSize_rS;
        SurfaceArea_3_Max  = stats3.SurfaceArea(idx3,:) * ImagePixelSize_rS^2;
        Centroid_3_Max = stats3.Centroid(idx3,:);
        ConvexVolume_3_Max = nan;
        PrincipalAxisLength_3_Max = stats3.PrincipalAxisLength(idx3,:).* ImagePixelSize_rS;

        Solidity_CH_3_Max  = stats3.Solidity(idx3,:);
        Extent_3_Max  = stats3.Extent(idx3,:);
        EquivDiameter_3_Max  = stats3.EquivDiameter(idx3,:).* ImagePixelSize_rS;
        Orientation_3_Max  = stats3.Orientation(idx3,:);
        EigenVectors_3_Max  = stats3.EigenVectors{idx3,:};
        EigenValues_3_Max  = stats3.EigenValues{idx3,:};

    catch ME %ME is an MException struct            
        Volume_3_Max = nan;
        BoundingBox_3_Max = nan(2,1);
        SurfaceArea_3_Max  = nan;
        Solidity_CH_3_Max  = nan;
        Extent_3_Max = nan;
        Centroid_3_Max = nan;
        ConvexVolume_3_Max = nan;
        PrincipalAxisLength_3_Max = nan(3,1);

        EquivDiameter_3_Max  = nan;
        Orientation_3_Max  = nan(1,3);
        EigenVectors_3_Max  = nan(3,3);
        EigenValues_3_Max  = nan(3,1);

        warning('%s regionprops3 ERROR, %s (%s)\n',Opt.Print_Str,ME.message,ME.identifier)
    end
end

R_stats.Area = Area;
R_stats.Surface = Surface;
R_stats.BoundingBox = BoundingBox;
R_stats.radius_spherefit = radius_spherefit;
R_stats.radii_elpsfit = radii_elpsfit;
R_stats.chi2_elpsfit = chi2_elpsfit;
R_stats.AspectRatio = AspectRatio;
R_stats.Centroid = Centroid;


R_stats3.Volume_3_Max = Volume_3_Max;
R_stats3.BoundingBox_3_Max = BoundingBox_3_Max;
R_stats3.SurfaceArea_3_Max = SurfaceArea_3_Max;
R_stats3.Solidity_CH_3_Max = Solidity_CH_3_Max;
R_stats3.Extent_3_Max = Extent_3_Max;
R_stats3.Centroid_3_Max = Centroid_3_Max;
R_stats3.ConvexVolume_3_Max = ConvexVolume_3_Max;
R_stats3.PrincipalAxisLength_3_Max = PrincipalAxisLength_3_Max;
R_stats3.EquivDiameter_3_Max = EquivDiameter_3_Max;
R_stats3.Orientation_3_Max = Orientation_3_Max;
R_stats3.EigenVectors_3_Max = EigenVectors_3_Max;
R_stats3.EigenValues_3_Max = EigenValues_3_Max;




