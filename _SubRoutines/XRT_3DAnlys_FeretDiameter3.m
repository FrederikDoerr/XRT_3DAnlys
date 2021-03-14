function[max_dist,XRT_FeretMax3_Results] = XRT_3DAnlys_FeretDiameter3(V,Opt)
%XRT_3DAnlys_FeretDiameter3 to analyse maximum Feret Diameter.
%   [R_stats,R_stats3,R] = XRT_3DAnlys_FeretDiameter3(V,Opt) function to  
%   claculate maximum Feret Diameter from a 3D VoxelList.
%
%   OPTIONS:
%      'Print_Log_On'      Enable fprintf outputs to monitor/record
%                           progress (default: false)
%      'ExpShorthand'       Sample ID (default: '')
%      'AppShorthand'       Application ID (default: mfilename())
% 
%   Links
%       https://blogs.mathworks.com/steve/2017/09/29/feret-diameter-introduction/
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


if numel(size(V)) > 2 % 3D or 2D Image
    try
        % Find XYZ LIst for all Voxels of V
        stats3 = regionprops3(V,'VoxelList');
        V_XYZ_reProp = stats3.VoxelList{1};

        k = convhull(V_XYZ_reProp,'Simplify',true);
        hull_voxel = V_XYZ_reProp(k,:);


        dx = hull_voxel(:,1) - hull_voxel(:,1)';
        dy = hull_voxel(:,2) - hull_voxel(:,2)';
        dz = hull_voxel(:,3) - hull_voxel(:,3)';

        pairwise_dist = sqrt(abs(dx).^2+abs(dy).^2 + abs(dz).^2);
        [max_dist,idx_dist] = max(pairwise_dist(:));
        
        XRT_FeretMax3_Results.CH_voxel = hull_voxel;
        
        [dx_idx1,dx_idx2] = ind2sub(size(dx),idx_dist);
        XRT_FeretMax3_Results.idx1_x = hull_voxel(dx_idx1,1);
        XRT_FeretMax3_Results.idx2_x = hull_voxel(dx_idx2,1)';
        XRT_FeretMax3_Results.dx = dx(idx_dist);
        
        [dy_idx1,dy_idx2] = ind2sub(size(dy),idx_dist);
        XRT_FeretMax3_Results.idx1_y = hull_voxel(dy_idx1,2);
        XRT_FeretMax3_Results.idx2_y = hull_voxel(dy_idx2,2)';
        XRT_FeretMax3_Results.dy = dy(idx_dist);
        
        [dz_idx1,dz_idx2] = ind2sub(size(dz),idx_dist);
        XRT_FeretMax3_Results.idx1_z = hull_voxel(dz_idx1,3);
        XRT_FeretMax3_Results.idx2_z = hull_voxel(dz_idx2,3)';
        XRT_FeretMax3_Results.dz = dz(idx_dist);
        
    catch ME %ME is an MException struct
        max_dist = nan;
        XRT_FeretMax3_Results.CH_voxel = nan;
        XRT_FeretMax3_Results.idx_dist = nan;
        XRT_FeretMax3_Results.idx1_x = nan;
        XRT_FeretMax3_Results.idx2_x = nan;
        XRT_FeretMax3_Results.dx = nan;
        XRT_FeretMax3_Results.idx1_y = nan;
        XRT_FeretMax3_Results.idx2_y = nan;
        XRT_FeretMax3_Results.dy = nan;
        XRT_FeretMax3_Results.idx1_z = nan;
        XRT_FeretMax3_Results.idx2_z = nan;
        XRT_FeretMax3_Results.dz = nan;

        fprintf('%s - %s: %s ERROR, %s (%s)\n',Opt.ExpShorthand, ...
            Opt.AppShorthand,mfilename(),ME.message,ME.identifier)

    end
else

    try
        % Find XYZ LIst for all Voxels of V
        stats = regionprops(V,'PixelList');

        V_XYZ_reProp = stats.PixelList;

        k = convhull(V_XYZ_reProp,'Simplify',true);
        hull_voxel = V_XYZ_reProp(k,:);

        dx = hull_voxel(:,1) - hull_voxel(:,1)';
        dy = hull_voxel(:,2) - hull_voxel(:,2)';
        
        % pairwise_dist = hypot(dx,dy,dz);
        pairwise_dist = sqrt(abs(dx).^2+abs(dy).^2);
        
        [max_dist,idx_dist] = max(pairwise_dist(:));
        
        
        XRT_FeretMax3_Results.CH_voxel = hull_voxel;
        
        [dx_idx1,dx_idx2] = ind2sub(size(dx),idx_dist);
        XRT_FeretMax3_Results.idx1_x = hull_voxel(dx_idx1,1);
        XRT_FeretMax3_Results.idx2_x = hull_voxel(dx_idx2,1)';
        XRT_FeretMax3_Results.dx = dx(idx_dist);
        
        [dy_idx1,dy_idx2] = ind2sub(size(dy),idx_dist);
        XRT_FeretMax3_Results.idx1_y = hull_voxel(dy_idx1,2);
        XRT_FeretMax3_Results.idx2_y = hull_voxel(dy_idx2,2)';
        XRT_FeretMax3_Results.dy = dy(idx_dist);
        
        XRT_FeretMax3_Results.idx1_z = nan;
        XRT_FeretMax3_Results.idx2_z = nan;
        XRT_FeretMax3_Results.dz = nan;

    catch ME %ME is an MException struct
        
        max_dist = nan;
        XRT_FeretMax3_Results.CH_voxel = nan;
        XRT_FeretMax3_Results.idx_dist = nan;
        XRT_FeretMax3_Results.idx1_x = nan;
        XRT_FeretMax3_Results.idx2_x = nan;
        XRT_FeretMax3_Results.dx = nan;
        XRT_FeretMax3_Results.idx1_y = nan;
        XRT_FeretMax3_Results.idx2_y = nan;
        XRT_FeretMax3_Results.dy = nan;
        XRT_FeretMax3_Results.idx1_z = nan;
        XRT_FeretMax3_Results.idx2_z = nan;
        XRT_FeretMax3_Results.dz = nan;

        fprintf('%s - %s: %s ERROR, %s (%s)\n',Opt.ExpShorthand, ...
            Opt.AppShorthand,mfilename(),ME.message,ME.identifier)
    end 
end
            
