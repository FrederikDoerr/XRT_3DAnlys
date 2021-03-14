function [V,R] = XRT_3DAnlys_Watershed(V,Opt)
%XRT_3DAnlys_Watershed Marker-Supported Watershed transformation.
%   [V,R] = XRT_3DAnlys_Watershed(V,Opt) computes a distance map and
%   local logic markers to guide watershed transformation the V.
%
%   OPTIONS:
%      'bwdist_method'      Distance map (default: 'chessboard')
%                           
%      'imext_Size'         Defined structuring element type for each
%                           imdilate step (default: 'cube')
%      'R_Ctrl_On'          Return image data of foreground marker 
%                           (marker_Img), distance map (D_Img) and 
%                           watershed lines (L_Img),(default: false)
%      'Print_Log_On'      Enable fprintf outputs to monitor/record
%                           progress (default: true)
%      'ExpShorthand'       Sample ID (default: '')
%      'AppShorthand'       Application ID (default: mfilename())
%
%   Example
%   -------------
%   Make a 3-D binary image containing two overlapping spheres.
%       center1 = -10;
%       center2 = -center1;
%       dist = sqrt(3*(2*center1)^2);
%       radius = dist/2 * 1.4;
%       lims = [floor(center1-1.2*radius) ceil(center2+1.2*radius)];
%       [x,y,z] = meshgrid(lims(1):lims(2));
%       bw1 = sqrt((x-center1).^2 + (y-center1).^2 + ...
%         (z-center1).^2) <= radius;
%       bw2 = sqrt((x-center2).^2 + (y-center2).^2 + ...
%         (z-center2).^2) <= radius;
%       bw = bw1 | bw2;
%       figure, isosurface(x,y,z,bw,0.5), axis equal, title('BW')
%       xlabel x, ylabel y, zlabel z
%       xlim(lims), ylim(lims), zlim(lims)
%       view(3), camlight, lighting gouraud
%
%   Watershed Segmentation
%       Opt_Watershed.Print_Log_On = true;
%       Opt_Watershed.bwdist_method = 'euclidean';
%       Opt_Watershed.imext_Size = 4;
%       Opt_Watershed.R_Ctrl_On = true;
%       [bw_S,bw_R] = XRT_3DAnlys_Watershed(bw,Opt_Watershed);
%
%   Visualise results
%       figure, isosurface(x,y,z,bw_S,0.5), axis equal, title('BW')
%       xlabel x, ylabel y, zlabel z
%       xlim(lims), ylim(lims), zlim(lims)
%       view(3), camlight, lighting gouraud
%
%
%   Author
%   --------
%   Frederik Doerr, Aug 2020 (MATLAB R2019b)
%   frederik.doerr(at)strath.ac.uk | CMAC (http://www.cmac.ac.uk/)
%   Github: https://github.com/frederik-d/XRT_3DAnlys


%% Setup
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


if ~isfield(Opt,'R_Ctrl_On')
    Opt.R_Ctrl_On = false;
end

if ~isfield(Opt,'bwdist_method')
    
    Opt.bwdist_method = 'chessboard';
    warning('%s Default bwdist_method (chessboard) \n',Opt.Print_Str)
end

if ~isfield(Opt,'imext_Size')
    Opt.imext_Size = 10;
    warning('%s Default bwdist_method (chessboard) \n',Opt.Print_Str)    
end


%% Start % % % % % % % % % % % % % % %
if Opt.Print_Log_On
    fprintf('%s INITIATED \n',Opt.Print_Str)
end



%% Distance Transformation % % % % % % % % % % % % % % % %
if Opt.Print_Log_On
    fprintf('%s (1) - bwdist (%s)\n',Opt.Print_Str,Opt.bwdist_method)
end

D = -bwdist(~V,Opt.bwdist_method);
% % % % % % % %


%% Foreground marker in local minima % % % % % % % % % % % % % % % %
if Opt.Print_Log_On
    fprintf('%s (2) - imextendedmin \n',Opt.Print_Str)
end

marker = imextendedmin(D,Opt.imext_Size);
% % % % % % % %

% % % % % % % % % % % % % % % %
if Opt.Print_Log_On
    fprintf('%s (3) - imimposemin\n',Opt.Print_Str)
end

D2 = imimposemin(D,marker);
D2(~V) = Inf;
% % % % % % % %

% % % % % % % % % % % % % % % %
if Opt.Print_Log_On
    fprintf('%s (4) - Watershed \n',Opt.Print_Str)
end

Ld = watershed(D2,26);
V(Ld == 0) = 0;



%% Covert D, Ld to int8 % % % % % % % % % % % % % % % %
% Converted to 8-bit/logical images for visualisation
if Opt.R_Ctrl_On
    min_D = min(D(:));
    max_D = 0;
    if Opt.Print_Log_On
        fprintf('%s Convert D (maxDist %.2f)\n',Opt.Print_Str,abs(min_D))
    end

    D_Img = (D + (0-min_D))./(max_D-min_D);
    D_Img = uint8(D_Img * 255);


    % Ld
    if Opt.Print_Log_On
        fprintf('%s Convert Ld \n',Opt.Print_Str)
    end
   
    Ld_Img = zeros(size(Ld),'logical');
    Ld_Img(Ld ~= 0) = 1;

    
    R.marker_Img = marker;
    R.D_Img = D_Img;
    R.Ld_Img = Ld_Img;
else
    R = struct();
end


% % % % % % % % % % % % % % % %
if Opt.Print_Log_On
    fprintf('%s COMPLETED \n',Opt.Print_Str)
end
