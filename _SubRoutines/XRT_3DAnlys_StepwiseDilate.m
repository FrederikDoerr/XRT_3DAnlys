function[V,R] = XRT_3DAnlys_StepwiseDilate(V,V_ROI,Opt)
%XRT_3DAnlys_StepwiseDilate Stepwise dilation operation on binary volume.
%   [V,R] = XRT_3DAnlys_StepwiseDilate(V,V_ROI,Opt) applies a stepwise
%   dilation to the binary volume V  by adding pixels to the exterior 
%   of objects. V_ROI is a binary 3D mask to guide the region growing
%   process. Select options to control dilation until doing so would  
%   result in previously unconnected objects being 26-connected.
% 
%   R: structure with process meta-data
%
%   OPTIONS:
%      'numDilationSteps'   Define number of dilation steps
%      'method'             '3d_stepwisedilate': Method for stepwise 
%                           imdilate in V_ROI logic mask
%                           
%                           '3d_bwmorph_thicken': Method for stepwise 
%                           imdilate in V_ROI logic mask without merging 
%                           previously isolated volumes
%                           
%      'strel_type'         Defined structuring element type for each
%                           imdilate step (default: 'cube')
%      'strel_size'         Defined structuring element size for each
%                           imdilate step (default: 3)
%      'Print_Log_On'      Enable fprintf outputs to monitor/record
%                           progress (default: true)
%      'ExpShorthand'       Sample ID (default: '')
%      'AppShorthand'       Application ID (default: mfilename())
%
%   Example
%   --------
%   % This example carries out the '3d_bwmorph_thicken' operation
%   % on a 3D binary volume with two spheres.
%     center1 = 40;
%     center2 = center1 + 40;
%     dist = sqrt(1*(2*center1)^2);
%     radius = dist/4 * 1.4;
%     lims = [0 ceil(center2+1.8*radius)];
%     [x,y,z] = meshgrid(lims(1):lims(2));
%     bw1 = sqrt((x-center1).^2 + (y-center1).^2 + ...
%        (z-center1).^2) <= radius;
%     bw2 = sqrt((x-center2).^2 + (y-center2).^2 + ...
%        (z-center2).^2) <= radius;
% 
%     V = bw1 | bw2;
% 
%     figure, isosurface(x,y,z,V,0.5), axis equal, title('V')
%     xlabel x, ylabel y, zlabel z
%     xlim(lims), ylim(lims), zlim(lims)
%     view(3), camlight, lighting gouraud
% 
%     % Define full image V_ROI
%     V_ROI = ones(size(V),'logical');
% 
%
%     % 3d_bwmorph_thicken
%     Opt.numDilationSteps = 10;
%     Opt.method = '3d_bwmorph_thicken';
%     Opt.ExpShorthand = 'Spheres (3d_bwmorph_thicken)';
%     Opt.Print_Log_On = true;
%     [V_thi,~] = XRT_3DAnlys_StepwiseDilate(V,V_ROI,Opt);
% 
%     figure, isosurface(x,y,z,V_thi,0.5), axis equal, title('V_thi')
%     xlabel x, ylabel y, zlabel z
%     xlim(lims), ylim(lims), zlim(lims)
%     view(3), camlight, lighting gouraud
%
%
%   Author
%   --------
%   Frederik Doerr, Aug 2020 (MATLAB R2019b)
%   frederik.doerr(at)strath.ac.uk | CMAC (http://www.cmac.ac.uk/)
%   Github: https://github.com/frederik-d/XRT_3DAnlys
%
%
%   See also IMERODE, IMDILATE, IMBOTHAT, IMTOPHAT, IMCLOSE, IMOPEN,
%   BWMORPH, BWMORPH3, BWSKEL.

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


if ~isfield(Opt,'strel_type') 
    strel_type = 'cube';
    if ~isfield(Opt,'strel_size') 
        Opt.strel_size = 3; % strel_size = 1 does not work
    end
else
    strel_type = Opt.strel_type;
    if ~isfield(Opt,'strel_size') 
        Opt.strel_size = 3; % strel_size = 1 does not work
    end
end

if ~isfield(Opt,'Print_Log_On') 
    Opt.Print_Log_On = true;
end

% Pass options to return structure
R = Opt;

switch lower(Opt.method)
    case {'3d_stepwisedilate'}
        % Stepwise dialte V within the V_ROI
        for k = 1: Opt.numDilationSteps
            se = strel(strel_type,Opt.strel_size);
            V = imdilate(V,se);

            % Re-evaluate Dilation Step with V_ROI
            V(~V_ROI) = 0;

            if Opt.Print_Log_On
                fprintf('%s %s - Step %.0f/%.0f \n',Opt.Print_Str,Opt.method,k, Opt.numDilationSteps)
            end
        end

    case {'3d_bwmorph_thicken'}
        % Thicken 3D V without fusing individual volume elements (gap detetcion)

         [L,numLabels] = bwlabeln(V);

        for k = 1: Opt.numDilationSteps
            V_Exp = zeros(size(V),'double'); 
            for i1 = 1:numLabels
                se = strel(strel_type,Opt.strel_size);
                V_obj = imdilate(L == i1,se);
                L(V_obj) = i1;
                V_Exp = V_Exp+V_obj;
            end

            % Find overlap between Region Labels (Gap_Detect)
            Gap_Detect = V_Exp > 1;
            Gap_Detect = imdilate(Gap_Detect,strel('cube',3));

            L (Gap_Detect) = 0;

            % Re-evaluate Dilation Step with V_ROI
            L(~V_ROI) = 0;

            if Opt.Print_Log_On
                fprintf('%s %s - Step %.0f/%.0f \n',Opt.Print_Str,Opt.method,k,Opt.numDilationSteps)   
            end
        end
        V = L > 0;
        R.Gap_Detect = Gap_Detect;

    otherwise
        warning('%s %s NOT defined!\n',Opt.Print_Str,Opt.method)   
end

