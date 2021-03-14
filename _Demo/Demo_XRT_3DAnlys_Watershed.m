%------------------------------------------------------------------------------------------------
% Code written by Frederik Doerr, Aug 2020 (MATLAB R2019b)
% Application: Demo for XRT_3DAnlys_Watershed
% 
% https://github.com/frederik-d
% Contact: frederik.doerr(at)strath.ac.uk | CMAC (http://www.cmac.ac.uk/)

clear all
close all
clc

%% Setup
set(0,'DefaultFigureVisible','on');

% Folder locations
path = matlab.desktop.editor.getActiveFilename;
[Opt.path_Mat_Demo,name,ext] = fileparts(path);

Opt.path_MAIN = fileparts(Opt.path_Mat_Demo);
Opt.path_Mat_SR = fullfile(Opt.path_MAIN,'_SubRoutines');
addpath(genpath(Opt.path_MAIN));

Opt.path_Mat_Demo_Data = fullfile(Opt.path_MAIN,'_Demo');

%% Demo 1: Two spheres
fprintf('\n - - -\nDemo 1: Two spheres\n - - -\n')

% Make a 3-D binary image containing two overlapping spheres.
center1 = -10;
center2 = -center1;
dist = sqrt(3*(2*center1)^2);
radius = dist/2 * 1.4;
lims = [floor(center1-1.2*radius) ceil(center2+1.2*radius)];
[x,y,z] = meshgrid(lims(1):lims(2));
bw1 = sqrt((x-center1).^2 + (y-center1).^2 + ...
  (z-center1).^2) <= radius;
bw2 = sqrt((x-center2).^2 + (y-center2).^2 + ...
  (z-center2).^2) <= radius;
bw = bw1 | bw2;
figure, isosurface(x,y,z,bw,0.5), axis equal, title('BW')
xlabel x, ylabel y, zlabel z
xlim(lims), ylim(lims), zlim(lims)
view(3), camlight, lighting gouraud

% Watershed Segmentation
Opt_Watershed.Print_Ctrl_On = false;
Opt_Watershed.bwdist_method = 'euclidean';
Opt_Watershed.imext_Size = 4;
Opt_Watershed.R_Ctrl_On = true;
[bw_S,bw_R] = XRT_3DAnlys_Watershed(bw,Opt_Watershed);

% Visualise results
figure, isosurface(x,y,z,bw_S,0.5), axis equal, title('BW')
xlabel x, ylabel y, zlabel z
xlim(lims), ylim(lims), zlim(lims)
view(3), camlight, lighting gouraud

% Include watershed lines
bw_S_Ld = bw_S;
bw_S_Ld(~bw_R.Ld_Img) = 1;
% Transform to label matrix
[L,~] = bwlabeln(bw_S);
ViewPnl = uipanel(figure,'Title','Labeled Volume');
h = labelvolshow(L,bw_S_Ld,'Parent',ViewPnl);

% Create an array of camera positions around the unit circle
vec = linspace(0,2*pi(),60)';
myPosition = [cos(vec) sin(vec) ones(size(vec))]*5;
for idx = 1:60
    h.CameraPosition = myPosition(idx,:);
    pause(0.05)
end

%% Demo 2: Capsule pellet segmentation
%  
%  % % Reference (open access):
%  Doerr, F. J. S., Florence, A. J. (2020)
%  A micro-XRT image analysis and machine learning methodology for the characterisation of multi-particulate capsule formulations. 
%  International Journal of Pharmaceutics: X. 
%  https://doi.org/10.1016/j.ijpx.2020.100041
%  Data repository: https://doi.org/10.15129/e5d22969-77d4-46a8-83b8-818b50d8ff45
%  Video Abstract: https://strathprints.strath.ac.uk/id/eprint/71463

fprintf('\n - - -\nDemo 2: Capsule pellet segmentation\n - - -\n')

% Load 3D XRT image data
ImgStack_Load_Options.ExpShorthand = 'Demo';
ImgStack_Load_Options.AppShorthand = 'Capsule pellet segmentation';
ImgStack_Load_Options.Print_Ctrl_On = true;
ImgStack_Load_Options.pool_mode = true;
ImgStack_Load_Options.ImgFormat = 'bmp';
ImgStack_Load_Options.path_ImgFolder = fullfile(Opt.path_Mat_Demo_Data,'XRT_Capsule\V_CS_R8');
[V_CS] = XRT_3DAnlys_ImgStack_Load(ImgStack_Load_Options);

ImgStack_Load_Options.ImgFormat = 'bmp';
ImgStack_Load_Options.path_ImgFolder = fullfile(Opt.path_Mat_Demo_Data,'XRT_Capsule\V_CP_ROI_R8');
[V_CP_ROI] = XRT_3DAnlys_ImgStack_Load(ImgStack_Load_Options);

% Watershed Segmentation
Opt_Watershed.Print_Ctrl_On = true;
Opt_Watershed.bwdist_method = 'euclidean';
Opt_Watershed.imext_Size = 4;
Opt_Watershed.R_Ctrl_On = true;
[V_CS_S,R] = XRT_3DAnlys_Watershed(V_CP_ROI,Opt_Watershed);

% Visualise results
[L,numObj] = bwlabeln(V_CS_S);
for k = 1:numObj
    L(L==k) = ceil(rand(1)*127);
end
ViewPnl = uipanel(figure,'Title','Labeled Volume');
h = labelvolshow(L,V_CS_S,'Parent',ViewPnl);


