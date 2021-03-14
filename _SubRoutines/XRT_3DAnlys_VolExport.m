function XRT_3DAnlys_VolExport(V,Opt)
%XRT_3DAnlys_VolExport Export Surface rendered object.
%   XRT_3DAnlys_VolExport(V,Opt) function to generate surface rendered
%   object and mesh with stl-export
%

%   OPTIONS:
%      'Export_Mesh_On'     Use vol2mesh
%      'Print_Log_On'       Enable fprintf outputs to monitor/record
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

if ~isfield(Opt,'Print_Log_On')
    Opt.Print_Log_On = false;
end
Opt.Print_Str = sprintf('%s%s:',Opt.Print_Str,mfilename());

AxisOff_On = Opt.AxisOff_On;
ParforComp = Opt.ParforComp;

if ~isfield(Opt,'Export_Mesh_On')
   Opt.Export_Mesh_On = false;
end

if ~exist(Opt.path_ImgFolder,'dir')
    mkdir(Opt.path_ImgFolder);
end

v = double(V);

%# visualize the volume
figure1 = figure;
axes1 = axes('Parent',figure1);
p = patch( isosurface(v,0) );                 %# create isosurface patch
isonormals(v, p)                              %# compute and set normals
set(p, 'FaceColor','k', 'EdgeColor','none')   %# set surface props
daspect([1 1 1])                              %# axes aspect ratio
view(3), axis vis3d tight, box on, grid on    %# set axes props
camproj perspective                           %# use perspective projection
camlight, lighting phong, alpha(.75)           %# enable light, set transparency
if AxisOff_On
    set(axes1,'XColor',[1 1 1],'YColor',[1 1 1],'ZColor',[1 1 1]);
end
xlabel('x [px]')
ylabel('y [px]')
zlabel('z [px]')

saveas(figure1,fullfile(Opt.path_ImgFolder,sprintf('%s_1',Opt.ImgName)),'png');   

AZ = 52.5; EL = 30;
view([AZ,EL])
saveas(figure1,fullfile(Opt.path_ImgFolder,sprintf('%s_2',Opt.ImgName)),'png');   

AZ = 142.5; EL = 30;
view([AZ,EL])
saveas(figure1,fullfile(Opt.path_ImgFolder,sprintf('%s_3',Opt.ImgName)),'png');   


if Opt.Export_Mesh_On
    addpath(fullfile(Opt.path_FileExchange,'iso2mesh-1.8-win32\iso2mesh'))
    ix = 1:size(V,1);
    iy = 1:size(V,2);
    iz = 1:size(V,3);
    opt = 2;
    maxvol = 2;
    dofix = 1;
%     method,isovalues
    [node,elem,face]=vol2mesh(V,ix,iy,iz,opt,maxvol,dofix);
    
    addpath(fullfile(Opt.path_FileExchange,'stlwrite'))
    filename = sprintf('%s.stl',Opt.ImgName);
    stlwrite(filename)
end
