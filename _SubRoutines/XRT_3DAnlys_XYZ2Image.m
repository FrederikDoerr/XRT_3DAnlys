function[V] = XRT_3DAnlys_XYZ2Image(XYZ,Opt)
%XRT_3DAnlys_XYZ2Image Convert list to 3D matrix.
%   [V] = XRT_3DAnlys_XYZ2Image(XYZ,Opt) uses values in list XYZ (3 or 4
%   columns) to generate 3D image data matrix.
%
%
%   Author
%   --------
%   Frederik Doerr, Aug 2020 (MATLAB R2019b)
%   frederik.doerr(at)strath.ac.uk | CMAC (http://www.cmac.ac.uk/)
%   Github: https://github.com/frederik-d/XRT_3DAnlys

if nargin < 2
    Opt = [];
end

[numRow,numCol] = size(XYZ);

if ~isfield(Opt,'Int_ScalingFactor') && numCol == 4
    Opt.Int_ScalingFactor = 1; %255/max(XYZ(:,4));
end


if numCol == 4
    % Additional intensity information
    max_x = max(XYZ(:,1));
    max_y = max(XYZ(:,2));
    max_z = max(XYZ(:,3));
    
    % 1) Create empty image volume
    V = uint8(zeros(max_y,max_x,max_z));
    
    % 2) Fill Image Volume
    for i = 1:numRow
        V(XYZ(i,2),XYZ(i,1),XYZ(i,3)) = XYZ(i,4)*(Opt.Int_ScalingFactor);
    end
    
elseif numCol == 3
    max_x = max(XYZ(:,1));
    max_y = max(XYZ(:,2));
    max_z = max(XYZ(:,3));
    
    % 1) Create empty image volume
    V = zeros(max_y,max_x,max_z,'logical');
    
    % 2) Fill Image Volume
    for i = 1:numRow
        V(XYZ(i,2),XYZ(i,1),XYZ(i,3)) = 1;
    end

end


