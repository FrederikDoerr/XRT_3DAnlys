function[V,R] = XRT_3DAnlys_NoiseReduction(V,Opt)
%XRT_3DAnlys_NoiseReduction Noise Reduction of V.
%   [V,R] = XRT_3DAnlys_NoiseReduction(V,Opt) for 3D noise reduction of 
%   binary 3D image dataset V.
%
%   OPTIONS:
%      'Method'             Selected noise reduction method:
%                                     'clearborder'
%                                     'imclose'
%                                     'imopen'
%                                     'sweep'
%
%   Author
%   --------
%   Frederik Doerr, Aug 2020 (MATLAB R2019b)
%   frederik.doerr(at)strath.ac.uk | CMAC (http://www.cmac.ac.uk/)
%   Github: https://github.com/frederik-d/XRT_3DAnlys



R = struct();
if isfield(Opt,'Processing_mode')
    switch lower(Opt.Processing_mode)
        case{'cpu'}
            Processing_mode_cpu = true;
    end
end


switch lower(Opt.Method)
    case {'clearborder'}
        addpath(fullfile(Opt.filesPath_FileExchange,'bwclearborder'))
        V = bwclearborder(V,26);
        
    case {'imclose'}
        se = strel('sphere',Opt.Close_Size);
        V = imclose(V,se); 
        
    case {'imopen'}
        se = strel('sphere',Opt.Open_Size);
        V = imopen(V,se);
        
    case {'sweep'}
        addpath(fullfile(Opt.filesPath_FileExchange,'bwareafilt2'))
        range = [];
        d = 'largest';
        p = 1;
        V = bwareafilt2(V, range, p, d);
end  

