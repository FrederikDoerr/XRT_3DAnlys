function [V] = XRT_3DAnlys_Shrinkwrap(V,Opt)
%XRT_3DAnlys_Shrinkwrap Shrinkwrap ROI of V.
%   [V] = XRT_3DAnlys_ROI_Reduction(V,Opt) function to establish a 3D ROI
%   using a shrinkwrap method for V
%
%   OPTIONS:
%      'Method'             Select method:
%                                   'shrinkwrap_op'
%                                   'stepwise_imclose'
%                                   'stepwise'
%      'strel_type'         Defined structuring element type for imclose
%                           step (default: 'cube')
%      'strel_size'         Defined structuring element size for each
%                           imdilate step (default: 3)
%      'Print_Log_On'       Enable fprintf outputs to monitor/record
%                           progress (default: false)
%      'ExpShorthand'       Sample ID (default: '')
%      'AppShorthand'       Application ID (default: mfilename())
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


if ~isfield(Opt,'imclose_strel_type')
    Opt.imclose_strel_type = 'disk';
    warning('%s strel_type NOT defined (default: disk) \n',Opt.Print_Str)
end


switch lower(Opt.Method)
    case {'shrinkwrap_op'}
        if Opt.Print_Log_On 
            fprintf('%s Enter %s, se_diskRadius >> %.0f << \n', ...
                Opt.Print_Str, Opt.Method,Opt.imclose_strel_size)
        end
        V = imfill(V,'holes');
        se = strel(Opt.imclose_strel_type,Opt.imclose_strel_size); % An offsetstrel object represents a nonflat morphological structuring element
        V_prc1 = imclose(V,se);
        if Opt.Print_Log_On 
            fprintf('%s %s 1) imclose COMPLETE \n',Opt.Print_Str, ...
                Opt.Method)
        end

        V_prc2 = imfill(V_prc1,'holes');      
        if Opt.Print_Log_On 
            fprintf('%s %s 2) imfill COMPLETE \n',Opt.Print_Str, ...
                Opt.Method)
        end
        
        % Match morph operation to internal chnages
        V_prc_Delta = imabsdiff(V_prc2,V);
                       
        % Detect areas with change (markers)
        V_processed_Delta_marker =  imabsdiff(V_prc1,V_prc2);
               
        % Use markers 
        V_prc_Delta = imreconstruct(V_processed_Delta_marker,V_prc_Delta);
        if Opt.Print_Log_On 
            fprintf('%s %s 3) match COMPLETE \n',Opt.Print_Str, ...
                Opt.Method)
        end
        
        V = imadd(V,V_prc_Delta);
        
        if ~islogical(V)
            V = logical(V);
        end
        if Opt.Print_Log_On 
            fprintf('%s %s, se_diskRadius >> %.0f << \n',Opt.Print_Str, ...
            Opt.Method,Opt.imclose_strel_size)
        end

    case {'stepwise_imclose'}        
        imclose_strel_size = 1:Opt.imclose_strel_size;
        for i = 1:Opt.imclose_strel_size
            se = strel(Opt.imclose_strel_type,imclose_strel_size(i)); % An offsetstrel object represents a nonflat morphological structuring element

            V_prc1 = imclose(V, se);
            if Opt.Print_Log_On 
                fprintf('%s %s (%.0f / %.0f) 1) imclose COMPLETE \n',Opt.Print_Str,Opt.Method,i,Opt.imclose_strel_size)
            end

            V_prc2 = imfill(V_prc1,'holes');
            V_prc2(V_prc1) = 0;
            if Opt.Print_Log_On 
                fprintf('%s %s (%.0f / %.0f) 2) imfill COMPLETE \n',Opt.Print_Str,Opt.Method,i,Opt.imclose_strel_size)
            end
            
            V(V_prc2) = 1;
            if Opt.Print_Log_On 
                fprintf('%s %s (%.0f / %.0f) \n',Opt.Print_Str,Opt.Method,i,Opt.imclose_strel_size)
            end
        end
         
    case {'shrinkwrap'}
        se = strel(Opt.imclose_strel_type,Opt.imclose_strel_size); % An offsetstrel object represents a nonflat morphological structuring element
        
        V = imclose(V, se);
        if Opt.Print_Log_On 
            fprintf('%s %s 1) imclose COMPLETE \n',Opt.Print_Str,Opt.Method)
        end

        V = imfill(V,'holes');
        if Opt.Print_Log_On 
            fprintf('%s %s 2) imfill COMPLETE \n',Opt.Print_Str,Opt.Method)
            fprintf('%s %s, se_diskRadius >> %.0f << \n',Opt.Print_Str, ...
            Opt.Method,Opt.imclose_strel_size)
        end
        
    otherwise
        fprintf('%s No Method Matched.\n',Opt.Print_Str)
end