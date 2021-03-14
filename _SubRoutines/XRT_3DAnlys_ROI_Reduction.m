function[V,R] = XRT_3DAnlys_ROI_Reduction(V,Opt)
%XRT_3DAnlys_ROI_Reduction ROI reduction of V.
%   [V,R] = XRT_3DAnlys_ROI_Reduction(V,Opt) function to reduce the images 
%   to their minimum size (bounding box). Various options to force ROI
%   reduction where pre-defined or add spacers to padarray with background
%   pixel.
%
%   OPTIONS:
%      'numSpacer'          Number of blank spacers.
%      'numImg_Start'      	Selected parameter z-dim start (dim = 3).
%      'num_Img_End'        Selected parameter z-dim end (dim = 3).
%      'row_Min'            Selected parameter y-dim start (dim = 1).
%      'row_Max'            Selected parameter y-dim end (dim = 1).
%      'column_Min'         Selected parameter x-dim start (dim = 2).
%      'column_Max'         Selected parameter x-dim start (dim = 2).
%      'Print_Log_On'      Enable fprintf outputs to monitor/record
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


if isfield(Opt,'numSpacer') && ~isnan(Opt.numSpacer)
    numSpacer = Opt.numSpacer;
else
    numSpacer = 0;
end


if ~any(V(:)) && ~(isfield(Opt,'V_any') && Opt.V_any)
    R.V_any = false;
    if Opt.Print_Log_On
        warning('%s FAILED (~any(V(:)))\n',Opt.Print_Str)
    end
    
    R.numImg_Start = nan;
    R.num_Img_End = nan;
    R.row_Min = nan;
    R.row_Max = nan;
    R.column_Min = nan;
    R.column_Max = nan;
    R.numSpacer = nan;
else
    R.V_any = true;
    [~,~,numImg] = size(V);
    if numImg > 1
        if isfield(Opt,'numImg_Start') && ...
            isfield(Opt,'num_Img_End') && ...
            isfield(Opt,'row_Min') && ...
            isfield(Opt,'row_Max') && ...
            isfield(Opt,'column_Min') && ...
            isfield(Opt,'column_Max')

            Opt.ROI_Space_FullyDefined = true;
        else
            Opt.ROI_Space_FullyDefined = false;
        end

        if ~Opt.ROI_Space_FullyDefined
            row_Min_List = [];
            row_Max_List = [];
            column_Min_List = [];
            column_Max_List = [];
            Control_ParticleDetect_list = zeros(1,numImg);
            for k = 1:numImg
                I = V(:,:,k);
                if max(max(I))
                    % Find limits of ROI region
                    vertical = any(I, 2);
                    horizontal = any(I, 1);

                    row_Min = find(vertical, 1, 'first'); % Y1
                    row_Min_List = [row_Min_List,row_Min]; %#ok<*AGROW>
                    row_Max = find(vertical, 1, 'last'); % Y2
                    row_Max_List = [row_Max_List,row_Max];
                    column_Min = find(horizontal, 1, 'first'); % X1
                    column_Min_List = [column_Min_List,column_Min];
                    column_Max = find(horizontal, 1, 'last'); % X2
                    column_Max_List = [column_Max_List,column_Max];

                    Control_ParticleDetect_list(k) = true;
                else
                    Control_ParticleDetect_list(k) = false;
                end
            end
            numImg_Start = find(Control_ParticleDetect_list, 1, 'first');
            num_Img_End = find(Control_ParticleDetect_list, 1, 'last');
        end

        if isfield(Opt,'numImg_Start') && ~isnan(Opt.numImg_Start)
            numImg_Start = Opt.numImg_Start;
        end

        if isfield(Opt,'num_Img_End') && ~isnan(Opt.num_Img_End)
            num_Img_End = Opt.num_Img_End;
        end

        if isfield(Opt,'row_Min') && ~isnan(Opt.row_Min)
            row_Min = Opt.row_Min;
        else
            row_Min = min(row_Min_List);
        end

        if isfield(Opt,'row_Max') && ~isnan(Opt.row_Max)
            row_Max = Opt.row_Max;
        else
            row_Max = max(row_Max_List);
        end

        if isfield(Opt,'column_Min') && ~isnan(Opt.column_Min)
            column_Min = Opt.column_Min;
        else
            column_Min = min(column_Min_List);
        end

        if isfield(Opt,'column_Max') && ~isnan(Opt.column_Max)
            column_Max = Opt.column_Max;
        else
            column_Max = max(column_Max_List);
        end
        
        V = V(:,:,numImg_Start:num_Img_End);
        V = V(row_Min:row_Max,column_Min:column_Max,:);
        if numSpacer > 0 
            for k = 1:numSpacer
                [numRow,numCol,~] = size(V);
                spacer_V = zeros(numRow,numCol,'logical');
                V = cat(3,spacer_V, V);
                V = cat(3, V, spacer_V);
                % Expand V_crop with zeros
                V = padarray(V,[1 1]);
            end
        end
        R.numImg_Start = numImg_Start;
        R.num_Img_End = num_Img_End;
        R.row_Min = row_Min;
        R.row_Max = row_Max;
        R.column_Min = column_Min;
        R.column_Max = column_Max;
        R.numSpacer = numSpacer;

        if Opt.Print_Log_On
            fprintf('%s COMPLETE\n',Opt.Print_Str)
        end
    else
        if  isfield(Opt,'row_Min') && ...
            isfield(Opt,'row_Max') && ...
            isfield(Opt,'column_Min') && ...
            isfield(Opt,'column_Max')

            Opt.ROI_Space_FullyDefined = true;
        else
            Opt.ROI_Space_FullyDefined = false;
        end

        if ~Opt.ROI_Space_FullyDefined
            % Find limits of ROI region
            vertical = any(V, 2);
            horizontal = any(V, 1);

            row_Min = find(vertical, 1, 'first'); % Y1
            row_Max = find(vertical, 1, 'last'); % Y2
            column_Min = find(horizontal, 1, 'first'); % X1
            column_Max = find(horizontal, 1, 'last'); % X2
        end


        if isfield(Opt,'row_Min') && ~isnan(Opt.row_Min)
            row_Min = Opt.row_Min;
        end

        if isfield(Opt,'row_Max') && ~isnan(Opt.row_Max)
            row_Max = Opt.row_Max;
        end

        if isfield(Opt,'column_Min') && ~isnan(Opt.column_Min)
            column_Min = Opt.column_Min;
        end

        if isfield(Opt,'column_Max') && ~isnan(Opt.column_Max)
            column_Max = Opt.column_Max;
        end

        V = V(row_Min:row_Max,column_Min:column_Max);
        if numSpacer > 0 
            for k = 1:numSpacer
                % Expand V_crop with zeros
                V = padarray(V,[1 1]);
            end
        end
        R.numImg_Start = nan;
        R.num_Img_End = nan;
        R.row_Min = row_Min;
        R.row_Max = row_Max;
        R.column_Min = column_Min;
        R.column_Max = column_Max;
        R.numSpacer = numSpacer;
    end
end
