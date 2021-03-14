%------------------------------------------------------------------------------------------------
% Code written by Frederik Doerr, Aug 2020 (MATLAB R2019b)
% Application: XRT Image and Data Analysis Framework
% 
% https://github.com/frederik-d
% Contact: frederik.doerr(at)strath.ac.uk | CMAC (http://www.cmac.ac.uk/)

fclose('all');
close all
clear all %#ok<CLALL>
clc

set(0,'DefaultFigureVisible','off');
% set(0,'DefaultFigureVisible','on');



%% Setup
fprintf('--- %s ---\n',mfilename())
fprintf('- Code written by Frederik Doerr, Aug 2020 (MATLAB R2019b)\n')
fprintf('- Application: XRT Image and Data Analysis Framework\n\n')
fprintf('- https://github.com/frederik-d\n')
fprintf('- Contact: info@cmac.ac.uk | CMAC (http://www.cmac.ac.uk/)\n')

% Main folder location
path = matlab.desktop.editor.getActiveFilename;
[XRT_Opts.filesPath_MAIN,~,~] = fileparts(path);

% Fixed File Folder Locations
XRT_Opts.filesPath_Mat_SR = fullfile(XRT_Opts.filesPath_MAIN,'_SubRoutines');
XRT_Opts.filesPath_FileExchange = fullfile(XRT_Opts.filesPath_MAIN,'_FileExchange');
XRT_Opts.filesPath_ImgPrc_Parameter = fullfile(XRT_Opts.filesPath_MAIN,'_ImgPrc_ParameterFiles');
XRT_Opts.filesPath_ImgPrc_Parameter_Archive = fullfile(XRT_Opts.filesPath_MAIN,'_Export\_Archive');
XRT_Opts.filesPath_Export = fullfile(XRT_Opts.filesPath_MAIN,'_Export');

% Excel File: XRT_ParameterImport
XRT_Opts.ImgPrc_Parameter_name = 'XRT_ParameterImport.xlsx';
XRT_Opts.ImgPrc_Parameter_book = fullfile(XRT_Opts.filesPath_ImgPrc_Parameter,XRT_Opts.ImgPrc_Parameter_name);
XRT_Opts.ImgPrc_3DAnlys = readtable(XRT_Opts.ImgPrc_Parameter_book,'Sheet','3DAnlys');

addpath(genpath(XRT_Opts.filesPath_FileExchange))
addpath(genpath(XRT_Opts.filesPath_Mat_SR))

% Get local computer ID
addpath(fullfile(XRT_Opts.filesPath_FileExchange,'getComputerName'))
XRT_Opts.localPC_name = getComputerName();
XRT_Opts.datestring = datestr(datetime('now'),'yyyy-mm-dd');
XRT_Opts.idx_MasterDB = length(dir(fullfile(XRT_Opts.filesPath_Export,sprintf('XRT_3DAnlys_MasterDB_%s_%s_*.csv',XRT_Opts.datestring,XRT_Opts.localPC_name))))+1;
XRT_Opts.filesName_MasterDB = sprintf('XRT_3DAnlys_MasterDB_%s_%s_%03.0f.csv',XRT_Opts.datestring,XRT_Opts.localPC_name,XRT_Opts.idx_MasterDB);
XRT_Opts.filesPath_MasterDB = fullfile(XRT_Opts.filesPath_Export,XRT_Opts.filesName_MasterDB);
XRT_Opts.numExp = find(cellfun(@isempty,XRT_Opts.ImgPrc_3DAnlys.ExpShorthand),1,'first')-1;
if isempty(XRT_Opts.numExp)
    XRT_Opts.numExp = size(XRT_Opts.ImgPrc_3DAnlys.ExpShorthand,1);
end

fprintf('--- %s. Setup Complete ---\n',mfilename())
pause(2)

%% Enter Analysis
for Batch3DAnlys_Iter = 1:XRT_Opts.numExp

    
    %% Initialize Anlysis Run
    clearvars -except XRT_3DAnlys_MasterDB XRT_Opts Batch3DAnlys_Iter pool
    
    %% Load Parameter and Create ExpShorthand Folder Directory
    XRT_3DAnlys_Opts = struct();
    XRT_ImgPrc_Opts = struct();
    [XRT_3DAnlys_Opts,XRT_ImgPrc_Opts] = XRT_3DAnlys_ParameterLoader(XRT_Opts,Batch3DAnlys_Iter,XRT_3DAnlys_Opts,XRT_ImgPrc_Opts);
    
    % Export Datatable - Create master-cell structure for outputs
    if ~exist('XRT_3DAnlys_MasterDB','var')
        XRT_3DAnlys_MasterDB = { ...
            'ExpShorthand', 'DateNumber', 'ImagePixelSize','ImagePixelSize_rS', 'DataResizeFactor',...
            'size_ImageBox_X','size_ImageBox_Y','size_ImageBox_Z', ...
            'V_BoundingBox_MaxLength','V_BoundingBox_MinLength', ...
            'V_ROI_BoundingBox_MaxLength','V_ROI_BoundingBox_MinLength',...
            'NumEl_Volume','V_Volume_um','V_ROI_Volume_um','V_ROI_CH_Volume_um', ...
            'V_Surface_um','V_ROI_Surface_um', 'V_max_FeretDiameter',...
            'V_radius_spherefit_um','V_ROI_radius_spherefit_um', ...
            'V_radii_elpsfit_1_um','V_radii_elpsfit_2_um','V_radii_elpsfit_3_um','V_chi2_elpsfit_um', ...
            'V_ROI_radii_elpsfit_1_um','V_ROI_radii_elpsfit_2_um','V_ROI_radii_elpsfit_3_um','V_ROI_chi2_elpsfit_um',...
            'V_AspectRatio','V_ROI_AspectRatio','V_ROI_Sphericity', ...
            'V_Poros_Volume_um','V_Poros_closed_Volume_um','V_Poros_open_Volume_um', ...
            'V_Solidity','V_Solidity_CH','V_ROI_Solidity_CH','V_Extent','V_ROI_Extent','V_EquivDiameter','V_ROI_EquivDiameter', ...
            'V_Orientation_Phi','V_Orientation_Theta','V_Orientation_Psi','V_ROI_Orientation_Phi','V_ROI_Orientation_Theta','V_ROI_Orientation_Psi', ...
            };
    end
    XRT_3DAnlys_Opts.idx_newRow = size(XRT_3DAnlys_MasterDB,1)+1;
    
    
    if isfield(XRT_Opts,'ImgPrc_3DAnlys')
        XRT_3DAnlys_Opts.ExpShorthand = char(XRT_Opts.ImgPrc_3DAnlys.ExpShorthand(Batch3DAnlys_Iter));
        XRT_3DAnlys_Opts.Info = char(XRT_Opts.ImgPrc_3DAnlys.Info(Batch3DAnlys_Iter));
    else
      	prompt = {'ExpShorthand'};
        prompttitle = 'ExpShorthand';
        num_lines = 1; 
        ExpShorthand = inputdlg(prompt,prompttitle,num_lines);
        options.Interpreter='tex';
        XRT_3DAnlys_Opts.ExpShorthand = char(ExpShorthand);
    end
    

    %% Start log-file
    diary(sprintf('LOG_XRT_3DAnlys_%s_%s.txt',XRT_3DAnlys_Opts.ExpShorthand,XRT_Opts.datestring))
    diary on
    
    % Get local computer ID
    addpath(XRT_Opts.filesPath_FileExchange)
    XRT_3DAnlys_Opts.localPC_name = XRT_Opts.localPC_name;
    
    XRT_3DAnlys_Opts.DateNumber = datenum(date);
    
    % TIC/ TOC: Set time zero
    XRT_3DAnlys_Opts.tic = tic;
    XRT_3DAnlys_Opts.PrcTime_Logger = {};
    XRT_3DAnlys_Opts.PrcTime_Logger{1,1} = 'MainScript';
    XRT_3DAnlys_Opts.PrcTime_Logger{1,2} = 'SubScript';
    XRT_3DAnlys_Opts.PrcTime_Logger{1,3} = 'Time [sec]';
    XRT_3DAnlys_Opts.PrcTime_Logger{end+1,1} = mfilename();
    XRT_3DAnlys_Opts.PrcTime_Logger{end,2} = '--- Start ----';
    XRT_3DAnlys_Opts.PrcTime_Logger{end,3} = toc(XRT_3DAnlys_Opts.tic);

	% log-file
    fprintf('START ---- %s ----  (PC: %s, Elapsed time: %.0f sec)\n',XRT_3DAnlys_Opts.ExpShorthand,XRT_3DAnlys_Opts.localPC_name,toc(XRT_3DAnlys_Opts.tic))
    fprintf('Date: %s \n',datestr(datetime('now')))
    fprintf('ExpShorthand: %s \n',XRT_3DAnlys_Opts.ExpShorthand)
    fprintf('Info: %s \n',XRT_3DAnlys_Opts.Info)

    % Update XRT_3DAnlys_MasterDB
    XRT_3DAnlys_MasterDB{XRT_3DAnlys_Opts.idx_newRow, strcmp('ExpShorthand',XRT_3DAnlys_MasterDB(1,:))} = XRT_3DAnlys_Opts.ExpShorthand; %#ok<*SAGROW>
    XRT_3DAnlys_MasterDB{XRT_3DAnlys_Opts.idx_newRow, strcmp('DateNumber',XRT_3DAnlys_MasterDB(1,:))} = XRT_3DAnlys_Opts.DateNumber;
    XRT_3DAnlys_MasterDB{XRT_3DAnlys_Opts.idx_newRow, strcmp('ImagePixelSize',XRT_3DAnlys_MasterDB(1,:))} = XRT_3DAnlys_Opts.ImagePixelSize;    
    XRT_3DAnlys_MasterDB{XRT_3DAnlys_Opts.idx_newRow, strcmp('ImagePixelSize_rS',XRT_3DAnlys_MasterDB(1,:))} = XRT_3DAnlys_Opts.ImagePixelSize*XRT_3DAnlys_Opts.DataResizeFactor;
    XRT_3DAnlys_MasterDB{XRT_3DAnlys_Opts.idx_newRow, strcmp('DataResizeFactor',XRT_3DAnlys_MasterDB(1,:))} = XRT_3DAnlys_Opts.DataResizeFactor;

    % Save Process Parameter
    XRT_3DAnlys_ParameterSave(XRT_Opts,XRT_3DAnlys_Opts,XRT_ImgPrc_Opts,Batch3DAnlys_Iter)

    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    %% XRT_3DAnlys_VolumeAnlys
    run XRT_3DAnlys_VolumeAnlys
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    %% XRT_3DAnlys_PorosAnlys
    run XRT_3DAnlys_PorosAnlys
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    %% XRT_3DAnlys_BinSegment
    run XRT_3DAnlys_BinSegment

    %% Export latest XRT_3DAnlys_MasterDB to csv report-file

    % Single Analysis Table
    XRT_3DAnlys_MasterDB_Table_Exp = cell2table(XRT_3DAnlys_MasterDB(end,:),...
        'VariableNames',XRT_3DAnlys_MasterDB(1,:));
    
    formatOut = 'yyyy-mm-dd_HH-MM-SS';
    datestring = datestr(datetime('now'),formatOut);
    filename = fullfile(pwd,sprintf('MasterDB_%s_%s.csv',XRT_3DAnlys_Opts.ExpShorthand,datestring));
    writetable(XRT_3DAnlys_MasterDB_Table_Exp,filename)
    
    % Full Table
    XRT_3DAnlys_MasterDB_Table = cell2table(XRT_3DAnlys_MasterDB(2:end,:),...
        'VariableNames',XRT_3DAnlys_MasterDB(1,:));
    writetable(XRT_3DAnlys_MasterDB_Table,XRT_Opts.filesPath_MasterDB)
    
    
    XRT_3DAnlys_Opts.PrcTime_Logger{end+1,1} = mfilename();
    XRT_3DAnlys_Opts.PrcTime_Logger{end,2} = '--- End ----';
    XRT_3DAnlys_Opts.PrcTime_Logger{end,3} = toc(XRT_3DAnlys_Opts.tic);
    
    % Save XRT_3DAnlys_Options.PrcTime_Logger
    datestring = datestr(date,'yyyy-mm-dd');
    filename = fullfile(pwd,sprintf('PrcTime_Logger_%s_%s.csv',XRT_3DAnlys_Opts.ExpShorthand,datestring));
    fid_PrcTime_Logger = fopen(filename,'w+'); 
    fprintf(fid_PrcTime_Logger,'%s, %s, %s \n',XRT_3DAnlys_Opts.PrcTime_Logger{1,:});
    for i_PrcTime_Logger = 2:size(XRT_3DAnlys_Opts.PrcTime_Logger,1)
        fprintf(fid_PrcTime_Logger,'%s, %s, %f \n',XRT_3DAnlys_Opts.PrcTime_Logger{i_PrcTime_Logger,:});
    end
    fclose(fid_PrcTime_Logger);
    
    fprintf('\n- - - - - - - - - - - -   - - - - - - - - - - - -   - - - - - - - - - - - -\n')
   	fprintf('XRT_3DAnlys: COMPLETE %d of %d (%s)\n',Batch3DAnlys_Iter,XRT_Opts.numExp,XRT_3DAnlys_Opts.ExpShorthand)
    fprintf('- - - - - - - - - - - -   - - - - - - - - - - - -   - - - - - - - - - - - -\n\n')
    pause(10)
    diary off
 	clc
    
end

set(0,'DefaultFigureVisible','on');
fprintf('--- Export XRT_3DAnlys_MasterDB. Exit %s ---\n',mfilename())


% MathWorks Licence
% Copyright (c) 2016-2020, Frederik Doerr
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
%     * Neither the name of the University of Strathclyde, CMAC nor the names
%       of its contributors may be used to endorse or promote products derived
%       from this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

