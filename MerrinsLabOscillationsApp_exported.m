classdef MerrinsLabOscillationsApp_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                  matlab.ui.Figure
        GridLayout                matlab.ui.container.GridLayout
        ChooseFiletoImportButton  matlab.ui.control.Button
        FileNameEditField         matlab.ui.control.EditField
        UITableImported           matlab.ui.control.Table
        UITableTimeWindows        matlab.ui.control.Table
        DataSummaryPanel          matlab.ui.container.Panel
        UITableOutput             matlab.ui.control.Table
        SaveFileEditField         matlab.ui.control.EditField
        SaveOutputButton          matlab.ui.control.Button
        TabGroup                  matlab.ui.container.TabGroup
        PlotsTab                  matlab.ui.container.Tab
        UIAxesAnal                matlab.ui.control.UIAxes
        WaveletviewTab            matlab.ui.container.Tab
        UIAxes4                   matlab.ui.control.UIAxes
        UIAxes3                   matlab.ui.control.UIAxes
        UIAxesWVLT                matlab.ui.control.UIAxes
        ControlsPanel             matlab.ui.container.Panel
        UpdateOutputButton        matlab.ui.control.Button
        DetrendCheckBox           matlab.ui.control.CheckBox
        UITableAnalVals           matlab.ui.control.Table
        Instructions              matlab.ui.control.Label
        SensorTypeDropDown        matlab.ui.control.DropDown
        SensorTypeDropDownLabel   matlab.ui.control.Label
        UIAxesSelector            matlab.ui.control.UIAxes
        ContextMenu               matlab.ui.container.ContextMenu
        AddRow                    matlab.ui.container.Menu
        DelRow                    matlab.ui.container.Menu
    end

    
    % This app written by J.D. Rogers <rogersjd@gmail.com> on Oct 11, 2025
    % Based on Merrins lab code, converted to GUI
    % last updated: 20251016


    % shared variables across functions
    properties (Access = private)
        cols2preview; % which columns in the imported table to plot in the preview
        rows2plot=[1]; % which row in the AnalVal table to plot in the anal window 
        alltimes;  % time in seconds for the entire data set
        timewindows=table(); % table containing the start and end times for each time window
        datapath; % path where data is located and output is saved
        OutTable = table() %store all the output info in this table
    end


    methods (Access = private)

        function updatePreviewPlot(app)
            % replot the full time sequence for the given selected columns
            plot(app.UIAxesSelector,app.UITableImported.Data{:,1},app.UITableImported.Data{:,app.cols2preview});
            legend(app.UIAxesSelector,app.UITableImported.Data.Properties.VariableDescriptions{app.cols2preview});
            title(app.UIAxesSelector,'select column(s) on left to plot, click to select time windows');
            ylim(app.UIAxesSelector,'auto');
            yl = app.UIAxesSelector.YLim(); % get the ylims so we can make time window patches spanning the y axis
            
            % grab the time windows and convert to index vals
            tmin = cell2mat(app.UITableTimeWindows.Data{:,2});
            tmax = cell2mat(app.UITableTimeWindows.Data{:,3});
            tminind = [];
            tmaxind = [];

            % add color patches for time windows in the UIAxesSelector
            hold(app.UIAxesSelector,"on")
            colors = colororder;
            for ii = 1:size(app.UITableTimeWindows.Data,1)
                tminind(ii) = find(app.alltimes>=tmin(ii),1,'first');
                tmaxind(ii) = find(app.alltimes<=tmax(ii),1,'last');
                
                patch(app.UIAxesSelector,[tmin(ii) tmin(ii) tmax(ii) tmax(ii)], [yl(1) yl(2) yl(2) yl(1)], ...
                      colors(ii,:)*0.25, 'FaceAlpha',0.2, 'EdgeColor',colors(ii,:)*0.5, 'DisplayName',app.UITableTimeWindows.Data{ii,1}{1});
            end
            ylim(app.UIAxesSelector,yl);
            hold(app.UIAxesSelector,"off")
        end

        function updateUITableAnalVals(app)
            %TableAnalVals =  app.UITableAnalVals.Data;
            nTWindows = size(app.UITableTimeWindows.Data,1);
            % For each time window row, make a set of nRegions of rows to
            % the analysis table
            nRegions = size(app.UITableImported.Data,2)-1; % (-1) since first column is time
            for ii=1:nTWindows
                % Add nRegions of rows for each time window to the tale
                app.UITableAnalVals.Data.TWindow((ii-1)*nRegions+[1:nRegions])=app.UITableTimeWindows.Data{ii,1};
                
                % set the start and end time (and indicesfor each table row
                tmin = app.UITableTimeWindows.Data{ii,2}{1}; % index {1} to return value instead of cell array
                tmax = app.UITableTimeWindows.Data{ii,3}{1};
                app.UITableAnalVals.Data.tmin((ii-1)*nRegions+[1:nRegions])=tmin; % store these vals in each row
                app.UITableAnalVals.Data.tmax((ii-1)*nRegions+[1:nRegions])=tmax;
                tminind = find(app.alltimes>=tmin,1,'first');
                tmaxind = find(app.alltimes<=tmax,1,'last');
                app.UITableAnalVals.Data.tminind((ii-1)*nRegions+[1:nRegions])= tminind; % store these vals in each row
                app.UITableAnalVals.Data.tmaxind((ii-1)*nRegions+[1:nRegions])= tmaxind;
                
                xdata=app.UITableImported.Data{tminind:tmaxind,1};
                xdata=xdata-xdata(1); 
                app.UITableAnalVals.Data.xdata((ii-1)*nRegions+[1:nRegions])={xdata};
                for row=1:nRegions
                    ydata = app.UITableImported.Data{tminind:tmaxind,row+1}; % row+1 because first row of imported data is time
                    threshold = app.UITableAnalVals.Data.Threshold((ii-1)*nRegions+row);

                    % Find peaks and troughs
                    [peaks,troughs] = app.peakdetect(ydata,threshold,xdata);
                    if (app.DetrendCheckBox.Value) && (size(troughs,1)>1) % detrend the data, use JDR method for now until I better understand Sophie's
                        p = polyfit(troughs(:,1),troughs(:,2),1);
                        ydata=ydata-p(1)*xdata;
                        peaks(:,2)=peaks(:,2)-p(1)*peaks(:,1);
                        troughs(:,2)=troughs(:,2)-p(1)*troughs(:,1);
                    end

                    app.UITableAnalVals.Data.ydata((ii-1)*nRegions+row)={ydata};
                    app.UITableAnalVals.Data.peaks((ii-1)*nRegions+row)={peaks};
                    app.UITableAnalVals.Data.troughs((ii-1)*nRegions+row)={troughs};
                    
                end

            end

            app.updateAnalPlot()

        end

        function updateOutTable(app)
            % app.UITableAnalVals.Data.Ignore
            regions2output = find(app.UITableAnalVals.Data.Ignore==false);
            app.OutTable = table(); % reinitialize each time so we only get the rows we are not ignoring
            app.OutTable.Region = app.UITableAnalVals.Data.Region(regions2output); 
            % | Region | NPulses | Baseline | Peak | Amplitude | Period | Plateau Fraction | ... 
            % Active Area | Average Platwidth | Average Basewidth | Silent Phase | ...
            % Ave Y Value | Notes | TimeWindow | Threshold | Ignore | start time | end time | ...
            % peaks | troughs | trendline | 
            %additional vars are hidden by making columns zero width
            % app.OutTable = app.UITableAnalVals.Data.Region;
            
            tmp=cell2mat(cellfun(@size,app.UITableAnalVals.Data.peaks(regions2output),'UniformOutput',false)); 
            app.OutTable.NPulses(:)=tmp(:,1);

            app.OutTable.Baseline(:)=0;
            app.OutTable.maxPeak(:)=0;
            app.OutTable.peakAmplitude(:)=0;
            app.OutTable.Period(:)=0;
            app.OutTable.Threshold(:)=app.UITableAnalVals.Data.Threshold(regions2output);
            app.OutTable.PlatFraction(:)=0;
            app.OutTable.ActiveArea(:)=0;
            app.OutTable.AvePlatWidth(:)=0;
            app.OutTable.AveBaseWidth(:)=0;
            app.OutTable.SilentPhase(:)=0;
            app.OutTable.AverageYval(:)=0;
            app.OutTable.Notes(:)={''};
            app.UITableOutput.Data = app.OutTable;
            assignin('base','outtable',app.OutTable);
            
        end

        function updateAnalPlot(app)
            
            % colorsAnal = ones([size(app.rows2plot,1) 4]);
            % colors=colororder;
            % for ii=1:size(app.rows2plot,1); colorsAnal(ii,1:3)=colors(1,:); end
            % colorsAnal(app.UITableAnalVals.Data(app.rows2plot,:).Ignore,4)=0.15;

            %h = plot(app.UIAxesAnal,app.UITableAnalVals.Data.xdata(app.rows2plot),app.UITableAnalVals.Data.ydata(app.rows2plot));
            app.rows2plot;
            xs=app.UITableAnalVals.Data.xdata(app.rows2plot);
            ys=app.UITableAnalVals.Data.ydata(app.rows2plot);
            peaks = app.UITableAnalVals.Data.peaks(app.rows2plot);
            troughs = app.UITableAnalVals.Data.troughs(app.rows2plot);

            h = plot(app.UIAxesAnal,xs{1},ys{1});
            if app.UITableAnalVals.Data(app.rows2plot,:).Ignore(1); h.Color(4)=0.15; end % set alpha low for ignored rows
            hold(app.UIAxesAnal,'on')
            if size(app.rows2plot,1)>1 % start with 2 so colors are correct
                for ii=2:size(app.rows2plot,1)
                    h(ii) = plot(app.UIAxesAnal,xs{ii},ys{ii});
                    if app.UITableAnalVals.Data(app.rows2plot,:).Ignore(ii); h(ii).Color(4)=0.15; end % set alpha low for ignored rows
                end
            end
                
            % size(h)
            % for ii=1:size(h,2)
            %     h(ii).Color=colorsAnal(ii,:);
            % end
            
            size(h)
            try
                for ii=1:size(app.rows2plot,1)
                    if app.UITableAnalVals.Data(app.rows2plot,:).Ignore(ii) 
                        scatter(app.UIAxesAnal,peaks{ii}(:,1),peaks{ii}(:,2),'ro',MarkerEdgeAlpha=0.35);
                        scatter(app.UIAxesAnal,troughs{ii}(:,1),troughs{ii}(:,2),'co',MarkerEdgeAlpha=0.35);
                    else
                        scatter(app.UIAxesAnal,peaks{ii}(:,1),peaks{ii}(:,2),'ro');
                        scatter(app.UIAxesAnal,troughs{ii}(:,1),troughs{ii}(:,2),'co');
                    end
                end
            catch
                disp('failed to plot peaks and troughs')
            end
            hold(app.UIAxesAnal,'off')
            size(h)
        end


        % *****************************************************************
        % code / functions from previous version 
        % *****************************************************************
        function [maxtab, mintab]=peakdetect(~, v, delta, x) % add ~ for functions in appdesigner to account for callbacks - JDR
            %PEAKDET Detect peaks in a vector
            %        [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
            %        maxima and minima ("peaks") in the vector V.
            %        MAXTAB and MINTAB consists of two columns. Column 1
            %        contains indices in V, and column 2 the found values.
            %      
            %        With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
            %        in MAXTAB and MINTAB are replaced with the corresponding
            %        X-values.
            %
            %        A point is considered a maximum peak if it has the maximal
            %        value, and was preceded (to the left) by a value lower by
            %        DELTA.
            
            maxtab = [];
            mintab = [];
            
            % v = v(:); % Just in case this wasn't a proper vector
            
            nargin
            if nargin < 3
              x = (1:length(v))';
            else 
              x = x(:);
              if length(v)~= length(x)
                error('Input vectors v and x must have same length');
              end
            end
            
            if (length(delta(:)))>1
              error('Input argument DELTA must be a scalar');
            end
            
            if delta <= 0
              error('Input argument DELTA must be positive');
            end
            
            mn = Inf; mx = -Inf;
            mnpos = NaN; mxpos = NaN;
            
            % delta === "significantly"
            % mx = max to beat
            % mn = min to beat
            lookformax = 1; %look for maximums first
            
            for i=1:length(v)
              this = v(i);
              if this > mx, mx = this; mxpos = x(i); end % if we're (still) going up, assign this as mx (the mx point to beat)
              if this < mn, mn = this; mnpos = x(i); end % if we're (still) going down, assign this as mn (the mn point to beat)
              
              if lookformax %if looking for maximums
                if this < mx-delta  % only if we're significantly lower than mx- we beat it
                  maxtab = [maxtab ; mxpos mx]; %then assign mx as a "maximum"
                  mn = this; mnpos = x(i); %then assign this as mn (the new mn point to beat)
                  lookformax = 0; %now start looking for minimums.
                end  
              else % if looking for minimums
                if this > mn+delta % only if we're significantly higher than mn- we beat it
                 % mintab = [mintab ; mnpos mn]; % then assign mn as a "minimum"
                  mx = this; mxpos = x(i); % then assign this as mx- the new mx point to beat
                  %lookformax = 1; %now start looking for maximums.
                  ii = 1;
                  for ii=1:length(v)
                      if lookformax == 0
                          currentpos = i-ii;
                          current = v(currentpos);
                          if current < v(currentpos+1) %& current > mn+delta   %if we're going down and there isn't a huge jump to mn
                              mn = current; mnpos = x(currentpos);
                          else
                              mintab = [mintab ; mnpos mn];
                              lookformax = 1;
                          end
                      end
                  end
                end
              end
            end
        end





    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: ChooseFiletoImportButton
        function ChooseFiletoImportButtonPushed(app, event)
            [filename, path]=uigetfile({'*.xlsx';'*.xls'}, 'Select Excel File');
            fullfilename = fullfile(path, filename);
            app.FileNameEditField.Value=filename;
            if isequal(filename, 0) % Check if user canceled
                return;
            end        
            % read in the excel data to a table and asign to a UITable
            TableData = readtable(fullfilename,"Sheet",1);
            app.UITableImported.Data = TableData;
            app.UITableImported.ColumnName = TableData.Properties.VariableDescriptions;
            
            % set output file name
            [~,filenamesansext,~] = fileparts(filename);
            app.datapath = path;
            app.SaveFileEditField.Value = [filenamesansext  '_analyzed-'  char(datetime('now','Format','y-MMM-d'))  '.xlsx'];
            
            % initialize the preview window
            app.cols2preview = 2; % start with first data column, update by selecting one or more collumns in the table later
            app.alltimes = app.UITableImported.Data{:,1};
            plot(app.UIAxesSelector,app.alltimes,app.UITableImported.Data{:,app.cols2preview});
            legend(app.UIAxesSelector,app.UITableImported.Data.Properties.VariableDescriptions{app.cols2preview});
            title(app.UIAxesSelector,'Data Preview, select column(s) on left to plot');
            xlabel(app.UIAxesSelector,app.UITableImported.Data.Properties.VariableDescriptions{1});
            ylabel(app.UIAxesSelector,'Val');

            % initialize a time window data table with the first time window add more later by adding rows
            startT = min(app.UITableImported.Data{:,1}); % use min and max instead of 1 and end to handle NaNs
            endT = max(app.UITableImported.Data{:,1});
            app.UITableTimeWindows.Data = table({'Window 1'},{startT},{endT},'VariableNames',{'Time Window','start [s]','end [s]'});

            % initialize UITableAnalVals table
            % Note that this table includes region, time window, threshold and ignore flag. 
            
            app.UITableAnalVals.Data = table(); % init a table and then add variables
            app.UITableAnalVals.Data.Region = app.UITableImported.Data.Properties.VariableDescriptions(2:end)'; % transpose cell array to get iselts as rows
            app.UITableAnalVals.Data.TWindow(:) = app.UITableTimeWindows.Data{1,1};
            app.UITableAnalVals.Data.Threshold(:) = 0.04;
            app.UITableAnalVals.Data.Ignore(:) = false;
            % window, peaks, and troughs. Much is redundant but needed to
            % make plotting multiple rows in live view easier
            app.UITableAnalVals.Data.tmin(:)=0;
            app.UITableAnalVals.Data.tmax(:)=0;
            app.UITableAnalVals.Data.tminind(:)=0;
            app.UITableAnalVals.Data.tmaxind(:)=0;
            app.UITableAnalVals.Data.xdata(:)={0}; % since the rest will take array vals, initialize them as cell arrays
            app.UITableAnalVals.Data.ydata(:)={0};
            app.UITableAnalVals.Data.peaks(:)={0};
            app.UITableAnalVals.Data.troughs(:)={0};
            % app.UITableAnalVals.Data.trendline(:)=0;

            
            % initialize the data analysis plot 
            tmin = cell2mat(app.UITableTimeWindows.Data{1,2});
            tmax = cell2mat(app.UITableTimeWindows.Data{1,3});
            
            tminind = find(app.alltimes>=tmin,1,'first');
            tmaxind = find(app.alltimes<=tmax,1,'last');
            
            ydata = app.UITableImported.Data{tminind:tmaxind,1};
            xdata = app.UITableImported.Data{tminind:tmaxind,2};

            app.updateUITableAnalVals();
            plot(app.UIAxesAnal,ydata,xdata);

            % initialize output table to hold all values to be exported
            % from analysis. These will add variable numbers of columns to
            % the exported table:
            % | Region | NPulses | Baseline | Peak | Amplitude | Period | Plateau Fraction | ... 
            % Active Area | Average Platwidth | Average Basewidth | Silent Phase | ...
            % Ave Y Value | Notes | TimeWindow | Threshold | Ignore | start time | end time | ...
            % peaks | troughs | trendline | 
            %additional vars are hidden by making columns zero width
            % app.OutTable = app.UITableAnalVals.Data.Region;
            % app.OutTable.NPulses(:)=0;
            % app.OutTable.Baseline(:)=0;
            % app.OutTable.Peak(:)=0;
            % app.OutTable.Amplitude(:)=0;
            % app.OutTable.Period(:)=0;
            % app.OutTable.PlatFraction(:)=0;
            % app.OutTable.ActiveArea(:)=0;
            % app.OutTable.AvePlatWidth(:)=0;
            % app.OutTable.AveBaseWidth(:)=0;
            % app.OutTable.SilentPhase(:)=0;
            % app.OutTable.AverageYval(:)=0;
            % app.OutTable.Notes(:)=0;

        end

        % Value changed function: FileNameEditField
        function FileNameEditFieldValueChanged(app, event)
            value = app.FileNameEditField.Value;
            % TODO: change file if filename is updated by typing in this
            % field, currently only displays filename selected via GUI
        end

        % Selection changed function: UITableImported
        function UITableImportedSelectionChanged(app, event)
            app.cols2preview = app.UITableImported.Selection; % since selection mode is set to column, this will only return col values
            app.updatePreviewPlot()
        end

        % Value changed function: SensorTypeDropDown
        function SensorTypeDropDownValueChanged(app, event)
            sensorType = app.SensorTypeDropDown.Value;
            
        end

        % Cell edit callback: UITableTimeWindows
        function UITableTimeWindowsCellEdit(app, event)
            indices = event.Indices;
            newData = event.NewData;
            app.updatePreviewPlot()
            app.updateUITableAnalVals()
            app.updateAnalPlot()
        end

        % Menu selected function: AddRow
        function AddRowMenuSelected(app, event)
            nRegions = size(app.UITableImported.Data,2)-1; % -1 so we don't count the time column
            app.UITableTimeWindows.Data(end+1,2:end)=app.UITableTimeWindows.Data(end,2:end);
            app.UITableTimeWindows.Data(end,1)={' '}; % must assign a string containering cell array so legend works for patches
            app.UITableAnalVals.Data = [app.UITableAnalVals.Data(:,:); app.UITableAnalVals.Data(1:nRegions,:)];
            app.UITableAnalVals.Data.TWindow(end-nRegions:end)={' '};
            % app.updateUITableAnalVals()
        end

        % Menu selected function: DelRow
        function DelRowMenuSelected(app, event)
            if size(app.UITableTimeWindows.Data,1)<2
                disp('can not delete only row')
            else
                % TODO: delete based on time window name
                nRegions = size(app.UITableImported.Data,2)-1; % -1 so we don't count the time column

                app.UITableAnalVals.Data(end-nRegions:end,:)=[];
                app.UITableTimeWindows.Data(end,:)=[];
            end
            % app.updateUITableAnalVals()
        end

        % Cell edit callback: UITableAnalVals
        function UITableAnalValsCellEdit(app, event)
            indices = event.Indices;
            newData = event.NewData;
            app.updateUITableAnalVals();
            app.updateAnalPlot();
        end

        % Selection changed function: UITableAnalVals
        function UITableAnalValsSelectionChanged(app, event)
            AnalValsSelection = app.UITableAnalVals.Selection;
            app.rows2plot = unique(AnalValsSelection(:,1)); % use unique to avoid double counting when multiple cols are selected too
            app.updateAnalPlot();
        end

        % Value changed function: DetrendCheckBox
        function DetrendCheckBoxValueChanged(app, event)
            app.updateUITableAnalVals();
            app.updateAnalPlot();
        end

        % Button pushed function: SaveOutputButton
        function SaveOutputButtonPushed(app, event)
            if isequal(app.SaveFileEditField.Value, 0) % Check if user removed filename
                return;
            end
            fullsavefile = fullfile(app.datapath,app.SaveFileEditField.Value);
            if isfile(fullsavefile)
               % file already exists, ask if you want to overwrite
               msg = "File exists, overwrite? Otherwise cancel and update filename before saving";
               title = "Output file exists";
               selection = uiconfirm(app.UIFigure,msg,title, ...
                    "Options",["Overwrite","Cancel"], ...
                    "DefaultOption",2);
               if selection == "Overwrite"
                   writetable(app.OutTable,fullsavefile)
               end
            else
                writetable(app.OutTable,fullsavefile)
            end
        end

        % Value changed function: SaveFileEditField
        function SaveFileEditFieldValueChanged2(app, event)
            %savefilename = app.SaveFileEditField.Value;
               
        end

        % Button pushed function: UpdateOutputButton
        function UpdateOutputButtonPushed(app, event)
            app.updateOutTable()
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1624 1271];
            app.UIFigure.Name = 'MATLAB App';

            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {'5x', '10x'};
            app.GridLayout.RowHeight = {'fit', '2x', '2x', '8x', '7x'};

            % Create UIAxesSelector
            app.UIAxesSelector = uiaxes(app.GridLayout);
            title(app.UIAxesSelector, 'Title')
            xlabel(app.UIAxesSelector, 'X')
            ylabel(app.UIAxesSelector, 'Y')
            zlabel(app.UIAxesSelector, 'Z')
            app.UIAxesSelector.Layout.Row = [2 3];
            app.UIAxesSelector.Layout.Column = 2;

            % Create ControlsPanel
            app.ControlsPanel = uipanel(app.GridLayout);
            app.ControlsPanel.Title = 'Controls';
            app.ControlsPanel.Layout.Row = 4;
            app.ControlsPanel.Layout.Column = 1;

            % Create SensorTypeDropDownLabel
            app.SensorTypeDropDownLabel = uilabel(app.ControlsPanel);
            app.SensorTypeDropDownLabel.HorizontalAlignment = 'right';
            app.SensorTypeDropDownLabel.Position = [13 454 71 22];
            app.SensorTypeDropDownLabel.Text = 'Sensor Type';

            % Create SensorTypeDropDown
            app.SensorTypeDropDown = uidropdown(app.ControlsPanel);
            app.SensorTypeDropDown.Items = {'None', 'Fura', 'Perceval', 'Laconic'};
            app.SensorTypeDropDown.ValueChangedFcn = createCallbackFcn(app, @SensorTypeDropDownValueChanged, true);
            app.SensorTypeDropDown.Position = [99 454 100 22];
            app.SensorTypeDropDown.Value = 'None';

            % Create Instructions
            app.Instructions = uilabel(app.ControlsPanel);
            app.Instructions.Position = [251 387 280 89];
            app.Instructions.Text = {'Instructions: '; '1. Update time window start and end above'; '2. Right click table to add or remove time windows'; '3. Select row(s) below to plot, update params'; '4. Click Update Ouput to refresh output table'; '5. Save to output file'};

            % Create UITableAnalVals
            app.UITableAnalVals = uitable(app.ControlsPanel);
            app.UITableAnalVals.ColumnName = {'Region'; 'TWindow'; 'Threshold'; 'Ignore'; 'tmin'; 'tmax'};
            app.UITableAnalVals.ColumnWidth = {'auto', 'auto', 'auto', 'auto', 0, 0};
            app.UITableAnalVals.RowName = {};
            app.UITableAnalVals.ColumnEditable = [false false true true false false];
            app.UITableAnalVals.CellEditCallback = createCallbackFcn(app, @UITableAnalValsCellEdit, true);
            app.UITableAnalVals.SelectionChangedFcn = createCallbackFcn(app, @UITableAnalValsSelectionChanged, true);
            app.UITableAnalVals.Position = [1 51 530 325];

            % Create DetrendCheckBox
            app.DetrendCheckBox = uicheckbox(app.ControlsPanel);
            app.DetrendCheckBox.ValueChangedFcn = createCallbackFcn(app, @DetrendCheckBoxValueChanged, true);
            app.DetrendCheckBox.Text = 'Detrend';
            app.DetrendCheckBox.Position = [19 433 65 22];
            app.DetrendCheckBox.Value = true;

            % Create UpdateOutputButton
            app.UpdateOutputButton = uibutton(app.ControlsPanel, 'push');
            app.UpdateOutputButton.ButtonPushedFcn = createCallbackFcn(app, @UpdateOutputButtonPushed, true);
            app.UpdateOutputButton.BackgroundColor = [0.0667 0.4431 0.7451];
            app.UpdateOutputButton.Position = [5 15 259 23];
            app.UpdateOutputButton.Text = 'Update Output with Current Analysis Params';

            % Create TabGroup
            app.TabGroup = uitabgroup(app.GridLayout);
            app.TabGroup.Layout.Row = [4 5];
            app.TabGroup.Layout.Column = 2;

            % Create PlotsTab
            app.PlotsTab = uitab(app.TabGroup);
            app.PlotsTab.Title = 'Plots';

            % Create UIAxesAnal
            app.UIAxesAnal = uiaxes(app.PlotsTab);
            title(app.UIAxesAnal, 'Title')
            xlabel(app.UIAxesAnal, 'X')
            ylabel(app.UIAxesAnal, 'Y')
            zlabel(app.UIAxesAnal, 'Z')
            app.UIAxesAnal.Position = [0 1 1061 923];

            % Create WaveletviewTab
            app.WaveletviewTab = uitab(app.TabGroup);
            app.WaveletviewTab.Title = 'Wavelet view';

            % Create UIAxesWVLT
            app.UIAxesWVLT = uiaxes(app.WaveletviewTab);
            title(app.UIAxesWVLT, 'Wavelet Analysis')
            xlabel(app.UIAxesWVLT, 'X')
            ylabel(app.UIAxesWVLT, 'Y')
            zlabel(app.UIAxesWVLT, 'Z')
            app.UIAxesWVLT.Position = [1 40 793 685];

            % Create UIAxes3
            app.UIAxes3 = uiaxes(app.WaveletviewTab);
            app.UIAxes3.Position = [1 724 793 200];

            % Create UIAxes4
            app.UIAxes4 = uiaxes(app.WaveletviewTab);
            app.UIAxes4.Position = [793 40 269 685];

            % Create DataSummaryPanel
            app.DataSummaryPanel = uipanel(app.GridLayout);
            app.DataSummaryPanel.Title = 'Data Summary';
            app.DataSummaryPanel.Layout.Row = 5;
            app.DataSummaryPanel.Layout.Column = 1;

            % Create SaveOutputButton
            app.SaveOutputButton = uibutton(app.DataSummaryPanel, 'push');
            app.SaveOutputButton.ButtonPushedFcn = createCallbackFcn(app, @SaveOutputButtonPushed, true);
            app.SaveOutputButton.BackgroundColor = [0.0667 0.4431 0.7451];
            app.SaveOutputButton.Position = [13 18 100 23];
            app.SaveOutputButton.Text = 'Save Output';

            % Create SaveFileEditField
            app.SaveFileEditField = uieditfield(app.DataSummaryPanel, 'text');
            app.SaveFileEditField.ValueChangedFcn = createCallbackFcn(app, @SaveFileEditFieldValueChanged2, true);
            app.SaveFileEditField.Position = [124 18 397 22];

            % Create UITableOutput
            app.UITableOutput = uitable(app.DataSummaryPanel);
            app.UITableOutput.ColumnName = {'Region'; 'nPeaks'; 'Baseline'; 'maxPeak'; 'PeakAmplitude'; 'Period'; 'Threshold'; 'PlatFrac'; 'ActiveArea'; 'AvePlatWidth'; 'AveBaseWidth'; 'SlientPhase'; 'AveYval'; 'Notes'};
            app.UITableOutput.RowName = {};
            app.UITableOutput.ColumnEditable = [false false false false false false false false false false false false false true];
            app.UITableOutput.Position = [5 57 516 248];

            % Create UITableTimeWindows
            app.UITableTimeWindows = uitable(app.GridLayout);
            app.UITableTimeWindows.ColumnName = {'Time Window'; 'Start [s]'; 'End [s]'};
            app.UITableTimeWindows.RowName = {};
            app.UITableTimeWindows.ColumnEditable = true;
            app.UITableTimeWindows.RowStriping = 'off';
            app.UITableTimeWindows.CellEditCallback = createCallbackFcn(app, @UITableTimeWindowsCellEdit, true);
            app.UITableTimeWindows.Layout.Row = 3;
            app.UITableTimeWindows.Layout.Column = 1;

            % Create UITableImported
            app.UITableImported = uitable(app.GridLayout);
            app.UITableImported.ColumnName = {'Column 1'; 'Column 2'; 'Column 3'; 'Column 4'};
            app.UITableImported.RowName = {};
            app.UITableImported.SelectionType = 'column';
            app.UITableImported.SelectionChangedFcn = createCallbackFcn(app, @UITableImportedSelectionChanged, true);
            app.UITableImported.Layout.Row = 2;
            app.UITableImported.Layout.Column = 1;

            % Create FileNameEditField
            app.FileNameEditField = uieditfield(app.GridLayout, 'text');
            app.FileNameEditField.ValueChangedFcn = createCallbackFcn(app, @FileNameEditFieldValueChanged, true);
            app.FileNameEditField.Editable = 'off';
            app.FileNameEditField.Layout.Row = 1;
            app.FileNameEditField.Layout.Column = 2;

            % Create ChooseFiletoImportButton
            app.ChooseFiletoImportButton = uibutton(app.GridLayout, 'push');
            app.ChooseFiletoImportButton.ButtonPushedFcn = createCallbackFcn(app, @ChooseFiletoImportButtonPushed, true);
            app.ChooseFiletoImportButton.BackgroundColor = [0.0667 0.4431 0.7451];
            app.ChooseFiletoImportButton.Layout.Row = 1;
            app.ChooseFiletoImportButton.Layout.Column = 1;
            app.ChooseFiletoImportButton.Text = 'Choose File to Import';

            % Create ContextMenu
            app.ContextMenu = uicontextmenu(app.UIFigure);

            % Create AddRow
            app.AddRow = uimenu(app.ContextMenu);
            app.AddRow.MenuSelectedFcn = createCallbackFcn(app, @AddRowMenuSelected, true);
            app.AddRow.Text = 'AddRow';

            % Create DelRow
            app.DelRow = uimenu(app.ContextMenu);
            app.DelRow.MenuSelectedFcn = createCallbackFcn(app, @DelRowMenuSelected, true);
            app.DelRow.Text = 'DelRow';
            
            % Assign app.ContextMenu
            app.UITableTimeWindows.ContextMenu = app.ContextMenu;

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = MerrinsLabOscillationsApp_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end