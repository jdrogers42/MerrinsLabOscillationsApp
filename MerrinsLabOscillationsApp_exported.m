classdef MerrinsLabOscillationsApp_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                  matlab.ui.Figure
        GridLayout                matlab.ui.container.GridLayout
        LoadstateButton           matlab.ui.control.Button
        SaveFileEditField         matlab.ui.control.EditField
        SaveOutputButton          matlab.ui.control.Button
        ChooseFiletoImportButton  matlab.ui.control.Button
        FileNameEditField         matlab.ui.control.EditField
        UITableImported           matlab.ui.control.Table
        UITableTimeWindows        matlab.ui.control.Table
        DataSummaryPanel          matlab.ui.container.Panel
        UITableOutput             matlab.ui.control.Table
        TabGroup                  matlab.ui.container.TabGroup
        PlotsTab                  matlab.ui.container.Tab
        UIAxesAnal                matlab.ui.control.UIAxes
        WaveletviewTab            matlab.ui.container.Tab
        UIAxesWVLTy               matlab.ui.control.UIAxes
        UIAxesWVLTx               matlab.ui.control.UIAxes
        UIAxesWVLT                matlab.ui.control.UIAxes
        ControlsPanel             matlab.ui.container.Panel
        UpdateOutputButton        matlab.ui.control.Button
        ShowPlatsCheckBox         matlab.ui.control.CheckBox
        SheetDropDown             matlab.ui.control.DropDown
        SheetLabel                matlab.ui.control.Label
        UsefindpeaksCheckBox      matlab.ui.control.CheckBox
        IgnoreOutliersCheckBox    matlab.ui.control.CheckBox
        DetrendCheckBox           matlab.ui.control.CheckBox
        Instructions              matlab.ui.control.Label
        SensorTypeDropDown        matlab.ui.control.DropDown
        SensorTypeLabel           matlab.ui.control.Label
        UITableAnalVals           matlab.ui.control.Table
        UIAxesSelector            matlab.ui.control.UIAxes
        ContextMenu               matlab.ui.container.ContextMenu
        AddRow                    matlab.ui.container.Menu
        DelRow                    matlab.ui.container.Menu
    end

    
    % app written by J.D. Rogers <rogersjd@gmail.com> on Oct 11, 2025
    % Based on Merrins lab code, converted to GUI
    % last updated: 20251216


    % shared variables across functions
    properties (Access = private)
        fullfilename; % file selected to import
        sheet; % which sheet of the selected file to work with
        DialogChooseSheet;
        cols2preview; % which columns in the imported table to plot in the preview
        rows2plot=[1]; % which row in the AnalVal table to plot in the anal window 
        alltimes;  % time for the entire data set
        timewindows=table(); % table containing the start and end times for each time window
        datapath; % path where data is located and output is saved
        OutTable = table(); %store all the output info in this table
        TableData; % Table to hold values directly read from file
        artifacts; % 2 column vector of [locations peaks] of artifacts in the data to be ignored
        timecol; % the column containing the time values (not always the first col)
        timeunits; 
        firstdatacol; % first column with data, assume it starts with #
    end


    methods (Access = private)
        
        function initializeTables(app)
            % read in the excel data to a table and asign to a UITable
            app.TableData = readtable(app.fullfilename,"Sheet",app.SheetDropDown.Value);
            app.UITableImported.Data = app.TableData;
            app.UITableImported.ColumnName = app.TableData.Properties.VariableDescriptions;
            % identify columns that are data
            app.timecol = find(contains(app.TableData.Properties.VariableDescriptions,'Time'));
            timeunitsstartindex = strfind(app.TableData.Properties.VariableDescriptions{app.timecol},'[');
            timeunitsendindex = strfind(app.TableData.Properties.VariableDescriptions{app.timecol},']');
            app.timeunits = app.TableData.Properties.VariableDescriptions{app.timecol}(timeunitsstartindex:timeunitsendindex);
            app.firstdatacol = app.timecol+1; %find(contains(app.TableData.Properties.VariableDescriptions,'#'),1);

            % Identify artifacts. Assume that artifacts are present
            % across all regions, so use the mean of all regions. This
            % greatly reduces the standard deviation and smoothes the
            % curves, but not the artifact making it more prominent. Then
            % use findpeaks() to identify artifacts that are less than 2.5
            % samples in width and have a prominance more than 2x the
            % standard deviation of the average signal. Set the values to
            % NaN but plot them in the preview window so we know what was
            % removed.

            meandata = mean(app.TableData{:,app.firstdatacol:end}');
            [pks,locs]=findpeaks(meandata,'MinPeakProminence',2*std(meandata),'MaxPeakWidth',2.5);
            [pksinv,locsinv]=findpeaks(-meandata,'MinPeakProminence',2*std(meandata),'MaxPeakWidth',2.5);
            % assignin("base","pks",pks);assignin("base","pksinv",pksinv)
            pks = [pks,pksinv]; locs = [locs,locsinv];
            vals = app.TableData{locs,app.firstdatacol:end}; % values of each region atthe artifact location
            app.artifacts = [locs' vals]; % artifacts to ignore if box is checke % assignin("base","artifacts",app.artifacts)
            % initialize the preview window
            app.cols2preview = app.firstdatacol; % start with first data column, update by selecting one or more collumns in the table later
            % app.alltimes = app.UITableImported.Data{:,1};
            app.alltimes = app.UITableImported.Data{:,app.timecol};
            % assignin("base","tabledata",app.TableData)
            plot(app.UIAxesSelector,app.alltimes,app.UITableImported.Data{:,app.cols2preview});
            if app.artifacts % only plot the artifacts if they exist
            hold(app.UIAxesSelector,"on");
                plot(app.UIAxesSelector,app.UITableImported.Data{app.artifacts(:,1),1}, ...
                                    mean(app.artifacts(:,app.firstdatacol:end),2), ...
                                    "rx",'DisplayName','Artifact');
                hold(app.UIAxesSelector,"off");
            end
            legend(app.UIAxesSelector,app.UITableImported.Data.Properties.VariableDescriptions{app.cols2preview},Location="east");
            title(app.UIAxesSelector,'Data Preview, select column(s) on left to plot');
            xlabel(app.UIAxesSelector,app.UITableImported.Data.Properties.VariableDescriptions{1});
            ylabel(app.UIAxesSelector,'Val');

            % initialize a time window data table with the first time window add more later by adding rows
            startT = min(app.UITableImported.Data{:,app.timecol}); % use min and max instead of 1 and end to handle NaNs
            endT = max(app.UITableImported.Data{:,app.timecol});
            app.UITableTimeWindows.Data = table({'Window 1'},{startT},{endT},'VariableNames',{'Time Window',['start ' app.timeunits],['end ' app.timeunits]});

            % initialize UITableAnalVals table
            % Note that this table includes region, time window, threshold and ignore flag. 
            
            app.UITableAnalVals.Data = table(); % init a table and then add variables
            app.UITableAnalVals.Data.Region = app.UITableImported.Data.Properties.VariableDescriptions(app.firstdatacol:end)'; % transpose cell array to get iselts as rows
            app.UITableAnalVals.Data.TWindow(:) = app.UITableTimeWindows.Data{1,1};
            app.UITableAnalVals.Data.Threshold(:) =  std(app.UITableImported.Data{:,app.firstdatacol:end}); % 0.04; % old code used 0.04 as default threshhold for peak finding, trying std() as a starting point instead
            app.UITableAnalVals.Data.Ignore(:) = false;
            % window, peaks, and troughs. Much is redundant but needed to
            % make plotting multiple rows in live view easier
            app.UITableAnalVals.Data.tmin(:)=0;
            app.UITableAnalVals.Data.tmax(:)=0;
            app.UITableAnalVals.Data.tminind(:)=0;
            app.UITableAnalVals.Data.tmaxind(:)=0;
            % since the rest will take array vals, initialize them as cell arrays
            app.UITableAnalVals.Data.xdata(:)={0}; 
            app.UITableAnalVals.Data.ydata(:)={0};
            app.UITableAnalVals.Data.peaks(:)={0};
            app.UITableAnalVals.Data.troughs(:)={0};

            
            % initialize the data analysis plot 
            tmin = cell2mat(app.UITableTimeWindows.Data{1,2});
            tmax = cell2mat(app.UITableTimeWindows.Data{1,3});
            
            tminind = find(app.alltimes>=tmin,1,'first');
            tmaxind = find(app.alltimes<=tmax,1,'last');
            
            xdata = app.UITableImported.Data{tminind:tmaxind,app.timecol};
            ydata = app.UITableImported.Data{tminind:tmaxind,app.firstdatacol};

            app.updateUITableAnalVals();
            plot(app.UIAxesAnal,xdata,ydata);
        end

        function updatePreviewPlot(app)
            % replot the full time sequence for the given selected columns
            plot(app.UIAxesSelector,app.UITableImported.Data{:,app.timecol},app.UITableImported.Data{:,app.cols2preview});
            legend(app.UIAxesSelector,app.UITableImported.Data.Properties.VariableDescriptions{app.cols2preview},Location="east");
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
            %if app.IgnoreOutliersCheckBox.Value; plot(app.UIAxesSelector,app.artifacts(:,1),mean(app.artifacts(:,2:end),2),"rx",'DisplayName','Artifact');end
            plot(app.UIAxesSelector,app.UITableImported.Data{app.artifacts(:,1),app.timecol},mean(app.artifacts(:,2:end),2),"rx",'DisplayName','Artifact')
            hold(app.UIAxesSelector,"off");
        end

        function updateUITableAnalVals(app)
                        %TableAnalVals =  app.UITableAnalVals.Data;
            nTWindows = size(app.UITableTimeWindows.Data,1);
            % For each time window row, make a set of nRegions of rows to
            % the analysis table
            nRegions = size(app.UITableImported.Data,2)+1-app.firstdatacol; % only count the number of data cols
            for ii=1:nTWindows
                % Add nRegions of rows for each time window to the tale
                app.UITableAnalVals.Data.TWindow((ii-1)*nRegions+(1:nRegions))=app.UITableTimeWindows.Data{ii,1};
                % set the start and end time (and indicesfor each table row
                tmin = app.UITableTimeWindows.Data{ii,2}{1}; % index {1} to return value instead of cell array
                tmax = app.UITableTimeWindows.Data{ii,3}{1};
                app.UITableAnalVals.Data.tmin((ii-1)*nRegions+(1:nRegions))=tmin; % store these vals in each row
                app.UITableAnalVals.Data.tmax((ii-1)*nRegions+[1:nRegions])=tmax;
                tminind = find(app.alltimes>=tmin,1,'first');
                tmaxind = find(app.alltimes<=tmax,1,'last');
                app.UITableAnalVals.Data.tminind((ii-1)*nRegions+[1:nRegions])= tminind; % store these vals in each row
                app.UITableAnalVals.Data.tmaxind((ii-1)*nRegions+[1:nRegions])= tmaxind;
                
                if tminind<tmaxind % if the window start is after the end, just chose the end value
                    xdata=app.UITableImported.Data{tminind:tmaxind,app.timecol};
                else
                    xdata=app.UITableImported.Data{tmaxind:end,app.timecol};
                end
                xdata=xdata-xdata(1); 
                
                % wavelet filter bank
                if contains(app.timeunits, '[s]')
                    fb = cwtfilterbank('wavelet','morse','SignalLength',length(xdata),...
                        'WaveletParameters',[3 60],'VoicesPerOctave',10,'SamplingPeriod',seconds((xdata(2)-xdata(1))), 'PeriodLimits',[seconds(10) seconds(600)]);
                elseif contains(app.timeunits, '[min]')
                    fb = cwtfilterbank('wavelet','morse','SignalLength',length(xdata),...
                        'WaveletParameters',[3 60],'VoicesPerOctave',10,'SamplingPeriod',minutes((xdata(2)-xdata(1))), 'PeriodLimits',[minutes(10/60) minutes(600/60)]);
                end                    
                    
                app.UITableAnalVals.Data.xdata((ii-1)*nRegions+[1:nRegions])={xdata};
                for row=1:nRegions
                    ydata = app.UITableImported.Data{tminind:tmaxind,row-1+app.firstdatacol}; % skip first non-data rows
                    threshold = app.UITableAnalVals.Data.Threshold((ii-1)*nRegions+row);
                    
                    % Find peaks and troughs
                    if (app.UsefindpeaksCheckBox.Value)
                        [pks,locs,widths,proms] = findpeaks(ydata,'MinPeakProminence',threshold,'MinPeakWidth',4);
                        peaks = [xdata(locs) pks];
                        [pks,locs,widths,proms] = findpeaks(-ydata,'MinPeakProminence',threshold,'MinPeakWidth',4);
                        troughs = [xdata(locs) -pks];
                    else
                        [peaks,troughs] = app.peakdetect(ydata,threshold,xdata);
                        
                        TroughSize = size(troughs,1); %grab the first value, the number of troughs.
                        PeakSize = size(peaks,1);

                        if(PeakSize > 1) % changed to >1 to ensure we have room to subtract the end peaks below -- JDR
                            if (peaks(1,1) < troughs(1,1))
                                peaks = peaks(2:PeakSize,:); %remove first peak point if it is before the first trough point
                                PeakSize = PeakSize - 1; %adjust peakSize by -1
                            end
                            if (peaks(PeakSize,1) > troughs(TroughSize,1)) %remove last peak point if it is after the last trough point
                                peaks = peaks(1:end-1,:);
                                PeakSize = PeakSize - 1;
                            end
                        end
                    end
                    
                    if (app.DetrendCheckBox.Value) && (size(troughs,1)>1) % detrend the data
                        p = polyfit(troughs(:,1),troughs(:,2),1);
                        ydata=ydata-p(1)*xdata;
                        peaks(:,2)=peaks(:,2)-p(1)*peaks(:,1);
                        troughs(:,2)=troughs(:,2)-p(1)*troughs(:,1);
                    end

                    % Find plateau area, s

                    app.UITableAnalVals.Data.ydata((ii-1)*nRegions+row)={ydata};
                    app.UITableAnalVals.Data.peaks((ii-1)*nRegions+row)={peaks};
                    app.UITableAnalVals.Data.troughs((ii-1)*nRegions+row)={troughs};

                    % calculate average wavelet over time
                    sig=fillmissing(ydata,'nearest');
                    [cfs,frq,coi] = wt(fb,sig);
                    frq = 1./seconds(frq);
                    coi = 1./seconds(coi);
                    p   = 1./ frq;
                    % assignin("base","mcfs",mean(abs(cfs)'))
                    % assignin("base","frq",frq)
                    app.UITableAnalVals.Data.aveWavelet((ii-1)*nRegions+row)={[frq mean(abs(cfs)')']};


                end

            end

            app.updateAnalPlot();

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
            
            pks = app.UITableAnalVals.Data.peaks(regions2output);
            trs = app.UITableAnalVals.Data.troughs(regions2output);
           

            tmp=cell2mat(cellfun(@size,pks,'UniformOutput',false)); 
            app.OutTable.NPulses(:)=tmp(:,1);

            tmp=cell2mat(cellfun(@mean,trs,'UniformOutput',false));
            app.OutTable.Baseline(:)=tmp(:,2);

            tmp=cell2mat(cellfun(@mean,pks,'UniformOutput',false));
            app.OutTable.avgPeak(:)=tmp(:,2);

            app.OutTable.peakAmplitude(:)=app.OutTable.avgPeak(:)-app.OutTable.Baseline(:);

            tmp=cell2mat(cellfun(@size,pks,'UniformOutput',false));
            app.OutTable.Period(:)=cellfun(@(x) x(end,1)-x(1,1),trs)./tmp(:,1);

            app.OutTable.Threshold(:)=app.UITableAnalVals.Data.Threshold(regions2output);

            % ideal to identify key points (peak, trough, rising edge,
            % trailing edge. Then the fractions and rise/fall time can be
            % calculated for different sensors from these vals. 

            app.OutTable.PlatFraction(:)=0;
            app.OutTable.ActiveArea(:)=0;
            app.OutTable.AvePlatWidth(:)=0;
            app.OutTable.AveBaseWidth(:)=0;
            app.OutTable.SilentPhase(:)=0;
            app.OutTable.AverageYval(:)=0;
            app.OutTable.Notes(:)={''};
            
            
            % app.SensorTypeDropDown.Value
            % if app.SensorTypeDropDown.Value=='Calcium'
                xdata = app.UITableAnalVals.Data.xdata(regions2output);
                ydata = app.UITableAnalVals.Data.ydata(regions2output);
                threshold = app.UITableAnalVals.Data.Threshold(regions2output);
                for row=1:size(app.OutTable,1)
                    [avgpf avgpa avgpw avgbw avgsp] = app.platfunction([xdata{row}, ydata{row}], 1, pks{row}, trs{row}, threshold(row));
                    app.OutTable.PlatFraction(row)=avgpf;
                    app.OutTable.ActiveArea(row)=avgpa;
                    app.OutTable.AvePlatWidth(row)=avgpw;
                    app.OutTable.AveBaseWidth(row)=avgbw;
                    app.OutTable.SilentPhase(row)=avgsp;
                end
            % end
            app.UITableOutput.Data = app.OutTable;

            % save vars into worspace for debugging purposes
            % assignin('base','trs',trs)
            % assignin('base','outtable',app.OutTable);
            
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
            hold(app.UIAxesAnal,'on');
            if size(app.rows2plot,1)>1 % start with 2 so colors are correct
                for ii=2:size(app.rows2plot,1)
                    h(ii) = plot(app.UIAxesAnal,xs{ii},ys{ii});
                    if app.UITableAnalVals.Data(app.rows2plot,:).Ignore(ii); h(ii).Color(4)=0.15; end % set alpha low for ignored rows
                end
            end
            
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
            hold(app.UIAxesAnal,'off');
            if app.ShowPlatsCheckBox.Value
                [avgpf avgpa avgpw avgbw avgsp] = app.platfunction([xs{1}, ys{1}], 1, peaks{1}, troughs{1}, 0); % this will recalculate the plateaus, but also plot them in the anal window for the first selected row only
            end
        end


        function updateWaveletPlot(app)
            % Create Time-Frequency Representations, and seperate the 10-60s signal and 60-600s signal and the cwt coeffs
            % written by Huixia Ren, 2025/10/07
            % Updated by J.D. Rogers 2025/10/31

            % must handle NaN values before sending to filterbank, for now
            % set NaN to nearest val
            t=app.UITableAnalVals.Data.xdata{app.rows2plot(1)};
            sig=fillmissing(app.UITableAnalVals.Data.ydata{app.rows2plot(1)},'nearest');
            
            % cwt parameters
            % we used the analytic Morse(3,60) wavelet with period limits from 10s to 600s as the continuous wavelet
            % transform filter bank. The number of voices per octave was 10. So the
            % continuous wavelet transform return a matrix of 60 x T.  The mean Row
            % 1-26 27-60
            
            if contains(app.timeunits, '[s]')
                    fb = cwtfilterbank('wavelet','morse','SignalLength',length(sig),...
                        'WaveletParameters',[3 60],'VoicesPerOctave',10,'SamplingPeriod',seconds((t(2)-t(1))), 'PeriodLimits',[seconds(10) seconds(600)]);
                    
                    [cfs,frq,coi] = wt(fb,sig);
                    frq = 1./seconds(frq);
                    coi = 1./seconds(coi);
                    p   = 1./ frq;
            
            elseif contains(app.timeunits, '[min]')
                    fb = cwtfilterbank('wavelet','morse','SignalLength',length(sig),...
                        'WaveletParameters',[3 60],'VoicesPerOctave',10,'SamplingPeriod',minutes((t(2)-t(1))), 'PeriodLimits',[minutes(10/60) minutes(600/60)]);
            
                    [cfs,frq,coi] = wt(fb,sig);
                    frq = 1./minutes(frq);
                    coi = 1./minutes(coi);
                    p   = 1./ frq;
            
            end

            
            
            
            % for i=1:size(cfs,2)
            %    cfs( frq < coi(i), i) = 0;
            % end

            plot(app.UIAxesWVLTx, t,sig);
            %ylabel('Ca^2^+')

            plot(app.UIAxesWVLTy,mean(abs(cfs(:,:))'),frq);
            set(app.UIAxesWVLTy,'yscale','log');
            %xlabel('cwt')

            pcolor(app.UIAxesWVLT,t,frq,abs(cfs),LineStyle='none');
            % shading interp;
            %axis tight;
            set(app.UIAxesWVLT,'yscale','log');
            %xlabel('Time (s)');ylabel('Frequency (Hz)')
            linkaxes([app.UIAxesWVLT app.UIAxesWVLTx],'x');
            linkaxes([app.UIAxesWVLT app.UIAxesWVLTy],'y');
            
        end

        % *****************************************************************
        % code / functions from previous version 
        % *****************************************************************
        function [maxtab, mintab]=peakdetect(app, v, delta, x) % add ~ for functions in appdesigner to account for callbacks - JDR
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
            
            nargin;
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


        function [pArea] = platarea(~,baseheight,platheight,troughs,peaks,curve) % added ~ argument for callback -JDR
            
            pArea = [];
            pSize = size(peaks);
            pSize = pSize(1);
            
            for pulse=1:pSize
                tmp = abs(curve(:,1)-troughs(pulse,1));
                [ida ida] = min(tmp);
                tmp = abs(curve(:,1)-troughs((pulse+1),1));
                [idc idc] = min(tmp);
                tmp = abs(curve(:,1)-peaks(pulse,1));
                [idb idb] = min(tmp);
            
                y1 = baseheight(pulse);
                y2 = [];
            
                %find the index (x value) where pulse begins- where curve dips below
                %baseheight
                for i=idb:-1:ida
                    if (curve(i,2) < y1)
                        break
                    end
                end
                lbound = i+1;
            
                %find the index (x value) where pulse ends- where curve dips below
                %baseheight
                for i=idb:idc
                    if (curve(i,2) < y1)
                        break
                    end
                end
                rbound = i-1;
            
                xscale = curve(lbound:rbound,1); % the x values of the pulse
                y2 = curve(lbound:rbound,2); % the y values of the pulse
            
                %instead of making y values all the way up to peak, plateau it off when
                %y values are greater than platheight
                for i=1:length(y2) 
                    if (y2(i) > platheight(pulse))
                        y2(i) = platheight(pulse);
                    end
                end
            
                % area(xscale,y2,y1);
                % colormap cool;
                
                if ~isempty(xscale)
                    yv = [baseheight(pulse); y2; baseheight(pulse); baseheight(pulse)];
                    xv = [xscale(1); xscale; xscale(end); xscale(1)];
                    pArea(pulse) = polyarea(xv,yv);
                end
            end
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [avgpf avgpa avgpw avgbw avgsp] = platfunction(app, curve, curveNum, peaks, troughs, threshold)
        %      Determines left and right edges of plateau based off of user input.
        %      Graphically represents plateau and base, and calculates plateau
        %      width over the base width (aka the Plateau Fraction)
        %
        % Modified to take only one curve
        doplots = app.ShowPlatsCheckBox.Value;

        %Call global variable used to hold spreadsheet
        
        %Initialize arrays for later use
        basewidth = [];
        platwidth = [];
        baseheight = [];
        amplitudes = [];
        minIndicies =[];
        maxIndicies = [];
        left_edge = [];
        right_edge = [];
        expoData = [];
        pArea = [];
        %Plateau Fractions
        pf = [];
        
        minSize = length(troughs(:,1));
        maxSize = length(peaks(:,1));
        
        %Determine parameters before messing with data
        avgBaseline = mean(troughs(:,2));
        avgPeak = mean(peaks(:,2));
        deltaR = avgPeak - avgBaseline;
        
        PeakSize = size(peaks);
        PeakSize = PeakSize(1);
        period = (troughs(end,1) - troughs(1,1))/PeakSize; % calculate period!!!
        
        %period = (curve(end,1) - curve(1,1))/length(peaks(:,2));
        
        
        %Determine average height between adjacent bases
        for i = 2:minSize
            basewidth(i-1) = (troughs(i,1)-troughs((i-1),1));
            baseheight(i-1) = (troughs(i,2)+troughs((i-1),2))/2;
        end
        
        %Prepare initial parameters for export
        % expoData1 = [curveNum maxSize avgBaseline(1) avgPeak(1) deltaR(1) period(1) threshold];
        
        %Determine amplitude for each pulse using base midpoint
        amplitudes = (peaks(:,2) - baseheight(:));
        amplitudes = transpose(amplitudes);
        
        %Ask user for plateau definition
        % ampperc = input('At what percent of the amplitude would you like \nto analyze the plateau fraction? (default 50) ');
        % if isempty(ampperc)
        %     ampperc = 50;
        % elseif (ampperc > 95)
        %     ampperc = 95;
        % elseif (ampperc < 5)
        %     ampperc = 5;
        % end
        ampperc=50;
        fprintf('Analyzing at %d%% of the amplitude ...\n',ampperc);
        ampperc = (ampperc/100);
        
        %Create array of plateau heights
        platheight = (ampperc.*amplitudes)+baseheight;
        
        %Determine length of curve
        curveSize = length(curve(:,1));
        
        %Scan curve for min and peak indicies
        curMin = 1;
        curMax = 1;
        for i=1:curveSize
            if (curMin <= minSize)
                if (curve(i,1) == troughs(curMin,1))
                    minIndicies = [minIndicies i];
                    curMin = curMin + 1;
                end
            end
            if (curMax <= maxSize)
                if (curve(i,1) == peaks(curMax,1))
                    maxIndicies = [maxIndicies i];
                    curMax = curMax + 1;
                end
            end
        end
        
        if doplots
            % f2 = figure(2);
            % clf;
            % set(f2, 'Position', [660 50 600 370]);
            % assignin("base","curve",curve)
            plot(app.UIAxesAnal,curve(:,1),curve(:,2));
            hold(app.UIAxesAnal,"on")
            % title(['Perceval curve ', num2str(curveNum)]);
            % xlabel('Time (min)');
            % ylabel(strcat('Perceval Ratio'));
            % grid on;
            
            % hold on;
        end

        pArea = app.platarea(baseheight,platheight,troughs,peaks,curve);
        
        if doplots
            % hold(app.UIAxesAnal,'on')
            plot(app.UIAxesAnal,peaks(:,1),peaks(:,2), 'Color', 'r', 'Marker', '*', 'LineStyle', 'none');
            plot(app.UIAxesAnal,troughs(:,1),troughs(:,2), 'Color', 'g', 'Marker', '*', 'LineStyle', 'none');
        end

        for pulse=1:maxSize
            inc_high = minIndicies(pulse);
            for i=minIndicies(pulse):maxIndicies(pulse);
                inc_low = inc_high;
                low_t = i - 1;
                inc_high = curve(i,2);
                if (inc_high >= platheight(pulse))
                    break
                end
            end
            
            %calculate slope
            m = (inc_high - inc_low)*10;
            %calculate difference between lower bound and platheight
            p_diff = platheight(pulse) - inc_low;
            %calculate difference between lower bound time and platheight time
            t_diff = p_diff/m;
            %estimate location of edge
            left_edge(pulse) = curve(low_t,1) + t_diff;
            if doplots
                hb = line(app.UIAxesAnal,[troughs(pulse,1) troughs((pulse+1),1)], [baseheight(pulse) baseheight(pulse)], 'Color', 'g','LineWidth',2);
                plot(app.UIAxesAnal,troughs(pulse,1),baseheight(pulse), 'Color', 'g', 'Marker', 'O', 'LineStyle', 'none');
                plot(app.UIAxesAnal,troughs(pulse,1),baseheight(pulse), 'Color', 'g', 'Marker', '+', 'LineStyle', 'none');
                plot(app.UIAxesAnal,troughs((pulse+1),1),baseheight(pulse), 'Color', 'g', 'Marker', 'O', 'LineStyle', 'none');
                plot(app.UIAxesAnal,troughs((pulse+1),1),baseheight(pulse), 'Color', 'g', 'Marker', '+', 'LineStyle', 'none');
                line(app.UIAxesAnal,[troughs(pulse,1) troughs(pulse,1)], [troughs(pulse,2) baseheight(pulse)], 'Color', 'c');
                line(app.UIAxesAnal,[troughs((pulse+1),1) troughs((pulse+1),1)], [baseheight(pulse) troughs((pulse+1),2)], 'Color', 'c');
                plot(app.UIAxesAnal,left_edge(pulse),platheight(pulse), 'Color', 'm', 'Marker', 'O', 'LineStyle', 'none');
                plot(app.UIAxesAnal,left_edge(pulse),platheight(pulse), 'Color', 'm', 'Marker', '+', 'LineStyle', 'none');
            end
        end
        
        for pulse=maxSize:-1:1
            inc_high = minIndicies(pulse+1);
            for i=minIndicies(pulse+1):-1:maxIndicies(pulse);
                inc_low = inc_high;
                low_t = i + 1;
                inc_high = curve(i,2);
                if (inc_high >= platheight(pulse))
                    break
                end
            end
            
            %calculate slope
            m = (inc_high - inc_low)*10;
            %calculate difference between lower bound and 40percent
            p_diff = platheight(pulse) - inc_low;
            %calculate difference between lower bound time and 40percent time
            t_diff = p_diff/m;
            %estimate location of edge
            right_edge(pulse) = curve(low_t,1) - t_diff;
            if doplots
                plot(app.UIAxesAnal,right_edge(pulse),platheight(pulse), 'Color', 'm', 'Marker', 'O', 'LineStyle', 'none');
                plot(app.UIAxesAnal,right_edge(pulse),platheight(pulse), 'Color', 'm', 'Marker', '+', 'LineStyle', 'none');
                hp = line(app.UIAxesAnal,[left_edge(pulse) right_edge(pulse)], [platheight(pulse) platheight(pulse)], 'Color', 'm','LineWidth',2);
            end
        end

        if doplots
            legend(app.UIAxesAnal,[hp,hb],'Plateau Width','Base Width','Location','Northwest');
            % title(['Curve ', curveNum]);
            % xlabel('Time (min)');
            % ylabel('Fura-2  Ratio');
            hold(app.UIAxesAnal,'off');
        end

        %Create an array of plateau widths and plateau fractions
        platwidth = right_edge - left_edge;
        pf = platwidth ./ basewidth;
        
        % %Determine average plateau width
        % avgpf = mean(pf);
        % avgpa = mean(pArea);
        
        %Determine average plateau width added by ss 10-9-19
        avgpf = mean(pf);
        avgpa = mean(pArea);
        avgpw = mean(platwidth);
        avgbw = mean(basewidth);
        
        for m =1:PeakSize
            Silencep(1:m)=(basewidth(1:m)-platwidth(1:m));
        end
        
        avgsp = mean(Silencep);
        
        % expoData1 = [expoData1 avgpf(1) avgpa(1) avgpw(1) avgbw(1) avgsp(1)];

        end
                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: ChooseFiletoImportButton
        function ChooseFiletoImportButtonPushed(app, event)
            [filename, path]=uigetfile({'*.xlsx';'*.xls'}, 'Select Excel File');
            if isequal(filename, 0) % Check if user canceled
                return;
            else
                app.fullfilename = fullfile(path, filename);
                app.FileNameEditField.Value=filename;
            end
            %if length(sheetnames(fullfilename))>1 % file contains multiple sheets, ask user to select one
            app.SheetDropDown.Items = sheetnames(app.fullfilename);
            %app.SheetDropDown.Items{1}
            app.SheetDropDown.Value = app.SheetDropDown.Items{1};

            % set output file name
            [~,filenamesansext,~] = fileparts(app.fullfilename);
            app.datapath = path;
            app.SaveFileEditField.Value = [filenamesansext  '_analyzed-'  char(datetime('now','Format','y-MMM-d'))  '.xlsx'];

            app.initializeTables();

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
            numericcols = zeros(size(app.cols2preview));
            for ii=1:length(numericcols) % remove non-numeric cols from preview list
                numericcols(ii)=isnumeric(app.UITableImported.Data{:,app.cols2preview(ii)});
            end
            app.cols2preview=app.cols2preview(logical(numericcols)); 
            if isempty(app.cols2preview), app.cols2preview = app.firstdatacol; end % fall back to first col if no valid data cols selected
            app.updatePreviewPlot();
        end

        % Value changed function: SensorTypeDropDown
        function SensorTypeDropDownValueChanged(app, event)
            sensorType = app.SensorTypeDropDown.Value;
            
        end

        % Cell edit callback: UITableTimeWindows
        function UITableTimeWindowsCellEdit(app, event)
            indices = event.Indices;
            newData = event.NewData;
            % assignin("base","timewinsa",app.UITableTimeWindows)
            if app.UITableTimeWindows.Data{indices(1),2}{1}>=app.UITableTimeWindows.Data{indices(1),3}{1}
                app.UITableTimeWindows.Data{indices(1),3}={app.alltimes(end)}; % the values in the time window are cell arrays, so assign it that way
            end
            app.updatePreviewPlot();
            app.updateUITableAnalVals();
            app.updateAnalPlot();
        end

        % Menu selected function: AddRow
        function AddRowMenuSelected(app, event)
            nRegions = size(app.UITableImported.Data,2)+1-app.firstdatacol; % -1 so we don't count the time column
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
                nRegions = size(app.UITableImported.Data,2)+1-app.firstdatacol; % -1 so we don't count the time column

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
            app.updateWaveletPlot();
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
            [~,fullsavematfile,~]=fileparts(fullsavefile)
            % fullsavematfile = fullsavematfile+'.mat'
            
            % temp table to store region data in excel sheets
            t = table()
            regions2output = find(app.UITableAnalVals.Data.Ignore==false);
            
            if isfile(fullsavefile)
               % file already exists, ask if you want to overwrite
               msg = "File exists, overwrite? Otherwise cancel and update filename before saving";
               title = "Output file exists";
               selection = uiconfirm(app.UIFigure,msg,title, ...
                    "Options",["Overwrite","Cancel"], ...
                    "DefaultOption",2);
               if selection == "Overwrite"
                   writetable(app.OutTable,fullsavefile)
                   for row=1:size(app.OutTable,1)
                    t=table()
                    % t.peaksX = app.UITableAnalVals.Data.peaks{row}(:,1);
                    % t.peaksY = app.UITableAnalVals.Data.peaks{row}(:,2);
                    % t.troughsX = app.UITableAnalVals.Data.troughs{row}(:,1);
                    % t.troughsY = app.UITableAnalVals.Data.troughs{row}(:,2);
                    t.freguqncy = app.UITableAnalVals.Data.aveWavelet{row}(:,1);
                    t.aveWavelet = app.UITableAnalVals.Data.aveWavelet{row}(:,2);
                    writetable(t,fullsavefile,Sheet=app.OutTable.Region{row}(1:3));
                   end
                   % save((fullsavematfile+".mat"))
               end
            else
                writetable(app.OutTable,fullsavefile)
                % assignin('base','OutTable',app.OutTable);
                % assignin('base','AnalTable',app.UITableAnalVals.Data);
                for row=1:size(app.OutTable,1)
                    t=table()
                    % t.peaksX = app.UITableAnalVals.Data.peaks{row}(:,1);
                    % t.peaksY = app.UITableAnalVals.Data.peaks{row}(:,2);
                    % t.troughsX = app.UITableAnalVals.Data.troughs{row}(:,1);
                    % t.troughsY = app.UITableAnalVals.Data.troughs{row}(:,2);
                    t.freguqncy = app.UITableAnalVals.Data.aveWavelet{row}(:,1);
                    t.aveWavelet = app.UITableAnalVals.Data.aveWavelet{row}(:,2);
                    writetable(t,fullsavefile,Sheet=app.OutTable.Region{row}(1:3));
                end
                % save((fullsavematfile+".mat"))
            end
            % assignin('base','xdata',app.UITableAnalVals.Data.xdata);
            % assignin('base','ydata',app.UITableAnalVals.Data.ydata);
            % assignin('base','peaks',app.UITableAnalVals.Data.peaks);
            % assignin('base','troughs',app.UITableAnalVals.Data.troughs);
            
        end

        % Value changed function: SaveFileEditField
        function SaveFileEditFieldValueChanged2(app, event)
            %savefilename = app.SaveFileEditField.Value;
               
        end

        % Button pushed function: UpdateOutputButton
        function UpdateOutputButtonPushed(app, event)
            app.updateOutTable();
        end

        % Value changed function: IgnoreOutliersCheckBox
        function IgnoreOutliersCheckBoxValueChanged(app, event)
            value = app.IgnoreOutliersCheckBox.Value;
            % if the checkbox is changed, toggle the vals in the table
            % between NaN and the original vals stored in app.artifacts
            if value
                app.UITableImported.Data{app.artifacts(:,1),2:end}=NaN;
            else
                app.UITableImported.Data{app.artifacts(:,1),2:end}=app.artifacts(:,2:end);
            end
            app.updatePreviewPlot();
            app.updateUITableAnalVals();
            app.updateAnalPlot();
        end

        % Value changed function: UsefindpeaksCheckBox
        function UsefindpeaksCheckBoxValueChanged(app, event)
            value = app.UsefindpeaksCheckBox.Value;
            app.updateUITableAnalVals();
        end

        % Value changed function: SheetDropDown
        function SheetDropDownValueChanged(app, event)
            value = app.SheetDropDown.Value;
            app.initializeTables()
        end

        % Button pushed function: LoadstateButton
        function LoadstateButtonPushed(app, event)
            % [filename, path]=uigetfile({'*.mat'}, 'Select MAT File');
            % load(fullfile(path,filename));
        end

        % Value changed function: ShowPlatsCheckBox
        function ShowPlatsCheckBoxValueChanged(app, event)
            value = app.ShowPlatsCheckBox.Value;
            
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.AutoResizeChildren = 'off';
            app.UIFigure.Position = [100 100 1600 1280];
            app.UIFigure.Name = 'MATLAB App';

            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {'2x', '3x', '5x', '4x', '1x'};
            app.GridLayout.RowHeight = {'fit', '3x', '2x', 180, '6x', '6x', 'fit'};

            % Create UIAxesSelector
            app.UIAxesSelector = uiaxes(app.GridLayout);
            title(app.UIAxesSelector, 'Title')
            xlabel(app.UIAxesSelector, 'X')
            ylabel(app.UIAxesSelector, 'Y')
            zlabel(app.UIAxesSelector, 'Z')
            app.UIAxesSelector.Layout.Row = [1 3];
            app.UIAxesSelector.Layout.Column = [3 5];

            % Create UITableAnalVals
            app.UITableAnalVals = uitable(app.GridLayout);
            app.UITableAnalVals.ColumnName = {'Region'; 'TWindow'; 'Threshold'; 'Ignore'; 'tmin'; 'tmax'; 'tminind'; 'tmaxind'; 'xdata'; 'ydata'; 'peaks'; 'troughs'; 'aveWavelet'};
            app.UITableAnalVals.RowName = {};
            app.UITableAnalVals.ColumnEditable = [false false true true false false false false];
            app.UITableAnalVals.CellEditCallback = createCallbackFcn(app, @UITableAnalValsCellEdit, true);
            app.UITableAnalVals.SelectionChangedFcn = createCallbackFcn(app, @UITableAnalValsSelectionChanged, true);
            app.UITableAnalVals.Layout.Row = 5;
            app.UITableAnalVals.Layout.Column = [1 2];

            % Create ControlsPanel
            app.ControlsPanel = uipanel(app.GridLayout);
            app.ControlsPanel.AutoResizeChildren = 'off';
            app.ControlsPanel.Title = 'Controls';
            app.ControlsPanel.Layout.Row = 4;
            app.ControlsPanel.Layout.Column = [1 2];

            % Create SensorTypeLabel
            app.SensorTypeLabel = uilabel(app.ControlsPanel);
            app.SensorTypeLabel.HorizontalAlignment = 'right';
            app.SensorTypeLabel.Position = [19 104 75 22];
            app.SensorTypeLabel.Text = 'Sensor Type:';

            % Create SensorTypeDropDown
            app.SensorTypeDropDown = uidropdown(app.ControlsPanel);
            app.SensorTypeDropDown.Items = {'None', 'Calcium', 'Lactate', 'ATP/ADP'};
            app.SensorTypeDropDown.ValueChangedFcn = createCallbackFcn(app, @SensorTypeDropDownValueChanged, true);
            app.SensorTypeDropDown.Tooltip = {'Select the sensor'};
            app.SensorTypeDropDown.Position = [109 104 100 22];
            app.SensorTypeDropDown.Value = 'None';

            % Create Instructions
            app.Instructions = uilabel(app.ControlsPanel);
            app.Instructions.Position = [240 50 280 103];
            app.Instructions.Text = {'Instructions:'; '1. Choose file, select sheet from dropdown '; '2. Update time window start and end'; '3. Right click table to add or remove time windows'; '4. Select row(s) below to plot, update params'; '6. Click Update Ouput to refresh output table'; '7. Save to output file'};

            % Create DetrendCheckBox
            app.DetrendCheckBox = uicheckbox(app.ControlsPanel);
            app.DetrendCheckBox.ValueChangedFcn = createCallbackFcn(app, @DetrendCheckBoxValueChanged, true);
            app.DetrendCheckBox.Tooltip = {'Flatten the curves by removing the slope.'};
            app.DetrendCheckBox.Text = 'Detrend';
            app.DetrendCheckBox.Position = [19 56 65 22];
            app.DetrendCheckBox.Value = true;

            % Create IgnoreOutliersCheckBox
            app.IgnoreOutliersCheckBox = uicheckbox(app.ControlsPanel);
            app.IgnoreOutliersCheckBox.ValueChangedFcn = createCallbackFcn(app, @IgnoreOutliersCheckBoxValueChanged, true);
            app.IgnoreOutliersCheckBox.Tooltip = {'Remove artifacts (shown as x in preview window) by identifying peaks that occur across all regions, are less than 1 timepoint wide, and are more than 2 standard deviations in prominence'};
            app.IgnoreOutliersCheckBox.Text = 'Ignore Outliers';
            app.IgnoreOutliersCheckBox.Position = [19 79 100 22];

            % Create UsefindpeaksCheckBox
            app.UsefindpeaksCheckBox = uicheckbox(app.ControlsPanel);
            app.UsefindpeaksCheckBox.ValueChangedFcn = createCallbackFcn(app, @UsefindpeaksCheckBoxValueChanged, true);
            app.UsefindpeaksCheckBox.Tooltip = {'Peak and trough detection method:'; 'Default is to use the code developed by S. Sdao, but if checked, use the matlab findpeaks() method instead (JDR)'};
            app.UsefindpeaksCheckBox.Text = 'Use findpeaks()';
            app.UsefindpeaksCheckBox.Position = [19 33 105 22];

            % Create SheetLabel
            app.SheetLabel = uilabel(app.ControlsPanel);
            app.SheetLabel.HorizontalAlignment = 'right';
            app.SheetLabel.Position = [19 131 39 22];
            app.SheetLabel.Text = 'Sheet:';

            % Create SheetDropDown
            app.SheetDropDown = uidropdown(app.ControlsPanel);
            app.SheetDropDown.Items = {'Sheet1'};
            app.SheetDropDown.ValueChangedFcn = createCallbackFcn(app, @SheetDropDownValueChanged, true);
            app.SheetDropDown.Position = [73 131 100 22];
            app.SheetDropDown.Value = 'Sheet1';

            % Create ShowPlatsCheckBox
            app.ShowPlatsCheckBox = uicheckbox(app.ControlsPanel);
            app.ShowPlatsCheckBox.ValueChangedFcn = createCallbackFcn(app, @ShowPlatsCheckBoxValueChanged, true);
            app.ShowPlatsCheckBox.Text = 'Show Plats';
            app.ShowPlatsCheckBox.Position = [134 33 82 22];

            % Create UpdateOutputButton
            app.UpdateOutputButton = uibutton(app.ControlsPanel, 'push');
            app.UpdateOutputButton.ButtonPushedFcn = createCallbackFcn(app, @UpdateOutputButtonPushed, true);
            app.UpdateOutputButton.BackgroundColor = [0.0667 0.4431 0.7451];
            app.UpdateOutputButton.Position = [19 8 259 23];
            app.UpdateOutputButton.Text = 'Update Output with Current Analysis Params';

            % Create TabGroup
            app.TabGroup = uitabgroup(app.GridLayout);
            app.TabGroup.AutoResizeChildren = 'off';
            app.TabGroup.Layout.Row = [4 6];
            app.TabGroup.Layout.Column = [3 5];

            % Create PlotsTab
            app.PlotsTab = uitab(app.TabGroup);
            app.PlotsTab.AutoResizeChildren = 'off';
            app.PlotsTab.Title = 'Plots';

            % Create UIAxesAnal
            app.UIAxesAnal = uiaxes(app.PlotsTab);
            title(app.UIAxesAnal, 'Title')
            xlabel(app.UIAxesAnal, 'X')
            ylabel(app.UIAxesAnal, 'Y')
            zlabel(app.UIAxesAnal, 'Z')
            app.UIAxesAnal.Position = [0 31 944 833];

            % Create WaveletviewTab
            app.WaveletviewTab = uitab(app.TabGroup);
            app.WaveletviewTab.AutoResizeChildren = 'off';
            app.WaveletviewTab.Title = 'Wavelet view';

            % Create UIAxesWVLT
            app.UIAxesWVLT = uiaxes(app.WaveletviewTab);
            title(app.UIAxesWVLT, 'Wavelet Analysis')
            xlabel(app.UIAxesWVLT, 'X')
            ylabel(app.UIAxesWVLT, 'Y')
            zlabel(app.UIAxesWVLT, 'Z')
            app.UIAxesWVLT.Position = [1 31 731 634];

            % Create UIAxesWVLTx
            app.UIAxesWVLTx = uiaxes(app.WaveletviewTab);
            app.UIAxesWVLTx.Position = [1 664 731 200];

            % Create UIAxesWVLTy
            app.UIAxesWVLTy = uiaxes(app.WaveletviewTab);
            app.UIAxesWVLTy.Position = [731 31 212 634];

            % Create DataSummaryPanel
            app.DataSummaryPanel = uipanel(app.GridLayout);
            app.DataSummaryPanel.AutoResizeChildren = 'off';
            app.DataSummaryPanel.Title = 'Data Summary';
            app.DataSummaryPanel.Layout.Row = 6;
            app.DataSummaryPanel.Layout.Column = [1 2];

            % Create UITableOutput
            app.UITableOutput = uitable(app.DataSummaryPanel);
            app.UITableOutput.ColumnName = {'Region'; 'nPeaks'; 'avgBaseline'; 'avgPeak'; 'avgPeakAmplitude'; 'Period'; 'Threshold'; 'PlatFrac'; 'ActiveArea'; 'AvePlatWidth'; 'AveBaseWidth'; 'SlientPhase'; 'AveYval'; 'Notes'};
            app.UITableOutput.RowName = {};
            app.UITableOutput.ColumnEditable = [false false false false false false false false false false false false false true];
            app.UITableOutput.Position = [4 -1 516 248];

            % Create UITableTimeWindows
            app.UITableTimeWindows = uitable(app.GridLayout);
            app.UITableTimeWindows.ColumnName = {'Time Window Label'; 'Start time'; 'End time'};
            app.UITableTimeWindows.RowName = {};
            app.UITableTimeWindows.ColumnEditable = true;
            app.UITableTimeWindows.RowStriping = 'off';
            app.UITableTimeWindows.CellEditCallback = createCallbackFcn(app, @UITableTimeWindowsCellEdit, true);
            app.UITableTimeWindows.Layout.Row = 3;
            app.UITableTimeWindows.Layout.Column = [1 2];

            % Create UITableImported
            app.UITableImported = uitable(app.GridLayout);
            app.UITableImported.ColumnName = {'Column 1'; 'Column 2'; 'Column 3'; 'Column 4'};
            app.UITableImported.RowName = {};
            app.UITableImported.SelectionType = 'column';
            app.UITableImported.SelectionChangedFcn = createCallbackFcn(app, @UITableImportedSelectionChanged, true);
            app.UITableImported.Layout.Row = 2;
            app.UITableImported.Layout.Column = [1 2];

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

            % Create SaveOutputButton
            app.SaveOutputButton = uibutton(app.GridLayout, 'push');
            app.SaveOutputButton.ButtonPushedFcn = createCallbackFcn(app, @SaveOutputButtonPushed, true);
            app.SaveOutputButton.BackgroundColor = [0.0667 0.4431 0.7451];
            app.SaveOutputButton.Layout.Row = 7;
            app.SaveOutputButton.Layout.Column = 1;
            app.SaveOutputButton.Text = 'Save Output';

            % Create SaveFileEditField
            app.SaveFileEditField = uieditfield(app.GridLayout, 'text');
            app.SaveFileEditField.ValueChangedFcn = createCallbackFcn(app, @SaveFileEditFieldValueChanged2, true);
            app.SaveFileEditField.Layout.Row = 7;
            app.SaveFileEditField.Layout.Column = 2;

            % Create LoadstateButton
            app.LoadstateButton = uibutton(app.GridLayout, 'push');
            app.LoadstateButton.ButtonPushedFcn = createCallbackFcn(app, @LoadstateButtonPushed, true);
            app.LoadstateButton.BackgroundColor = [0.149 0.149 0.149];
            app.LoadstateButton.Layout.Row = 7;
            app.LoadstateButton.Layout.Column = 5;
            app.LoadstateButton.Text = 'Load state';

            % Create ContextMenu
            app.ContextMenu = uicontextmenu(app.UIFigure);

            % Create AddRow
            app.AddRow = uimenu(app.ContextMenu);
            app.AddRow.MenuSelectedFcn = createCallbackFcn(app, @AddRowMenuSelected, true);
            app.AddRow.Text = 'AddRow';

            % Create DelRow
            app.DelRow = uimenu(app.ContextMenu);
            app.DelRow.MenuSelectedFcn = createCallbackFcn(app, @DelRowMenuSelected, true);
            app.DelRow.Text = 'DelLastRow';
            
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