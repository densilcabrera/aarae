% The tree-node selection callback
% This function is almost the brains behind the tree's response to
% selection commands, it controls also visibility of some functions
% depending on what type of signal is selected and its contents.
function mySelectFcn(tree, ~)
    pause on
    % Get handles of main window
    aarae_fig = findobj('Tag','aarae');
    mainHandles = guidata(aarae_fig);
    selectedNodes = tree.getSelectedNodes; % Get selected leaf
    if ~isempty(selectedNodes)
        % Call the 'desktop'
        hMain = getappdata(0,'hMain');
        contents = cell(length(selectedNodes),1);
        for n = 1:length(selectedNodes)
            contents{n,1} = selectedNodes(n).handle.Userdata;
        end
        selaudio = cellfun(@isfield,contents,repmat({'audio'},length(selectedNodes),1));
        seldata = cellfun(@isfield,contents,repmat({'data'},length(selectedNodes),1));
        if all(selaudio) || (length(seldata) <= 2 && all(seldata))
            set(mainHandles.compare_btn,'Enable','on')
            if all(selaudio), mainHandles.compareaudio = 1; end
            if all(seldata), mainHandles.compareaudio = 0; end
        else
            set(mainHandles.compare_btn,'Enable','off')
        end 
        selectedNodes = selectedNodes(1);
        audiodata = selectedNodes.handle.UserData;
        if ~isempty(audiodata) && strcmp(audiodata.datatype,'syscal')
            mainHandles.syscalstats = audiodata;
            set(mainHandles.signaltypetext,'String',[selectedNodes.getName.char ' selected']);
        end
        if ~isempty(audiodata) && isfield(audiodata,'audio')% If there's audio data saved in the leaf...

%             Details = audiodata;
%             if isfield(Details,'datatype'), Details = rmfield(Details,'datatype'); end
%             if isfield(Details,'funcallback'), Details = rmfield(Details,'funcallback'); end
%             if isfield(Details,'properties'), Details = rmfield(Details,'properties'); end %#ok : used in lines 30 and 31
%             audiodatatext = evalc('Details');
%             clear('Details')

            % method of getting audiodatatext (replaces code above)
            audiodatatext = [char(10),' audio [',num2str(size(audiodata.audio,1)),'x',...
                num2str(size(audiodata.audio,2)),'x',...
                num2str(size(audiodata.audio,3)),'x',...
                num2str(size(audiodata.audio,4)),'x',...
                num2str(size(audiodata.audio,5)),'x',...
                num2str(size(audiodata.audio,6)),']'];
            if isfield(audiodata, 'audio2')
                audiodatatext = [audiodatatext, char(10), ' audio2 [',num2str(size(audiodata.audio2,1)),'x',...
                num2str(size(audiodata.audio2,2)),']'];
            end
            if isfield(audiodata,'chanID')
                if length(audiodata.chanID) == 1
                    audiodatatext = [audiodatatext, char(10), ' chanID {', audiodata.chanID{1,1}, '}'];
                elseif length(audiodata.chanID) == 2
                    audiodatatext = [audiodatatext, char(10), ' chanID {', audiodata.chanID{1,1}, ';', audiodata.chanID{2,1}, '}'];
                else
                    audiodatatext = [audiodatatext, char(10), ' chanID {', audiodata.chanID{1,1}, ';', audiodata.chanID{2,1}, '...}'];
                end
            end
            if isfield(audiodata,'bandID')
                if length(audiodata.bandID) > 1
                audiodatatext = [audiodatatext,char(10), ' bandID [',...
                    num2str(min(audiodata.bandID)), ' to ',num2str(max(audiodata.bandID)) ']'];
                else
                    audiodatatext = [audiodatatext,char(10), ' bandID [',...
                    num2str(audiodata.bandID), ']'];
                end
            end
            if isfield(audiodata,'cal')
                if length(audiodata.cal) == 1
                    audiodatatext = [audiodatatext, char(10), ' cal ', num2str(audiodata.cal)];
                elseif length(audiodata.cal) == 2
                    audiodatatext = [audiodatatext, char(10), ' cal [', num2str(audiodata.cal(1)), ',', num2str(audiodata.cal(2)), ']'];
                elseif length(audiodata.cal) == 3
                    audiodatatext = [audiodatatext, char(10), ' cal [', num2str(audiodata.cal(1)), ';', num2str(audiodata.cal(2)), ';', num2str(audiodata.cal(3)), ']'];
                elseif length(audiodata.cal) == 4
                    audiodatatext = [audiodatatext, char(10), ' cal [', num2str(audiodata.cal(1)), ';', num2str(audiodata.cal(2)), ';', num2str(audiodata.cal(3)), ';', num2str(audiodata.cal(4)), ']'];
                else 
                    audiodatatext = [audiodatatext, char(10), ' cal [', num2str(audiodata.cal(1)), ';', num2str(audiodata.cal(2)), ';', num2str(audiodata.cal(3)), ';', num2str(audiodata.cal(4)), '...]'];
                end
            end
            
            set(mainHandles.audiodatatext,'String',['Selected: ' selectedNodes.getName.char audiodatatext]); % Output contents in textbox below the tree
            set(mainHandles.datatext,'Visible','off');
            set(mainHandles.datatext,'String',[]);
            set(mainHandles.data_panel1,'Visible','off');
            set(mainHandles.data_panel2,'Visible','off'); 
            set(mainHandles.tools_panel,'Visible','on');
            if ~strcmp(audiodata.datatype,'syscal')
                set([mainHandles.edit_btn mainHandles.cal_btn],'Enable','on');
            else
                set(mainHandles.edit_btn,'Enable','on');
                set(mainHandles.cal_btn,'Enable','off');
            end
            set(mainHandles.axestime,'Visible','on');
            set(mainHandles.axesfreq,'Visible','on');
            cla(mainHandles.axesdata)
            set(mainHandles.axesdata,'Visible','off');
            set(mainHandles.time_popup,'Visible','on');
            set(mainHandles.freq_popup,'Visible','on');
            set([mainHandles.text16,mainHandles.text17,mainHandles.text18,mainHandles.text19,mainHandles.text20,mainHandles.text21],'Visible','on')
            set([mainHandles.To_time,mainHandles.To_freq],'String','0')
            mainHandles.To_time_IN = 0;
            mainHandles.To_freq_IN = 0;
            if length(audiodata.audio) <= mainHandles.Settings.maxtimetodisplay*audiodata.fs
                set([mainHandles.Tf_time,mainHandles.Tf_freq],'String',num2str(length(audiodata.audio)/audiodata.fs))
                mainHandles.Tf_time_IN = length(audiodata.audio)/audiodata.fs;
                mainHandles.Tf_freq_IN = length(audiodata.audio)/audiodata.fs;
            else
                set([mainHandles.Tf_time,mainHandles.Tf_freq],'String',num2str(mainHandles.Settings.maxtimetodisplay))
                mainHandles.Tf_time_IN = mainHandles.Settings.maxtimetodisplay;
                mainHandles.Tf_freq_IN = mainHandles.Settings.maxtimetodisplay;
            end
            set([mainHandles.text20,mainHandles.text21],'String',[num2str(length(audiodata.audio)/audiodata.fs) ' s'])
            set(mainHandles.process_panel,'Visible','on');
            set(mainHandles.playback_panel,'Visible','on');
            set(mainHandles.channel_panel,'Visible','off');
            set(mainHandles.procat_box,'Value',1);
            set(mainHandles.proc_box,'Visible','off');
            set(mainHandles.proc_btn,'Visible','off');
            set(mainHandles.proc_help_btn,'Visible','off');
            set(mainHandles.analysis_panel,'Visible','on');
            set(mainHandles.funcat_box,'Value',1);
            set(mainHandles.fun_box,'Visible','off');
            set(mainHandles.analyze_btn,'Visible','off');
            set(mainHandles.analyser_help_btn,'Visible','off');
            if isfield(audiodata,'properties'), set(mainHandles.properties_btn,'Visible','on'); else set(mainHandles.properties_btn,'Visible','off'); end
            setappdata(hMain,'testsignal', audiodata); % Set leaf contents in the 'desktop'
            if ~ismatrix(audiodata.audio)
                set(mainHandles.channel_panel,'Visible','on');
                set(mainHandles.IN_nchannel,'String','1');
                set(mainHandles.tchannels,'String',['/ ' num2str(size(audiodata.audio,2))]);
                if ndims(audiodata.audio) == 3, cmap = colormap(hsv(size(audiodata.audio,3))); end
                if ndims(audiodata.audio) == 4, cmap = colormap(copper(size(audiodata.audio,4))); end
                if ndims(audiodata.audio) >= 5, cmap = colormap(cool(size(audiodata.audio,5))); end
                set(mainHandles.aarae,'DefaultAxesColorOrder',cmap)
            else
                cmap = colormap(lines(size(audiodata.audio,2)));
                set(mainHandles.aarae,'DefaultAxesColorOrder',cmap)
            end
            pause(0.001)
            refreshplots(mainHandles,'time')
            pause(0.001)
            refreshplots(mainHandles,'freq')
            %if isfield(audiodata,'audio2') && ~isempty(audiodata.audio2)% && ismatrix(audiodata.audio)
                set(mainHandles.IR_btn,'Enable','on');
            %else
            %    set(mainHandles.IR_btn,'Enable','off');% Display process IR button if selection is a measurement based on a sine sweep
            %end
            pause(0.001)
        elseif ~isempty(audiodata) && ~isfield(audiodata,'audio')% If there's data saved in the leaf but not audio...
            plot(mainHandles.axestime,0,0)
            plot(mainHandles.axesfreq,0,0)
            set(mainHandles.axestime,'Visible','off');
            set(mainHandles.axesfreq,'Visible','off');
            set(mainHandles.time_units_display,'Visible','off');
            set(mainHandles.freq_units_display,'Visible','off');
            set([mainHandles.text16,mainHandles.text17,mainHandles.text18,mainHandles.text19,mainHandles.text20,mainHandles.text21],'Visible','off')
            if isfield(audiodata,'data') || isfield(audiodata,'tables')
                set(mainHandles.axesdata,'Visible','on');
            else
                set(mainHandles.axesdata,'Visible','off');
            end
            set(mainHandles.audiodatatext,'String',[]);
            datatext = evalc('audiodata');
            set(mainHandles.datatext,'Visible','on');
            set(mainHandles.datatext,'String',['Selected: ' selectedNodes.getName.char datatext]); % Output contents in textbox below the tree
            cla(mainHandles.axesdata)
        % Enable chart data visualization
            if isfield(audiodata,'data')
                set(mainHandles.data_panel1,'Visible','on');
                filltable(audiodata,mainHandles)
                doresultplot(mainHandles,mainHandles.axesdata)
                mainHandles.tabledata = get(mainHandles.cattable,'Data');
            else
                set(mainHandles.data_panel1,'Visible','off');
            end
        % Enable table data visualization
            if isfield(audiodata,'tables')
                set(mainHandles.data_panel2,'Visible','on');
                set(mainHandles.ntable_popup,'String',{audiodata.tables(:).Name})%cellstr(num2str((1:length(audiodata.tables))')));
                ntable = get(mainHandles.ntable_popup,'Value');
                Xvalues = get(mainHandles.Xvalues_sel,'SelectedObject');
                Xvalues = get(Xvalues,'tag');
                switch Xvalues
                    case 'radiobutton1'
                        if size(audiodata.tables(ntable).Data,2) < get(mainHandles.Yvalues_box,'Value'), set(mainHandles.Yvalues_box,'Value',1); end
                        bar(mainHandles.axesdata,audiodata.tables(ntable).Data(:,get(mainHandles.Yvalues_box,'Value')),'FaceColor',[0 0.5 0.5])
                        set(mainHandles.axesdata,'Xtick',1:length(audiodata.tables(ntable).RowName),'XTickLabel',audiodata.tables(ntable).RowName)
                        set(mainHandles.Xvalues_box,'String',audiodata.tables(ntable).RowName,'Value',1)
                        set(mainHandles.Yvalues_box,'String',audiodata.tables(ntable).ColumnName)
                    case 'radiobutton2'
                        if size(audiodata.tables(ntable).Data,1) < get(mainHandles.Yvalues_box,'Value'), set(mainHandles.Yvalues_box,'Value',1); end
                        bar(mainHandles.axesdata,audiodata.tables(ntable).Data(get(mainHandles.Yvalues_box,'Value'),:),'FaceColor',[0 0.5 0.5])
                        set(mainHandles.axesdata,'Xtick',1:length(audiodata.tables(ntable).ColumnName),'XTickLabel',audiodata.tables(ntable).ColumnName)
                        set(mainHandles.Xvalues_box,'String',audiodata.tables(ntable).ColumnName,'Value',1)
                        set(mainHandles.Yvalues_box,'String',audiodata.tables(ntable).RowName)
                end
                ycontents = cellstr(get(mainHandles.Yvalues_box,'String'));
                ylabel(mainHandles.axesdata,ycontents{get(mainHandles.Yvalues_box,'Value')})
            else
                set(mainHandles.data_panel2,'Visible','off');
            end
            pause(0.001)
            set(mainHandles.tools_panel,'Visible','on');
            set([mainHandles.edit_btn mainHandles.cal_btn],'Enable','off')
            set([mainHandles.time_popup,mainHandles.freq_popup,mainHandles.smoothtime_popup,mainHandles.smoothfreq_popup,mainHandles.To_freq,mainHandles.Tf_freq,mainHandles.To_time,mainHandles.Tf_time],'Visible','off');
            set(mainHandles.logtime_chk,'Visible','off');
            set(mainHandles.logfreq_chk,'Visible','off');
            set(mainHandles.process_panel,'Visible','off');
            set(mainHandles.analysis_panel,'Visible','off');
            set(mainHandles.playback_panel,'Visible','off');
            set(mainHandles.channel_panel,'Visible','off');
            set([mainHandles.complextime mainHandles.complexfreq],'Visible','off')
            if isfield(audiodata,'properties'), set(mainHandles.properties_btn,'Visible','on'); else set(mainHandles.properties_btn,'Visible','off'); end
            setappdata(hMain,'testsignal', []);
            pause(0.001)
        else
            % If selection has no data, hide everything and don't display
            % data
            plot(mainHandles.axestime,0,0)
            semilogx(mainHandles.axesfreq,0,0)
            plot(mainHandles.axesdata,0,0)
            set(mainHandles.axestime,'Visible','off');
            set(mainHandles.axesfreq,'Visible','off');
            set(mainHandles.axesdata,'Visible','off');
            set([mainHandles.text16,mainHandles.text17,mainHandles.text18,mainHandles.text19,mainHandles.text20,mainHandles.text21],'Visible','off')
            set(mainHandles.audiodatatext,'String',[]);
            set(mainHandles.datatext,'Visible','off');
            set(mainHandles.datatext,'String',[]);
            set(mainHandles.data_panel1,'Visible','off');
            set(mainHandles.data_panel2,'Visible','off');
            set(mainHandles.tools_panel,'Visible','off');
            set(mainHandles.process_panel,'Visible','off');
            set(mainHandles.analysis_panel,'Visible','off');
            set(mainHandles.playback_panel,'Visible','off');
            set(mainHandles.channel_panel,'Visible','off');
            set([mainHandles.time_popup,mainHandles.freq_popup,mainHandles.smoothtime_popup,mainHandles.smoothfreq_popup,mainHandles.To_freq,mainHandles.Tf_freq,mainHandles.To_time,mainHandles.Tf_time],'Visible','off');
            set(mainHandles.logtime_chk,'Visible','off');
            set(mainHandles.logfreq_chk,'Visible','off');
            set([mainHandles.complextime mainHandles.complexfreq],'Visible','off')
            set(mainHandles.properties_btn,'Visible','off');
            setappdata(hMain,'testsignal', []);
            pause(0.001)
        end
    end
    pause off
    guidata(aarae_fig,mainHandles);
end  % mySelectFcn

function filltable(audiodata,handles)
fields = fieldnames(audiodata);
fields = fields(3:end-1);
categories = fields(mod(1:length(fields),2) == 1);
catdata = cell(size(categories));
catunits = cell(size(categories));
catorcont = cell(size(categories));
for n = 1:length(categories)
    catunits{n,1} = audiodata.(matlab.lang.makeValidName([categories{n,1} 'info'])).units;
    catorcont{n,1} = audiodata.(matlab.lang.makeValidName([categories{n,1} 'info'])).axistype;
    if islogical(catorcont{n,1}) && catorcont{n,1} == true
        catdata{n,1} = ':';
    else
        catdata{n,1} = '[1]';
    end
end
dat = [categories,catdata,catunits,catorcont];
set(handles.cattable, 'Data', dat);
naxis = length(find([catorcont{:}] == true));
switch naxis
    case 0
        set(handles.chartfunc_popup,'String',{'distributionPlot','boxplot','notboxplot'},'Value',1)
    case 1
        set(handles.chartfunc_popup,'String',{'plot','semilogx','semilogy','loglog','distributionPlot','boxplot','notboxplot'},'Value',1)
    case 2
        set(handles.chartfunc_popup,'String',{'mesh','surf','imagesc'},'Value',1)
    otherwise
        set(handles.chartfunc_popup,'String',{[]},'Value',1)
end
end % End of function filltable