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
        signaldata = selectedNodes.handle.UserData;
        if ~isempty(signaldata) && strcmp(signaldata.datatype,'syscal')
            mainHandles.syscalstats = signaldata;
            set(mainHandles.signaltypetext,'String',[selectedNodes.getName.char ' selected']);
        end
        if ~isempty(signaldata) && isfield(signaldata,'audio')% If there's audio data saved in the leaf...

%             Details = signaldata;
%             if isfield(Details,'datatype'), Details = rmfield(Details,'datatype'); end
%             if isfield(Details,'funcallback'), Details = rmfield(Details,'funcallback'); end
%             if isfield(Details,'properties'), Details = rmfield(Details,'properties'); end %#ok : used in lines 30 and 31
%             audiodatatext = evalc('Details');
%             clear('Details')

            % method of getting audiodatatext (replaces code above)
            audiodatatext = [char(10),' audio [',num2str(size(signaldata.audio,1)),'x',...
                num2str(size(signaldata.audio,2)),'x',...
                num2str(size(signaldata.audio,3)),'x',...
                num2str(size(signaldata.audio,4)),'x',...
                num2str(size(signaldata.audio,5)),'x',...
                num2str(size(signaldata.audio,6)),']'];
            if isfield(signaldata, 'audio2')
                audiodatatext = [audiodatatext, char(10), ' audio2 [',num2str(size(signaldata.audio2,1)),'x',...
                num2str(size(signaldata.audio2,2)),']'];
            end
            if isfield(signaldata, 'fs')
                audiodatatext = [audiodatatext, char(10), ' fs ',num2str(signaldata.fs), ' Hz'];
            end
            if isfield(signaldata,'chanID')
                if length(signaldata.chanID) == 1
                    audiodatatext = [audiodatatext, char(10), ' chanID {', signaldata.chanID{1,1}, '}'];
                elseif length(signaldata.chanID) == 2
                    audiodatatext = [audiodatatext, char(10), ' chanID {', signaldata.chanID{1,1}, ';', signaldata.chanID{2,1}, '}'];
                else
                    audiodatatext = [audiodatatext, char(10), ' chanID {', signaldata.chanID{1,1}, ';', signaldata.chanID{2,1}, '...}'];
                end
            end
            if isfield(signaldata,'bandID')
                if length(signaldata.bandID) > 1
                audiodatatext = [audiodatatext,char(10), ' bandID [',...
                    num2str(min(signaldata.bandID)), ' to ',num2str(max(signaldata.bandID)) ']'];
                else
                    audiodatatext = [audiodatatext,char(10), ' bandID [',...
                    num2str(signaldata.bandID), ']'];
                end
            end
            if isfield(signaldata,'cal')
                if length(signaldata.cal) == 1
                    audiodatatext = [audiodatatext, char(10), ' cal ', num2str(signaldata.cal)];
                elseif length(signaldata.cal) == 2
                    audiodatatext = [audiodatatext, char(10), ' cal [', num2str(signaldata.cal(1)), ',', num2str(signaldata.cal(2)), ']'];
                elseif length(signaldata.cal) == 3
                    audiodatatext = [audiodatatext, char(10), ' cal [', num2str(signaldata.cal(1)), ';', num2str(signaldata.cal(2)), ';', num2str(signaldata.cal(3)), ']'];
                elseif length(signaldata.cal) == 4
                    audiodatatext = [audiodatatext, char(10), ' cal [', num2str(signaldata.cal(1)), ';', num2str(signaldata.cal(2)), ';', num2str(signaldata.cal(3)), ';', num2str(signaldata.cal(4)), ']'];
                else 
                    audiodatatext = [audiodatatext, char(10), ' cal [', num2str(signaldata.cal(1)), ';', num2str(signaldata.cal(2)), ';', num2str(signaldata.cal(3)), ';', num2str(signaldata.cal(4)), '...]'];
                end
            end
            
            set(mainHandles.audiodatatext,'String',['Selected: ' selectedNodes.getName.char audiodatatext]); % Output contents in textbox below the tree
            set(mainHandles.datatext,'Visible','off');
            set(mainHandles.datatext,'String',[]);
            set(mainHandles.data_panel1,'Visible','off');
            set(mainHandles.data_panel2,'Visible','off'); 
            set(mainHandles.tools_panel,'Visible','on');
            if ~strcmp(signaldata.datatype,'syscal')
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
            if size(signaldata.audio,1) <= mainHandles.Settings.maxtimetodisplay*signaldata.fs
                set([mainHandles.Tf_time,mainHandles.Tf_freq],'String',num2str(size(signaldata.audio,1)/signaldata.fs))
                mainHandles.Tf_time_IN = size(signaldata.audio,1)/signaldata.fs;
                mainHandles.Tf_freq_IN = size(signaldata.audio,1)/signaldata.fs;
            else
                set([mainHandles.Tf_time,mainHandles.Tf_freq],'String',num2str(mainHandles.Settings.maxtimetodisplay))
                mainHandles.Tf_time_IN = mainHandles.Settings.maxtimetodisplay;
                mainHandles.Tf_freq_IN = mainHandles.Settings.maxtimetodisplay;
            end
            set([mainHandles.text20,mainHandles.text21],'String',[num2str(size(signaldata.audio,1)/signaldata.fs) ' s'])
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
            if isfield(signaldata,'properties'), set(mainHandles.properties_btn,'Visible','on'); else set(mainHandles.properties_btn,'Visible','off'); end
            setappdata(hMain,'testsignal', signaldata); % Set leaf contents in the 'desktop'
            if ~ismatrix(signaldata.audio)
                set(mainHandles.channel_panel,'Visible','on');
                set(mainHandles.IN_nchannel,'String','1');
                if size(signaldata.audio,5) >1
                    if get(mainHandles.chandisp, 'Value') == 2
                        set(mainHandles.tchannels,'String',['/ ' num2str(size(signaldata.audio,5))]);
                    else
                        set(mainHandles.tchannels,'String',['/ ' num2str(size(signaldata.audio,2))]);
                    end
                    set(mainHandles.chandisp,'String',{'Input Channels'; 'Output Channels'})
                else
                    set(mainHandles.tchannels,'String',['/ ' num2str(size(signaldata.audio,2))]);
                    set(mainHandles.chandisp,'String',{'Input Channels'})
                end
                if ndims(signaldata.audio) == 3, cmap = colormap(hsv(size(signaldata.audio,3))); end
                if ndims(signaldata.audio) == 4, cmap = colormap(copper(size(signaldata.audio,4))); end
                if ndims(signaldata.audio) >= 5, cmap = colormap(cool(size(signaldata.audio,5))); end
                set(mainHandles.aarae,'DefaultAxesColorOrder',cmap)
            else
                cmap = colormap(lines(size(signaldata.audio,2)));
                set(mainHandles.aarae,'DefaultAxesColorOrder',cmap)
            end
            pause(0.001)
            refreshplots(mainHandles,'time')
            pause(0.001)
            refreshplots(mainHandles,'freq')
            %if isfield(signaldata,'audio2') && ~isempty(signaldata.audio2)% && ismatrix(signaldata.audio)
                set(mainHandles.IR_btn,'Enable','on');
            %else
            %    set(mainHandles.IR_btn,'Enable','off');% Display process IR button if selection is a measurement based on a sine sweep
            %end
            pause(0.001)
        elseif ~isempty(signaldata) && ~isfield(signaldata,'audio')% If there's data saved in the leaf but not audio...
            plot(mainHandles.axestime,0,0)
            plot(mainHandles.axesfreq,0,0)
            set(mainHandles.axestime,'Visible','off');
            set(mainHandles.axesfreq,'Visible','off');
            set(mainHandles.time_units_display,'Visible','off');
            set(mainHandles.freq_units_display,'Visible','off');
            set([mainHandles.text16,mainHandles.text17,mainHandles.text18,mainHandles.text19,mainHandles.text20,mainHandles.text21],'Visible','off')
            if isfield(signaldata,'data') || isfield(signaldata,'tables')
                set(mainHandles.axesdata,'Visible','on');
            else
                set(mainHandles.axesdata,'Visible','off');
            end
            set(mainHandles.audiodatatext,'String',[]);
            datatext = evalc('signaldata');
            set(mainHandles.datatext,'Visible','on');
            set(mainHandles.datatext,'String',['Selected: ' selectedNodes.getName.char datatext]); % Output contents in textbox below the tree
            cla(mainHandles.axesdata)
        % Enable chart data visualization
            if isfield(signaldata,'data')
                set(mainHandles.data_panel1,'Visible','on');
                filltable(signaldata,mainHandles)
                doresultplot(mainHandles,mainHandles.axesdata)
                mainHandles.tabledata = get(mainHandles.cattable,'Data');
            else
                set(mainHandles.data_panel1,'Visible','off');
            end
        % Enable table data visualization
            if isfield(signaldata,'tables')
                set(mainHandles.data_panel2,'Visible','on');
                set(mainHandles.ntable_popup,'String',{signaldata.tables(:).Name})%cellstr(num2str((1:length(signaldata.tables))')));
                ntable = get(mainHandles.ntable_popup,'Value');
                Xvalues = get(mainHandles.Xvalues_sel,'SelectedObject');
                Xvalues = get(Xvalues,'tag');
                switch Xvalues
                    case 'radiobutton1'
                        if size(signaldata.tables(ntable).Data,2) < get(mainHandles.Yvalues_box,'Value'), set(mainHandles.Yvalues_box,'Value',1); end
                        bar(mainHandles.axesdata,signaldata.tables(ntable).Data(:,get(mainHandles.Yvalues_box,'Value')),'FaceColor',[0 0.5 0.5])
                        set(mainHandles.axesdata,'Xtick',1:length(signaldata.tables(ntable).RowName),'XTickLabel',signaldata.tables(ntable).RowName)
                        set(mainHandles.Xvalues_box,'String',signaldata.tables(ntable).RowName,'Value',1)
                        set(mainHandles.Yvalues_box,'String',signaldata.tables(ntable).ColumnName)
                    case 'radiobutton2'
                        if size(signaldata.tables(ntable).Data,1) < get(mainHandles.Yvalues_box,'Value'), set(mainHandles.Yvalues_box,'Value',1); end
                        bar(mainHandles.axesdata,signaldata.tables(ntable).Data(get(mainHandles.Yvalues_box,'Value'),:),'FaceColor',[0 0.5 0.5])
                        set(mainHandles.axesdata,'Xtick',1:length(signaldata.tables(ntable).ColumnName),'XTickLabel',signaldata.tables(ntable).ColumnName)
                        set(mainHandles.Xvalues_box,'String',signaldata.tables(ntable).ColumnName,'Value',1)
                        set(mainHandles.Yvalues_box,'String',signaldata.tables(ntable).RowName)
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
            if isfield(signaldata,'properties'), set(mainHandles.properties_btn,'Visible','on'); else set(mainHandles.properties_btn,'Visible','off'); end
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

function filltable(signaldata,handles)
fields = fieldnames(signaldata);
fields = fields(3:end-1);
categories = fields(mod(1:length(fields),2) == 1);
catdata = cell(size(categories));
catunits = cell(size(categories));
catorcont = cell(size(categories));
for n = 1:length(categories)
    catunits{n,1} = signaldata.(matlab.lang.makeValidName([categories{n,1} 'info'])).units;
    catorcont{n,1} = signaldata.(matlab.lang.makeValidName([categories{n,1} 'info'])).axistype;
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