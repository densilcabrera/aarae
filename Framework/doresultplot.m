function doresultplot(handles,haxes)

selectedNodes = handles.mytree.getSelectedNodes;
audiodata = selectedNodes(1).handle.UserData;
chartmenu = cellstr(get(handles.chartfunc_popup,'String'));
chartfunc = chartmenu{get(handles.chartfunc_popup,'Value')};
cattable = get(handles.cattable);
sel = strjoin(cattable.Data(:,2).',',');
if isempty(sel), sel = '[1]'; end
if ~isempty(findobj('Tag','tempaxes')) && get(haxes,'Parent') == handles.aarae
    delete(findobj('Tag','tempaxes'));
end
try
    tabledata = get(handles.cattable,'Data');
    catorcont = tabledata(:,4);
    if any(cellfun(@isempty,catorcont)), catorcont(cellfun(@isempty,catorcont)) = {false}; end
    naxis = length(find([catorcont{:}] == true));
    eval(['data = squeeze(audiodata.data(' sel '));'])
    if naxis < 2 % Visualization of 2D data
        cmap = colormap(hsv(size(data,2))); %#ok : data defined in line 15
        set(get(haxes,'Parent'),'DefaultAxesColorOrder',cmap)
        % Line plotting (plot, semilogx, semilogy, loglog)
        if ~strcmp(chartfunc,'distributionPlot') && ~strcmp(chartfunc,'boxplot') && ~strcmp(chartfunc,'notboxplot')
            cla(haxes,'reset')
            axis = find([catorcont{:}] == true);
            Xdata = audiodata.(genvarname(tabledata{axis(1,1),1}));
            if ~isnumeric(Xdata)
                if iscell(Xdata), Xdata = cell2mat(Xdata); end
            end
            if ~isequal(size(Xdata),size(audiodata.data))
                eval([chartfunc '(haxes,Xdata,data)'])
            else
                eval([chartfunc '(haxes,squeeze(Xdata(' sel ')),data)'])
            end
            xlabel(haxes,strrep([tabledata{axis(1,1),1} ' [' audiodata.(genvarname([tabledata{axis(1,1),1} 'info'])).units ']'],'_',' '),'HandleVisibility','on')
            ylabel(haxes,strrep(audiodata.datainfo.units,'_',' '),'HandleVisibility','on')
        else % Distribution, box plotting and not box plotting
            cla(haxes,'reset')
            catdim = find(cellfun(@isempty,tabledata(:,4)),1);
            if isempty(catdim), catdim = 1; end
            eval(['catdim = length(size(audiodata.data(' sel ')));'])
            if strcmp(chartfunc,'distributionPlot')
                eval([chartfunc '(haxes,data,''xNames'',audiodata.(genvarname(cattable.Data{catdim,1}))(' cattable.Data{catdim,2} '))'])
            end
            if strcmp(chartfunc,'boxplot')
                eval([chartfunc '(haxes,data,''labels'',audiodata.(genvarname(cattable.Data{catdim,1}))(' cattable.Data{catdim,2} '))'])
            end
            if strcmp(chartfunc,'notboxplot')
                notBoxPlot(data,[],[],[],haxes)
                haxes2 = findobj('Tag','tempaxes','Parent',get(haxes,'Parent'));
                if isvector(data), set(haxes2,'XLim',[0.5 1.5],'Xtick',1); end
                set(haxes,'XLim',get(haxes2,'XLim'))
                set(haxes2,'XTickLabel',eval(['audiodata.(genvarname(cattable.Data{catdim,1}))(' cattable.Data{catdim,2} ')']))
            end
            xlabel(haxes,strrep([cattable.Data{catdim,1} ' [' audiodata.(genvarname([cattable.Data{catdim,1} 'info'])).units ']'],'_',' '))
            ylabel(haxes,strrep(audiodata.datainfo.units,'_',' '))
        end
        handles.tabledata = tabledata;
        guidata(handles.aarae,handles)
    elseif naxis == 2 % Visualization of 3D data
        cla(haxes,'reset')
        [~,~,ext] = fileparts(handles.Settings.colormap);
        if isempty(ext)
            eval(['cmap = colormap(' lower(handles.Settings.colormap) '(128));']);
        else
            custom_colormap = importdata([cd '/Utilities/Custom_colormaps/' handles.Settings.colormap]);
            cmap = colormap(custom_colormap);
        end
        set(get(haxes,'Parent'),'DefaultAxesColorOrder',cmap)
        axis = find([catorcont{:}] == true);
        Xdata = audiodata.(genvarname(tabledata{axis(1,1),1})); 
        Ydata = audiodata.(genvarname(tabledata{axis(1,2),1}));
        if ~isnumeric(Xdata)
            if iscell(Xdata), Xdata = cell2mat(Xdata); end
        end
        if ~isnumeric(Ydata)
            if iscell(Ydata), Ydata = cell2mat(Ydata); end
        end
        %aaraecmap = importdata([cd '/Utilities/aaraecmap.mat']);
        if ~isequal([length(Ydata),length(Xdata)],size(data)), data = data'; end %#ok : Used in line above
        % Mesh and surf plotting
        if ~strcmp(chartfunc,'imagesc')
            eval([chartfunc '(haxes,Xdata,Ydata,data)'])
            %colormap(handles.Settings.colormap)
            xlabel(haxes,strrep([tabledata{axis(1,1),1} ' [' audiodata.(genvarname([tabledata{axis(1,1),1} 'info'])).units ']'],'_',' '),'HandleVisibility','on')
            ylabel(haxes,strrep([tabledata{axis(1,2),1} ' [' audiodata.(genvarname([tabledata{axis(1,2),1} 'info'])).units ']'],'_',' '),'HandleVisibility','on')
            zlabel(haxes,strrep(audiodata.datainfo.units,'_',' '),'HandleVisibility','on')
        else
        % Imagesc plotting
            eval([chartfunc '(Xdata,1:length(Ydata),data,''Parent'',haxes)'])
            set(haxes,'YTickLabel',num2str(Ydata(get(haxes,'Ytick')).'))
            set(haxes,'YDir','normal')
            %colormap(handles.Settings.colormap)
            xlabel(haxes,strrep([tabledata{axis(1,1),1} ' [' audiodata.(genvarname([tabledata{axis(1,1),1} 'info'])).units ']'],'_',' '))
            ylabel(haxes,strrep([tabledata{axis(1,2),1} ' [' audiodata.(genvarname([tabledata{axis(1,2),1} 'info'])).units ']'],'_',' '))
        end
        handles.tabledata = tabledata;
        guidata(handles.aarae,handles)
    else
        set(handles.cattable,'Data',handles.tabledata)
        warndlg('Cannot display plots with more than 2 main axis defined!','AARAE info','modal')
    end
catch error
    if length(find(cellfun(@length,cellfun(@str2num,tabledata(cell2mat(catorcont) == 0,2),'UniformOutput',false)) ~= 1)) == 2 && (strcmp(chartfunc,'plot') || strcmp(chartfunc,'semilogx') || strcmp(chartfunc,'semilogy') || strcmp(chartfunc,'loglog'))
    % 2D data subplotting
        selections = cellfun(@str2num,tabledata(:,2),'UniformOutput',false);
        if any(cellfun(@isempty,selections)), selections(cellfun(@isempty,selections)) = {1}; end
        cats = tabledata(cellfun(@length,selections) ~= 1,1);
        choice = questdlg('Would you like to display data subplots?',...
                          'AARAE subplots',...
                          cats{:},'Both',...
                          cats{1,1});
        switch choice
            case cats{1,1}
                dims = cellfun(@str2num,tabledata(cellfun(@length,selections) ~= 1,2),'UniformOutput',false);
                indexes = find(cellfun(@length,selections) ~= 1);
                dim2 = dims{1,1};
                figure;
                k = 1;
                [r,c] = subplotpositions(size(data,2), 0.5);
                for i = 1:size(data,2)
                    hsub = subplot(r,c,k);
                    newsel = strsplit(sel,',');
                    newsel{1,indexes(1,1)} = num2str(i);
                    newsel = strjoin(newsel,',');
                    if ~isequal(size(Xdata),size(audiodata.data))
                        eval([chartfunc '(hsub,Xdata,squeeze(audiodata.data(' newsel ')))'])
                    else
                        eval([chartfunc '(hsub,squeeze(Xdata(' newsel ')),squeeze(audiodata.data(' newsel ')))'])
                    end
                    %eval([chartfunc '(hsub,Xdata,squeeze(audiodata.data(' newsel ')))'])
                    xlabel(hsub,strrep([tabledata{axis(1,1),1} ' [' audiodata.(genvarname([tabledata{axis(1,1),1} 'info'])).units ']'],'_',' '))
                    ylabel(hsub,strrep(audiodata.datainfo.units,'_',' '))
                    if isnumeric(audiodata.(genvarname(cats{1,1}))), audiodata.(genvarname(cats{1,1})) = num2cell(audiodata.(genvarname(cats{1,1}))); end
                    dim2str = audiodata.(genvarname(cats{1,1})){dim2(i)};
                    if isnumeric(dim2str), dim2str = num2str(dim2str); end
                    title(hsub,dim2str)
                    k = k + 1;
                end
            case cats{2,1}
                dims = cellfun(@str2num,tabledata(cellfun(@length,selections) ~= 1,2),'UniformOutput',false);
                indexes = find(cellfun(@length,selections) ~= 1);
                dim3 = dims{2,1};
                figure;
                k = 1;
                [r,c] = subplotpositions(size(data,3), 0.5);
                for i = 1:size(data,3)
                    hsub = subplot(r,c,k);
                    newsel = strsplit(sel,',');
                    newsel{1,indexes(2,1)} = num2str(i);
                    newsel = strjoin(newsel,',');
                    if ~isequal(size(Xdata),size(audiodata.data))
                        eval([chartfunc '(hsub,Xdata,squeeze(audiodata.data(' newsel ')))'])
                    else
                        eval([chartfunc '(hsub,squeeze(Xdata(' newsel ')),squeeze(audiodata.data(' newsel ')))'])
                    end
                    %eval([chartfunc '(hsub,Xdata,squeeze(audiodata.data(' newsel ')))'])
                    xlabel(hsub,strrep([tabledata{axis(1,1),1} ' [' audiodata.(genvarname([tabledata{axis(1,1),1} 'info'])).units ']'],'_',' '))
                    ylabel(hsub,strrep(audiodata.datainfo.units,'_',' '))
                    if isnumeric(audiodata.(genvarname(cats{2,1}))), audiodata.(genvarname(cats{2,1})) = num2cell(audiodata.(genvarname(cats{2,1}))); end
                    dim3str = audiodata.(genvarname(cats{2,1})){dim3(i)};
                    if isnumeric(dim3str), dim3str = num2str(dim3str); end
                    title(hsub,dim3str)
                    k = k + 1;
                end
            case 'Both'
                dims = cellfun(@str2num,tabledata(cellfun(@length,selections) ~= 1,2),'UniformOutput',false);
                indexes = find(cellfun(@length,selections) ~= 1);
                dim2 = dims{1,1};
                dim3 = dims{2,1};
                figure;
                k = 1;
                [r,c] = subplotpositions(size(data,2)*size(data,3), 0.5);
                for i = 1:size(data,2)
                    for j = 1:size(data,3)
                        hsub = subplot(r,c,k);
                        newsel = strsplit(sel,',');
                        newsel{1,indexes(1,1)} = num2str(i);
                        newsel{1,indexes(2,1)} = num2str(j);
                        newsel = strjoin(newsel,',');
                        if ~isequal(size(Xdata),size(audiodata.data))
                            eval([chartfunc '(hsub,Xdata,squeeze(audiodata.data(' newsel ')))'])
                        else
                            eval([chartfunc '(hsub,squeeze(Xdata(' newsel ')),squeeze(audiodata.data(' newsel ')))'])
                        end
                        %eval([chartfunc '(hsub,Xdata,squeeze(audiodata.data(' newsel ')))'])
                        xlabel(hsub,strrep([tabledata{axis(1,1),1} ' [' audiodata.(genvarname([tabledata{axis(1,1),1} 'info'])).units ']'],'_',' '))
                        ylabel(hsub,strrep(audiodata.datainfo.units,'_',' '))
                        if isnumeric(audiodata.(genvarname(cats{1,1}))), audiodata.(genvarname(cats{1,1})) = num2cell(audiodata.(genvarname(cats{1,1}))); end
                        dim2str = audiodata.(genvarname(cats{1,1})){dim2(i)};
                        if isnumeric(dim2str), dim2str = num2str(dim2str); end
                        if isnumeric(audiodata.(genvarname(cats{2,1}))), audiodata.(genvarname(cats{2,1})) = num2cell(audiodata.(genvarname(cats{2,1}))); end
                        dim3str = audiodata.(genvarname(cats{2,1})){dim3(j)};
                        if isnumeric(dim3str), dim3str = num2str(dim3str); end
                        title(hsub,[dim2str ',' dim3str])
                        k = k + 1;
                    end
                end
            case ''
                set(handles.cattable,'Data',handles.tabledata)
                catorcont = handles.tabledata(:,4);
                if any(cellfun(@isempty,catorcont)), catorcont(cellfun(@isempty,catorcont)) = {false}; end
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
                warndlg(error.message,'AARAE info','modal')
        end
    elseif length(find(cellfun(@length,cellfun(@str2num,tabledata(cell2mat(catorcont) == 0,2),'UniformOutput',false)) ~= 1)) == 1 && (strcmp(chartfunc,'mesh') || strcmp(chartfunc,'surf') || strcmp(chartfunc,'imagesc'))
    % 3D data subplotting
        selections = cellfun(@str2num,tabledata(:,2),'UniformOutput',false);
        if any(cellfun(@isempty,selections)), selections(cellfun(@isempty,selections)) = {1}; end
        cats = tabledata(cellfun(@length,selections) ~= 1,1);
        choice = questdlg('Would you like to display data subplots?',...
                          'AARAE subplots',...
                          'Yes','No',...
                          'Yes');
        switch choice
            case 'Yes'
                dims = cellfun(@str2num,tabledata(cellfun(@length,selections) ~= 1,2),'UniformOutput',false);
                indexes = find(cellfun(@length,selections) ~= 1);
                dim2 = dims{1,1};
                figure;
                [~,~,ext] = fileparts(handles.Settings.colormap);
                if isempty(ext)
                    eval(['cmap = colormap(' lower(handles.Settings.colormap) '(128));']);
                else
                    custom_colormap = importdata([cd '/Utilities/Custom_colormaps/' handles.Settings.colormap]);
                    cmap = colormap(custom_colormap);
                end
                set(get(haxes,'Parent'),'DefaultAxesColorOrder',cmap)
                k = 1;
                [r,c] = subplotpositions(length(dim2), 0.5);
                for i = 1:length(dim2)
                    hsub = subplot(r,c,k);
                    newsel = strsplit(sel,',');
                    newsel{1,indexes(1,1)} = num2str(i);
                    newsel = strjoin(newsel,',');
                    eval(['data = squeeze(audiodata.data(' newsel '));'])
                    if ~isequal([length(Ydata),length(Xdata)],size(data)), data = data'; end
                    if ~strcmp(chartfunc,'imagesc')
                        eval([chartfunc '(hsub,Xdata,Ydata,data)'])
                        %colormap(aaraecmap)
                    else
                        eval([chartfunc '(Xdata,1:length(Ydata),data,''Parent'',hsub)'])
                        set(hsub,'YTickLabel',num2str(Ydata'))
                        set(hsub,'YDir','normal')
                        %colormap(aaraecmap)
                    end
                    xlabel(hsub,strrep([tabledata{axis(1,1),1} ' [' audiodata.(genvarname([tabledata{axis(1,1),1} 'info'])).units ']'],'_',' '))
                    ylabel(hsub,strrep([tabledata{axis(1,2),1} ' [' audiodata.(genvarname([tabledata{axis(1,2),1} 'info'])).units ']'],'_',' '))
                    zlabel(hsub,strrep(audiodata.datainfo.units,'_',' '))
                    if isnumeric(audiodata.(genvarname(cats{1,1}))), audiodata.(genvarname(cats{1,1})) = num2cell(audiodata.(genvarname(cats{1,1}))); end
                    dim2str = audiodata.(genvarname(cats{1,1})){dim2(i)};
                    if isnumeric(dim2str), dim2str = num2str(dim2str); end
                    title(hsub,dim2str)
                    k = k + 1;
                end
            otherwise
                set(handles.cattable,'Data',handles.tabledata)
                catorcont = handles.tabledata(:,4);
                if any(cellfun(@isempty,catorcont)), catorcont(cellfun(@isempty,catorcont)) = {false}; end
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
                warndlg(error.message,'AARAE info','modal')
        end
    else
        set(handles.cattable,'Data',handles.tabledata)
        catorcont = handles.tabledata(:,4);
        if any(cellfun(@isempty,catorcont)), catorcont(cellfun(@isempty,catorcont)) = {false}; end
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
        doresultplot(handles,haxes);
        warndlg(error.message,'AARAE info','modal')
    end
end