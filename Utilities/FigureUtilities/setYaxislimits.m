function setYaxislimits
% This function can be used as a figure callback when subplots are used (so
% that all subplots can be rescaled identically).
compplot = gcf;
iplots = findobj(compplot,'Type','axes');
ylims = cell2mat(get(iplots,'Ylim'));
prompt = {'Ymin','Ymax','Y-scale: change between lin and log [0 | 1]'};
dlg_title = 'Y axis limits';
num_lines = 1;
def = {num2str(min(ylims(:,1))),num2str(max(ylims(:,2))),'0'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
if ~isempty(answer)
    ymin = str2num(answer{1});
    ymax = str2num(answer{2});
    ylinlog = str2num(answer{3});
    if ~isempty(ymin) && ~isempty(ymax) && ~isempty(ylinlog)
        set(iplots,'Ylim',[ymin ymax])
        if ylinlog == 1
            yscale = get(iplots(1),'YScale');
            if strcmp(yscale,'linear')
                set(iplots,'YScale','log','YTickLabelMode','auto')
            else
                set(iplots,'YScale','linear')
            end
        end
    end
else
    warndlg('Invalid entry','AARAE info')
end
end
