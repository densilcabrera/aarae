function setXaxislimits
% This function can be used as a figure callback when subplots are used (so
% that all subplots can be rescaled identically).
compplot = gcf;
iplots = findobj(compplot,'Type','axes');
xlims = cell2mat(get(iplots,'Xlim'));
prompt = {'Xmin','Xmax','X-scale: change between lin and log [0 | 1]'};
dlg_title = 'X axis limits';
num_lines = 1;
def = {num2str(min(xlims(:,1))),num2str(max(xlims(:,2))),'0'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
if ~isempty(answer)
    xmin = str2num(answer{1});
    xmax = str2num(answer{2});
    xlinlog = str2num(answer{3});
    if ~isempty(xmin) && ~isempty(xmax) && ~isempty(xlinlog)
        set(iplots,'Ylim',[xmin xmax])
        if xlinlog == 1
            xscale = get(iplots(1),'XScale');
            if strcmp(xscale,'linear')
                set(iplots,'XScale','log','XTickLabelMode','auto')
            else
                set(iplots,'XScale','linear')
            end
        end
    end
else
    warndlg('Invalid entry','AARAE info')
end
end
