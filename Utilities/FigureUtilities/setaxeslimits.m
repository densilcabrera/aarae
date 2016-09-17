function setaxeslimits
% For the main audio comparison plots created by AARAE's '=' button
compplot = gcf;
iplots = findobj(compplot,'Type','axes');
xlims = cell2mat(get(iplots,'Xlim'));
ylims = cell2mat(get(iplots,'Ylim'));
prompt = {'Xmin','Xmax','X-scale: change between lin and log [0 | 1]',...
    'Ymin','Ymax','Y-scale: change between lin and log [0 | 1]'};
dlg_title = 'Axes limits';
num_lines = 1;
def = {num2str(min(xlims(:,1))),num2str(max(xlims(:,2))),'0',num2str(min(ylims(:,1))),num2str(max(ylims(:,2))),'0'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
if ~isempty(answer)
    xmin = str2num(answer{1});
    xmax = str2num(answer{2});
    xlinlog = str2num(answer{3});
    ymin = str2num(answer{4});
    ymax = str2num(answer{5});
    ylinlog = str2num(answer{6});
    if ~isempty(xmin) && ~isempty(xmax) && ~isempty(xlinlog)...
            && ~isempty(ymin) && ~isempty(ymax) && ~isempty(ylinlog)
        set(iplots,'Xlim',[xmin xmax])
        set(iplots,'Ylim',[ymin ymax])
        if xlinlog == 1
            xscale = get(iplots(1),'XScale');
            if strcmp(xscale,'linear')
                set(iplots,'XScale','log','XTickLabelMode','auto')
            else
                set(iplots,'XScale','linear')
            end
        end
        if ylinlog == 1
            yscale = get(iplots(1),'YScale');
            if strcmp(yscale,'linear')
                set(iplots,'YScale','log','YTickLabelMode','auto')
            else
                set(iplots,'YScale','linear')
            end
        end
    else
        warndlg('Invalid entry','AARAE info')
    end
end
