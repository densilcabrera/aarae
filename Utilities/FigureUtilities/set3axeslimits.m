function set3axeslimits
% For the 3-axis audio comparison scatter plots created by AARAE's '=' button
compplot = gcf;

iplots = findobj(compplot,'Type','axes');
xlims = cell2mat(get(iplots,'Xlim'));
ylims = cell2mat(get(iplots,'Ylim'));
Levellims = ylims(2,:);
Timelims = xlims(4,:);
Freqlims = xlims(2,:);
prompt = {'Level min','Level max','Freq min',...
    'Freq max','Time min','Time max'};
dlg_title = 'Axes limits';
num_lines = 1;
def = {num2str(Levellims(:,1)),num2str(Levellims(:,2)),...
    num2str(Freqlims(:,1)),num2str(Freqlims(:,2)),...
    num2str(Timelims(:,1)),num2str(Timelims(:,2))};
answer = inputdlg(prompt,dlg_title,num_lines,def);
if ~isempty(answer)
    Lmin = str2num(answer{1});
    Lmax = str2num(answer{2});
    fmin = str2num(answer{3});
    fmax = str2num(answer{4});
    tmin = str2num(answer{5});
    tmax = str2num(answer{6});
    if ~isempty(Lmin) && ~isempty(Lmin) && ~isempty(fmin)...
            && ~isempty(fmin) && ~isempty(tmin) && ~isempty(tmax)
        set(iplots(2),'Xlim',[fmin,fmax])
        set(iplots(3),'Xlim',[tmin,tmax])
        set(iplots(4),'Xlim',[tmin,tmax])
        set(iplots(2),'Ylim',[Lmin,Lmax])
        set(iplots(3),'Ylim',[Lmin,Lmax])
        set(iplots(4),'Ylim',[fmin,fmax])
    else
        warndlg('Invalid entry','AARAE info')
    end
end
