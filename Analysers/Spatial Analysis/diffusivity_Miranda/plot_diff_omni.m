function plot_diff_omni(diffusivity, omniparameters, plot_setup, FileName, PathName)

if plot_setup.rt == 1;
    figure1 = figure;

    Tickarray = [];

    for I = 1:length(omniparameters.bands);
        labelarray{I} = num2str(omniparameters.bands(I));
        Tickarray = cat(2,Tickarray,I);
    end

    % Create axes
    axes1 = axes('Parent',figure1,'YGrid','on','XTickLabel',labelarray,...
        'XTick',Tickarray,...
        'FontWeight','bold',...
        'FontSize',16);
    box(axes1,'on');
    hold(axes1,'all');

    % Create multiple lines using matrix input to bar
    bar1 = bar([omniparameters.T30;omniparameters.T20]');
    set(bar1(1),'FaceColor',[0 0 0],'DisplayName','T30');
    set(bar1(2),'DisplayName','T20');

    % Create xlabel
    xlabel('Frequencies','FontWeight','bold','FontSize',16);

    % Create ylabel
    ylabel('Seconds','FontWeight','bold','FontSize',16);

    % Create title
    title('Reverb Time','FontWeight','bold','FontSize',16);

    legend(axes1,'show');
    set(legend,'Location','NorthEastOutside');

end



saveas(gcf, strcat(PathName, 'RT', FileName), 'fig')
saveas(gcf, strcat(PathName, 'RT', FileName), 'jpg')



if plot_setup.st == 1;
    figure1 = figure;

    Tickarray = [];

    for I = 1:length(omniparameters.bands)+1;
        if I <= length(omniparameters.bands); 
            labelarray{I} = num2str(omniparameters.bands(I));
            Tickarray = cat(2,Tickarray,I);
        else
            labelarray{I} = 'Mean';
            Tickarray = cat(2,Tickarray,I);
        end
    end

    % Create axes
    axes1 = axes('Parent',figure1,'YGrid','on','XTickLabel',labelarray,...
        'XTick',Tickarray,...
        'FontWeight','bold',...
        'FontSize',16);
    box(axes1,'on');
    hold(axes1,'all');

    % Create multiple lines using matrix input to bar
    bar1 = bar([omniparameters.STearly_band,omniparameters.STearly;omniparameters.STlate_band,omniparameters.STlate]');
    set(bar1(1),'FaceColor',[0 0 0],'DisplayName','STearly');
    set(bar1(2),'DisplayName','STlate');

    % Create xlabel
    xlabel('Frequencies','FontWeight','bold','FontSize',16);

    % Create ylabel
    ylabel('dB','FontWeight','bold','FontSize',16);

    % Create title
    title('Stage Support','FontWeight','bold','FontSize',16);

    legend(axes1,'show');
    set(legend,'Location','NorthEastOutside');

end

saveas(gcf, strcat(PathName, 'ST', FileName), 'jpg')
saveas(gcf, strcat(PathName, 'ST', FileName), 'fig')



[~, diff_win_time] = max(omniparameters.magnitudedB);

magnitudedBforplot = omniparameters.magnitudedB(diff_win_time:end);

diffusivityforplot = diffusivity.diff_wideband(diff_win_time:end);

analysis_times_for_plot = omniparameters.analysis_times(diff_win_time:end);
analysis_times_for_plot = analysis_times_for_plot - min(analysis_times_for_plot);

[~, sampleindex_for_diff_plot] = min(abs((analysis_times_for_plot - plot_setup.diffusivityplot_length)));

[fitresult, ~] = createFit(analysis_times_for_plot(1:sampleindex_for_diff_plot), diffusivityforplot(1:sampleindex_for_diff_plot));

h = gcf;

close (h);

fitcurve = feval(fitresult, analysis_times_for_plot(1:sampleindex_for_diff_plot));

figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,...
    'YTickLabel',{'0','0.2','0.4','0.6','0.8','1'},...
    'YTick',[0 0.2 0.4 0.6 0.8 1],...
    'YGrid','on',...
    'XGrid','on',...
    'GridLineStyle','--',...
    'FontWeight','bold',...
    'FontSize',16);
xlim(axes1,[0 plot_setup.diffusivityplot_length]);
ylim(axes1,[0 1]);
box(axes1,'on');
hold(axes1,'all');

% Create plot
plot1 = plot(analysis_times_for_plot(1:sampleindex_for_diff_plot),[diffusivityforplot(1:sampleindex_for_diff_plot) fitcurve],'Parent',axes1);
set(plot1(1),'LineWidth',2,'Color',[0 0 0]);
set(plot1(2),'LineWidth',1.5,'LineStyle',':',...
    'Color',[0.24705882370472 0.24705882370472 0.24705882370472]);


% Create xlabel
xlabel('Time (s)','FontWeight','bold','FontSize',16);

% Create ylabel
ylabel('Diffusivity','FontWeight','bold','FontSize',16);

% Create axes
axes2 = axes('Parent',figure1,'YTick',[-50 -40 -30 -20 -10 0],...
    'YAxisLocation','right',...
    'YColor',[1 0 0],...
    'GridLineStyle','--',...
    'FontWeight','bold',...
    'FontSize',16,...
    'ColorOrder',[0 0.5 0;1 0 0;0 0.75 0.75;0.75 0 0.75;0.75 0.75 0;0.25 0.25 0.25;0 0 1],...
    'Color','none');
xlim(axes2,[0 plot_setup.diffusivityplot_length]);
ylim(axes2,[-50 0]);
hold(axes2,'all');

% Create plot
plot(analysis_times_for_plot,magnitudedBforplot,'Parent',axes2,'LineWidth',1.5,'LineStyle','--','Color',[1 0 0]);

% Create ylabel
ylabel('Level (dB)','VerticalAlignment','cap','FontWeight','bold',...
    'FontSize',16,...
    'Color',[1 0 0]);
h=gcf;

set(h,'Position',[1 1 758 504]);

saveas(gcf, strcat(PathName, 'Diff', FileName), 'jpg')
saveas(gcf, strcat(PathName, 'Diff', FileName), 'fig')

