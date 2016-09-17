function OUT = circular_rigid_piston_infbaffle(a,c,rho,v0,method,val1,val2,val3,val4,val5)
% This function calculates the sound radiation from a rigid circular piston
% in an infinite baffle.
%
% SOURCES:
% Morse, P. M., & Ingard, K. U. (1968). Theoretical acoustics. Princeton
% University Press.
%
% Pierce, A. D. (1981). Acoustics: an introduction to its physical
% principles and applications (Vol. 20). New York: McGraw-Hill.
%
% Jacobsen, F., Poulsen, T., Rindel, J. H., Gade, A. C., & Ohlrich, M.
% (2011). Fundamentals of acoustics and noise control. Department of
% Electrical Engineering, Technical University of Denmark.
%
% Jacobsen, F., & Juhl, P. (2010). Radiation of sound. DTU, Lyngby.


if nargin == 0
    str = {'1. One frequency only'
        '2. Full spectrum on-axis';
        '3. Full spectrum around arc';
        '4. Full spectrum at one arbitrary location (IR output only)';
        '5. Radiation summary parameters (e.g. impedance)'};
    % 6. Visualisation of filtered pulse propagation in space??????
    [method,ok] = listdlg('PromptString','Select the calculation method',...
        'SelectionMode','single',...
        'ListString',str,...
        'ListSize', [400,400]);
    if ~ok
        OUT = [];
        return
    end
end

switch method
    case 1 % ONE FREQUENCY ONLY
        if nargin == 0
            param = inputdlg({'Frequency to plot (Hz)';...
                'Speed of sound in medium (m/s)';...
                'Density of medium (kg/m^3)';...
                'Volume velocity of source (m^3/s)';...
                'Piston radius (m)';...
                'Maximum distance (m)';...
                'On-axis divisions per metre';...
                '2-d wave field divisions per metre';
                'Far-field analytic approximation [0] or numeric calculation (slower!) [1]'},...
                'Circular Piston Radiation Input Parameters',...
                [1 80],...
                {'1000', '342', '1.225', '0.0001', '0.2', '2', '500','100','0'});
            if length(param) < 9, param = []; end
            if ~isempty(param)
                [val1, f] = deal(str2num(char(param(1))));
                c = str2num(char(param(2)));
                rho = str2num(char(param(3)));
                v0 = str2num(char(param(4)));
                a = str2num(char(param(5)));
                [val2, size_of_calc] = deal(str2num(char(param(6))));
                [val3, divisions_per_m] = deal(str2num(char(param(7))));
                [val4, griddivisions_per_m] = deal(str2num(char(param(8))));
                [val5, numeric] = deal(str2num(char(param(9))));
            else
                OUT = [];
                return
            end
        end
        % *****************************************************************
        % ON-AXIS PLOT OF SOUND PRESSURE, PARTICLE VELOCITY AND
        % INTENSITY
        t=0; % time
        omega = 2.*pi.*f; % angular frequency
        k = (omega./c); % wave number
        X = 0:1/divisions_per_m:size_of_calc; % distance
        beta = X ./ (X.^2 + a.^2).^0.5;
        gamma = k./2 .*((X.^2 + a.^2).^0.5 - X);
        p = rho.*c.*v0.*exp(1i*(omega.*t-k.*X)).*(1-exp(-2i*gamma));
        u = v0.*exp(1i*(omega.*t-k.*X)).*(1-beta.*exp(-2i*gamma));
        Iactive = (p.*conj(u));
        Ireact = abs(imag(conj(u).*p));
        
        figure1 = figure('Name','Circular piston in infinite baffle');
        subplot1 = subplot(2,1,1,'Parent',figure1);
        hold(subplot1,'on');
        plot(X,mag2db(abs(p)./2e-5),'r','DisplayName','Pressure')
        plot(X,mag2db(abs(u)./5e-8),'b','DisplayName','Particle velocity')
        plot(X,pow2db(abs(Iactive)./1e-12),'k','DisplayName','Intensity (active)')
        plot(X,pow2db(abs(Ireact)./1e-12),'Color',[0.5 0.7 0.5],'DisplayName','Intensity (reactive)')
        xlabel('Distance (m)')
        ylabel('Level (dB)')
        legend('show')
        title([num2str(f) ' Hz, source radius = ' num2str(a) ' m, volume velocity = ' num2str(v0) ' m^3/s, surface velocity = ' num2str(v0./(pi*a.^2)) ' m/s'])
        
        subplot2=subplot(2,1,2,'Parent',figure1);
        hold(subplot2,'on');
        plot(X,180.*angle(p)./pi,'r','DisplayName','Pressure')
        plot(X,180.*angle(u)./pi,'b','DisplayName','Particle velocity')
        dp = 180.*(angle(p)-angle(u))./pi;
        dp(dp<-180) = dp(dp<-180)+180;
        plot(X,dp,':k','LineWidth',1,...
            'DisplayName','P-U phase difference')
        xlabel('Distance (m)')
        ylabel('Phase (deg)')
        ylim([-180 180])
        set(subplot2,'YTick',[-180 -90 0 90 180]);
        legend('show')
        
        
        % *****************************************************************
        % PRESSURE SOUND FIELD SURFACE CHART
        if numeric == 1
            % CREATE SOURCE ARRAY (upper half)
            subsourcespacing = 0.125*c/24000; % 1/8 wavelength @ Nyquist freq
            Ys=-a:subsourcespacing:a;
            Yslen = length(Ys);
            Zs = (0:subsourcespacing:a)';
            Zslen = length(Zs);
            Ys = repmat(Ys,[Zslen,1]);
            Zs = repmat(Zs,[1,Yslen]);
            semicircle = abs(Ys+1i*Zs)<= a;
            Ys = Ys(semicircle); % Y subsource coordinates
            Zs = Zs(semicircle); % Z subsource coordinates
            N = length(Ys);
            
            % DATA AT ONE FREQUENCY OVER THE SPACE
            omega = 2.*pi.*f; % angular frequency
            k = (omega./c);
            X = 0:1/griddivisions_per_m:size_of_calc;
            X = repmat(X,[length(X),1]);
            Y = X';
            p = zeros(size(Y));
            Qs = subsourcespacing.^2 .* v0 ./ (pi*a.^2);
            % add monopoles in half-space
            % double values except when Z == 0 to make semicircle into
            % circle
            for J = 1:N
                r = (X.^2 + (Ys(J)-Y).^2 + Zs(J).^2).^0.5;
                if Zs(J) == 0
                    p = p + 2.*((-1i.*rho.*c.*k)./(4.*pi)).*Qs.*(exp(1i.*k.*r)./r);
                else
                    p = p + 4.*((-1i.*rho.*c.*k)./(4.*pi)).*Qs.*(exp(1i.*k.*r)./r);
                end
            end
        else % analytic far-field approximation
            omega = 2.*pi.*f; % angular frequency
            k = (omega./c);
            X = 0:1/griddivisions_per_m:size_of_calc;
            X = repmat(X,[length(X),1]);
            Y = X';
            r = (X.^2 + Y.^2).^0.5; % distance from source centre
            theta = atan2(Y,X); % angle from source centre
            p = zeros(size(Y));
            p(theta==0)=-1i.*omega.*rho.*c.*v0.*(a.^2).*0.5.*(exp(1i.*k.*r(theta==0))./r(theta==0));
            p(theta~=0)=-1i.*omega.*rho.*c.*v0.*(a.^2).*(besselj(1,(k.*a.*sin(theta(theta~=0))))./(k.*a.*sin(theta(theta~=0)))).*(exp(1i.*k.*r(theta~=0))./r(theta~=0));
        end
        % level and phase 3D mesh figure
        xsize = 1024;
        x = (0:xsize)'./xsize;
        Sat = abs(sin(2*pi*x));
        Hue = circshift(x,-round(xsize*30/360));
        Val = abs(x-0.5).*2;
        circularcolormap = hsv2rgb([Hue Sat Val]);
        fig0 = figure('ColorMap',circularcolormap,...
            'Name', 'Circular piston in infinite baffle');
        ax0 = axes('Parent',fig0,'CLim', [-180 180]);
        
        ax0.DataAspectRatio = [1 1 1];
        ax0.PlotBoxAspectRatio = [1 1 1];
        
        hold(ax0,'on');
        h=surf(pow2db(abs(p)./2e-5),'Parent',ax0);
        h.CData = 180.*angle(p)./pi;
        h.FaceColor = 'flat'; % 'interp' does not work for circular colormaps
        h.EdgeLighting = 'gouraud';
        h.FaceAlpha = 1;
        h.FaceLighting = 'gouraud';
        h.AmbientStrength = 0.6;
        h.DiffuseStrength = 1;
        h.EdgeAlpha = 0.1;
        title([num2str(f) ' Hz, source radius = ' num2str(a) ' m'])
        camlight left
        %     set(gca, 'DataAspectRatio', [1 1 1])
        %     set(gca, 'PlotBoxAspectRatio', [1 1 1])
        view(ax0,[76.8 50.8]);
        grid(ax0,'on');
        xlabel('Axial distance (m)');
        ylabel('Lateral distance (m)');
        zlabel('Sound pressure level (dB)');
        %zlim([max(max(p_db))-40 max(max(p_db))]);
        colorbar('peer',ax0,'Ticks',[-180 -90 0 90 180],...
            'TickLabels',{'-180','-90','0 deg','90','180'});
        hold off
        
        % *****************************************************************
        % INTENSITY VECTOR CHART (IF NUMERIC)
        if numeric == 1
            X = 0:size_of_calc/50:size_of_calc;
            X = repmat(X,[length(X),1]);
            Y = X';
            [p,u,v] = deal(zeros(size(Y)));
            %Qs = subsourcespacing.^2 .* v0 ./ (pi*a.^2);
            % add monopoles in half-space
            % double values except when Z == 0 to make semicircle into
            % circle
            for J = 1:N
                r = (X.^2 + (Ys(J)-Y).^2 + Zs(J).^2).^0.5;
                if Zs(J) == 0
                    p1 = 2.*((-1i.*rho.*c.*k)./(4.*pi)).*Qs.*(exp(1i.*k.*r)./r);
                    p = p + p1;
                    %verticalangle = 0;
                    horizontalangle = atan2(repmat(Ys(J),size(X)),X);
                    u1 = p1./(rho*c).*(1+1./(1i.*k.*r));
                    u = u + u1 .*cos(horizontalangle); % axial velocity
                    v = v + u1 .*sin(horizontalangle); % lateral velocity
                else
                    p1 = 4.*((-1i.*rho.*c.*k)./(4.*pi)).*Qs.*(exp(1i.*k.*r)./r);
                    p = p + p1;
                    verticalangle = atan2(repmat(Zs(J),size(X)),X);
                    horizontalangle = atan2(repmat(Ys(J),size(X)),X);
                    u1 = p1./(rho*c).*(1+1./(1i.*k.*r));
                    u = u + u1 .* cos(verticalangle).*cos(horizontalangle); % axial velocity
                    v = v + u1 .* cos(verticalangle).*sin(horizontalangle); % lateral velocity
                end
            end
            figure('Name','Circular piston in infinite baffle')
            %quiver(X,Y,abs(u),abs(v))
            streamslice(X,Y,abs(u),abs(v))
        end
        
        % *****************************************************************
        % PRESSURE FAR-FIELD DIRECTIVITY CHART
        step = 0.1; % resolution in degrees
        theta = pi.*(0:step:90)./180;
        step = pi*step/180; % step in rad
        r = 1000; % large distance!
        p = zeros(size(theta));
        p(1) = -1i.*omega.*rho.*c.*v0.*(a.^2).*0.5.*(exp(1i.*k.*r)./r);
        p(2:end) = -1i.*omega.*rho.*c.*v0.*(a.^2).*(besselj(1,(k.*a.*sin(theta(2:end))))./(k.*a.*sin(theta(2:end)))).*(exp(1i.*k.*r)./r);
        p = abs(p)./abs(p(1)); % normalise to on-axis
        L = mag2db(p); % in dB
        Lmin = -40;
        L(L<Lmin)=Lmin;
        Q = 4 * p(1).^2 ./ ...
            (p(1).^2 .* sin(theta(1)-step/2) + ...
            sum(p(2:end-1).^2 .* (sin(theta(2:end-1)+step/2)-sin(theta(2:end-1)-step/2))) + ...
            p(end).^2 .* sin(theta(end)+step/2));
        
        % polar plot
        try % the function polarplot was introduced to Matlab in 2016
            polarfig = figure('Name','Sound level re 0 deg');
            polarplot([theta pi+1e-99 3*pi/2-1e99 theta+3*pi/2],...
                [L Lmin Lmin flip(L)],...
                'Color', [0 0.7 0],'LineWidth',1)
            rlim([10*floor(Lmin/10) max(L)])
        catch
            delete(polarfig)
            figure('Name','Sound magnitude re 0 deg');
            p=polar([theta pi+1e-99 3*pi/2-1e99 theta+3*pi/2],...
                [p 0 0 flip(p)]);
            set(p,'Color',[0 0.7 0])
        end
        title([num2str(f) ' Hz, Source radius: ' num2str(a) ' m, Directivity factor: ' num2str(Q) ', DI: ' num2str(pow2db(Q)) ' dB'])
        OUT = [];
        
        
        
    case 2 % ON-AXIS FULL SPECTRUM
        if nargin == 0
            param = inputdlg({'Speed of sound in medium (m/s)';...
                'Density of medium (kg/m^3)';...
                'Volume velocity of source (m^3/s)';...
                'Piston radius (m)';...
                'Maximum distance (m)';...
                'Chart data divisions per metre';...
                'Generate IRs? No [0], Pressure [1], Velocity [2]';
                'IR divisions per metre'},...
                'Circular Piston Radiation Input Parameters',...
                [1 80],...
                {'342', '1.225', '0.0001', '0.2', '2', '100','1','2'});
            if length(param) < 8, param = []; end
            if ~isempty(param)
                c = str2num(char(param(1)));
                rho = str2num(char(param(2)));
                v0 = str2num(char(param(3)));
                a = str2num(char(param(4)));
                [val1, size_of_calc] = deal(str2num(char(param(5))));
                [val2, divisions_per_m] = deal(str2num(char(param(6))));
                [val3, makeIR] = deal(str2num(char(param(7))));
                [val4, IRdivisions_per_m] = deal(str2num(char(param(8))));
                val5 = [];
            else
                OUT = [];
                return
            end
        end
        
        % Spatial decay rate analysis
        X = 10.^((-20:0.25:20)./10);
        farray = [16; 31.5; 63; 125; 250; 500; 1000; 2000; 4000; 8000; 16000; 31500];
        omegaarray = 2.*pi.*farray;
        karray = omegaarray./c;
        beta = X ./ (X.^2 + a.^2).^0.5;
        gamma = repmat(karray,[1,length(X)])./2 ...
            .*((repmat(X,[length(karray),1]).^2 + a.^2).^0.5...
            - repmat(X,[length(karray),1]));
        Pspectrum = rho.*c.*v0.*...
            exp(1i*(-repmat(karray,[1,length(X)])...
            .*repmat(X,[length(karray),1])))...
            .*(1-exp(-2i*gamma));
        Uspectrum = v0.*exp(1i*(-repmat(karray,[1,length(X)])...
            .*repmat(X,[length(karray),1])))...
            .*(1-repmat(beta,[length(karray),1])...
            .*exp(-2i*gamma));
        Lp = mag2db(abs(Pspectrum)./2e-5);
        Lu = mag2db(abs(Uspectrum)./5e-8);
        LIact = pow2db(real(conj(Uspectrum).*Pspectrum)./1e-12);
        LIreact = pow2db(abs(imag(conj(Uspectrum).*Pspectrum)./1e-12));

        
        
        % Spatial decay rate impedance & inv sq deviation plot
        figure('Name','On-axis impedance and deviation from inverse (square) dispersion 0.01 m to 100 m');
        for fcount = 1:12
            axes1 = subplot(3,4,fcount);
            semilogx(X,mag2db(real(Pspectrum(fcount,:)./Uspectrum(fcount,:))./(rho*c)),'g','LineWidth',1,'DisplayName','Re(Z)/(rho.c)');
            hold on
            semilogx(X,mag2db(imag(Pspectrum(fcount,:)./Uspectrum(fcount,:))./(rho*c)),'c','LineWidth',1,'DisplayName','Im(Z)/(rho.c)');
            invsq = mag2db(max(X)./X)-Lu(fcount,:)+Lu(fcount,end);
            semilogx(X,-invsq,'b','LineWidth',1,'DisplayName','Lu deviation');
            invsq = mag2db(max(X)./X)-LIact(fcount,:)+LIact(fcount,end);
            semilogx(X,-invsq,'k','LineWidth',1,'DisplayName','LI active dev');
            invsq = mag2db(max(X)./X)-Lp(fcount,:)+Lp(fcount,end);
            semilogx(X,-invsq,'r','LineWidth',1,'DisplayName','Lp deviation');
            ylim([-50 30])
            title([num2str(farray(fcount)) ' Hz'])
            if fcount>=9
                xlabel('Distance (m)')
            end
            if mod(fcount,4)==1
                ylabel('Value (dB)')
            end
            if fcount==12;
                legend('show','Location','southeast');
            end
            set(axes1,'XGrid','on','XMinorTick','on','XScale','log','YGrid','on')
        end
                
        
        % Spatial decay rate sound level plot
        figure('Name','On-axis spatial decay 0.01 m to 100 m');
        for fcount = 1:12
            axes1 = subplot(3,4,fcount);
            h(1)=semilogx(X,Lu(fcount,:),'b','LineWidth',1,'DisplayName','Lu');
            hold on
            h(2)=semilogx(X,LIact(fcount,:),'k','LineWidth',1,'DisplayName','LI active');
            h(3)=semilogx(X,LIreact(fcount,:),'Color',[0.5 0.7 0.5],'LineWidth',1,'DisplayName','LI reactive');
            h(4)=semilogx(X,Lp(fcount,:),'r','LineWidth',1,'DisplayName','Lp');
            offset = 80+10.*round(min(min(Lp))./10);
            for lineoffsset = -300:10:100
                semilogx(X,offset+lineoffsset+mag2db(1./(X)),'Color',[0.8 0.8 0.8],'LineWidth',0.5)
            end
            miny = min(min([Lp Lu]));
            maxy = max(max([Lp Lu]));
            ylim([miny maxy])
            title([num2str(farray(fcount)) ' Hz'])
            if fcount>=9
                xlabel('Distance (m)')
            end
            if mod(fcount,4)==1
                ylabel('Level (dB)')
            end
            if fcount==12;
                legend(h,'Location','southwest');
            end
            set(axes1,'XGrid','on','XMinorTick','on','XScale','log','YGrid','on')
        end
        

        
        

        
        % main analyis
        X = 0:1/divisions_per_m:size_of_calc;
        farray = (0:20:24000)';
        omegaarray = 2.*pi.*farray;
        karray = omegaarray./c;
        beta = X ./ (X.^2 + a.^2).^0.5;
        gamma = repmat(karray,[1,length(X)])./2 ...
            .*((repmat(X,[length(karray),1]).^2 + a.^2).^0.5...
            - repmat(X,[length(karray),1]));
        Pspectrum = rho.*c.*v0.*...
            exp(1i*(-repmat(karray,[1,length(X)])...
            .*repmat(X,[length(karray),1])))...
            .*(1-exp(-2i*gamma));
        Uspectrum = v0.*exp(1i*(-repmat(karray,[1,length(X)])...
            .*repmat(X,[length(karray),1])))...
            .*(1-repmat(beta,[length(karray),1])...
            .*exp(-2i*gamma));
        phasediff = 180.*(mod(pi+angle(Pspectrum)-angle(Uspectrum),2*pi)-pi)./pi;
        
        % level and phase 3D mesh figure
        %HSVmap = hsv(1024);
        phasediffmap = [(1:-1/1023:0)', (0:1/1023:1)', [(0:1/1023:0.5)';(0.5:-1/1023:0)']];
        fig0 = figure('Name','On-Axis Pressure','Colormap',phasediffmap);
        ax0 = axes('Parent',fig0,'CLim', [-90 90]);
        hold(ax0,'on');
        h=surf(X,farray,mag2db(abs(Pspectrum)./2e-5),'Parent',ax0);
        h.CData = phasediff;
        h.FaceColor = 'flat'; % 'interp' does not work for circular colormaps
        h.EdgeLighting = 'gouraud';
        h.FaceAlpha = 1;
        h.FaceLighting = 'gouraud';
        h.AmbientStrength = 0.6;
        h.DiffuseStrength = 1;
        h.EdgeAlpha = 0.1;
        zlim([max(max(mag2db(abs(Pspectrum)./2e-5)))-40, max(max(mag2db(abs(Pspectrum)./2e-5)))]);
        xlabel('Distance (m)');
        ylabel('Frequency (Hz)');
        zlabel('Sound Pressure Level (dB)');
        title(['Source radius = ' num2str(a) ' m (colormap shows p-u phase difference)'])
        camlight left
        view(ax0,[-26.8 61.2]);
        grid(ax0,'on');
        colorbar('peer',ax0,'Ticks',[-180 -90 0 90 180],...
            'TickLabels',{'-180','-90','0 deg','90','180'});
        hold off
        
        
        fig0 = figure('Name','On-Axis Particle Velocity','Colormap',phasediffmap);
        ax0 = axes('Parent',fig0,'CLim', [-90 90]);
        hold(ax0,'on');
        h=surf(X,farray,mag2db(abs(Uspectrum)./5e-8),'Parent',ax0);
        h.CData = phasediff;
        h.FaceColor = 'flat'; % 'interp' does not work for circular colormaps
        h.EdgeLighting = 'gouraud';
        h.FaceAlpha = 1;
        h.FaceLighting = 'gouraud';
        h.AmbientStrength = 0.6;
        h.DiffuseStrength = 1;
        h.EdgeAlpha = 0.1;
        zlim([max(max(mag2db(abs(Uspectrum)./5e-8)))-40, max(max(mag2db(abs(Uspectrum)./5e-8)))]);
        xlabel('Distance (m)');
        ylabel('Frequency (Hz)');
        zlabel('Particle Velocity Level (dB)');
        title(['Source radius = ' num2str(a) ' m (colormap shows p-u phase difference)'])
        camlight left
        view(ax0,[-26.8 61.2]);
        grid(ax0,'on');
        colorbar('peer',ax0,'Ticks',[-180 -90 0 90 180],...
            'TickLabels',{'-180','-90','0 deg','90','180'});
        hold off
        
        % plot of pressure, intensity and particle velocity at furthest
        % distance
        figure1 = figure('Name','Rigid circular piston in infinite baffle on-axis');
        axes1 = subplot(3,1,1:2,'Parent',figure1);
        %axes1 = axes('Parent',figure1);
        hold(axes1,'on');
        semilogx1 = semilogx(farray,...
            [mag2db(abs(Pspectrum(:,end))./2e-5),...
            mag2db(abs(Uspectrum(:,end))./5e-8),...
            pow2db(abs(real(Pspectrum(:,end).*conj(Uspectrum(:,end))))./1e-12),...
            pow2db(abs(imag(Pspectrum(:,end).*conj(Uspectrum(:,end))))./1e-12)],'LineWidth',1);
        set(semilogx1(1),'DisplayName','Sound pressure','Color',[1 0 0]);
        set(semilogx1(2),'DisplayName','Particle velocity','Color',[0 0 1]);
        set(semilogx1(3),'DisplayName','Active intensity','Color',[0 0 0]);
        set(semilogx1(4),'DisplayName','Reactive intensity','Color',[0.5 0.7 0.5]);
        xlabel('Frequency (Hz)');
        ylabel(['Level (dB) at maximum distance (' num2str(X(end)) ' m)']);
        xlim([20 20000])
        box(axes1,'on');
        set(axes1,'XGrid','on','XMinorTick','on','XScale','log','YGrid','on');
        legend1 = legend(axes1,'show');
        set(legend1,'Location','northwest');
        title(['Piston radius = ' num2str(a) ' m, volume velocity = ' num2str(v0) ' m^3/s, surface velocity = ' num2str(v0./(pi*a.^2)) ' m/s'])
        hold off
        
        subplot(3,1,3)
        axes2 = subplot(3,1,3,'Parent',figure1);
        semilogx(farray,phasediff(:,end),...
            'LineWidth',1,'DisplayName','P-U phase difference','Color',[0 0.7 0]);
        xlabel('Frequency (Hz)');
        xlim([20 20000])
        ylabel('P-U phase difference (deg)');
        ylim([-90 90])
        set(axes2,'XGrid','on','XMinorTick','on','XScale','log','YGrid','on',...
            'YTick',[-90 -45 0 45 90],'YTickLabel',{'-90','-45','0','45','90'});
        
        % plot of pressure, intensity and particle velocity at 1 m
        beta = 1 ./ (1 + a.^2).^0.5;
        gamma = karray./2 ...
            .*((ones([length(karray),1]).^2 + a.^2).^0.5...
            - ones([length(karray),1]));
        P1spectrum = rho.*c.*v0.*...
            exp(1i*(-karray))...
            .*(1-exp(-2i*gamma));
        U1spectrum = v0.*exp(1i*(-karray))...
            .*(1-repmat(beta,[length(karray),1])...
            .*exp(-2i*gamma));
        
        figure1 = figure('Name','Rigid circular piston in infinite baffle 1 m on-axis');
        axes1 = subplot(3,1,1:2,'Parent',figure1);
        %axes1 = axes('Parent',figure1);
        hold(axes1,'on');
        semilogx1 = semilogx(farray,...
            [mag2db(abs(P1spectrum)./2e-5),...
            mag2db(abs(U1spectrum)./5e-8),...
            pow2db(abs(real(P1spectrum.*conj(U1spectrum)))./1e-12),...
            pow2db(abs(imag(P1spectrum.*conj(U1spectrum)))./1e-12)],'LineWidth',1);
        set(semilogx1(1),'DisplayName','Sound pressure','Color',[1 0 0]);
        set(semilogx1(2),'DisplayName','Particle velocity','Color',[0 0 1]);
        set(semilogx1(3),'DisplayName','Active intensity','Color',[0 0 0]);
        set(semilogx1(4),'DisplayName','Reactive intensity','Color',[0.5 0.7 0.5]);
        xlabel('Frequency (Hz)');
        ylabel('Level (dB) at 1 m');
        xlim([20 20000])
        box(axes1,'on');
        set(axes1,'XGrid','on','XMinorTick','on','XScale','log','YGrid','on');
        legend1 = legend(axes1,'show');
        set(legend1,'Location','northwest');
        title(['Piston radius = ' num2str(a) ' m, volume velocity = ' num2str(v0) ' m^3/s, surface velocity = ' num2str(v0./(pi*a.^2)) ' m/s'])
        hold off
        
        axes2 = subplot(3,1,3,'Parent',figure1);
        phasediff = 180.*(mod(pi+angle(P1spectrum)-angle(U1spectrum),2*pi)-pi)./pi;
        semilogx(farray,phasediff,...
            'LineWidth',1,'DisplayName','P-U phase difference','Color',[0 0.7 0]);
        xlabel('Frequency (Hz)');
        xlim([20 20000])
        ylabel('P-U phase difference (deg)');
        ylim([-90 90])
        set(axes2,'XGrid','on','XMinorTick','on','XScale','log','YGrid','on',...
            'YTick',[-90 -45 0 45 90],'YTickLabel',{'-90','-45','0','45','90'});
        

        
        if makeIR == 1
            X = 0:1/IRdivisions_per_m:size_of_calc;
            farray = (0:1:24000)';
            omegaarray = 2.*pi.*farray;
            karray = omegaarray./c;
            %beta = X ./ (X.^2 + a.^2).^0.5;
            gamma = repmat(karray,[1,length(X)])./2 ...
                .*((repmat(X,[length(karray),1]).^2 + a.^2).^0.5...
                - repmat(X,[length(karray),1]));
            Pspectrum = rho.*c.*v0.*...
                exp(1i*(-repmat(karray,[1,length(X)])...
                .*repmat(X,[length(karray),1])))...
                .*(1-exp(-2i*gamma));
            OUT.audio = ifftshift(ifft([Pspectrum;flip(conj(Pspectrum(2:end,:)))]));
            OUT.chanID = makechanID(length(X),2,[X',zeros(length(X),2)]);
            OUT.fs = 48000;
            OUT.cal = zeros(1,size(OUT.audio,2));
            OUT.properties.units = 'Pa';
            OUT.properties.units_ref = 2e-5;
            OUT.properties.units_type = 1;
        elseif makeIR == 2
            X = 0:1/IRdivisions_per_m:size_of_calc;
            farray = (0:1:24000)';
            omegaarray = 2.*pi.*farray;
            karray = omegaarray./c;
            beta = X ./ (X.^2 + a.^2).^0.5;
            gamma = repmat(karray,[1,length(X)])./2 ...
                .*((repmat(X,[length(karray),1]).^2 + a.^2).^0.5...
                - repmat(X,[length(karray),1]));
            Uspectrum = v0.*exp(1i*(-repmat(karray,[1,length(X)])...
                .*repmat(X,[length(karray),1])))...
                .*(1-repmat(beta,[length(karray),1])...
                .*exp(-2i*gamma));
            OUT.audio = ifftshift(ifft([Uspectrum;flip(conj(Uspectrum(2:end,:)))]));
            OUT.chanID = makechanID(length(X),2,[X',zeros(length(X),2)]);
            OUT.fs = 48000;
            OUT.cal = zeros(1,size(OUT.audio,2));
            OUT.properties.units = 'm/s';
            OUT.properties.units_ref = 5e-8;
            OUT.properties.units_type = 1;
        end
        
        
        
    case 3 % ARC FULL SPECTRUM
        if nargin == 0
            param = inputdlg({'Speed of sound in medium (m/s)';...
                'Density of medium (kg/m^3)';...
                'Volume velocity of source (m^3/s)';...
                'Piston radius (m)';...
                'Arc distance from source centre (m)';...
                'Angular resolution of arc (deg) for chart';...
                'Generate IRs? No [0], Yes (pressure) [1]';
                'IR angular resolution (deg)';...
                'Far-field analytic approximation [0] or numeric calculation (slow!) [1]'},...
                'Circular Piston Radiation Input Parameters',...
                [1 80],...
                {'342', '1.225', '0.0001', '0.2', '1', '1','1','30','0'});
            if length(param) < 9, param = []; end
            if ~isempty(param)
                c = str2num(char(param(1)));
                rho = str2num(char(param(2)));
                v0 = str2num(char(param(3)));
                a = str2num(char(param(4)));
                [val1, dist] = deal(str2num(char(param(5))));
                [val2, chartanglestep] = deal(str2num(char(param(6))));
                [val3, makeIR] = deal(str2num(char(param(7))));
                [val4, IRanglestep] = deal(str2num(char(param(8))));
                [val5, numeric] = deal(str2num(char(param(9))));
            else
                OUT = [];
                return
            end
        end
        if numeric == 1
            % CREATE SOURCE ARRAY (upper half)
            subsourcespacing = 0.125*c/24000; % 1/8 wavelength @ Nyquist freq
            Ys=-a:subsourcespacing:a;
            Yslen = length(Ys);
            Zs = (0:subsourcespacing:a)';
            Zslen = length(Zs);
            Ys = repmat(Ys,[Zslen,1]);
            Zs = repmat(Zs,[1,Yslen]);
            semicircle = abs(Ys+1i*Zs)<= a;
            Ys = Ys(semicircle); % Y subsource coordinates
            Zs = Zs(semicircle); % Z subsource coordinates
            N = length(Ys);
            
            if makeIR
                angledeg = 0:IRanglestep:90;
                karray = 2.*pi.*(1:24000)'./c;
                OUT.audio = zeros(48001,length(angledeg));
                Pspectrum = zeros(length(karray),length(angledeg));
                X = dist.*cos(pi.*angledeg(:)'./180);
                Y = dist.*sin(pi.*angledeg(:)'./180);
                X = repmat(X,[length(karray),1]);
                Y = repmat(Y,[length(karray),1]);
                karray = repmat(karray(:),[1,length(angledeg)]);
                Qs = subsourcespacing.^2 .* v0 ./ (pi*a.^2);
                Z0term = 2.*((-1i.*rho.*c.*karray)./(4.*pi)).*Qs; % pre-calculate for speed
                Z1term = 2.*Z0term;
                for J = 1:N
                    r = (X.^2 + (Ys(J)-Y).^2 + Zs(J).^2).^0.5;
                    if Zs(J) == 0
                        Pspectrum = Pspectrum + Z0term.*(exp(1i.*karray.*r)./r);
                    else
                        Pspectrum = Pspectrum + Z1term.*(exp(1i.*karray.*r)./r);
                    end
                end
                Pspectrum(isnan(Pspectrum)) = 0;
                OUT.audio = ifftshift(ifft([zeros(1,length(angledeg));Pspectrum;flip(conj(Pspectrum))]));
                [Xval,Yval,Zval]=sph2cart(pi.*angledeg(:)./180, zeros(length(angledeg),1), dist.*ones(length(angledeg),1));
                OUT.chanID = makechanID(length(angledeg),2,[Xval,Yval,Zval]);
                OUT.fs = 48000;
                OUT.cal = zeros(1,size(OUT.audio,2));
                OUT.properties.units = 'Pa';
                OUT.properties.units_ref = 2e-5;
                OUT.properties.units_type = 1;
            end
            
            angledeg = 0:chartanglestep:90;
            farray = (0:100:24000)';
            karray = 2.*pi.*farray'./c;
            Pspectrum = zeros(length(karray),length(angledeg));
            X = dist.*cos(pi.*angledeg(:)'./180);
            Y = dist.*sin(pi.*angledeg(:)'./180);
            X = repmat(X,[length(karray),1]);
            Y = repmat(Y,[length(karray),1]);
            karray = repmat(karray(:),[1,length(angledeg)]);
            Qs = subsourcespacing.^2 .* v0 ./ (pi*a.^2);
            for J = 1:N
                r = (X.^2 + (Ys(J)-Y).^2 + Zs(J).^2).^0.5;
                if Zs(J) == 0
                    Pspectrum = Pspectrum + 2.*((-1i.*rho.*c.*karray)./(4.*pi)).*Qs.*(exp(1i.*karray.*r)./r);
                else
                    Pspectrum = Pspectrum + 4.*((-1i.*rho.*c.*karray)./(4.*pi)).*Qs.*(exp(1i.*karray.*r)./r);
                end  
            end
            Pspectrum(isnan(Pspectrum)) = 0;
            
        else
            if makeIR
                angledeg = 0:IRanglestep:90;
                omega = 2.*pi.*(1:24000)'; % angular frequency
                karray = omega./c;
                OUT.audio = zeros(48001,length(angledeg));
                Pspectrum = zeros(length(karray),length(angledeg));
                X = dist.*cos(pi.*angledeg(:)'./180);
                Y = dist.*sin(pi.*angledeg(:)'./180);
                X = repmat(X,[length(karray),1]);
                Y = repmat(Y,[length(karray),1]);
                karray = repmat(karray(:),[1,length(angledeg)]);
                omega = repmat(omega(:),[1,length(angledeg)]);
                r = (X.^2 + Y.^2).^0.5; % distance from source centre
                theta = atan2(Y,X); % angle from source centre
                Pspectrum(theta==0)=-1i.*omega(theta==0).*rho.*c.*v0.*(a.^2).*0.5.*(exp(1i.*karray(theta==0).*r(theta==0))./r(theta==0));
                Pspectrum(theta~=0)=-1i.*omega(theta~=0).*rho.*c.*v0.*(a.^2).*(besselj(1,(karray(theta~=0).*a.*sin(theta(theta~=0))))./(karray(theta~=0).*a.*sin(theta(theta~=0)))).*(exp(1i.*karray(theta~=0).*r(theta~=0))./r(theta~=0));
                Pspectrum(isnan(Pspectrum)) = 0;
                OUT.audio = ifftshift(ifft([zeros(1,length(angledeg));Pspectrum;flip(conj(Pspectrum))]));
                [Xval,Yval,Zval]=sph2cart(pi.*angledeg(:)./180, zeros(length(angledeg),1), dist.*ones(length(angledeg),1));
                OUT.chanID = makechanID(length(angledeg),2,[Xval,Yval,Zval]);
                OUT.fs = 48000;
                OUT.cal = zeros(1,size(OUT.audio,2));
                OUT.properties.units = 'Pa';
                OUT.properties.units_ref = 2e-5;
                OUT.properties.units_type = 1;
            end
            angledeg = 0:chartanglestep:90;
            farray = (0:10:24000)';
            omega = 2.*pi.*farray; % angular frequency
            karray = omega./c;
            Pspectrum = zeros(length(karray),length(angledeg));
            X = dist.*cos(pi.*angledeg(:)'./180);
            Y = dist.*sin(pi.*angledeg(:)'./180);
            X = repmat(X,[length(karray),1]);
            Y = repmat(Y,[length(karray),1]);
            r = (X.^2 + Y.^2).^0.5; % distance from source centre
            theta = atan2(Y,X); % angle from source centre
            karray = repmat(karray(:),[1,length(angledeg)]);
            omega = repmat(omega(:),[1,length(angledeg)]);
            Pspectrum(theta==0)=-1i.*omega(theta==0).*rho.*c.*v0.*(a.^2).*0.5.*(exp(1i.*karray(theta==0).*r(theta==0))./r(theta==0));
            Pspectrum(theta~=0)=-1i.*omega(theta~=0).*rho.*c.*v0.*(a.^2).*(besselj(1,(karray(theta~=0).*a.*sin(theta(theta~=0))))./(karray(theta~=0).*a.*sin(theta(theta~=0)))).*(exp(1i.*karray(theta~=0).*r(theta~=0))./r(theta~=0));
            Pspectrum(isnan(Pspectrum)) = 0;
        end
        
        
        %         xsize = 1024;
        %         x = (0:xsize)'./xsize;
        %         Sat = abs(sin(2*pi*x));
        %         Hue = circshift(x,-round(xsize*30/360));
        %         Val = abs(x-0.5).*2;
        %         circularcolormap = hsv2rgb([Hue Sat Val]);
        
        HSVmap = hsv(1024);
        fig0 = figure('Name','Pressure as a function of angle','Colormap',HSVmap);
        ax0 = axes('Parent',fig0,'CLim', [-180 180],'XTick',[0 15 30 45 60 75 90]);
        hold(ax0,'on');
        h=surf(angledeg,farray,mag2db(abs(Pspectrum)./2e-5),'Parent',ax0);
        h.CData = 180.*angle(Pspectrum)./pi;
        h.FaceColor = 'texturemap'; % 'interp' does not work for circular colormaps
        h.EdgeLighting = 'gouraud';
        h.FaceAlpha = 1;
        h.FaceLighting = 'gouraud';
        h.AmbientStrength = 0.6;
        h.DiffuseStrength = 1;
        h.EdgeAlpha = 0.1;
        zlim([max(max(mag2db(abs(Pspectrum)./2e-5)))-60, max(max(mag2db(abs(Pspectrum)./2e-5)))]);
        xlabel('Angle (deg)')
        ylabel('Frequency (Hz)')
        zlabel('Sound Pressure Level (dB)')
        title(['Source radius = ' num2str(a) ' m (colormap shows phase)'])
        xlim([0 90])
        camlight left
        view(ax0,[71 44]);
        grid(ax0,'on');
        colorbar('peer',ax0,'Ticks',[-180 -90 0 90 180],...
            'TickLabels',{'-180','-90','0 deg','90','180'});
        hold off
        
        % 3D polar plot surface (freq is z axis)
        
        % Q and DI
        % half-space circular harmonic analysis???
    case 4
        if nargin == 0
            param = inputdlg({'Speed of sound in medium (m/s)';...
                'Density of medium (kg/m^3)';...
                'Volume velocity of source (m^3/s)';...
                'Piston radius (m)';...
                'Pressure only [1]; Pressure and radial particle velocity [2];Pressure and XY particle velocity [3]';...
                'X coordinate (axial distance) (m)';...
                'Y coordinate (lateral distance) (m)';...
                'Z coordinate (vertical distance) (m)'},...
                'Circular Piston Radiation Input Parameters',...
                [1 80],...
                {'342', '1.225', '0.0001', '0.2', '2'});
            if length(param) < 5, param = []; end
            if ~isempty(param)
                c = str2num(char(param(1)));
                rho = str2num(char(param(2)));
                v0 = str2num(char(param(3)));
                a = str2num(char(param(4)));
                [val1] = str2num(char(param(5)));
                [val2,X] = deal(str2num(char(param(6))));
                [val3,Y] = deal(str2num(char(param(7))));
                [val4,Z] = deal(str2num(char(param(8))));
                val5 = [];

            else
                OUT = [];
                return
            end
        end
        
        warndlg('SORRY - CODE FOR THIS OPTION IS NOT WRITTEN YET')
        OUT = [];
        return
    case 5
        if nargin == 0
            param = inputdlg({'Speed of sound in medium (m/s)';...
                'Density of medium (kg/m^3)';...
                'Volume velocity of source (m^3/s)';...
                'Piston radius (m)';},...
                'Circular Piston Radiation Input Parameters',...
                [1 80],...
                {'342', '1.225', '0.0001', '0.2'});
            if length(param) < 4, param = []; end
            if ~isempty(param)
                c = str2num(char(param(1)));
                rho = str2num(char(param(2)));
                v0 = str2num(char(param(3)));
                a = str2num(char(param(4)));
                [val1,val2,val3,val4,val5] = deal([]);
            else
                OUT = [];
                return
            end
        end
        f = 20:20000;
        k = 2*pi*f./c;
        Q = (k.*a).^2 ./ (1 -besselj(1,2.*k.*a)./(k.*a)); % directivity factor
        Force = rho.*c.*v0.*(1-besselj(1,2.*k.*a)./(k*a) + 1i*StruveH1(2.*k.*a)./(k.*a));
        Zrad = Force ./ (pi*a.^2 .* v0); % radiation impedance
        Power = 0.5 .* abs(v0).^2 .* rho.*c./(pi.*a.^2).*(1 -besselj(1,2.*k.*a)./(k.*a));
        
        figure1 = figure('Name','Rigid circular piston in infinite baffle');
        axes1 = axes('Parent',figure1);
        hold(axes1,'on');
        plot1 = plot(f,[real(Zrad)./(rho*c);imag(Zrad)./(rho*c)],'LineWidth',1);
        set(plot1(1),'DisplayName','Real(Z)','Color',[1 0 0]);
        set(plot1(2),'DisplayName','Imag(Z)','Color',[0 0 1]);
        xlabel('Frequency (Hz)');
        ylabel('Radiation impedance relative to rho*c');
        xlim([20 20000])
        box(axes1,'on');
        set(axes1,'XGrid','on','XMinorTick','on','XScale','log','YGrid','on','YScale','log');
        legend1 = legend(axes1,'show');
        set(legend1,'Location','southeast');
        title(['Piston radius = ' num2str(a) ' m'])
        
        figure1 = figure('Name','Rigid circular piston in infinite baffle');
        axes1 = axes('Parent',figure1);
        hold(axes1,'on');
        semilogx1 = semilogx(f,pow2db(Q),'LineWidth',1);
        set(semilogx1,'DisplayName','Directivity Index','Color',[0 0.7 0]);
        xlabel('Frequency (Hz)');
        ylabel('Directivity Index (dB)');
        xlim([20 20000])
        ylim([0 max(pow2db(Q))])
        set(axes1,'YColor',[0 0 0]);
        yyaxis right
        set(axes1,'YColor',[0 0 0]);
        ylim([1 max(Q)])
        axes1.YScale = 'log';
        ylabel('Directivity Factor')
        box(axes1,'on');
        set(axes1,'XColor',[0 0 0],'XGrid','on','XMinorTick','on','XScale','log',...
        'YGrid','on','ZColor',[0 0 0]);
        title(['Piston radius = ' num2str(a) ' m'])
        
        figure1 = figure('Name','Rigid circular piston in infinite baffle');
        axes1 = axes('Parent',figure1);
        hold(axes1,'on');
        semilogx1 = semilogx(f,pow2db(Power./1e-12),'LineWidth',1);
        set(semilogx1,'DisplayName','Sound Power Level','Color',[0.5 0.5 0]);
        xlabel('Frequency (Hz)');
        ylabel('Sound power level (dB)');
        xlim([20 20000])
        box(axes1,'on');
        set(axes1,'XGrid','on','XMinorTick','on','XScale','log','YGrid','on');
        title(['Piston radius = ' num2str(a) ' m, volume velocity = ' num2str(v0) ' m^3/s, surface velocity = ' num2str(v0./(pi*a.^2)) ' m/s'])
    otherwise
        OUT = [];
        return
end
% Function callbacks
OUT.funcallback.name = 'circular_rigid_piston_infbaffle.m';
OUT.funcallback.inarg = {a,c,rho,v0,method,val1,val2,val3,val4,val5};
end %eof

%**************************************************************************
% Copyright (c) 2016, Densil Cabrera
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%  * Redistributions of source code must retain the above copyright notice,
%    this list of conditions and the following disclaimer.
%  * Redistributions in binary form must reproduce the above copyright
%    notice, this list of conditions and the following disclaimer in the
%    documentation and/or other materials provided with the distribution.
%  * Neither the name of the University of Sydney nor the names of its contributors
%    may be used to endorse or promote products derived from this software
%    without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
% TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
% PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER
% OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%**************************************************************************