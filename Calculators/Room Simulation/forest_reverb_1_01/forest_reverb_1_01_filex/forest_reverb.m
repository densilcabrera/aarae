function audio = forest_reverb(N_trees,N_sc_max,forestwidth,minradius,maxradius,fs)

%This script calculates the impulse response for a forest of 
%hard cylinders. There are a number of parameters at the top of the file
%that you can adjust.

%This is currently copyright Travis Wiens 2008 but I intend to GPL it. 
%Please contact me if you want me to hurry this process up or you want
%a commerical license. email: travis.mlfx@nutaksas.com

% Slightly modified by Densil Cabrera 2013 for the AARAE project
% Basically the script was converted to a function, returning audio to the
% calling function. Also random number seeding was updated due to change in
% Matlab language. Modification and inclusion in AARAE based on the BSD
% license (from Matlab Central).

save_data=false;%should we save data?
save_wav=false;%should we write a wav file?
plot_data=true;%should we plot data?
zero_delay=true;%if true, chops pre-echo so that the direct path impulse arrives 
%   at sample 1, otherwise it will arrive at sample k_prewindow

if nargin < 1
    N_trees=10;%number of trees in forest
end

if nargin < 2
    N_sc_max=5;%maximum number of scatterings in each path (calculation time increases a lot)
end

if nargin <3
    X_limits=100*[-1 -1;1 1];%(m) limits of forest [x y]
else
    X_limits = [-forestwidth -forestwidth; forestwidth forestwidth];
end

if nargin <4
    a_tree_limits=[0.1 0.2];%(m) limtis to tree radius
else
    a_tree_limits=[minradius maxradius];
end

x_sp=1/2;%exponential on sound propagation equation 
% (1 for spherical spreading, 1/2 for cylindrical, 0 for plane wave)
d0=1;%(m) distance for gain G0;
G0=1;%gain at d0

if nargin <6
    F_s=44100;%(Hz) Sampling frequency
else
    F_s = fs;
end
N_fft=512;%number of points in fft and individual impulse responses
N_m=0;%number of terms in scattering equation series (or use 0 to set automatically)
f_f=17e3;%(Hz) filter cutoff frequency
minphase=false;%change phase to minimum phase filter form?
c=340;%(m/s) speed of sound
N_bits=24;%number of bits to save wav (some older software may require 16)

k_prewindow=100;%number of samples to take from end of impulse (pre echo)
N_window=300;%total number of samples to use from impulse

%end of parameters

%seed = rand('seed');%get seed from PRNG (in case we want to re-run something)
%rand('seed',seed);%set PRNG to an old method (since we don't need security)
%rng('default');

X_tree=[X_limits(1,1)+diff(X_limits(:,1))*rand(1,N_trees); X_limits(1,2)+diff(X_limits(:,2))*rand(1,N_trees)];
%(m) tree position vectors ([x;y])
X_t=[X_limits(1,1)+diff(X_limits(:,1))*rand(1,1); X_limits(1,2)+diff(X_limits(:,2))*rand(1,1)];
%(m) transmitter position [x;y]
X_r=[X_limits(1,1)+diff(X_limits(:,1))*rand(1,1); X_limits(1,2)+diff(X_limits(:,2))*rand(1,1)];
%(m) receiver position [x;y]

a_tree=a_tree_limits(1)*ones(1,N_trees)+diff(a_tree_limits)*rand(1,N_trees);%(m) tree radii

if plot_data
    figure(1);%plot forest layout
    h1=plot(X_t(1,:),X_t(2,:),'kx');
    hold on
    h2=plot(X_r(1,:),X_r(2,:),'k+');
    h3=plot(X_tree(1,:),X_tree(2,:),'go');
    legend([h1 h2 h3],'transmitter','reciever','tree','Location','Best')
    axis([X_limits(1,1) X_limits(2,1) X_limits(1,2) X_limits(2,2)])
    hold off
    xlabel('X_1 (m)')
    ylabel('X_2 (m)')
end



tic;

X_all=[X_t X_tree X_r];%aggregate all the location vectors

Y_sc=cell(N_trees+1,N_trees,N_trees+1);%cells to store scattering transfer functions

count=0;%counter
N_sc_calc=N_trees*((N_trees+1)+ (N_trees+3)*(N_trees)/2);%total number of scatterings to calculate
dt_i=.45;%(s) estimated calculation time per scattering (this is for a Core2Duo 6400 @ 2.13 GHz with 2 Gb RAM)
fprintf('Total impulses = %d, estimated t1 = %0.1f s = %0.1f m = %0.1f h\n',N_sc_calc,N_sc_calc*dt_i,N_sc_calc*dt_i/60,N_sc_calc*dt_i/3600)

N_sc_max_paths=sum(N_trees.^(1:N_sc_max).*(1:N_sc_max));
%maximum possible number of scatterings in the paths. This includes self
%reflection so it will be overestimated. If anyone knows an analytical
%expression for the proper number, I'd appreciate it.

dt_s=1e-4;%(s) estimated time to add an impulse
fprintf('Max scatterings = %d, t2 = %0.1f s = %0.1f m = %0.1f h (overestimate)\n',N_sc_max_paths,N_sc_max_paths*dt_s,N_sc_max_paths*dt_s/60,N_sc_max_paths*dt_s/3600)

for i=1:(N_trees+1);%select from transmitter plus trees
    X1=X_all(:,i);%incident ray origin
    t=toc;
    fprintf('Origin %d/%d. %d/%d impulses done. %0.1f s elapsed %0.1f s remaining\n',i,(N_trees+1),count,N_sc_calc, t, t*(N_sc_calc/count-1))
    for j=(1:(N_trees))+1;%select from trees only
        X2=X_all(:,j);%scattering point
        if i==1
            k_idx=2:(N_trees+2);
            %select from trees that haven't been starting points plus
            %receiver. Note that the path ABC has the same transfer
            %function as path CBA, so we can save a bunch of calculations.
        else
            k_idx=i:(N_trees+2);
        end
        for k=k_idx;
          X3=X_all(:,k); %scattered ray endpoint
          theta1=atan2(X2(2,:)-X1(2,:),X2(1,:)-X1(1,:));%angle of incident ray
          theta2=atan2(X3(2,:)-X2(2,:),X3(1,:)-X2(1,:));%angle of scattered ray
          phi=theta2-theta1;%scattering angle (0 for no scattering, pi for full reflection)
          [y_tmp, Y_sc{i,j,k}]=scatter_impulse(a_tree(j-1),phi,F_s,N_fft,N_m,f_f, minphase);
          %calculate complex transfer function of scattering Y_sc
          count=count+1;%increment counter
        end
    end
end
fprintf('Total impulses generated = %d\n',count)
t1=toc;
fprintf('Impulse generation time=%f s (%f s/scatter)\n',t, t/count)

N_path_max=sum(N_trees.^(1:N_sc_max));%maximum possible number of paths (includes self reflection)




d1=sqrt(sum((X_t-X_r).^2));%(m) distance for direct path
k_d1=round(d1/c*F_s)+k_prewindow;%(samples) delay for direct path. 
%Note zero is actually at sample number k_prewindow.  This is to prevent any
%of the prewindowed impulses from falling outside the array if there any
%trees or near the line better transmitter and receiver.

G1=G0*(d0./d1).^x_sp;%gain for direct path

N_imp=N_fft+k_prewindow;%initial size of impulse response
y_imp=zeros(1,N_imp);
y_imp(k_prewindow)=G1;%perfect impulse response of direct path
count=1;%reset counter
count_ref=0;%reset counter
for N_r=1:N_sc_max
    fprintf('%d scatterer(s): ',N_r)
    idx_tree=ones(1,N_r);%this matrix includes the index of the trees in this path
    while idx_tree(end)<=N_trees %stop when all paths have been cycled through
        if all(diff(idx_tree)~=0) 
            %only calculate if there are no consecutive identical trees, 
            %ie ABCA and ABAB are ok but AABA and ABBB aren't
            X_path=[X_t X_tree(:,idx_tree) X_r];%locations of points on sound path
            Y_path=ones(1,N_fft);%path transfer function
            d_path=0;%(m) total distance along path
            for j=1:N_r
                d=sqrt(sum((X_path(:,j)-X_path(:,j+1)).^2));%(m) length of incident ray
                d_path=d_path+d;%add to cumulative path distance
                if j==1%sound source
                    idx1=1;%index of incident ray origin
                else
                    idx1=idx_tree(j-1)+1;
                end
                idx2=idx_tree(j)+1;%index of scattering point
                if j==N_r%receiver
                    idx3=N_trees+2;%index of scattered ray endpoint
                else
                    idx3=idx_tree(j+1)+1;
                end
                if (idx3<idx1);%switch direction because path ABC=CBA
                    idx_tmp=idx1;
                    idx1=idx3;
                    idx3=idx_tmp;
                end
                Y_path=Y_path.*Y_sc{idx1,idx2,idx3};%apply scattered TF
                count_ref=count_ref+1;%increment counter
            end
            d=sqrt(sum((X_path(:,j+1)-X_path(:,j+2)).^2));%(m) lenght of last scattered ray
            d_path=d_path+d;%total path length
            k_d=round(d_path/c*F_s)-k_d1+k_prewindow;%(samples) delay for each path
            y_path=real(ifft(Y_path));%get path's impulse response (use real part due to rounding errors)
            G_path=G0*(d0./d_path).^x_sp;%gain for path due to spreading
            
            y_path=circshift(y_path',k_prewindow)';%shift impulse to use preecho
            if numel(y_imp)<(k_d+N_window)%check if path impulse will fit into y_imp
                y_imp=[y_imp zeros(1,k_d+N_window-numel(y_imp))];%expand size
            end
            y_imp(k_d+(1:N_window))=y_imp(k_d+(1:N_window))+G_path*y_path(1:N_window);%add path impulse
            
            count=count+1;
        
        end
        %increment idx_tree (like it's an N_tree-base number) for next path
        idx_tree(1)=idx_tree(1)+1;
        for j=1:(N_r-1)
            if idx_tree(j)>N_trees
                idx_tree(j)=1;
                idx_tree(j+1)=idx_tree(j+1)+1;
            end
        end
        
    end
    fprintf('%d cumulative scatterings, %d paths\n',count_ref,count-1)
end

t2=toc;
fprintf('Inpulse addition time=%f s (%e s/scattering)\n',t2-t1,(t2-t1)/count_ref)
fprintf('Total time=%f s\n',t2)

N_sc=count_ref;%total number of scatterings
N_paths=count-1;%total number of paths

if zero_delay
    y_imp(1:(k_prewindow-1))=[];%remove all pre-echo (possible if there is 
    %a tree near the line between transmitter and receiver) and shift
    %response
end

N_imp=numel(y_imp);%size of impulse response

if plot_data
    figure(2)
    plot((1:N_imp)/F_s,y_imp)
    xlabel('t (s)')
    ylabel('y_{imp}')

    figure(3)
    semilogy((1:N_imp)/F_s,abs(y_imp))
    axis([0 N_imp/F_s  1e-10 1])
    xlabel('t (s)')
    ylabel('y_{imp}')
end

dstring=datestr(now,30);%date/time string
if save_wav
    fname=['forest_impulse_' dstring '.wav'];
    wavwrite(y_imp'./max(abs(y_imp))*(1-2^(-N_bits+1)),F_s,N_bits,fname);
 end
if save_data
    fname=['forest_impulse_' dstring '.mat'];
    save(fname,'y_imp','k_prewindow','N_window','X_r','X_tree','X_t','N_m','N_fft','N_trees','N_sc_max','a_tree','N_sc','N_paths','Y_sc','seed','N_sc','N_paths','t1','t2','f_f','c','minphase','F_s','x_sp','G0','d0','d1','G1')
end

audio = y_imp' ./ max(abs(y_imp));