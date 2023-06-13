function [x,xm,x_opt,y,ym,y_opt,sta] = fn_proc_filen_(filename,prefs)

% load raw data
dta = fn_dta_load_(filename);

if isempty(prefs.plot_raw)
    plot_raw = input('\nTo plot raw coordinates press 1,\nor 0 if not to, and then ENTER!\n');
else
    plot_raw = prefs.plot_raw;
end
if plot_raw
    
    % plot raw data
    fn_dta_plot_raw_(dta);
    
end

% truncate raw data
[x,y,x_opt,y_opt,xm,ym,sta] = fn_dta_trun_(dta);

if isempty(prefs.plot_trun)
    plot_trun = input('\nTo plot truncated coordinates press 1,\nor 0 if not to, and then ENTER!\n');
else
    plot_trun = prefs.plot_trun;
end
if plot_trun

    % plot truncated data
    fn_dta_plot_trun_(x,y,x_opt,y_opt,filename);
    
end

end
%% fn_dta_load_
function dta = fn_dta_load_(filename)

eval(['load ',filename]);
dta = fn_dta_(output);

end
%% fn_dta
function output = fn_dta_(input)

    N = length(input.workingFrames);
    xy = zeros(N,2);
    cc = xy;
    cascade = xy;
    points = xy;
    foreground = xy;
    output = struct('xy', xy, 'cc', cc, 'cascade', cascade, 'points', ...
        points, 'foreground',foreground);
    for n = 1:N
        % xy
        if ~isempty(input.xy{n})
            output.xy(n,1) = input.xy{n}(1);
            output.xy(n,2) = input.xy{n}(2);
        else
            if n == 1
                output.xy(n,1) = 0;
                output.xy(n,2) = 0;
            else
                output.xy(n,1) = output.xy(n-1,1);
                output.xy(n,2) = output.xy(n-1,2);
            end
        end
        % cc
        if ~isempty(input.coords_cc{n})
            output.cc(n,1) = input.coords_cc{n}(1);
            output.cc(n,2) = input.coords_cc{n}(2);
        else
            if n == 1
                output.cc(n,1) = 0;
                output.cc(n,2) = 0;
            else
                output.cc(n,1) = output.cc(n-1,1);
                output.cc(n,2) = output.cc(n-1,2);                
            end
        end
        % cascade
        if ~isempty(input.coords_cascade{n})
            output.cascade(n,1) = input.coords_cascade{n}(1);
            output.cascade(n,2) = input.coords_cascade{n}(2);
        else
            if n == 1
                output.cascade(n,1) = 0;
                output.cascade(n,2) = 0;
            else
                output.cascade(n,1) = output.cascade(n-1,1);
                output.cascade(n,2) = output.cascade(n-1,2);
            end
        end
        % points
        if ~isempty(input.coords_points{n})
            output.points(n,1) = input.coords_points{n}(1);
            output.points(n,2) = input.coords_points{n}(2);
        else
            if n == 1
                output.points(n,1) = 0;
                output.points(n,2) = 0;
            else
                output.points(n,1) = output.points(n-1,1);
                output.points(n,2) = output.points(n-1,2);
            end
        end
        % foreground
        if ~isempty(input.coords_foreground{n})
            output.foreground(n,1) = input.coords_foreground{n}(1);
            output.foreground(n,2) = input.coords_foreground{n}(2);
        else
            if n == 1
                output.foreground(n,1) = 0;
                output.foreground(n,2) = 0;
            else
                output.foreground(n,1) = output.foreground(n-1,1);
                output.foreground(n,2) = output.foreground(n-1,2);
            end
        end
    end
    
end
%% fn_dta_plot_raw_
function fn_dta_plot_raw_(dta)

if isempty(dta)
    
    error('No data!');
    
else
    
    figure
    plot(dta.xy(:,1),'c','LineWidth',1.5)
    hold on
    plot(dta.cc(:,1),'r')
    plot(dta.cascade(:,1),'m')
    plot(dta.foreground(:,1),'g')
    title('\bfRaw X-coordinates')
    legd = legend('Optimal X&Y','X-correll','Cascade','Foreground');
    set(legd,'FontAngle','italic','FontSize',8,'Location','SouthEast');
    axis auto
    hold off
    
    figure
    plot(dta.xy(:,2),'c','LineWidth',1.5)
    hold on
    plot(dta.cc(:,2),'r')
    plot(dta.cascade(:,2),'m')
    plot(dta.foreground(:,2),'g')
    title('\bfRaw Y-coordinates')
    legd = legend('Optimal X&Y','X-correll','Cascade','Foreground');
    set(legd,'FontAngle','italic','FontSize',8,'Location','SouthEast');
    axis auto
    hold off
    
end

end
%% fn_dta_trun_
function [x,y,x_opt,y_opt,xm,ym,sta] = fn_dta_trun_(dta)

% collect all the coordinates
dta_prov = [dta.xy(:,1),dta.cc(:,1),dta.cascade(:,1),dta.foreground(:,1),...
    dta.xy(:,2),dta.cc(:,2),dta.cascade(:,2), dta.foreground(:,2)];

% search for zeros in all the coordinates
idx0 = find(any(dta_prov==0),1,'first');

% first non-zero point used for truncation
sta = find(dta_prov(:,idx0)~=0,1,'first');

% optimal coordinates
x_opt = dta.xy(sta:end,1);
y_opt = dta.xy(sta:end,2);

% all other coordinates
x = [dta.cc(sta:end,1), dta.cascade(sta:end,1), dta.foreground(sta:end,1)]; 
y = [dta.cc(sta:end,2), dta.cascade(sta:end,2), dta.foreground(sta:end,2)]; 

% all mean coordinates
xm = mean(x,2);
ym = mean(y,2);

end
%% fn_dta_plot_trun_
function fn_dta_plot_trun_(x,y,x_opt,y_opt,filename)

if isempty(x) || isempty(y) || isempty(x_opt) || isempty(y_opt)

    error('No data!')
    
else
    
    % general figure settings -------------------------------
    % -------------------------------------------------------
    set(0, 'defaultFigurePaperType', 'A4')
    set(0, 'defaultFigurePaperUnits', 'centimeters')
    set(0, 'defaultFigurePaperPositionMode', 'auto')
    set(0, 'defaultFigurePaperOrientation', 'landscape')
    scrsz = get(0,'ScreenSize');
    % -------------------------------------------------------
    
    % -------------------------------------------------------
    hf1=figure;
    % -------------------------------------------------------
    plot(x_opt,'c','LineWidth',1.5)
    hold on
    plot(x)
    title('\bfTruncated X-coordinates')
    legd = legend('Optimal X&Y','X-correll','Cascade','Foreground');
    set(legd,'FontAngle','italic','FontSize',8,'Location','SouthEast');
    axis auto
    hold off
    % -------------------------------------------------------
    set(gcf,'Visible','on');
    set(gcf,'Position',[1 scrsz(2) 0.75*scrsz(3) 0.75*scrsz(4)]);
    print(hf1, '-dpdf',[filename,'_truncated_COORD_X.pdf']);
    % -------------------------------------------------------
    close(hf1);
    % -------------------------------------------------------
    
    % -------------------------------------------------------
    hf2=figure;
    % -------------------------------------------------------
    plot(y_opt,'c','LineWidth',1.5)
    hold on
    plot(y)
    title('\bfTruncated Y-coordinates')
    legd = legend('Optimal X&Y','X-correll','Cascade','Foreground');
    set(legd,'FontAngle','italic','FontSize',8,'Location','SouthEast');
    axis auto
    hold off
    % -------------------------------------------------------
    set(gcf,'Visible','on');
    set(gcf,'Position',[1 scrsz(2) 0.75*scrsz(3) 0.75*scrsz(4)]);
    print(hf2, '-dpdf',[filename,'_truncated_COORD_Y.pdf']);
    % -------------------------------------------------------
    close(hf2);
    % -------------------------------------------------------
    
end

end