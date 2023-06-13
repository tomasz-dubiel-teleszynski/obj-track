function [X_fin,X_diff] = fn_fin_cord_(filename,labl,x,xm,x_opt,finests,names,varargin)

I = numel(varargin);
if numel(names) ~= I
    error('Different number of names than state estimates!')
end
linesty = [{'g'},{'r'},{'b'},{'--m'},{'--k'}];

% general figure settings -------------------------------
% -------------------------------------------------------
set(0, 'defaultFigurePaperType', 'A4')
set(0, 'defaultFigurePaperUnits', 'centimeters')
set(0, 'defaultFigurePaperPositionMode', 'auto')
set(0, 'defaultFigurePaperOrientation', 'landscape')
scrsz = get(0,'ScreenSize');
% -------------------------------------------------------

% -------------------------------------------------------
hf=figure;
% -------------------------------------------------------
plot(x_opt,'c','LineWidth',1.5)
hold on
plot(xm,'y','LineWidth',1.5)
X_fin = zeros(size(x,1),I);
for i = 1:I
    x_fin = fn_fin_cord_calc_(x,finests,varargin{i});
    plot(x_fin,linesty{i},'LineWidth',1)
    X_fin(:,i) = x_fin;
end
title(['\bf',filename])
xlabel(['\bf Optimal, mean and final ',labl,' coordinate estimates'])
legd = legend(['OPT','MEAN',names]);
set(legd,'FontAngle','italic','FontSize',8,'Location','SouthEast');
axis auto
hold off
% -------------------------------------------------------
set(gcf,'Visible','on');
set(gcf,'Position',[1 scrsz(2) 0.75*scrsz(3) 0.75*scrsz(4)]);
print(hf, '-dpdf',[filename,'_final_COORD_',labl,'.pdf']); 
% -------------------------------------------------------
close(hf);
% -------------------------------------------------------

% -------------------------------------------------------
hf=figure;
% -------------------------------------------------------
plot(x_opt-xm,'y','LineWidth',1.5)
hold on
X_diff = zeros(size(x,1),I);
for i = 1:I
    X_diff(:,i) = x_opt-X_fin(:,i);
    plot(X_diff(:,i),linesty{i},'LineWidth',1)
end
title(['\bf',filename])
xlabel(['\bf Differences between optimal and mean and final ',labl,' coordinate estimates'])
legd = legend(['MEAN',names]);
set(legd,'FontAngle','italic','FontSize',8,'Location','SouthEast');
axis auto
hold off
% -------------------------------------------------------
set(gcf,'Visible','on');
set(gcf,'Position',[1 scrsz(2) 0.75*scrsz(3) 0.75*scrsz(4)]);
print(hf, '-dpdf',[filename,'_final_DIFF_',labl,'.pdf']); 
% -------------------------------------------------------
close(hf);
% -------------------------------------------------------

end
%% fn_fin_cord_calc_
function fin_cord = fn_fin_cord_calc_(x,finests,state)

idx = fn_idx_(x);
wgts = fn_wgts_(finests);

fin_cord = sum(kron(ones(size(x,1),1),wgts').*state(:,logical(idx)),2);

end
%% fn_wgts_
function wgts = fn_wgts_(finests)

wgts = 1./sqrt(exp(finests(2:end))); % inverse standard deviation of noise
wgts = wgts/sum(wgts);

end
%% fn_idx_
function idx = fn_idx_(x)

idx = zeros(2*size(x,2)+1,1);
for n = 1:numel(idx)-1
    if mod(n,2) ~= 0
        idx(n)=1;
    end
end

end