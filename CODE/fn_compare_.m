function fn_compare_(path2frames,filename,sta,x_opt,xm,X_fin,y_opt,ym,Y_fin)

[~,remain] = strtok(filename,'_');

% load results of alignment/filtering
eval(['load ',filename,'_results'])
% load frames used
eval(['load ',path2frames,remain,'.mat'])

nframes = numel(frames_used);

[N,M] = size(X_fin);

circsty = [{'og'},{'or'},{'ob'},{'om'},{'ok'}];

start = nframes - (sta + N - 1) ;

for i = start+sta:sta+N-1
    imshow(frames_used(1,i).cdata)
    hold on
    idx = i-sta-start+1;
    plot(x_opt(idx),y_opt(idx),'xc','MarkerSize',10)
    plot(xm(idx),ym(idx),'xy','MarkerSize',10)
    for m = 1:M
        plot(X_fin(idx,m),Y_fin(idx,m),circsty{m})
    end
    legd = legend(['OPT','MEAN',X_names]);
    set(legd,'FontAngle','italic','FontSize',8,'Location','SouthEast');
    pause(0.5)
    hold off
end


end

