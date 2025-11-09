function [T, P] = groupplot(X, id, opt)
%[T, P] = groupplot(X, id, opt);
%
%INPUT:
%  X    Data to be plotted. Either 2D or 3D
%  id   Grouping. Can have any number of columns. Each unique combination
%        will be grouped and an ellipse will be drawn (not advisable for
%        more than 10 unique combinations)
%  opt  Option structure. Type 'opt = groupplot' for default values
%
%OUTPUT:
%  T    Rotated and centered score values for each group
%  P    Loading for each group
%
% See also: codeplot

% AAR 150714 Added some more possibilities in the options
% AAR 250713 Now giving the local scores and loadings as outputs
% AAR 180713 It's now possible to give the color code to be used as an
%             input variable
% AAR 010710 You can now have it calculate the confidence interval based on
%             bootstrapping. Seems quite handy
% AAR 100407 Works better with 2D
% AAR 200207 Sort the groups
% AAR 021106

if nargin == 0
    T.unc = [0 1];
    T.col = numcol( 10);
    T.alpha = [.7 1];
    T.lw = .5;
    T.new = true;
    T.view = true;
    T.Info = {
        'unc - Two numbers. First indicates, 0 - std, 1 - s.e. and 2 - confidence intervals';
        '       If given as a fraction this will be used as the correction factor to calculate';
        '       the standard error. The second number is the of errors added to each treatment';
        'col - A matrix with three columns indicating the colors which are to be used. Don''t';
        '       worry if it''s too large, it will only use the number of colors you need';
        'alpha - alpha for the face and edge of the ellipses (expert use)';
        'lw - Linewidth of the circles (expert use)';
        'new - True if you want a new figure, false if not';
        'view - True if you want to make the figure'};
    return
end

if nargin < 3
    opt = groupplot;
end

if isempty( opt.unc)
    opt.unc = [0 1];
end

if length( opt.unc) ~= 2
    opt.unc = [opt.unc(1) 1];
end

[rx,cx]=size(X);
%Check dimenstion
if cx<2 || cx>3
    X=X';
    [rx,cx]=size(X);
    if cx<2 || cx>3
        error('X should be 2D or 3D')
    end
end
if size(X,1)~=size(id,1)
    id=id';
    if size(X,1)~=size(id,1)
        error('X and id should have the same number of samples')
    end
end

[uid,nid]=fjernlike(id);
if size( opt.col, 1) ~= length( uid)
    if length( uid) > 10 
        if opt.view
            disp( 'Since there are more than 10 groups, and the correct colormap has been given')
            disp( ' everything will be grey-scaled')
        end
        opt.col = ones( length( uid), 1) * [.7 .7 .7];        
    end
end

%Sort the groups so that they are in increasing order
[uid,p]=sortrows(uid);
nid=nid(p,:);

% figure
if opt.view
    if opt.new
        cla
    end
    hold on
end
if cx==2
    an=0:pi/16:2*pi;
    x{1}=cos(an);
    x{2}=sin(an);
else
    [x{1},x{2},x{3}]=sphere(36);
end
for i=1:length(uid)
    tempid=nid(i,:);
    tempid(isnan(tempid))=[];
    if length(tempid)>=cx
        [xn,xm]=cen_std(X(tempid,:));
        [t,p]=pca(xn,cx);
        if opt.view
            for j=1:cx
                ax{j}=ones(size(x{j}))*xm(j);
                for k=1:cx
                    switch opt.unc(1)
                        case 0
                            ax{j}=ax{j}+x{k}*p(j,k)*std(t(:,k)) * opt.unc(2);
                        case 1
                            ax{j}=ax{j}+x{k}*p(j,k)*std(t(:,k))/sqrt( size( t, 1) ) * opt.unc(2);
                        case 2
                            ci = ciboot( t( :, k), 'mean', 5, opt.unc(2), 1e3);
                            ax{j} = ax{j} + x{k} * p( j, k) * max( abs( ci) );
                        otherwise
                            se = sqrt( sum( cen_std( t( :, k) ).^2) * opt.unc(1) );
                            ax{j} = ax{j} + x{k} * p( j, k) * se * opt.unc(2);
                    end
                end
            end
            if cx==2
                g(i)=patch(ax{1},ax{2},'-');
            else
                g(i)=surf(ax{1},ax{2},ax{3});
            end
        end
    else
        disp('One (or more) of the groups have fewer members than the wanted dimension')
        uid(i)=NaN;
    end
    T{i} = t;
    P{i} = p;
end

if opt.view
    if cx>2
        shading interp
        for i=1:length(uid)
            if ~isnan(uid(i))
                set(g(i),'FaceColor',opt.col(i,:))%,'FaceAlpha',.7)
            end
        end
        light('Position',max(X)+range(X)*.1)
        lighting phong
    else
        for i=1:length(uid)
            if ~isnan(uid(i))
                set(g(i),'FaceColor', opt.col(i,:),'EdgeColor', opt.col(i,:), ...
                    'FaceAlpha',opt.alpha(1), 'EdgeAlpha', opt.alpha(2), 'LineWidth', opt.lw)
            end
        end
    end
end