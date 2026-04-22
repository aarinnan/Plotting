function opt = codeplot( X, id, opt, leg, incol)
% opt = codeplot( X, id, opt, leg, incol);
%
% Makes a coded plot with 2 or 3 axes (size(X,2)) with as many different
% groups as identified in id. Can take up to four columns of id (different
% way of groups, e.g. {'Location','Type','Time'}
%
%INPUT
% X       Data to be plotted (max 3 columns)
% id      Identifier. Same number of rows as X (max 4 columns). (Please
%          note that the function sorts the IDs, using 'sortrows')
% opt     Options for 'codeplot'
% leg     If added, this will be used for the legends. Should be a cell matrix
%          with as many cells as columns in id. Each cell should hold the right
%          number of rows (equal number of group in corresponding column)
% incol  If added, this decides the coloring used during the mapping
% See also: GROUPPLOT

% AAR 280414 Added more colors in the first grouping possible (new input)
% AAR 130613 If there's only one coding the circles get the same color in
%             the outline as the interior
% AAR 271006

%280920 AAR There's "always" a discrepency between the colors in the plot
%            and the legends. This is because the color codes are sorted,
%            the legends are not. I think I want to remove the legend
%            option, and just always use the input values for the legends.
%            SO, one ALWAYS get the legend (unless one opt this out in an
%            'option' structure) 

if nargin == 0
    opt = struct( 'AxisNames', { ''}, 'Legends', 'on');
    opt.Colors = struct( 'FaceColor', numcol( 10), 'EdgeColor', flipud( numcol( 10) ));
    opt.MarkerSize = 10;
    opt.Info = { 'AxisNames  - The names of the axis to be used in the plot'
        '          given as a cell array'
        'Legends    - Either on or off'
        'Colors     - A structure with facecolors (1st coloring) and edgecolors'
        'MarkerSize - Size of the markers'};
    return    
end

if nargin < 3
    opt = codeplot;
end

[ xr, xc] = size( X);
shap = [ 'd', 's', 'o', '^', 'v', 'p', '<', '>', 'h'];
[ ir, ic] = size(id);
if ic > 4
    error( '''id'' can maximum have four columns')
end
if ir ~= xr
    error( '''X'' and ''id'' must have the same number of rows')
end

%Find the unique inputs of id
[ uid, nid] = fjernlike( id);
[ uid, j] = sortrows( uid);
nid=nid(j,:); %Why do I sort these?
for i = 1:size( id, 2)
    mid( i) = length( fjernlike( uid( :, i)));
end

%280920 AAR I should sort the ID columns in case some are with very many
%            different values. This, as I can only have 10 different
%            symbols

%Making the legends, now within the code itself
if nargin < 4 || isempty( leg)
    if strcmp( lower( opt.Legends), 'on')
        for cl = 1:size( uid, 2)
            leg{ cl} = cellstr( num2str( uid( :, cl) ) );
        end
    end
end

if mid(1) > size( opt.Colors.FaceColor, 1)
    error( ['This function cannot handle more than ' num2str( size( incol, 1) ) ' groups in the first column'])
end
if max( mid( 2:end) ) > 10
    error('This function cannot handle more than 10 groups. Sorry')
else
end

%Recode in case the groups are not 1:1:end
newuid=uid;
for i=1:size(uid,2)
    temp=sort(fjernlike(uid(:,i)));
    for j=1:length(temp)
        newuid(uid(:,i)==temp(j),i)=j;
    end
end
uid=newuid;
clear newuid

if length(mid)==1
    uid(:,2)=3;
elseif length(mid)==3
    uid(:,4)=9;
end

figure
subplot(4,4,[1:3 5:7 9:11 13:15])
hold on
for i=1:size(uid,1) %Number of unique groups
    tempid=nid(i,~isnan(nid(i,:)));
    if xc==2
        h = plot(X(tempid,1),X(tempid,2),['k' shap(uid(i,2))],'MarkerFaceColor', opt.Colors.FaceColor(uid(i,1),:), 'MarkerSize', opt.MarkerSize);
        if length(mid)>2
            set( h, 'MarkerEdgeColor', opt.Colors.EdgeColor( uid( i, 3), :), 'LineWidth', 2, 'MarkerSize', opt.MarkerSize);
        else
            set( h, 'MarkerEdgeColor', opt.Colors.FaceColor( uid( i, 1), :) )
        end
%             plot(X(tempid,1),X(tempid,2),shap(uid(i,4)),'MarkerEdgeColor',col(uid(i,3),:),'LineWidth',2,'MarkerSize',10)
%         end
    else
        plot3(X(tempid,1),X(tempid,2),X(tempid,3),['k' shap(uid(i,2))],'MarkerFaceColor',opt.Colors.FaceColor(uid(i,1),:), 'MarkerSize', opt.MarkerSize)
        if length(mid)>2
            plot3(X(tempid,1),X(tempid,2),X(tempid,3),shap(uid(i,4)),'MarkerEdgeColor',opt.Colors.EdgeColor(uid(i,3),:),'LineWidth',2,'MarkerSize',opt.MarkerSize)
        end
    end
end
if nargin > 3
    if isempty( opt.AxisNames)
        ax{1} = '';
        ax{2} = '';
        ax{3} = '';
    else
        ax=cellstr( opt.AxisNames);
    end
    if length(ax)==xc
        xlabel(ax{1},'FontWeight','Demi')
        ylabel(ax{2},'FontWeight','Demi')
        if xc == 3
            zlabel(ax{3},'FontWeight','Demi')
        end
    end
    if nargin > 3
        if ~isempty( leg)
            for i=1:length(leg)
                if ~iscell( leg{i})
                    templeg = cellstr( leg{i});
                else
                    templeg = leg{i};
                end
                subplot(4,4,i*4)
                hold on
                if i-floor(i/2)*2==1
                    for j=1:length( templeg)
                        if i == 1
                            tempcol = opt.Colors.FaceColor;
                        else
                            tempcol = opt.Colors.EdgeColor;
                        end
                        plot(0,length(templeg)+1-j,'ko','MarkerFaceColor',tempcol(j,:), 'MarkerSize', opt.MarkerSize)
                        text(.5,length(templeg)+1-j,templeg{j})
                    end
                else
                    for j=1:length(templeg)
                        plot(0,length( templeg)+1-j,['k' shap(j)],'MarkerFaceColor','k', 'MarkerSize', opt.MarkerSize)
                        text(.5,length( templeg)+1-j, templeg{j})
                    end
                end
                axis([0 2 .5 length(leg{i})+.5])
                axis off
            end
            subplot(4,4,4)
            if length(mid)<3
                title('Legend','FontWeight','Demi')
            else
                title('Legend - Inner symbol','FontWeight','Demi')
                subplot(4,4,12)
                title('Legend - Outer symbol','FontWeight','Demi')
            end
        end
    end %if nargin == 4
end
g = get( gcf, 'Children');
if length( g) == 2
    pos = get( g(1), 'Position');
    pos( [2 4]) = [.5 .4];
    set( g(1), 'Position', pos)
end