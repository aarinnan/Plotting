function op = icplot( X, opt)
%
%icplot( X, opt)
% Makes a plot to view the number of misclassifications in iECVA/ iPLS-DA
%
%INPUT:
% X    Dataset-struct or X-matrix with pre-processed data. Add all data in
%       order to get the misclassification in %
% opt  Struct-array. Type 'opt = icplot' for default

% GNU-license, 2010 - 
% This M-file and the code in it may be changed and modified as seen fit by
% the user. However, any redistribution of this code should follow the
% GNU-license, both if you give it away gratis or for a fee.
%
% Åsmund Rinnan
% Quality and Technology
% Department of Food Science
% Faculty of Sciences
% University of Copenhagen
% Rolighedsvej 30, 1958 Frederiksberg C, Denmark
% e-mail: aar@life.ku.dk

if nargin == 0
    op.Ax = NaN;
    op.axlim = [NaN NaN];
    op.mis = NaN;
    op.lim = NaN;
    op.ind = NaN;
    op.iset = 1;
    op.nLV = NaN;
    op.plot = [3 1];
    op.dir = 1;
    op.info = {'Ax = the axis for your spectra'; ...
        'axlim = minimum and maximum of y-axis'; ...
        'mis = Number of misclassifications in each interval + full (last number)'; ...
        'lim = Limits for the number of misclassifications used in the legend (5 numbers)'; ...
        'ind = Vector indicating intervals. Can be given in 3 ways, see below'; ...
        'iset = 1 - vector with same number for same interval, length = length of axis'; ...
        '       2 - vector specifying start and end of intervals (according to axis)'; ...
        '       3 - as above but with variable number instead'; ...
        'nLV = Number of LV''s used in each interval + full (last number)'; ...
        'plot = Two numbers. First number is plot version. 1-3 can all be given'; ...
        '       Second number indicates with or without number of LV''s indicated'; ...
        'dir = direction of x-axis. 1 - normal, other - reverse'};
    
    return
end

if isstruct( X)
    %Allow for input as given as dataset
    xspec = mean( X.data);
    num = size( X.data, 1)/ 100;
    if isnan( opt.Ax)
        ax = X.axisscale{ 2, 1};
    end
    clear X
else %Other input type
    if size( X, 1) > 1
        xspec = mean( X);
        num = size( X, 1)/ 100;
    else %If mean already given as input, there is no need to calculate the mean
        xspec = X;
        num = 1;
    end
    clear X
    ax = opt.Ax;
end

ax = ax(:)';

%X should be given in a way so that xspec can be calculated as the mean
%spectrum
if size( xspec, 1) ~= 1 || size( xspec, 2) == 1
    error( '''X'' not given correctly')
end

%Check if the 'ind' is given "normally" 
if length( opt.ind) ~= length( ax) && opt.iset == 1
    error( '''ind'' and ''ax'' should have the same length')
else %If it is given as the beginning and end of intervals, the total number should be divisible with 2
    if length( opt.ind)/2 ~= floor( length( opt.ind)/2)
        error( '''ind'' should be a set of paired numbers')
    end
end

if opt.iset < 1 || opt.iset > 3
    error( '''iset'' can only have the value 1, 2 or 3')
else
    %Transform the information into variable numbers
    if opt.iset == 2    
        for ci = 1:length( opt.ind)
            [~, j] = min( abs( ax - opt.ind( ci) ) );
            opt.ind( ci) = j;
        end
    end
    %Transform the information further into the "normal" format
    if opt.iset > 1
        temp = reshape( opt.ind, [2 length( opt.ind)/2] )';        
        temp = sortrows( temp, 1);
        for ci = 1:size( temp, 1)
            id = temp( :, 1) - temp( ci, 2) == 0;
            if sum( id) == 1
                temp( id, 1) = temp( temp( id, 1) - temp( :, 2) == 0, 2) + 1;
            end
        end
        %Make the "normal" format
        opt.ind = ones( length( ax), 1) * NaN;
        for ci = 1:size( temp, 1)
            opt.ind( temp( ci, 1):temp( ci, 2) ) = ci;
        end
    else
    opt.ind = opt.ind(:);
    end
end

if length( opt.nLV) ~= length( opt.mis)
    if ~isnan( opt.nLV)
        error( '''nLV'' and ''mis'' should have the same length')
    end
end

if (length( fjernlike( vec( opt.ind( ~isnan( opt.ind) ) ) ) ) + 1) ~= length( opt.mis)
    error( '''mis'' should be one longer than the number of intervals')
end

%In case these are data from iPLS-DA
if size( opt.nLV, 1) < size( opt.nLV, 2)
    opt.nLV = opt.nLV';
end
if size( opt.nLV, 2) > 1
    opt.nLV = [min( opt.nLV, [], 2)' max( opt.nLV, [], 2)']; %Range
    tempLV = mean( opt.nLV, 2);
    txtfor = '%2.0f-%2.0f';
else
    txtfor = '%2.0f';
    tempLV = opt.nLV;
end

%Set the ax limits if they haven't been set
if any( isnan( opt.axlim) ) || length( opt.axlim) ~= 2
    if any( opt.mis == 0)
        yax = - max( opt.mis)/ (num * 20);
    else
        yax = 0;
    end
    opt.axlim = [yax max( opt.mis)/ num];
    opt.axlim( 2) = opt.axlim( 2) + (opt.axlim( 2) - opt.axlim( 1))/ 20;
end
opt.axlim = sort( opt.axlim);

%Find the intervals and their centers
[i, j] = fjernlike( opt.ind);
j( isnan( i), :) = [];
tot = sum( ~isnan( j), 2);
tot = [j( :, 1) j( :, 1) + tot - 1];
id = round(mean( tot, 2))';
tot( tot( :, 1) > 1, 1) = tot( tot( :, 1) > 1, 1) - 1;
axid = ax(id);

fs = 10 - floor( size( tot, 1)/ 20);

if opt.plot( 1) == 1 || opt.plot( 1) == 4
    %Find the coloring range. Can have 5 number of colors
    if isnan( opt.lim)
        lim = [min( tempLV) median( tempLV)];
        if diff( lim) < 4
            lim( 2) = lim( 1) + 4;
        end
        opt.lim = floor( lim(1):diff( lim)/4:lim(2) );
        limfac = lim;
    else
        limfac = opt.lim;
    end
    axnLV = ones( length( tempLV) + 1, 1) * 5;
    for cm = 4:-1:1
        axnLV( tempLV <= opt.lim( cm) ) = cm;
    end
end

if opt.plot( 1) > 1
    %Find the coloring range. Can have 5 number of colors
    temp = opt.mis;
    if isnan( sum( opt.lim) ) || length( opt.lim) < 5 %If the limits are not given
        %0 is reserved to the green line/ white patch
        lim = [max( [min( temp), 1]) median( temp)];
        if diff( lim) < 3
            lim( 2) = lim( 1) + 3;
        end
        lim = [0 floor( lim(1):diff( lim)/3:lim(2) )];
    else
        lim = opt.lim(:)';
    end
    %Index the intervals according to color
    axmis = ones( length( temp), 1) * 5;
    for cm = 4:-1:1
        axmis( temp(:) <= lim( cm) ) = cm;
    end
end

%Need to find the intervals which have not been tested
id = find( tot( 2:end, 1) - tot( 1:end-1, 2) > 0);
if ~isempty( id)
    not = tot( sort( [id id+1], 2), :);
    not = [not( 1:2:end, 2) not( 2:2:end, 1)];
    %If the beginning of the spectra is missing
    if tot( 1, 1) > 1
        not = [1 tot( 1, 2); not];
    end
    %If the end of the spectra is missing
    if tot( end, 2) < length( ax)
        not = [not; tot( end, 2) length( ax)];
    end
else 
    not = [];
end

clf
switch opt.plot(1)
    case 1 %"Standard"-plot
        %Check if there are some classifications which are zero and set the
        %lower axis limit accordingly
        if sum( opt.mis == 0)
            yax = - max( opt.mis)/ (num * 20);
        else
            yax = 0;
        end
        
        xspec = xspec - min( xspec);
        
        %The spectra are either shown from 5-80% of the plot, or 10-85%
        ylim = ( max( opt.mis)/num - yax)/40;
        ylim = [ylim ylim + 20/3];
        numhist = hist( xspec, 50);
        if sum( numhist(1:3) )/ sum( numhist) > .25
            xspec = xspec * ( (max( opt.mis)/num) / range( xspec) * .75) + max( opt.mis) / (num * 10);
        else
            xspec = xspec * ( (max( opt.mis)/num) / range( xspec) * .75) + max( opt.mis) / (num * 20);
        end
        
        subplot( 2, 1, 1)
        %The colors
        if opt.plot(2) == 0
            col = flipud( [.3 .45 .6 .8 1]' * ones( 1, 3) );
        else %If you don't want the numbers, you get colors
            col = ones( 5, 3);
        end

        hold on
        if ~isempty( not)
            %Color all missing intervals
            for cb = 1:size( not, 1)
                g = patch( ax( not( cb, [1 2 2 1]) ), opt.axlim( [1 1 2 2]), [.95 .95 .95]);
                set( g, 'EdgeColor', [.95 .95 .95])
                plot( ax( not( cb, 1):not( cb, 2) ), xspec( not( cb, 1):not( cb, 2)), '-', 'Color', [.9 .9 .9]) %The line of the spectra
            end
        end
            
        for cb = 1:size( tot, 1)
            patch( ax( tot( cb, [1 2 2 1])), [yax yax opt.mis( cb)/num opt.mis( cb)/num], col( axnLV( cb), :))
            plot( ax( tot( cb, 1):tot( cb, 2) ), xspec( tot( cb, 1):tot( cb, 2) ), 'k-', 'LineWidth', 1)
        end
        %Make a line above and below the graph
        plot( ax( [1 end] ), ones( 2, 1) * opt.axlim, 'k-')
                                
        plot( ax( [1 end]), opt.mis( [end end])/num, 'r--', 'LineWidth', 1.5 )
        if ~isnan( opt.nLV( end, 1) )
            if length( fjernlike( opt.nLV( end, :)') ) == 1
                tempLV = opt.nLV( end, 1);
                txtf = '%2.0f';
            else
                tempLV = opt.nLV( end, :);
                txtf = txtfor;
            end
            if ax(1) > ax(end)
                text( ax( end) - (max(ax) - min( ax)) * .025, opt.mis( end)/num, num2str( tempLV, txtf) )
            else
                text( ax( end) + (max(ax) - min( ax)) * .025, opt.mis( end)/num, num2str( tempLV, txtf) )
            end
        end        
        
        axis( [min( ax) max( ax) opt.axlim ] )
        if num == 1
            ylabel( '# of misclassifications', 'FontSize', 12, 'FontWeight', 'Demi')
        else
            ylabel( '% of misclassifications', 'FontSize', 12, 'FontWeight', 'Demi')
        end
        
        set( gca, 'Box', 'on')
        
        if opt.plot(2) == 0
            subplot( 2, 1, 2)
            for cg = 1:5
                patch( [1 2 2 1], 6 - [cg cg cg + .7 cg + .7], col( cg, :) )
            end
            leg{ 5} = [num2str( limfac(end)) '+'];
            for cg = 1:4
                leg{cg} = num2str( limfac( cg) );
                if limfac( cg) + 1 < limfac( cg + 1)
                    leg{ cg} = [leg{ cg} '-' num2str( limfac( cg + 1) - 1)];
                end
            end
            text( ones( 5, 1) * 2.5, 6 - ( (1:5)' + .35), leg)
            axis( [.6 3.6 0 5.3] )
            title( '# of factors', 'FontSize', 12, 'FontWeight', 'Demi')
            set( gca, 'XTick', [], 'YTick', [], 'Box', 'on')
            g = get( gcf, 'Children');
            pos = [.82 .7 .15 .2; .1 .1 .7 .8];
            for cg = 1:2
                set( g(cg), 'Position', pos( cg, :) )
            end
        else
            g = get( gcf, 'Children');
            set( g, 'Position', [.1 .1 .8 .8])
        end
        
    case {2, 3, 4} %"Munck"-like plot
        
        %The legend should be a from-to legend
        temp = [0 lim( 1:3) + 1; lim( 1:4)];
        temp( 1, 1) = 0;
        temp( :, end + 1) = [lim( 4) + 1; lim(5)];
%         temp( 2, end) = 
%         id = [false temp( 1, 2:end) - temp( 2, 1:end-1) == 0];
%         temp( 1, id) = temp( 1, id) + 1;
%         temp( 1, end + 1) = lim( 4) + 1;
%         temp( 2, end) = lim( end);
        leg = cell( 5, 1);
        if num == 1
            numfor = '%2.0f';
            titname = '# of misclass';
        else
            numfor = '%2.1f';
            titname = '% of misclass';
        end
        leg{ 5} = [num2str( temp( 1, 5)/num, numfor ) '+'];
        for cl = 1:4
            if temp( 1, cl) == temp( 2,cl)
                leg{ cl} = num2str( temp( 1, cl)/num, numfor );
            else
                leg{ cl} = [num2str( temp( 1, cl)/num, numfor ) '-' num2str( temp( 2, cl)/num, numfor )];
            end
        end

        ylim = [min( xspec) max( xspec)];
        ylim = [ylim(1) - range( ylim) * .05 ylim(2) + range( ylim) * .2];

        %Check if the spectra should be shifted a bit up
        xb = hist( vec( xspec), 40);
        if sum( xb( 1:3) ) > length( xspec) * .25
            xspec = xspec + range( ylim) * .05;
        end

        if opt.plot(1) == 2
            col = flipud( (0:.25:1)' * ones( 1, 3) );%zeros( 5, 1) (0:.25:1)' zeros( 5, 1) ] );
            
            subplot( 4, 1, 1)
            for cb = 1:size( tot, 1)
                g = patch( ax( tot( cb, [1 2 2 1])), ylim( [1 1 2 2] ), col( axmis( cb), :));
            end
            
            hold on
            plot( ax, xspec, 'r', 'LineWidth', 2)
            
            axis( [min( ax) max( ax) ylim] )
            
            %Hardcode the legends as the legend command does not really
            %work :(
            subplot( 4, 1, 2)
            for cg = 1:5
                patch( [1 2 2 1], 6 - [cg cg cg + .7 cg + .7], col( cg, :) )
            end
            text( ones( 5, 1) * 2.5, 6 - ( (1:5)' + .35), leg)
            if num == 1
                axis( [.6 3.6 0 5.3] )
            else
                axis( [.6 5.1 0 5.3] )
            end
            
            set( gca, 'XTick', [], 'YTick', [], 'Box', 'on')
            title( titname, 'FontSize', 12, 'FontWeight', 'Demi')
            
        else
            
            %Change color code from green through yellow to red to
            %green to grey
%             col = [0 .6 0; .6 .9 .1; 1 .9 0; .9 0 0; .9 0 0];
            col = [0 .6 0; .6 .9 .1; .6 .6 .6; .8 .8 .8; .8 .8 .8];
            mr = {'-', '-', '-', '-', ':'};
%             col = [0 .8 0; .9 .9 0; .9 .6 0; .9 0 0; .75 .75 .75];
            lw = [3 3 2 2 1];
            
            if opt.plot(1) == 4
                cback = flipud( (.6:.1:1)' * ones( 1, 3) );
            end

            subplot( 4, 1, 1)
            hold on
            for cb = 1:size( tot, 1)                
                if opt.plot(1) == 4
                    g = patch( ax( tot( cb, [1 2 2 1])), ylim( [1 1 2 2]), cback( axnLV( cb), :));
                    set( g, 'EdgeColor', cback( axnLV( cb), :) )
                end
                if opt.plot(2) == 1 && opt.plot(1) == 3
                    plot( ones( 2, 1) * ax( tot( cb, 1:2)), ylim' * ones( 1, 2), 'k-')
                end
                g = plot( ax( tot( cb, 1):tot( cb, 2) ), xspec( tot( cb, 1):tot( cb, 2) ), mr{axmis( cb)});
                set( g, 'Color', col( axmis( cb), :), 'LineWidth', lw( axmis( cb) ) )
            end

            if ~isempty( not)
                %Color all missing intervals black
                for cb = 1:size( not, 1)
                    q = patch( ax( not( cb, [1 2 2 1]) ), ylim( [1 1 2 2]), [.98 .98 .98]);
                    plot( ax( not( cb, 1):not( cb, 2) ), xspec( not( cb, 1):not( cb, 2) ), '-', 'Color', [.85 .85 .85])
                end
%                 plot( ax( [1 end]), ones( 2, 1) * ylim, 'k-')
            end
            
            %set( gca, 'Color', [.95 .95 .95], 'Box', 'on')

            axis( [min( ax) max( ax) ylim ylim])                       
            
            %Add the legends manually :( (see reasoning above)
            subplot( 4, 1, 2)
            cla
            hold on
            for cg = 1:5
                plot( [1 2], 6 - [cg cg], mr{cg}, 'LineWidth', 2, 'Color', col( cg, :) )
            end
            text( ones( 5, 1) * 2.5, 6 - ( (1:5)'), leg)
            set( gca, 'XTick', [], 'YTick', [], 'Box', 'on')
            title( titname, 'FontSize', 12, 'FontWeight', 'Demi')
            if num == 1
                axis( [.6 3.5 .3 5.7] )
            else
                axis( [.6 5.0 .3 5.7] )
            end
            
            if opt.plot(1) == 4
                %Add # of fac legend
                subplot( 4, 1, 4)
                for cg = 1:5
                    patch( [1 2 2 1], 6 - [cg cg cg + .7 cg + .7], cback( cg, :) )
                end
                leg{ 5} = [num2str( limfac(end)) '+'];
                for cg = 1:4
                    leg{cg} = num2str( limfac( cg) );
                    if limfac( cg) + 1 < limfac( cg + 1)
                        leg{ cg} = [leg{ cg} '-' num2str( limfac( cg + 1) - 1)];
                    end
                end
                text( ones( 5, 1) * 2.5, 6 - ( (1:5)' + .35), leg)
                axis( [.6 3.6 0 5.3] )
                title( '# of factors', 'FontSize', 10, 'FontWeight', 'Demi')
                set( gca, 'XTick', [], 'YTick', [], 'Box', 'on')
            end
            
        end
               
        %Add the number of misclassifiations for the whole spectra
        subplot( 4, 1, 3)
        patch( [0 0 1 1], [0 1 1 0], col(axmis( end), :) )
        if opt.plot( 2) == 1 && ~isnan( opt.nLV( end) )
            if length( fjernlike( opt.nLV( end, :)' ) ) == 1
                tempLV = opt.nLV( end, 1);
                txtf = '%2.0f';
            else
                tempLV = opt.nLV( end, :);
                txtf = txtfor;
            end
            g = text( .5, .5, num2str( tempLV, txtf), 'HorizontalAlignment', 'center');
            if opt.plot( 1) == 2 && axmis( end) == 5 %Then the box is black and the text has to be white
                set( g, 'Color', [1 1 1])
            end
        end
        set( gca, 'XTick', [], 'YTick', []) %Remove the ticks
        title( 'Full spectra', 'FontWeight', 'Demi')
        
        %Move the subplots into the right positions
        g = gcf;
        h = get( g, 'Children');
        if opt.plot(1) == 4
            pos = [.82 .2 .15 .05; .82 .4 .15 .2; .82 .7 .15 .2; .1 .1 .7 .8];
        else
            pos = [.82 .4 .15 .05; .82 .7 .15 .2; .1 .1 .7 .8];            
        end
        for cp = 1:size( pos, 1)
            set( h(cp), 'Position', pos( cp, :) )
        end
        
    otherwise
        error( 'That plot version has not yet been made')
end

g = gcf;
h = get( g, 'Children');
h = h(end);

%Need a text string at the bottom with nLV
if opt.plot(2) == 1 && ~isnan( sum( opt.nLV(:) ) )
    axes( h)
    g = text( axid, ones( length( axid), 1) * (ylim(1) + range( ylim) * .0375), ...
        num2str( opt.nLV( 1:end-1, :), txtfor), 'HorizontalAlignment', 'center', ...
        'FontSize', fs, 'Rotation', 90);
    if opt.plot( 1) == 2
        set( g( axmis(1:end-1) == 5), 'Color', [1 1 1], 'HorizontalAlignment', 'Center', 'FontSize', fs)
    end
end

%This part is needed in case the intervals are continous
if opt.plot(1) == 2
    set( h, 'Color', [0 0 0])
end

%Reverse the x-axis if it is NMR/ IR spectra
if ax(1) > ax(end) || opt.dir ~= 1
    set( h, 'XDir', 'reverse')
end
axes( h)

%--------------------------------------------------------------------------
function xnew=vec(x)
%xnew=vec(x)
%
% Vectorices the matrix x

xnew=x(:);

%--------------------------------------------------------------------------
function [svar,ind]=fjernlike(x,n,lim)
%[svar,ind]=fjernlike(x,n,lim)
%
% Removes equal rows with respect to a specific column. NaN-values will not
% be counted as a unique value
%
% INPUT:
%  x      The x-matrix, can also be a string matrix
%  n      The column to check by 
%  lim    Can be added to specify a %-deviation allowed for equality
%
% OUTPUT:
%  svar   The matrix consisting of the unique values in column n
%  ind    The indexes for the equal rows

if isempty(x)
    svar=[];
    ind=[];
    return
end
if iscell(x)
    x=char(x);
end
if nargin<3
    lim=0;
end
[~,k]=size(x);
if nargin==1
    n=1:k;
end
if ischar(x(1,:))
    [xnew,i]=sortrows(x,n);
    a=double(xnew);
    c=a(1:end-1,:)~=a(2:end,:);
%     d=a(1:end-1,:)==a(2:end,:);
    %If there is only a character vector, this will also check for those
    if size(c,2)>1
        b=[1 find(sum(c, 2)>0)+1]';
    else
        b=[1;find(c>0)+1]; 
    end
else
    [xnew,i]=sortrows(x(:,n));
    a=xnew(1:end-1,:)*(1+lim)-xnew(2:end,:)*(1-lim);
    a(a>0)=0;
    a=a.^2;
    if length(n)>1
        a = sum(a, 2)';
    end
    b=[1;find(a~=0)+1];
end
svar=x(sort(i(b)),:);

if nargout==2
    %Find the unique rows
    a=zeros(size(x,1),1);
    a(b)=1;
    c=cumsum(a);
    %Find the equal rows
    a(a==1)=b-[0;b(1:end-1)];
    d=(1:length(a))'-cumsum(a)+1;
    [~,p]=sort(i(b));
    %Convert that into a row index with original coordinates
    ind=full(spconvert([c d i]));
    ind(ind==0)=NaN;
    %It is utterly important to include the following line!
    %If not the index of row 12 will NOT(!) be the unique row given in svar
    ind=ind(p,:);
end

%--------------------------------------------------------------------------
function x=range(y)
% x=range(y)
%
% Find the range in the matrix/ vector y
%
% INPUT:
%  y   Matrix/ vector for investigation
%
% OUTPUT:
%  x   Range in y given as a row vector

% 280906 AAR

x=nanmax(y)-nanmin(y);