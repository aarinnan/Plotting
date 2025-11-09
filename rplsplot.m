function rplsplot( X, ax, evolve, lab)
%rplsplot( X, ax, evolve, lab)
%
%INPUT:
% X       Pre-processed X data
% ax      Axis of X data
% evolve  Results from 'rpls'
% lab     Label on the x-axis
%
%See also: rpls

% Copyright, 2017 - 
% This M-file and the code in it belongs to the holder of the
% copyrights and is made public under the following constraints:
% It must not be changed or modified and code cannot be added.
% The file must be regarded as read-only. 
%
% Åsmund Rinnan
% Chemometrics and Analytical Technology
% Department of Food Science
% Faculty of Science
% University of Copenhagen
% Rolighedsvej 30, 1958 Frederiksberg C, Denmark
% e-mail: aar@food.ku.dk

%AAR 010813

temp = evolve.w( 2:end, :);
temp( :, end + 1) = temp( :, end);
temp( end + 1, :) = temp( end, :);
templog = log10( temp + 1e-35);

figure
h = surf( ax( [1 1:end]), (1:size( templog, 1)), templog);
set( h, 'EdgeColor', 'flat')
%Define a new set of colors
jet = colormap;
col = [(.9:-.02:.71)' * ones( 1, 3); jet( 25:-1:1, :)];
colormap( col)
view( [0 90]);

%Include the performance axis (RMSECV)
rmsid = (1:size( evolve.detail.id, 1))' * ones( 1, size( evolve.detail.id, 2) );
rmsid( isnan( evolve.detail.id) ) = NaN;
rmsid = vec( rmsid');
rmsid( isnan( rmsid) ) = [];
%AAR 121114 Don't know what happens here, really
% rms = evolve.rms( rmsid( 2:end-1));
rms = evolve.rms( rmsid( 2:end) );

%If class
% rmstxt = num2str( rms');
%else
rmstxt = num2str( rms', '%4.4f');

%Get it into the same axis scale
rms = 2 - rms;
rms = rms - min( rms);
rms = rms/ max( rms);

%In order to get different colors for the rms-axis, I do a bit of a trick
[i, k] = sort( rms');
temp = [diff( i) zeros( size( i, 1) - 1, 1)];
for ce = 1:( 64 - length( temp) - 1)
    [~, m] = max( temp( :, 1));
    temp( m, 2) = temp( m, 2) + 1;
    temp( m, 1) = temp( m, 1)/ 2;
end
temp = [temp; zeros( 1, 2)];
temp( end, 3) = 64;
for ct = (length( rms) - 1):-1:1
    temp( ct, 3) = temp( ct + 1, 3) - 1 - temp( ct, 2);
end
[~, j] = fjernlike( i);
id = find( sum( ~isnan( j), 2) > 1);
for ci = 1:length( id)
    temp( j( id( ci), ~isnan( j( id( ci), :) ) ), 3) = temp( j( id( ci), 1), 3);
end
rms( k) = temp( :, 3);

%Need to get these values into the 64 bins which currently are in use
if ax(1) > ax(end)
    xlim = [ax(end) ax(end) - ( max( ax) - min( ax))/ 32];
else
    xlim = [ax(end) ax(end) + ( max( ax) - min( ax))/ 32];
end
ylim = [1 size( templog, 1)];
h = patch( xlim( [1 1 2 2]), ylim( [1 2 2 1]), [1 1 1]);
set( h, 'EdgeColor', [1 1 1])
g = colorbar( 'WestOutside');
set( g, 'YTick', -35:5:0, 'YTickLabel', num2str( 10.^( -35:5:0)' ) )
col = [linspace( .5, 1, 64)' linspace( 0, 1, 64)' * ones( 1, 2)];
if ax(1) > ax(end)
    xlim = [xlim(2) xlim(1) - ( max( ax) - min( ax))/ 32 * 3];
else
    xlim = [xlim(2) xlim(1) + ( max( ax) - min( ax))/ 32 * 3];
end
for cr = 1:length( rms)
    patch( xlim( [1 2 2 1]), [cr cr cr + 1 cr + 1], col( rms( cr), :) )
end

if ~iscell( lab)
    lab = {lab};
end

xlabel( lab{1}, 'FontSize', 14, 'FontWeight', 'Demi')
ylabel( 'Iteration number', 'FontSize', 14, 'FontWeight', 'Demi')
if ax(1) > ax(end) %Typical for NMR data at least
    text( ones( length( rmstxt), 1) * ( ax(end) - ( max( ax) - min( ax))/ 10), (1:length( rmstxt) )' + .5, rmstxt)
%     text( (ax(end) - range( ax)/ 32), ylim(2) + range( ylim)/ 20, 'MisClass', 'FontSize', 12, 'FontWeight', 'Demi', 'HorizontalAlignment', 'center')
    text( (ax(end) - ( max( ax) - min( ax))/ 32), ylim(2) + ( max( ylim) - min( ylim))/ 20, 'RMSECV', 'FontSize', 12, 'FontWeight', 'Demi', 'HorizontalAlignment', 'center')
else
    text( ones( length( rmstxt), 1) * ( ax(end) + ( max( ax) - min( ax))/ 10), (1:length( rmstxt) )' + .5, rmstxt)
%     text( (ax(end) + range( ax)/ 32), ylim(2) + range( ylim)/ 20, 'MisClass', 'FontSize', 12, 'FontWeight', 'Demi', 'HorizontalAlignment', 'center')
    text( (ax(end) + ( max( ax) - min( ax))/ 32), ylim(2) + ( max( ylim) - min( ylim))/ 20, 'RMSECV', 'FontSize', 12, 'FontWeight', 'Demi', 'HorizontalAlignment', 'center')
end

%Also need the average spectra on top
if min( size( X) ) > 1
    xmean = mean( X);
else
    xmean = X;
end
xmean = xmean - min( xmean);
xmean = xmean/ max( xmean);
const = find( sum( evolve.w > 1e-4, 2) < evolve.num(2)/ 4, 2);
const(1) = [];
const(2) = ( ylim(2) - const( 1) ) * .75;
xmean = xmean * const(2) + const(1);
hold on
plot( ax, xmean, '-', 'Color', [.8 0 0], 'LineWidth', 2)
if ax(1) > ax(end)
    axis( [ax(end) - ( max( ax) - min( ax))/ 32 * 3 ax(1) 1 ylim(2)])
else
    axis( [ax(1) ax(end) + ( max( ax) - min( ax))/ 32 * 3 1 ylim(2)])
end

%Add the optimal model as a dotted line

ytick = get( gca, 'YTick');
set( gca, 'YTick', ytick + .5, 'YTickLabel', num2str( ytick') )
[~, op] = min( evolve.rms);
plot( ax( [1 end]), evolve.detail.id( op, [1 1]) - .5, 'r--', 'LineWidth', 2) 
if ax(1) > ax(end)
    set( gca, 'XDir', 'reverse')
end