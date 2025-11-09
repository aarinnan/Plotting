function op = iplsplot( modi, modf, opt)
%iplsplot( modi, modf, opt)
% Make the "normal" i-PLS plot for calibration
%
%INPUT:
% modi  Model structure from 'ipls'
% modf  Model structure from 'cv', all interval model
% opt   Options, type opt = iplsplot for default values
%
% See also: icplot

%281212 AAR

if nargin == 0
    op.Ax = NaN;
    op.Xmean = NaN;
    op.ind = NaN;
    op.nLV = NaN;
    op.plot = [3 1];
    op.dir = 1;
    op.info = {'Ax = the axis for your spectra'; ...
        'Xmean = Average spectra of the samples'; ...
        'ind = Vector indicating intervals, same numbers = same interval'; ...
        'plot = Two numbers. First number is plot version. 1-3 can all be given'; ...
        '       Second number indicates with or without number of LV''s indicated'; ...
        'dir = direction of x-axis. 1 - normal, other - reverse'};    
    return
end

num = length( modi);
lim = modi(1).Val.rms(1);
if isnan( opt.nLV) || length( opt.nLV) ~= length( fjernlike( opt.ind) ) + 1
    for ci = 1:num
        nLV( ci) = estfac( modi( ci) );
        rms( ci) = min( lim, modi( ci).Val.rms( nLV( ci) + 1) );
    end
    
    nLV( ci + 1) = estfac( modf);
    rms( ci + 1) = min( lim, modf.Val.rms( nLV( ci + 1) ) );
else
    nLV = opt.nLV
    for ci = 1:num
        rms( ci) = modi( ci).Val.rms( opt.nLV( ci) );
    end
    rms( ci + 1) = modf.Val.rms( opt.nLV( ci + 1) );
end

nLV( rms == lim) = 0;

%Find the range of the two axis
temp = [min( rms) max( rms)];
ylim = temp + [-1 1] * range( temp)/ 20;
xlim = [min( opt.Ax) max( opt.Ax)];

%Find the center of each interval
[~, j] = fjernlike( opt.ind);
i = nanmedian( j, 2);
temp = [floor( i) ceil( i)];
tempax = mean( opt.Ax( temp), 2);

%Plotting of the actual figure
figure
g = bar( tempax, rms( 1:end-1));
set( g, 'FaceColor', [.7 .7 .7])

%Adjust the spectra to fit the current plot
xtemp = opt.Xmean - min( opt.Xmean);
xtemp = xtemp./ max( xtemp); %Now the spectra goes from 0-1
%The spectra should cover the middle 75% of the plot (i.e. 12.5% empty
%space on each side
xtemp = xtemp * (range( ylim) * .75) + ylim(1) + range( ylim)/ 8;

%Plot the average spectra
hold on
h = plot( opt.Ax, xtemp, 'k-');
f = plot( xlim, rms( [end end]), 'r-');
set( h, 'LineWidth', 2)
set( f, 'LineWidth', 2)

%Swap the direction of the x-axis if set by the user
if opt.dir ~= 1    
    set( gca, 'Xdir', 'reverse')
end
if opt.dir == 1
    text( xlim(2), rms(end), num2str( nLV( end) ) )
else
    text( xlim(1), rms( end), num2str( nLV( end) ) )
end
axis( [xlim ylim])
t = text( tempax, ones( length( tempax), 1) * (ylim(1) + range( ylim)/ 20), num2str( vec( nLV( 1:end-1) ) ) );
set( t, 'HorizontalAlignment', 'center')

%--------------------------------------------------------------------------
function x=range(y, dim)
% x=range(y, dim)
%
% Find the range in the matrix/ vector y
%
% INPUT:
%  y   Matrix/ vector for investigation
%  dim In what dimension the range should be calculated <default = 1>
%
% OUTPUT:
%  x   Range in y

% 190115 AAR Included dimensionality
% 280906 AAR

if nargin < 2
    dim = 1;
end

if dim ~= 1
    dim = 2;
end

x = nanmax( y, dim) - nanmin( y, dim);

%--------------------------------------------------------------------------
function [a,i]=nanmax(x, dim, opt)
% [a,i]=nanmax(x, dim, opt)
%
% Calculate the max without taking into account the missing values.
%
% INPUT
%  x    X-matrix/vector
%  dim  Dimension to search for the max-value (default = 1)
%  opt  If opt=0 Inf is the max-value. opt=1 indicates the max number not
%        taking Inf into account
%        Default: opt = 1
%
% OUTPUT
%  a  Max-values
%  i  The index
%
% See also: isnan, max

% 131211 AAR Included the direction as well
% 121005 AAR I am not normally interested in Inf as the maximum
% 010905 AAR

if nargin < 2
    dim = 1;
end
if size( x, 1) == 1
    dim = 2;
end

if nargin < 3
    opt=1;
end
x(isnan(x))=-Inf;
if opt==1
    x(isinf(x))=-Inf;
end
[a,i]=max(x, [], dim);
if opt==1
    a(isinf(a))=NaN;
end

%--------------------------------------------------------------------------
function [a,i]=nanmin(x, dim, opt)
% [a,i]=nanmin(x, dim, opt)
%
% Calculate the min without taking into account the missing values.
%
% INPUT
%  x    X-matrix/vector
%  dim  Direction, 1 - columnwise (default), 2 - rowwise
%  opt  opt=1: -Inf are removed.
%        Default: opt = 1
%
% OUTPUT
%  a    Min-values
%  i    Index
%
% See also: nanmax, nanmean, nanmedian, nanstd, nansum, nanvar

% 180111 AAR Can now also calculate the minimum value in a specific
%             direction
% 121005 AAR Normally not interested in -Inf values
% 010905 AAR

if nargin<2
    opt=1;
    dim = 1;
end
if nargin < 3
    opt = 1;
end
if size( x, 1) == 1
    dim = 2;
end

if dim < 1 || dim > 2
    error( '''dim'' can only have the value ''1'' or ''2''')
end

x(isnan(x))=Inf;
if opt==1
    x(isinf(x))=Inf;
end
[a,i]=min(x, [], dim);
if opt==1
    a(isinf(a))=NaN;
end