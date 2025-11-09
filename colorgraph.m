function out = colorgraph( X, Y, opt)
% out = colorgraph( X, Y, opt);
%Colors the samples (X) according to a value (given in Y)
%
%INPUT:
% X   Samples x Variables
% Y   Vector of values to color according to
% opt Options. Type 'opt = colorgraph;' for default values

%170915 AAR

if nargin == 0
    out.ax = NaN;
    out.axdir = true;
    out.map = 'jet';
    out.legend = 'none';
    return
end

if size( Y, 2) > 1
    error( 'Y can only have one column');
end

if size( X, 1) ~= size( Y, 1)
    error( 'X and Y must have the same number of rows')
end

if nargin < 3
    opt = colorgraph;
end

if isnan( opt.ax)
    opt.ax = 1:size( X, 2);
end

% pt = 1;
if length( opt.ax) ~= size( X, 2)
%     if length( opt.ax == size( X, 1) )
%         %Different plot being made
%         pt = 2;
%     else
        opt.ax = 1:size( X, 2);
%     end
end

%Define the colors
try
    col = eval( [ 'colormap(' opt.map '(100) );']);
catch
    if isnumeric( opt.map) && size( opt.map, 2) == 3
        col = opt.map;
    else
        error( 'opt.map has not been set to a valid value')
    end
end

%061118 AAR A different approach
% Y = repnum( Y);
Y = Y - min( Y);
Y = Y/ range( Y);
Y = floor( Y * ( size( col, 1) - 1) ) + 1;

% clf
if iscell( opt.legend)
    if length( fjernlike( Y) ) ~= length( opt.legend)
        error( 'The number of unique ''Y'' and entries in opt.legend is different')
    end
    uid = sort( fjernlike( Y) );
    hold on
    for cu = 1:length( uid)
        id = find( Y == uid( cu), 1);
        g( cu) = plot( opt.ax, X( find( Y == uid(cu), 1), :), '-', 'Color', col( uid( cu), :) );
    end    
end

h = plot( opt.ax, X);

for cs = 1:size( X, 1)
    set( h(cs), 'Color', col( Y( cs), :) );
end

if ~opt.axdir 
    set( gca, 'XDir', 'reverse')
end

if iscell( opt.legend)
    legend( opt.legend)
end

out = col;