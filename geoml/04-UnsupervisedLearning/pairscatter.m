function pairscatter(X,label,grp,varargin)
% Pairwise Scatter of Columns of Matrix
% This is similar to plotmatrix except that it only plots each pair of
% variables once, eliminates the histograms on the diagonals and enables
% the user to also specify labeling for each axes.  
%
% pairscatter(X, label, grp)
%    Draws a pairwise scatter of the columns of X. label is an optional
%    cell array of labels, one for each column of X. grp is an optional
%    grouping variable similar to that used in gscatter
%
% pairscatter(X, label, grp, 'plotArgs', {'name',value,...})
%    specifies a list of parameter-value pairs to be passed on to the
%    plotting routine to customize the look of the markers in the scatter
%    plot. Default = {'marker', '.', 'markersize', 3}
%
% Copyright 2014-2019 MathWorks, Inc.
figure
param = parsePVPairs({'plotArgs','showDensity'}, {{},false}, {@iscell, @islogical}, varargin); 

param.plotArgs = parsePVPairs({'marker','markersize','linestyle'}, {'o', 2,'none'}, [], param.plotArgs, 'keepUnmatched', true, 'pvOutput', true);


m = size(X,2);
nCharts = nchoosek(m,2);
nCol = ceil(sqrt(nCharts));
nRow = ceil(nCharts/nCol);
ctr = 0;

if nargin < 2 || isempty(label), label = 1:m; end
if nargin < 3 || isempty(grp), grp = []; end
    
if isnumeric(label)
    label = cellfun(@num2str,num2cell(label),'uniformoutput',false);
end
args = {};
if nRow > 5 || nCol > 5
    args = {'xtickmode','manual','ytickmode','manual','xtick',[],'ytick',[]};
end
for i = 1:m
    for j = i+1:m
        ctr = ctr + 1;
        ax(ctr) = subplot(nRow,nCol,ctr,'fontsize',8,args{:});
           if isempty(grp)
                plot(X(:,i), X(:,j), param.plotArgs{:});
            else
                h = gscatter(X(:,i), X(:,j), grp);
                legend location best
                set(h,param.plotArgs{:})
            end
        
        axis(ax(ctr),'tight');
        %xlabel(ax(ctr),label{i});
        %ylabel(ax(ctr),label{j});
        title(sprintf('%s vs %s',label{i},label{j}));
        drawnow
    end
end

function param = parsePVPairs(names, defaults, validators, pvs, varargin) 
% Parse Parameter Value inputs & return structure of parameters
% param = parsePVPairs(names, defaults, validators, pvs) 
%
% All inputs are cell arrays. names contains the names of the parameters,
% defaults - their default values, validators - the function handles each
% of which validate the input value, and pvs - the set of parameter value
% inputs. param is a structure with a field corresponding to each named
% parameter.
%
% param = parsePVPairs(..., 'caseSensitive', true) enforces a case sensitive
% match. Default = false
%
% param = parsePVPairs(..., 'partialMatch', false) enforces a full
% parameter match. Default = true
%
% param = parsePVPairs(..., 'keepUnmatched', false) enforces that all
% paramters match those listed in names. Default = false
%
% param = parsePVPairs(..., 'pvOutput', true) returns output as a cell
% array of parameter value pairs. Default = false 

% Copyright 2014 MathWorks, Inc.

if ~isempty(varargin) % To prevent infinite recursion
    localParam = parsePVPairs({'caseSensitive','partialMatch','keepUnmatched','pvOutput'},...
        {false, true, false, false},{@islogical, @islogical, @islogical, @islogical}, varargin);
else
    localParam.caseSensitive = false;
    localParam.partialMatch = true;
    localParam.keepUnmatched = true;
    localParam.pvOutput = false;
end

if nargin < 2 || isempty(defaults), defaults = repmat({[]},1,length(names)); end
if nargin < 3 || isempty(validators), validators = repmat({{}},1,length(names)); end    

p = inputParser;
for i = 1:length(names)
    if ~isempty(validators{i})
        p.addParamValue(names{i},defaults{i},validators{i});
    else
        p.addParamValue(names{i},defaults{i});
    end
end
p.CaseSensitive = localParam.caseSensitive;
p.PartialMatching = localParam.partialMatch;
p.KeepUnmatched = localParam.keepUnmatched;
p.parse(pvs{:});
param = p.Results;

if localParam.pvOutput
    param = [fieldnames(param)';struct2cell(param)'];
    param = param(:)';
end
