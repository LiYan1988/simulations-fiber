function test_input(varargin)

checkPositiveScalar = @(x, variableName) assert(isnumeric(x) && ...
    isscalar(x) && (x>0), ...
    sprintf('%s should be a positive number.', variableName));

p = inputParser;
checkInputCell = @(x, variableName) assert(iscell(x) || ...
    (isnumeric(x) && isscalar(x) && (x>0)),...
    sprintf('%s should be either a cell or a positive number.', variableName));
addParameter(p, 'inputCell', cell(1), @(x)checkInputCell(x, 'InputCell'))
parse(p, varargin{:})

fprintf('%d\n', size(p.Results.inputCell, 2))
fprintf('%d\n', isfield(p.Results, 'ok'))