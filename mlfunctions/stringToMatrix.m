function mat = stringToMatrix(string)

% NOTE: The following also changes fractions into doubles ...
% saves time to avoid unnecessary symbols
%if hasSymbols(string)
%    mat = str2sym(string)
%else
%    mat = str2num(string)
%end

mat = str2sym(string);

function val = hasSymbols(string)
    matches = regexp(string, "[a-zA-Z]");
    val = length(matches) > 0;

