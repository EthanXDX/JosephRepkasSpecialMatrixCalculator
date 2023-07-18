function stringRep = toString(mat);
% We need to manually remove the spaces from each entry,
% particularly when it comes to symbolic entries such as "a + 1"
matOfStrings = replace(string(mat), " ", "");
% This is a hacky but quick way of joining
% entries with spaces, joining lines with semicolons,
% and wrapping the whole things withs square bracket.
% It should probably be replaced% with more careful code ...
stringRep = replace(mat2str(matOfStrings), '"', '');
