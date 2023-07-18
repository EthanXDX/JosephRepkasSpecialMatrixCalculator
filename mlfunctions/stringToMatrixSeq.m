function matrixSeq = stringToMatrixSeq(matrixSeq)

matStrings = matrixSeq.split(newline);
matrixSeq = cell(0);
for i = 1:length(matStrings)
    matrixSeq{i} = stringToMatrix(matStrings{i});
end
