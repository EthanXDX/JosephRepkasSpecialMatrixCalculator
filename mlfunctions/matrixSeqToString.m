% matrixSeq: a row array wherein each entry is a matrix
function stringRep = matrixSeqToString(matrixSeq)

stringRep = "";

for i = 1:length(matrixSeq)
    if i > 1
        stringRep = stringRep + newline;
    end
    stringRep = stringRep + matrixToString(matrixSeq{i});
end
