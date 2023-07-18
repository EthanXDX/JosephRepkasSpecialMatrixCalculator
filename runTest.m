function returnValue = runTest()

addpath('./joemat')
addpath('./mlfunctions')

L5_1_generators = {[0.25 1 0 0;0 0.25 0 0;0 0 0.25 0;0 0 0 -0.75],[0 1 0 0;0 0 1 0;0 0 0 0;0 0 0 0],[0 0 1 0;0 0 0 0;0 0 0 0;0 0 0 0]}
L5_2_generators = {[0.25 0 0 0;0 0.25 0 1;0 0 -0.75 0;0 0 0 0.25],[0 1 0 0;0 0 0 0;0 0 0 0;0 0 0 0],[0 0 0 1;0 0 0 0;0 0 0 0;0 0 0 0]};
a=sym("a");
L5_4_generators = {[0 1 0 0;0 0 1 0;0 0 0 1;0 0 0 0],[0 0 a 0;0 0 0 a+1;0 0 0 0;0 0 0 0],[0 0 0 1;0 0 0 0;0 0 0 0;0 0 0 0]};
L5_5_generators = {[0 1 0 0;0 0 0 0;0 0 0 1;0 0 0 0],[0 0 a 0;0 0 0 a+1;0 0 0 0;0 0 0 0],[0 0 0 1;0 0 0 0;0 0 0 0;0 0 0 0]};

%generatorsString = "[0.25 1 0 0;0 0.25 0 0;0 0 0.25 0;0 0 0 -0.75]" + newline + "[0 1 0 0;0 0 1 0;0 0 0 0;0 0 0 0]" + newline + "[0 0 1 0;0 0 0 0;0 0 0 0;0 0 0 0]";
%clib.joemat.setLieAlgebra(generatorsString);
%disp("Generators:" + newline + generatorsString);
%disp("Dimension: " + clib.joemat.getLieAlgebraDim());
%disp("Basis: " + newline + clib.joemat.getLieAlgebraBasis());
%disp("Normalizer: " + newline + clib.joemat.getLieAlgebraNormalizer());
%disp("Centralizer: " + newline + clib.joemat.getLieAlgebraCentralizer());


generators = L5_4_generators;
disp("Generators:" + newline + matrixSeqToString(generators));
disp("Dimension: " + getLieAlgebraDim(generators));
disp("Basis: " + newline + matrixSeqToString(getLieAlgebraBasis(generators)));
disp("Normalizer: " + newline + matrixSeqToString(getLieAlgebraNormalizer(generators)));
disp("Centralizer: " + newline + matrixSeqToString(getLieAlgebraCentralizer(generators)));

returnValue = 0;
