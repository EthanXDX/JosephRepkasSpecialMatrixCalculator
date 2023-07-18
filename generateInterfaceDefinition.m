mex -setup cpp;
libName = "joemat";
libFile = "bin/libjoemat.dll";
switch computer('arch')
    case 'win64'
    case 'glnxa64'
        libFile = "bin/libjoemat.so";
    case 'maci64'
        libFile = "bin/libjoemat.dylib";
    otherwise
        warning('Unexpected architecture. Assuming Windows')
end
headerFile = "include/joemat.h";
clibgen.generateLibraryDefinition(headerFile, ...
    "Libraries", libFile, ...
    "PackageName", libName, ...
    "OverwriteExistingDefinition", true, ...
    "AdditionalCompilerFlags", "-std=c++17 -arch x86_64", ...
    "Verbose", true, ...
    "TreatConstCharPointerAsCString", true);
