#include "../lib/src/headers/lin_alg.h"
#include "../lib/src/headers/lie_algebra.h"
#include "../lib/src/headers/compat.h"

#include <iostream>
#include <fstream>

#ifndef MAIN_H
#define MAIN_H

string RawInputGetter(int i);
string BasisGetter(int i);
string DimGetter(int i);
string NormalizerGetter(int i);
string NormalizerDimGetter(int i);
string CentralizerGetter(int i);
string CentralizerDimGetter(int i);
void append_column(string heading, std::function<string(int)> AlgebraInvariantGetter);
string sanitize_csv_cell(string text);
string get_file_contents(string file_path);

class InputParser {
    public:
        InputParser (int &argc, char **argv);
        const string& get_cmd_option(const string &option) const;
        bool cmd_option_exists(const string &option) const;
    private:
        vector<string> tokens;
};

string get_cmd_option_fallback(InputParser input, string option, string fallback);

void print_usage_message();

#endif // MAIN_H