#include "main.h"

vector<lie_algebra*> algebraSeq;
vector<string> rawAlgebraStrings;
vector<string> lines;

int main(int argc, char* argv[]) {
    InputParser input(argc, argv);
    if (input.cmd_option_exists("-h")) {
        print_usage_message();
        return 0;
    }
    string input_path = get_cmd_option_fallback(input, "-i", "matrixInput.txt");
    string output_path = get_cmd_option_fallback(input, "-o", "output.csv");

    string matrixInput = get_file_contents(input_path);
    matrixInput = replace(matrixInput, "\r", "");
    compress_spaces_and_newlines(matrixInput);
    algebraSeq = toLieAlgebraSequence(matrixInput);
    rawAlgebraStrings = split(matrixInput, "\n@\n");

    // Initialize lines (makes later code simpler)
    for (int i = 0; i < algebraSeq.size() + 1; i++) {
        lines.push_back("");
    }

    append_column("Raw Input", RawInputGetter);
    append_column("Basis", BasisGetter);
    append_column("Dimension", DimGetter);
    append_column("Normalizer", NormalizerGetter);
    append_column("Normalizer Dimension", NormalizerDimGetter);
    append_column("Centralizer", CentralizerGetter);
    append_column("Centralizer Dimension", CentralizerDimGetter);

    std::ofstream output;
    output.open(output_path);
    for (int i = 0; i < lines.size(); i++) {
        output << lines[i] << '\n';
    }
    output.close();

    return 0;
}

string RawInputGetter(int i) {
    return rawAlgebraStrings[i];
}

string BasisGetter(int i) {
    return toString(algebraSeq[i]->get_basis());
}

string DimGetter(int i) {
    return std::to_string(algebraSeq[i]->get_dim());
}

string NormalizerGetter(int i) {
    return toString(algebraSeq[i]->compute_normalizer());
}

string NormalizerDimGetter(int i) {
    int dim = (algebraSeq[i]->compute_normalizer())->get_dim();
    return std::to_string(dim);
}

string CentralizerGetter(int i) {
    return toString(algebraSeq[i]->compute_centralizer());
}

string CentralizerDimGetter(int i) {
    int dim = (algebraSeq[i]->compute_centralizer())->get_dim();
    return std::to_string(dim);
}

void append_column(string heading, std::function<string(int)> AlgebraInvariantGetter) {
    if (!lines[0].empty()) lines[0] += ",";
    lines[0] += heading;
    for (int i = 0; i < algebraSeq.size(); i++) {
        if (!lines[i+1].empty()) lines[i+1] += ",";
        lines[i+1] += sanitize_csv_cell(AlgebraInvariantGetter(i));
    }
}

string sanitize_csv_cell(string text) {
    return replace(text, "\n", " ");
}

string get_file_contents(string file_path) {
    // see https://stackoverflow.com/a/2602258
    std::ifstream t(file_path);
    std::stringstream buffer;
    buffer << t.rdbuf();
    return buffer.str();
}

// see https://stackoverflow.com/a/868894
InputParser::InputParser (int &argc, char **argv) {
    for (int i=1; i < argc; ++i)
        this->tokens.push_back(std::string(argv[i]));
}

const string& InputParser::get_cmd_option(const string &option) const {
    std::vector<std::string>::const_iterator itr;
    itr =  std::find(this->tokens.begin(), this->tokens.end(), option);
    if (itr != this->tokens.end() && ++itr != this->tokens.end()) {
        return *itr;
    }
    static const std::string empty_string("");
    return empty_string;
}

bool InputParser::cmd_option_exists(const string &option) const {
    return std::find(this->tokens.begin(), this->tokens.end(), option)
        != this->tokens.end();
}

string get_cmd_option_fallback(InputParser input, string option, string fallback) {
    string value;
    const string &raw_value = input.get_cmd_option(option);
    if (raw_value.empty()) value = fallback;
    else value = raw_value;
    return value;
}

void print_usage_message() {
    std::cout << "usage: joemat -i matrixInput.txt -o output.csv\n";
}