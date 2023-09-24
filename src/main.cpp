#include "main.h"

vector<lie_algebra*> algebraSeq;
vector<string> algebra_names;
vector<string> algebra_inputs;
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

    vector<int> input_indices = {};

    //Look for flags that modify what is done
    bool derived = input.cmd_option_exists("--derived");
    bool lower = input.cmd_option_exists("--lower");
    bool include_basis = input.cmd_option_exists("--include-basis");
    bool include_inputs = input.cmd_option_exists("--include-inputs");
    bool compute_all_min_ranks = input.cmd_option_exists("--compute-all-min-rank");
    bool compute_top_min_rank = input.cmd_option_exists("--compute-min-rank");

    if (derived || lower) {
        vector<lie_algebra*> new_algebraSeq = vector<lie_algebra*>();
        for (int i = 0; i < algebraSeq.size(); i++) {//(lie_algebra* alg : algebraSeq) {
            lie_algebra* alg = algebraSeq[i];
            new_algebraSeq.push_back(alg);

            vector< lie_algebra* > derived_series = alg->compute_derived_series();
            if (derived) {
                new_algebraSeq.insert(new_algebraSeq.end(), derived_series.begin(), derived_series.end());
            }
            vector< lie_algebra* > lower_series = alg->compute_lower_central_series();
            if (lower) {
                new_algebraSeq.insert(new_algebraSeq.end(), lower_series.begin(),lower_series.end());
            }

            algebra_names.push_back(std::string ("Input ") +  std::to_string(i+1));
            algebra_inputs.push_back(rawAlgebraStrings[i]);
            if (derived) {
                for (int j = 0; j < alg->compute_derived_series().size(); j++) {
                    algebra_names.push_back(std::string ("Derived ") +  std::to_string(j+1));
                    algebra_inputs.push_back("");
                }
            }
            if (lower) {
                for (int j = 0; j < alg->compute_lower_central_series().size(); j++) {
                    algebra_names.push_back(std::string ("Lower ") +  std::to_string(j+1));
                    algebra_inputs.push_back("");
                }
            }
        }
            algebraSeq = new_algebraSeq;
    } else {
        for (int i = 0; i < algebraSeq.size(); i++) {
            algebra_names.push_back(std::string ("Input ") +  std::to_string(i+1));
            algebra_inputs.push_back(rawAlgebraStrings[i]);
        }
    }


    // Initialize lines (makes later code simpler)
    for (int i = 0; i < algebraSeq.size() + 1; i++) {
        lines.push_back("");
    }

    append_column("Names", NameGetter);
    if (include_inputs) append_column("Raw Input", RawInputGetter);
    if (include_basis) append_column("Basis", BasisGetter);
    append_column("Dimension", DimGetter);
    if (include_basis) append_column("Normalizer", NormalizerGetter);
    append_column("Normalizer Dimension", NormalizerDimGetter);
    if (include_basis) append_column("Centralizer", CentralizerGetter);
    append_column("Centralizer Dimension", CentralizerDimGetter);
    if (compute_all_min_ranks) append_column("Minimum rank", MinimumRankGetter);
    else if (compute_top_min_rank) append_column("Minimum rank", TopMinimumRankGetter);
    append_column("Maximum rank", MaximumRankGetter);

    std::ofstream output;
    output.open(output_path);
    for (int i = 0; i < lines.size(); i++) {
        output << lines[i] << '\n';
    }
    output.close();

    return 0;
}

string NameGetter(int i) {
    return algebra_names[i];
}

string RawInputGetter(int i) {
    return algebra_inputs[i];
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
    int dim = algebraSeq[i]->compute_normalizer()->get_dim();
    return std::to_string(dim);
}

string CentralizerGetter(int i) {
    return toString(algebraSeq[i]->compute_centralizer());
}

string CentralizerDimGetter(int i) {
    int dim = algebraSeq[i]->compute_centralizer()->get_dim();
    return std::to_string(dim);
}

string MinimumRankGetter(int i) {
    int min_rank = algebraSeq[i]->min_rank();
    return std::to_string(min_rank);
}

string TopMinimumRankGetter(int i) {
    if (algebra_names[i][0] != 'I') return "";
    int min_rank = algebraSeq[i]->min_rank();
    return std::to_string(min_rank);
}

string MaximumRankGetter(int i) {
    int max_rank = algebraSeq[i]->max_rank();
    return std::to_string(max_rank);
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