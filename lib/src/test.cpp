//
// Created by cole on 24/05/23.
//
#include "headers/utils.h"
#include "headers/lin_alg.h"
#include "headers/lie_algebra.h"
#include <random>

//TODO: Add tests for derived / central series, max_rank

lie_algebra* get_L5_1() {
    g::matrix a = {{1,4,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,-3}};
    g::matrix b = {{0,1,0,0},{0,0,1,0},{0,0,0,0},{0,0,0,0}};
    g::matrix c = {{0,0,1,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
    std::vector< g::matrix> matrices = {a,b,c};
    lie_algebra* out = new lie_algebra(matrices);
    return out;
}

lie_algebra* get_L5_2() {
    g::matrix a = {{1,0,0,0},{0,1,0,4},{0,0,-3,0},{0,0,0,1}};
    g::matrix b = {{0,1,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
    g::matrix c = {{0,0,0,1},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
    std::vector< g::matrix> matrices = {a,b,c};
    lie_algebra* out = new lie_algebra(matrices);
    return out;
}

lie_algebra* get_L5_4(g::symbol x) {
    g::matrix a = {{0,1,0,0},{0,0,1,0},{0,0,0,1},{0,0,0,0}};
    g::matrix b = {{0,0,x,0},{0,0,0,x+1},{0,0,0,0},{0,0,0,0}};
    g::matrix c = {{0,0,0,1},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
    std::vector< g::matrix> matrices = {a,b,c};
    lie_algebra* out = new lie_algebra(matrices);

    return out;
}

lie_algebra* get_L5_5(g::symbol x) {
    g::matrix a = {{0,1,0,0},{0,0,0,0},{0,0,0,1},{0,0,0,0}};
    g::matrix b = {{0,0,x,0},{0,0,0,x+1},{0,0,0,0},{0,0,0,0}};
    g::matrix c = {{0,0,0,1},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
    std::vector< g::matrix> matrices = {a,b,c};
    lie_algebra* out = new lie_algebra(matrices);
    return out;
}

lie_algebra* get_L3_025_5() {
    g::matrix m1 = {{-1, 0, 0, 0, 0, 2, 0, 0, 0, 0},
                    {0,  1, 0,  0, 0, 0, 0, 0, 0, 0},
                    {0,  0, -1, 0, 0, 0, 0, 0, 0, 0},
                    {0,  0, 0,  1, 0, 0, 0, 0, 0, 0},
                    {0,  0, 0,  0, 0, 0, 0, 1, 0, 0},
                    {0,  0, 1,  0, 0, -1, 0, 0, 0, 0},
                    {0,  0, 0,  0, 0, 0, 0, 0, 0, 1},
                    {0,  0, 0,  0, 0, 0, 0, 0, 0, 0},
                    {0,  0, 0,  0, 0, 0, 0, 0, 1, 0},
                    {0,  0, 0,  0, 0, 0, 0, 0, 0, 0}};
    g::matrix m2 = {{0, 0, 0, 0, 2, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 2, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 1, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 1, 1, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 1, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 1, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 1, 0, 0, 0, 0, 0, 0}};
    g::matrix m3 = {{0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 1, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
                    {0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
    mat_vec basis = mat_vec {m1,m2,m3};
    mat_vec brackets = {lin_alg::bracket(m1,m2), lin_alg::bracket(m2,m3), lin_alg::bracket(m1,m3)};
    g::matrix M = lin_alg::bracket(m1,m2).add(m2).sub(m3);
    lie_algebra* L3 = new lie_algebra(basis);
    return L3;
}

lie_algebra* get_L3_0_9a(g::symbol x){
    g::matrix m1 = {{0, 0, 0, 0, 0, 0, 2, 0, 0, 0}, {0, 2, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, -2, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 1, 0, 0, 0, 1, 0}, {0, 0, 0, 0, 0, -1, 0, 0, 0, 1}, {0, 0, 0, 1, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 1, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, -1}};
    g::matrix m2 = {{2*x, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 2*x+2, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, -6*x-2, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 2*x, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 2*x+1, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, -2*x-1, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 2*x, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, -2*x, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 2*x+1, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, -2*x-1}};
    g::matrix m3 = {{0, 0, 0, 0, 2, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0, 0, 0, 0, 0,}, {0, 0, 0, 0, 0, 0, 0, 1, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 1, 0},  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
    mat_vec basis = {m1,m2,m3};
    lie_algebra* L3_0_9a = new lie_algebra(basis);
    return L3_0_9a;
}

void test_spanning_subsequence(){
    g::symbol x("x");
    g::symbol y("y");

    g::matrix a = {{1,0}};
    g::matrix b = {{2,0}};
    g::matrix c = {{3,1}};
    mat_vec v = {a,b,c};

    mat_vec subsequence = lin_alg::spanning_subsequence(v);

    std::cout << "Printing matrices" << std::endl;
    utils::print_matrices(v);
    std::cout << "Printing spanning subsequence" << std::endl;
    utils::print_matrices(subsequence);

    a = {{1,0},{0,1}};
    b = {{2,0},{2,1}};
    c = {{3,1},{1,2}};
    g::matrix d = {{3,1},{3,1}};
    v = {a,b,c,d};

    subsequence = lin_alg::spanning_subsequence(v);

    std::cout << "Printing matrices" << std::endl;
    utils::print_matrices(v);
    std::cout << "Printing spanning subsequence" << std::endl;
    utils::print_matrices(subsequence);
}

void test_lie_alg_equals() {
    g::matrix a = {{1,0},{0,-1}};
    g::matrix b = {{0,1},{0,0}};
    g::matrix c = {{0,0},{1,0}};
    lie_algebra alg_1 = {{a,b,c}};
    lie_algebra alg_2 = {{b,c}};
    if(! alg_1.equals(&alg_2) ) {
        throw std::runtime_error("I goofed");
    }
    lie_algebra alg_3 = {{a}};
    if ( alg_1.equals(&alg_3)) {
        throw std::runtime_error("I goofed");
    }
}

void test_get_sl(int n) {
    lie_algebra* sl = lie_algebra::get_sl(n);
    if (sl->get_dim() != n*n-1) {
        throw std::runtime_error("I goofed");
    }
}

void test_bracket_algebra_sl_sl(int n) {
    lie_algebra* sl = lie_algebra::get_sl(n);
    if (! sl->equals(bracket_lie_algebras(sl, sl))) {
        throw std::runtime_error("I goofed");
    }
}

void test_get_normalizer_element() {
    g::matrix a = {{0,1,0}, {0,0,1}, {0,0,0}};
    g::matrix b = {{1,0,0}, {0,1,0}, {0,0,-2}};
    g::matrix c = {{1,0,0}, {0,1,0}, {0,0,1}};
    g::matrix d = {3,3};
    d(1,2) = 1;

    mat_vec L_basis = {d};
    lie_algebra* sl = lie_algebra::get_sl(3);
    lie_algebra* L = new lie_algebra(L_basis);

    mat_vec M = sl->get_basis();
    mat_vec normalizer = L->compute_normalizer_element(d, M);

    normalizer = L->compute_normalizer_element(a, M);
    std::cout << "Printing a" << std::endl;
    utils::print_matrix(a);
    std::cout << "Printing basis of L" << std::endl;
    utils::print_matrices(L->get_basis());
    std::cout << "Printing normalizer of a in M" << std::endl;
    utils::print_matrices(normalizer);

    mat_vec a_brackets;
    for (g::matrix v : normalizer){
        a_brackets.push_back(lin_alg::bracket(a,v));
    }
    std::cout << "Printing brackets" << std::endl;
    utils::print_matrices(a_brackets);

    std::cout << "Printing b" << std::endl;
    utils::print_matrix(b);
    std::cout << "Computing normalizer of b in M" << std::endl;
    normalizer = L->compute_normalizer_element(b, M);
    utils::print_matrices(normalizer);

    mat_vec b_brackets;
    for (g::matrix v : normalizer){
        b_brackets.push_back(lin_alg::bracket(b,v));
    }
    std::cout << "Printing brackets" << std::endl;
    utils::print_matrices(b_brackets);
}

void test_gaussian_elimination(){
    g::symbol x("x");
    g::symbol y("y");
    g::matrix a = {3,3};
    g::matrix b = {3,3};
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            b(i,j) = 1;
        }
    }
    g::matrix c = {{x, 0, 0}, {0, 0, y}};
    g::matrix d = {3,2};
    g::matrix e = {{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,-2,-1},{0,0,0,1,0,0,0,0},{-1,0,0,0,1,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,-1,0,0,0,0},{0,1,0,0,0,-1,0,0},{0,0,0,0,0,1,0,0}};

    g::matrix gaussian_a_computed = lin_alg::gaussian_elimination(a);
    g::matrix gaussian_b_computed = lin_alg::gaussian_elimination(b);
    g::matrix gaussian_c_computed = lin_alg::gaussian_elimination(c);
    g::matrix gaussian_d_computed = lin_alg::gaussian_elimination(d);
    g::matrix gaussian_e_computed = lin_alg::gaussian_elimination(e);

    g::matrix gaussian_a_actual = {{0,0,0},{0,0,0},{0,0,0}};
    g::matrix gaussian_b_actual = {{1,1,1},{0,0,0},{0,0,0}};
    g::matrix gaussian_c_actual = {{x,0,0},{0,0,y}};
    g::matrix gaussian_d_actual = {{0,0},{0,0},{0,0}};
    g::matrix gaussian_e_actual = {{-1,0,0,0,1,0,0,0},{0,1,0,0,0,0,0,0},{0,0,0,1,0,0,0,0},{0,0,0,0,0,1,0,0},
                                        {0,0,0,0,0,0,-2,-1},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}};
                                        

    if(!utils::matrix_eq(gaussian_a_computed, gaussian_a_actual)) {
        throw std::runtime_error("gaussian_elimination fails on a");
    }
    if(!utils::matrix_eq(gaussian_b_computed, gaussian_b_actual)) {
        throw std::runtime_error("gaussian_elimination fails on b");
    }
    if(!utils::matrix_eq(gaussian_c_computed, gaussian_c_actual)) {
        throw std::runtime_error("gaussian_elimination fails on c");
    }
    if(!utils::matrix_eq(gaussian_d_computed, gaussian_d_actual)) {
        throw std::runtime_error("gaussian_elimination fails on d");
    }
    if(!utils::matrix_eq(gaussian_e_computed, gaussian_e_actual)) {
        throw std::runtime_error("gaussian_elimination fails on e");
    }
}

void test_nullspace(){
    g::symbol x("x");
    g::symbol y("y");
    g::matrix a = {3,3};
    g::matrix b = {3,3};
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            b(i,j) = 1;
        }
    }
    g::matrix c = {{x, 0, 0}, {0, 0, y}};
    g::matrix d = {3,2};
    g::matrix e = {{x, 1, 0}, {0, 0, y}};
    g::matrix f = {{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,-2,-1},{0,0,0,1,0,0,0,0},{-1,0,0,0,1,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,-1,0,0,0,0},{0,1,0,0,0,-1,0,0},{0,0,0,0,0,1,0,0}};

    mat_vec null_a_computed = lin_alg::nullspace(a);
    mat_vec null_b_computed = lin_alg::nullspace(b);
    mat_vec null_c_computed = lin_alg::nullspace(c);
    mat_vec null_d_computed = lin_alg::nullspace(d);
    mat_vec null_e_computed = lin_alg::nullspace(e);
    mat_vec null_f_computed = lin_alg::nullspace(f);

    std::vector< g::exvector > null_a_actual = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    std::vector< g::exvector > null_b_actual = {{1, -1, 0}, {1, 0, -1}};
    std::vector< g::exvector > null_c_actual = {{0, 1, 0}};
    std::vector< g::exvector > null_d_actual = {{1, 0}, {0, 1}};
    std::vector< g::exvector > null_e_actual = {{1, -x, 0}};
    std::vector< g::exvector > null_f_actual = {{0,0,1,0,0,0,0,0},{1,0,0,0,1,0,0,0},{0,0,0,0,0,0,-1,2}};


    if(!lin_alg::equals(&null_a_computed, &null_a_actual)) {
        throw std::runtime_error("nullspace fails on a");
    }
    if(!lin_alg::equals(&null_b_computed, &null_b_actual)) {
        throw std::runtime_error("nullspace fails on b");
    }
    if(!lin_alg::equals(&null_c_computed, &null_c_actual)) {
        throw std::runtime_error("nullspace fails on c");
    }
    if(!lin_alg::equals(&null_d_computed, &null_d_actual)) {
        throw std::runtime_error("nullspace fails on d");
    }
    if(!lin_alg::equals(&null_e_computed, &null_e_actual)) {
        throw std::runtime_error("nullspace fails on e");
    }
    if(!lin_alg::equals(&null_f_computed, &null_f_actual)) {
        throw std::runtime_error("nullspace fails on f");
    }
}

void test_sl_ize(){
    lie_algebra* sl = lie_algebra::get_sl(3);
    mat_vec sl_basis = sl->get_basis();

    g::matrix a = {{-14,2,3},{4,5,6},{7,8,9}};
    g::matrix b = {3,3};
    g::matrix c = {{-2,0,0},{0,1,0},{0,0,1}};
    g::matrix d = {{1, 3, 0, -6}, {2, 1, 0, 4}, {0, 0, 1, 8}, {-2, 4, 0, -3}};

    g::exvector slize_a = lin_alg::sl_ize(a,3);
    g::exvector slize_b = lin_alg::sl_ize(b,3);
    g::exvector slize_c = lin_alg::sl_ize(c,3);
    g::exvector slize_d = lin_alg::sl_ize(d,4);

    g::exvector slize_a_real = {2,4,3,7,6,8,-5,-9};
    g::exvector slize_b_real = {0,0,0,0,0,0,0,0};
    g::exvector slize_c_real = {0,0,0,0,0,0,-1,-1};
    g::exvector slize_d_real = {3,2,0,0,-6,-2,0,0,4,4,8,0,-1,-1,3};

    g::matrix a_reconstruction = lin_alg::vector_to_matrix(slize_a, sl_basis);
    g::matrix b_reconstruction = lin_alg::vector_to_matrix(slize_b, sl_basis);
    g::matrix c_reconstruction = lin_alg::vector_to_matrix(slize_c, sl_basis);

    sl = lie_algebra::get_sl(4);
    sl_basis = sl->get_basis();
    g::matrix d_reconstruction = lin_alg::vector_to_matrix(slize_d, sl_basis);


    // test slize correct
    if (!utils::exvector_eq(slize_a, slize_a_real)) { 
        throw std::runtime_error("slize fails on a");
    }
    if (!utils::exvector_eq(slize_b, slize_b_real)) {
        throw std::runtime_error("slize fails on b");
    }
    if (!utils::exvector_eq(slize_c, slize_c_real)) {
        throw std::runtime_error("slize fails on c");
    }
    if (!utils::exvector_eq(slize_d, slize_d_real)) {
        throw std::runtime_error("slize fails on d");
    }

    // test slize is inverse of vector_to_matrix with sl_basis
    if (!utils::matrix_eq(a, a_reconstruction)) {
        throw std::runtime_error("slize fails to be an inverse to vector_to_matrix on a");
    }
    if (!utils::matrix_eq(b, b_reconstruction)) {
        throw std::runtime_error("slize fails to be an inverse to vector_to_matrix on b");
    }
    if (!utils::matrix_eq(c, c_reconstruction)) {
        throw std::runtime_error("slize fails to be an inverse to vector_to_matrix on c");
    }
    if (!utils::matrix_eq(d, d_reconstruction)) {
        throw std::runtime_error("slize fails to be an inverse to vector_to_matrix on d");
    }
}

void test_matricize(){
    g::symbol x("x");
    g::matrix a = {{1},{2},{3},{4},{5}};
    g::matrix b = {{x}, {1}, {0}, {2}};

    g::exvector c = {1,2,3,4,5,6,7,8,9};
    g::exvector d = {x,2,0,-x,1,2};
    
    g::matrix a1_matricize = lin_alg::matricize(a,1,5);
    g::matrix a2_matricize = lin_alg::matricize(a,5,1);
    g::matrix b1_matricize = lin_alg::matricize(b,2,2);
    g::matrix b2_matricize = lin_alg::matricize(b,2,2,true);
    
    g::matrix c1_matricize = lin_alg::matricize(c,3,3);
    g::matrix c2_matricize = lin_alg::matricize(c,3,3,true);
    g::matrix d1_matricize = lin_alg::matricize(d,2,3);
    g::matrix d2_matricize = lin_alg::matricize(d,2,3,true);

    g::matrix a1_matricize_true = {{1,2,3,4,5}};
    g::matrix a2_matricize_true = {{1},{2},{3},{4},{5}};
    g::matrix b1_matricize_true = {{x,1},{0,2}};
    g::matrix b2_matricize_true = {{x,0},{1,2}};
    
    g::matrix c1_matricize_true = {{1,2,3},{4,5,6},{7,8,9}};
    g::matrix c2_matricize_true = {{1,4,7},{2,5,8},{3,6,9}};
    g::matrix d1_matricize_true = {{x,2,0},{-x,1,2}};
    g::matrix d2_matricize_true = {{x,0,1},{2,-x,2}};
    
    if(!utils::matrix_eq(a1_matricize, a1_matricize_true)) {
        throw std::runtime_error("matricize(a,1,5) failed");
    }
    if(!utils::matrix_eq(a2_matricize, a2_matricize_true)) {
        throw std::runtime_error("matricize(a,5,1) failed");
    }
    if(!utils::matrix_eq(b1_matricize, b1_matricize_true)) {
        throw std::runtime_error("matricize(b,2,2) by rows failed");
    }
    if(!utils::matrix_eq(b2_matricize,b2_matricize_true)) {
        throw std::runtime_error("matricize(b,2,2) by cols failed");
    }
    if(!utils::matrix_eq(c1_matricize,c1_matricize_true)) {
        throw std::runtime_error("matricize(c,3,3) by rows failed");
    }
    if(!utils::matrix_eq(c2_matricize,c2_matricize_true)) {
        throw std::runtime_error("matricize(c,3,3) by cols failed");
    }
    if(!utils::matrix_eq(d1_matricize,d1_matricize_true)) {
        throw std::runtime_error("matricize(d,2,3) by rows failed");
    }
    if(!utils::matrix_eq(d2_matricize,d2_matricize_true)) {
        throw std::runtime_error("matricize(d,2,3) by cols failed");
    }
}

void test_normalizer() {
    g::symbol x("x");
    lie_algebra* L5_1 = get_L5_1();
    lie_algebra* L5_2 = get_L5_2();
    lie_algebra* L5_4 = get_L5_4(x);
    lie_algebra* L5_5 = get_L5_5(x);
    lie_algebra* L3_0_9a = get_L3_0_9a(x);
    // lie_algebra* test_1 = get_test_alg_1();

    lie_algebra* normalizer_L5_1 = L5_1->compute_normalizer();
    lie_algebra* normalizer_L5_2 = L5_2->compute_normalizer();
    lie_algebra* normalizer_L5_4 = L5_4->compute_normalizer();
    lie_algebra* normalizer_L5_5 = L5_5->compute_normalizer();
    lie_algebra* normalizer_L3_0_9a = L3_0_9a->compute_normalizer();
    // lie_algebra* normalizer_test_1 = test_1->compute_normalizer();
    
    mat_vec L5_1_normalizer_basis = {
        {{0,  1,  0,  0}, {0,  0,  0,  0}, {0,  0,  0,  0}, {0,  0,  0,  0}}, 
        {{0,  0,  1,  0}, {0,  0,  0,  0}, {0,  0,  0,  0}, {0,  0,  0,  0}},
        {{0,  0,  0,  0}, {0,  0,  1,  0}, {0,  0,  0,  0}, {0,  0,  0,  0}},
        {{1,  0,  0,  0}, {0,  1,  0,  0}, {0,  0,  1,  0}, {0,  0,  0,  -3}}
    };

    mat_vec L5_2_normalizer_basis = {
        {{1,  0,  0,  0}, {0,  0,  0,  0}, {0,  0,  -1,  0}, {0,  0,  0,  0}},
        {{0,  1,  0,  0}, {0,  0,  0,  0}, {0,  0,  0,  0}, {0,  0,  0,  0}},
        {{0,  0,  0,  1}, {0,  0,  0,  0}, {0,  0,  0,  0}, {0,  0,  0,  0}},
        {{0,  0,  0,  0}, {0,  0,  0,  1}, {0,  0,  0,  0}, {0,  0,  0,  0}},
        {{0,  0,  0,  0}, {0,  1,  0,  0}, {0,  0,  -2,  0}, {0,  0,  0,  1}}
    };

    mat_vec L5_4_normalizer_basis = {
        {{0, 0, 1, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}},                               
        {{0, (2*x + 1)/(x+1),  0,  0}, {0,  0,  1,  0}, {0,  0,  0,  0}, {0,  0,  0,  0}},      
        {{-3,  0,  0,  0}, {0, -1,  0,  0}, {0,  0,  1,  0}, {0,  0,  0,  3}},                  
        {{0,  0,  0,  1}, {0,  0,  0,  0}, {0,  0,  0,  0}, {0,  0,  0,  0}},                   
        {{0,  0,  0,  0}, {0,  0,  0,  1}, {0,  0,  0,  0}, {0,  0,  0,  0}},                   
        {{0,  -x/(x+1), 0,  0}, {0,  0,  0,  0}, {0,  0,  0,  1}, {0,  0,  0, 0}}               
    };

    mat_vec L5_5_normalizer_basis = {
        {{0 , 1 , 0 , 0}, {0 , 0 , 0 , 0}, {0,  0,  0,  0}, {0,  0,  0,  0}},
        {{0 , 0 , 0 , 0}, {0 , 1 , 0 , 0}, {0,  0,  -1,  0}, {0,  0,  0,  0}},
        {{0 , 0 , 1 , 0}, {0 , 0 , 0 , 0}, {0,  0,  0,  0}, {0,  0,  0,  0}},
        {{0 , 0 , 0 , 1}, {0 , 0 , 0 , 0}, {0,  0,  0,  0}, {0,  0,  0,  0}},
        {{0 , 0 , 0 , 0}, {0 , 0 , 0 , 1}, {0,  0,  0,  0}, {0,  0,  0,  0}}, 
        {{0 , 0 , 0 , 0}, {0 , 0 , 0 , 0}, {0,  0,  0,  1}, {0,  0,  0,  0}},
        {{-1,  0,  0, 0},{ 0,  0,  0,  0}, { 0,  0,  0,  0}, {0,  0,  0,  1}}
    };

    lie_algebra normalizer_actual_L5_1 = {L5_1_normalizer_basis, true};
    lie_algebra normalizer_actual_L5_2 = {L5_2_normalizer_basis, true};
    lie_algebra normalizer_actual_L5_4 = {L5_4_normalizer_basis, true};
    lie_algebra normalizer_actual_L5_5 = {L5_5_normalizer_basis, true};

    if(!normalizer_actual_L5_1.equals(normalizer_L5_1)) {
        throw std::runtime_error("I normalizer goofed with L5_1");
    }
    if(!normalizer_actual_L5_2.equals(normalizer_L5_2)) {
        throw std::runtime_error("I normalizer goofed with L5_2");
    }
    if(!normalizer_actual_L5_4.equals(normalizer_L5_4)) {
        throw std::runtime_error("I normalizer goofed with L5_4");
    }
    if(!normalizer_actual_L5_5.equals(normalizer_L5_5)) {
        throw std::runtime_error("I normalizer goofed with L5_5");
    }
}

void test_centralizer() {
    g::symbol x("x");
    lie_algebra* L5_1 = get_L5_1();
    lie_algebra* L5_2 = get_L5_2();
    lie_algebra* L5_4 = get_L5_4(x);
    lie_algebra* L5_5 = get_L5_5(x);
    
    lie_algebra* sl = lie_algebra::get_sl(4);

    lie_algebra* centralizer_L5_1 = L5_1->compute_centralizer();
    lie_algebra* centralizer_L5_2 = L5_2->compute_centralizer();
    lie_algebra* centralizer_L5_4 = L5_4->compute_centralizer();
    lie_algebra* centralizer_L5_5 = L5_5->compute_centralizer();
    
    mat_vec L5_1_cent_basis = {
        {{1,  0,  0,  0}, {0,  1,  0,  0}, {0,  0,  1,  0}, {0,  0,  0,  -3}}, 
        {{0,  0,  1,  0}, {0,  0,  0,  0}, {0,  0,  0,  0}, {0,  0,  0,  0}}
    };

    mat_vec L5_2_cent_basis = {
        {{1,  0,  0,  0}, {0,  1,  0,  0}, {0,  0,  -3,  0}, {0,  0,  0,  1}},
        {{0,  0,  0,  1}, {0,  0,  0,  0}, {0,  0,  0,  0}, {0,  0,  0,  0}}
    };
    
    mat_vec L5_4_cent_basis = {
        {{0,  0,  1,  0}, {0,  0,  0,  1}, {0,  0,  0,  0}, {0,  0,  0,  0}},
        {{0,  0,  0,  1}, {0,  0,  0,  0}, {0,  0,  0,  0}, {0,  0,  0,  0}}
    };

    mat_vec L5_5_cent_basis = {
        {{0,  x/(x+1),  0,  0}, {0,  0,  0,  0}, {0,  0,  0,  1}, {0,  0,  0,  0}},
        {{0,  0,  0,  1}, {0,  0,  0,  0}, {0,  0,  0,  0}, {0,  0,  0,  0}},
        {{0,  0,  1,  0}, {0,  0,  0,  1}, {0,  0,  0,  0}, {0,  0,  0,  0}},
    };

    lie_algebra centralizer_actual_L5_1 = {L5_1_cent_basis, true};
    lie_algebra centralizer_actual_L5_2 = {L5_2_cent_basis, true};
    lie_algebra centralizer_actual_L5_4 = {L5_4_cent_basis, true};
    lie_algebra centralizer_actual_L5_5 = {L5_5_cent_basis, true};
    
    if(! centralizer_actual_L5_1.equals(centralizer_L5_1)) {
        throw std::runtime_error("I centralizer goofed with L5_1");
    }
    if(! centralizer_actual_L5_2.equals(centralizer_L5_2)) {
        throw std::runtime_error("I centralizer goofed with L5_2");
    }

    if(! centralizer_actual_L5_4.equals(centralizer_L5_4)) {
        throw std::runtime_error("I centralizer goofed with L5_4");
    }
    if(! centralizer_actual_L5_5.equals(centralizer_L5_5)) {
        throw std::runtime_error("I centralizer goofed with L5_5");
    }
}

void test_minrank() {
    g::symbol x("x");
    lie_algebra* L5_1 = get_L5_1();
    lie_algebra* L5_2 = get_L5_2();
    lie_algebra* L5_4 = get_L5_4(x);
    lie_algebra* L5_5 = get_L5_5(x);

    if (L5_1->min_rank() != 1) {
        throw std::runtime_error("The minrank of L5_1 is wrong, it should be 1");
    }
    if (L5_2->min_rank() != 1) {
        throw std::runtime_error("The minrank of L5_2 is wrong, it should be 1");
    }
    if (L5_4->min_rank() != 1) {
        throw std::runtime_error("The minrank of L5_3 is wrong, it should be 1");
    }
    if (L5_5->min_rank() != 1) {
        throw std::runtime_error("The minrank of L5_4 is wrong, it should be 1");
    }
}

void test_get_leading_term() {
    g::symbol x("x");
    g::ex p = g::pow(x,3) + g::pow(x,2) - 2*x + 1;
    g::ex q = g::pow(x,2) - 1;
    g::ex h = -3 * g::pow(x,4) + 2*x;
    g::ex f = 2 * g::pow(x,3) + 3 * x - 1;
    g::symbol y("y");
    g::ex r = g::pow(y,2) + 2 * y + 1 + x*y + y*y*y*3 - x*y*y*2;
    g::ex s = y*y;

    std::vector< g::symbol > x_vec = {x};
    std::vector< g::symbol > y_vec = {x,y};

    if (!get_leading_term(p,x_vec).is_equal(g::pow(x,3))) {
        throw std::runtime_error("The leading term of p is wrong");
    }
    if (!get_leading_term(q,x_vec).is_equal(g::pow(x,2))) {
        throw std::runtime_error("The leading term of q is wrong");
    }
    if (!get_leading_term(h,x_vec).is_equal(-3 * g::pow(x,4))) {
        throw std::runtime_error("The leading term of h is wrong");
    }
    if (!get_leading_term(f,x_vec).is_equal(2 * g::pow(x,3))) {
        throw std::runtime_error("The leading term of f is wrong");
    }
    if (!get_leading_term(r,y_vec).is_equal(  -2*x*g::pow(y,2))) {
        throw std::runtime_error("The leading term of r is wrong");
    }
    if (!get_leading_term(s,y_vec).is_equal(  y*y)) {
        throw std::runtime_error("The leading term of s is wrong");
    }
}

void test_polynomial_divide_one_var() {
    g::symbol x("x");
    g::ex p = g::pow(x,3) + g::pow(x,2) - 2*x + 1;
    g::ex q = g::pow(x,2) - 1;
    g::exvector q_vec = {q};
    std::vector< g::symbol > x_vec = {x};
    g::exvector out = polynomial_divide(p,q_vec,x_vec);
    
    if(!out[0].is_equal(x+1)) {
        throw std::runtime_error("The univariate polynomial didn't divide properly");
    }

    if(!out[1].is_equal(2-x)) {
        throw std::runtime_error("The univariate polynomials didn't divide properly");
    }

    p = 1;
    q = 2;
    out = polynomial_divide(p,{q},{x});
    
    if(!out[0].is_equal(1.0/2.0)) {
        throw std::runtime_error("The univariate polynomial didn't divide properly");
    }

    if(!out[1].is_equal(0)) {
        throw std::runtime_error("The univariate polynomials didn't divide properly");
    }
}

void test_polynomial_divide_two_var() {
    g::symbol x("x"); 
    g::symbol y("y");
    
    // g::ex p = g::pow(x,3)*g::pow(y,2)-2*g::pow(x,3)+g::pow(x,2)*y+x*y+1;
    g::ex p = x*x + x - y*y + y;
    g::exvector q_vec = {x*y + 1, x + y};
    std::vector< g::symbol > var_vec = {x,y};
    g::exvector out = polynomial_divide(p,q_vec,var_vec);

    if(!out[0].is_equal(-1)) {
        throw std::runtime_error("The duovariate polynomial didn't divide properly");
    }
    if(!out[1].is_equal(x + 1)) {
        throw std::runtime_error("The duovariate polynomials didn't divide properly");
    }
    if(!out[2].is_equal(-y*y + 1)) {
        throw std::runtime_error("The duovariate polynomials didn't divide properly");
    }

    q_vec = {x + y, x*y + 1};
    var_vec = {x,y};
    out = polynomial_divide(p,q_vec,var_vec);

    if(!out[0].is_equal(x-y+1)) {
        throw std::runtime_error("The duovariate polynomial didn't divide properly");
    }
    if(!out[1].is_equal(0)) {
        throw std::runtime_error("The duovariate polynomials didn't divide properly");
    }
    if(!out[2].is_equal(0)) {
        throw std::runtime_error("The duovariate polynomials didn't divide properly");
    }
}

void test_polynomial_divide_three_var() {
    g::symbol x("x"); 
    g::symbol y("y"); 
    g::symbol z("z");
    
    g::ex p = z*y*x*x - 2*x*y*y + x*z*y - y*z + 3*x - z + 2;
    g::exvector q_vec = {x,y,z};
    std::vector< g::symbol > var_vec = {x,y,z};
    g::exvector out = polynomial_divide(p,q_vec,var_vec);

    if(!out[0].is_equal(x*z*y-2*y*y+z*y+3)) {
        throw std::runtime_error("The trivariate polynomial didn't divide properly");
    }
    if(!out[1].is_equal(-z)) {
        throw std::runtime_error("The trivariate polynomials didn't divide properly");
    }
    if(!out[2].is_equal(-1)) {
        throw std::runtime_error("The trivariate polynomials didn't divide properly");
    }
    if(!out[3].is_equal(2)) {
        throw std::runtime_error("The trivariate polynomials didn't divide properly");
    }    
}


void test_polynomial_divide_dividend_lower_deg() {
    g::symbol x("x"); 
    g::symbol y("y"); 
    
    g::ex p = x - y + 1;
    g::exvector q_vec = {x*x - y, y*y};

    g::exvector out = polynomial_divide(p,q_vec,{x,y});

    if(!out[0].is_equal(0)) {
        throw std::runtime_error("The polynomial didn't divide properly by the lower degree polynomial (1.1)");
    }
    if(!out[1].is_equal(0)) {
        throw std::runtime_error("The polynomial didn't divide properly by the lower degree polynomial (1.2)");
    }
    if(!out[2].is_equal(p)) {
        throw std::runtime_error("The polynomial didn't divide properly by the lower degree polynomial (1.3)");
    }

    p = 1;
    q_vec = {x*x - y, y*y, 3};
    
    out = polynomial_divide(p,q_vec,{x,y});

    if(!out[0].is_equal(0)) {
        throw std::runtime_error("The polynomial didn't divide properly by the lower degree polynomial (2.1)");
    }
    if(!out[1].is_equal(0)) {
        throw std::runtime_error("The polynomial didn't divide properly by the lower degree polynomial (2.2)");
    }
    g::ex l = 1;
    g::ex lll = 3;
    g::ex llll = l/lll;
    if(!out[2].is_equal(llll)) {
        throw std::runtime_error("The polynomial didn't divide properly by the lower degree polynomial (2.3)");
    }
    if(!out[3].is_equal(0)) {
        throw std::runtime_error("The polynomial didn't divide properly by the lower degree polynomial (2.4)");
    }
}

void test_polynomial_divide() {
    test_polynomial_divide_one_var();
    test_polynomial_divide_two_var();
    test_polynomial_divide_three_var();
    test_polynomial_divide_dividend_lower_deg();
}

void test_S() {
    g::symbol x("x"); 
    g::symbol y("y"); 
    g::symbol z("z");

    g::ex p = y*x*x - 2*x*y*y + x*z*y - y*z + 3*x - z + 2;
    g::ex q = z*x + x - y*y + y - 4;

    g::ex out = S_function_from_Dummit_foot(p,q,{x,y,z});
    g::ex answer = -x*x*y + x*y*y*y - 2*x*y*y*z - x*y*y + x*y*z*z + 4*x*y + 3*x*z - y*z*z - z*z + 2*z;
}

void test_get_groebner_basis(){
    g::symbol x("x"); 
    g::symbol y("y");

    g::ex f1 = x*x*x*y - x*y*y + 1;
    g::ex f2 = x*x*y*y - y*y*y - 1;

    std::vector< g::ex > out = get_groebner_basis({f1,f2},{x,y});

    if(size(out) != 4) {
        throw std::runtime_error("Groebner basis size is wrong");
    }
    if(!out[0].is_equal(f1)) {
        throw std::runtime_error("Groebner basis is wrong");
    }
    if(!out[1].is_equal(f2)) {
        throw std::runtime_error("Groebner basis is wrong");
    }
    if(!out[2].is_equal(x+y)) {
        throw std::runtime_error("Groebner basis is wrong");
    }
    if(!out[3].is_equal(y*y*y*y - y*y*y - 1)) {
        throw std::runtime_error("Groebner basis is wrong");
    }   
}

void test_grob(){
    g::symbol x("x");
    g::symbol y("y");
    g::symbol t_0("t_0");
    g::symbol t_1("t_1");
    g::symbol t_2("t_2");

    g::exvector I1 = {x*x + y - x*y + 2, x*y + y*y - 2, x*x*x - y};                             // Contains 1
    g::exvector I2 = {x*x + y - x*y + 2, x*x*x*x - x*y*y + 2*x*x + 2*y*y - 4*x + 4*y - 4};      // Does not contain 1
    g::exvector I3 = {1};                                                                       // Contains 1
    g::exvector I4 = {x * y + 1, x - 1, y - 1, x + y};                                          // Contains 1
    g::exvector I5 = {x};                                                                       // Does not contain 1
    g::exvector I6 = {x*x + 2*x - 1, x*x*x + 2*x + 1, x*x*x - x*x};                             // Contains 1
    g::exvector I7 = {0, t_0, x, t_2, 0, 0, 0, 1+x, 0, 0, 0, t_0, 0, 0, 0, 0};                  // Contains 1

    if(!grob(I1, {x,y})) {
        throw std::runtime_error("I1 is actually grob.");
    }
    
    if(grob(I2, {x,y})) {
        throw std::runtime_error("I2 is actually not grob.");
    }
    
    if(!grob(I3, {x,y})) {
        throw std::runtime_error("I3 is actually grob.");
    }
    
    if(!grob(I4, {x,y})) {
        throw std::runtime_error("I4 is actually grob.");
    }

    if(grob(I5, {x,y})) {
        throw std::runtime_error("I5 is actually not grob.");
    }

    if(!grob(I6, {x,y})) {
        throw std::runtime_error("I6 is actually grob.");
    }

    if(!grob(I7, {t_0,t_1,t_2,x})) {
        throw std::runtime_error("I7 is actually grob.");
    }
}

void test_get_minors() {
    g::symbol x("x");
    g::matrix m = {{x + 1, 2, 0},{x, 1, 1},{0, x - 2, x*x}};
    g::matrix m11 = {{1,1}, {x-2, x*x}};
    g::matrix m12 = {{x,1}, {0, x*x}};
    g::matrix m13 = {{x,1}, {0, x-2}};
    g::matrix m21 = {{2,0}, {x-2,x*x}};
    g::matrix m22 = {{x+1,0}, {0,x*x}};
    g::matrix m23 = {{x+1,2}, {0,x-2}};
    g::matrix m31 = {{2,0}, {1,1}};
    g::matrix m32 = {{x+1,0}, {x,1}};
    g::matrix m33 = {{x+1,2}, {x,1}};
    std::vector< g::matrix > m_minors_matrices = {m33, m32, m31, m23, m22, m21, m13, m12, m11};
    g::exvector m_minors = {};
    for (g::matrix mat : m_minors_matrices) {
        m_minors.push_back(mat.determinant());
    }
    g::exvector m_minors_out = lin_alg::get_minors(m, 2);
    g::exvector m_minors_out1 = lin_alg::get_minors(m, 1);
    g::exvector m_minors_out2 = lin_alg::get_minors(m, 3);
        
    if (size(m_minors_out) != 9) {
        throw std::runtime_error("wrong number of minors");
    }
    for (int i = 0; i < 9; i ++){
        if (m_minors_out[i].is_equal(m_minors[i]) == false) {
            throw std::runtime_error("wrong minor at index " + std::to_string(i));
        }
    }
}

int main() {
    test_lie_alg_equals();
    test_get_sl(6);
    test_bracket_algebra_sl_sl(6);
    test_gaussian_elimination();
    test_nullspace();
    test_sl_ize();
    test_matricize();
    test_normalizer();
    test_centralizer();
    test_minrank();
    test_get_leading_term();
    test_polynomial_divide();
    test_get_groebner_basis();
    test_S();
    test_grob();
    test_get_minors();

    // Do not run, only prints, no auto-tests for now
    // test_spanning_subsequence();
    // test_get_normalizer_element();

    return 0;
}
