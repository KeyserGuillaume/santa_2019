#include "day_cost_lower_bound.h"


double get_day_cost(const unsigned int &i, const unsigned int &ip1) {
    unsigned int gap = (i > ip1) ? i - ip1 : ip1 - i;
    return (i - MIN_NB_PEOPLE_PER_DAY)/400.*pow(i, 0.5 + gap/50.);
}

double get_day_cost_DP_lb(const std::vector<std::vector<unsigned int>> &possible_quantities) {
    std::vector<std::vector<double>> value_functions(NB_DAYS, std::vector<double>(0));
    std::vector<std::vector<unsigned int>> argmin_value_functions(NB_DAYS, std::vector<unsigned int>(0));
    for (unsigned int m = 0; m < possible_quantities[NB_DAYS - 1].size(); m++){
        value_functions[NB_DAYS - 1].push_back(get_day_cost(possible_quantities[NB_DAYS - 1][m], possible_quantities[NB_DAYS - 1][m]));
    }
    for (int i = NB_DAYS - 2; i >= 0; i--){
        for (unsigned int n = 0; n < possible_quantities[i].size(); n++) {
            double mini = 100000, current;
            unsigned int argmin;
            for (unsigned int m = 0; m < possible_quantities[i + 1].size(); m++) {
                current = value_functions[i + 1][m] + get_day_cost(possible_quantities[i][n], possible_quantities[i + 1][m]);
                if (current < mini) {
                    mini = current;
                    argmin = m;
                }
            }
            value_functions[i].push_back(mini);
            argmin_value_functions[i + 1].push_back(argmin);
        }
    }
    double mini = 100000;
    unsigned int argmin;
    for (unsigned int m = 0; m < possible_quantities[0].size(); m++) {
        if (value_functions[0][m] < mini){
            mini = value_functions[0][m];
            argmin = m;
        }
    }
//    std::vector<unsigned int> solution(NB_DAYS, 0);
//    solution[0] = possible_quantities[0][argmin];
//    for (unsigned int i = 1; i < NB_DAYS; i++){
//        argmin = argmin_value_functions[i][argmin];
//        solution[i] = possible_quantities[i][argmin];
//    }
//    print_nicely(solution);
    return mini;
}

double get_binary_day_cost_lb(const unsigned int &a_i, const unsigned int &b_i, const unsigned int &a_ip1,
                             const unsigned int &b_ip1) {
    unsigned int gap;
    if (b_i <= a_ip1)
        gap = a_ip1 - b_i;
    else if (b_ip1 <= a_i)
        gap = a_i - b_ip1;
    else
        gap = 0;
    return (a_i - MIN_NB_PEOPLE_PER_DAY)/400.*pow(a_i, 0.5 + gap/50.);
}

double get_ternary_day_cost_lb(const unsigned int &a_im1, const unsigned int &b_im1, const unsigned int &a_i,
                              const unsigned int &b_i, const unsigned int &a_ip1, const unsigned int &b_ip1) {
    double mini = 100000;
    for (unsigned int m = a_i; m <= b_i && mini; m++)
        mini = std::min(double(mini), 0.5*get_binary_day_cost_lb(a_im1, b_im1, m, m) + 0.5*get_binary_day_cost_lb(m, m, a_ip1, b_ip1));
    return mini;
}

double get_day_cost_DP_ultimate_lb(const std::vector<std::vector<unsigned int>> &possible_quantities, const unsigned int& sum_ni) {
    std::vector<unsigned int> Nmax(NB_DAYS, MIN_NB_PEOPLE_PER_DAY);
    std::vector<unsigned int> Nmin(NB_DAYS, MAX_NB_PEOPLE_PER_DAY);
    for (unsigned int i = 0; i < NB_DAYS; i++){
        for (unsigned int j = 0; j < possible_quantities[i].size(); j++){
            Nmax[i] = std::max(Nmax[i], possible_quantities[i][j]);
            Nmin[i] = std::min(Nmin[i], possible_quantities[i][j]);
        }
    }

    std::vector<unsigned int> Mmin(NB_DAYS);
    std::vector<unsigned int> Mmax(NB_DAYS);
    for (unsigned int i = 0; i < NB_DAYS; i++){
        unsigned int s1 = 0;
        unsigned int s2 = 0;
        unsigned int s3 = 0;
        unsigned int s4 = 0;
        for (unsigned int j = 0; j < i; j++) {
            s1 += Nmin[j];
            s2 += Nmax[j];
        }
        for (unsigned int j = i; j < NB_DAYS; j++) {
            s3 += Nmin[j];
            s4 += Nmax[j];
        }
        Mmin[i] = std::max(int(s3), int(sum_ni - s2));
        Mmax[i] = std::min(int(s4), int(sum_ni - s1));
    }

    // value_function[t][n][M - Mmin],
    // 0       <= t <= 99,  
    // corresponds to nth element in possible_quantities[t],
    // Mmin    <= M <= Mmax.
    std::vector<std::vector<std::vector<double>>> value_functions(NB_DAYS, std::vector<std::vector<double>>(0));
    std::vector<std::vector<std::vector<unsigned int>>> argmin_value_functions(NB_DAYS, std::vector<std::vector<unsigned int>>(0));
    for (unsigned int m = 0; m < possible_quantities[NB_DAYS - 1].size(); m++){
        unsigned int q = possible_quantities[NB_DAYS - 1][m];
        double cost = get_day_cost(q, q);
        value_functions[NB_DAYS - 1].push_back(std::vector<double>(Mmax[NB_DAYS - 1] - Mmin[NB_DAYS - 1] + 1, 100000));
        value_functions[NB_DAYS - 1][m][q - Mmin[NB_DAYS - 1]] = cost;
    }
    for (int i = NB_DAYS - 2; i >= 0; i--){
        for (unsigned int n = 0; n < possible_quantities[i].size(); n++) {
            value_functions[i].push_back(std::vector<double>(Mmax[i] - Mmin[i] + 1, 100000));
            argmin_value_functions[i + 1].push_back(std::vector<unsigned int>(Mmax[i] - Mmin[i] + 1));
            for (unsigned int M = Mmin[i]; M <= Mmax[i]; M++) {
                unsigned int q = possible_quantities[i][n];
                if (M - q < Mmin[i + 1] || M - q > Mmax[i + 1])
                    continue;
                double mini = 100000, current;
                unsigned int argmin;
                for (unsigned int m = 0; m < possible_quantities[i + 1].size(); m++) {
                    current = value_functions[i + 1][m][M - q - Mmin[i + 1]] +
                              get_day_cost(q, possible_quantities[i + 1][m]);
                    if (current < mini) {
                        mini = current;
                        argmin = m;
                    }
                }
                value_functions[i][n][M - Mmin[i]] = mini;
                argmin_value_functions[i + 1][n][M - Mmin[i]] = argmin;
            }
        }
    }

    double mini = 100000;
    unsigned int argmin;
    for (unsigned int m = 0; m < possible_quantities[0].size(); m++) {
        if (value_functions[0][m][0] < mini){
            mini = value_functions[0][m][0];
            argmin = m;
        }
    }

    std::vector<unsigned int> solution(NB_DAYS, 0);
    solution[0] = possible_quantities[0][argmin];
    unsigned int M = sum_ni;
    for (unsigned int i = 1; i < NB_DAYS; i++){
        unsigned int prev_argmin = argmin;
        argmin = argmin_value_functions[i][argmin][M - Mmin[i - 1]];
        M -= possible_quantities[i - 1][prev_argmin];
        solution[i] = possible_quantities[i][argmin];
    }
    print_nicely(solution);
    unsigned int s = 0; float c = 0;
    for (unsigned int i = 0; i < NB_DAYS; i++){
        s += solution[i];
    }
    for (unsigned int i = 0; i < NB_DAYS - 1; i++){
        c += get_day_cost(solution[i], solution[i + 1]);
    }
    c += get_day_cost(solution[NB_DAYS - 1], solution[NB_DAYS - 1]);

    std::cout << "the ultimate DP gives a day costs lower bound of "  << mini << " " << s << " " << c << std::endl;
    
    return 0;
}

