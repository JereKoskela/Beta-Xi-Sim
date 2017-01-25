/*
** Copyright (C) 2017 Jere Koskela <koskela@math.tu-berlin.de>
**  
** This file is part of Beta-Xi-Sim.
** 
** Beta-Xi-Sim is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
** 
** Beta-Xi-Sim is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with Beta-Xi-Sim.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <cassert>
#include <cstdlib>
#include <gsl/gsl_matrix.h>
#include <iostream>
#include <unistd.h>
#include "ancestry.hh"

void print_ancestry(const Ancestry &ancestry) {
    for (int i = 0; i < ancestry.site_no; i++) {
        for (int j = 0; j < 2 * ancestry.sample_size - 1; j++) {
            std::cout << gsl_matrix_get(ancestry.anc, i, j) << " ";
        }
        std::cout << std::endl;
    }
    return;
}

void print_types(const Ancestry &ancestry) {
    for (int i = 0; i < ancestry.site_no; i++) {
        for (int j = 0; j < 2 * ancestry.sample_size - 1; j++) {
            std::cout << gsl_matrix_get(ancestry.type, i, j) << " ";
        }
        std::cout << std::endl;
    }
    return;
}

void print_coalescence_times(const Ancestry &ancestry) {
    for (int i = 0; i < ancestry.site_no; i++) {
        for (int j = 0; j < 2 * ancestry.sample_size - 1; j++) {
            std::cout << gsl_matrix_get(ancestry.t, i, j) << " ";
        }
        std::cout << std::endl;
    }
    return;
}

void print_singletons(const Ancestry &ancestry) {
    for (int i = 0; i < ancestry.site_no; i++) {
        for (int j = 0; j < ancestry.sample_size; j++) {
            std::cout << gsl_matrix_get(ancestry.type, i, j) << " ";
        }
        std::cout << std::endl;
    }
    return;
}

void print_pairwise_differences(const Ancestry &ancestry) {
    int ancestor_1, ancestor_2, diff_1 = 0, diff_2 = 0;
    for (int i = 0; i < ancestry.sample_size - 1; i++) {
        for (int j = i + 1; j < ancestry.sample_size; j++) {
            for (int k = 0; k < ancestry.site_no - 1; k++) {
                ancestor_1 = gsl_matrix_get(ancestry.anc, k, i);
                ancestor_2 = gsl_matrix_get(ancestry.anc, k, j);
                diff_1 = gsl_matrix_get(ancestry.type, k, i) 
                            + gsl_matrix_get(ancestry.type, k, j);
                while (ancestor_1 != ancestor_2) {
                    if (ancestor_1 < ancestor_2) {
                        diff_1 += gsl_matrix_get(ancestry.type, k, ancestor_1);
                        ancestor_1 = gsl_matrix_get(ancestry.anc, k, ancestor_1);
                    } else {
                        diff_1 += gsl_matrix_get(ancestry.type, k, ancestor_2);
                        ancestor_2 = gsl_matrix_get(ancestry.anc, k, ancestor_2);
                    }
                }
                for (int q = k + 1; q < ancestry.site_no; q++) {
                    ancestor_1 = gsl_matrix_get(ancestry.anc, q, i);
                    ancestor_2 = gsl_matrix_get(ancestry.anc, q, j);
                    diff_2 = gsl_matrix_get(ancestry.type, q, i) 
                                + gsl_matrix_get(ancestry.type, q, j);
                    while (ancestor_1 != ancestor_2) {
                        if (ancestor_1 < ancestor_2) {
                            diff_2 += gsl_matrix_get(ancestry.type, q, ancestor_1);
                            ancestor_1 = gsl_matrix_get(ancestry.anc, q, ancestor_1);
                        } else {
                            diff_2 += gsl_matrix_get(ancestry.type, q, ancestor_2);
                            ancestor_2 = gsl_matrix_get(ancestry.anc, q, ancestor_2);
                        }
                    }
                    std::cout << diff_1 << " " << diff_2 << std::endl;
                }
            }
        }
    }
    return;
}

void print_nonoverlapping_pairwise_differences(const Ancestry &ancestry) {
    int ancestor_1, ancestor_2, diff_1 = 0, diff_2 = 0;
    for (int i = 0; i < ancestry.sample_size - 1; i += 2) {
        for (int k = 0; k < ancestry.site_no - 1; k += 2) {
            ancestor_1 = gsl_matrix_get(ancestry.anc, k, i);
            ancestor_2 = gsl_matrix_get(ancestry.anc, k, i + 1);
            diff_1 = gsl_matrix_get(ancestry.type, k, i) 
                        + gsl_matrix_get(ancestry.type, k, i + 1);
            while (ancestor_1 != ancestor_2) {
                if (ancestor_1 < ancestor_2) {
                    diff_1 += gsl_matrix_get(ancestry.type, k, ancestor_1);
                    ancestor_1 = gsl_matrix_get(ancestry.anc, k, ancestor_1);
                } else {
                    diff_1 += gsl_matrix_get(ancestry.type, k, ancestor_2);
                    ancestor_2 = gsl_matrix_get(ancestry.anc, k, ancestor_2);
                }
            }
            ancestor_1 = gsl_matrix_get(ancestry.anc, k + 1, i);
            ancestor_2 = gsl_matrix_get(ancestry.anc, k + 1, i + 1);
            diff_2 = gsl_matrix_get(ancestry.type, k + 1, i) 
                    + gsl_matrix_get(ancestry.type, k + 1, i + 1);
            while (ancestor_1 != ancestor_2) {
                if (ancestor_1 < ancestor_2) {
                    diff_2 += gsl_matrix_get(ancestry.type, k + 1, ancestor_1);
                    ancestor_1 = gsl_matrix_get(ancestry.anc, k + 1, ancestor_1);
                } else {
                    diff_2 += gsl_matrix_get(ancestry.type, k + 1, ancestor_2);
                    ancestor_2 = gsl_matrix_get(ancestry.anc, k + 1, ancestor_2);
                }
            }
            std::cout << diff_1 << " " << diff_2 << std::endl;
        }
    }
    return;
}

void print_pairwise_coalescence_times(const Ancestry &ancestry) {
    int ancestor_1, ancestor_2;
    double diff_1 = 0.0, diff_2 = 0.0;
    for (int i = 0; i < ancestry.sample_size - 1; i++) {
        for (int j = i + 1; j < ancestry.sample_size; j++) {
            for (int k = 0; k < ancestry.site_no - 1; k++) {
                ancestor_1 = gsl_matrix_get(ancestry.anc, k, i);
                ancestor_2 = gsl_matrix_get(ancestry.anc, k, j);
                while (ancestor_1 != ancestor_2) {
                    if (ancestor_1 < ancestor_2) {
                        ancestor_1 = gsl_matrix_get(ancestry.anc, k, ancestor_1);
                    } else {
                        ancestor_2 = gsl_matrix_get(ancestry.anc, k, ancestor_2);
                    }
                }
                diff_1 = gsl_matrix_get(ancestry.t, k, ancestor_1);
                for (int q = k + 1; q < ancestry.site_no; q++) {
                    ancestor_1 = gsl_matrix_get(ancestry.anc, q, i);
                    ancestor_2 = gsl_matrix_get(ancestry.anc, q, j);
                    while (ancestor_1 != ancestor_2) {
                        if (ancestor_1 < ancestor_2) {
                            ancestor_1 = gsl_matrix_get(ancestry.anc, q, ancestor_1);
                        } else {
                            ancestor_2 = gsl_matrix_get(ancestry.anc, q, ancestor_2);
                        }
                    }
                    diff_2 = gsl_matrix_get(ancestry.t, q, ancestor_1);
                    std::cout << diff_1 << " " << diff_2 << std::endl;
                }
            }
        }
    }
    return;
}

void print_nonoverlapping_pairwise_coalescence_times(const Ancestry &ancestry) {
    int ancestor_1, ancestor_2;
    double diff_1 = 0.0, diff_2 = 0.0;
    for (int i = 0; i < ancestry.sample_size - 1; i += 2) {
        for (int k = 0; k < ancestry.site_no - 1; k++) {
            ancestor_1 = gsl_matrix_get(ancestry.anc, k, i);
            ancestor_2 = gsl_matrix_get(ancestry.anc, k, i + 1);
            while (ancestor_1 != ancestor_2) {
                if (ancestor_1 < ancestor_2) {
                    ancestor_1 = gsl_matrix_get(ancestry.anc, k, ancestor_1);
                } else {
                    ancestor_2 = gsl_matrix_get(ancestry.anc, k, ancestor_2);
                }
            }
            diff_1 = gsl_matrix_get(ancestry.t, k, ancestor_1);
            ancestor_1 = gsl_matrix_get(ancestry.anc, k + 1, i);
            ancestor_2 = gsl_matrix_get(ancestry.anc, k + 1, i + 1);
            while (ancestor_1 != ancestor_2) {
                if (ancestor_1 < ancestor_2) {
                    ancestor_1 = gsl_matrix_get(ancestry.anc, k + 1, ancestor_1);
                } else {
                    ancestor_2 = gsl_matrix_get(ancestry.anc, k + 1, ancestor_2);
                }
            }
            diff_2 = gsl_matrix_get(ancestry.t, k + 1, ancestor_1);
            std::cout << diff_1 << " " << diff_2 << std::endl;
        }
    }
    return;
}

void print_singleton_pairs(const Ancestry &ancestry) {
    for (int i = 0; i < ancestry.site_no; i++) {
        for (int j = 0; j < ancestry.sample_size - 1; j += 2) {
            std::cout << gsl_matrix_get(ancestry.type, i, j) + gsl_matrix_get(ancestry.type, i, j + 1) << " ";
        }
        std::cout << std::endl;
    }
    return;
}

void print_normalised_sfs(const Ancestry &ancestry) {
    std::vector<double> sfs(ancestry.sample_size - 1, 0.0);
    double total = 0.0;
    int index = 0;
    for (int i = 0; i < ancestry.site_no; i++) {
        total = 0.0;
        for (int k = 0; k < ancestry.sample_size - 1; k++) {
            sfs[k] = 0.0;
        }
        index = 0;
        while (gsl_matrix_get(ancestry.anc, i, index) != -1) {
            sfs[gsl_matrix_get(ancestry.ord, i, index) - 1] += (double)(gsl_matrix_get(ancestry.type, i, index));
            total += (double)(gsl_matrix_get(ancestry.type, i, index));
            index++;
        }
        for (int k = 0; k < ancestry.sample_size - 1; k++) {
            if (total > 0.0) {
                sfs[k] /= total;
            }
            std::cout << sfs[k] << " ";
        }
        std::cout << std::endl;
    }
    return;
}

void print_sfs(const Ancestry &ancestry) {
    std::vector<double> sfs(ancestry.sample_size - 1, 0.0);
    double total = 0.0;
    int index = 0;
    for (int i = 0; i < ancestry.site_no; i++) {
        total = 0.0;
        for (int k = 0; k < ancestry.sample_size - 1; k++) {
            sfs[k] = 0.0;
        }
        index = 0;
        while (gsl_matrix_get(ancestry.anc, i, index) != -1) {
            sfs[gsl_matrix_get(ancestry.ord, i, index) - 1] += (double)(gsl_matrix_get(ancestry.type, i, index));
            total += (double)(gsl_matrix_get(ancestry.type, i, index));
            index++;
        }
        for (int k = 0; k < ancestry.sample_size - 1; k++) {
            std::cout << sfs[k] << " ";
        }
        std::cout << std::endl;
    }
    return;
}

void print_average_branch_lengths(const Ancestry &ancestry) {
    std::vector<double> branches(ancestry.sample_size - 1, 0.0);
    double total = 0.0, t;
    int index = 0;
    for (int i = 0; i < ancestry.site_no; i++) {
        index = 0;
        while (gsl_matrix_get(ancestry.anc, i, index) != -1) {
            t = gsl_matrix_get(ancestry.t, i, gsl_matrix_get(ancestry.anc, i, index))
                - gsl_matrix_get(ancestry.t, i, index);
            branches[gsl_matrix_get(ancestry.ord, i, index) - 1] += t;
            total += t;
            index++;
        }
    }
    for (int k = 0; k < ancestry.sample_size - 1; k++) {
        if (total > 0.0) {
            branches[k] /= total;
        }
        std::cout << branches[k] << " ";
    }
    std::cout << std::endl;
    return;
}

void fixed_s(const Ancestry &ancestry, const int s) {
    std::vector<double> branches(ancestry.sample_size - 1, 0.0);
    std::vector<int> sfs(ancestry.sample_size - 1, 0.0);
    double total = 0.0, t;
    int index = 0;
    for (int i = 0; i < ancestry.site_no; i++) {
        total = 0.0;
        for (int k = 0; k < ancestry.sample_size - 1; k++) {
            branches[k] = 0.0;
            sfs[k] = 0;
        }
        index = 0;
        while (gsl_matrix_get(ancestry.anc, i, index) != -1) {
            t = (double)(gsl_matrix_get(ancestry.t, i, gsl_matrix_get(ancestry.anc, i, index)) - gsl_matrix_get(ancestry.t, i, index));
            branches[gsl_matrix_get(ancestry.ord, i, index) - 1] += t;
            total += t;
            index++;
        }
        for (int j = 0; j < ancestry.sample_size - 1; j++) {
            branches[j] /= total;
        }
        for (int j = 0; j < s; j++) {
            double coin = gsl_rng_uniform(ancestry.gen);
            double lb = 0.0;
            double ub = branches[0];
            int index = 0;
            while(!(lb <= coin && coin < ub)) {
                lb = ub;
                ub += branches[index + 1];
                index++;
            }
            sfs[index]++;
        }
        for (int j = 0; j < ancestry.sample_size - 1; j++) {
            std::cout << sfs[j] << " ";
        }
        std::cout << std::endl;
    }
    return;
}

int main(int argc, char **argv) {
    if (argc != 8) {
        std::cout << "Call " << argv[0] << " <sample size> <site number> <alpha> <growth rate> <model> <diploid> <replicate number>" << std::endl;
        std::cout << "<model> = 0 => Beta(2 - alpha, alpha)-coalescent (growth rate is ignored)." << std::endl;
        std::cout << "<model> = 1 => Exponential growth coalescent (alpha is ignored)." << std::endl;
        std::cout << "<model> = 2 => Algebraic growth coalescent (alpha is ignored)." << std::endl;
        std::cout << "<diploid> = 0/1 => Haploid/diploid model." << std::endl;
        return 1;
    }

    int sample_size = atoi(argv[1]);
    int site_no = atoi(argv[2]);
    double alpha = atof(argv[3]);
    double growth_rate = atof(argv[4]);
    int model = atoi(argv[5]);
    int diploid = atoi(argv[6]);
    int replicate_number = atoi(argv[7]);

    assert(sample_size > 1);
    assert(site_no > 0);
    assert(model == 0 || model == 1 || model == 2);
    if (model == 0) {
        assert(alpha > 1.0);
        assert(alpha < 2.0);
    } else if (model == 1 || model == 2) {
        assert(growth_rate > 0.0);
    }
    assert(diploid == 0 || diploid == 1);
    assert(replicate_number > 0);
    
    double mu = 1e-8; // Rough per-site-per-generation mutation probability: test with both 1e-8 and 1e-7
    double locus_length = 26e6; // Average length of Atlantic cod linkage groups
    double Ne = 1000.0;
    double mutation_rate = 0.0;
    if (model == 0) {
        mutation_rate = 2.0 * mu * locus_length * pow(Ne, alpha - 1.0);
    } else {
        mutation_rate = 2.0 * mu * locus_length * Ne;
    }
    Ancestry ancestry(sample_size, site_no, model, diploid, mutation_rate, alpha, growth_rate);

    for (int i = 0; i < replicate_number; i++) {
        ancestry.clear();
        ancestry.simulate();
        //print_types(ancestry);
        //print_ancestry(ancestry);
        //print_coalescence_times(ancestry);
        //print_singletons(ancestry);
        //print_pairwise_differences(ancestry);
        //print_nonoverlapping_pairwise_differences(ancestry);
        //print_pairwise_coalescence_times(ancestry);
        //print_nonoverlapping_pairwise_coalescence_times(ancestry);
        //print_singleton_pairs(ancestry);
        print_normalised_sfs(ancestry);
        //print_sfs(ancestry);
        //print_average_branch_lengths(ancestry);
        //fixed_s(ancestry, 50);
    }
    
    return 1;
}
