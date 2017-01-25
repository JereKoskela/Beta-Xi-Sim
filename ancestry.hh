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

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_matrix.h>
#include <iostream>
#include <limits>
#include <unistd.h>
#include <vector>

struct Ancestry {

    Ancestry(const int sample, const int site, const int type_flag,
        const int diploid_flag, const double mut_rate, const double a, 
        const double g)
    : sample_size(sample), site_no(site), coalescent_type(type_flag),	diploid(
    diploid_flag), mutation_rate(mut_rate), beta_param_a(2.0 - a), beta_param_b(a),
    growth_rate(g),next_parent(site_no, sample_size), coalescers(site_no, 0),
    sample(sample_size, -1), families(4, 0), n_length(site_no, sample_size),
    coalescence_probabilities(site_no * sample_size, 0.0)
    {
        type = gsl_matrix_alloc(site_no, 2 * sample_size - 1);
        anc = gsl_matrix_alloc(site_no, 2 * sample_size - 1);
        ord = gsl_matrix_alloc(site_no, 2 * sample_size - 1);
        t = gsl_matrix_alloc(site_no, 2 * sample_size - 1);
        n = gsl_matrix_alloc(site_no, sample_size);
        for (int i = 0; i < site_no; i++) {
            for (int j = 0; j < 2 * sample_size - 1; j++) {
                gsl_matrix_set(type, i, j, 0);
                gsl_matrix_set(anc, i, j, -1);
                gsl_matrix_set(t, i, j, 0.0);
                if (j < sample_size) {
                    gsl_matrix_set(ord, i, j, 1);
                    gsl_matrix_set(n, i, j, j);
                } else {
                    gsl_matrix_set(ord, i, j, 0);
                }
            }
        }
        gen = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(gen, time(NULL) * getpid());
    }

    ~Ancestry() {
        gsl_matrix_free(type);
        gsl_matrix_free(anc);
        gsl_matrix_free(ord);
        gsl_matrix_free(t);
        gsl_matrix_free(n);
        gsl_rng_free(gen);
    }

    void clear() {
        for (int i = 0; i < site_no; i++) {
            for (int j = 0; j < 2 * sample_size - 1; j++) {
                gsl_matrix_set(type, i, j, 0);
                gsl_matrix_set(anc, i, j, -1);
                gsl_matrix_set(t, i, j, 0.0);
                if (j < sample_size) {
                    gsl_matrix_set(ord, i, j, 1);
                    gsl_matrix_set(n, i, j, j);
                } else {
                    gsl_matrix_set(ord, i, j, 0);
                }
            }
            next_parent[i] = sample_size;
            n_length[i] = sample_size;
        }
        return;
    }
    
    double acceptance_prob(const double x, const double c) {
        double ret = 0.0, k;
        double m = (double)(site_no);
        if (x > 1e-4) {
            for (int i = 0; i < site_no; i++) {
                k = (double)(n_length[i]);
                if (k * log(1.0 - x) > log(k) + log(x) + (k - 1.0) * log(1.0 - x)) {
                    ret += k * log(1.0 - x) - 2.0 * log(x) / m + log(1.0 + exp(log(k) 
                            + log(x) - log(1.0 - x)));
                } else {
                    ret += log(k) + log(x) + (k - 1.0) * log(1.0 - x) - 2.0 * log(x) / 
                    m + log(exp(log(1.0 - x) - log(k) - log(x)) + 1.0);
                }
            }
            ret = exp(-log(c) - 2.0 * log(x)) - exp(ret - log(c));
        } else {
            std::vector<double> tmp_hap(site_no, -1.0);
            double tmp_max = -std::numeric_limits<double>::max();
            for (int i = 0; i < site_no; i++) {
                if (n_length[i] > 2) {
                    tmp_hap[i] = gsl_sf_lnchoose(n_length[i], 3);
                    tmp_max = fmax(tmp_max, tmp_hap[i]);
                }
            }
            for (int i = 0; i < site_no; i++) {
                if (tmp_hap[i] > -0.5) {
                    ret += exp(tmp_hap[i] - tmp_max);
                }
            }
            ret = tmp_max + log(ret);
            ret = 1.0 - exp(log(2.0) + ret - log(c) + log(x));
        }
        assert(ret <= 1.01);
        assert(ret >= -0.01);
        return ret;
    }
    
    void sample_families(const int i) {
        if (diploid == 0) {
            families[0] = coalescers[i];
            families[1] = 0;
            families[2] = 0;
            families[3] = 0;
        } else {
            families[0] = gsl_ran_binomial(gen, 1.0 / 4.0, coalescers[i]);
            families[1] = gsl_ran_binomial(gen, 1.0 / 3.0, coalescers[i] 
                            - families[0]);
            families[2] = gsl_ran_binomial(gen, 1.0 / 2.0, coalescers[i] 
                            - families[0] - families[1]);
            families[3] = coalescers[i] - families[0] - families[1] 
                            - families[2];
        }
        return;
    }
    
    void simulate_lambda_event(int &pop_size, double &time) {
        double coal_no = 0.0;
        for (int i = 0; i < site_no; i++) {
            coal_no += (double)(n_length[i] * (n_length[i] - 1) / 2);
        }
        time += gsl_ran_exponential(gen, 1.0 / coal_no);
        double impact = gsl_ran_beta(gen, beta_param_a, beta_param_b);
        double coin = gsl_rng_uniform(gen);
        if (coin < acceptance_prob(impact, coal_no)) {
            int index = 0;
            coin = gsl_rng_uniform(gen);
            double lb = 0.0;
            double ub = (double)(n_length[0] * (n_length[0] - 1)) / (2.0 * coal_no);
            while (!(lb <= coin && coin < ub)) {
                index++;
                lb = ub;
                ub += (double)(n_length[index] * (n_length[index] - 1)) / (2.0 * coal_no);
            }
            for (int i = 0; i < site_no; i++) {
                if (i == index) {
                    coalescers[i] = 2 + gsl_ran_binomial(gen, impact, n_length[i] - 2);
                } else {
                    coalescers[i] = gsl_ran_binomial(gen, impact, n_length[i]);
                }
            }
            for (int i = 0; i < site_no; i++) {
                if (coalescers[i] > 1) {
                    sample_families(i);
                    index = 0;
                    int total_records = 0, total_sampled = 0;
                    for (int j = 0; j < 4; j++) {
                        if (families[j] > 1) {
                            index++;
                            total_records = 0;
                            total_sampled = 0;
                            while (total_sampled < families[j]) {
                                coin = gsl_rng_uniform(gen);
                                if ((double)(n_length[i] - total_records) * coin >= 
                                    (double)(families[j] - total_sampled)) {
                                    total_records++;
                                } else {
                                    sample[total_sampled] = total_records;
                                    total_records++;
                                    total_sampled++;
                                }
                            }
                            std::sort(sample.begin(), sample.begin() + families[j]);
                            gsl_matrix_set(t, i, next_parent[i], time);
                            int child;
                            for (int k = families[j] - 1; k > -1; k--) {
                                child = gsl_matrix_get(n, i, sample[k]);
                                gsl_matrix_set(anc, i, child, next_parent[i]);
                                gsl_matrix_set(ord, i, next_parent[i], gsl_matrix_get(ord, i,
                                    next_parent[i]) + gsl_matrix_get(ord, i, child));
                                gsl_matrix_set(type, i, child, gsl_ran_poisson(gen, 
                                mutation_rate * (time - gsl_matrix_get(t, i, child)) / 2.0));
                                gsl_matrix_set(n, i, sample[k], gsl_matrix_get(n, i, n_length[i] - 1));
                                n_length[i]--;
                            }
                            next_parent[i]++;
                            pop_size -= families[j] - 1;
                        }
                    }
                    for (int j = index; j > 0; j--) {
                        gsl_matrix_set(n, i, n_length[i], next_parent[i] - j);
                        n_length[i]++;
                    }
                }
            }
        }
        return;
    }
    
    void simulate_exp_growth_event(int &pop_size, double &time) {
        double diploidy = 1.0;
        if (diploid == 1) {
            diploidy = 4.0;
        }
        int coal_no = 0; // number of possible binary mergers
        for (int i = 0; i < site_no; i++) {
            coal_no += n_length[i] * (n_length[i] - 1) / 2;
        }
        double tmp = log(diploidy) + log(growth_rate) + log(-log(gsl_rng_uniform(gen)))
                    - log(coal_no) - growth_rate * time;
        if (tmp > 0.0) {
            tmp = tmp + log(1.0 + exp(-tmp));
        } else {
            tmp = log(1.0 + exp(tmp));
        }
        time += exp(log(tmp) - log(growth_rate));
        int index = 0;
        double lb = 0.0;
        double ub = (double)(n_length[0] * (n_length[0] - 1)) / (double)(2 * coal_no);
        double coin = gsl_rng_uniform(gen);
        while (!(lb <= coin && coin < ub)) {
            index++;
            lb = ub;
            ub += (double)(n_length[index] * (n_length[index] - 1)) / (double)(2 * coal_no);
        }
        int child_ind_1, child_ind_2;
        do {
            child_ind_1 = floor(gsl_rng_uniform(gen) * (double)(n_length[index]));
            child_ind_2 = floor(gsl_rng_uniform(gen) * (double)(n_length[index]));
        } while(child_ind_1 == child_ind_2);
        int child_1 = gsl_matrix_get(n, index, child_ind_1);
        int child_2 = gsl_matrix_get(n, index, child_ind_2);
        gsl_matrix_set(anc, index, child_1, next_parent[index]);
        gsl_matrix_set(anc, index, child_2, next_parent[index]);
        gsl_matrix_set(t, index, next_parent[index], time);
        gsl_matrix_set(ord, index, next_parent[index], gsl_matrix_get(ord, index,
                        child_1) + gsl_matrix_get(ord, index, child_2));
        gsl_matrix_set(type, index, child_1, gsl_ran_poisson(gen, mutation_rate
                        * (time - gsl_matrix_get(t, index, child_1)) / 2.0));
        gsl_matrix_set(type, index, child_2, gsl_ran_poisson(gen, mutation_rate
                        * (time - gsl_matrix_get(t, index, child_2)) / 2.0));
        gsl_matrix_set(n, index, std::max(child_ind_1, child_ind_2), 
                        gsl_matrix_get(n, index, n_length[index] - 1));
        gsl_matrix_set(n, index, std::min(child_ind_1, child_ind_2), 
                        gsl_matrix_get(n, index, n_length[index] - 2));
        gsl_matrix_set(n, index, n_length[index] - 2, next_parent[index]);
        next_parent[index]++;
        n_length[index]--;
        pop_size--;
        return;
    }
    
    void simulate_alg_growth_event(int &pop_size, double &time) {
        double diploidy = 1.0;
        if (diploid == 1) {
            diploidy = 4.0;
        }
        int coal_no = 0; // number of possible binary mergers
        for (int i = 0; i < site_no; i++) {
            coal_no += n_length[i] * (n_length[i] - 1) / 2;
        }
        // We use a log-sum-exp trick here to ensure a stable time increment computation
        double tmp_1 = (growth_rate + 1.0) * log(time);
        double tmp_2 = log(diploidy) + log(-log(gsl_rng_uniform(gen))) 
                        + log(growth_rate + 1.0) - log((double)(coal_no));
        if (tmp_1 > tmp_2) {
            tmp_1 = tmp_1 + log(1.0 + exp(tmp_2 - tmp_1));
        } else {
            tmp_1 = tmp_2 + log(exp(tmp_1 - tmp_2) + 1.0);
        }
        time = exp(tmp_1 / (growth_rate + 1.0));
        int index = 0;
        double lb = 0.0;
        double ub = (double)(n_length[0] * (n_length[0] - 1)) / (double)(2 * coal_no);
        double coin = gsl_rng_uniform(gen);
        while (!(lb <= coin && coin < ub)) {
            index++;
            lb = ub;
            ub += (double)(n_length[index] * (n_length[index] - 1)) / (double)(2 * coal_no);
        }
        int child_ind_1, child_ind_2;
        do {
            child_ind_1 = floor(gsl_rng_uniform(gen) * (double)(n_length[index]));
            child_ind_2 = floor(gsl_rng_uniform(gen) * (double)(n_length[index]));
        } while(child_ind_1 == child_ind_2);
        int child_1 = gsl_matrix_get(n, index, child_ind_1);
        int child_2 = gsl_matrix_get(n, index, child_ind_2);
        gsl_matrix_set(anc, index, child_1, next_parent[index]);
        gsl_matrix_set(anc, index, child_2, next_parent[index]);
        gsl_matrix_set(t, index, next_parent[index], time);
        gsl_matrix_set(ord, index, next_parent[index], gsl_matrix_get(ord, index,
                        child_1) + gsl_matrix_get(ord, index, child_2));
        gsl_matrix_set(type, index, child_1, gsl_ran_poisson(gen, mutation_rate
                        * (time - gsl_matrix_get(t, index, child_1)) / 2.0));
        gsl_matrix_set(type, index, child_2, gsl_ran_poisson(gen, mutation_rate
                        * (time - gsl_matrix_get(t, index, child_2)) / 2.0));
        gsl_matrix_set(n, index, std::max(child_ind_1, child_ind_2), 
                        gsl_matrix_get(n, index, n_length[index] - 1));
        gsl_matrix_set(n, index, std::min(child_ind_1, child_ind_2), 
                        gsl_matrix_get(n, index, n_length[index] - 2));
        gsl_matrix_set(n, index, n_length[index] - 2, next_parent[index]);
        next_parent[index]++;
        n_length[index]--;
        pop_size--;
        return;
    }
    
    void simulate() {
        int pop_size = site_no * sample_size;
        double time = 0.0;
        while (pop_size > site_no) {
            switch(coalescent_type) {
                case 0: simulate_lambda_event(pop_size, time); 
                    break;
                case 1: simulate_exp_growth_event(pop_size, time); 
                    break;
                case 2: simulate_alg_growth_event(pop_size, time);
                    break;
                default: std::cout << "Coalescent type " << coalescent_type << " not recognised." << std::endl;
                    std::cout << "Call ./unlinked for a list of implemented types." << std::endl;
                    gsl_matrix_free(type);
                    gsl_matrix_free(anc);
                    gsl_matrix_free(ord);
                    gsl_matrix_free(t);
                    gsl_matrix_free(n);
                    abort();
            }
        }
        return;
    }
    
    int sample_size, site_no, coalescent_type, diploid;
    double mutation_rate, beta_param_a, beta_param_b, growth_rate;
    std::vector<int> next_parent, coalescers, sample, families, n_length;
    std::vector<double> coalescence_probabilities;
    gsl_matrix *type;
    gsl_matrix *anc;
    gsl_matrix *ord;
    gsl_matrix *t;
    gsl_matrix *n;
    gsl_rng *gen;
};