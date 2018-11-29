#ifndef ANCESTRY
#define ANCESTRY

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_gamma.h>
#include <iostream>
#include <libconfig.h++>
#include <limits>
#include <unistd.h>
#include <vector>
#include "branch.hh"

struct Ancestry {
    // ===============================================
    // Functions to be defined by the user begin here.
    // See function descriptions for specifications.
    // ===============================================
    
    double max_fitness(const int island) const {
        // Returns the maximum reproductive fitness that a set of parents can 
        // have on the given island. 
        
        // ** Must coincide with the maximum value of the fitness function! **
        // ** This requirement is not checked by the code, but without it   **
        // ** the algorithm does not correspond to any model!               **
        
        double ret = (double)number_of_loci * selection_rates[island];
        return ret;
    }
    
    double fitness(const std::vector<int> &inc, 
                   const std::vector<int> &che) const {
        // inc and che specify the branch indices of every locus on every 
        // chromosome in a selective branch event. They correspond in the 
        // obvious way to the incoming and checking branches. The pattern
        // within each vector is 
        // (chr1 loc1, chr1 loc2, ..., chr1 locn, chr2 loc1, ... chr2 locn).
                       
        // ** The maximum value must coincide with max_fitness!             **
        // ** This requirement is not checked by the code, but without it   **
        // ** the algorithm does not correspond to any model!               **
        int island = branches[inc[0]].island;
        double ret = 0.0;
        for (int i = 0; i < number_of_loci; i++) {
            int check_inc = 1;
            int check_che = 1;
            if (branches[inc[i]].selective_type == 0 
                && branches[inc[number_of_loci + i]].selective_type == 0) {
                check_inc = 0;
            }
            if (branches[che[i]].selective_type == 0 
                && branches[che[number_of_loci + i]].selective_type == 0) {
                check_inc = 0;
            }
            if (check_inc == 1 || check_che == 1) {
                ret += selection_rates[island];
            }
        }
        return ret;
    }
    
    int sample_initial_selective_type(const int locus) const {
        // Return a type from the desired distribution for selective types of
        // roots. Under stationarity this should be the stationary distribution
        // of mutate_selective_type(), but this is not enforced in the code and
        // other initial distributions can be used if desired.
        
        // The selective type space is implicit in the definitions of this 
        // function as well as mutate_selective_type(). It is stored as an 
        // integer in the branch struct.
        int ret = 1;
        if (locus > -1 && locus < number_of_loci 
            && gsl_rng_uniform(gen) < 0.9) {
            ret = 0;
        }
        return ret;
    }
    
    int mutate_selective_type(const int child, const int parent) const {
        // Sample a change in selective type given a mutation event.
        // The selective type space is implicit in the definitions of this 
        // function as well as sample_initial_selective_type(). It is stored 
        // as an integer in the branch struct.
        int locus = branches[child].locus;
        int mut_count = gsl_ran_poisson(gen, (branches[parent].leaf_time 
            - branches[child].leaf_time) * selective_mutation_rates[locus]);
        // Mutation matrix = {{0.9, 0.1}, {0.9, 0.1}}
        int ret = branches[parent].selective_type;
        for (int i = 0; i < mut_count; i++) {
            if (ret == 0) {
                if (gsl_rng_uniform(gen) < 0.1) {
                    ret = 1;
                }
            } else {
                if (gsl_rng_uniform(gen) < 0.9) {
                    ret = 0;
                }
            }   
        }
        return ret;
    }
    
    double sample_neutral_mutation(const int i, double tot_anc_length) const {
        // Sample the location of a neutral mutation on branch i given that
        // one has happened. Location needs to lie within ancestral material
        // on that branch, the length of which is provided as an argument.
        int block = 0;
        double coin = gsl_rng_uniform(gen);
        double ub = (branches[i].ancestral_blocks[block + 1] 
                    - branches[i].ancestral_blocks[block]) / tot_anc_length;
        while (ub < coin) {
           block += 2;
           ub += (branches[i].ancestral_blocks[block + 1] 
                    - branches[i].ancestral_blocks[block]) / tot_anc_length;
        }
        double ret = branches[i].ancestral_blocks[block] 
                + (branches[i].ancestral_blocks[block + 1] 
                - branches[i].ancestral_blocks[block]) * gsl_rng_uniform(gen);
        return ret;
    }
    
    double recombination_rate(const int i) const {
        // Take the index of a branch, and return the total recombination rate
        // between all links within ancestral material on that branch.
        
        // ** This needs to be consistent with                               **
        // ** sample_recombination_point()! The code does not check for this **
        // ** requirement, it is left to the user to ensure a flexible       **
        // ** algorithm.                                                     **
        double ret = (branches[i].ancestral_blocks.back() 
                        - branches[i].ancestral_blocks[0]) 
                    * recombination_rates[branches[i].locus];
        return ret;
    }
    
    double sample_recombination_point(const int i) const {
        // Sample a recombination point among all links surrounded by 
        // ancestral material on branch with index i, given that such a
        // recombination event is happening on that branch.
        
        // ** This needs to be consistent with recombination_rate()!      **
        // ** The code does not check for this requirement, it is left to **
        // ** the user to ensure a flexible algorithm.                    **
        double ret = branches[i].ancestral_blocks[0] + gsl_rng_uniform(gen) 
                    * (branches[i].ancestral_blocks.back() 
                    - branches[i].ancestral_blocks[0]);
        return ret;
    }
    
    double migration_probability(const int i, const int j) const {
        // Returns the forwards-in-time probability of migrating from island i
        // to island j given a migration event from island i. The 
        // coalescent (reverse time) migration rate from i to j is 
        // migration_rates[j] * migration_probability(j, i), up to factors
        // based on different island population sizes and growth rates.
        double ret = 0.0;
        if (number_of_islands > 1 && j != i) {
            ret = 1.0 / (double)(number_of_islands - 1);
        }
        return ret;
    }
    
    double sample_impact() const {
        // Return the impact of a particular event in the Poisson construction
        // of a Lambda-coalescent on the specified island.
        double ret = 0.0;
        if (alpha < 2.0) {
            ret = gsl_ran_beta(gen, 2.0 - alpha, alpha);
        }
        return ret;
    }
    
    // =======================================================================
    // User-defined functions end here. Functions below this point do not need
    // to be changed unless the scope of the algorithm is being altered.
    // =======================================================================
    
    Ancestry(const char *filename)
    : number_of_loci(), number_of_islands(), sample_size(0), alpha(), 
    growth_rate(), sim_time(0.0), selection_rates(), mutation_rates(), 
    selective_mutation_rates(), migration_rates(), recombination_rates(), 
    relative_population_sizes(), branches(0), active_branches(0), 
    break_points() {
        read_config(filename);
        std::vector<double> tmp_break_points(2, 0.0);
        tmp_break_points[1] = 1.0;
        for (int i = 0; i < number_of_loci; i++) {
            break_points.push_back(tmp_break_points);
        }
        gen = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(gen, time(NULL) * getpid());
        per_island_locus_recomb_rates = gsl_matrix_alloc(number_of_islands, 
                                                         number_of_loci);
        gsl_matrix_set_zero(per_island_locus_recomb_rates);
        for (int i = 0; i < sample_size * number_of_loci; i++) {
            gsl_matrix_set(per_island_locus_recomb_rates, branches[i].island, 
                           branches[i].locus, gsl_matrix_get(
                               per_island_locus_recomb_rates, 
                               branches[i].island, branches[i].locus) 
                           + recombination_rate(i));
        }
    }

    ~Ancestry() {
        gsl_rng_free(gen);
        gsl_matrix_free(per_island_locus_recomb_rates);
    }
    
    void read_config(const char *filename) {
        double effective_population_size, per_site_mutation_probability;
        libconfig::Config cfg;
        cfg.readFile(filename);
        cfg.lookupValue("alpha", alpha);
        cfg.lookupValue("effective_population_size", 
                        effective_population_size);
        cfg.lookupValue("growth_rate", growth_rate);
        growth_rate *= pow(effective_population_size, alpha - 1.0);
        cfg.lookupValue("per_site_mutation_probability", 
                        per_site_mutation_probability);
        number_of_islands = cfg.getRoot()["sample_sizes"].getLength();
        number_of_loci = cfg.getRoot()["locus_lengths"].getLength();
        Branch tmp(0, 0, 0.0, 1.0, 0.0);
        double tmp_double;
        for (int i = 0; i < number_of_islands; i++) {
            int count = cfg.getRoot()["sample_sizes"][i];
            sample_size += count;
            std::vector<std::vector<int> > within_island_tmp(number_of_loci);
            std::vector<int> within_locus_tmp(count);
            tmp.island = i;
            for (int l = 0; l < number_of_loci; l++) {
                for (int k = 0; k < count; k++) {
                    within_locus_tmp[k] = (int)branches.size() + k;
                }
                tmp.locus = l;
                branches.insert(branches.end(), count, tmp);
                within_island_tmp[l] = within_locus_tmp;
            }
            active_branches.push_back(within_island_tmp);
            relative_population_sizes.push_back(
                cfg.getRoot()["relative_population_sizes"][i]);
            tmp_double = cfg.getRoot()["selection_strengths"][i];
            selection_rates.push_back(tmp_double 
                * pow(effective_population_size, alpha - 1.0));
            tmp_double = cfg.getRoot()["migrant_fractions"][i];
            migration_rates.push_back(tmp_double 
                * pow(effective_population_size, alpha - 1.0));
        }
        for (int i = 0; i < number_of_loci; i++) {
            tmp_double = cfg.getRoot()["locus_lengths"][i];
            mutation_rates.push_back(per_site_mutation_probability 
                * pow(effective_population_size, alpha - 1.0) * tmp_double);
            tmp_double = cfg.getRoot()["selective_window_lengths"][i];
            selective_mutation_rates.push_back(per_site_mutation_probability 
                * pow(effective_population_size, alpha - 1.0) * tmp_double);
            tmp_double = cfg.getRoot()["recombination_probabilities"][i];
            recombination_rates.push_back(pow(effective_population_size, 
                                              alpha - 1.0) * tmp_double);
        }
        return;
    }
    
    void reset(const char *filename) {
        libconfig::Config cfg; 
        cfg.readFile(filename);
        Branch tmp(0, 0, 0.0, 1.0, 0.0);
        branches.clear();
        active_branches.clear();
        for (int i = 0; i < number_of_islands; i++) {
            int count = cfg.getRoot()["sample_sizes"][i];
            std::vector<std::vector<int> > within_island_tmp(number_of_loci);
            std::vector<int> within_locus_tmp(count);
            tmp.island = i;
            for (int l = 0; l < number_of_loci; l++) {
                for (int k = 0; k < count; k++) {
                    within_locus_tmp[k] = (int)branches.size() + k;
                }
                tmp.locus = l;
                branches.insert(branches.end(), count, tmp);
                within_island_tmp[l] = within_locus_tmp;
            }
            active_branches.push_back(within_island_tmp);
        }
        for (int i = 0; i < number_of_loci; i++) {
            break_points[i].clear();
            break_points[i].push_back(0.0);
            break_points[i].push_back(1.0);
        }
        gsl_matrix_set_zero(per_island_locus_recomb_rates);
        for (int i = 0; i < sample_size * number_of_loci; i++) {
            gsl_matrix_set(per_island_locus_recomb_rates, branches[i].island, 
                           branches[i].locus, gsl_matrix_get(
                               per_island_locus_recomb_rates, 
                               branches[i].island, branches[i].locus) 
                           + recombination_rate(i));
        }
        sim_time = 0.0;
        return;
    }
    
    double acceptance_prob(const double impact, const int island) {
        double ret = 1.0;
        if (alpha < 2.0) {
            double n = 0.0;
            double pairs = (double)(active_branches[island][0].size()
                * (active_branches[island][0].size() - 1)) / 2.0;
            for (int i = 1; i < number_of_loci; i++) {
                pairs += (double)(active_branches[island][i].size() 
                    * (active_branches[island][i].size() - 1)) / 2.0;
            }
            if (impact > 1e-4) {
                ret = 0.0;
                for (int i = 0; i < number_of_loci; i++) {
                    n = (double)(active_branches[island][i].size() - 1);
                    ret += n * log(1.0 - impact) + log(1.0 + n * impact);
                }
                ret = exp(log(1.0 - exp(ret)) - 2.0 * log(impact) 
                    - log(pairs));
            } else {
                ret = impact * (double)(active_branches[island][0].size()
                    * (active_branches[island][0].size() - 1)
                    * (active_branches[island][0].size() - 2)) 
                    / (3.0 * pairs);
                for (int i = 1; i < number_of_loci; i++) {
                    ret += impact * (double)(active_branches[island][i].size()
                        * (active_branches[island][i].size() - 1)
                        * (active_branches[island][i].size() - 2)) 
                        / (3.0 * pairs);
                }
                ret = 1.0 - ret;
            }
        }
        return ret;
    }
    
    void simulate_migration_event(const int island, const int target_island, 
                                  const int locus) {
        int child_index = floor((double)(active_branches[island][locus].size())
            * gsl_rng_uniform(gen));
        int child = active_branches[island][locus][child_index];
        Branch tmp(locus, target_island, 0.0, 1.0, sim_time);
        tmp.children.push_back(child);
        tmp.ancestral_blocks = branches[child].ancestral_blocks;
        tmp.virtual_flag = branches[child].virtual_flag;
        branches[child].parents.push_back(branches.size());
        active_branches[island][locus].erase(
            active_branches[island][locus].begin() + child_index);
        active_branches[target_island][locus].push_back(branches.size());
        branches.push_back(tmp);
        gsl_matrix_set(per_island_locus_recomb_rates, island, locus, 
                       gsl_matrix_get(per_island_locus_recomb_rates, island, 
                                      locus) - recombination_rate(child));
        gsl_matrix_set(per_island_locus_recomb_rates, target_island, locus, 
                       gsl_matrix_get(per_island_locus_recomb_rates, 
                                      target_island, locus) 
                       + recombination_rate(branches[child].parents[0]));
        return;
    }
    
    void simulate_selection_event(const int island, const int locus) {
        int child_index = floor((double)(active_branches[island][locus].size())
            * gsl_rng_uniform(gen));
        int child = active_branches[island][locus][child_index];
        branches[child].parents.push_back(branches.size());
        Branch tmp(locus, island, 0.0, 1.0, sim_time);
        tmp.children.push_back(child);
        tmp.ancestral_blocks = branches[child].ancestral_blocks;
        tmp.virtual_flag = branches[child].virtual_flag;
        active_branches[island][locus].erase(active_branches[island]
            [locus].begin() + child_index);
        active_branches[island][locus].push_back(branches.size());
        branches.push_back(tmp);
        // Incoming and checking branch indices consist of 
        // 2 * number_of_loci entries to track both chromosomes of both
        // parents. The first number_of_loci entries correspond to the first 
        // chromosome, then the second. WLOG inheritance is always assumed
        // to be from the first chromosome when an actual selection event 
        // takes place.
        tmp.children.clear();
        tmp.virtual_flag = 1;
        for (int j = 0; j < 2; j++) {
            for (int i = 0; i < number_of_loci; i++) {
                tmp.locus = i;
                branches[child].incoming.push_back(branches.size());
                branches[child].checking.push_back(branches.size() + 1);
                active_branches[island][i].push_back(branches.size());
                active_branches[island][i].push_back(branches.size() + 1);
                branches.insert(branches.end(), 2, tmp);
                gsl_matrix_set(per_island_locus_recomb_rates, island, i, 
                               gsl_matrix_get(per_island_locus_recomb_rates, 
                                              island, locus) 
                               + recombination_rate(branches.size() - 1) 
                               + recombination_rate(branches.size() - 2));
            }
        }
        return;
    }
    
    void simulate_recombination_event(const int island, const int locus) {
        int child_ind = 0;
        double upper_bound = recombination_rate(active_branches[island][locus]
            [child_ind]) / gsl_matrix_get(per_island_locus_recomb_rates, 
                                          island, locus);
        double coin = gsl_rng_uniform(gen);
        while (upper_bound < coin) {
            child_ind++;
            upper_bound += recombination_rate(active_branches[island][locus]
                [child_ind]) / gsl_matrix_get(per_island_locus_recomb_rates, 
                                              island, locus);
        }
        int child = active_branches[island][locus][child_ind];
        double split_point = sample_recombination_point(child);
        break_points[locus].insert(std::upper_bound(
            break_points[locus].begin(), break_points[locus].end(), 
            split_point), split_point);
        int index = 0;
        while (branches[child].ancestral_blocks[index] < split_point) {
            index++;
        }
        std::vector<double> points_1(0);
        std::vector<double> points_2(0);
        points_1.insert(points_1.begin(), 
                        branches[child].ancestral_blocks.begin(), 
                        branches[child].ancestral_blocks.begin() + index);
        points_2.insert(points_2.begin(), 
                        branches[child].ancestral_blocks.begin() + index, 
                        branches[child].ancestral_blocks.end());
        if (index % 2 == 1) {
            points_1.push_back(split_point);
            points_2.insert(points_2.begin(), split_point);
        }
        Branch tmp_1(locus, island, 0.0, 1.0, sim_time);
        Branch tmp_2(locus, island, 0.0, 1.0, sim_time);
        tmp_1.ancestral_blocks = points_1;
        tmp_2.ancestral_blocks = points_2;
        tmp_1.children.push_back(child);
        tmp_2.children.push_back(child);
        tmp_1.virtual_flag = branches[child].virtual_flag;
        tmp_2.virtual_flag = branches[child].virtual_flag;
        branches[child].parents.push_back(branches.size());
        branches[child].parents.push_back(branches.size() + 1);
        active_branches[island][locus].erase(active_branches[island]
                                            [locus].begin() + child_ind);
        active_branches[island][locus].push_back(branches.size());
        active_branches[island][locus].push_back(branches.size() + 1);
        branches.push_back(tmp_1);
        branches.push_back(tmp_2);
        gsl_matrix_set(per_island_locus_recomb_rates, island, locus, 
                       gsl_matrix_get(per_island_locus_recomb_rates, island,
                                      locus) - recombination_rate(child) 
                            + recombination_rate(branches[child].parents[0]) 
                            + recombination_rate(branches[child].parents[1]));
        return;
    }

    void simulate_coalescence_event(const int island, const int locus) {
        double impact = sample_impact();
        if (gsl_rng_uniform(gen) < acceptance_prob(impact, island)) {
            int size_at_locus = 0;
            for (int l = 0; l < number_of_loci; l++) {
                if (l == locus) {
                    size_at_locus = 2 + gsl_ran_binomial(gen, impact, 
                        active_branches[island][l].size() - 2);
                } else {
                    size_at_locus = gsl_ran_binomial(gen, impact, 
                        active_branches[island][l].size());
                }
                if (size_at_locus > 1) {
                    int mergers, running_total = 0;
                    int actual_mergers = 0;
                    int max_families = 4;
                    for (int i = 0; i < max_families; i++) {
                        mergers = gsl_ran_binomial(gen, 1.0 
                            / (double)(max_families - i), size_at_locus 
                            - running_total);
                        running_total += mergers;
                        if (mergers > 1) {
                            actual_mergers++;
                            std::vector<int> children(mergers, -1);
                            int total_records = 0;
                            int total_sampled = 0;
                            while (total_sampled < mergers) {
                                if ((double)(active_branches[island][l].size()
                                    - total_records) * gsl_rng_uniform(gen) 
                                    >= (double)(mergers - total_sampled)) {
                                    total_records++;
                                } else {
                                    children[total_sampled] = total_records;
                                    total_records++;
                                    total_sampled++;
                                }
                            }
                            Branch tmp(l, island, 0.0, 1.0, sim_time);
                            tmp.virtual_flag = 1;
                            for (int j = 0; j < (int)children.size(); j++) {
                                tmp.children.push_back(active_branches
                                    [island][l][children[j]]);
                                gsl_matrix_set(per_island_locus_recomb_rates, 
                                    island, l, gsl_matrix_get(
                                    per_island_locus_recomb_rates, island, l)
                                    - recombination_rate(tmp.children.back()));
                                if (branches[tmp.children[j]].virtual_flag
                                    == 0) {
                                    tmp.virtual_flag = 0;
                                }
                            }
                            std::vector<double> tmp_blocks = branches
                                [tmp.children[0]].ancestral_blocks;
                            int ind;
                            for (int j = 1; j < mergers; j++) {
                                ind = 0;
                                for (int k = 0; k <= (int)branches[tmp.
                                    children[j]].ancestral_blocks.size() 
                                    / 2; k += 2) {
                                    if (ind < (int)tmp_blocks.size()) {
                                        while (tmp_blocks[ind] < branches[
                                            tmp.children[j]].ancestral_blocks
                                            [k]) {
                                            ind += 2;
                                            if (ind == 
                                                (int)tmp_blocks.size()) {
                                                break;
                                            }
                                        }
                                    }
                                    tmp_blocks.insert(tmp_blocks.begin() 
                                        + ind, branches[tmp.children[j]]
                                        .ancestral_blocks.begin() + k, 
                                        branches[tmp.children[j]]
                                        .ancestral_blocks.begin() + k + 2);
                                }
                            }
                            ind = 0;
                            while (ind < (int)tmp_blocks.size() - 2) {
                                while (tmp_blocks[ind + 2] 
                                        < tmp_blocks[ind + 1]) {
                                    tmp_blocks[ind + 1] = fmax(tmp_blocks
                                        [ind + 1], tmp_blocks[ind + 3]);
                                    tmp_blocks.erase(tmp_blocks.begin() 
                                        + ind + 2, tmp_blocks.begin() + ind 
                                        + 4);
                                    if (ind + 2 > (int)tmp_blocks.size() - 1) {
                                        break;
                                    }
                                }
                                ind += 2;
                            }
                            tmp.ancestral_blocks = tmp_blocks;
                            for (int j = mergers - 1; j > -1; j--) {
                                branches[tmp.children[j]].parents
                                    .push_back(branches.size());
                                active_branches[island][l].erase(
                                    active_branches[island][l].begin() 
                                    + children[j]);
                            }
                            branches.push_back(tmp);
                            gsl_matrix_set(per_island_locus_recomb_rates, 
                                island, l, gsl_matrix_get(
                                per_island_locus_recomb_rates, island, l)
                                + recombination_rate(branches.size() - 1));
                        }
                    }
                    for (int i = actual_mergers; i > 0; i--) {
                        active_branches[island][l].push_back(
                            branches.size() - i);
                    }
                }
            }
        }
        return;
    }
    
    double time_increment(int &event_type, int &island, int &target_island,
                          int &locus) const {
        double ret = std::numeric_limits<double>::max();
        double tmp = 0.0;
        double helpful_factor = pow(relative_population_sizes[0], 3.0 - alpha);
        for (int i = 0; i < number_of_islands; i++) {
            helpful_factor += pow(relative_population_sizes[i], 3.0 - alpha);
        }
        for (int i = 0; i < number_of_islands; i++) {
            for (int j = 0; j < number_of_loci; j++) {
                if (active_branches[i][j].size() > 1) {
                    if (growth_rate * (alpha - 1.0) > 0.0) {
                        tmp = log(1.0 - 2.0 * growth_rate * (alpha - 1.0)
                            * pow(relative_population_sizes[i], alpha - 1.0)
                            * exp(-growth_rate * (alpha - 1.0) * sim_time) 
                            * helpful_factor * log(gsl_rng_uniform_pos(gen)) 
                            / (double)(active_branches[i][j].size() 
                            * (active_branches[i][j].size() - 1))) 
                            / (growth_rate * (alpha - 1.0));
                    } else {
                        tmp = gsl_ran_exponential(gen, 2.0 * helpful_factor * 
                            pow(relative_population_sizes[i], alpha - 1.0)
                            / (double)(active_branches[i][j].size() 
                            * (active_branches[i][j].size() - 1)));
                    }
                    if (tmp < ret) {
                        ret = tmp;
                        island = i;
                        event_type = 0;
                        locus = j;
                    }
                }
                if (active_branches[i][j].size() > 0) {
                    for (int k = 0; k < number_of_islands; k++) {
                        if (k != i && migration_rates[k] 
                            * migration_probability(k, i) > 0.0) {
                            tmp = gsl_ran_exponential(gen, 
                                relative_population_sizes[i] 
                                / ((double)active_branches[i][j].size() 
                                * relative_population_sizes[k] 
                                * migration_rates[k] 
                                * migration_probability(k, i)));
                            if (tmp < ret) {
                                ret = tmp;
                                island = i;
                                event_type = 1;
                                target_island = k;
                                locus = j;
                            }
                        }
                    }
                    if (max_fitness(i) > 0.0) {
                        tmp = gsl_ran_exponential(gen, 1.0 
                            / ((double)active_branches[i][j].size() 
                            * max_fitness(i)));
                        if (tmp < ret) {
                            ret = tmp;
                            island = i;
                            event_type = 2;
                            locus = j;
                        }
                    }
                    if (gsl_matrix_get(per_island_locus_recomb_rates, i, j) 
                        > 0.0) {
                        tmp = gsl_ran_exponential(gen, 1.0 / gsl_matrix_get(
                            per_island_locus_recomb_rates, i, j));
                        if (tmp < ret) {
                            ret = tmp;
                            island = i;
                            event_type = 3;
                            locus = j;
                        }
                    }
                }
            }
        }
        return ret;
    }
    
    void simulate_event() {
        int event_type = -1;
        int island = -1;
        int target_island = -1;
        int locus = -1;
        sim_time += time_increment(event_type, island, target_island, locus);
        switch (event_type) {
            case 0: simulate_coalescence_event(island, locus);
                    break;
            case 1: simulate_migration_event(island, target_island, locus);
                    break;
            case 2: simulate_selection_event(island, locus);
                    break;
            case 3: simulate_recombination_event(island, locus);
                    break;
            default: std::cout << "unrecognised event type" << std::endl;
                    abort();
        }
        return;
    }
    
    void extract_trees() {
        // ==============================================
        // This function assumes simulate() has been run!
        // ==============================================
        for (int i = (int)branches.size() - 1; i > -1; i--) {
            if (branches[i].parents.size() == 0) {
                branches[i].selective_type = 
                    sample_initial_selective_type(branches[i].locus);
            } else {
                if ((int)branches[i].incoming.size() == 2 * number_of_loci) {
                    if (gsl_rng_uniform(gen) < fitness(branches[i].incoming, 
                        branches[i].checking) 
                        / max_fitness(branches[i].island)) {
                        if (branches[i].virtual_flag == 1) {
                            branches[branches[i].parents[0]].children.clear();
                            branches[i].parents[0] 
                                = branches[i].incoming[branches[i].locus];
                            branches[branches[i].parents[0]].children
                                .push_back(i);
                        } else {
                            branches[branches[i].parents[0]].children.clear();
                            branches[branches[i].parents[0]].virtual_flag = 1;
                            std::vector<int> to_check 
                                = branches[branches[i].parents[0]].parents;
                            while (to_check.size() > 0) {
                                int virtual_q = 1;
                                for (unsigned int j = 0; j 
                                    < branches[to_check[0]].children.size(); 
                                    j++) {
                                    if (branches[branches[to_check[0]].children
                                        [j]].virtual_flag == 0) {
                                        virtual_q = 0;
                                        break;
                                    }
                                }
                                if (virtual_q == 1) {
                                    branches[to_check[0]].virtual_flag = 1;
                                    to_check.insert(to_check.end(), 
                                        branches[to_check[0]].parents.begin(), 
                                        branches[to_check[0]].parents.end());
                                }
                                to_check.erase(to_check.begin());
                            }
                            branches[i].parents[0] 
                                = branches[i].incoming[branches[i].locus];
                            branches[branches[i].parents[0]].children
                                .push_back(i);
                            branches[branches[i].parents[0]].virtual_flag = 0;
                            to_check 
                                = branches[branches[i].parents[0]].parents;
                            while (to_check.size() > 0) {
                                if (branches[to_check[0]].virtual_flag == 1) {
                                    branches[to_check[0]].virtual_flag = 0;
                                    to_check.insert(to_check.end(), 
                                        branches[to_check[0]].parents.begin(), 
                                        branches[to_check[0]].parents.end());
                                }
                                to_check.erase(to_check.begin());
                            }
                        }
                    }
                }
                // We assume the selective window is fully linked to the left
                // edge of each observed, neutral locus
                branches[i].selective_type = mutate_selective_type(i, 
                                                    branches[i].parents[0]);
            }
        }
        return;
    }
    
    void fill_orders() {
        // ===================================================
        // This function assumes extract_trees() has been run!
        // ===================================================
        for (int i = 0; i < sample_size * number_of_loci; i++) {
            branches[i].order.insert(branches[i].order.begin(), 
                            break_points[branches[i].locus].size() - 1, 1);
        }
        for (unsigned int i = sample_size * number_of_loci; 
             i < branches.size(); i++) {
            if (branches[i].virtual_flag == 0 
                && branches[i].parents.size() > 0) {
                branches[i].order.insert(branches[i].order.begin(), 
                    break_points[branches[i].locus].size() - 1, 0);
                for (int j = 0; j 
                    < (int)break_points[branches[i].locus].size() - 1; j++) {
                    for (unsigned int k = 0; k < branches[i].children.size(); 
                         k++) {
                        if (branches[branches[i].children[k]].virtual_flag 
                            == 0) {
                            if (branches[branches[i].children[k]].contains(
                                (break_points[branches[i].locus][j + 1] 
                                + break_points[branches[i].locus][j]) / 2.0) 
                                == 1) {
                                branches[i].order[j] 
                                    += branches[branches[i].children[k]]
                                    .order[j];
                            }
                        }
                    }
                }
            }
        }
        return;
    }
    
    void print_normalised_sfs() const {
        // =================================================
        // This function assumes fill_orders() has been run!
        // =================================================
        std::vector<int> tmp(sample_size - 1, 0);
        std::vector<std::vector<int> > sfs(number_of_loci, tmp);
        for (unsigned int i = 0; i < branches.size(); i++) {
            if (branches[i].virtual_flag == 0 
                && branches[i].parents.size() > 0) {
                for (int j = 0; j < (int)break_points[branches[i].locus].size()
                    - 1; j++) {
                    if (branches[i].contains((break_points[branches[i].locus]
                        [j + 1] + break_points[branches[i].locus][j])/ 2.0) 
                        == 1 && branches[i].order[j] > 0) {
                        // zero intervals can arise if a partial chromosome 
                        // branches due to selection - the virtual branches
                        // carry tracked material which nevertheless has zero
                        // descendants in the leaves
                        sfs[branches[i].locus][branches[i].order[j] - 1] 
                            += gsl_ran_poisson(gen, 
                            (break_points[branches[i].locus][j + 1] 
                            - break_points[branches[i].locus][j]) 
                            * (branches[branches[i].parents[0]].leaf_time 
                            - branches[i].leaf_time) 
                            * mutation_rates[branches[i].locus]);
                    }
                }
            }
        }
        double total = 0.0;
        for (int i = 0; i < number_of_loci; i++) {
            total = 0.0;
            for (int j = 0; j < sample_size - 1; j++) {
                total += (double)sfs[i][j];
            }
            if (total > 0.0) {
                for (int j = 0; j < sample_size - 1; j++) {
                    std::cout << (double)sfs[i][j] / total << " ";
                }
            } else {
                for (int j = 0; j < sample_size - 1; j++) {
                    std::cout << "0 ";
                }
            }
            std::cout << std::endl;
        }
        return;
    }
    
    void print_normalised_branch_lengths() const {
        // ========================================================
        // This function assumes fill_orders() has been run!
        // ========================================================
        gsl_matrix *branch_lengths = gsl_matrix_alloc(number_of_loci, 
                                                      sample_size - 1);
        gsl_matrix_set_zero(branch_lengths);
        std::vector<double> total(number_of_loci, 0.0);
        for (unsigned int i = 0; i < branches.size(); i++) {
            if (branches[i].virtual_flag == 0 && branches[i].parents.size()
                > 0) {
                for (int j = 0; j < (int)break_points[branches[i].locus].size()
                    - 1; j++) {
                    if (branches[i].contains((break_points[branches[i].locus]
                        [j + 1] + break_points[branches[i].locus][j])/ 2.0) 
                        == 1 && branches[i].order[j] > 0) {
                        // zero intervals can arise if a partial chromosome 
                        // branches due to selection - the virtual branches
                        // carry tracked material which nevertheless has zero
                        // descendants in the leaves
                        gsl_matrix_set(branch_lengths, branches[i].locus, 
                            branches[i].order[j] - 1, 
                            gsl_matrix_get(branch_lengths, branches[i].locus, 
                            branches[i].order[j] - 1) 
                            + (branches[branches[i].parents[0]].leaf_time 
                            - branches[i].leaf_time) 
                            * (break_points[branches[i].locus][j + 1] 
                            - break_points[branches[i].locus][j]));
                        total[branches[i].locus] 
                            += (branches[branches[i].parents[0]].leaf_time 
                            - branches[i].leaf_time) 
                            * (break_points[branches[i].locus][j + 1] 
                            - break_points[branches[i].locus][j]);
                    }
                }
            }
        }
        for (int i = 0; i < number_of_loci; i++) {
            for (int j = 0; j < sample_size - 1; j++) {
                if (total[i] > 0.0) {
                    std::cout << gsl_matrix_get(branch_lengths, i, j) 
                    / total[i] << " ";
                } else {
                    std::cout << 0 << " ";
                }
            }
            std::cout << std::endl;
        }
        return;
    }
    
    void print_fixed_s(const int s) const {
        // ========================================================
        // This function assumes fill_orders() has been run!
        // ========================================================
        gsl_matrix *branch_lengths = gsl_matrix_alloc(number_of_loci, 
                                                      sample_size - 1);
        gsl_matrix *sfs = gsl_matrix_alloc(number_of_loci, sample_size - 1);
        gsl_matrix_set_zero(branch_lengths);
        gsl_matrix_set_zero(sfs);
        std::vector<double> total(number_of_loci, 0.0);
        for (unsigned int i = 0; i < branches.size(); i++) {
            if (branches[i].virtual_flag == 0 && branches[i].parents.size() 
                > 0) {
                for (int j = 0; j < (int)break_points[branches[i].locus].size() 
                    - 1; j++) {
                    if (branches[i].contains((break_points[branches[i].locus]
                        [j + 1] + break_points[branches[i].locus][j])/ 2.0) 
                        == 1 && branches[i].order[j] > 0) {
                        // zero intervals can arise if a partial chromosome 
                        // branches due to selection - the virtual branches
                        // carry tracked material which nevertheless has zero
                        // descendants in the leaves
                        gsl_matrix_set(branch_lengths, branches[i].locus, 
                            branches[i].order[j] - 1, 
                            gsl_matrix_get(branch_lengths, branches[i].locus, 
                            branches[i].order[j] - 1) + 
                            (branches[branches[i].parents[0]].leaf_time 
                            - branches[i].leaf_time) 
                            * (break_points[branches[i].locus][j + 1] 
                            - break_points[branches[i].locus][j]));
                        total[branches[i].locus] 
                            += (branches[branches[i].parents[0]].leaf_time 
                            - branches[i].leaf_time) 
                            * (break_points[branches[i].locus][j + 1] 
                            - break_points[branches[i].locus][j]);
                    }
                }
            }
        }
        for (int i = 0; i < number_of_loci; i++) {
            for (int j = 0; j < s; j++) {
                double coin = gsl_rng_uniform(gen);
                double ub = gsl_matrix_get(branch_lengths, i, 0) / total[i];
                int index = 0;
                while (coin > ub) {
                    index++;
                    ub += gsl_matrix_get(branch_lengths, i, index) / total[i];
                }
                gsl_matrix_set(sfs, i, index, 
                               gsl_matrix_get(sfs, i, index) + 1);
            }
            for (int j = 0; j < sample_size - 1; j++) {
                std::cout << gsl_matrix_get(sfs, i, j) << " ";
            }
            std::cout << std::endl;
        }
        return;
    }
    
    void simulate() {
        int total_active_branches = sample_size * number_of_loci;
        while (total_active_branches > number_of_loci) {
            simulate_event();
            total_active_branches = 0;
            for (int i = 0; i < number_of_islands; i++) {
                for (int j = 0; j < number_of_loci; j++) {
                    total_active_branches += (int)active_branches[i][j].size();
                }
            }
        }
        extract_trees();
        fill_orders();
        return;
    }
    
    int number_of_loci, number_of_islands, sample_size;
    double alpha, growth_rate, sim_time; 
    std::vector<double> selection_rates, mutation_rates;
    std::vector<double> selective_mutation_rates, migration_rates;
    std::vector<double> recombination_rates, relative_population_sizes;
    std::vector<Branch> branches;
    std::vector<std::vector<std::vector<int> > > active_branches;
    std::vector<std::vector<double> > break_points;
    gsl_rng *gen;
    gsl_matrix *per_island_locus_recomb_rates;
};

#endif
