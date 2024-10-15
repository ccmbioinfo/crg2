#!/usr/bin/env python3

#
# Expansion Hunter Denovo
# Copyright (c) 2017 Illumina, Inc.
#
# Author: Egor Dolzhenko <edolzhenko@illumina.com>,
#         Sai Chen <schen6@illumina.com>
# Concept: Michael Eberle <meberle@illumina.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

import numpy as np
import scipy.special
import scipy.stats


def calculate_approximate_pvalue(cases, controls):
    all_counts = cases + controls
    ranks = scipy.stats.rankdata(all_counts)
    case_rank_sum = np.sum(ranks[:len(cases)])

    mu_cases = len(cases) * (len(cases) + len(controls) + 1) / 2
    sigma_cases = np.sqrt(len(cases) * len(controls) *
                          (len(cases) + len(controls) + 1) / 12)
    z_cases = (case_rank_sum - mu_cases) / sigma_cases

    return 1 - scipy.stats.norm.cdf(z_cases)


def calculate_permutation_pvalue(cases, controls, num_permutations):
    all_counts = cases + controls
    ranks = scipy.stats.rankdata(all_counts)
    num_cases = len(cases)
    true_case_rank_sum = np.sum(ranks[:num_cases])

    permuted_case_ranks = np.random.choice(
        ranks, size=(num_permutations, num_cases))
    permuted_case_rank_sums = np.sum(permuted_case_ranks, axis=1)

    num_case_rank_sums_as_extreme_as_true = np.sum(
        permuted_case_rank_sums >= true_case_rank_sum)

    return (num_case_rank_sums_as_extreme_as_true + 1) / (num_permutations + 1)


def wilcoxon_rank_sum_test(test_params, cases, controls):
    test_method = test_params['method']
    if test_method == 'normal':
        return calculate_approximate_pvalue(cases, controls)
    elif test_method == 'resample':
        return calculate_permutation_pvalue(cases, controls, test_params['num_resamples'])
    else:
        assert False, '{} is an unknown method type'.format(test_method)
