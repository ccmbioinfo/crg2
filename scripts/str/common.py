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

import collections
import logging
import json
import scipy.stats as stats

from . import regiontools

from .wilcoxontest import wilcoxon_rank_sum_test


def init_logger():
    logging.basicConfig(format='%(asctime)s: %(message)s', level=logging.INFO)


def load_manifest(manifest_path):
    '''Extract sample information from a manifest file.

    '''

    # pylint: disable=I0011,C0103
    Sample = collections.namedtuple('Sample', 'id status path')

    samples = []
    with open(manifest_path, 'r') as manifest_file:
        for line in manifest_file:
            sample_id, status, path = line.split()
            sample = Sample(id=sample_id, status=status, path=path)
            samples.append(sample)
    return samples


def load_combined_json(json_path):
    logging.info('Loading %s', json_path)
    with open(json_path, 'r') as json_file:
        combined_json = json.load(json_file)
    return combined_json


def filter_counts_by_magnitude(count_table, count_cutoff):
    filtered_count_table = []
    for row in count_table:
        max_count = max(count for _, count in row['sample_counts'].items())
        if max_count >= count_cutoff:
            filtered_count_table.append(row)

    return filtered_count_table


def filter_counts_by_region(count_table, target_regions):
    filtered_count_table = []
    for row in count_table:
        chrom, start, end = row['region'].replace(':', '-').split('-')
        start, end = int(start), int(end)
        region = regiontools.Region(chrom, start, end)
        overlaps_target_region = any(regiontools.compute_distance(
            region, target) == 0 for target in target_regions)

        if overlaps_target_region:
            filtered_count_table.append(row)

    return filtered_count_table


def filter_counts(count_table, count_cutoff, target_regions=None):
    filtered_count_table = filter_counts_by_magnitude(
        count_table, count_cutoff)
    if target_regions:
        filtered_count_table = filter_counts_by_region(
            filtered_count_table, target_regions)
    return filtered_count_table


def extract_case_control_assignments(samples):
    sample_status = {}
    for sample in samples:
        sample_status[sample.id] = sample.status
    return sample_status


def test_samples(test_params, sample_status, sample_counts):
    control_samples = [sample for sample,
                       status in sample_status.items() if status == 'control']
    case_samples = [sample for sample,
                    status in sample_status.items() if status == 'case']

    control_counts = [sample_counts[s]
                      if s in sample_counts else 0 for s in control_samples]

    case_counts = [sample_counts[s]
                   if s in sample_counts else 0 for s in case_samples]

    pvalue = wilcoxon_rank_sum_test(test_params, case_counts, control_counts)

    return pvalue


def compare_counts(test_params, sample_status, count_table):
    for row in count_table:
        # Generate counts before testing
        pvalue = test_samples(test_params, sample_status, row['sample_counts'])
        row['pvalue'] = pvalue


def correct_pvalues(count_table):
    num_tests = len(count_table)
    for row in count_table:
        row['bonf_pvalue'] = min(row['pvalue'] * num_tests, 1.0)
