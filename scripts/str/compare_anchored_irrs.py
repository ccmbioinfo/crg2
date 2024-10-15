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

import argparse
import logging
import sys
import json

from core import regiontools, common


def load_parameters():
    '''Capture command-line parameters.

    '''
    parser = argparse.ArgumentParser(
        description='compare counts of anchored in-repeat reads')

    # Add required arguments.
    manifest_help = 'TSV file with id, case/control status, and path of each sample'
    input_help = 'JSON file with combined counts of anchored in-repeat reads'
    output_help = ('BED file with regions containing anchored in-repeat reads and associated'
                   ' p-values')
    required_arg_group = parser.add_argument_group('required arguments')
    required_arg_group.add_argument(
        '--manifest', help=manifest_help, required=True)
    required_arg_group.add_argument(
        '--inputCounts', help=input_help, required=True)
    required_arg_group.add_argument(
        '--outputRegions', help=output_help, required=True)

    # Add optional arguments.
    def_min_count = 5
    min_count_help = ('minimum number reads in a region for downstream analysis'
                      ' (default: {})').format(def_min_count)
    target_regions_help = 'BED file with regions to which analysis should be restricted'
    def_test_method = 'normal'
    test_method_help = 'Method of calculating Wilcoxon Rank-Sum Test p-value (default: {})'.format(
        def_test_method)
    def_num_resamples = 1000000
    num_resamples_help = 'Number of iterations for the resampling test (default: {})'.format(
        def_num_resamples)
    parser.add_argument(
        '--minCount', help=min_count_help, default=def_min_count, type=int)
    parser.add_argument('--targetRegions',
                        help=target_regions_help, default=None)
    parser.add_argument('--testMethod', help=test_method_help,
                        choices=['normal', 'resample'], default=def_test_method)
    parser.add_argument(
        '--numResamples', help=num_resamples_help, default=def_num_resamples, type=int)

    args = parser.parse_args()

    parameter_encoding = ' '.join(sys.argv[1:])
    logging.info('Starting with these parameters: %s', parameter_encoding)

    test_params = {'method': args.testMethod,
                   'num_resamples': args.numResamples}
    return {'manifest_path': args.manifest, 'counts_path': args.inputCounts,
            'min_count': args.minCount, 'output_path': args.outputRegions,
            'target_regions': args.targetRegions, 'test_params': test_params}


def generate_count_table(combined_counts):
    count_table = []
    for unit, rec in combined_counts.items():
        if 'RegionsWithIrrAnchors' not in rec:
            continue

        for region, sample_counts in rec['RegionsWithIrrAnchors'].items():
            table_row = {'region': region, 'unit': unit}
            table_row['sample_counts'] = sample_counts
            count_table.append(table_row)

    logging.info('Loaded %i regions', len(count_table))
    return count_table


def load_target_regions(fname):
    logging.info('Loading target regions from %s', fname)
    regions = []
    with open(fname, 'r') as bed_file:
        for line in bed_file:
            chrom, start, end, *_ = line.split()
            start, end = int(start), int(end)
            region = regiontools.Region(chrom, start, end)
            regions.append(region)
    return regions


def output_results(count_table, output_path):
    with open(output_path, 'w') as output_file:
        for row in count_table:
            chrom, start, end = row['region'].replace(':', '-').split('-')
            unit = row['unit']
            pvalue, bonf_pvalue = row['pvalue'], row['bonf_pvalue']

            sample_counts = row['sample_counts']
            encoded_counts = ['{}:{}'.format(s, c)
                              for s, c in sample_counts.items()]
            encoded_counts = ','.join(encoded_counts)
            print(chrom, start, end, unit, pvalue, bonf_pvalue, encoded_counts,
                  sep='\t', file=output_file)


def main():
    common.init_logger()
    parameters = load_parameters()
    combined_counts = common.load_combined_json(parameters['counts_path'])
    samples = common.load_manifest(parameters['manifest_path'])
    target_regions = None
    if parameters['target_regions']:
        target_regions = load_target_regions(parameters['target_regions'])
    count_table = generate_count_table(combined_counts)
    count_table = common.filter_counts(
        count_table, parameters['min_count'], target_regions)
    logging.info('%i regions left after initial filtering', len(count_table))
    sample_status = common.extract_case_control_assignments(samples)
    common.compare_counts(
        parameters['test_params'], sample_status, count_table)
    common.correct_pvalues(count_table)
    output_results(count_table, parameters['output_path'])
    logging.info('Done')


if __name__ == '__main__':
    main()
