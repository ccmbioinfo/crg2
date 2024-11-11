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
        description='combine counts of in-repeat reads across multiple samples')
    parser.add_argument(
        '--manifest',
        help='TSV file with id, case/control status, and path of each sample',
        required=True)
    parser.add_argument(
        '--combinedCounts',
        help='file with combined counts of anchored in-repeat reads',
        required=True)
    parser.add_argument(
        '--minUnitSize',
        help='size of smallest repeat unit to process (in bp)',
        default=2, type=int, required=False)
    parser.add_argument(
        '--maxUnitSize',
        help='size of the largest repeat unit to process (in bp)',
        default=20, type=int, required=False)
    args = parser.parse_args()

    parameter_encoding = ' '.join(sys.argv[1:])
    logging.info('Starting with these parameters: %s', parameter_encoding)

    return {'manifest_path': args.manifest, 'output_path': args.combinedCounts,
            'min_unit_size': args.minUnitSize,
            'max_unit_size': args.maxUnitSize}


def merge_combined_regions(combined_air_regions):
    for unit in combined_air_regions:
        combined_air_regions[unit].merge()


def process_sample(parameters, denovo_json, sample_id, combined_air_regions,
                   combined_pir_regions):
    min_unit_size = parameters['min_unit_size']
    max_unit_size = parameters['max_unit_size']

    for unit in denovo_json:
        if not min_unit_size <= len(unit) <= max_unit_size:
            continue

        if unit == 'Depth' or unit == 'ReadLength':
            continue

        if 'RegionsWithIrrAnchors' in denovo_json[unit]:
            anc_irr_record = denovo_json[unit]['RegionsWithIrrAnchors']
            current_regions = regiontools.create_region_collection_from_denovo_record(
                sample_id, anc_irr_record)

            if unit not in combined_air_regions:
                combined_air_regions[unit] = regiontools.RegionCollection()
            combined_air_regions[unit].extend(current_regions)

        if denovo_json[unit]['IrrPairCount']:
            num_irr_pairs = denovo_json[unit]['IrrPairCount']
            if unit not in combined_pir_regions:
                combined_pir_regions[unit] = {}
            combined_pir_regions[unit][sample_id] = num_irr_pairs


def process_samples(parameters, samples, sample_depths):
    combined_air_regions = {}
    combined_pir_regions = {}

    for sample_num, sample in enumerate(samples):
        logging.info('Processing sample %s', sample.id)
        with open(sample.path, 'r') as json_file:
            denovo_json = json.load(json_file)
            sample_depths[sample.id] = denovo_json['Depth']
            process_sample(parameters, denovo_json, sample.id,
                           combined_air_regions,
                           combined_pir_regions)

        if sample_num > 0 and sample_num % 50 == 0:
            logging.info('Merging regions')
            merge_combined_regions(combined_air_regions)

    logging.info('Performing final merge')
    merge_combined_regions(combined_air_regions)

    return {'air': combined_air_regions, 'pir': combined_pir_regions}


def normalize_count(sample_depth, count):
    target_depth = 40
    return target_depth * count / sample_depth


def normalize_counts(sample_depths, combined_regions):
    for unit, regions in combined_regions['air'].items():
        for region in regions:
            for sample_id, count in region.feature_counts._count_dict.items():
                sample_depth = sample_depths[sample_id]
                region.feature_counts._count_dict[sample_id] = normalize_count(
                    sample_depth, count)
    for unit, sample_counts in combined_regions['pir'].items():
        for sample_id, count in sample_counts.items():
            sample_depth = sample_depths[sample_id]
            combined_regions['pir'][unit][sample_id] = normalize_count(
                sample_depth, count)


def create_json(combined_regions):
    output_json = {}
    for unit, regions in combined_regions['air'].items():
        output_json[unit] = {}
        output_json[unit]['RegionsWithIrrAnchors'] = regions.as_dict()

    for unit, counts in combined_regions['pir'].items():
        if unit not in output_json:
            output_json[unit] = {}
        output_json[unit]['IrrPairCounts'] = counts

    return output_json


def write_json(parameters, output_json):
    with open(parameters['output_path'], 'w') as output_file:
        json_encoding = json.dumps(output_json, sort_keys=True,
                                   indent=4, separators=(',', ': '))
        print(json_encoding, file=output_file)


def main():
    common.init_logger()
    parameters = load_parameters()
    samples = common.load_manifest(parameters['manifest_path'])
    sample_depths = {}
    combined_regions = process_samples(parameters, samples, sample_depths)
    normalize_counts(sample_depths, combined_regions)
    output_json = create_json(combined_regions)
    write_json(parameters, output_json)


if __name__ == '__main__':
    main()
