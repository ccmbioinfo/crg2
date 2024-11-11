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

import sys
import copy


class FeatureCounts(object):
    def __init__(self, count_dict):
        self._count_dict = count_dict

    def __len__(self):
        return len(self._count_dict)

    def __getitem__(self, sample_id):
        return self._count_dict[sample_id]

    def __repr__(self):
        encoding = ('{}={}'.format(s, self._count_dict[s])
                    for s in sorted(self._count_dict))
        encoding = ','.join(encoding)
        return encoding

    def __iter__(self):
        for sample, count in self._count_dict.items():
            yield (sample, count)

    def __eq__(self, other):
        # pylint: disable=I0011,W0212
        return self._count_dict == other._count_dict

    def combine(self, other):
        # pylint: disable=I0011,W0212
        for sample, count in other._count_dict.items():
            if sample in self._count_dict:
                self._count_dict[sample] += count
            else:
                self._count_dict[sample] = count


class Region(object):
    def __init__(self, chrom=None, start=0, end=0, feature_counts=None):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.feature_counts = feature_counts

    def is_unset(self):
        return self.chrom is None

    def __repr__(self):
        region_repr = '{}:{}-{}'.format(self.chrom, self.start, self.end)
        if self.feature_counts:
            region_repr += '({})'.format(self.feature_counts)
        return region_repr

    def __eq__(self, other):
        return (self.chrom == other.chrom and self.start == other.start and
                self.end == other.end and
                self.feature_counts == other.feature_counts)

    def __lt__(self, other):
        self_tuple = (self.chrom, self.start, self.end)
        other_tuple = (other.chrom, other.start, other.end)
        return self_tuple < other_tuple


def compute_distance(region_a, region_b):
    distance = sys.maxsize
    if region_a.is_unset() or region_b.is_unset():
        raise Exception('Cannot compute distance between unset regions')

    if region_a.chrom == region_b.chrom:
        if region_a.end < region_b.start:  # Is region_a left of region_b?
            distance = region_b.start - region_a.end
        elif region_b.end < region_a.start:  # Is region_b left or region_a?
            return region_a.start - region_b.end
        else:  # If neither is true, the regions must overlap.
            distance = 0

    return distance


class RegionCollection(object):
    def __init__(self, regions=None):
        self._regions = regions
        if self._regions is None:
            self._regions = []

    def __len__(self):
        return len(self._regions)

    def __getitem__(self, index):
        return self._regions[index]

    def __repr__(self):
        encoding = '\n'.join(repr(r) for r in self._regions)
        return encoding

    def __eq__(self, other):
        # pylint: disable=I0011,W0212
        for region_a, region_b in zip(self._regions, other._regions):
            if region_a != region_b:
                return False
        return True

    def __iter__(self):
        for region in self._regions:
            yield region

    def as_dict(self):
        region_dict = {}
        for region in self._regions:
            region_encoding = '{}:{}-{}'.format(region.chrom,
                                                region.start, region.end)
            region_dict[region_encoding] = dict(region.feature_counts)
        return region_dict

    def merge(self, max_dist=500):
        self._regions.sort()

        merged_regions = []
        aggregate_region = Region()

        for region in self._regions:
            if aggregate_region.is_unset():
                aggregate_region = region
                continue

            if compute_distance(aggregate_region, region) <= max_dist:
                aggregate_region.end = max(aggregate_region.end, region.end)
                aggregate_region.feature_counts.combine(region.feature_counts)
            else:
                merged_regions.append(aggregate_region)
                aggregate_region = region

        if not aggregate_region.is_unset():
            merged_regions.append(aggregate_region)

        self._regions = merged_regions

    def extend(self, other):
        # pylint: disable=I0011,W0212
        self._regions += other._regions


def create_region_from_denovo_record(sample_id, region_count):
    region_encoding, count = region_count
    region_encoding = region_encoding.replace(':', '-').split('-')
    chrom, start, end = region_encoding
    start, end = int(start), int(end)
    feture_counts = FeatureCounts({sample_id: count})
    return Region(chrom, start, end, feture_counts)


def create_region_collection_from_denovo_record(sample_id, denovo_record):
    regions = []
    for region_count in denovo_record.items():
        region = create_region_from_denovo_record(sample_id, region_count)
        regions.append(region)
    return RegionCollection(regions)
