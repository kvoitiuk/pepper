//
// Created by Kishwar Shafin on 05/03/19.
//

#ifndef HELEN_PILEUP_GENERATOR_H
#define HELEN_PILEUP_GENERATOR_H

#include <math.h>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <cstring>
#include <assert.h>
#include <string.h>
#include <stdio.h>
using namespace std;
#include "../dataio/bam_handler.h"


namespace ImageOptions {
    static constexpr uint8_t MAX_COLOR_VALUE = 254;
    static constexpr uint8_t BASE_QUALITY_CAP = 40;
    static constexpr uint8_t MAP_QUALITY_CAP = 60;
    static constexpr int IMAGE_HEIGHT = 100;
    static constexpr int IMAGE_CHANNELS = 3;
};

struct ImagePixel{
    uint8_t base_color;
    uint8_t base_quality_color;
    uint8_t map_quality_color;
    uint8_t strand_color;

    ImagePixel(char base, int mapping_quality, int base_quality, bool is_reverse) {
        map<char, uint8_t> global_base_color = {{'C', 50}, {'T', 100}, {'G', 150},
                {'A', 200}, {'*', 250}, {'.', 10}, {'N', 10}};
        base_color = global_base_color[base];

        base_quality_color = (double) ImageOptions::MAX_COLOR_VALUE *
                          ((double) min(base_quality, ImageOptions::BASE_QUALITY_CAP)
                           / (double) ImageOptions::BASE_QUALITY_CAP);


        map_quality_color = (double) ImageOptions::MAX_COLOR_VALUE *
                         ((double) min(mapping_quality, ImageOptions::MAP_QUALITY_CAP) /
                          (double) ImageOptions::MAP_QUALITY_CAP);

        strand_color = is_reverse ? ImageOptions::MAX_COLOR_VALUE : 70;
    }
    ImagePixel() {
        base_color = 0;
        base_quality_color = 0;
        map_quality_color = 0;
        strand_color = 0;
    }
    void operator=(const ImagePixel& that) {
        this->base_color = that.base_color;
        this->base_quality_color = that.base_quality_color;
        this->map_quality_color = that.map_quality_color;
        this->strand_color = that.strand_color;
    }
};

struct PositionMap {
    long long ref_pos;
    int insert_pos;

    PositionMap(long long rp, int ip) {
        ref_pos = rp;
        insert_pos = ip;
    }
};

struct ReadPositionMap {
    string read_id;
    long long ref_pos;
    int insert_pos;

    ReadPositionMap(string rid, long long rp, int ip) {
        ref_pos = rp;
        insert_pos = ip;
        read_id = rid;
    }
};

struct ReadPositionMapComparator {
    bool operator() ( ReadPositionMap a, ReadPositionMap b ) const {
        if(a.read_id == b.read_id) {
            if (a.ref_pos == b.ref_pos) {
                return a.insert_pos < b.insert_pos;
            } else {
                return a.ref_pos < b.ref_pos;
            }
        } else {
            return a.read_id < b.read_id;
        }
    }
};

struct PositionMapComparator {
    bool operator() ( PositionMap a, PositionMap b ) const {
        if (a.ref_pos == b.ref_pos) {
            return a.insert_pos < b.insert_pos;
        } else {
            return a.ref_pos < b.ref_pos;
        }
    }
};


class PileupGenerator {
    long long ref_start;
    long long ref_end;
    string chromosome_name;
    string reference_sequence;
    bool train_mode;

    map<ReadPositionMap, ImagePixel, ReadPositionMapComparator> pixel_summaries;
    map<PositionMap, int, PositionMapComparator> position_to_image_index;
    map<string, pair<long long, long long> > read_start_end_map;
    map<long long, int> longest_insert_count;
    map<long long, double> coverage;

    map< pair<long long, int>, char> insert_labels;
    map< long long, char> base_labels;
    map<int, long long> row_end_map;
public:
    vector< vector< vector<uint8_t> > > image;
    vector<uint8_t> labels;
    vector<pair<long long, int> > genomic_pos;
    vector<int> bad_label_positions;

    PileupGenerator(string reference_sequence,
                     string chromosome_name,
                     long long ref_start,
                     long long ref_end,
                     bool train_mode);
    void generate_summary(vector <type_read> &reads,
                          long long start_pos,
                          long long end_pos);

    void iterate_over_read(type_read read, long long region_start, long long region_end);
    void generate_labels(type_read truth_reads, long long region_start, long long region_end);
    void generate_ref_features();
    void debug_print();
    int get_image_row(long long start_pos, long long end_pos);
    void generate_image(vector<string> read_id_list,
                        int image_length,
                        long long region_start,
                        long long region_end);
};


#endif //HELEN_PILEUP_GENERATOR_H
