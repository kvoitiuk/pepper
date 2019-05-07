//
// Created by Kishwar Shafin on 11/1/18.
//

#include "../../headers/pileup_summary/pileup_generator.h"

PileupGenerator::PileupGenerator(string reference_sequence, string chromosome_name, long long ref_start,
                                 long long ref_end, bool train_mode) {
    this->reference_sequence = reference_sequence;
    this->ref_start = ref_start;
    this->ref_end = ref_end;
    this->chromosome_name = chromosome_name;
    this->train_mode = train_mode;
}


int get_feature_index(char base, bool is_reverse) {
    base = toupper(base);
    if (is_reverse) {
        if (base == 'A') return 0;
        if (base == 'C') return 1;
        if (base == 'G') return 2;
        if (base == 'T') return 3;
        return 8;
    } else {
        // tagged and forward
        if (base == 'A') return 4;
        if (base == 'C') return 5;
        if (base == 'G') return 6;
        if (base == 'T') return 7;
        return 9;
    }
}


uint8_t get_labels(char base) {
    base = toupper(base);
    if (base == 'A') return 1;
    if (base == 'C') return 2;
    if (base == 'G') return 3;
    if (base == 'T') return 4;
    if (base == '*') return 0; // this is for deleted bases, but the number is so small that it creates confusion
    if (base == '#') return 0;
    return 0;
}


void PileupGenerator::iterate_over_read(type_read read, long long region_start, long long region_end) {
    int read_index = 0;
    long long ref_position = read.pos;
    int cigar_index = 0;
    int base_quality = 0;
    long long reference_index;
    int mapping_quality = read.mapping_quality;
    string read_id = read.read_id;
    long long read_start = -1;
    long long read_end = -1;

    for (auto &cigar: read.cigar_tuples) {
        if (ref_position > region_end) break;
        switch (cigar.operation) {
            case CIGAR_OPERATIONS::EQUAL:
            case CIGAR_OPERATIONS::DIFF:
            case CIGAR_OPERATIONS::MATCH:
                cigar_index = 0;
                if (ref_position < ref_start) {
                    cigar_index = min(ref_start - ref_position, (long long) cigar.length);
                    read_index += cigar_index;
                    ref_position += cigar_index;
                }
                for (int i = cigar_index; i < cigar.length; i++) {
                    reference_index = ref_position - ref_start;
                    //read.base_qualities[read_index] base quality
                    if (ref_position >= ref_start && ref_position <= ref_end) {
                        if(read_start == -1) {
                            read_start = ref_position;
                            read_end = ref_position;
                        }
                        read_start = min(read_start, ref_position);
                        read_end = max(read_end, ref_position);

                        char base = read.sequence[read_index];
                        int base_quality = read.base_qualities[read_index];

                        // create a pixel map for this read
                        ReadPositionMap read_position(read_id, ref_position, 0);
                        ImagePixel pixel_value(base, base_quality, mapping_quality, read.flags.is_reverse);
                        pixel_summaries[read_position] = pixel_value;
                        coverage[ref_position] += 1.0;
                    }
                    read_index += 1;
                    ref_position += 1;
                }
                break;
            case CIGAR_OPERATIONS::IN:
//                base_qualities = read.base_qualities.begin() + read_index, read.base_qualities.begin() + (read_index + cigar.length);
                reference_index = ref_position - ref_start - 1;

                if (ref_position - 1 >= ref_start &&
                    ref_position - 1 <= ref_end) {

                    if(read_start == -1) {
                        read_start = ref_position;
                        read_end = ref_position;
                    }
                    read_start = min(read_start, ref_position);
                    read_end = max(read_end, ref_position);
                    // process insert allele here
                    string alt;
                    alt = read.sequence.substr(read_index, cigar.length);
                    for (int i = 0; i < cigar.length; i++) {
                        // create a pixel map for this read and the insert position
                        char base = alt[i];
                        int base_quality = read.base_qualities[read_index + i];
                        ReadPositionMap read_position(read_id, ref_position - 1, i + 1);
                        ImagePixel pixel_value(base, base_quality, mapping_quality, read.flags.is_reverse);
                        pixel_summaries[read_position] = pixel_value;
                    }
                    longest_insert_count[ref_position - 1] = std::max(longest_insert_count[ref_position - 1],
                                                                      (int) alt.length());
                }
                read_index += cigar.length;
                break;
            case CIGAR_OPERATIONS::REF_SKIP:
            case CIGAR_OPERATIONS::PAD:
            case CIGAR_OPERATIONS::DEL:
//                base_quality = read.base_qualities[max(0, read_index)];
                reference_index = ref_position - ref_start - 1;
                if(read_start == -1) {
                    read_start = ref_position;
                    read_end = ref_position;
                }
                // process delete allele here
                for (int i = 0; i < cigar.length; i++) {
                    if (ref_position + i >= ref_start && ref_position + i <= ref_end) {
                        read_start = min(read_start, ref_position + i);
                        read_end = max(read_end, ref_position + i);

                        // update the summary of base
                        char base = '*';
                        int base_quality = 15;
                        ReadPositionMap read_position(read_id, ref_position + i, 0);
                        ImagePixel pixel_value(base, base_quality, mapping_quality, read.flags.is_reverse);
                        pixel_summaries[read_position] = pixel_value;
                    }
                }
                ref_position += cigar.length;
                break;
            case CIGAR_OPERATIONS::SOFT_CLIP:
                read_index += cigar.length;
                break;
            case CIGAR_OPERATIONS::HARD_CLIP:
                break;
        }
    }
    read_start_end_map[read.read_id] = make_pair(read_start, read_end);
}

bool check_base(char base) {
    if(base=='A' || base=='a' ||
       base=='C' || base=='c' ||
       base=='T' || base=='t' ||
       base =='G' || base=='g' || base == '*' || base == '#') return true;
    return false;
}

void PileupGenerator::generate_labels(type_read read, long long region_start, long long region_end) {
    int read_index = 0;
    long long ref_position = read.pos;
    int cigar_index = 0;
    int base_quality = 0;
    long long reference_index;

    for (auto &cigar: read.cigar_tuples) {
        if (ref_position > region_end) break;
        switch (cigar.operation) {
            case CIGAR_OPERATIONS::EQUAL:
            case CIGAR_OPERATIONS::DIFF:
            case CIGAR_OPERATIONS::MATCH:
                cigar_index = 0;
                if (ref_position < ref_start) {
                    cigar_index = min(ref_start - ref_position, (long long) cigar.length);
                    read_index += cigar_index;
                    ref_position += cigar_index;
                }
                for (int i = cigar_index; i < cigar.length; i++) {
                    reference_index = ref_position - ref_start;
//                    cout<<ref_position<<" "<<ref_end<<" "<<region_end<<endl;
                    if (ref_position >= ref_start && ref_position <= ref_end) {
                        char base = read.sequence[read_index];
                        base_labels[ref_position] = base;
                    }
                    read_index += 1;
                    ref_position += 1;
                }
                break;
            case CIGAR_OPERATIONS::IN:
                reference_index = ref_position - ref_start - 1;
                if (ref_position - 1 >= ref_start &&
                    ref_position - 1 <= ref_end) {
                    // process insert allele here
                    string alt;
                    alt = read.sequence.substr(read_index, cigar.length);

                    for (int i = 0; i < longest_insert_count[ref_position - 1]; i++) {
                        char base = '#';
                        if (i < alt.length()) {
                            base = alt[i];
                        }
                        pair<long long, int> position_pair = make_pair(ref_position - 1, i);
                        insert_labels[position_pair] = base;
                    }
                }
                read_index += cigar.length;
                break;
            case CIGAR_OPERATIONS::REF_SKIP:
            case CIGAR_OPERATIONS::PAD:
            case CIGAR_OPERATIONS::DEL:
                reference_index = ref_position - ref_start - 1;

                if (ref_position >= ref_start && ref_position <= ref_end) {
                    // process delete allele here
                    for (int i = 0; i < cigar.length; i++) {
                        if (ref_position + i >= ref_start && ref_position + i <= ref_end) {
                            // DELETE
                            char base = '*';
                            base_labels[ref_position + i] = '*';
                        }
                    }
                }
                ref_position += cigar.length;
                break;
            case CIGAR_OPERATIONS::SOFT_CLIP:
                read_index += cigar.length;
                break;
            case CIGAR_OPERATIONS::HARD_CLIP:
                break;
        }
    }

    for (long long pos = region_start; pos <= region_end; pos++) {
        if(coverage[pos] > 0) {
            labels.push_back(get_labels(base_labels[pos]));
        } else {
            labels.push_back(get_labels('*'));
        }

        // if the label contains anything but ACTG
        if(!check_base(base_labels[pos])) {
            bad_label_positions.push_back(labels.size());
        }

        genomic_pos.push_back(make_pair(pos, 0));
        if (longest_insert_count[pos] > 0) {
            for (int ii = 0; ii < longest_insert_count[pos]; ii++) {
                if (insert_labels[make_pair(pos, ii)]) {
                    labels.push_back(get_labels(insert_labels[make_pair(pos, ii)]));

                    // if the label contains anything but ACTG
                    if(!check_base(insert_labels[make_pair(pos, ii)])) {
                        bad_label_positions.push_back(labels.size());
                    }
                }
                else labels.push_back(get_labels('#'));
            }
        }
    }
    bad_label_positions.push_back(labels.size());
}


void PileupGenerator::debug_print() {

    if(train_mode) {
        for(int i=0;i<labels.size();i++) {
            if((int)labels[i] == 0) cout<<"*";
            if((int)labels[i] == 1) cout<<"A";
            if((int)labels[i] == 2) cout<<"C";
            if((int)labels[i] == 3) cout<<"G";
            if((int)labels[i] == 4) cout<<"T";
        }
        cout<<endl;
    }


    map<uint8_t, char> global_base_decoder = {{50, 'C'}, {100, 'T'}, {150, 'G'},
                                              {200, 'A'}, {250, '*'}, {10, '.'}, {10, 'N'}};
    for(int i=0; i<image.size(); i++) {
        for(int j=0; j<image[i].size(); j++) {
            int base_color = image[i][j][0];
            if(base_color == 0)cout<<"-";
            else cout<<global_base_decoder[base_color];
        }
        cout<<endl;
    }

    for(int i=0; i<image.size(); i++) {
        for(int j=0; j<image[i].size(); j++) {
            int strand_color = image[i][j][2];
            if(strand_color == 0)cout<<"-";
            if(strand_color == 70)cout<<"0";
            if(strand_color == 254)cout<<"1";
        }
        cout<<endl;
    }
}

int PileupGenerator::get_image_row(long long start_pos, long long end_pos) {
    for(int i=0; i < ImageOptions::IMAGE_HEIGHT; i++) {
        if (row_end_map[i] < start_pos) {
            row_end_map[i] = end_pos;
            return i;
        }
    }
    return -1;
}

void PileupGenerator::generate_image(vector<string> read_id_list, int image_length, long long region_start,
                                     long long region_end) {
    vector< vector< vector<uint8_t> > > image_matrix(ImageOptions::IMAGE_HEIGHT,
                                                     vector<vector<uint8_t> >(image_length,
                                                     vector<uint8_t>(ImageOptions::IMAGE_CHANNELS)));

    for (auto &read_id:read_id_list) {
        // this probably means the read don't belong to this region, but this needs to be debugged and confirmed
        if(read_start_end_map[read_id].first < 0 || read_start_end_map[read_id].second < 0) continue;

        long long read_start = max(read_start_end_map[read_id].first, region_start);
        long long read_end = min(read_start_end_map[read_id].second, region_end);

        int image_row = get_image_row(read_start, read_end);
        if(image_row == -1) {
            continue;
        }

        for(int base_pos = read_start; base_pos<= read_end; base_pos ++) {
            // first the base position
            ReadPositionMap read_position_map(read_id, base_pos, 0);
            PositionMap position_map(base_pos, 0);
            ImagePixel pixel_values = pixel_summaries[read_position_map];
            int image_column_index = position_to_image_index[position_map];

            image_matrix[image_row][image_column_index][0] = pixel_values.base_color;
            image_matrix[image_row][image_column_index][1] = pixel_values.base_quality_color;
            image_matrix[image_row][image_column_index][2] = pixel_values.strand_color;

            for(int ins_pos = 0; ins_pos < longest_insert_count[base_pos]; ins_pos++) {
                ReadPositionMap read_position_map(read_id, base_pos, ins_pos + 1);
                PositionMap position_map(base_pos, ins_pos + 1);
                int image_column_index = position_to_image_index[position_map];

                if(pixel_summaries.find(read_position_map) != pixel_summaries.end()) {
                    ImagePixel pixel_values = pixel_summaries[read_position_map];
                    image_matrix[image_row][image_column_index][0] = pixel_values.base_color;
                    image_matrix[image_row][image_column_index][1] = pixel_values.base_quality_color;
                    image_matrix[image_row][image_column_index][2] = pixel_values.strand_color;
                } else {
                    image_matrix[image_row][image_column_index][0] = 10;
                    image_matrix[image_row][image_column_index][1] = 0;
                    image_matrix[image_row][image_column_index][2] = 0;
                }
            }
        }
//        cout<<read_start<<" "<<read_end<<" "<<image_row<<endl;
    }
    image = image_matrix;
}


bool ReadCompareFunction ( type_read a, type_read b ) {
    if(a.pos == b.pos) {
        return a.pos_end < b.pos_end;
    } else {
        return a.pos < b.pos;
    }
}


void PileupGenerator::generate_summary(vector <type_read> &reads,
                                        long long start_pos,
                                        long long end_pos) {

//    sort(reads.begin(), reads.end(), ReadCompareFunction);
    vector<string> read_id_list;
    for (auto &read:reads) {
        // this populates base_summaries and insert_summaries dictionaries
        if(read.mapping_quality > 0) {
            iterate_over_read(read, start_pos, end_pos);
            read_id_list.push_back(read.read_id);
        }
    }

    int image_length = 0;
    // after all the dictionaries are populated, we can simply walk through the region and generate a sequence
    for (long long pos = start_pos; pos <= end_pos; pos++) {
        genomic_pos.push_back(make_pair(pos, 0));
        PositionMap position(pos, 0);
        position_to_image_index[position] = image_length;
        image_length += 1;
        if (longest_insert_count[pos] > 0) {
            for (int ii = 0; ii < longest_insert_count[pos]; ii++) {
                genomic_pos.push_back(make_pair(pos, ii + 1));
                PositionMap position(pos, ii + 1);
                position_to_image_index[position] = image_length;
                image_length += 1;
            }
        }
    }

    generate_image(read_id_list, image_length, start_pos, end_pos);
//     at this point everything should be generated
//    debug_print();
}