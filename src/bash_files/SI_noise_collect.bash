#!/usr/bin/bash

sinter collect\
    --processes 64\
    --circuits circuits/SI/*/*.stim\
    --metadata_func "sinter.comma_separated_key_values(path)"\
    --decoders pymatching\
    --max_shots 1_000_000\
    --quiet\
    --max_errors 100_000\
    --save_resume_filepath collected_stats/SI_stats_v2.csv;\
sinter collect\
    --processes 64\
    --circuits circuits/SI/*/*.stim\
    --metadata_func "sinter.comma_separated_key_values(path)"\
    --decoders pymatching\
    --max_shots 10_000_000\
    --quiet\
    --max_errors 1_000\
    --save_resume_filepath collected_stats/SI_stats_v2.csv;\
sinter collect\
    --processes 64\
    --circuits circuits/SI/*/*.stim\
    --metadata_func "sinter.comma_separated_key_values(path)"\
    --decoders pymatching\
    --max_shots 100_000_000\
    --quiet\
    --max_errors 1_000\
    --save_resume_filepath collected_stats/SI_stats_v2.csv;\
sinter collect\
    --processes 64\
    --circuits circuits/SI/*/*.stim\
    --metadata_func "sinter.comma_separated_key_values(path)"\
    --decoders pymatching\
    --max_shots 1_000_000_000\
    --quiet\
    --max_errors 500\
    --save_resume_filepath collected_stats/SI_stats_v2.csv;\
sinter collect\
    --processes 64\
    --circuits circuits/SI/*/*.stim\
    --metadata_func "sinter.comma_separated_key_values(path)"\
    --decoders pymatching\
    --quiet\
    --max_shots 10_000_000_000\
    --max_errors 100\
    --save_resume_filepath collected_stats/SI_stats_v2.csv;\
sinter collect\
    --processes 64\
    --circuits circuits/SI/*/*.stim\
    --metadata_func "sinter.comma_separated_key_values(path)"\
    --decoders pymatching\
    --quiet\
    --max_shots 10_000_000_000\
    --max_errors 200\
    --save_resume_filepath collected_stats/SI_stats_v2.csv;\
sinter collect\
    --processes 64\
    --circuits circuits/SI/*/*.stim\
    --metadata_func "sinter.comma_separated_key_values(path)"\
    --decoders pymatching\
    --quiet\
    --max_shots 20_000_000_000\
    --max_errors 200\
    --save_resume_filepath collected_stats/SI_stats_v2.csv;\
sinter collect\
    --processes 64\
    --circuits circuits/SI/*/*.stim\
    --metadata_func "sinter.comma_separated_key_values(path)"\
    --decoders pymatching\
    --quiet\
    --max_shots 30_000_000_000\
    --max_errors 100\
    --save_resume_filepath collected_stats/SI_stats_v2.csv;\
sinter collect\
    --processes 64\
    --circuits circuits/SI/*/*.stim\
    --metadata_func "sinter.comma_separated_key_values(path)"\
    --decoders pymatching\
    --quiet\
    --max_shots 40_000_000_000\
    --max_errors 100\
    --save_resume_filepath collected_stats/SI_stats_v2.csv;\
sinter collect\
    --processes 64\
    --circuits circuits/SI/*/*.stim\
    --metadata_func "sinter.comma_separated_key_values(path)"\
    --decoders pymatching\
    --quiet\
    --max_shots 50_000_000_000\
    --max_errors 50\
    --save_resume_filepath collected_stats/SI_stats_v2.csv;\
sinter collect\
    --processes 64\
    --circuits circuits/SI/*/*.stim\
    --metadata_func "sinter.comma_separated_key_values(path)"\
    --decoders pymatching\
    --quiet\
    --max_shots 50_000_000_000\
    --max_errors 100\
    --save_resume_filepath collected_stats/SI_stats_v2.csv;\
sinter collect\
    --processes 64\
    --quiet\
    --circuits circuits/SI/*/*.stim\
    --metadata_func "sinter.comma_separated_key_values(path)"\
    --decoders pymatching\
    --max_shots 70_000_000_000\
    --max_errors 100\
    --save_resume_filepath collected_stats/SI_stats_v2.csv;\
sinter collect\
    --processes 64\
    --quiet\
    --circuits circuits/SI/*/*.stim\
    --metadata_func "sinter.comma_separated_key_values(path)"\
    --decoders pymatching\
    --max_shots 100_000_000_000\
    --max_errors 100\
    --save_resume_filepath collected_stats/SI_stats_v2.csv;\
sinter collect\
    --processes 64\
    --quiet\
    --circuits circuits/SI/*/*.stim\
    --metadata_func "sinter.comma_separated_key_values(path)"\
    --decoders pymatching\
    --max_shots 250_000_000_000\
    --max_errors 100\
    --save_resume_filepath collected_stats/SI_stats_v2.csv;\
sinter collect\
    --processes 64\
    --quiet\
    --circuits circuits/SI/*/*.stim\
    --metadata_func "sinter.comma_separated_key_values(path)"\
    --decoders pymatching\
    --max_shots 300_000_000_000\
    --max_errors 100\
    --save_resume_filepath collected_stats/SI_stats_v2.csv;\
sinter collect\
    --processes 64\
    --quiet\
    --circuits circuits/SI/*/*.stim\
    --metadata_func "sinter.comma_separated_key_values(path)"\
    --decoders pymatching\
    --max_shots 500_000_000_000\
    --max_errors 50\
    --save_resume_filepath collected_stats/SI_stats_v2.csv;\
sinter collect\
    --processes 64\
    --circuits circuits/SI/*/*.stim\
    --metadata_func "sinter.comma_separated_key_values(path)"\
    --decoders pymatching\
    --quiet\
    --max_shots 1_000_000_000_000\
    --max_errors 50\
    --save_resume_filepath collected_stats/SI_stats_v2.csv;\
