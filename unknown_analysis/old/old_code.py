# collected_dir = os.path.join(root_dir, 'output', 'collected')
# ranks_dir = os.path.join(root_dir, 'output', 'ranks')


# # Make the collect_fastq script
# for ext in ['.sh', '.ps1']:
#   sep = '\\' if ext == '.ps1' else '/'
#   with open('run_collect_fastq' + ext, 'w') as out:
#     for info in library_info:
#       i = join_path(sep, fastq_dir, get_lib_name(info) + '.fastq')
#       o = join_path(sep, collected_dir, get_lib_name(info) + '.tsv')
#       out.write(f'python collect_fastq.py -i {i} -o {o}\n')
# # yjl090_WT_sgCD_R1_antisense.tsv

# # Make the get_ranks script
# for ext in ['.sh', '.ps1']:
#   sep = '\\' if ext == '.ps1' else '/'
#   with open('run_get_ranks' + ext, 'w') as out:
#     for info in library_info:
#       i = join_path(sep, args.n, get_lib_name_long(info) + '.tsv')
#       r = join_path(sep, collected_dir, get_lib_name(info) + '.tsv')
#       o = join_path(sep, ranks_dir, get_lib_name(info) + '_ranks.tsv')
#       out.write(f'python get_ranks.py -i {i} -r {r} -o {o}\n')
