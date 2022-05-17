x = pd.read_csv(r"C:\Users\tchan\Code\SCMB_Project\DistanceGraph\files_input\1DSB_WT_R1_sgA_sense.tsv", sep="\t")
y = pd.read_csv(r"C:\Users\tchan\Code\SCMB_Project\DistanceGraph\files_input\1DSB_WT_R1_sgA_sense_0mut.tsv", sep="\t")
y['CIGAR'] = y['Sequence'].str.len().astype(str) + 'M'
x = pd.concat([x, y], axis='index')
x = x[['Sequence', 'CIGAR', 'yjl255_freq']]
x = x.sort_values('yjl255_freq', ascending=False)
x['Sequence'] = x['Sequence'].str.replace('-', '')
x.to_csv('1_tj_edit/yjl255_orig.tsv', sep='\t', index=False)

x[['CIGAR']].to_csv('1_tj_edit/cigar_orig.tsv', sep='\t', index=False)
x[['Sequence']].to_csv('1_tj_edit/sequence_orig.tsv', sep='\t', index=False)
shape0 = x.shape[0]

x = pd.read_csv('1_tj_edit/yjl255_new.tsv', sep='\t')
x[['CIGAR']][0:shape0].to_csv('1_tj_edit/cigar_new.tsv', sep='\t', index=False)
x[['Sequence']][0:shape0].to_csv('1_tj_edit/sequence_new.tsv', sep='\t', index=False)