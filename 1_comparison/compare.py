import pandas as pd


############### Combine the indel and 0mut ###############
orig_indel = pd.read_csv('1_comparison/orig_indel_old.tsv', sep = '\t')
orig_indel['Num_Subst'] = 0
orig_0mut = pd.read_csv('1_comparison/orig_0mut_old.tsv', sep = '\t')
orig_0mut['CIGAR'] = orig_0mut['Sequence'].str.len().astype(str) + 'M'

orig_indel_0mut = pd.concat([orig_0mut, orig_indel], axis = 'index')
orig_indel_0mut = orig_indel_0mut.sort_values('Count', ascending = False).reset_index(drop = True)
orig_indel_0mut.to_csv('1_comparison/orig_indel_0mut_old.tsv', sep = '\t', index = False)
############################################################

########## Find the difference of the orig and new ##########
orig_indel_0mut = pd.read_csv('1_comparison/orig_indel_0mut_old.tsv', sep = '\t')
new_indel_0mut = pd.read_csv('1_comparison/new_indel_0mut_4.tsv', sep = '\t')

# orig_indel_0mut = orig_indel_0mut.set_index(['Sequence', 'CIGAR'])['Count']
# new_indel_0mut = new_indel_0mut.set_index(['Sequence', 'CIGAR'])['Count']
orig_indel_0mut = orig_indel_0mut.groupby('Sequence').agg(
  Count = ('Count', 'sum'),
  Num_Subst = ('Num_Subst', 'first'),
  CIGAR = ('CIGAR', 'first'),
)
new_indel_0mut = new_indel_0mut.groupby('Sequence').agg(
  Count = ('Count', 'sum'),
  Num_Subst = ('Num_Subst', 'first'),
  CIGAR = ('CIGAR', 'first'),
)

new_indel_0mut, orig_indel_0mut = new_indel_0mut.align(
  orig_indel_0mut,
  join = 'outer',
  axis = 'index',
)

new_indel_0mut = new_indel_0mut.reset_index()
orig_indel_0mut = orig_indel_0mut.reset_index()

new_indel_0mut.loc[new_indel_0mut['Count'].isna(), 'Count'] = 0
orig_indel_0mut.loc[orig_indel_0mut['Count'].isna(), 'Count'] = 0
new_indel_0mut.loc[new_indel_0mut['Num_Subst'].isna(), 'Num_Subst'] = orig_indel_0mut['Num_Subst']
orig_indel_0mut.loc[orig_indel_0mut['Num_Subst'].isna(), 'Num_Subst'] = new_indel_0mut['Num_Subst']
new_indel_0mut.loc[new_indel_0mut['CIGAR'].isna(), 'CIGAR'] = orig_indel_0mut['CIGAR']
orig_indel_0mut.loc[orig_indel_0mut['CIGAR'].isna(), 'CIGAR'] = new_indel_0mut['CIGAR']

diff_indel_0mut = pd.DataFrame({
  'Sequence': new_indel_0mut['Sequence'],
  'CIGAR': new_indel_0mut['CIGAR'],
  'Count_orig': orig_indel_0mut['Count'],
  'Count_new': new_indel_0mut['Count'],
  'Count_diff': new_indel_0mut['Count'] - orig_indel_0mut['Count'],
  'Num_Subst' : new_indel_0mut['Num_Subst'],
})
diff_indel_0mut = diff_indel_0mut.loc[diff_indel_0mut['Count_diff'] != 0]
diff_indel_0mut['Count_orig'] = diff_indel_0mut['Count_orig'].astype(int)
diff_indel_0mut['Count_new'] = diff_indel_0mut['Count_new'].astype(int)
diff_indel_0mut['Count_diff'] = diff_indel_0mut['Count_diff'].astype(int)
diff_indel_0mut['Num_Subst'] = diff_indel_0mut['Num_Subst'].astype(int)

diff_indel_0mut['Count_diff_abs'] = diff_indel_0mut['Count_diff'].abs()
diff_indel_0mut = diff_indel_0mut.sort_values('Count_diff_abs', ascending = False)
diff_indel_0mut = diff_indel_0mut.drop('Count_diff_abs', axis = 'columns')
diff_indel_0mut = diff_indel_0mut.drop('Num_Subst', axis = 'columns')
diff_indel_0mut.to_csv('1_comparison/diff_indel_0mut_old.tsv', sep = '\t', index = False)
