import cigar_utils
import re

def get_alignment(ref_seq, read_seq, ref_pos, cigar):
  """
    Get the alignment represented by the CIGAR string in "alignment matrix" format.
    Example:
      ref_seq    = "ACTCAT"
      read_seq   = "TCAGG"
      cigar      = "2D3M2I1D"

               =>
       
      ref_align  = "ACTCA--T"
      read_align = "--TCAGG-"

    Parameters
    ----------
    ref_seq  : the reference nucleotide sequence
    ref_pos  : position where the alignment starts on the reference sequence (1-based).
    read_seq : the read nucleotide sequence
    cigar    : the CIGAR string (from Bowtie2) representing the alignment

    Returns
    -------
    A tuple (ref_align, read_align) :
      ref_align  : the ref alignment string
      read_align : the read alignment string
    These are strings of the same length which contain letters "A", "C", "G", "T", "-".
  """
  cigar_parsed = cigar_utils.parse_cigar(cigar)

  ref_align = ''
  read_align = '-' * (ref_pos - 1)
  read_pos = 1

  for variation_info in cigar_parsed:
    var_type = variation_info['type']
    var_count = variation_info['count']
    if var_type == 'M':
      for _ in range(var_count):
        ref_align += ref_seq[ref_pos - 1]
        read_align += read_seq[read_pos - 1]
        ref_pos += 1
        read_pos += 1
    elif var_type == 'I':
      for _ in range(var_count):
        ref_align += '-'
        read_align += read_seq[read_pos - 1]
        read_pos += 1
    elif var_type == 'D':
      for _ in range(var_count):
        ref_align += ref_seq[ref_pos - 1]
        read_align += '-'
        ref_pos += 1
  return ref_align, read_align

def count_variations(ref_align, read_align):
  """
    Return the number of insertions, deletions, and substitions in the alignment.

    Parameters
    ----------
    ref_align  : the alignment string for the reference sequence.
    read_align : the alignment string for the read sequence.

    Returns
    -------
    A tuple (num_ins, num_del, num_subst) :
      num_ins   : the number of insertions
      num_del   : the number of deletions
      num_subst : the number of substitutions
  """
  assert len(ref_align) == len(read_align), "Alignment strings must be the same length"
  num_ins = 0
  num_del = 0
  num_subst = 0
  for i in range(len(ref_align)):
    if ref_align[i] != read_align[i]:
      if ref_align[i] == '-':
        num_ins += 1
      elif read_align[i] == '-':
        num_del += 1
      else:
        # both ref and read are nucleotides but not equal
        num_subst += 1
  return num_ins, num_del, num_subst
      

def get_variation_pos(ref_align, read_align):
  assert len(ref_align) == len(read_align), "Alignment strings must be the same length"

  ref_pos = 1
  ins_pos = []
  del_pos = []
  subst_pos = []
  for i in range(len(ref_align)):
    if ref_align[i] == '-':
      ins_pos.append(ref_pos - 1) # insertions are mapped to the previous reference positition
    elif read_align[i] == '-':
      del_pos.append(ref_pos)
      ref_pos += 1
    elif ref_align[i] != read_align[i]:
      subst_pos.append(ref_pos)
      ref_pos += 1
    else: # match
      ref_pos += 1
  return ins_pos, del_pos, subst_pos

def get_orig_seq(align_str):
  """
    Get the original sequence of the alignment string.
    Equivalently, remove the "-" characters from align_str.
  """
  return ''.join(x for x in align_str if x != '-')

def get_cigar(ref_align, read_align):
  """
    Construct a CIGAR string from an alignment.

    Parameters
    ----------
    ref_align  : the alignment string for the reference sequence.
    read_align : the alignment string for the read sequence.

    Returns
    -------
    cigar : the CIGAR string of the alignment in the same format output by Bowtie2.

    Example:
      ref_align  : ATCGA--TGCT
      read_align : --CGATTTGAT
               =>
      cigar : 2D3M2I4M

    Note: substitutions and matches are both converted to "M".
  """
  assert len(ref_align) > 0, 'Empty alignment string'
  assert len(ref_align) == len(read_align), 'Unequal lengths alignment strings'

  cigar = ''
  last_type = None
  curr_count = None
  for i in range(len(ref_align)):
    curr_type = None
    if ref_align[i] == '-': # insertion
      curr_type = 'I'
    elif read_align[i] == '-': # deletion
      curr_type = 'D'
    else: # substitution or match
      curr_type = 'M'
    if curr_type != last_type:
      if last_type is not None:
        cigar += str(curr_count) + last_type
      last_type = curr_type
      curr_count = 1
    else:
      curr_count += 1
  cigar += str(curr_count) + last_type
  return cigar