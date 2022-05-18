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
        # both are nucleotides but not equal
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


def get_indel_pos(ref_align, read_align):
  """
    Get the insertion and deletion positions of the alignment.
    The returned list of positions are 0-based indices that index ref_align and seq_align.

    Parameters
    ----------
    ref_align  : the alignment string for the reference sequence.
    read_align : the alignment string for the read sequence.

    Returns
    -------
    positions : a list of 0-based indices, with respect to the input alignment strings,
      which give the positions of the insertions, deletions, and substitions.
  """
  assert len(ref_align) == len(read_align), "Alignment strings must be the same length"
  positions = []
  for i in range(len(ref_align)):
    if (ref_align[i] == '-') or (read_align[i] == '-'):
      positions.append(i)
  return positions

def get_ref_positions(ref_align, read_align, ref_pos):
  """
    Get the map from alignment indices to reference indices.

    Example:
      input:
        ref_pos    = 1
        ref_align  = A C G - - T C A
        read_align = - - G G T - - C
      output:
        pos_map    = 1 2 3 3 3 4 5 6

    Parameters
    ----------
    ref_align  : the alignment string for the reference sequence.
    read_align : the alignment string for the read sequence.
    ref_pos    : left-most position on the reference sequence that the read aligns with.

    Returns
    -------
    pos_map : a list with same length as ref_align.
      pos_map[i] gives the 1-based position on reference sequence corresponding
      with ref_align[i] and seq_align[i].
  """
  assert len(ref_align) == len(read_align), "Alignment strings must be the same length"
  pos_map = []
  for i in range(len(ref_align)):
    if ref_align[i] == '-':
      pos_map.append(ref_pos - 1) # insertions are mapped to the previous reference positition
    else:
      pos_map.append(ref_pos)
      ref_pos += 1
  return pos_map

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
  cigar_long = ''
  for i in range(len(ref_align)):
    if ref_align[i] == '-': # insertion
      cigar_long += 'I'
    elif read_align[i] == '-': # deletion
      cigar_long += 'D'
    else: # substitution or match
      cigar_long += 'M'
  cigar = ''
  for match in re.findall(r'I+|D+|M+', cigar_long):
    cigar += str(int(len(match))) + match[0]
  return cigar