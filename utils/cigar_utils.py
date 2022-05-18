import re

def parse_cigar(cigar):
  """
    Parse a CIGAR string into runs of variations.
    The CIGAR is assumed to be from the output of Bowtie2 and only contain
    only integers and the alphabets "I" (insertion), "D" (deletion), "M" (mismatch).
    Please refer to the SAM format spec or Bowtie2 documentation for more info.

    Parameters
    ----------
    cigar : the CIGAR string

    Returns
    -------
    variations : a list of dictionaries with items:
      "type" -> one of "I", "D", "M"
      "count" -> a positive integer

    Example
    -------
    cigar = "4I40M3D"
       ->
    variations = [
      {"type": "I", "count": 4},
      {"type": "M", "count": 41},
      {"type": "D", "count": 3},
    ]
  """
  variations = re.findall(r'(\d+)([IMD])', cigar)
  variations = [
    {'type': x[1], 'count': int(x[0])}
    for x in variations
  ]
  return variations