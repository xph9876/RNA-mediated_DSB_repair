import common

def get_alignment_window(
  ref_align,
  read_align,
  dsb_pos,
  window_size,
  anchor_size,
  anchor_mismatch_limit,
):
  ref_i = 1 # index on the original reference sequence (without "-")
  read_i = 1 # index on the original read sequence (without "-")

  ref_align_window = ''
  read_align_window = ''

  window_start, window_end = common.get_window_range(dsb_pos, window_size)
  left_anchor_start = window_start - anchor_size
  left_anchor_end = window_start - 1
  right_anchor_start = window_end + 1
  right_anchor_end = window_end + anchor_size

  left_anchor_mismatches = 0
  right_anchor_mismatches = 0
  for i in range(min(len(ref_align))):
    if ref_i in range(left_anchor_start, left_anchor_end + 1):
      # Check the mismatches/in/dels on the left anchor
      if ref_align[i] != read_align[i]:
        left_anchor_mismatches += 1
      elif (ref_align[i] == '-') or (read_align[i] == '-'):
        left_anchor_mismatches = np.inf
    elif ref_i in range(window_start, window_end + 1):
      # extract the window around the DSB
      ref_align_window += ref_align[i]
      read_align_window += read_align[i]
    elif ref_i in range(right_anchor_start, right_anchor_end + 1):
      # Check the anchor/in/dels on the right anchor
      if ref_align[i] != read_align[i]:
        right_anchor_mismatches += 1
      elif (ref_align[i] == '-') or (read_align[i] == '-'):
        right_anchor_mismatches = np.inf

    # increment counters
    if ref_align[i] != '-':
      ref_i += 1
    if read_align[i] != '-':
      read_i += 1

    if ref_i > right_anchor_end:
      break

  if ref_i <= right_anchor_end:
    print(
      'Warning: read did not align across window and anchor:\n' +
      ref_align + '\n' +
      read_align + '\n'
    )
    return None

  if (
    (left_anchor_mismatches > anchor_mismatch_limit) or
    (right_anchor_mismatches > anchor_mismatch_limit)
  ):
    return None

  return {
    'ref_align': ref_align_window,
    'read_align': read_align_window,
  }
