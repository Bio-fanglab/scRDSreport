# Conservative sample grouping algorithm

## Goal

Given one sample ID per cell, propose a sample-level `group` and `replicate`
without silently treating technical naming artifacts as biological truth.

## Inputs and outputs

- Input: cell metadata, an optional explicit sample column, and an optional user
  sample map.
- Output: one row per sample with `group`, `replicate`, `confidence`,
  `needs_review`, `grouping_rule`, and `n_cells`.

## Decision rule

1. A user map always wins and receives confidence `user`.
2. An explicit terminal `rep1`, `replicate2`, or `r3` token is stripped from
   the group and receives confidence `high`.
3. A donor/subject/patient/mouse token can serve as a replicate block, but the
   proposed remaining condition receives confidence `medium` and needs review.
4. A bare terminal number is accepted only when the same base occurs at least
   twice. It receives confidence `medium` and needs review. Time/dose-like
   prefixes are excluded.
5. Otherwise every sample remains its own group, replicate is missing, and the
   result receives confidence `low`.

Pseudocode:

```
if user_map: return validated user_map
for sample in samples:
    if explicit_rep_suffix(sample): split group and replicate; high confidence
    else if subject_token(sample): remove subject token; medium confidence
    else if repeated_base_plus_number(sample): split; medium confidence
    else: group = sample; replicate = NA; low confidence
```

The algorithm is linear in the total number of characters across sample names,
apart from a small hash-table count of candidate bases. It deliberately does
not infer differential-expression contrasts. Biological replication and
contrast direction must be confirmed by the user.
