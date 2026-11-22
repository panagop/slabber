"""Source-CSV shape, grid, coverage and benchmark integrity checks.

This is NOT a transcription audit against printed sources. Only built-in
Python modules are used to permit quick testing before numerical solvers.
"""
import csv
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'czerny_dataset'


def read(name):
    with (DATA/name).open(newline='', encoding='utf-8-sig') as f:
        return list(csv.DictReader(f))


def verify():
    four = read('czerny_coefficients.csv')
    three = read('czerny_three_edge_coefficients.csv')
    coverage = read('coverage.csv')
    assert len(four) == 1701, f'Four-edge entries changed: {len(four)}'
    assert len(three) == 456, f'Three-edge entries changed: {len(three)}'
    assert len(coverage) == 15, f'Cases changed: {len(coverage)}'
    assert len({r['case_id'] for r in four}) == 9
    assert len({r['case_id'] for r in three}) == 6
    assert {r['case_id'] for r in four}.isdisjoint({r['case_id'] for r in three})
    keys=set()
    for label, rows in [('four',four),('three',three)]:
        for row in rows:
            key=row['case_id'],row['quantity'],float(row['ratio_ly_lx'])
            assert key not in keys, f'Duplicate key ({label}): {key}'
            keys.add(key)
            val=float(row['table_value'])
            assert val > 0, (label,row)
            assert row['source_url'].startswith('https://'), (label,row)
    ratios = defaultdict(set)
    counts = Counter()
    for row in four:
        ratios[row['case_id']].add(round(float(row['ratio_ly_lx']),2))
        counts[row['case_id']]+=1
    expected={round(1+i*.05,2) for i in range(21)}
    for case,nodes in ratios.items():
        assert nodes == expected, f'Incomplete 0.05 grid: {case}'
    for row in coverage:
        assert row['edges_x0_y0_xa_yb'].count(',')==3
    print('PASS: 9 four-edge cases, 1,701 coefficients and 21-ratio grids')
    print('PASS: 6 three-edge cases and 456 Hahn coefficients')
    print('PASS: 15 coverage entries, positive values, unique keys, source URLs')
    print('NOTE: source-to-scan transcription not independently audited')


if __name__ == '__main__':
    verify()
