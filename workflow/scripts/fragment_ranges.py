"""ENCODE's fragment length ranges in bp: NFR below 150, mono-nucleosome above 150 up to 300.

fraglen_qc.py scores with them, and qc_figures.py shades them. This module
imports nothing, so the report env can load it.
"""

NFR_UPPER = 150
MONO_LOWER = 150
MONO_UPPER = 300
