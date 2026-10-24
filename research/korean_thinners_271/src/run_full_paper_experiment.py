# -*- coding: utf-8 -*-
"""Full pipeline: 00–09 genuine raw/evidence workflow + three manuscript RQs.

Examples:
 python research/korean_thinners_271/src/run_full_paper_experiment.py --offline
 python research/korean_thinners_271/src/run_full_paper_experiment.py --offline \
   --noru-zip PATH --samhwa-zip PATH --samples 200
"""
import argparse
from pipeline_270_core import cli as run_00_to_09
from paper_rq_270 import main as run_rqs

if __name__ == '__main__':
    # The existing pipeline CLI already owns all flags; forward unchanged.
    run_00_to_09('all')
    run_rqs([])