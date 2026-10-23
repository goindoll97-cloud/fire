# -*- coding: utf-8 -*-
"""Execute evidence-only 02a, 02b, optional 02c, and 02d in one run.

Never creates regulatory, molecular-property or product-compatibility approvals.
"""
import argparse
from pathlib import Path
from auto_enrich_270 import property_enrich, official_twa_evidence, product_evidence, report_queue


def main(argv=None):
    p=argparse.ArgumentParser(description='Automatic 270 MSDS research evidence enrichment; approvals remain manual')
    p.add_argument('--noru-zip',type=Path)
    p.add_argument('--samhwa-zip',type=Path)
    p.add_argument('--limit',type=int,default=54)
    p.add_argument('--offline',action='store_true')
    p.add_argument('--no-apply',action='store_true')
    a=p.parse_args(argv)
    if not 1<=a.limit<=54:p.error('--limit must be 1..54')
    if bool(a.noru_zip)!=bool(a.samhwa_zip):p.error('Both ZIP arguments are required together')
    print('[02a] PubChem evidence' + (' (offline: omitted)' if a.offline else ''))
    if not a.offline:
        property_enrich(limit=a.limit,apply=not a.no_apply)
    print('[02b] KOSHA and Korean Law official evidence')
    official_twa_evidence(limit=a.limit,offline=a.offline)
    if a.noru_zip and a.samhwa_zip:
        print('[02c] MSDS original section-1 use evidence')
        product_evidence(a.noru_zip,a.samhwa_zip)
    else:
        print('[02c] SKIPPED: two original PDF archives not provided')
    print('[02d] Remaining review queue')
    report_queue()
    print('AUTO EVIDENCE FINISHED: no HANNA/TWA/same-use approval has been granted')
if __name__=='__main__':main()