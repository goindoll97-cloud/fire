# -*- coding: utf-8 -*-
"""Numbered follow-up after stage 02: pubchem evidence without automatic approval."""
from auto_enrich_270 import main
if __name__ == '__main__':
    import sys
    main(['pubchem', *sys.argv[1:]])