# -*- coding: utf-8 -*-
"""Numbered follow-up after stage 02: products evidence without automatic approval."""
from auto_enrich_270 import main
if __name__ == '__main__':
    import sys
    main(['products', *sys.argv[1:]])