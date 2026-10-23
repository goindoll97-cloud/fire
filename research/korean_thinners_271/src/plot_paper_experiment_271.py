# -*- coding: utf-8 -*-
"""Optional manuscript figures from validated experiment summary; needs matplotlib."""
from __future__ import annotations
import argparse
from pathlib import Path
from paper_experiment_271 import read_csv


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--reports',type=Path,required=True)
    args=ap.parse_args()
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows=read_csv(args.reports/'paper_제품별_안정성.csv')
    by={}
    for r in rows:
        by.setdefault(r['제품ID'],{})[r['표본추출']]=r
    pids=sorted(by)
    labels=[];legacy=[];gibbs=[]
    for pid in pids:
        if 'legacy' not in by[pid] or 'gibbs' not in by[pid]:
            continue
        labels.append(pid)
        legacy.append((float(by[pid]['legacy']['시드별_압력비중앙값_중앙값'])-1)*100)
        gibbs.append((float(by[pid]['gibbs']['시드별_압력비중앙값_중앙값'])-1)*100)
    if not labels:
        raise ValueError('Both legacy and gibbs sampling required for paired figure')
    xs=list(range(len(labels)))
    width=.38
    fig,ax=plt.subplots(figsize=(9,4.6))
    ax.bar([x-width/2 for x in xs],legacy,width,label='Legacy sampler')
    ax.bar([x+width/2 for x in xs],gibbs,width,label='Bounded-simplex Gibbs sampler')
    ax.set_xticks(xs,labels)
    ax.set_ylabel('Median HANNA/Raoult pressure ratio - 1 (%)')
    ax.set_title('Sampling-method sensitivity of equilibrium-model discrepancy')
    ax.legend()
    fig.tight_layout()
    fig.savefig(args.reports/'Figure_1_sampling_sensitivity.png',dpi=220)
    plt.close(fig)
    print('Created Figure_1_sampling_sensitivity.png; figures are not model validation')

if __name__=='__main__':
    main()