# pDm8 and yDm8 in the optic-lobe multiome

**Truong Nguyen, CUNY School of Medicine**

Primary data: GSE305940, a 10x single-nucleus RNA + ATAC multiome of the developing *Drosophila melanogaster* optic lobe. Expression was checked in an independent optic-lobe scRNA-seq atlas, GSE156455. No other accession was analyzed.

## Question

Dm8 neurons are one published class. beat-IIIc (CG15138) and DIP-gamma (CG14521) were used to ask whether pale and yellow cells can be separated inside that class, and whether the two groups differ in other genes and in chromatin accessibility.

Labels use raw RNA counts:

- **pDm8:** beat-IIIc > 0 and DIP-gamma = 0
- **yDm8:** DIP-gamma > 0 and beat-IIIc = 0
- Cells with both markers, or neither, were held out of the differential tests

The two markers are largely mutually exclusive among marker-positive nuclei: 86.6% at P24 and 86.8% at P48. Adult and P0 were not tested. Adult beat-IIIc is too rare in this nuclear RNA assay, and P0 has no genotype labels to pair on.

## Design

For each stage, counts were summed within each DGRP line and subtype. A line was kept only if it had at least 5 nuclei of each subtype: 10 lines at P24 (198 pDm8 vs 353 yDm8) and 11 lines at P48 (159 vs 224). The model is `~ line + subtype` on TMM log-counts per million, fit with limma-trend. limma-voom, edgeR quasi-likelihood, a female-line fit, and leave-one-line-out were used as checks. Positive log2 fold-change means higher in pDm8.

A gene was tested if it was detected in at least 10% of nuclei of one subtype. A peak was tested at a 5% detection threshold. The reporting cut is FDR < 0.05 and absolute log2 fold-change ≥ 0.5.

## Expression

beat-IIIc and DIP-gamma are the labeling genes, so their differential expression is a control.

Three other genes stay higher in yDm8 at both P24 and P48 in limma-trend, voom, edgeR, the female-line fit, and every leave-one-line-out:

| Gene | P24 log2FC (p/y) | P48 log2FC (p/y) |
|---|---:|---:|
| nAChRalpha6 | −2.61 | −1.32 |
| CG44153 | −3.04 | −2.01 |
| Tl | −3.63 | −2.18 |

nAChRalpha6 is detected in most cells of both subtypes, so it is a shift in level. Sulf1, Drl-2, and kek6 agree in direction but do not hold in every check.

In GSE156455, the same raw-count rule gives 122 pDm8 vs 167 yDm8 cells at 24 h and 207 vs 220 at 48 h. Grouped by genotype, with at least 3 nuclei of each subtype, nAChRalpha6, CG44153, and Tl are higher in yDm8 in every genotype at both times.

## Chromatin

Six peaks inside Vsx2, on chromosome X, are more open in pDm8 at both stages (line-level log2 fold-change about +2 to +4). They remain significant in voom, edgeR, and the female-line fit. Vsx2 RNA itself was too rare to test, so this is a chromatin result, not an mRNA result.

At DIP-gamma, accessibility does not point one way: some peaks are higher in yDm8 and some in pDm8. No tested peak at beat-IIIc separates the subtypes.

Coverage tracks at these three loci count each fragment once and are in `figures/09_coverage_tracks.png`. The y-axis is mean coverage per base, per million ATAC peak counts.

## Peak–gene links

A link requires a significant peak, a differentially expressed gene at the same stage, the same direction of effect, and either a transcription start site within 50 kb or a published peak–gene link from GSE305940 at p < 0.05. That definition returns 21 associations, all higher in yDm8. Fifteen involve DIP-gamma. These are candidates. They are not a demonstration that the peaks are enhancers.

## What this does not show

- Enhancer function was not tested.
- Each DGRP line is one animal at a stage. Ten or eleven lines is a real replicate design, and a modest one.
- Sex is constant within each line, so it cannot be added on top of line. A female-line fit is the sex check. Nuclei labeled male express roX1 and roX2; nuclei labeled female rarely do.
