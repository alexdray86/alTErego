# alTErego

[Read alTErego documentation](https://alexdray86.github.io/alTErego/build/html/index.html)

Compute transposable element (TE)-derived regulons from gene expression at the single-cell level. We use a similar approach to SCENIC method, where we replace transcription factor (TF)-specific motifs by transposable element sequences. The goal is to predict triplets of (transcription factor, transposable element family, putative targets) by using clues from co-expression in scRNA-seq integrated with TF-TE known interaction derived from a ChIP-seq database (by default, we use ENCODE + Imbeault et al. 2018 datasets). 
