.. pyTEnrich documentation master file, created by
   sphinx-quickstart on Wed Dec 30 11:29:02 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. image:: images/alTErego_logo_new.png

.. toctree::
   :maxdepth: 2
   :caption: Contents: 
   
   usage/installation.rst
   usage/execution.rst
   source/modules.rst

Overview of the methods
=======================

This program aims at defining regulons of putative targets for Transcription Factors (TFs) targetting Transposable Elements (TEs) sequences. The whole method is built on the hypothesis that TEs found in close proximity of gene promoters are potentially behaving as enhancers. 

By taking advantage of GRNBoost2 from `arboreto <https://arboreto.readthedocs.io/en/latest/>`_  to compute modules of co-expressed TFs and putative targets. We then used TE abundance around TSS with a similar approach that `SCENIC method <https://www.nature.com/articles/nmeth.4463>`_ uses with motifs. Using enrichment of TF->TE interactions from a large panel of ChIP-seq data, we then prune the results to keep only regulons from known TF->TE interactions. 

.. image:: images/alTErego_design.png

