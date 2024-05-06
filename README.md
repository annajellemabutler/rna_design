# RNA Design

## Project Outline
This script identifies primary RNA sequences which meet very specific primary and secondary structure constraints. It was used to design synthetic RNA molecules on which a novel high-throughput RNA structure determination method was tested and validated. The project formed part of a five-week laboratory rotation I completed in Prof. Sid Dey's lab at the University of California, Santa Barbara. 

## Script Overview 
### Constraints
`rna_design.py` generates a `.csv` file containing RNA sequences which match both:
* Primary structure constraints
  * 68 bp length
  * ~50% overall GC content
  * 3' tail of 15 adenosines
  * Reasonable distribution of cytosines, including
    * ≥ 1 cytosine in the loop region
    * 2-3 cytosines in the 5' tail, including at least a few with preceding adenosine
    * 5-6 cytosines in the stem, but 0 in the first 2 bp
* Energy & folding constraints
  * Predicted minimum free energy (MFE) structure is a single stem-loop (with no bulges or interior loops)
  * Predicted suboptimal structures have absolute Gibbs free energy values ≤ 2 kcal/mole smaller than the MFE structure

### Pipeline
This is completed in the following steps:
1. An RNA sequence with the specified length and GC content is randomly generated. These sequences contain higher GC content in the stem region to introduce bias towards forming stable stem loops. 
2. The minimum free energy secondary structure is predicted using `RNAfold` and compared to the desired secondary structure (by default, a perfect stem-loop). 
3. If the structure is a match, the free energy of the optimal structure is compared with the free energy of the nextmost suboptimal secondary structure (generated using `RNAsubopt`). Only if the absolute difference in dG is greater than the specified value is the sequence kept.
4. The sequence set is filtered for those which satisfy the primary structure constraints (nucleotide counts, etc.).


