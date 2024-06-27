## RNA Design

Script generates semi-random primary RNA sequences under specific primary and secondary structure constraints. Originally used in the design of synthetic RNA molecules on which a novel high-throughput RNA structure determination method was tested and validated. 

### Project-specific constraints
* Primary structure
  * 68 bp length
  * 50% overall GC content
  * Long (≥ 10 nucleotides) 3' poly-A tail
  * Reasonable distribution of cytosines, including
    * ≥ 1 cytosine in the loop region
    * 2-3 cytosines in the 5' tail, including at least a few with preceding adenosine
    * 5-6 cytosines in the stem, but 0 in the first 2 bp
* Energy & folding constraints
  * Predicted minimum free energy (MFE) structure is a single stem-loop with no bulges or interior loops
  * Predicted suboptimal structures have absolute Gibbs free energy values ≤ 2 kcal/mole smaller than the MFE structure

### Overview
1. An RNA sequence with the specified length and GC content is randomly generated. GC content is increased in the stem region to bias the region towards forming stable stem loops. 
2. The minimum free energy secondary structure is predicted using `RNAfold` and compared to the desired secondary structure (a perfect stem-loop, by default). 
3. If the structure is a match, the free energy (dG) of the optimal structure is compared with the dG of the nextmost suboptimal secondary structure (generated using `RNAsubopt`), and the sequence is retained only if the absolute difference in dG is greater than a specified value. 
4. The sequence set is filtered for those which satisfy the remaining primary structure constraints.
5. A `.csv` file containing the matching sequences is saved. 

