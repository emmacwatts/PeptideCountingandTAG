# Tasks: 

* within the secretome (subset of proteome) calculate all potential peptides √
* Per protein
    * Identify all K - these are potential label sites √

* Identify all R and K as potential cut sites √

* Cut for peptides - allowing for 2 missed cleavages √

* Identify those that are outside of detectable range (X-Y) and the ratio of indetectables across the proteome √

* Identify how many of the potential peptides are repeat signals - this depends on the input
* So exclude any that are covering the same peptide based on mapped position

* Labelled peptides will come back as a list per protein

# Questions
* Assumed when labelled, still account for two missed cleavages without counting labelled K
* Is the theoretical maximum still useful if it doesn't consider inactive lysines
* Secretome NbD IDs list?
* Rough cut-off for size max and min
* The other side to consider was how many peptides that you get back determine a unique signal - although the way this code is written determines all peptides that should come back so this might not be useful.

# Outcomes
* For the full proteome, the number of peptides if fully labelled lysines (i.e. theoretical maximum) is 5,185,205 allowing 2 missed cleavages (not yet accounting for a detectable range). Total peptides if unlabelled is 11,549,494.