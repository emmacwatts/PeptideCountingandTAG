#This file is used in case functions will not run in ipynb due to size.

import theoreticalLabellingFunctions as t
import time

def testerReturnAllPeptideCounts(inputFile):
    start = time.perf_counter()
    totalUnlabelled, totalLabelled = t.returnAllPeptideCounts(inputFile)
    stop = time.perf_counter()
    print(f"Ran functions in {stop - start:0.4f} seconds. Calculated unlabelled peptide total of {totalUnlabelled}, with labelled peptide total {totalLabelled}")

testerReturnAllPeptideCounts("Mock_NbD_Proteome.fasta")
#Ran functions in 0.3215 seconds. Calculated unlabelled peptide total of 657, with labelled peptide total 288

testerReturnAllPeptideCounts("LongMock_NbD_Proteome.fasta")
#Ran functions in 1.3450 seconds. Calculated unlabelled peptide total of 503461, with labelled peptide total 223186

testerReturnAllPeptideCounts("NbD_Proteome.fasta")
#Ran functions in 31.0720 seconds. Calculated unlabelled peptide total of 11549494, with labelled peptide total 5185205

