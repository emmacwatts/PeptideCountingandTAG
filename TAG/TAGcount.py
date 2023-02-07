def maptoGenomeScaffold(transcript):
    """
    Maps an NbD transcript sequence to the longest scaffold in the NbD genome.

    """

    from Bio import SeqIO
    import pandas as pd

    scaffoldSequence = ""

    for fasta in SeqIO.parse("TAG/NbDE_Genome.fasta", "fasta"):
        sequence = fasta.seq

        if transcript in str(sequence):
            if len(sequence) > len(scaffoldSequence):
                scaffoldSequence = sequence
    
    return scaffoldSequence

def nextStop(transcript, scaffoldSequence, stopCodons = ["TAA", "TGA"]):

    transcriptStartInd = scaffoldSequence.find(transcript)

    for ind, base in enumerate(scaffoldSequence[transcriptStartInd:]):
        if ind % 3 == 0 and scaffoldSequence[ind: ind+3] in stopCodons:
            extendedTranscript = scaffoldSequence[transcriptStartInd: ind]
            break

    return extendedTranscript

def TAGextensions(transcriptomeFile = "TAG/NbDE_Transcriptome.fasta"):
    """
    Evaluates TAGs in the transcriptome - total number and sequences if extended.

    """
    from Bio import SeqIO
    import pandas as pd

    TAGcount = 0

    extendedTAGsdf = pd.DataFrame(columns=["NbdID", "OriginalSeq"])

    for fasta in SeqIO.parse(transcriptomeFile, "fasta"):
        sequence = fasta.seq

        for codonInd in range(0, len(sequence), 3):
            codon = sequence[codonInd:codonInd+3]
            if codon == 'TAG':
                TAGcount+= 1
                extendedTAGsdf.loc[len(extendedTAGsdf)] = [fasta.id, str(fasta.seq)]
            
    extendedTAGSdf["ScaffoldSeq"] = maptoGenomeScaffold(extendedTAGsdf["OriginalSeq"])

    return TAGcount, extendedTAGsdf 

#Test/runner
def testmaptoGenomeScaffold():
    import time
    
    start = time.perf_counter()
    scaffoldSequence = maptoGenomeScaffold("ATGGATCTTGATCTTACACCAAAGTTGGCAAAGCAAGTGTACGGAGGAGATGGTGGTTCTTATCATGCATGGTGTCCTAATGATCTGCCTATGTTAAAAGAAGGAAACATTGGTGGTGCTAAACTTGCTCTTTCTAAGAATGGTTTTGCTTTGCCTCGTTACTCTGATTCTGCTAAAGTTGCCTATGTTCTTCAAGGTTGTGGAGTTGCTGGAATTGTTCTTCCAGAGAAAGAAGAGAAGGTGCTAGCGATTAAGACAGGAGATGCTATAGCCCTTCCTTTTGGTGTTGTAACATGGTGGTACAACAAAGAGGACACTGAACTTGTGATTCTGTTCCTTGGTGATACCAAAACTGCAAACAAAGCTGGTTCATTCACAGACATGTACTTAACTGGCTCGAATGGCATTTTCACTGGTTTTTCAACCGAGTTTGTTAGCAGGGCATGGAATGTAGAAGAGAGTGTTGCTAAAACTCTTGTTAGCTCCCAAACTGCTCAGGGTATCGTGAAGCTTGACGCAGGCTTCCAAATGCCCGAGCCTAAGCAAGGCCACAGGGATGGTATGGTTCTCAACTGTTTGGAAGCCCCATTGGATGTTGACATTAAGGGTGGTGGAAAGGTTGTTGTTTTGAATACCAAGAACTTGCCTTTGGTCGGTGAAGTTGGACTTGGTGCTGATCTTGTGAGGTTGAATGGAAGTGCTATGTGCTCTCCTGGTTTTTCATGTGACTCAGCTCTTCAGGTCACTTACATTGTTAGAGGCAGTGGCCGAGTTCAAGTCGTTGGCCCTGATGGTAAGCGCGTTCTTGAAACACACATCAAGGCGGGCAATCTCTTCATCGTTCCAAGG")    
    stop = time.perf_counter()

    print(f"Ran functions in {stop - start:0.4f} seconds.")

    return scaffoldSequence

def testTAGextensions():
    import time
    
    start = time.perf_counter()
    TAGcount, extendedTAGSdf = TAGextensions("TAG/Mock_NbDE_Transcriptome.fasta")
    stop = time.perf_counter()
    print(f"Ran functions in {stop - start:0.4f} seconds. Calculated total TAG count as {TAGcount}")

    return extendedTAGSdf

scaffoldSequence = testmaptoGenomeScaffold()

extendedTAGSdf = testTAGextensions()
#Ran functions in 23.5676 seconds. Calculated total TAG count as 17001

