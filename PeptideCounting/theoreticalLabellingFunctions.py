def aminoacidDetect(proteinsFASTA, aminoacids):

    """
    Detects specified amino acids in all sequences from a fasta file.

    params:
        proteinsFASTA: fasta file of protein sequences
        aminoacids: Amino acids to detect as a list of strings
    output:
        aaCounts: dictionary of the counts of specified amino acids the sequences.
    """

    from Bio import SeqIO

    #start dictionary where amino acids will be counted
    aaCounts = {}

    #For each amino acid specified, start a dictionary key that is currently set to 0.
    for aa in aminoacids:
        aaCounts[aa] = 0

    #Iterate through sequences in the fasta file 
    for fasta in SeqIO.parse(proteinsFASTA, "fasta"):
        sequence = fasta.seq

        #Iterate through amino acid in the protein
        for letter in sequence:
            if letter in aminoacids:
                #If the letter is in aaCounts dictionary, add 1 to that value
                aaCounts[letter] += 1
        
    return aaCounts

def cutSiteSearchRegions(proteinSequence, cutsites = ['R', 'K']):
    """
    Determines search regions for peptides - proteinSequences cleaved to start at each possible cutsite.

    params:
        proteinSequence: the protein sequence of interest (str)
        cutSites: list of amino acids at which the enzyme cuts (list of str)

    returns:
        searchRegions: proteinSequence cleaved to begin at each cutsite in the sequence (list of str)
        
     """
    #This will be a list of all cut sites in the sequence, with start 0. start sites of peptides will iterate through this list.
    cutSitePositions = []

    #Iterate through the sequence and collect all cut site positions
    for index, aa in enumerate(proteinSequence):
        if aa in cutsites:
            cutSitePositions.append(index)
        
    searchRegions = [proteinSequence] #starting region will be the whole protein (starts at 0)
    for cutSite in cutSitePositions:
            searchRegions.append(proteinSequence[cutSite+1:]) #note that this is +1 as the cut peptide itself should be left of the cut site
            
    return searchRegions

def potentialPeptides(proteinSequence, cutSites = ['R', 'K'], missedcleavagemax = 2):
    """
    Calculates all potential peptides for a protein sequence given the amino acids at which the enzyme cuts.

    params:
        proteinSequence: string of the protein sequenc of interest.
        cutSites: list of strings of amino acids at which the enzyme cuts. By default, this is set for trypsin.
    
    output:
        allPeptidesList: all potential peptides for the given inputs.
    """
    #Start list for the output - all potential peptides for this protein
    allPeptidesList = []

    #Counter for missed cleavages (allows for two)

    searchRegions = cutSiteSearchRegions(proteinSequence, cutsites = cutSites)
    
    #Iterate through the string and add potential peptides according to allowed missed cleavages.
    for searchRegion in searchRegions:
        peptideCount = 0
        for index, aa in enumerate(searchRegion):
        #If cut site is encountered, save sequence
            if aa in cutSites:
                stopaa = index 
                allPeptidesList.append(searchRegion[:stopaa+1]) #again to keep the cut protein on the left side
                peptideCount +=1
                #Once we've accounted for two missed cleavages, stop iterating and move to the next start site
                if peptideCount == missedcleavagemax+1: #missed cleavage allowance +1 = expected number of peptides output
                    break
            elif index == len(searchRegion)-1:
                allPeptidesList.append(searchRegion[:index+1])
                break                

    #Remove peptides which are out of range from the full list
    #allPeptidesList = inRangePeptideCalc(allPeptidesList, minSize = 5, maxSize = 10) hashed out for now as the min/max is nonsensical

    allPeptidesLength = len(allPeptidesList)
    
    return allPeptidesLength

def inRangePeptideCalc(peptideList, minSize = 5, maxSize = 10):
    """
    Determine the proportion of peptides that are undetectable given a size range.count

    params:
        peptideList: peptides to be considered
        minSize: minimum detectable peptide size.
        maxSize: maximum detectable peptide size.

    output:
        detectableRatio: ratio of peptides that are in range/total
        inRangePeptides: truncated peptideList for only proteins within the range
    """
    #note that if the min/max needs to be in kDa, need a kDa calculator function here

    inRangePeptides = []

    for peptide in peptideList:
        if len(peptide) > minSize and len(peptide) < maxSize:
            inRangePeptides.append(peptide)

    detectableRatio = len(inRangePeptides)/len(peptideList)

    return detectableRatio, inRangePeptides

def returnAllPeptideCounts(proteinsfasta, maxmissedcleave = 2):
    """
    for a set of proteins, the total peptide count is given for when the protein is unlabelled and labelled at all lysines. 

    params:
        proteinsfasta: a fasta file of all protein sequences of interest
        maxmissedcleave: number of missed cleavages allowed (how many cut sites the enzyme might skip at most)
    
    output:
        totalUnlabelled: total count of peptides if the sequences are unlabelled.
        totalLabelled: total count of peptides if the sequences are labelled at all lysines.
    """
    from Bio import SeqIO

    totalUnlabelled = 0
    totalLabelled = 0

    #Iterate through sequences in the fasta file 
    for fasta in SeqIO.parse(proteinsfasta, "fasta"):
        sequence = fasta.seq

        totalUnlabelled += potentialPeptides(sequence, cutSites = ['R', 'K'], missedcleavagemax= maxmissedcleave)
        totalLabelled += potentialPeptides(sequence, cutSites = ['R'], missedcleavagemax= maxmissedcleave)
        
    return totalUnlabelled, totalLabelled
