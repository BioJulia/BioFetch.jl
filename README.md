# BioFetch.jl
 Easily fetch biological sequences from online sources

BioFetch provides a higher-level interface to retrieve data from sequence databases and provide them in
a readily manipulatable form via FASTX.jl and GenomeAnnotations.jl.

Currently only supports Entrez (NCBI) Nucleotide ("nuccore") and Protein ("protein") databases, but
the intent is to generalize it.

Examples:
```julia
fetchseq(nuccore, "AH002844")                             # retrive one FASTA nucleotide record
fetchseq(protein, "CAA41295.1", "NP_000176"; format = gb) # retrieve two GenBank flatfile protein records
```
