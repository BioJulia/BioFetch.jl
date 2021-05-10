# BioFetch.jl
 Easily fetch biological sequences from online sources

BioFetch provides a higher-level interface to retrieve data from sequence databases and provide them in
a readily manipulatable form via FASTX.jl and GenomeAnnotations.jl.

Currently supports Entrez (NCBI) Nucleotide and Protein databases, as well as UniProt and Ensembl.

Examples:
```julia
fetchseq("AH002844")                # retrive one FASTA NCBI nucleotide record
fetchseq("CAA41295.1", "NP_000176") # retrieve two FASTA NCBI protein records
fetchseq("Q00987")                  # retrieve one UniProt protein record
fetchseq("ENSG00000141510")         # retrieve one Ensembl gene record
```
