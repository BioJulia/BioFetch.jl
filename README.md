# BioFetch.jl
 Easily fetch biological sequences from online sources

BioFetch provides a higher-level interface to retrieve data from sequence databases and provide them in
a readily manipulatable form via FASTX.jl and GenomeAnnotations.jl.

Currently supports Entrez (NCBI) Nucleotide and Protein databases, as well as UniProt and Ensembl.

Examples:
```julia
fetchseq("AH002844")                             # retrive one NCBI nucleotide record as FASTA
fetchseq("CAA41295.1", "NP_000176", format = gb) # retrieve two NCBI protein records as GenBank Flat File
fetchseq("Q00987")                               # retrieve one UniProt protein record as FASTA
fetchseq("ENSG00000141510")                      # retrieve one Ensembl gene record's proteins as FASTA
fetchseq("NC_036893.1", range = 81_775_230 .+ (1:1_000_000))         # retrive a 1 Mb segment of a FASTA NCBI genomic record
fetchseq("NC_036893.1", range = 81000000:81999999, revstrand = true) # retrive a 1 Mb segment of a FASTA NCBI genomic record on the reverse strand
```
