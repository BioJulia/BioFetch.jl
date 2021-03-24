module BioFetch

using FASTX
using GenomicAnnotations
using BioServices.EUtils

@enum Format fasta gb

@enum Database nuccore protein

"""
    fetchseq(db::Database, ids::AbstractString...; format::Format = fasta)

Fetches sequence data from a database by accession number in either FASTA format or GenBank flatfile format.

Currently only supports `nuccore` and `protein` Entrez databases.

```julia
fetchseq(nuccore, "AH002844")                             # retrive one FASTA nucleotide record
fetchseq(protein, "CAA41295.1", "NP_000176"; format = gb) # retrieve two GenBank flatfile protein records
```
"""
function fetchseq(db::Database, ids::AbstractString...; format::Format = fasta)
    response = efetch(db = String(Symbol(db)), id = ids, rettype = String(Symbol(format)), retmode="text")
    body = IOBuffer(response.body)
    if format == fasta
        reader = FASTA.Reader(body)
        records = [record for record ∈ reader]
    else
        records = readgbk(body)
    end
    return records
end

end