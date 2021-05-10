module BioFetch

using FASTX
using GenomicAnnotations
using BioServices.EUtils

@enum Format fasta gb

"""
    fetchseq(ids::AbstractString...; format::Format = fasta)

Fetches sequence data from a database by accession number in either FASTA format or GenBank flatfile format.
Nucleotide and protein records may be mixed. Results will be returned in the order provided.

```julia
fetchseq("AH002844")                             # retrive one FASTA NCBI nucleotide record
fetchseq("CAA41295.1", "NP_000176"; format = gb) # retrieve two GenBank flatfile NCBI protein records
```
"""
function fetchseq(ids::AbstractString...; format::Format = fasta)
    ncbinucleotides = []
    ncbiproteins = []
    ncbinucleotidesorder = []
    ncbiproteinsorder = []
    for (i, id) ∈ enumerate(ids)
        ncbiprotein = startswith(id, r"[NX]P_|[A-Z]{3}[0-9]")
        ncbinucleotide = startswith(id, r"[NX][CGRMW]_|[A-Z]{2}[0-9]|[A-Z]{4,6}[0-9]")
        ncbiprotein && ncbinucleotide && throw("ambiguous identifier $id, cannot infer database")
        ncbiprotein || ncbinucleotide || throw("could not infer database for $id")

        if ncbinucleotide
            push!(ncbinucleotides, id)
            push!(ncbinucleotidesorder, i)
        elseif ncbiprotein
            push!(ncbiproteins, id)
            push!(ncbiproteinsorder, i)
        end
    end

    results = [fetchseq_ncbi(ncbinucleotides, "nuccore", format);
               fetchseq_ncbi(ncbiproteins, "protein", format)]

    return [results[i] for i ∈ [ncbinucleotidesorder; ncbiproteinsorder]]
end

function fetchseq_ncbi(ids::AbstractVector{<:AbstractString}, db::AbstractString, format::Format)
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

function fetchseq_ebi(ids::AbstractString...; format::Format = fasta)
    response = 
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