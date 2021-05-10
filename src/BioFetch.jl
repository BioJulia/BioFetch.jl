module BioFetch

using FASTX
using GenomicAnnotations
using BioServices.EUtils

export fetchseq

@enum Format fasta gb

"""
    fetchseq(ids::AbstractString...; format::Format = fasta)

Fetches sequence data from a database by accession number in either FASTA format or GenBank flatfile format.
Nucleotide and protein records may be mixed. Results will be returned in the order provided.
Supports NCBI, Ensembl, and UniProt accession numbers.

```julia
fetchseq("AH002844")                             # retrive one FASTA NCBI nucleotide record
fetchseq("CAA41295.1", "NP_000176"; format = gb) # retrieve two GenBank flatfile NCBI protein records
```
"""
function fetchseq(ids::AbstractString...; format::Format = fasta)
    ncbinucleotides = String[]
    ncbiproteins = String[]
    ebiuniprots = String[]
    ebiensembls = String[]
    order = IdDict(ncbinucleotides => [], ncbiproteins => [], ebiuniprots => [], ebiensembls => [])
    for (i, id) ∈ enumerate(ids)
        ebiensembl = startswith(id, r"ENS[A-Z][0-9]{11}")
        ncbiprotein = startswith(id, r"[NX]P_|[A-Z]{3}[0-9]")
        ncbinucleotide = startswith(id, r"[NX][CGRMW]_|[A-Z]{2}[0-9]|[A-Z]{4,6}[0-9]") && !ebiensembl
        ebiuniprot = startswith(id, r"[A-Z][0-9][A-Z0-9]{4}")
        ncbiprotein || ncbinucleotide || ebiuniprot || ebiensembl || throw("could not infer database for $id")

        array = ncbinucleotide ? ncbinucleotides :
                ncbiproteins ? ncbiproteins :
                ebiuniprot ? ebiuniprots :
                ebiensembl ? ebiensembls :
                error()
        push!(array, id)
        push!(order[array], id)
    end

    results = [fetchseq_ncbi(ncbinucleotides, "nuccore"; format);
               fetchseq_ncbi(ncbiproteins, "protein"; format);
               fetchseq_uniprot(ebiuniprots; format);
               fetchseq_ensembl(ebiensembls; format)]
    order = [order[ncbinucleotides]; order[ncbiproteins]; order[ebiuniprots]; order[ebiensembls]]
            
    return results[order]
end

function fetchseq_ncbi(ids::AbstractVector{<:AbstractString}, db::AbstractString; format::Format = fasta)
    isempty(ids) && return []
    response = efetch(; db, id = ids, rettype = String(Symbol(format)), retmode="text")
    println(response.body)
    body = IOBuffer(response.body)
    if format == fasta
        reader = FASTA.Reader(body)
        records = [record for record ∈ reader]
    else
        records = readgbk(body)
    end
    return records
end

function fetchseq_uniprot(ids::AbstractVector{<:AbstractString}; format::Format = fasta)
    isempty(ids) && return []
    records = []
    contenttype = format == fasta ? "text/x-fasta" : format == gb ? "text/flatfile" : error("unknown format")
    for id ∈ ids
        response = ebiproteins(; accession = id, contenttype)
        body = IOBuffer(response.body)
        if format == fasta
            reader = FASTA.Reader(body)
            records = first(reader)
        else
            records = readgbk(body)
        end
        push!(records, record)
    end
    return records
end

function fetchseq_ensembl(ids::AbstractVector{<:AbstractString}; format::Format = fasta)
    isempty(ids) && return []
    records = []
    contenttype = format == fasta ? "text/x-fasta" : format == gb ? "text/flatfile" : error("unknown format")
    for id ∈ ids
        response = ebiproteins(; dbtype = "Ensembl", dbid = id, contenttype)
        body = IOBuffer(response.body)
        if format == fasta
            reader = FASTA.Reader(body)
            records = first(reader)
        else
            records = readgbk(body)
        end
        push!(records, record)
    end
    return records
end

end