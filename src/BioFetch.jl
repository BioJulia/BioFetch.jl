module BioFetch

using FASTX
using GenomicAnnotations
using BioServices.EUtils
using BioServices.EBIProteins

export fetchseq, SeqFormat, fasta, gb

@enum SeqFormat fasta gb

"""
    fetchseq(ids::AbstractString...; format::SeqFormat = fasta, range::AbstractUnitRange = nothing, revstrand = false)

Fetches sequence data from a database by accession number in either FASTA format or GenBank flatfile format.
Nucleotide and protein records may be mixed. Results will be returned in the order provided.
Supports NCBI, Ensembl, and UniProt accession numbers.

GenBank format may not work right now.

```julia
fetchseq("AH002844")                # retrive one FASTA NCBI nucleotide record
fetchseq(["CAA41295.1", "NP_000176"]) # retrieve two FASTA NCBI protein records
fetchseq("NC_036893.1", 81_775_230 .+ (1:1_000_000))         # retrive a 1 Mb segment of a FASTA NCBI genomic record
fetchseq("NC_036893.1", 81000000:81999999, true) # retrive a 1 Mb segment of a FASTA NCBI genomic record on the reverse strand
```

`range` and `revstrand` are only valid for nucleotide queries
"""
function fetchseq(ids::AbstractVector{<:AbstractString}, range::Union{Nothing, AbstractRange} = nothing, revstrand = false; format::SeqFormat = fasta)
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
                ncbiprotein ? ncbiproteins :
                ebiuniprot ? ebiuniprots :
                ebiensembl ? ebiensembls :
                error()
        push!(array, id)
        push!(order[array], i)
    end

    results = [fetchseq_ncbi(ncbinucleotides, "nuccore"; format, range, revstrand);
               fetchseq_ncbi(ncbiproteins, "protein"; format);
               fetchseq_uniprot(ebiuniprots; format);
               fetchseq_ensembl(ebiensembls; format)]
    order = [order[ncbinucleotides]; order[ncbiproteins]; order[ebiuniprots]; order[ebiensembls]]

    return results[order]
end

function fetchseq(id::AbstractString, range::Union{Nothing, AbstractRange} = nothing, revstrand = false; format::SeqFormat = fasta)
    ebiensembl = startswith(id, r"ENS[A-Z][0-9]{11}")
    ncbiprotein = startswith(id, r"[NX]P_|[A-Z]{3}[0-9]")
    ncbinucleotide = startswith(id, r"[NX][CGRMW]_|[A-Z]{2}[0-9]|[A-Z]{4,6}[0-9]") && !ebiensembl
    ebiuniprot = startswith(id, r"[A-Z][0-9][A-Z0-9]{4}")

    result = ncbinucleotide ? fetchseq_ncbi(id, "nuccore"; format, range, revstrand) :
           ncbiprotein ? fetchseq_ncbi(id, "protein"; format) :
           ebiuniprot ? fetchseq_uniprot([id]; format) :
           ebiensembl ? fetchseq_ensembl([id]; format) :
           error("could not infer database for $id")
    
    return result
end

function fetchseq_ncbi(ids, db::AbstractString; format::SeqFormat = fasta, range::Union{Nothing, AbstractRange} = nothing, revstrand = false)
    isempty(ids) && return []
    if isnothing(range)
        response = efetch(; db, id = ids, rettype = String(Symbol(format)), retmode="text", strand = revstrand ? 2 : 1)
    else
        response = efetch(; db, id = ids, rettype = String(Symbol(format)), retmode="text", seq_start = first(range), seq_stop = last(range), strand = revstrand ? 2 : 1)
    end
    body = IOBuffer(response.body)
    reader = format == fasta ? FASTA.Reader(body) :
             format == gb ? GenBank.Reader(body) :
             nothing
    records = [record for record ∈ reader]
    return records
end

function fetchseq_uniprot(ids::AbstractVector; format::SeqFormat = fasta)
    isempty(ids) && return []
    records = []
    contenttype = format == fasta ? "text/x-fasta" : format == gb ? "text/flatfile" : error("unknown format")
    for id ∈ ids
        response = ebiproteins(; accession = id, contenttype)
        body = IOBuffer(response.body)
        reader = format == fasta ? FASTA.Reader(body) :
                 format == gb ? GenBank.Reader(body) :
                 nothing
        records = [record for record ∈ reader]
        append!(records, records)
    end
    return records
end

function fetchseq_ensembl(ids::AbstractVector; format::SeqFormat = fasta)
    isempty(ids) && return []
    records = []
    contenttype = format == fasta ? "text/x-fasta" : format == gb ? "text/flatfile" : error("unknown format")
    for id ∈ ids
        response = ebiproteins(; dbtype = "Ensembl", dbid = id, contenttype)
        body = IOBuffer(response.body)
        reader = format == fasta ? FASTA.Reader(body) :
                 format == gb ? GenBank.Reader(body) :
                 nothing
        records = [record for record ∈ reader]
        append!(records, records)
    end
    return records
end

end