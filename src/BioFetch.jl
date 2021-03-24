module BioFetch

using FASTX
using GenomicAnnotations
using BioServices.EUtils

@enum Format fasta gb


function fetchseq(ids::AbstractString...; format::Format = :fasta)
    response = efetch(db = "nuccore", id = ids, rettype = String(format), retmode="text")
    body = IOBuffer(response.body)
    if format == :fasta
        reader = FASTA.Reader(body)
        records = [record for record ∈ reader]
    else
        records = readgbk(body)
    end
    return records
end



end