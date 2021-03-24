module BioFetch

using FASTX
using BioServices.EUtils


function fetchseq(ids::AbstractString...)
    response = efetch(db="nuccore", id=ids, rettype="fasta", retmode="text")
    reader = FASTA.Reader(IOBuffer(response.body))
    records = [record for record ∈ reader]
    return records
end



end