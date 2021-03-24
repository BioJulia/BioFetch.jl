module BioFetch

using Downloads
using FASTX

const ENTREZ_BASE_URL = raw"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

const ENTREZ_FETCH = "efetch.fcgi"

function fetch(ids::AbstractString...)
    url = ENTREZ_BASE_URL * ENTREZ_FETCH * "?db=nuccore&id=" * join(ids, ",") * "&rettype=fasta&retmode=text"
    path = Downloads.download(url)
    records = FASTA.Record[]
    open(FASTA.Reader, path) do reader
        record = FASTA.Record()
        while !eof(f)
            read!(f, record)
            push!(records, record)
        end
    end
    return records
end


end