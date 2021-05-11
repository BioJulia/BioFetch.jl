using BioFetch, FASTX, GenomicAnnotations, Test

nr002726 = fetchseq("NR_002726.2")
@test length(nr002726) == 1
@test typeof(nr002726[1]) == FASTA.Record

nr002726gb = fetchseq("NR_002726.2", format = gb)
@test length(nr002726gb) == 1
@test typeof(nr002726gb[1]) == GenomicAnnotations.Record{Gene}

seqs = fetchseq("CAA41295.1", "NP_000176", format = fasta)
@test length(seqs) == 2
@test typeof(seqs[1]) == FASTA.Record

@test length(fetchseq("ENSG00000141510")) == 42
