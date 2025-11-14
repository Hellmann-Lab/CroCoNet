# List of genomic annotations

List of genomic annotations per species as GRanges objects. The human
annotation is the GTF of Hg38 GENCODE release 32 primary assembly, while
the gorilla and cynomolgus macaque annotations were created by
transferring the human annotation onto the gorGor6 and macFas6 genomes
via the tool Liftoff (https://github.com/agshumate/Liftoff). Each
annotation was subsetted for the 300 genes that feature in this example
dataset.

## Usage

``` r
gtf_list
```

## Format

A named list of 3 GRanges objects. Each object has 8 columns.

Attributes:

- seqnames:

  Chromosome or contig name.

- start:

  Genomic start location.

- end:

  Genomic end location.

- width:

  Width of the feature in base pairs.

- strand:

  Genomic strand ("+" or "-").

- source:

  The prediction program or public database where the annotations came
  from.

- type:

  Feature type (gene, transcript, exon, CDS, UTR, start_codon,
  stop_codon or Selenocysteine).

- score:

  The degree of confidence in the feature's existence and coordinates.

- phase:

  One of '0', '1' or '2'. '0' means that the first base of the feature
  is the first base of a codon, '1' that the second base is the first
  base of a codon, and so on.

- gene_id:

  Unique identifier of the gene.

- gene_name:

  Name of the gene.

- transcript_id:

  Unique identifier of the transcript.

- transcript_name:

  Name of the transcript.

## Source

\<ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.gtf.gz\>,
\<https://hgdownload.soe.ucsc.edu/goldenPath/gorGor6/bigZips/gorGor6.fa.gz\>,
\<https://ftp.ensembl.org/pub/release-109/fasta/macaca_fascicularis/dna/Macaca_fascicularis.Macaca_fascicularis_6.0.dna_sm.toplevel.fa.gz\>
