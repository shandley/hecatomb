# Example config files

These example files are provided for your convenience. You will certailnly need to edit them, but they will clue you in on what you need.

Our database structure looks like this:

```
databases/
├── bac_giant_unique_species
│   ├── bac_uniquespecies_giant.masked_Ns_removed.fasta
│   └── ref
├── contaminants
│   ├── nebnext_adapters.fa
│   ├── primerB.fa
│   ├── rc_primerB_ad6.fa
│   ├── ref
│   └── vector_contaminants.fa.gz
└── human_masked
    ├── human_virus_masked.fasta
        └── ref
```

The *fastq* directory points to where your reads are!

Note that you can easily add a directory for (for e.g.) `mouse_masked` data, and then add that to the config file instead of the `human_masked` data.
