# snpeff.sh
Script to prepare file for input to snpeff, and then run snpeff, select sites of interest (0/4 fold degenrate, UTR etc).

Briefly first part of the script edits VCF file to have NAs specified as the alternate allele for *all* sites (variant and invarint).
Then when you run SNPeff on this file, it spits out annotations on all possile substitutions for this site.

The script then matches these annotations for each site to pull sites of interest. 

