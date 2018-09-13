```Canvas_master.sh``` runs Canvas on both tumor and normal WGS samples.
What is yet unclear is how to remove germline calls that pop-up in the somatic VCF, see:
```https://github.com/Illumina/canvas/issues/67``` &
```https://github.com/Illumina/canvas/issues/60```

The ```filter13```exclusion file in the ```Canvas_master.sh``` contains the centromeres as provided by Illumina (see Canvas GitHub page) plus all non-primary chromosomes of GRCh38 (```grep -E 'chrU|_random|chrM'```).
