#!/bin/bash

################################################################################################
################################################################################################

## Get the files required for Canvas, as instructed on the manual page:

################################################################################################
################################################################################################

## ChromSizes:
wget http://canvas-cnv-public.s3.amazonaws.com/GRCh38/WholeGenomeFasta/GenomeSize.xml

## genome:
wget http://canvas-cnv-public.s3.amazonaws.com/GRCh38/WholeGenomeFasta/genome.fa

## index:
wget http://canvas-cnv-public.s3.amazonaws.com/GRCh38/WholeGenomeFasta/genome.fa.fai

## kmers:
wget http://canvas-cnv-public.s3.amazonaws.com/GRCh38/kmer.fa
