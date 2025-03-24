# DIMPLE quick QC app

This is a simple Shiny-python app to quickly assess a mutational variant library using whole-plasmid sequencing data. Although it is designed for our DIMPLE pipeline, it should be usable for many kinds of variant libraries.

## Description

Variant libraries, by design, contain a mixture of bases in mutated regions, but should be constant elsewhere. A consensus sequence generated from a library should show less certainty in the variant regions as a result. Although this does not provide variant-specific information, and assumes that variants are relatively close to the reference sequence, this will indicate general library quality. We provide the following data:

* Number of reads: the total number of reads at each position. This is a measure of the coverage at each position.

* Number of variant reads: the number of reads that are not the reference base (not including indels) at each position.

* Variant fraction: the fraction of reads that are not the reference base (not including indels) at each position.

* Insertion count: the number of insertions at each position. This may not be reliable due to the consensus calling pipeline.

* Deletion count: the number of deletions at each position.

* Indel fraction: the fraction of reads that are insertions or deletions at each position.

* Indel to substitution ratio: the ratio of indels to substitutions at each position.

* Entropy: the Shannon entropy of the base calls (but not indels) at each position. This is a measure of the uncertainty of the base calls at each position. A high entropy indicates a mix of bases, while a low entropy indicates a consensus.

* Effective entropy: to account for the fact that even in a variant region the reference base will be the most common, we calculate the entropy of the non-reference bases only. This should be much more sensitive, and the variation in this value across a mutated region should indicate the evenness of the variant distribution.

* Percent of max entropy: the maximum effective entropy at a position is simply ln(3), as there are three non-reference bases. This is the percentage of maximum possible entropy at each position.

* Max counts of non-reference base: the number of reads of the most common non-reference base at each position.

### Selecting specific ranges

Normally, only a portion of the reference sequence will be mutated. By entering a range, you can focus on a specific region of the sequence. The mean values of the selected and unselected ranges are then shown on the left.

This will also allow one more metric to be calculated:

* Variant fraction % expected: the expected fraction of variant reads at each position, assuming a random distribution of bases in the variants. This can indicate problems with large amounts of starting template in the library.

## Input

The app is designed to parse the results of whole-plasmid sequencing from [Plasmidsaurus](https://plasmidsaurus.com/), and thus uses a non-standard format. The app uses the per-base csv file provided in the "(sample name)_per-base-data" folder.

## Usage

This is meant to be deployed as a Shiny app. It can be installed locally using the requirements.txt file and running the app.py script. A basic Dockerfile is also provided.


## Caveats

This is intended to be a cost-effective and quick QC tool, not a comprehensive analysis of library. While it has been able to identify problematic libraries in our hands, you should always check the results manually.