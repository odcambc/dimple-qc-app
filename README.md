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

## Installation

This is meant to be deployed as a Shiny app. It can be installed locally by either adding the dependencies to your environment manager of choice, or by using the provided Dockerfile, then running the app.py script.

### Install with uv (recommended)

```bash
uv sync
uv run shiny run
```

### Install with poetry

```bash
poetry install
poetry run shiny run
```

### Install with pip

```bash
pip install -r requirements.txt
python app.py
```

### Install with Docker

```bash
docker build -t dimple-qc-app .
docker run -p 8080:8080 dimple-qc-app
```

## Example workflow
We have found the following to work quite well for us
* After sub-pool cloning, we send the entire subpool to Plasmidsaurous for sequencing, rather than picking colonies in step 13.3 in the [dimple protocol](https://www.protocols.io/view/dimple-library-generation-and-assembly-protocol-rm7vzy7k8lx1/v6?step=8&version_warning=no)
* This can either be done via whole-plasmid sequencing of the subpool itself, or with premium PCR sequencing following a brief (5-10 cycle) PCR of the relevant portion of the subpool. Premium PCR sequencing will generate substantially more reads, which makes interpretation easier, but requires an additional PCR step.
* For an individual subpool, we expect the variant fraction to be around 0.12. For a properly constructed subpool, we don't expect to see any shift over the length of the subpool, and we don't expect any individual variants to be much higher than any others. It should look something like this: <img width="894" height="262" alt="image" src="https://github.com/user-attachments/assets/7f6b4da9-0442-44af-895f-ff791d861c25" />
* One advantage of this approach is that multiple subpools can be combined together and sequenced in the same reaction. You should see a proportional decrease in variant fraction as you add more subpools. e.g. expect an average around 0.06 for two subpools, and 0.04 for three subpools. As more subpools are added, the cost decreases, but it also becomes more challenging to distinguish signal from noise. <img width="925" height="314" alt="image" src="https://github.com/user-attachments/assets/c5f55cde-a518-49fb-91ea-09866510beef" />




## Caveats

This is intended to be a cost-effective and quick QC tool, not a comprehensive analysis of library. While it has been able to identify problematic libraries in our hands, you should always check the results manually.
