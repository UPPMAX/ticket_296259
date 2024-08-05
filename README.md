# ticket_296259


Discuss:

- [ ] Goal: how to be able to fix the problem,
      maybe the VarTrix developers are better to ask
- [ ] Validate BAM file

```bash
module load bioinfo-tools
module load vartrix/1.1.22
module load picard/3.1.1
java -jar $PICARD ValidateSamFile --INPUT user_filename.bam
```

- [ ] What is the history of the BAM file? Goal: reproduce problem on public BAM file
- [ ] I've creates an Issue at VarTrix: [Help find problematic file with error 'Failed to seek to offset'](https://github.com/10XGenomics/vartrix/issues/124)
- [ ] Decide upon plan

### Solution

I predict there is something wrong with the bam input file.
The user thinks so to ([see communication of 2024-08-02](#2024-08-02)).
The bam input file cannot be shared.

However, it can be validated (code from [the UPPMAX documentation on Picard](https://docs.uppmax.uu.se/software/picard/)):

```bash
module load bioinfo-tools
module load vartrix/1.1.22
module load picard/3.1.1
java -jar $PICARD ValidateSamFile --INPUT $VARTRIX_TEST/test_dna.bam
```

The user tried to delete the analysis files
and started the analysis from scratch again,
but I am unsure if this was correct ([see communication of 2024-08-02](#2024-08-02))

- [ ] I've creates an Issue at VarTrix: [Help find problematic file with error 'Failed to seek to offset'](https://github.com/10XGenomics/vartrix/issues/124)

### Problem

> I am trying to run vartrix in Bianca.
> The program seems to work fine when I try the test run (files in $VARTRIX_TEST) but crashes when I use my input files.
> The errors I get:
> "Failed to seek to offset ((some number))"
> "Segmentation fault. Core dump" and I get a core.XXXX file. When I check the core.XXXX with gdb it says that it "Failed to read a valid object file image from memory" and "Program terminated with signal 6"
> I have no idea what it is trying to tell me. Could you maybe help me?
> (attached screenshots)

![](screenshot_2024_07_09_18_28.png)

![](screenshot_2024_07_09_18_38.png)

![](screenshot_2024_07_09_18_43.png)

## Notes

I assume this is about the Rust program called vartrix at [https://github.com/10XGenomics/vartrix](https://github.com/10XGenomics/vartrix).

A GitHub seach on `Failed to seek to offset`:

![](github_search.png)

 * [UNRESOLVED] [https://github.com/mflamand/Bullseye/issues/2#issuecomment-1067969449](https://github.com/mflamand/Bullseye/issues/2#issuecomment-1067969449)
 * [UNRESOLVED] [https://github.com/GMOD/jbrowse-components/issues/1990#issuecomment-844406099](https://github.com/GMOD/jbrowse-components/issues/1990#issuecomment-844406099)

I feel the problem is in [https://github.com/samtools/htslib](https://github.com/samtools/htslib), due to this PR:

 * [https://github.com/samtools/htslib/pull/1504](https://github.com/samtools/htslib/pull/1504)

Lines of code are at:

 * [https://github.com/samtools/htslib/blob/master/hts.c#L4182](https://github.com/samtools/htslib/blob/master/hts.c#L4182)
 * [https://github.com/samtools/htslib/blob/master/hts.c#L4203](https://github.com/samtools/htslib/blob/master/hts.c#L4203)

Both work on a `BGZF` pointer, a structure 
defined at [https://github.com/samtools/htslib/blob/develop/htslib/bgzf.h#L68](https://github.com/samtools/htslib/blob/develop/htslib/bgzf.h#L68).

As the user can get to work testing files, I feel it is most likely that the user's input files are in an invalid format.

I contacted the VarTrix maintainers to help me and the user diagnose the faulty file at [this Issue](https://github.com/10XGenomics/vartrix/issues/124).

## Communication



### 2024-08-02

From the user:

> The command I used to run a test in Vartrix:

```console
$ module load bioinfo-tools
$ module load vartrix/1.1.22
$ ls -ltr $VARTRIX_TEST/*
$ vartrix --bam $VARTRIX_TEST/test_dna.bam \
  --cell-barcodes $VARTRIX_TEST/dna_barcodes.tsv \
  --fasta $VARTRIX_TEST/test_dna.fa \
  --vcf $VARTRIX_TEST/test_dna.vcf
```

> Trying to remove all analysis files:
> there are none created, I checked with '''ls -altr''' and I see nothing
> 
> Input files:
> I agree in that it is, most likely, my input files.
> My input bams have been tinkered with, because my bam files are output coming from RSEM and they are missing the CB tag needed by Vartrix (as Vartrix takes bam files with CellRanger format)
> I used java/jvarkit to add the CB tag to my RSEM bams
> I also modified the vcf file to tailor it to what we need, but when I compare it to the test vcf file they seem to have the same format

> Sharing my inputs:
> I'm not sure I can share my bam as it is sensitive information, but I have it in Bianca, so maybe you can access that directory?
> For the remaining files I could send you a link to a tar in OneDrive, would that work?
```
