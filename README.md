# ticket_296259

## Solution

I predict there is something wrong with the input files.
These could be of invalid type or incorrectly called.

- Could you share the input files?

As I assume the data is private, I'll ask:

- Could you share the command and scripts you used to call vatrix for the testing data?
- Could you share the command and scripts you used to call vatrix for the real data?

Besides that, I need to know

- What happens if you delete the analysis files (keeping the data!) 
  and start from scratch again? It seems like things are out of sync

## Problem

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
