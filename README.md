# ticket_296259

Abbreviation | In full
-------------|---------------------
`BAM`        | Binary Alignment Map
`CIGAR`      | Compact Idiosyncratic Gapped Alignment Report
`SAM`        | Sequence Alignment Map

## 2024-08-06

Bash script from user [run.sh](run.sh).

The problem is in an early line:

```bash
samtools view -b ${bamO} chr2 > SS2_19_037-H13_chr2.bam
```

This is a problem as a samtools view is not in the same format as a BAM
file. Using this line instead brings us closer to the truth:

```bash
samtools view -b ${bamO} chr2 > SS2_19_037-H13_chr2.sam
```

This is because a SAM file is a plaintext format that samtools
renders to. However, the user does use the samtools `-b` flag,
which, according to [the documentation](https://github.com/samtools/samtools/blob/master/doc/samtools-view.1)
states `-b, --bam Output in the BAM format`

I hypothesize that `-b` is ignored when an output file is not 
specified (with the `-o` flag),
which can be tested by seeing if this would indeed create a BAM file:

```bash
samtools view --bam ${bamO} chr2 --output SS2_19_037-H13_chr2.bam
```

- [ ] Test hypothesis

This is, however, not the solution. Instead, I predict this will work:

```bash
samtools view ${bamO} chr2 > SS2_19_037-H13_chr2.sam
# Work on the SAM file
samtools view SS2_19_037-H13_chr2.sam --bam --output SS2_19_037-H13_chr2.bam
```

- [ ] Test hypothesis

User text on how to get a FASTA:

> Fasta (too large to send as attachment) downloaded from https://www.gencodegenes.org/human/release_37.html (Fasta files > genome sequence)

Indeed, the FASTA file can be downloaded from [https://www.gencodegenes.org/human/release_37.html](https://www.gencodegenes.org/human/release_37.html),
which links to [https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/GRCh38.p13.genome.fa.gz](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/GRCh38.p13.genome.fa.gz).
I put that in a script [get_fasta.sh](get_fasta.sh).

- [ ] Verify `SS2_19_037-H13_chr2.bam` is a valid BAM file
    - `cat SS2_19_037-H13_chr2.bam` indeed shows a binary file
    - `samtools quickcheck SS2_19_037-H13_chr2.bam` gives no error
    - `samtools quickcheck -vvv SS2_19_037-H13_chr2.bam` gives no error and some good output:

```
richel@richel-N141CU:~/GitHubs/ticket_296259$ samtools quickcheck -vvv SS2_19_037-H13_chr2.bam
verbosity set to 3
checking SS2_19_037-H13_chr2.bam
opened SS2_19_037-H13_chr2.bam
SS2_19_037-H13_chr2.bam is sequence data
SS2_19_037-H13_chr2.bam has 116 targets in header.
SS2_19_037-H13_chr2.bam has good EOF block.
```

- [ ] Run use script, reproduce problem

Does not work yet:

```bash
[richel@rackham1 ticket_296259]$ cat slurm-49075340.out
[DEBUG][SamJdk] Compiling :
         1  import java.util.*;
         2  import java.util.stream.*;
         3  import java.util.function.*;
         4  import htsjdk.samtools.*;
         5  import htsjdk.samtools.util.*;
         6  public class SamJdkCustom173338094 extends com.github.lindenb.jvarkit.tools.samjs.SamJdk.AbstractFilter {
         7    public SamJdkCustom173338094(final SAMFileHeader header) {
         8    super(header);
         9    }
        10    @Override
        11    public Object apply(final SAMRecord record) {
        12     /** user's code starts here */
        13  String c=record.getReadName(); int h=0; int s=21; record.setAttribute(CB,c.substring(h,s));return record;
        14  /** user's code ends here */
        15     }
        16  }
/tmp/jvarkit.tmp3660225531974008139/SamJdkCustom173338094.java:13: error: cannot find symbol
String c=record.getReadName(); int h=0; int s=21; record.setAttribute(CB,c.substring(h,s));return record;
                                                                      ^
  symbol:   variable CB
  location: class SamJdkCustom173338094
1 error
[SEVERE][SamJdk]java.lang.RuntimeException: Cannot compile
java.lang.RuntimeException: java.lang.RuntimeException: Cannot compile
	at com.github.lindenb.jvarkit.lang.OpenJdkCompiler$DefaultOpenJdkCompiler.compileClass(OpenJdkCompiler.java:250)
	at com.github.lindenb.jvarkit.tools.samjs.SamJdk.doWork(SamJdk.java:584)
	at com.github.lindenb.jvarkit.util.jcommander.Launcher.instanceMain(Launcher.java:819)
	at com.github.lindenb.jvarkit.util.jcommander.Launcher.instanceMainWithExit(Launcher.java:982)
	at com.github.lindenb.jvarkit.tools.samjs.SamJdk.main(SamJdk.java:809)
	at java.base/jdk.internal.reflect.NativeMethodAccessorImpl.invoke0(Native Method)
	at java.base/jdk.internal.reflect.NativeMethodAccessorImpl.invoke(NativeMethodAccessorImpl.java:77)
	at java.base/jdk.internal.reflect.DelegatingMethodAccessorImpl.invoke(DelegatingMethodAccessorImpl.java:43)
	at java.base/java.lang.reflect.Method.invoke(Method.java:568)
	at com.github.lindenb.jvarkit.tools.jvarkit.JvarkitCentral$Command.execute(JvarkitCentral.java:260)
	at com.github.lindenb.jvarkit.tools.jvarkit.JvarkitCentral.run(JvarkitCentral.java:716)
	at com.github.lindenb.jvarkit.tools.jvarkit.JvarkitCentral.main(JvarkitCentral.java:727)
Caused by: java.lang.RuntimeException: Cannot compile
	at com.github.lindenb.jvarkit.lang.OpenJdkCompiler$DefaultOpenJdkCompiler.exec(OpenJdkCompiler.java:164)
	at com.github.lindenb.jvarkit.lang.OpenJdkCompiler$DefaultOpenJdkCompiler.compileClass(OpenJdkCompiler.java:227)
	... 11 more
[INFO][Launcher]samjdk Exited with failure (-1)
[main_samview] fail to read the header from "./SS2_19_037-H13-CB_chr2.sam".
samtools index: "SS2_19_037-H13_chr2_bam.bam" is in a format that cannot be usefully indexed
```

Problem seems that the text `CB` disappears. Let's fix it by `"CB"` to `\"CB\"`.

Works:

```
[richel@rackham1 ticket_296259]$ cat slurm-49075347.out
[DEBUG][SamJdk] Compiling :
         1  import java.util.*;
         2  import java.util.stream.*;
         3  import java.util.function.*;
         4  import htsjdk.samtools.*;
         5  import htsjdk.samtools.util.*;
         6  public class SamJdkCustom1568608340 extends com.github.lindenb.jvarkit.tools.samjs.SamJdk.AbstractFilter {
         7    public SamJdkCustom1568608340(final SAMFileHeader header) {
         8    super(header);
         9    }
        10    @Override
        11    public Object apply(final SAMRecord record) {
        12     /** user's code starts here */
        13  String c=record.getReadName(); int h=0; int s=21; record.setAttribute("CB",c.substring(h,s));return record;
        14  /** user's code ends here */
        15     }
        16  }
[richel@rackham1 ticket_296259]$ 
```

- [x] [FAILS] Confirm that `samtools quickcheck` on the `SS2_19_037-H13_chr2.bam` created by the original script
  indicates a problem

```bash
[richel@rackham1 ticket_296259]$ samtools quickcheck SS2_19_037-H13_chr2.bam
[richel@rackham1 ticket_296259]$ 
```

- [x] Send a feature request to samtools to detect this. Done at [https://github.com/samtools/samtools/issues/2094](https://github.com/samtools/samtools/issues/2094)

`samtools` determines the type of file from the input file, not from its filename

- [x] Confirm that `java -jar $PICARD ValidateSamFile --INPUT java -jar $PICARD ValidateSamFile --INPUT user_filename.bam` 
      indicates a problem

```bash
[richel@rackham1 ticket_296259]$ java -jar $PICARD ValidateSamFile --INPUT SS2_19_037-H13_chr2.bam 
Aug 06, 2024 1:14:05 PM com.intel.gkl.NativeLibraryLoader load
INFO: Loading libgkl_compression.so from jar:file:/sw/bioinfo/picard/3.1.1/rackham/picard.jar!/com/intel/gkl/native/libgkl_compression.so
[Tue Aug 06 13:14:05 CEST 2024] ValidateSamFile --INPUT SS2_19_037-H13_chr2.bam --MODE VERBOSE --MAX_OUTPUT 100 --IGNORE_WARNINGS false --VALIDATE_INDEX true --INDEX_VALIDATION_STRINGENCY EXHAUSTIVE --IS_BISULFITE_SEQUENCED false --MAX_OPEN_TEMP_FILES 8000 --SKIP_MATE_VALIDATION false --VERBOSITY INFO --QUIET false --VALIDATION_STRINGENCY STRICT --COMPRESSION_LEVEL 5 --MAX_RECORDS_IN_RAM 500000 --CREATE_INDEX false --CREATE_MD5_FILE false --help false --version false --showHidden false --USE_JDK_DEFLATER false --USE_JDK_INFLATER false
[Tue Aug 06 13:14:05 CEST 2024] Executing as richel@rackham1.uppmax.uu.se on Linux 3.10.0-1160.119.1.el7.x86_64 amd64; OpenJDK 64-Bit Server VM 17+35-2724; Deflater: Intel; Inflater: Intel; Provider GCS is available; Picard version: Version:3.1.1
WARNING	2024-08-06 13:14:05	ValidateSamFile	NM validation cannot be performed without the reference. All other validations will still occur.
ERROR::MISSING_READ_GROUP:Read groups is empty
WARNING::RECORD_MISSING_READ_GROUP:Read name Run0484_ACDWKAANXX_L4_R1_T2103_C1510673, A record is missing a read group
WARNING::MISSING_TAG_NM:Record 1, Read name Run0484_ACDWKAANXX_L4_R1_T2103_C1510673, NM tag (nucleotide differences) is missing
WARNING::RECORD_MISSING_READ_GROUP:Read name Run0484_ACDWKAANXX_L4_R1_T1209_C1324596, A record is missing a read group
[...]
```

## 2024-08-05

### Meeting Monday 2024-08-05 10:00

Discuss:

- [x] Goal: how to be able to fix the problem, maybe the VarTrix developers are better to ask

It seems fixable: we concluded that the user's script converted a BAM file
to a text file. I can fix that script :-)

- [x] Validate BAM file

```bash
module load bioinfo-tools
module load vartrix/1.1.22
module load picard/3.1.1
java -jar $PICARD ValidateSamFile --INPUT user_filename.bam
```

We found out the BAM file was not a binary file at all!
Instead, it was a text file in the same format as `samtools view` gives.

- [/] What is the history of the BAM file? Goal: reproduce problem on public BAM file

The user will send me her script

- [x] I've creates an Issue at VarTrix: [Help find problematic file with error 'Failed to seek to offset'](https://github.com/10XGenomics/vartrix/issues/124)

Not needed.

- [x] Decide upon plan

    - User sends script and subsetted BAM file before 14:00
    - Richel works on the fix
    - Next meeting: Tuesday Aug 6th 14:00

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

## 2024-08-02

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

### Notes

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

### Communication

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


