# ticket_296259


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