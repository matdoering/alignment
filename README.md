# alignment
Alignment of genomic sequences using the Needleman-Wunsch algorithm (global alignment) or the Smith-Waterman algorithm (local alignment).

## How to build and run?

To build: Just call `make`

For a global DNA alignment:

 ./Alignment TACTACCAGA TACTAGGGGG

For a global RNA alignment:

 ./Alignment UACUAGCCAGA UACUAGGGG

For a local alignment:
 ./Alignment -l TAC AAAAAAAAAAGGGTACTTTT


