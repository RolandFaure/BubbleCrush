# BubbleCrush
Crushes simple bubbles in GFA

## Usage
```
usage: bubble_crush.py [-h] -i INPUT_ASSEMBLY -o OUTPUT_ASSEMBLY
                       [-a ABSOLUTE_THRESHOLD] [-r RELATIVE_THRESHOLD] [-m]
                       [-t]

Parse a genomic assembly in GFA format and pop the bubbles

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_ASSEMBLY, --input-assembly INPUT_ASSEMBLY
                        Original assembly in GFA format (required)
  -o OUTPUT_ASSEMBLY, --output-assembly OUTPUT_ASSEMBLY
                        Output assembly with bubbles crushed (required)
  -a ABSOLUTE_THRESHOLD, --absolute-threshold ABSOLUTE_THRESHOLD
                        Branch of a bubble less covered than this will be
                        deleted. 0 to turn off. [5]
  -r RELATIVE_THRESHOLD, --relative-threshold RELATIVE_THRESHOLD
                        Branch of a bubble less covered than this compared to
                        other will be deleted. 0 to turn off. 1 for crushing
                        the less covered branch. [0]
  -m, --dont_merge      Do not merge the contigs after bubble popping
  -t, --transfer-coverage
                        Transfer the coverage of the deleted branch to the
                        other branch. [False]
```
