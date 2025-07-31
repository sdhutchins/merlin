# MERLIN - Multipoint Engine for Rapid Likelihood Inference

<!-- markdown-link-check-disable -->
[![Build Status](https://img.shields.io/badge/build-passing-brightgreen)](https://github.com/your-org/merlin)
<!-- markdown-link-check-enable -->

**MERLIN** is a powerful tool for rapid analysis of dense genetic maps using sparse gene flow trees. This version is provided free of charge and comes with no guarantees.

**Copyright (c) 2000 - 2002 Goncalo Abecasis**

## Table of Contents

- [Project Background](#project-background)
- [Install & Setup](#install--setup)
- [Usage](#usage)
- [Examples](#examples)
- [Directory Structure](#directory-structure)
- [Contributing](#contributing)
- [License](#license)
- [Authors](#authors)

## Project Background

MERLIN (Multipoint Engine for Rapid Likelihood Inference) is designed for genetic linkage analysis and haplotype inference. It provides:

- **Non-parametric linkage analysis** using NPL pairs and NPL all scoring functions
- **Variance components analysis** for quantitative trait loci (QTL)
- **Haplotyping** with inheritance vector optimization
- **Simulation** of marker genotypes
- **Kinship matrix** probability calculations

The tool is particularly effective for analyzing dense genetic maps and can handle complex pedigree structures efficiently.

## Install & Setup

### Prerequisites

- **C++ Compiler**: GCC 4.8+ or compatible
- **Make**: GNU Make 3.8+
- **Memory**: 256MB+ RAM (configurable)

### Instructions

1. **Clone the repository**:
   ```bash
   git clone <repository-url>
   cd merlin
   ```

2. **Compile MERLIN**:
   ```bash
   make
   ```

3. **Verify installation**:
   ```bash
   ls executables/
   # Should show: merlin, pedstats, and other executables
   ```

## Usage

### Basic Command Structure

```bash
merlin -d <datafile> -p <pedfile> [options]
```

### Input File Formats

MERLIN supports two main input formats:

1. **Linkage format** - Standard pedigree analysis format
2. **QTDT format** - Quantitative trait analysis format

For QTDT format, you'll need a separate map file with 3 columns:
```
CHROMOSOME    MARKER     POSITION
1             D1S123     134.0
```

### Common Analysis Types

#### Basic Analyses
```bash
# Error checking
merlin -d data.dat -p pedigree.ped --error

# Information content
merlin -d data.dat -p pedigree.ped --info

# IBD and kinship matrices
merlin -d data.dat -p pedigree.ped --ibd --kinship
```

#### Linkage Analyses
```bash
# Non-parametric linkage (affecteds only)
merlin -d data.dat -p pedigree.ped --pairs
merlin -d data.dat -p pedigree.ped --npl

# QTL analysis
merlin -d data.dat -p pedigree.ped --qtl
merlin -d data.dat -p pedigree.ped --deviate
```

#### Variance Components
```bash
# Variance component linkage analysis
merlin -d data.dat -p pedigree.ped --vc
```

#### Haplotyping
```bash
# Find best inheritance vectors
merlin -d data.dat -p pedigree.ped --best

# Sample random inheritance vectors
merlin -d data.dat -p pedigree.ped --sample

# List all non-recombinant haplotypes
merlin -d data.dat -p pedigree.ped --all --zero --founders
```

## Examples

### Minimal Pedstats Example

The `pedstats` component provides pedigree statistics and visualization:

```bash
# Basic pedigree statistics
./executables/pedstats -d examples/basic2.dat -p examples/basic2.ped

# Hardy-Weinberg testing
./executables/pedstats -d examples/basic2.dat -p examples/basic2.ped --hardyWeinberg

# Generate PDF output
./executables/pedstats -d examples/basic2.dat -p examples/basic2.ped --pdf
```

### Advanced Examples

```bash
# Large dataset analysis
./executables/pedstats -d examples/assoc.dat -p examples/assoc.ped --verbose

# Cluster-based SNP analysis
merlin -d examples/snp-scan.dat -p examples/snp-scan.ped -m examples/snp-scan.map \
       --clusters examples/snp-scan.clusters

# Parametric linkage analysis
merlin -d examples/parametric.dat -p examples/parametric.ped -m examples/parametric.map \
       --model examples/parametric.model --freq examples/parametric.freq
```

## Directory Structure

```
merlin/
├── README.md                    <- This file
├── PEDSTATS_FIXES.md           <- C++11 compilation fixes documentation
├── Makefile                     <- Main build configuration
├── LICENSE.twister              <- Mersenne Twister license
├── merlin.pdf                   <- Detailed documentation
├── .gitignore                   <- Git ignore rules
├── examples/                    <- Example data files
│   ├── HapMap.genotypes        <- HapMap format genotype data
│   ├── HapMap.template         <- Pedigree template
│   ├── assoc.tbl               <- Association analysis table
│   ├── basic2.freq             <- Basic frequency file
│   ├── parametric.freq         <- Parametric frequency file
│   ├── parametric.model        <- Disease model specification
│   ├── snp-scan.clusters       <- SNP clustering data
│   └── snp-scan.clusters-only  <- Clusters-only data
├── libsrc/                      <- Core library source code
├── merlin/                      <- Main MERLIN source code
├── clusters/                    <- Clustering algorithms
├── regress/                     <- Regression analysis
├── pdf/                         <- PDF generation utilities
├── offline/                     <- Offline analysis tools
├── extras/                      <- Additional utilities
└── executables/                 <- Compiled binaries
```

## Resource Management

### Memory Control

```bash
# Limit pedigree complexity (bits)
merlin -d data.dat -p pedigree.ped --bits:16

# Restrict memory allocation (MB)
merlin -d data.dat -p pedigree.ped --megabytes:256

# Enable memory swapping
merlin -d data.dat -p pedigree.ped --swap
```

### Performance Tips

- Set `--megabytes` to available physical memory size
- Use `--bits` to skip overly complex pedigrees
- For large pedigrees on 32-bit systems, limit to 2048MB

## Contributing

We welcome contributions! Please:

1. Report bugs and provide helpful comments
2. Follow the existing code style and conventions
3. Test your changes thoroughly
4. Update documentation as needed

**Note**: This version includes some debug code. Performance improvements are planned for future releases.

## License

- **MERLIN**: Free for academic use
- **Mersenne Twister**: See `LICENSE.twister` for specific conditions

**Redistribution**: Redistribution of MERLIN in source or compiled formats is not allowed.

## Authors

- **Goncalo Abecasis** - Primary Developer
  - Email: goncalo@umich.edu
  - Institution: University of Michigan

## Citation

**GR Abecasis, SS Cherny, WOC Cookson and LR Cardon (2002)**
MERLIN - Rapid analysis of dense genetic maps using sparse gene flow trees.
*Nature Genetics* 30:97-101

## Support

- **Latest Version**: http://www.sph.umich.edu/csg/abecasis/Merlin
- **Tutorial**: http://www.sph.umich.edu/csg/abecasis/Merlin/tour
- **Reference**: http://www.sph.umich.edu/csg/abecasis/Merlin/reference.html

**Registration**: Please register by emailing goncalo@umich.edu or filling out the web registration form. 