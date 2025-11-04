# PepQueryMHC

---
- [About](#about)
- [Usage](#usage)
  - [Quick start](#quick-start)
  - [Parameters](#parameters)
  - [Scan mode](#scan-mode)
  - [Target mode](#target-mode)
  - [FASTQ mode](#fastq-mode)
  - [Annotate mode](#annotate-mode)
  - [White list](#white-list)
- [Output](#output)
  - [scan.tsv](#scan.tsv)
  - [target.tsv](#target.tsv)
  - [gloc.tsv](#gloc.tsv)
  - [peptide.tsv](#peptide.tsv)
  - [annotate.tsv](#annotate.tsv)
- [Figures for the reviewers](#figures-for-the-reviewers)
- [License]($license)
---

## About

The accurate prioritization of tumor antigens, including aberrant translational products, is critical for the development of personalized cancer immunotherapies. PepQueryMHC estimates a comprehensive repertoire of local RNA expression of tumor antigens within minutes per sample.
<br>

## Usage
PepQueryMHC provides three main functions such as 1) scan mode, 2) target mode 3) FASTQ mode and 4) annotate mode. <br>
When you use FASTQ mode, please make sure that FASTQ files do not contain artifical sequences such as adaptors, barcodes and so on.

### Quick start
Scan mode
```bash
java -Xmx2G -jar PepQueryMHC.jar \
--mode scan \
--input peptides.tsv \
--bam sample.sorted.bam \
--output sample \
--thread 16
```
Target mode
```bash
java -Xmx2G -jar PepQueryMHC.jar \
--mode target \
--input peptides_locations_strands.tsv \
--bam sample.sorted.bam \
--output sample \
--thread 16
```
FASTQ mode (for single-end)
```bash
java -Xmx2G -jar PepQueryMHC.jar \
--mode fastq \
--input peptides.tsv \
--0 sample.trimmed.fastq.gz \
--output sample \
--strand f \
--thread 16
```
FASTQ mode (for paired-end)
```bash
java -Xmx2G -jar PepQueryMHC.jar \
--mode fastq \
--input peptides.tsv \
--1 sample.trimmed.fastq.1.gz \
--2 sample.trimmed.fastq.2.gz \
--output sample \
--strand rf \
--thread 16
```
Annotate mode
```bash
java -Xmx2G -jar PepQueryMHC.jar \
--mode annotate \
--input locations_strands.tsv \
--gtf reference_annotation.gtf \
--output sample
```

### Parameters
Y+: mandatory, Y: optional, N: none
|Option    | Description    | Type   | Default | Scan mode   | Target mode   | FASTQ mode | Annotate mode   |
| :---:    | :---:          | :---:   | :---:       | :---:       | :---:         | :---:           | :---:           |
| m/mode   | mode to use| scan\|target\|fastq\|annotate  | | Y+          | Y+            | Y+              | Y+              |
| i/input  | input file path| string  || Y+          | Y+            | Y+             | Y+             |
| o/output  | output base name path| string  || Y+          | Y+           | Y+             | Y+             |
| b/bam  | sorted bam/sam file path | bam\|sam  || Y+          | Y+            | N              | N              |
| 0/fastq_single  | fastq file path | fastq\|fastq.gz  || N          | N            | Y+              | N              |
| 1/fastq_paired_1  | fastq file path | fastq\|fastq.gz  || N          | N            | Y+              | N              |
| 2/fastq_paried_2  | fastq file path | fastq\|fastq.gz  || N          | N            | Y+              | N              |
| g/gtf  | gtf file path | string  || N          | N            | N            | Y+              |
| @/thread  | the number of threads | int  |4| Y          | Y            | Y            | N              |
| c/count  | tpye of reads being processed | primary\|all  | primary | Y          | Y            | N              | N              |
| l/lib_size  | tsv file including library size information | string |  | Y          | Y            | Y            | N              |
| w/white_list  | cell brcode list (tsv), only available in single-cell RNA-seq | string |  | Y          | Y            | Y            | N              |
| p/prob  | ignore region of interests with error > p| [0,1] | 0.05 | Y          | Y           | Y             | N              |
| e/equal  | specify isoleucine = leucine | none |  | Y          | Y           | Y            | N              |
| u/union  | specify the unit of the peptide read count | sum\|max | sum | Y          | Y            | N              | N              |
| s/strand  | specify strandedness. non: non-stranded, fr: fr-second strand, rf: fr-first strand, f: forward strand for single-end, r: reverse strand for single-end, auto: auto-detection. Auto-detection is only available if there is XS tag in a given bam file | non\|fr\|rf\|f\|r\|auto | auto | Y          | Y             | Y+            | N              |
| s/stretch  | output single line per annotation | none |  | N          | N            | N            | Y              |
| v/verbose  | print every messages being processed | none |  | Y          | Y            | Y              | Y              |

### Scan mode
**Input format**
|Sequence| User-defined column 1| ...   | User-defined column N |
| :---:    | :---:          | :---:   | :---:   |
|AACTKLAKKM| any value | ... | any value |


### Target mode
**Input format**
|Sequence| Location | Strand |User-defined column 1| ...   | User-defined column N |
| :---:    | :---:          | :---:          | :---:          | :---:   | :---:   |
|AACTKLAKKM| chr1:1-30 | + | any value | ... | any value |
|TKMQEPPALY| chr1:31-50\|chr1:81-90 | - | any value | ... | any value |
|KEKRKAPPR| . | . |  any value | ... | any value |

### FASTQ mode
**Input format**
|Sequence| User-defined column 1| ...   | User-defined column N |
| :---:    | :---:          | :---:   | :---:   |
|AACTKLAKKM| any value | ... | any value |
* input format is exactly the same as what used in scan mode.

### Annotate mode
**Input format**
| Location | Strand |User-defined column 1| ...   | User-defined column N |
| :---:          | :---:          | :---:          | :---:   | :---:   |
| chr1:1-30 | + |  any value | ... | any value |
| chr1:31-50\|chr1:81-90 | -  | any value | ... | any value |
| chr1:21-40\|chr1:87-90 | .  | any value | ... | any value |

### White list
A white list is a set of barcodes selected for inclusion in the analysis of single-cell RNA-seq data.<br>
**Input format**
| Barcode | 
| :---:   |
| AAACCTGAGCAATCTC-1 |
| AAACCTGAGCGTTTAC-1 |
| AAACCTGAGCTGCAAG-1 |
| AAACCTGCAAACGTGG-1 |
| AAACCTGCAAACTGCT-1 |
| AAACCTGCAACTGGCC-1 |
| AAACCTGCAAGCCCAC-1 |
| AAACCTGCACTTAAGC-1 |
| AAACCTGCAGCCACCA-1 |

## Output
PepQueryMHC provides four modes (scan, fastq, target, and annotate modes) and generates different levels of transcriptomic annotations.

### scan.tsv
This format is generated in both scan and fastq modes and contains detailed level of transcriptomic features.<br>
| Column                            | Description                                                                       | Example                     |  Scan mode   | FASTQ mode |
| :---:                             | :---:                                                                             | :---:                       | :---:        | :---:      |
| Matched_location                  | Matched genomic location                                                          | chr16:84101848-84101874     |  Y           | N          |
| Matched_mutations                 | Nucleotide variants including SNV and INDEL                                       | chr16:84101848A>T           |  Y           | N          |
| Matched_strand                    | RNA strand                                                                        | -                           |  Y           | N          |
| Matched_peptide                   | Translated nucleotide sequence directly from the input NGS                        | LLAETKIHL                   |  Y           | Y          |
| Matched_nucleotide                | Nucleotide sequence directly from the input NGS                                   | TAAGTGAATTTTTGTTTCAGCTAAAAG |  Y           | Y          |
| Matched_reference_nucleotide      | Reference nucleotide sequence of a given genomic region                           | aAAGTGAATTTTTGTTTCAGCTAAAAG |  Y           | N          |
| Matched_read_count                | The number of reads supporting the match                                          | 91                          |  Y           | Y          |
| Matched_RPHM                      | Normalized read count supporting tha match                                        | 32.6392389935464            |  Y           | Y          |
| Proportion                        | Proportion of this annotation among all reads matching the input peptide sequence | 0.947916666666666           |  Y           | Y          |

### target.tsv
This format is generated in target mode.<br>
| Column                            | Description                                                                       | Example                     |
| :---:                             | :---:                                                                             | :---:                       |
| Matched_location                  | Matched genomic location                                                          | chr16:84101848-84101874     |
| Matched_mutations                 | Nucleotide variants including SNV and INDEL                                       | chr16:84101848A>T           |
| Matched_strand                    | RNA strand                                                                        | -                           |
| Matched_peptide                   | Translated nucleotide sequence directly from the input NGS                        | LLAETKIHL                   |
| Matched_nucleotide                | Nucleotide sequence directly from the input NGS                                   | TAAGTGAATTTTTGTTTCAGCTAAAAG |
| Matched_reference_nucleotide      | Reference nucleotide sequence of a given genomic region                           | aAAGTGAATTTTTGTTTCAGCTAAAAG |
| Matched_read_count                | The number of reads supporting the match                                          | 91                          |
| Matched_RPHM                      | Normalized read count supporting tha match                                        | 32.6392389935464            |
| Proportion                        | Proportion of this annotation among all reads matching the input peptide sequence | 0.978494623655914           |


### gloc.tsv
This format is generated in both scan and target modes. It summarizes matches at the genomic location level, regardless of mutations<br>
| Column                            | Description                                                                       | Example                     |  Scan mode   | Target mode |
| :---:                             | :---:                                                                             | :---:                       | :---:        | :---:       |
| Matched_peptide                   | Translated nucleotide sequence directly from the input NGS                        | LLAETKIHL                   |  Y           | Y           |
| Matched_location                  | Matched genomic location                                                          | chr16:84101848-84101874     |  Y           | Y           |
| Matched_strand                    | RNA strand                                                                        | -                           |  Y           | Y           |
| Matched_read_count                | The number of reads supporting the match regardless of mutations                  | 93                          |  Y           | Y           |
| Matched_RPHM                      | Normalized read count supporting tha match regradless of mutations                | 33.3565849054925            |  Y           | Y           |
| Proportion                        | Proportion of this annotation among all reads matching the input peptide sequence | 0.96875                     |  Y           | Y           |

### peptide.tsv
This format is generated in both scan and target modes. It summarizes matches at the peptide level.<br>
| Column                            | Description                                                                       | Example                     |  Scan mode   | Target mode |
| :---:                             | :---:                                                                             | :---:                       | :---:        | :---:       |
| Matched_peptide(sum)              | Translated nucleotide sequence directly from the input NGS                        | LLAETKIHL                   |  Y           | Y           |
| Abundant_location                 | The most abundant genomic location among all matched locations                    | chr16:84101848-84101874     |  Y           | Y           |
| Abundant_strand                   | RNA strand of the most abundant genomic location                                  | -                           |  Y           | Y           |
| Proportion                        | Proportion of this annotation among all reads matching the input peptide sequence | 0.96875                     |  Y           | Y           |
| Matched_num_locations             | The number of matched genomic locations                                           | 2                           |  Y           | Y           |
| Matched_read_count                | Sum of the number of reads supporting the matches                                 | 96                          |  Y           | Y           |
| Matched_RPHM                      | Normalized read count supporting tha sum of matches                               | 34.4326037734116            |  Y           | Y           |


### annotate.tsv
This format is generated in annotate mode.<br>
| Column                            | Description                                                                       | Example                     |
| :---:                             | :---:                                                                             | :---:                       |
| Annotation_count                  | The number of possible annotations based on a given gene model                    | 1                           |
| Gene_id                           | Gene id of a given genomic location                                               | ENSG00000140943.18          |
| Gene_name                         | Gene name of a given genomic location                                             | MBTPS1                      |
| Gene_strand                       | Gene strand                                                                       | -                           |
| Gene_type                         | Gene type                                                                         | protein_coding              |
| Class_code                        | It assigns events, detailed in below                                              | 5`-UTR                      |
| Unique_class_code                 | It only reports unique class code by removing duplications                        | 5`-UTR                      |
| Warning_tag                       | It reports a warning if the gene annotation is incomplete                         | .                           |

Class_code<br>
-- In-frame (IF) <br>
-- Out-of-frame (OOF) <br>
-- Non-coding RNA (ncRNA) <br>
-- 5\`-untranslated-region (5\`-UTR) <br>
-- 3\`-untranslated-region (3\`-UTR) <br>
-- Intron-retention (IR) <br>
-- Antisense RNA (asRNA) <br>
-- Intergenic-region (IGR) <br>
-- Exon-skipping (ES) <br>
-- Alternative-splicing-site (ASS) <br>
-- Unknown <br>


## Figures for the reviewers
Figures in the paper can be generated using: 1) R code in figR folder, 2) Supplementary Tables 1,2,3, and 5, and 3) the meta dataset available at https://doi.org/10.5281/zenodo.17429717. 

## License
All code is available as under the <a href="https://creativecommons.org/licenses/by-nc/4.0/">Attribution-NonCommercial (CC BY-NC) 4.0 license</a>.