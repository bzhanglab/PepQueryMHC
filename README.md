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
- [Paper](#paper)
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

## Paper
Figures in the paper can be generated using: 1) R code in figR folder, 2) Supplementary Tables S1,S2,S3, and S5, and 3) the meta dataset available at https://doi.org/10.5281/zenodo.14984543. 

## License
All code is available as under the <a href="https://creativecommons.org/licenses/by-nc/4.0/">Attribution-NonCommercial (CC BY-NC) 4.0 license</a>.