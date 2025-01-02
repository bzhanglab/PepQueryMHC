# PepQueryHLA

---
- [About](#about)
- [Usage](#usage)
  - [Quick start](#quick-start)
  - [Parameters](#parameters)
  - [Scan mode](#scan-mode)
  - [Target mode](#target-mode)
  - [Annotate mode](#annotate-mode)
---

## About

The accurate prioritization of tumor antigens, including aberrant translational products, is critical for the development of personalized cancer immunotherapies. PepQueryHLA estimates a comprehensive repertoire of local RNA expression of tumor antigens within minutes per sample.
<br>

## Usage
PepQueryHLA provides three main functions such as 1) scan mode, 2) target mode and 3) annotate mode.

### Quick start
Scan mode
```bash
java -Xmx2G -jar PepQueryHLA.jar \
--mode scan \
--input peptides.tsv \
--bam sample.sorted.bam \
--output sample \
--thread 16
```
Target mode
```bash
java -Xmx1G -jar PepQueryHLA.jar \
--mode target \
--input peptides_locations_strands.tsv \
--bam sample.sorted.bam \
--output sample \
--thread 16
```
Annotate mode
```bash
java -Xmx2G -jar PepQueryHLA.jar \
--mode annotate \
--input locations_strands.tsv \
--gtf reference_annotation.gtf \
--output sample
```

### Parameters
|Option    | Description    | Value   | Scan mode   | Target mode   | Annotate mode   |
| :---:    | :---:          | :---:   | :---:       | :---:         | :---:           |
| m/mode   | mode to use| scan\|target\|annotate\|  | YY          | YY            | YY              |
| i/input  | input file path| string  | YY          | YY            | YY              |
