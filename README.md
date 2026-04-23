# **A Telomere-To-Telomere Genome Assembly Of Native Desert-Adapted Sheep *Ovis Aries* Breed “Najdi” From Saudi Arabia**

We assembled telomere-to-telomere (T2T) genome assemblies for Najdi sheep using PacBio HiFi, ONT, and Hi-C data.

This repository contains scripts and workflows used for:

- Genome assembly
- Genome polishing
- Assembly assessment
- Chromosome reconstruction
- Gap closure
- Telomere rescue

---

## Repository Structure

```text
.
├── 01_Assembly/
├── 02_Polish/
├── 03_Assessment/
├── 04_Chromosome_Reconstruction/
├── 05_Gap_Closure/
├── 06_Telomere_Rescue/
├── LICENSE
└── README.md
```

---

## Description

## 01_Assembly

Contains scripts for:

- PacBio HiFi de novo assembly using Hifiasm  
- Hi-C scaffolding  
- SALSA (Male assembly)  
- YaHS (Female assembly)  
- Initial chromosome-scale draft assemblies

---

## 02_Polish

Contains polishing workflows including:

- Racon  
- TGS-GapCloser  
- Funannotate clean

### Output

- Polished genome assemblies

---

## 03_Assessment

Contains quality assessment workflows for:

- BUSCO  
- QUAST  
- Merqury  
- NUCmer synteny analysis

### Metrics

- Completeness  
- Contiguity  
- QV  
- Structural consistency

---

## 04_Chromosome_Reconstruction

Contains workflows for:

- Synteny-guided chromosome assignment  
- Hi-C validation  
- Reference-guided chromosome reconstruction using RagTag

---

## 05_Gap_Closure

Contains scripts for:

- Gap detection  
- Minimap2 alignment  
- MAFFT validation  
- MEGA12 inspection  
- Custom genome completion scripts

---

## 06_Telomere_Rescue

Contains workflows for:

- Telomere detection using TIDK  
- Anchor-based telomere recovery  
- Chromosome end extension  
- Final telomere rescue

---

## Dependencies

- Python 3  
- Hifiasm  
- SALSA  
- YaHS  
- Racon  
- TGS-GapCloser  
- RagTag  
- Minimap2  
- BUSCO  
- QUAST  
- Merqury  
- TIDK

---

## Usage

Configure:

```bash
edit config.yaml
```

Run:

```bash
python driver.py --config config.yaml
```

