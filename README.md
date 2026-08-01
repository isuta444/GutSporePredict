# GutSporePredict

> Genome-scale prediction of bacterial sporulation, germination and complete lifecycle potential.

---

## Overview

GutSporePredict is an open-source comparative genomics platform for predicting the sporulation and germination capabilities of bacteria directly from genome sequences.

Unlike conventional annotation pipelines that simply detect homologous genes, GutSporePredict reconstructs the biological sporulation program by integrating gene evidence, competitive orthology assignment, biological modules, developmental stages, and lifecycle inference into a unified framework.

The project was originally designed for human gut-associated Firmicutes, particularly members of the families Lachnospiraceae, Oscillospiraceae, Clostridiaceae, Peptostreptococcaceae and related lineages. However, the overall framework is applicable to diverse spore-forming bacteria.

GutSporePredict performs the complete analysis automatically from assembled genome FASTA files and produces biologically interpretable lifecycle predictions.

---

## Why GutSporePredict?

Bacterial sporulation is among the most complex developmental programs known in prokaryotes.

Predicting whether a bacterium can truly complete sporulation cannot be achieved reliably by detecting only a few marker genes such as **spo0A**.

Different bacterial lineages possess distinct evolutionary solutions for sporulation and germination.

Many genes belong to homologous protein families.

Some proteins are lineage-specific.

Some functions are fulfilled by mutually exclusive orthologs.

Simple BLAST or HMM hit counting therefore frequently produces biologically inconsistent predictions.

GutSporePredict addresses these limitations by combining:

- curated biological knowledge
- competitive ortholog assignment
- hierarchical module evaluation
- developmental stage reconstruction
- confidence estimation
- complete lifecycle inference

rather than relying solely on gene presence.

---

## Biological Scope

GutSporePredict evaluates two interconnected biological systems.

### Sporulation

The sporulation program includes

- initiation
- asymmetric division
- sigma factor activation
- engulfment
- cortex synthesis
- coat assembly
- maturation

The framework reconstructs the progression of these developmental processes using curated biological knowledge.

---

### Germination

The germination program evaluates

- germinant receptors
- cortex lytic enzymes
- Csp-mediated pathways
- Ger-mediated pathways
- cortex degradation
- emergence from dormancy

Different germination strategies are evaluated independently before integration into the final lifecycle prediction.

---

## Philosophy

GutSporePredict is **not** intended to be another genome annotation pipeline.

Instead, it aims to function as a biological reasoning engine.

Rather than asking

> "Is this gene present?"

GutSporePredict asks

> "Can this organism realistically complete this biological process?"

This distinction forms the central philosophy of the project.

Every prediction is therefore based on multiple layers of biological evidence rather than single-gene detection.

---

## Key Features

- End-to-end genome analysis from FASTA files
- Automatic protein prediction using Prodigal
- HMMER3-based homology detection
- Competitive ortholog assignment
- Three-state gene calling
- Biological module reconstruction
- Developmental stage inference
- Lifecycle prediction
- Confidence scoring
- Reproducible command-line workflow
- Human-readable TSV outputs
- Modular biological knowledge base
- Extensible YAML-based rule system

---

## Pipeline Overview

The complete workflow consists of eight major steps.

```
Genome FASTA
      │
      ▼
Protein prediction
      │
      ▼
HMM profile search
      │
      ▼
Competitive ortholog assignment
      │
      ▼
Gene call generation
      │
      ▼
Module evaluation
      │
      ▼
Developmental stage evaluation
      │
      ▼
Lifecycle prediction
```

Each stage produces independent output files that can be inspected separately.

---

## Scientific Motivation

Sporulation is a critical survival strategy that enables bacteria to tolerate environmental stress through the production of metabolically dormant endospores.

Within the human gut microbiome, sporulation contributes to

- host-to-host transmission
- environmental persistence
- resistance to oxygen exposure
- resilience during antibiotic treatment
- long-term ecological stability

Despite its importance, experimentally determining sporulation capacity remains labor-intensive.

GutSporePredict was developed to provide a reproducible computational framework capable of predicting sporulation and germination potential directly from genome sequences while maintaining biological interpretability.

---

## Design Principles

GutSporePredict follows several guiding principles.

### Biological correctness

Predictions should reflect biological mechanisms rather than simple sequence similarity.

### Reproducibility

All analyses should be deterministic and fully reproducible.

### Transparency

Every prediction can be traced back to individual gene calls, biological modules, developmental stages, and supporting evidence.

### Extensibility

The knowledge base can be expanded without modifying source code.

### Comparative genomics

The software is designed for both individual genomes and large comparative genomic datasets.

---

## Current Status

Current version

```
4.0.0-beta1
```

Current capabilities

- Complete command-line interface
- End-to-end analysis
- Gene prediction
- Sporulation prediction
- Germination prediction
- Lifecycle prediction
- Module inference
- Developmental stage inference
- Automated testing
- Static type checking

Future releases will further extend evolutionary analyses, phylogenetic reconstruction, ancestral state estimation, and interactive reports.


---

# Installation

## System Requirements

GutSporePredict has been developed and tested primarily on Linux and macOS systems.

Windows users are encouraged to use Windows Subsystem for Linux (WSL2).

### Supported operating systems

- macOS (Apple Silicon and Intel)
- Linux (Ubuntu recommended)
- Windows (WSL2)

---

## Software Requirements

The following external programs are required.

| Software | Minimum Version | Purpose |
|-----------|-----------------|----------|
| Python | 3.12 | Runtime |
| Prodigal | 2.6.3 | Gene prediction |
| HMMER | 3.4 | HMM profile search |

---

## Python Dependencies

The core package has intentionally few dependencies.

Current required packages include

- PyYAML

Development dependencies

- pytest
- mypy
- ruff

---

## Installation using Micromamba

Clone the repository.

```bash
git clone https://github.com/YOUR_USERNAME/GutSporePredict.git

cd GutSporePredict
```

Create the environment.

```bash
micromamba env create \
    -f environments/gutsporepredict.yml
```

Activate the environment.

```bash
micromamba activate gutsporepredict
```

Install GutSporePredict.

```bash
pip install -e .
```

Verify the installation.

```bash
gutsporepredict doctor
```

---

## Installation using Conda

```bash
conda env create \
    -f environments/gutsporepredict.yml

conda activate gutsporepredict

pip install -e .
```

---

## Installation using pip

If all external dependencies are already available,

```bash
pip install -e .
```

is sufficient.

---

# Quick Start

Create a directory containing genome FASTA files.

```
genomes/

    genome1.fna

    genome2.fna

    genome3.fna
```

Run GutSporePredict.

```bash
gutsporepredict run \
    --genomes genomes \
    --output results
```

The analysis automatically performs

1. Protein prediction
2. HMM search
3. Gene assignment
4. Module evaluation
5. Stage evaluation
6. Lifecycle prediction

---

# Expected Runtime

Approximate runtimes for bacterial isolate genomes.

| Number of genomes | Runtime |
|-------------------|----------|
| 1 | <1 minute |
| 10 | Several minutes |
| 100 | Tens of minutes |
| 1000 | Depends on available CPU resources |

Actual runtime depends on

- genome size
- number of genes
- CPU threads
- storage performance

---

# Command Line Interface

Display general help.

```bash
gutsporepredict --help
```

Display run help.

```bash
gutsporepredict run --help
```

Display environment information.

```bash
gutsporepredict doctor
```

---

# Example

Example analysis.

```bash
gutsporepredict run \
    --genomes example_genomes \
    --output example_results
```

Example output.

```
results/

    01_proteins/

    02_hmmsearch_sporulation/

    03_hmmsearch_germination/

    04_gene_calls_sporulation/

    05_gene_calls_germination/

    06_gene_calls/

    07_modules/

    08_stages/

    lifecycle_summary.tsv

    logs/
```

---

# Input

GutSporePredict currently accepts assembled bacterial genomes.

Supported extensions

- .fa
- .fna
- .fasta
- .fas

Each FASTA file should contain one assembled genome.

Draft genomes are supported.

Complete genomes are recommended.

---

# Output

The final prediction is stored in

```
results/lifecycle_summary.tsv
```

Intermediate files are preserved to allow full inspection of every analysis step.

---

# Repository Structure

```
GutSporePredict/

├── config/
├── database/
├── docs/
├── environments/
├── examples/
├── knowledge/
├── scripts/
├── src/
├── tests/
├── README.md
├── LICENSE
└── pyproject.toml
```

---

# Testing

Run all automated tests.

```bash
pytest
```

Run static type checking.

```bash
mypy src tests
```

Run linting.

```bash
ruff check src tests
```

Run all quality checks.

```bash
ruff check src tests

mypy src tests

pytest
```

Current validation

- 75 unit tests
- strict mypy
- Ruff linting
- end-to-end CLI testing

---

# Troubleshooting

## Prodigal not found

Ensure Prodigal is installed and available in PATH.

```
prodigal -v
```

---

## HMMER not found

Verify installation.

```
hmmsearch -h
```

---

## No genomes detected

Verify that

- the input directory exists
- files have supported extensions
- FASTA files are not empty

---

## Empty prediction

Check

- Prodigal logs
- HMM search outputs
- lifecycle_summary.tsv

for diagnostic information.


---

# Biological Knowledge Base

One of the major characteristics of GutSporePredict is that biological knowledge is separated from the analysis engine.

Rather than embedding biological rules directly in Python code, GutSporePredict stores genes, modules, developmental stages, and inference rules as human-readable YAML files.

This design allows the biological knowledge base to evolve independently of the software implementation.

Consequently, updates to sporulation biology generally require modifications only to the knowledge base rather than changes to the source code.

---

# Knowledge Architecture

The biological knowledge is organized into four hierarchical layers.

```
Genes
   │
   ▼
Biological Modules
   │
   ▼
Developmental Stages
   │
   ▼
Lifecycle Prediction
```

Each layer summarizes the biological evidence from the previous layer.

---

# Gene Knowledge

Individual genes are described independently.

Each gene definition includes biological metadata such as

- canonical name
- aliases
- biological function
- developmental role
- pathway assignment
- literature references
- HMM profile

Example

```
knowledge/

    genes/

        spo0A.yaml

        sigH.yaml

        sigF.yaml

        sigE.yaml

        sigG.yaml

        sigK.yaml

        spoIIE.yaml

        gerAA.yaml

        gerAB.yaml

        ...
```

Gene definitions are intentionally independent of prediction algorithms.

---

# Biological Modules

Genes are grouped into biologically meaningful functional modules.

Examples include

- Spo0 phosphorelay
- Sigma factor activation
- Engulfment
- Cortex synthesis
- Coat assembly
- SASP production
- Ger receptor complex
- Cortex lytic enzyme system

Each module contains

- required genes
- optional genes
- alternative genes
- logical relationships
- confidence rules

Modules therefore represent biological processes rather than individual proteins.

---

# Developmental Stages

Modules are integrated into developmental stages.

Current stages include

| Stage | Biological Process |
|--------|--------------------|
| ST001 | Sporulation initiation |
| ST002 | Asymmetric septation |
| ST003 | Sigma activation |
| ST004 | Engulfment |
| ST005 | Intercellular signaling |
| ST006 | Cortex formation |
| ST007 | Coat assembly |
| ST008 | Spore maturation |
| ST009 | Germination competence |

Stage predictions summarize whether an organism is capable of completing each developmental transition.

---

# Lifecycle Inference

The final prediction is generated after integrating

- gene evidence
- module completeness
- developmental stages
- assessment scores
- confidence estimates

This hierarchical approach reduces false positive predictions caused by isolated homologous genes.

---

# Competitive Ortholog Assignment

Many sporulation proteins belong to homologous protein families.

A single HMM hit therefore does not necessarily identify the biologically correct ortholog.

GutSporePredict resolves these situations using competitive ortholog assignment.

Within each ortholog group,

1. all candidate HMM hits are collected,
2. scores are compared,
3. the biologically most likely ortholog is selected,
4. remaining candidates are rejected.

This procedure substantially reduces ambiguous annotations.

---

# Three-State Gene Calls

Each gene is assigned one of three states.

| State | Meaning |
|--------|---------|
| Present | High-confidence ortholog detected |
| Absent | No evidence detected |
| Ambiguous | Evidence exists but is insufficient for confident assignment |

Three-state classification preserves uncertainty instead of forcing binary decisions.

---

# Confidence Estimation

GutSporePredict reports confidence separately from prediction scores.

Typical confidence levels include

- High
- Moderate
- Low

Confidence depends on

- fraction of assessed genes
- module completeness
- competing evidence
- biological consistency

Confidence therefore reflects prediction reliability rather than biological importance.

---

# Knowledge Files

Current knowledge directories

```
knowledge/

    genes/

    modules/

    rules/

    stages.yaml
```

Configuration files

```
config/

    hmm/

    gtdb_targets/

    benchmarks/
```

These files can be modified without changing the software source code.

---

# Extending the Knowledge Base

New biological pathways can be incorporated by

1. adding new gene definitions,
2. defining biological modules,
3. updating developmental stages,
4. creating inference rules.

No modifications to the prediction engine are required.

This architecture makes GutSporePredict readily extensible to future discoveries in sporulation and germination biology.

---

# Design Philosophy

The prediction engine intentionally contains very little biological knowledge.

Instead,

- biological knowledge resides in YAML files,
- computational logic resides in Python,
- external tools perform sequence analyses.

This separation improves

- maintainability,
- reproducibility,
- transparency,
- extensibility.

It also enables independent peer review of both biological assumptions and software implementation.

---

# Current Biological Coverage

Current knowledge base includes

- Sporulation initiation
- Sigma factor cascade
- Engulfment machinery
- Cortex synthesis
- Coat formation
- Small acid-soluble proteins
- Germinant receptors
- Cortex lytic enzymes
- Csp-mediated germination
- Ger-mediated germination

Additional pathways will be incorporated in future releases as new biological evidence becomes available.


---

# Prediction Pipeline

GutSporePredict performs a complete end-to-end analysis beginning with assembled bacterial genomes and ending with biologically interpretable lifecycle predictions.

The pipeline consists of eight sequential stages.

```
Genome FASTA
      │
      ▼
Protein Prediction
      │
      ▼
HMM Search
      │
      ▼
Competitive Ortholog Assignment
      │
      ▼
Three-State Gene Calls
      │
      ▼
Biological Module Evaluation
      │
      ▼
Developmental Stage Evaluation
      │
      ▼
Lifecycle Prediction
```

Each stage produces independent output files, allowing every prediction to be inspected and reproduced.

---

# Step 1 — Protein Prediction

Input genomes are processed using Prodigal.

For each genome the following files are generated.

| File | Description |
|------|-------------|
| *.faa | Predicted protein sequences |
| *.ffn | Coding nucleotide sequences |
| *.gff | Gene coordinates |

Gene prediction is intentionally performed before all downstream analyses to ensure a consistent annotation strategy across all genomes.

---

# Step 2 — HMM Profile Search

Each predicted protein dataset is searched against curated HMM profiles using HMMER3.

Every reference gene has an associated profile HMM.

For each genome, GutSporePredict records

- bit score
- E-value
- domain E-value
- alignment coverage
- matched HMM
- matched protein

Only statistically significant hits proceed to the next stage.

---

# Step 3 — Competitive Ortholog Assignment

Many sporulation genes belong to homologous protein families.

A protein may therefore produce significant matches against multiple HMM profiles.

Instead of accepting all matches, GutSporePredict performs competitive assignment.

Within each competition group

1. all candidate hits are collected,
2. scores are ranked,
3. biological constraints are evaluated,
4. the most probable ortholog is retained,
5. competing assignments are rejected.

This substantially reduces false positive annotations.

---

# Step 4 — Three-State Gene Calling

Following competitive assignment, each reference gene receives one of three possible states.

| State | Interpretation |
|--------|----------------|
| Present | High-confidence ortholog detected |
| Ambiguous | Evidence exists but remains inconclusive |
| Absent | No convincing evidence detected |

Unlike binary annotation systems, ambiguous evidence is preserved throughout the analysis.

This prevents premature biological conclusions.

---

# Step 5 — Biological Module Evaluation

Individual genes are integrated into biological modules.

Examples include

- Spo0 phosphorelay
- Sigma factor cascade
- Engulfment machinery
- Cortex synthesis
- Coat formation
- SASP production
- Ger receptor complex
- Cortex lytic enzyme system

Each module is evaluated independently.

Outputs include

- completeness
- score
- confidence
- assessment fraction

Missing optional genes do not necessarily invalidate a module.

Conversely, absence of essential genes prevents module completion.

---

# Step 6 — Developmental Stage Evaluation

Modules are integrated into developmental stages.

Each stage represents a biologically meaningful transition during sporulation.

Current stages include

| Stage | Process |
|--------|---------|
| ST001 | Initiation |
| ST002 | Septation |
| ST003 | Sigma activation |
| ST004 | Engulfment |
| ST005 | Intercellular signaling |
| ST006 | Cortex formation |
| ST007 | Coat assembly |
| ST008 | Spore maturation |
| ST009 | Germination competence |

Stages summarize higher-order biological functions rather than individual genes.

---

# Step 7 — Lifecycle Prediction

The final prediction combines

- gene evidence
- module evidence
- developmental stages
- assessment fractions
- confidence estimates

Both sporulation and germination are evaluated independently before integration.

Possible lifecycle predictions include

| Prediction | Interpretation |
|------------|----------------|
| complete_lifecycle | Sporulation and germination appear complete |
| sporulation_only | Sporulation complete but germination incomplete |
| germination_only | Germination genes present without complete sporulation |
| incomplete_lifecycle | Partial evidence for both systems |
| non_spore_former | No convincing sporulation program detected |

---

# Scoring System

Each prediction produces a numerical score between 0 and 1.

General interpretation

| Score | Interpretation |
|--------|----------------|
| 0.90–1.00 | Nearly complete pathway |
| 0.70–0.90 | Strong biological evidence |
| 0.40–0.70 | Partial pathway |
| 0.10–0.40 | Weak evidence |
| 0.00–0.10 | No meaningful evidence |

Scores alone should not be interpreted without accompanying confidence estimates.

---

# Assessment Fraction

Not every biological pathway can always be evaluated completely.

GutSporePredict therefore reports an assessment fraction.

Assessment reflects the proportion of pathway components that could be evaluated reliably.

Example

```
sporulation_score       0.95

sporulation_assessment  0.82

confidence              high
```

High scores with poor assessment should be interpreted cautiously.

---

# Confidence Estimation

Confidence reflects the reliability of a prediction rather than biological significance.

Current confidence levels include

- High
- Moderate
- Low

Confidence depends upon

- available biological evidence
- completeness of assessment
- internal consistency
- competing ortholog assignments

---

# Intermediate Results

Unlike many prediction tools, GutSporePredict intentionally preserves every intermediate result.

Researchers can inspect

- HMM hits
- competitive assignments
- gene calls
- module scores
- stage scores
- lifecycle predictions

This design facilitates debugging, benchmarking, and biological interpretation.

---

# Reproducibility

Given identical

- genome assemblies,
- software versions,
- HMM profiles,
- knowledge base,

GutSporePredict produces deterministic and fully reproducible results.

No stochastic algorithms are used during prediction.

---

# Parallelization

HMM searches may be parallelized using

```
--threads
```

Increasing thread count accelerates sequence searches while leaving prediction logic unchanged.

Therefore, results remain identical regardless of CPU count.

---

# Error Handling

The pipeline validates

- genome inputs,
- external dependencies,
- required resources,
- intermediate outputs,
- final prediction files.

Whenever possible, informative error messages identify the failing analysis stage.

---

# Pipeline Outputs

The final prediction is only one component of the generated results.

The complete output directory documents every decision made during the analysis.

This enables

- reproducibility,
- manual inspection,
- downstream comparative analyses,
- method development,
- benchmarking against experimental data.


# Developer Guide

GutSporePredict was designed as a modular prediction framework.

Rather than embedding biological knowledge directly into Python source code, nearly all biological logic is stored as editable resource files.

This architecture allows biological knowledge to evolve independently from the software implementation.

---

# Project Architecture

```
gutsporepredict/

├── annotation/
├── assignment/
├── evidence/
├── gene_prediction/
├── io/
├── knowledge/
├── models/
├── pipeline/
├── prediction/
├── qc/
├── reference/
├── search/
├── validation/
└── visualization/
```

Each package has a single responsibility.

---

# Biological Knowledge Hierarchy

The prediction engine operates using four biological layers.

```
Genes
      │
      ▼
Modules
      │
      ▼
Developmental Stages
      │
      ▼
Lifecycle Prediction
```

This hierarchy mirrors biological organization.

Individual genes are rarely interpreted independently.

Instead, genes form functional modules, modules form developmental stages, and stages determine overall phenotype.

---

# Knowledge Files

All biological knowledge is stored as YAML.

Examples include

```
knowledge/

genes/
modules/
rules/
stages.yaml
```

No Python modification is required when biological definitions change.

---

# Gene Definitions

Each gene has an individual YAML file.

Example

```
knowledge/genes/spo0A.yaml
```

Typical fields include

- gene identifier
- aliases
- biological description
- essentiality
- competition group
- literature references

---

# Module Definitions

Modules define biological functions.

Example

```
knowledge/modules/SP004.yaml
```

Typical contents include

```
required_genes

optional_genes

minimum_score

assessment_rules
```

Modules are intentionally independent.

A new biological pathway can be added without modifying the prediction engine.

---

# Rule Engine

Rules define logical relationships.

Example

```
knowledge/rules/R0001.yaml
```

Rules may specify

```
IF

required modules

minimum score

minimum assessment

THEN

stage becomes COMPLETE
```

Future releases may support more sophisticated logical operators.

---

# Stage Definitions

Developmental stages are described in

```
knowledge/stages.yaml
```

Each stage contains

- identifier
- description
- required modules
- optional modules
- scoring strategy

The pipeline automatically evaluates stages using these definitions.

---

# HMM Profiles

Reference HMM profiles are stored separately from source code.

Each profile corresponds to one biological target.

```
spo0A.hmm

sigE.hmm

gerAA.hmm

...
```

New HMMs can be incorporated without changing Python code.

---

# Competition Groups

Some proteins belong to highly similar protein families.

Competition groups prevent multiple annotations for the same protein.

Current definitions are stored in

```
competition_groups.tsv

germination_competition_groups.tsv
```

These files determine which HMMs compete with one another.

---

# Resource Independence

The prediction engine never hardcodes

- gene names
- stage definitions
- module definitions
- HMM locations
- competition groups

Instead, these are discovered dynamically from package resources.

This separation greatly simplifies maintenance.

---

# Adding a New Gene

Adding support for a new gene typically requires four steps.

1.

Create

```
knowledge/genes/new_gene.yaml
```

2.

Add a reference HMM

```
new_gene.hmm
```

3.

Register the gene within the target configuration.

4.

Update any affected module definitions.

No Python code should require modification.

---

# Adding a New Module

Create a YAML file

```
knowledge/modules/SP010.yaml
```

Define

- required genes
- optional genes
- scoring rules

The module evaluator automatically incorporates the new module.

---

# Adding a New Developmental Stage

Edit

```
knowledge/stages.yaml
```

Define

- stage identifier
- participating modules
- assessment strategy

No pipeline modification is required.

---

# Testing Philosophy

Every major software component includes unit tests.

Current coverage includes

- knowledge loading
- reference parsing
- HMM parsing
- assignment logic
- evidence integration
- rule evaluation
- CLI behaviour

Regression tests ensure biological predictions remain stable across software revisions.

---

# Code Quality

Development follows strict quality control.

Before every commit

```
ruff check

mypy

pytest
```

should all complete successfully.

---

# Design Principles

GutSporePredict follows several core principles.

- deterministic predictions

- biologically interpretable outputs

- modular architecture

- reproducible analyses

- editable biological knowledge

- transparent intermediate results

- minimal hidden heuristics

These principles guide future development.

---

# Versioning

GutSporePredict follows semantic versioning.

```
MAJOR.MINOR.PATCH
```

Examples

```
4.0.0

4.1.0

4.1.2
```

Major releases introduce new prediction capabilities.

Minor releases extend biological knowledge.

Patch releases fix software defects.

---

# Citation

If GutSporePredict contributes to published work, please cite the associated publication once available.

Until then, please cite the GitHub repository.

```
Hisatomi A.

GutSporePredict:
A comparative genomics platform for predicting bacterial sporulation and germination potential.

GitHub repository.
```


# Roadmap

GutSporePredict is under active development.

The long-term goal is to provide a comprehensive evolutionary framework for bacterial sporulation and germination.

---

## Version 4.x

Current development focuses on improving prediction accuracy and usability.

Planned features include

- additional sporulation genes
- expanded germination pathways
- improved competition-group handling
- improved confidence scoring
- enhanced HTML reports
- automatic benchmark generation

---

## Version 5.x

Future releases will incorporate comparative genomics.

Planned features include

- GTDB-wide genome screening
- phylogenetic reconstruction
- ancestral state reconstruction
- gain/loss analysis
- evolutionary module conservation
- family-specific prediction models

---

## Version 6.x

Long-term development aims to integrate ecological information.

Potential features include

- habitat-aware prediction
- microbiome-wide lifecycle profiling
- metagenome support
- machine-learning assisted confidence estimation
- phenotype prediction from incomplete MAGs

---

# Limitations

GutSporePredict predicts biological potential from genome sequences.

It does **not** demonstrate biological activity.

Presence of sporulation or germination genes should not be interpreted as proof that an organism sporulates under laboratory conditions.

Experimental validation remains essential.

Current predictions may also be affected by

- incomplete genomes
- fragmented MAGs
- sequencing errors
- annotation quality
- currently available biological knowledge

Future releases will continue improving robustness as additional reference genomes and experimentally validated datasets become available.

---

# Frequently Asked Questions

## Which genomes can be analysed?

Any bacterial genome in FASTA format.

Supported extensions include

```
.fa
.fna
.fasta
.fas
```

---

## Are MAGs supported?

Yes.

However, fragmented assemblies may reduce prediction confidence because genes may be absent due to incomplete assembly rather than true biological absence.

---

## Does GutSporePredict require annotated genomes?

No.

Protein prediction is performed automatically using Prodigal.

---

## Can I replace Prodigal with Prokka?

Not yet.

Support for alternative annotation pipelines is planned for future releases.

---

## Can I add my own HMM profiles?

Yes.

New HMMs can be added together with corresponding knowledge definitions.

The prediction engine automatically incorporates newly registered targets.

---

## Is Linux required?

No.

GutSporePredict currently supports

- macOS
- Linux

Windows users are encouraged to use WSL2.

---

## Does the software require internet access?

No.

Once installed, all analyses are performed locally.

---

# Contributing

Contributions are welcome.

Examples include

- bug reports
- feature requests
- documentation improvements
- biological knowledge updates
- additional HMM profiles
- code improvements

Before submitting a pull request, please ensure that

```
ruff check

mypy

pytest
```

all complete successfully.

---

# License

GutSporePredict is released under the MIT License.

See

```
LICENSE
```

for details.

---

# Acknowledgements

Development of GutSporePredict has benefited from many discussions with collaborators working in

- comparative genomics
- bacterial sporulation
- bacterial germination
- gut microbiology
- microbial cultivation

The project also builds upon numerous community-developed bioinformatics tools, including

- Prodigal
- HMMER
- GTDB
- NCBI Datasets
- Python
- PyYAML

The authors gratefully acknowledge these open-source projects.

---

# Contact

**Atsushi Hisatomi**

RIKEN BioResource Research Center

Microbe Division

GitHub:
https://github.com/<YOUR_ACCOUNT>/GutSporePredict

Questions, bug reports, and feature requests are welcome through GitHub Issues.

---

# Citation

If GutSporePredict contributes to published work, please cite the associated publication.

Until a manuscript becomes available, please cite the GitHub repository.

```
Hisatomi A.

GutSporePredict:
Comparative genomics framework for predicting bacterial sporulation and germination potential.

GitHub repository.
```

---

**GutSporePredict**

*Genome → Genes → Modules → Developmental Stages → Lifecycle Prediction*