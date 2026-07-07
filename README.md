# Full Track Evaluator — Enformer/Borzoi Test Set Overlap

This Evaluator benchmarks track predictions against measured DNase-seq signal. It requests chromatin accessibility track predictions for K562 cell type (*homo_sapiens*) and calculates Pearson correlation and Jensen-Shannon divergence as metrics of Predictor accuracy.

## Important Links

- Main GAME Repository: [de-Boer-Lab/Genomic-API-for-Model-Evaluation](https://github.com/de-Boer-Lab/Genomic-API-for-Model-Evaluation)
- GAME Documentation: [ReadTheDocs](https://genomic-api-for-model-evaluation-documentation.readthedocs.io)
- Pre-built Evaluator container image and data: [Hugging Face](https://huggingface.co/datasets/deBoerLab/Track_Evaluator_Borzoi_Enformer_TestSet)
- List of all [GAME Modules](https://github.com/de-Boer-Lab/GAME_modules)

## Overview

This evaluator:
- Sends genomic sequences drawn from the overlapping test regions of Enformer and Borzoi to Predictors and requests accessibility track predictions
- Compares predictions to a measured DNase-seq signal track (K562) (.bigWig)
- Calculates Pearson r and Jensen-Shannon divergence per prediction task
- Supports both JSON and MessagePack data formats


To get started with running the pre-built container, download the container and data from Hugging Face (see [Important Links](#important-links) above).

### Running with Container

```bash
apptainer run --containall \
  --B /path/to/evaluator_data:/evaluator_data \
  --B /path/to/output:/predictions \
  track_evaluator.sif <predictor_ip> <predictor_port> /predictions
```

**Example:**
```bash
apptainer run --containall \
  --B ./evaluator_data:/evaluator_data \
  --B ./results:/predictions \
  full_track_evaluator.sif 192.168.1.100 8000 /predictions
```

If you would like to re-build the container, edit the corresponding paths in the `.def` file and build using the command below. Please make sure to follow the directory structure shown below.

### Container Build

```bash
apptainer build full_track_evaluator.sif evaluator.def
```

## Data Preparation

The input parquet file (`enformer_borzoi_test_seqs.parquet`) must be generated before running the evaluator. See `evaluator_data/README` for details on required input files and how to run `sequence_design.py` to produce this file. The measured data was `ENCFF972GVB.bigWig` was downloaded from here: https://www.encodeproject.org/experiments/ENCSR000EKS/

## Output Files

The evaluator generates two output files:

### 1. Raw Predictions JSON

`<EVALUATOR_NAME>_predictions_<input_file>_from_<PREDICTOR_NAME>.json`

Contains:
- Predictor name and metadata
- All prediction tasks and their parameters
- Predictions for each sequence

### 2. Evaluation Metrics CSV

`evaluation_summary_<EVALUATOR_NAME>.csv`

Contains:
- Evaluator name and description
- Predictor name and timestamp
- Pearson r and Jensen-Shannon divergence per task
- Prediction task metadata

**Example output:**
```tsv
Evaluator	Description	Predictor_name	Time_stamp	Metric	Value	Prediction_task(s)_data
Enformer_Borzoi_TestSetOverlap	DNase Track Request (K562)	MyPredictor	20250125-143022.123456	pearson_r	0.872	[{...}]
Enformer_Borzoi_TestSetOverlap	DNase Track Request (K562)	MyPredictor	20250125-143022.123456	Jensen-Shannon divergence	0.134	[{...}]
```

## API Specification

The Evaluator communicates with Eredictor services via REST API:

### Format Negotiation

```
GET /formats
```

Returns supported request and response formats.

### Prediction Request

```
POST /predict
REQUEST_FORMAT = "application/json"
RESPONSE_FORMAT = "application/msgpack"
```

**Request Body:**
```json
{
  "readout": "track",
  "prediction_tasks": [
    {
      "name": "TrackEvaluation",
      "type": "accessibility",
      "cell_type": "K562",
      "scale": "linear",
      "species": "homo_sapiens"
    }
  ],
  "sequences": {
    "chr1:100000-296608": "ATCGATCG...",
    "chr2:500000-696608": "CGATCGAT..."
  }
}
```

## Evaluation Metrics

**Pearson Correlation Coefficient (r)**
Measures the correlation between predicted and measured signal at base-pair resolution. Binned predictions are expanded to bp-level using the Predictor's reported `bin_size` and  `trim_upstream`.

**Jensen-Shannon Divergence**
Measures the distributional similarity between predicted and measured signal tracks.

## Directory Structure

```
Full_track_evaluator/
├── config.py                        # Configuration settings
├── evaluator_RestAPI.py             # Main evaluator script
├── data_loader.py                   # Input data loading and validation
├── evaluator_content_handler.py     # API communication logic
├── evaluator_metrics_calculator.py  # Metrics computation
├── evaluator.def                    # Apptainer container definition
├── sequence_design.py               # One-time input data preparation
├── evaluator_data/                  # Input data directory (see evaluator_data/README)
└── README.md                        # This file
```


## Development
This codebase serve as a templates for creating other custom track evaluators. The modular design makes it easy to adapt for different evaluation tasks. 
