# Custom Cell Segmentation with Cellpose + Vizgen Post-Processing Tool

This repository contains the training data, the final Cellpose model, and command line scripts for using the **[Vizgen Post-Processing Tool (VPT)](https://github.com/Vizgen/vizgen-postprocessing)**

---

## Overview
We developed a **custom Cellpose 2.0 model** using a human-in-the-loop approach to re-segment MERSCOPE data.  
This repo shares the training data, annotations, trained model, and workflows used for the results in the paper.

---

## Methods Summary
- **Segmentation input:** Maximum intensity projection of three cell boundary stains at z = 3 (`Cellbound_MAX`) + DAPI  
- **Annotation:** 53 patches extracted with Vizgen Post-Processing Tool (VPT), of which 48 were manually annotated and 5 reserved as a test set  
- **Preprocessing:** CLAHE normalisation (clip limit 0.01) applied to all images  
- **Training:** Pretrained *cyto2* model fine-tuned uwith:  
  - Epochs: 100  
  - Learning rate: 0.1  
  - Weight decay: 1×10⁻⁴  

- **VPT** with the cellpose2 plugin (vpt-plugin-cellpose2):  
  - Cell diameter = 93.11 px  
  - Flow threshold = 0.6  
  - Cell probability threshold = -4  
  - Minimum mask size = 100 px
  - The VPT workflow was implemented using an **[automated NextFlow pipeline](https://github.com/WEHI-SODA-Hub/spatialvpt)**
     
- **Downstream use:** VPT resegmented the data and regenerated the cell-by-gene matrix for downstream analysis.

---

## Repository Contents
- `/training_images` – image patches and annotations
- `/cellpose_model` – final trained Cellpose model weights  
- `/scripts` – example VPT commands and cellpose2_custom.json file used for the NextFlow pipeline 
