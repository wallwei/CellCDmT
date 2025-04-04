# CellCDmT
A cell-cell communication inference computing framwork.

## Overview of CellCDmT

<img src="https://github.com/wallwei/CellCDmT/blob/9fed3e871cfdff2ddd9a7656da6dfc9c7f634d9f/simple_flowchat.png" width = 50%>

## Classification for LRPs

<img src="https://github.com/wallwei/CellCDmT/blob/9fed3e871cfdff2ddd9a7656da6dfc9c7f634d9f/classifier.png" width = 50%>

## CCC inference

<img src="https://github.com/wallwei/CellCDmT/blob/a5e6f1d3e6af7dcede6cdc52e968ca5ec1df562d/flowchat_CCC.jpg" width = 50%>


# Environment
Python == 3.8.20

numpy == 1.19.5

pandas == 1.4.4

scikit-learn == 1.2.2

deep-forest == 0.1.7

catboost == 1.2.3

xgboost == 2.1.1

# Data
1. scRNA-seq can obtain from [GEO datasets](https://www.ncbi.nlm.nih.gov/gds)
2. Feature extraction tools([PyFeat](https://github.com/mrzResearchArena/PyFeat))

# Usage

## 1.We obtained LRP features by [PyFeat](https://github.com/mrzResearchArena/PyFeat)

## 2.Run CellCDmT to obtain LRIs,for example, run case_study_1:
```python
python code/case_study_1.py
```

## 3.Combining single-cell expression matrix to filter LRIs:
```python
python code/breast_case_study/LRI_filter.py
```

## 4.To get CCC score matrix and visualize CCC network
```python
python code/breast_case_study/3point-method.py
```
