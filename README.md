# CellCDmT
A Cell-cell communication inference computing framwork.

## Overview of CellCDmT
[simple_flowchat.png](https://github.com/wallwei/CellCDmT/blob/5330dbf01420a5110fdb7cd820aca78b7b83d7b6/simple_flowchat.png)

## Classification for LRPs
[classifier.png](https://github.com/wallwei/CellCDmT/blob/b2f3e5c0fe0b539395480feb709cb2ca105101a8/classifier.png)

## CCC inference
[flowchart_CCC.png](https://github.com/wallwei/CellCDmT/blob/b2f3e5c0fe0b539395480feb709cb2ca105101a8/flowchart_CCC.png)

# Environment
Python == 3.8.20

numpy=1.19.5

pandas=1.4.4

scikit-learn=1.2.2

deep-forest=0.1.7

catboost=1.2.3

ngboost=0.3.6

xgboost=2.1.1

# Usage

### 1.We obtained LRP features from [PyFeat](https://github.com/mrzResearchArena/PyFeat)

### 2.Run CellCDmT to obtain LRIs,for example, run case_study_1:
```python
python code/case_study_1.py
```

### 3.To get CCC score matrix
```python
python code/breast_case_study/3point-method.py
```
