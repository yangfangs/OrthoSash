# OrthoSash
Encoding sequence by simHash for infer Ortholog based machine learning method.
# Step 1: Encoding sequence

## Parameters

| Parameter |  Description                                                                                |
|:------- |:-------------------------------------------------------------------------------------------|
|  kmer_k    | Setting kmer length.                                            |
|  min_size     |  Hash size.                                                          |
|  genome     |  The path of protein.                                                             |

# Step 2: ML training
## Parameters
| Parameter |  Description                                                                                |
|:------- |:-------------------------------------------------------------------------------------------|
|  pos_path   | Positive pairs.                                            |
|  train_dir     |  The file to store training data.                                                          |
|  res_dir     |  Result file.                                                             |
|  files     |  The path of 78 species proteome.                                                         |


# Step 3: Prediction

## Parameters

| Parameter |  Description                                                                                |
|:------- |:-------------------------------------------------------------------------------------------|
|  all_trains   | The path of training data.                                            |
|  feature_path     |  The encoding feature path.                                                          |
|  res_dir     |  Result file.                                                             |
|  files     |  The path of 78 species proteome.                                                         |

# Python package (need install)
* numpy
* pandas
* sklearn
* mmh3
* Bio