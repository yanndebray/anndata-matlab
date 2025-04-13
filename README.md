
# Read AnnData file format in MATLAB Datastore

anndata is a Python package for handling annotated data - this is a thin MATLAB wrapper on it.

<a name="beginToc"></a>

## Table of Contents
&emsp;[Setup python env](#setup-python-env)
 
&emsp;[Create AnnData file in Python](#create-anndata-file-in-python)
 
&emsp;[Create a Datastore in MATLAB](#create-a-datastore-in-matlab)
 
&emsp;[Cast into a MATLAB datatype](#cast-into-a-matlab-datatype)
 
<a name="endToc"></a>

# Setup python env
```matlab
!python get-pip.py
```

```matlabTextOutput
python: can't open file '/OneDrive/MATLAB/anndata/get-pip.py': [Errno 2] No such file or directory
```

```matlab
!python -m pip install anndata
```

```matlabTextOutput
Defaulting to user installation because normal site-packages is not writeable
Requirement already satisfied: anndata in /home/matlab/.local/lib/python3.10/site-packages (0.11.3)
Requirement already satisfied: array-api-compat!=1.5,>1.4 in /home/matlab/.local/lib/python3.10/site-packages (from anndata) (1.11.2)
Requirement already satisfied: exceptiongroup in /home/matlab/.local/lib/python3.10/site-packages (from anndata) (1.2.2)
Requirement already satisfied: h5py>=3.7 in /home/matlab/.local/lib/python3.10/site-packages (from anndata) (3.13.0)
Requirement already satisfied: natsort in /home/matlab/.local/lib/python3.10/site-packages (from anndata) (8.4.0)
Requirement already satisfied: numpy>=1.23 in /home/matlab/.local/lib/python3.10/site-packages (from anndata) (2.2.4)
Requirement already satisfied: packaging>=20.0 in /home/matlab/.local/lib/python3.10/site-packages (from anndata) (24.2)
Requirement already satisfied: pandas!=2.1.0rc0,!=2.1.2,>=1.4 in /home/matlab/.local/lib/python3.10/site-packages (from anndata) (2.2.3)
Requirement already satisfied: scipy>1.8 in /home/matlab/.local/lib/python3.10/site-packages (from anndata) (1.15.2)
Requirement already satisfied: python-dateutil>=2.8.2 in /usr/lib/python3/dist-packages (from pandas!=2.1.0rc0,!=2.1.2,>=1.4->anndata) (2.8.2)
Requirement already satisfied: pytz>=2020.1 in /home/matlab/.local/lib/python3.10/site-packages (from pandas!=2.1.0rc0,!=2.1.2,>=1.4->anndata) (2025.1)
Requirement already satisfied: tzdata>=2022.7 in /home/matlab/.local/lib/python3.10/site-packages (from pandas!=2.1.0rc0,!=2.1.2,>=1.4->anndata) (2025.1)
```

# Create AnnData file in Python
```matlab
pyrunfile("create_anndata.py")
```

```matlabTextOutput
AnnData file 'example_anndata.h5ad' has been created.
```

# Create a Datastore in MATLAB
```matlab
ds = fileDatastore("*.h5ad","ReadFcn",@py.anndata.read_h5ad)
```

```matlabTextOutput
ds = 
  FileDatastore with properties:

                       Files: {
                              '/OneDrive/MATLAB/anndata/example_anndata.h5ad'
                              }
                     Folders: {
                              '/OneDrive/MATLAB/anndata'
                              }
                 UniformRead: 0
                    ReadMode: 'file'
                   BlockSize: Inf
                  PreviewFcn: @py.anndata.read_h5ad
      SupportedOutputFormats: ["txt"    "csv"    "dat"    "asc"    "xlsx"    "xls"    "parquet"    "parq"    "png"    "jpg"    "jpeg"    "tif"    "tiff"    "wav"    "flac"    "ogg"    "opus"    "mp3"    "mp4"    "m4a"]
                     ReadFcn: @py.anndata.read_h5ad
    AlternateFileSystemRoots: {}

```

Process the data in MATLAB

```matlab
data = readall(ds)
```

```matlabTextOutput
data = 1x1 cell array
    {1x1 py.anndata._core.anndata.AnnData}

```

```matlab
data{1}
```

```matlabTextOutput
ans = 
  Python AnnData with properties:

            T: [1x1 py.anndata._core.anndata.AnnData]
            X: [1x1 py.numpy.ndarray]
     filename: [1x1 py.NoneType]
      is_view: 0
     isbacked: 0
       layers: [1x1 py.anndata._core.aligned_mapping.Layers]
        n_obs: [1x1 py.int]
       n_vars: [1x1 py.int]
          obs: [1x1 py.pandas.core.frame.DataFrame]
    obs_names: [1x1 py.pandas.core.indexes.base.Index]
         obsm: [1x1 py.anndata._core.aligned_mapping.AxisArrays]
         obsp: [1x1 py.anndata._core.aligned_mapping.PairwiseArrays]
          raw: [1x1 py.NoneType]
        shape: [1x2 py.tuple]
          uns: [1x1 py.collections.OrderedDict]
          var: [1x1 py.pandas.core.frame.DataFrame]
    var_names: [1x1 py.pandas.core.indexes.base.Index]
         varm: [1x1 py.anndata._core.aligned_mapping.AxisArrays]
         varp: [1x1 py.anndata._core.aligned_mapping.PairwiseArrays]
         file: [1x1 py.anndata._core.file_backing.AnnDataFileManager]

    AnnData object with n_obs x n_vars = 100 x 10
        obs: 'cell_type', 'sample'
        var: 'gene_symbols'

```

# Cast into a MATLAB datatype
```matlab
double(data{1}.X)
```

```matlabTextOutput
ans = 100x10
    0.7588    0.1857    0.1584    0.3342    0.6175    0.2583    0.3866    0.3526    0.4029    0.5127
    0.0003    0.5374    0.6164    0.9802    0.6877    0.1799    0.7184    0.8926    0.9330    0.1712
    0.3698    0.1750    0.3881    0.2430    0.8898    0.3881    0.4184    0.8182    0.7536    0.6670
    0.2544    0.3308    0.5912    0.2956    0.5890    0.6597    0.4327    0.6635    0.7783    0.1073
    0.0528    0.9274    0.4190    0.2508    0.2715    0.3196    0.1004    0.4111    0.7852    0.5738
    0.8798    0.4915    0.7152    0.4133    0.4345    0.5131    0.2989    0.0315    0.1183    0.2182
    0.7209    0.4909    0.3734    0.9942    0.2845    0.1972    0.6129    0.6037    0.7029    0.0847
    0.8532    0.1064    0.6125    0.4100    0.3122    0.3772    0.5875    0.2830    0.5387    0.1962
    0.5707    0.7673    0.7089    0.7311    0.4006    0.0094    0.8223    0.3179    0.8020    0.7084
    0.7709    0.9276    0.1386    0.4550    0.3837    0.8726    0.4717    0.6807    0.7017    0.0160

```
