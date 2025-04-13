
# Read AnnData file format in MATLAB Datastore

[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=yanndebray/anndata-matlab)

anndata is a Python package for handling annotated data \- this is a thin MATLAB wrapper on it.

<a name="beginToc"></a>

## Table of Contents
&emsp;&emsp;[Setup python env](#setup-python-env)
 
&emsp;&emsp;[Create AnnData file in Python](#create-anndata-file-in-python)
 
&emsp;&emsp;[Create a Datastore in MATLAB](#create-a-datastore-in-matlab)
 
&emsp;&emsp;[Cast into a MATLAB datatype](#cast-into-a-matlab-datatype)
 
&emsp;&emsp;[Utils](#utils)
 
<a name="endToc"></a>

## Setup python env

Inspire from [matlab\-with\-python\-book/8\_Resources.md at main · yanndebray/matlab\-with\-python\-book](https://github.com/yanndebray/matlab-with-python-book/blob/main/8_Resources.md)

```matlab
websave("/tmp/get-pip.py","https://bootstrap.pypa.io/get-pip.py");
!python /tmp/get-pip.py
```

```matlabTextOutput
Defaulting to user installation because normal site-packages is not writeable
Collecting pip
  Downloading pip-25.0.1-py3-none-any.whl.metadata (3.7 kB)
Collecting setuptools
  Downloading setuptools-78.1.0-py3-none-any.whl.metadata (6.6 kB)
Collecting wheel
  Downloading wheel-0.45.1-py3-none-any.whl.metadata (2.3 kB)
Downloading pip-25.0.1-py3-none-any.whl (1.8 MB)
□[?25l   □[90m━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━□[0m □[32m0.0/1.8 MB□[0m □[31m?□[0m eta □[36m-:--:--□[0m
□[2K   □[90m━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━□[0m □[32m1.8/1.8 MB□[0m □[31m98.8 MB/s□[0m eta □[36m0:00:00□[0m
□[?25hDownloading setuptools-78.1.0-py3-none-any.whl (1.3 MB)
□[?25l   □[90m━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━□[0m □[32m0.0/1.3 MB□[0m □[31m?□[0m eta □[36m-:--:--□[0m
□[2K   □[90m━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━□[0m □[32m1.3/1.3 MB□[0m □[31m109.0 MB/s□[0m eta □[36m0:00:00□[0m
□[?25hDownloading wheel-0.45.1-py3-none-any.whl (72 kB)
Installing collected packages: wheel, setuptools, pip
□[33m  WARNING: The script wheel is installed in '/home/matlab/.local/bin' which is not on PATH.
  Consider adding this directory to PATH or, if you prefer to suppress this warning, use --no-warn-script-location.□[0m□[33m
□[0m□[33m  WARNING: The scripts pip, pip3 and pip3.10 are installed in '/home/matlab/.local/bin' which is not on PATH.
  Consider adding this directory to PATH or, if you prefer to suppress this warning, use --no-warn-script-location.□[0m□[33m
□[0mSuccessfully installed pip-25.0.1 setuptools-78.1.0 wheel-0.45.1
```

```matlab
!python -m pip install anndata
```

```matlabTextOutput
Defaulting to user installation because normal site-packages is not writeable
Collecting anndata
  Downloading anndata-0.11.4-py3-none-any.whl.metadata (9.3 kB)
Collecting array-api-compat!=1.5,>1.4 (from anndata)
  Downloading array_api_compat-1.11.2-py3-none-any.whl.metadata (1.9 kB)
Collecting exceptiongroup (from anndata)
  Downloading exceptiongroup-1.2.2-py3-none-any.whl.metadata (6.6 kB)
Collecting h5py>=3.7 (from anndata)
  Downloading h5py-3.13.0-cp310-cp310-manylinux_2_17_x86_64.manylinux2014_x86_64.whl.metadata (2.5 kB)
Collecting natsort (from anndata)
  Downloading natsort-8.4.0-py3-none-any.whl.metadata (21 kB)
Collecting numpy>=1.23 (from anndata)
  Downloading numpy-2.2.4-cp310-cp310-manylinux_2_17_x86_64.manylinux2014_x86_64.whl.metadata (62 kB)
Collecting packaging>=24.2 (from anndata)
  Downloading packaging-24.2-py3-none-any.whl.metadata (3.2 kB)
Collecting pandas!=2.1.0rc0,!=2.1.2,>=1.4 (from anndata)
  Downloading pandas-2.2.3-cp310-cp310-manylinux_2_17_x86_64.manylinux2014_x86_64.whl.metadata (89 kB)
Collecting scipy>1.8 (from anndata)
  Downloading scipy-1.15.2-cp310-cp310-manylinux_2_17_x86_64.manylinux2014_x86_64.whl.metadata (61 kB)
Collecting python-dateutil>=2.8.2 (from pandas!=2.1.0rc0,!=2.1.2,>=1.4->anndata)
  Downloading python_dateutil-2.9.0.post0-py2.py3-none-any.whl.metadata (8.4 kB)
Collecting pytz>=2020.1 (from pandas!=2.1.0rc0,!=2.1.2,>=1.4->anndata)
  Downloading pytz-2025.2-py2.py3-none-any.whl.metadata (22 kB)
Collecting tzdata>=2022.7 (from pandas!=2.1.0rc0,!=2.1.2,>=1.4->anndata)
  Downloading tzdata-2025.2-py2.py3-none-any.whl.metadata (1.4 kB)
Requirement already satisfied: six>=1.5 in /usr/lib/python3/dist-packages (from python-dateutil>=2.8.2->pandas!=2.1.0rc0,!=2.1.2,>=1.4->anndata) (1.16.0)
Downloading anndata-0.11.4-py3-none-any.whl (144 kB)
Downloading array_api_compat-1.11.2-py3-none-any.whl (53 kB)
Downloading h5py-3.13.0-cp310-cp310-manylinux_2_17_x86_64.manylinux2014_x86_64.whl (4.5 MB)
□[?25l   □[90m━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━□[0m □[32m0.0/4.5 MB□[0m □[31m?□[0m eta □[36m-:--:--□[0m
□[2K   □[90m━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━□[0m □[32m4.5/4.5 MB□[0m □[31m168.9 MB/s□[0m eta □[36m0:00:00□[0m
□[?25hDownloading numpy-2.2.4-cp310-cp310-manylinux_2_17_x86_64.manylinux2014_x86_64.whl (16.4 MB)
□[?25l   □[90m━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━□[0m □[32m0.0/16.4 MB□[0m □[31m?□[0m eta □[36m-:--:--□[0m
□[2K   □[90m━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━□[0m □[32m16.4/16.4 MB□[0m □[31m268.2 MB/s□[0m eta □[36m0:00:00□[0m
□[?25hDownloading packaging-24.2-py3-none-any.whl (65 kB)
Downloading pandas-2.2.3-cp310-cp310-manylinux_2_17_x86_64.manylinux2014_x86_64.whl (13.1 MB)
□[?25l   □[90m━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━□[0m □[32m0.0/13.1 MB□[0m □[31m?□[0m eta □[36m-:--:--□[0m
□[2K   □[90m━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━□[0m □[32m13.1/13.1 MB□[0m □[31m218.0 MB/s□[0m eta □[36m0:00:00□[0m
□[?25hDownloading scipy-1.15.2-cp310-cp310-manylinux_2_17_x86_64.manylinux2014_x86_64.whl (37.6 MB)
□[?25l   □[90m━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━□[0m □[32m0.0/37.6 MB□[0m □[31m?□[0m eta □[36m-:--:--□[0m
□[2K   □[90m━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━□[0m □[32m37.6/37.6 MB□[0m □[31m228.6 MB/s□[0m eta □[36m0:00:00□[0m
□[?25hDownloading exceptiongroup-1.2.2-py3-none-any.whl (16 kB)
Downloading natsort-8.4.0-py3-none-any.whl (38 kB)
Downloading python_dateutil-2.9.0.post0-py2.py3-none-any.whl (229 kB)
Downloading pytz-2025.2-py2.py3-none-any.whl (509 kB)
Downloading tzdata-2025.2-py2.py3-none-any.whl (347 kB)
Installing collected packages: pytz, tzdata, python-dateutil, packaging, numpy, natsort, exceptiongroup, array-api-compat, scipy, pandas, h5py, anndata
□[33m  WARNING: The scripts f2py and numpy-config are installed in '/home/matlab/.local/bin' which is not on PATH.
  Consider adding this directory to PATH or, if you prefer to suppress this warning, use --no-warn-script-location.□[0m□[33m
□[0m□[33m  WARNING: The script natsort is installed in '/home/matlab/.local/bin' which is not on PATH.
  Consider adding this directory to PATH or, if you prefer to suppress this warning, use --no-warn-script-location.□[0m□[33m
□[0mSuccessfully installed anndata-0.11.4 array-api-compat-1.11.2 exceptiongroup-1.2.2 h5py-3.13.0 natsort-8.4.0 numpy-2.2.4 packaging-24.2 pandas-2.2.3 python-dateutil-2.9.0.post0 pytz-2025.2 scipy-1.15.2 tzdata-2025.2
```

## Create AnnData file in Python
```matlab
pyrunfile("create_anndata.py")
```

```matlabTextOutput
AnnData file 'example_anndata.h5ad' has been created.
```

## Create a Datastore in MATLAB
```matlab
ds = fileDatastore("*.h5ad","ReadFcn",@py.anndata.read_h5ad)
```

```matlabTextOutput
ds = 
  FileDatastore with properties:

                       Files: {
                              '/MATLAB Drive/anndata-matlab/example_anndata.h5ad'
                              }
                     Folders: {
                              '/MATLAB Drive/anndata-matlab'
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

## Cast into a MATLAB datatype
```matlab
double(data{1}.X)
```

```matlabTextOutput
ans = 100x10
    0.3902    0.4530    0.2303    0.1644    0.1387    0.0044    0.5164    0.7349    0.9660    0.8410
    0.4711    0.7978    0.3014    0.0239    0.4256    0.7271    0.8429    0.4177    0.8212    0.5306
    0.7812    0.0382    0.7872    0.6018    0.3788    0.9143    0.8811    0.8398    0.1012    0.0655
    0.2274    0.4132    0.6263    0.3336    0.8878    0.1644    0.9931    0.9012    0.3154    0.1782
    0.3972    0.8497    0.0596    0.3916    0.3861    0.5077    0.0408    0.5115    0.6194    0.8040
    0.7103    0.5851    0.0006    0.2181    0.8362    0.1073    0.2539    0.0303    0.0563    0.5896
    0.4740    0.2416    0.3401    0.1660    0.4830    0.6881    0.7436    0.9198    0.6605    0.2095
    0.7402    0.6688    0.3366    0.9501    0.7017    0.8497    0.3403    0.9303    0.7928    0.3628
    0.3244    0.2612    0.2610    0.7588    0.8610    0.8312    0.8603    0.9579    0.8278    0.0482
    0.4528    0.3654    0.3428    0.1934    0.1149    0.1027    0.9573    0.6976    0.3066    0.8612

```

## Utils
-  Add badge to Open in MATLAB Online 

`[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=yanndebray/anndata-matlab)`

-  Export Livescript to Markdown (don't forget to add back the badge with the line above) 
```matlab
export livescript.mlx README.md;
```
