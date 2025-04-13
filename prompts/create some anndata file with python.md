Below is an example of how you can create a simple AnnData object and save it to an H5AD file using Python. First, make sure you have the necessary packages installed:

```bash
pip install anndata numpy pandas
```

Then, you can use the following Python script:

```python
import numpy as np
import pandas as pd
import anndata

# Create some random data: 100 cells and 10 genes
X = np.random.rand(100, 10)

# Create dummy observation metadata (for cells)
obs = pd.DataFrame({
    'cell_type': ['B' if i % 2 == 0 else 'T' for i in range(100)],
    'sample': ['sample1'] * 100
}, index=[f'cell_{i}' for i in range(100)])

# Create dummy variable metadata (for genes)
var = pd.DataFrame({
    'gene_symbols': [f'gene_{i}' for i in range(10)]
}, index=[f'gene_{i}' for i in range(10)])

# Create the AnnData object
adata = anndata.AnnData(X=X, obs=obs, var=var)

# Save the AnnData object to an .h5ad file
adata.write_h5ad("example_anndata.h5ad")

print("AnnData file 'example_anndata.h5ad' has been created.")
```

### Explanation

- **Data Matrix (X):** A 100×10 matrix of random values representing expression levels.
- **Observation Metadata (obs):** A DataFrame with cell-related metadata, including a simple cell type classification and sample name. The index labels each cell.
- **Variable Metadata (var):** A DataFrame with gene-related metadata, here just providing gene symbols, with each gene indexed.
- **AnnData Object:** Combines the data matrix, observations, and variable metadata.
- **File Saving:** The AnnData object is saved in H5AD format, which is commonly used for single-cell data.

This script should help you get started with creating and saving AnnData files using Python.