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