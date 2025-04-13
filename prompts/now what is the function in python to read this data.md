To read an AnnData file in Python, you can use the `read_h5ad()` function from the `anndata` package. Here’s a quick example:

```python
import anndata

# Read the AnnData file
adata = anndata.read_h5ad("example_anndata.h5ad")

# Verify by printing a summary of the object
print(adata)
```

This function loads the `.h5ad` file (which is a common format for single-cell data) and returns an AnnData object that you can work with in your analysis.