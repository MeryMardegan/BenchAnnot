# BenchAnnot environment setup

Create the locked environment:

```bash
conda-lock install --name benchannot 1_setup/conda-lock.yml
```

Activate it and install BenchAnnot in editable mode:

```bash
conda activate benchannot
python -m pip install -e .
```

Register the environment as a Jupyter kernel:

```bash
python -m ipykernel install \
    --user \
    --name benchannot \
    --display-name "Python (BenchAnnot)"
```

Test the installation and inspect the available kernels:

```bash
python -c "import benchannot; print(benchannot.__file__)"
jupyter kernelspec list
```

Start Jupyter Lab:

```bash
jupyter lab
```

Use the `Python (BenchAnnot)` kernel when opening the notebooks. Execute the
eukaryotic notebooks in this order:

1. `1_audit_prepare_reference.ipynb`;
2. `2_prepare_output_tools.ipynb`;
3. `3_functional_analysis.ipynb`.
