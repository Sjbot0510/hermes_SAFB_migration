# hermes_SAFB_migration

Migration project managed via Hermes (OpenShell sandbox environment).

## Environment

- **Build toolchain**: `x86_64-conda-linux-gnu-gcc` / `x86_64-conda-linux-gnu-g++`
- **Python**: `/sandbox/miniforge3/envs/build/bin/python`
- **Compiler vars**: `$CC`, `$CXX` (use these; do not call compilers directly)
- **Conda**: solver = `classic`, prefer `conda-forge`, `CONDA_NO_PLUGINS=true`, `CONDA_OVERRIDE_CUDA=""`

## Quick Start

```bash
# Initialize environment before any build/compile task
source /sandbox/setup_build.sh

# Git setup
cd /sandbox/hermes_SAFB_migration
git init
git remote add origin <your-remote-url>
git branch -M main
git add .
git commit -m "initial commit"
git push -u origin main
```

## Pre-installed Packages

### Build / Dev Tools
`gcc_linux-64`, `gxx_linux-64`, `cmake`, `make`, `gh`, `vim`, `ruff`, `black`, `mypy`, `pytest`

### Scientific / Data
`numpy`, `scipy`, `sympy`, `pandas`, `matplotlib`

### Processing / Parsing
`pymupdf`, `pylatexenc`, `latexcodec`, `python-docx`, `openpyxl`, `markdown`, `beautifulsoup4`, `lxml`

Do not reinstall these packages.

## Contributing

1. Create a branch: `git checkout -b <branch-name>`
2. Commit changes: `git commit -m "description"`
3. Push and open a PR via `gh`.
