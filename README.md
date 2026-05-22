# hermes_SAFB_migration

Hermes Sandbox for SAFB (Single Agent Framework Builder) Migration.

## Environment

This sandbox runs on NVIDIA OpenShell with Nemoclaw, providing:

- **Conda/Miniforge3** — package management at `/sandbox/miniforge3`
- **build-env** — Conda environment with Python 3.13, C/C++ tooling, and research libraries
- **GCC/G++ 15.2** — via Conda target-prefixed compilers (`$CC`, `$CXX`)
- **CMake 4.3.2** — for C/C++ project builds
- **Research stack** — NumPy, SciPy, SymPy, Pandas, Matplotlib, BeautifulSoup4, LXML, python-docx, OpenPyXL
- **Code quality** — pytest, Ruff, Black, MyPy

### Setup

Always initialize the build environment before any work:

```bash
source /sandbox/setup_build.sh
```

This configures:
- Conda activation (`build-env`)
- C/C++ compiler paths (`$CC`, `$CXX`)
- GitHub CLI authentication
- Git user configuration

### Verification

```bash
source /sandbox/setup_build.sh
which python
python --version
echo "$CC"
echo "$CXX"
$CC --version
$CXX --version
cmake --version
pytest --version
```

## GitHub Authentication

- **User:** Sjbot0510
- **Token:** Managed via `/sandbox/secrets/sjbot0510_github_token`
- **Config Dir:** `/sandbox/gh-config-sjbot0510`
- **Git Config:** `/sandbox/gitconfig-sjbot0510`
- **Protocol:** HTTPS

## Constraints

- No `sudo apt install` available — use Conda (`conda-forge`) or pip only
- Conda requires `CONDA_NO_PLUGINS=true`, `CONDA_OVERRIDE_CUDA=""`, and `solver classic` (due to sandbox restrictions)
- For CMake projects, always pass `-DCMAKE_C_COMPILER="$CC" -DCMAKE_CXX_COMPILER="$CXX"`
