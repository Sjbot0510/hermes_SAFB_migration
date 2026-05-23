# hermes_SAFB_migration

Hermes Sandbox — Single Agent Framework Builder Migration

## Overview

This repository hosts the Hermes OpenShell sandbox environment configured for
research, prototyping, and migration work. It provides a fully containerized
development setup with Python, C/C++, and research tooling.

## Environment

| Component | Path / Details |
|---|---|
| Conda / Miniforge3 | `/sandbox/miniforge3` |
| Conda Environment | `build` |
| Python | 3.13 |
| GCC / G++ | 15.2 (Conda target-prefixed) |
| CMake | 4.3.2 |
| Git Config | `/sandbox/gitconfig-sjbot0510` |
| GH Config | `/sandbox/gh-config-sjbot0510` |
| Token Source | `/sandbox/secrets/sjbot0510_github_token` |

### Pre-installed Packages

- **Build/Dev:** gcc, gxx, cmake, make, gh, vim, ruff, black, mypy, pytest
- **Scientific:** NumPy, SciPy, SymPy, Pandas, Matplotlib
- **Web scraping:** BeautifulSoup4, LXML
- **Document handling:** python-docx, OpenPyXL, PyMuPDF, LatexCodec
- **Parsing:** pymupdf, pylatexenc, latexcodec, markdown, beautifulsoup4, lxml
- **Testing & linting:** pytest, Ruff, Black, MyPy

## Quick Start

### 1. Initialize the build environment

```bash
source /sandbox/setup_build.sh
```

> **Important:** Always use `source`, never `bash`, so shell variables persist.

### 2. Verify the environment

```bash
which python
python --version
echo "$CC"
echo "$CXX"
$CC --version
$CXX --version
cmake --version
pytest --version
```

### 3. Build a C/C++ project

```bash
mkdir build && cd build
cmake .. \
  -DCMAKE_C_COMPILER="$CC" \
  -DCMAKE_CXX_COMPILER="$CXX"
make
```

## GitHub Authentication

- **Username:** `Sjbot0510`
- **Email:** `sjbotchen@gmail.com`
- **Token:** Loaded from `/sandbox/secrets/sjbot0510_github_token`
- **Protocol:** HTTPS (`git_protocol https`)
- **Credential helper:** `gh auth git-credential`

```bash
source /sandbox/setup_build.sh
gh auth status          # verify authentication
gh repo create          # create a new repo
gh pr create            # open a PR
```

## Sandbox Constraints

1. **No `sudo apt install`** — use Conda (`conda install -c conda-forge`) or pip.
2. **Conda safety flags** must be set:
   - `CONDA_NO_PLUGINS=true`
   - `CONDA_OVERRIDE_CUDA=""`
   - `solver classic`
   
   These prevent crashes from CUDA auto-detection and libmamba in the protected
   sandbox environment.
3. **Compiler paths** — use `$CC` and `$CXX`, not bare `gcc`/`g++`.
4. **CMake** — always pass both `-DCMAKE_C_COMPILER="$CC"` and `-DCMAKE_CXX_COMPILER="$CXX"`.

## Directory Structure

```
/sandbox/
├── miniforge3/              # Conda installation
├── secrets/                 # Token files (gitignored)
├── setup_build.sh           # Environment initialization script
├── gh-config-sjbot0510/     # GitHub CLI config
└── gitconfig-sjbot0510      # Git global config
```

## Hermes Agent Integration

This sandbox runs Hermes Agent with the following capabilities:

- **Self-improving skills** — persistent procedural knowledge that accumulates over sessions
- **Multi-platform gateways** — Telegram, Discord, Slack, WhatsApp, Signal, Matrix, email, SMS, and 10+ more
- **Provider-agnostic** — 20+ LLM providers (OpenRouter, Anthropic, OpenAI, DeepSeek, local models, etc.)
- **Persistent memory** — cross-session recall with pluggable backends (built-in, Honcho, Mem0)
- **Profiles** — multiple independent instances with isolated configs, sessions, and memory
- **Cron jobs** — scheduled task execution with delivery to messaging platforms
- **Extensible toolsets** — web, browser, terminal, file, vision, image_gen, tts, MCP, delegation, and more

### Quick Start with Hermes

```bash
# Initialize environment first
source /sandbox/setup_build.sh

# Interactive chat
hermes

# Single query
hermes chat -q "What is the capital of France?"

# Check health
hermes doctor

# View/configure tools
hermes tools
hermes config

# Skills management
hermes skills list
hermes skills search <query>
hermes skills install <id>

# Models and providers
hermes model
hermes auth add
```

### Key Paths

| Path | Purpose |
|------|---------|
| `~/.hermes/config.yaml` | Main configuration |
| `~/.hermes/.env` | API keys and secrets |
| `~/.hermes/sessions/` | Session transcripts |
| `~/.hermes/logs/` | Gateway and error logs |
| `~/.hermes/auth.json` | OAuth tokens and credential pools |
| `$HERMES_HOME/skills/` | Installed skills |

## Troubleshooting

| Problem | Solution |
|---|---|
| `gcc: command not found` | Use `$CC --version` instead |
| Conda crashes on install | Verify `CONDA_NO_PLUGINS=true` and `solver classic` |
| GitHub auth fails | Run `source /sandbox/setup_build.sh` first |
| CUDA-related errors | Ensure `CONDA_OVERRIDE_CUDA=""` is set |
| Hermes tools not available | Run `hermes tools`, enable required toolsets, then `/reset` |
| Model/provider issues | Run `hermes doctor` and check `.env` |
| Gateway dies on SSH logout | Enable linger: `sudo loginctl enable-linger $USER` |
| Gateway dies on WSL2 close | Add `systemd=true` to `/etc/wsl.conf` |
