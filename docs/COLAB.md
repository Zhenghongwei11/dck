# Colab Workflow (Recommended)

Goal: run memory‑heavy steps in Colab while keeping the repository fully reproducible for third‑party reviewers.

## Key principle

Google Drive can be used as a **personal cache** (persistence between sessions), but the **public reproduction path must not depend on your Drive**. Auditors should be able to run everything from public sources listed in `data/manifest.tsv`.

## Option A (simplest): Colab + Google Drive cache

1. Create a new Colab notebook (or use VS Code + Colab kernel).
2. Mount Drive:

```python
from google.colab import drive
drive.mount("/content/drive")
```

3. Choose a working directory on Drive:

```python
import os
WORKDIR = "/content/drive/MyDrive/ezhu"
os.makedirs(WORKDIR, exist_ok=True)
os.chdir(WORKDIR)
```

4. Clone the repo (or `git pull` if already cloned):

```python
!git clone https://github.com/<ORG>/<REPO>.git
%cd <REPO>
```

5. Run the pipeline:

```python
!make setup
!make reproduce
```

Notes:
- If you do not want to store large raw data in Drive, keep `data/raw/` under `/content` (ephemeral) and only copy `results/` + `plots/` back to Drive at the end.
- Colab VM storage is temporary and will reset; Drive persists but has your account quota.

## Option B (auditor‑friendly): Colab without Drive dependency

This is the recommended public reproduction approach:
- Do not require Drive at all.
- Download everything from authoritative public URLs in `data/manifest.tsv` and validate checksums.
- Run `make reproduce`.

Reviewers can reproduce from scratch without needing access to your Drive.

