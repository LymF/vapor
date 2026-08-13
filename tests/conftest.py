import sys
from pathlib import Path

# Repo root (parent of tests/) on sys.path so `import pipeline_config` / `import vapor` work.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
