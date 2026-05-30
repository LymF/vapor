# Contributing to VAPOR

## Releasing a new version

### 1. Update tool versions

Edit `containers.yaml` to bump tool versions, then regenerate the lock file:

```bash
conda activate snakemake
python3 scripts/pin_containers.py
git add containers.yaml containers.lock.yaml
```

### 2. Build and push custom Docker images

Three tools require custom images hosted on GHCR (GitHub Container Registry).
These must be rebuilt and pushed whenever their versions change.

**Genome maps** — pycirclize + matplotlib + biopython for circular genome visualisation:

```bash
docker build -f docker/Dockerfile.genome-map \
    -t ghcr.io/LymF/vapor-genome-map:1.0 .
docker push ghcr.io/LymF/vapor-genome-map:1.0
```

**GPU medaka** — ONT-official medaka image with CUDA support (optional):

```bash
docker build -f docker/Dockerfile.medaka-gpu \
    -t ghcr.io/LymF/vapor-medaka-gpu:2.2.0 .
docker push ghcr.io/LymF/vapor-medaka-gpu:2.2.0
```

**GPU COMEBin** — PyTorch CUDA image for GPU-accelerated binning (optional):

```bash
docker build -f docker/Dockerfile.comebin-gpu \
    -t ghcr.io/LymF/vapor-comebin-gpu:1.0.4 .
docker push ghcr.io/LymF/vapor-comebin-gpu:1.0.4
```

After pushing, update the version in `containers.yaml` under the `genome_map`,
`medaka`, or `comebin` key and re-run `pin_containers.py`.

### 3. Commit

```bash
git add containers.yaml containers.lock.yaml
git commit -m "Release vX.Y.Z: update container versions"
git tag vX.Y.Z
git push origin master --tags
```
