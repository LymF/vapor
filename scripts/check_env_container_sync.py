#!/usr/bin/env python3
"""
check_env_container_sync.py — garante que os dois caminhos de execução rodam
a mesma versão de cada ferramenta.

A vapor resolve conda e container SEPARADAMENTE:

    rule checkv:
        conda:     "../envs/env_viral.yaml"      # pin do env
        container: CONTAINERS.get("checkv")      # pin do containers.yaml

Nada no Snakemake garante que os dois pins concordem. Se divergirem, o mesmo
commit produz resultados diferentes conforme a flag usada
(`--sdm conda --use-conda` vs `--sdm apptainer`) — sem erro, sem aviso.

Este script casa cada regra que declara AMBOS e compara a versão do pacote
no env conda com a versão declarada em containers.yaml.

Uso:
    python3 scripts/check_env_container_sync.py            # relatório
    python3 scripts/check_env_container_sync.py --check    # exit 1 se divergir (CI)

O que é aceito sem reclamar:
  - entradas `custom: true` em containers.yaml (imagens próprias no GHCR, sem
    equivalente bioconda — genome_map, medaka-gpu)
  - pacotes instalados por `pip:` no env (comparados normalmente)
"""

import argparse
import glob
import os
import re
import sys

import yaml

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def env_pins(path):
    """{pacote: versão} de um env yaml, incluindo a seção pip."""
    data = yaml.safe_load(open(path)) or {}
    pins = {}
    unpinned = set()
    for dep in data.get("dependencies") or []:
        if isinstance(dep, str):
            _collect(dep, pins, unpinned)
        elif isinstance(dep, dict):
            for pdep in dep.get("pip") or []:
                if isinstance(pdep, str):
                    _collect(pdep.split("#")[0].strip(), pins, unpinned)
    return pins, unpinned


def _collect(dep, pins, unpinned):
    dep = dep.split("#")[0].strip()
    if not dep:
        return
    # pip usa '=='; conda usa '='. Faixas (>=, <=, ~=) NÃO são pin exato.
    m = re.match(r"^([A-Za-z0-9._-]+)\s*==?\s*([0-9][^\s=]*)$", dep)
    if m:
        pins[m.group(1)] = m.group(2)
        return
    m = re.match(r"^([A-Za-z0-9._-]+)\s*[><~!]=", dep)
    if m:
        unpinned.add(m.group(1))


def rule_pairs():
    """(regra, arquivo_env, chave_container) para regras que declaram os dois."""
    pairs = set()
    for f in sorted(glob.glob(os.path.join(ROOT, "rules", "*.smk"))):
        txt = open(f, "rb").read().decode()
        for block in re.split(r"\n(?=\s*(?:use )?rule [A-Za-z0-9_]+:)", txt):
            m = re.match(r"\s*(?:use )?rule ([A-Za-z0-9_]+):", block)
            if not m:
                continue
            head = re.split(r"\n\s*(?:shell|run):", block)[0]
            env = re.search(r'conda:\s*"\.\./envs/([^"]+)"', head)
            key = re.search(r'container:\s*CONTAINERS\.get\("([^"]+)"\)', head)
            if env and key:
                pairs.add((m.group(1), env.group(1), key.group(1)))
    return sorted(pairs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true",
                    help="exit 1 se houver divergência (para CI)")
    args = ap.parse_args()

    containers = yaml.safe_load(open(os.path.join(ROOT, "containers.yaml"))) or {}
    envs = {}
    for f in glob.glob(os.path.join(ROOT, "envs", "*.yaml")):
        envs[os.path.basename(f)] = env_pins(f)

    mismatched, unpinned, missing = [], [], []
    checked = 0

    for rule, envf, key in rule_pairs():
        spec = containers.get(key) or {}
        if spec.get("custom"):
            continue  # imagem própria, sem par bioconda
        pkg = spec.get("image", key)
        cver = str(spec.get("version", ""))
        pins, unp = envs.get(envf, ({}, set()))
        if pkg in pins:
            checked += 1
            if pins[pkg] != cver:
                mismatched.append((rule, envf, pkg, cver, pins[pkg]))
        elif pkg in unp:
            unpinned.append((rule, envf, pkg, cver))
        else:
            missing.append((rule, envf, pkg, cver))

    print(f"{checked} pares regra/ferramenta comparados\n")

    if mismatched:
        print("DIVERGENTES (mesma ferramenta, versões diferentes por caminho):")
        print(f"  {'regra':<26}{'env':<26}{'pacote':<18}{'container':<12}{'env'}")
        for r, e, p, cv, ev in mismatched:
            print(f"  {r:<26}{e:<26}{p:<18}{cv:<12}{ev}")
        print()
    if unpinned:
        print("SEM PIN EXATO no env (o caminho conda pode resolver outra versão):")
        for r, e, p, cv in unpinned:
            print(f"  {r:<26}{e:<26}{p:<18}container={cv}")
        print()
    if missing:
        print("Pacote ausente do env (informativo — pode vir por dependência):")
        for r, e, p, cv in missing:
            print(f"  {r:<26}{e:<26}{p:<18}container={cv}")
        print()

    problems = len(mismatched) + len(unpinned)
    if not problems:
        print("OK — conda e container pinados na mesma versão em todas as regras.")
    if args.check and problems:
        print(f"FALHOU: {len(mismatched)} divergência(s), {len(unpinned)} sem pin exato.",
              file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
