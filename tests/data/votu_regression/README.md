# Baseline de regressão do catálogo de vOTU

O teste `tests/test_votu_catalog_regression.py` roda sobre um recorte real da
corrida de junho/2026 e é pulado quando os dados não estão presentes.

Para gerar o insumo (requer skani >= 0.3.2 e os resultados em
`~/amazon/results`):

```bash
R=~/amazon/results
for s in P01_RNG_08_947 P01_RNG_08_948 P01_RNG_3_947 P01_RNG_3_948 \
         P06_TAP_3_957 P06_TAP_3_958; do
  awk -v S="$s" '/^>/{print ">"S"|"substr($1,2); next}{print}' \
      "$R/$s/viral/consensus/${s}_viral_nonredundant.fasta"
done > pool.fasta

skani triangle -i pool.fasta -o pool_sparse.tsv -t 8 --slow --sparse
```

Valores medidos em 2026-08-13 com skani 0.3.2:

| métrica | valor |
|---|---|
| sequências no pool | 9653 |
| vOTUs | 5524 |
| redução | 42,8 % |
| vOTUs em >1 amostra | 577 (10,4 %) |
