# Report React + D3 — abas procarióticas, pangenoma, diversidade e vista por amostra (plano 3 de 4)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** completar o report v2 até a paridade com o antigo — as abas Catálogo de MAGs, Defesa/AMR/Plasmídeos, Pangenoma, Diversidade e Leituras, mais a vista por amostra nos dois modos com a marca de procedência — e então remover o report antigo.

**Architecture:** fatia vertical, como no plano 2. Cada aba nasce com o bloco de dados que a alimenta e com as formas novas que ela exige, nunca antes. Nada de forma especulativa.

**Tech Stack:** React 18, D3 v7, esbuild, vitest + @testing-library/react, pytest.

**Spec:** `docs/superpowers/specs/2026-08-22-report-react-d3-design.md` (§5.4 a §5.8)

## Global Constraints

Valem inteiras as do plano 2, e mais três que são específicas do lado procariótico:

- **Herança é visível ou é mentira.** Sob o princípio (h) do `ROADMAP_SIMPLIFICACAO.md`, taxonomia GTDB, defesa, AMR, plasmídeo, KEGG e CAZy de um MAG foram computados no **representante** do cluster, que pode ter vindo de outra amostra. Todo valor herdado carrega marca e o representante no tooltip.
- **Nunca cortar um ID do catálogo no primeiro `__`.** `S1__binette_bin1__k141_1_5` cortado ali vira `S1` e atribui o achado à AMOSTRA. As vistas do Snakemake já reescrevem o prefixo; o report herda essa garantia e não refaz o corte.
- **O `?` do pangenoma não é ausência.** Três estados categóricos, `?` hachurado, denominador de frequência excluindo `?` e declarado no tooltip.

---

### Task 1: Aba Catálogo de MAGs — CONCLUÍDA (a16451b)

**Files:**
- Modify: `scripts/report/renderer_v2.py`
- Create: `src/report-ui/src/charts/Scatter.jsx`
- Create: `src/report-ui/src/panels/MagCatalog.jsx`
- Create: `src/report-ui/test/mag.test.jsx`
- Create: `tests/test_report_v2_prokaryotic.py`
- Modify: `src/report-ui/src/App.jsx`, `src/report-ui/src/styles.css`

**Blocos:** qualidade (CheckM2 × GUNC com zonas MIMAG), estrutura do catálogo (tamanho dos clusters + proveniência), taxonomia GTDB com o seletor de rank, metabolismo (KEGG por módulo com `missing_ko` no tooltip; CAZy por classe).

- [x] Testes primeiro (pytest para o bloco `prokaryotic`, vitest para o painel)
- [x] `build_prokaryotic()` lendo `mag_catalog/` — nunca as vistas por amostra
- [x] `Scatter.jsx` com as linhas de corte MIMAG desenhadas
- [x] Painel, aba condicionada a `data.prokaryotic`, bundle recompilado e commitado

### Task 2: Aba Defesa, AMR e plasmídeos — CONCLUÍDA (4cfa63b)

Heatmap MAG × tipo de sistema ordenado pela taxonomia; ARGs por classe de droga a partir do consenso (`n_tools >= 2`) com o número de ferramentas como encoding secundário; a ligação plasmídeo–AMR–defesa em três formas (UpSet de MAGs, trilha genômica do contig, stat tile da fração), com a ressalva de que replicon em contig de MAG é evidência de origem, não prova de plasmídeo intacto.

### Task 3: Aba Pangenoma

Matriz gene × membro em três estados; core/shell/cloud por cluster; tabela de candidatos com o critério que qualificou cada um, PlasmidFinder marcado como sinal de mobilidade e não como critério.

### Task 4: Abas Diversidade e Leituras

Alfa em strip/ridgeline conforme n; PCoA com % de variância nos eixos; Simpson e Chao1 só quando vieram de contagens; trilha sylph com `host_source` visível e o aviso de que o espaço de IDs do sylph não conversa com o dos contigs montados.

### Task 5: Vista por amostra (modos individual e comparação)

Os seis blocos da §5.8 no modo individual; pequenos múltiplos alinhados no modo comparação. Marca de procedência em tudo que é herdado.

### Task 6: Inventário da matriz de status

Hoje só três regras per-sample entram em `load_tool_status`. Ampliar para o conjunto real de regras da rodada, mantendo `ok`/`skipped`/`failed`/`unknown` distintos.

### Task 7: Remoção do report antigo

Só depois da paridade, em commit separado: `scripts/report/components/*.js`, `assets/echarts.min.js`, `rule generate_report`, e a troca do alvo em `_t_report()`.

---

## Verificação final do plano 3

- [ ] `python -m pytest tests/ -q` e `npx vitest run` passam inteiros
- [ ] Report da rodada real abre com todas as abas
- [ ] Modo escuro conferido em cada aba nova
- [ ] Nenhum valor herdado exibido sem marca de procedência
- [ ] Trilha desligada em `config.yaml` → estado vazio, nunca eixo quebrado
