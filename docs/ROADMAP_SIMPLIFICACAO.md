# Roadmap — simplificação da vapor

**Criado:** 2026-08-18
**Motivação:** a vapor acumulou ferramentas redundantes. O objetivo desta linha de
trabalho é **remover complexidade**, não adicionar, e só adicionar onde houver
lacuna funcional real.

Este documento existe para uma sessão nova conseguir continuar sem redescobrir as
decisões. Leia as seções "Decisões já tomadas" e "Achados que não devem ser
redescobertos" ANTES de propor alternativas — várias já foram avaliadas e
descartadas com motivo.

---

## Estado atual

| item | estado |
|---|---|
| status real no `done.txt` (15 regras) | **feito** — commit `4ff4507`, já no master |
| Remoção do VIBRANT | **feito** — commit `d8dbf7d` |
| (f) `bacphlip_votu` + (g) `eggnog_viral` | **feito** — commit `bf2d358`, regras globais |
| Relatório religado (Lifestyle + Putative AMGs) | **feito** — commit `1c66c6c` |
| (b) metaSPAdes + metaviralSPAdes | **feito** — só o MEGAHIT nas short reads |
| (c) um montador só para long reads | **feito** — Flye+Medaka (ONT) / metaMDBG (HiFi) |
| (h) computar no representante | **feito** — commit `94a9bf9`, mais 4 bugs de namespace |
| (d) remover merge_contigs / mmseqs2 | **feito** — o hub agora é a montagem |
| (e) mover o limiar de comprimento | **feito** — portão composto + unificação da cadeia skani de grupo + FASTA/TSV de descarte |
| Restante | nenhum item aberto nesta lista |

Branch de trabalho: `master` (o `refactor/unit-wildcard` já foi mergeado por
fast-forward e pode ser apagado).

---

## Sequência acordada

**As letras são identidade, não ordem.** Elas foram atribuídas na ordem em que
os itens surgiram na discussão; a ordem de execução foi revisada depois e está
abaixo. Não execute em ordem alfabética.

### Ordem de execução

```
(a) VIBRANT ........................ FEITO  d8dbf7d
(f) bacphlip_votu .................. FEITO  bf2d358
(g) eggnog_viral ................... FEITO  bf2d358
        ↓
(b) remover metaSPAdes + metaviralSPAdes ... FEITO
(c) um montador só para long reads ....... FEITO
        ↓
(h) migrar regras para o representante ... FEITO
        ↓
(d) remover merge_contigs / merge_lr / mmseqs2 ... FEITO
(e) mover o limiar de comprimento .......... FEITO
```

**Por que (h) vem antes de (d).** As duas mexem em *quem consome o quê*, e
incidem sobre o mesmo conjunto de pontos de referência — as 22 referências ao
`rep_seq.fasta` e as 7 regras que consomem `votu_catalog_reps` se sobrepõem.
Fazer as duas em sequência sobre os mesmos consumidores dobra a chance de erro
silencioso, que é exatamente a classe de falha que já custou caro aqui (ver
"Lição registrada sobre método").

Fazendo (h) primeiro, as regras que rodam no representante global **saem do
caminho** do hub `rep_seq.fasta` antes de ele ser desmontado. O (d) fica com um
conjunto menor e mais homogêneo de consumidores para reapontar.

(b) e (c) vêm antes das duas porque são pré-requisito de (d): o dedup só pode
sair depois que sobrar **um montador por trilha**, que é o que remove a razão de
existir dele.

Cada passo é verificado antes do seguinte (ver "Como verificar").

### (a) Remover VIBRANT — FEITO (commit `d8dbf7d`)

Motivo: parte relevante do falso positivo viral vem dele.

Superfície: `rules/viral_detection.smk` (regra + parsing no `viral_consensus`),
`rules/coassembly.smk` (gêmea), `rules/viral_binning.smk` (`make_votu_table`),
`scripts/make_votu_table.py`, `scripts/report/{data_loaders,renderer}.py`,
`scripts/report/components/{overview,viral,annotation,about}.js`, `Snakefile`
(`_VIBRANT_BASE` + alvo no `rule all`), `config.yaml`, `vapor.py`,
`envs/phage_vibrant.yaml` (apagar), `envs/env_viral.yaml`, `containers.yaml`,
`containers.lock.yaml`, docs.

**Perdas funcionais aceitas** (o VIBRANT era fonte única das duas):
- detecção de AMG → reposta em (g)
- lifestyle lítico/lisogênico dos vOTUs montados → reposta em (f)

Decisão de schema: as colunas de lifestyle e AMG do `make_votu_table.py`
**continuam existindo, preenchidas como vazio**, para não quebrar quem consome o
arquivo e para religar a fonte depois sem migração.

### (b) Remover metaSPAdes e metaviralSPAdes — FEITO

Fica só o MEGAHIT nas short reads.

Consequência: `rule cobra_spades` morre junto (o COBRA tem uma regra por
montador). `rules/quast.smk` e o relatório referenciam os montadores.

Precedente: o metaFun testou MEGAHIT vs metaSPAdes em 3 ambientes CAMI e ficou
só com o MEGAHIT — venceu em comprimento total e fração de bp alinhados, com
desempate por eficiência de recursos.

**O que saiu:** regras `metaspades`, `metaviral_spades` e `cobra_spades`; chaves
`use_spades`, `spades_mem`, `spades_kmers`, `spades_kmer_list`,
`cobra_spades_mink/maxk`; `spades=4.2.0` do `env_assembly`; entrada `spades` do
`containers.yaml`/`.lock.yaml`; rótulos do QUAST (5 → 3 conjuntos) e do relatório.

**O que ficou de propósito:** o `merge_contigs` continua existindo e continua
prefixando `MEGAHIT_` — ele só morre em (d), e manter o prefixo evita renomear
contig em todo o resultado já produzido. O `cobra_merge` também ficou, agora com
um diretório de corrida só.

**Verificação:** dry-run PE 195 → 189 jobs, 331 → 323 outputs. Delta = 2 amostras
× 4 arquivos (metaspades `contigs.fasta` + `.gfa`, metaviral `contigs.fasta`,
`cobra/spades/done.txt`) e 2 × 3 regras. Nada mais mudou. Dry-run também limpo
em SE (189 jobs) e long reads (167 jobs).

### (c) Um montador só para long reads — FEITO

Escolha por `lr_tech`: **Flye para ONT, metaMDBG para HiFi**. O `hifiasm_lr` sai
sempre; o `merge_lr` morre.

**O que saiu:** regras `hifiasm_lr` e `merge_lr`; o ramo hifiasm do `medaka_lr`
(output `hifiasm_pol`); o guard `USE_METAMDBG` e o ramo `--in-ont` do
`metaMDBG_lr`; constantes `LR_HIFIASM_HOM`, `LR_METAMDBG` e `USE_METAMDBG` do
`Snakefile`; chaves `lr_hifiasm_hom` e `lr_metaMDBG` do `config.yaml`;
`hifiasm`/`hifiasm_meta` do `env_flye`; entrada `hifiasm` dos containers;
`scripts/merge_lr_assemblies.py` (apagado, sem consumidor).

**Estrutura nova:** dentro de `if LONG_READS:` as regras passam a ser
instanciadas por `if LR_TECH == "ont":` / `else:`, então nunca existe montador
no DAG que não vai ser usado. O `mmseqs2` consome:
- ONT → `{sample}/assembly/lr/flye_polished/assembly.fasta`
- HiFi → `{sample}/assembly/lr/metaMDBG/assembly.fasta`

**Prefixo de header:** o `merge_lr` prefixava `FLYE_`/`HIFIASM_`/`MDBG_`. Com um
montador só o prefixo some para long reads. Verificado por grep que **ninguém
faz parsing dele** em `scripts/`, `rules/` ou `scripts/report/` — e o
`votu_catalog_pool` já prefixa por origem no catálogo global.

**Alvo do `rule all` reapontado:** o `Snakefile` ainda pedia
`assembly/lr/merged/done.txt` (o sentinel do `merge_lr`). Passou a pedir o
sentinel do montador único — `medaka_done.txt` (ONT) ou `metaMDBG/done.txt`
(HiFi). Esse era o único ponto que o grep por "merge_lr" não pegava, porque a
referência era um caminho literal.

**Verificação:** dry-run limpo nos quatro modos.

| modo | jobs antes | jobs depois |
|---|---|---|
| PE  | 189 | 189 (diff de outputs **vazio**) |
| SE  | 189 | 189 |
| LR ONT | 167 | 161 |
| LR HiFi | — | 159 |

Delta LR-ONT = 6 = 2 amostras × 3 regras (`hifiasm_lr`, `merge_lr`,
`metaMDBG_lr`); o conjunto de regras do DAG não mudou em mais nada. Short reads
intactas, como exigido.

**Bug pré-existente encontrado e corrigido (fora de (c)):** a trilha **HiFi não
rodava desde antes deste roadmap**. Em `rules/qc.smk` o `filtlong_lr` recebia,
quando `lr_tech != "ont"`, o caminho
`{sample}/lr_filtered/{sample}_placeholder.fastq.gz` — **arquivo que nenhuma
regra produz** → `MissingInputException` antes mesmo da montagem. Como o
`porechop_lr` já trata HiFi internamente (faz `cat` das reads cruas), o
`filtlong_lr` passou a consumir sempre a saída dele: uma trilha só para as duas
tecnologias, sem condicional. O condicional antigo tinha ainda um segundo
problema — o `and LONG_READS` era morto, já que a regra só existe dentro do
`if LONG_READS:`.

Depois do conserto o HiFi monta o DAG: **159 jobs**, com
`porechop_lr` (pass-through) → `filtlong_lr` → `metaMDBG_lr`, sem Flye nem
Medaka. ONT e short reads não mudaram nada.

### (d) Remover `merge_contigs`, `merge_lr` e `mmseqs2` (dedup) — FEITO

**Fechamento (2026-08-18).** O hub passou a ser a montagem, resolvida pelo
helper `_sample_contigs()` do `Snakefile` (mesmo estilo dos `_clean_r1`/
`_clean_lr` que já existiam) — e **não** por uma regra que copiasse o FASTA
para um caminho canônico, que duplicaria um arquivo enorme em disco por
amostra.

Eram **13** pontos, não 22: o (c) já havia eliminado parte quando o `merge_lr`
morreu. O `{sample}_cluster.tsv` do MMseqs2 **não tinha nenhum consumidor** —
morreu sem substituto.

O prefixo `MEGAHIT_` caiu junto com o `merge_contigs`; os contigs passam a ter
o nome cru do montador (`k141_10`), como a trilha de co-assembly sempre teve.
Auditado ponto a ponto: o único lugar que re-deriva chave a partir do header da
montagem é o `multisplit_catalog`, que já re-namespaceia tudo por conta própria
(`>S{sample}C{header}`). O resto ou lê o catálogo global (namespaceado por
origem, indiferente ao nome cru) ou trabalha dentro do espaço de nomes de uma
amostra só.

`min_seq_id`/`MIN_SEQ_ID` ficaram órfãos e saíram. O MMseqs2 **permanece** no
`env_assembly` e no `containers.yaml` — as regras de taxonomia usam.

**O "ganho colateral" do COBRA não existia.** O roadmap supunha que o COBRA
pudesse estar recebendo montagem diferente daquela contra a qual o BAM foi
mapeado. Verificado: antes do (d) o `bwa_index` e o `cobra_megahit` já
apontavam para o mesmo `rep_seq`; depois, ambos apontam para a mesma montagem.
O pareamento já estava correto — o (d) preserva, não conserta.

**Ressalva sobre comparabilidade (achado da auditoria).** A justificativa do
roadmap — *"o `easy-linclust` existia para colapsar redundância entre
montadores"* — é **parcial**. Com `--min-seq-id 0.95 -c 0.9 --cov-mode 2` ele
também removia contigs *contidas* dentro de outras da **mesma** montagem
(variantes de cepa, bolhas não resolvidas). Efeito: mais contigs redundantes na
detecção viral e no mapeamento, com multi-mapping um pouco maior e abundância
levemente diluída. No viral o catálogo global de vOTU (skani 95/85, mais
estrito) absorve; no binning, MetaBAT2/SemiBin2 absorvem. Impacto prático
pequeno, **mas resultados pós-(d) não são diretamente comparáveis aos
pré-(d)** — considerar ao reprocessar dataset já analisado.

**Bug que só a auditoria pegaria.** O `scripts/report/renderer.py` lia
`quast_data[s].get("deduplicated")`, e o rótulo do QUAST virou `"assembly"`.
O `parse_quast_all` chaveia pelo header do `report.tsv`, então `# contigs` e
`N50` de todas as amostras sairiam como "N/A" no card de Overview — sem erro,
sem log. **Nenhum dry-run pegaria isso: o relatório não roda em `-n`.** Vale
como lembrete de que o invariante do DAG não cobre o conteúdo do relatório.

**Verificação:** dry-run limpo nos 4 modos. PE 177 → 173, SE 177 → 173,
LR-ONT 149 → 147, LR-HiFi 147 → 145. Delta = `merge_contigs` (só SR) +
`mmseqs2` (todos), 2 amostras cada. Um grep por `mmseqs/` no dry-run inteiro
devolve só caminhos de taxonomia — nenhum `_rep_seq.fasta` em nenhum modo.

> **Pré-requisitos: (b), (c) e (h).** (b)/(c) porque o dedup só pode sair
> depois de sobrar um montador por trilha. (h) porque as regras que rodam
> no representante global precisam sair do caminho do hub `rep_seq.fasta`
> antes dele ser desmontado — ver "Ordem de execução".

**Esta é a mudança de maior risco.** O hub deixa de ser
`{sample}_rep_seq.fasta` e passa a ser o contig do MEGAHIT — 22 referências a
reapontar (`viral_detection`, `mapping`, `prok_binning`, `cobra`, `quast`,
`coassembly`, `Snakefile`).

Justificativa: o `easy-linclust` existia para colapsar redundância **entre
montadores**. Com um montador por trilha essa fonte deixa de existir.

**Precedente decisivo:** a trilha de co-assembly **já roda assim** — MEGAHIT →
`contigs.fa` → tudo, sem merge e sem dedup — e está em produção
(`results/coassembly/RNG/` tem corrida completa com bins, checkm2, gtdbtk e
viral). O padrão de destino já existe e está testado; é reapontar, não inventar.

Ganho colateral: sem dedup, o COBRA passa a receber a mesma montagem contra a
qual o BAM foi mapeado, que é o que ele exige.

### (e) Mover o limiar de comprimento em vez de apagá-lo — FEITO

**A tabela original desta seção (abaixo, riscada) propunha um `viral_min_contig`
de 3000 aplicado ANTES da detecção viral. Substituída em 2026-08-18 por um
desenho melhor, decidido com base em literatura — não reabrir sem novo motivo.**

> ~~| onde | hoje | destino |~~
> ~~|---|---|---|~~
> ~~| MEGAHIT | `--min-contig-len 3000` | **1000** |~~
> ~~| `merge_contigs` | filtra 3000 de novo | removido em (d) |~~
> ~~| detecção viral | **sem piso** | **`viral_min_contig`, 3000** (parâmetro novo) |~~
> ~~| binners | clampam sozinhos | inalterado |~~
> ~~| COBRA | recebe filtrado | recebe as curtas |~~

**O desenho final:** o piso de comprimento só vale onde não há evidência
independente de que a sequência é um genoma real. O filtro roda **DEPOIS do
vRhyme**, não antes, e é **composto** — uma sequência é mantida se qualquer
condição valer:

1. está num bin do vRhyme (o bin sustenta o contig curto), **ou**
2. o CheckV a classifica como **Complete / High-quality / Medium-quality**
   (ou completeness ≥ 50%). **Atenção:** os tiers aqui são o conjunto FIXO
   `MQ_TIERS` do `scripts/viral_length_gate.py`, e deliberadamente **não** o
   `VIRAL_KEEP_TIERS` configurável que o `is_mq()` do `votu_catalog.smk` usa.
   Os dois só coincidem com `viral_min_quality: medium`; com o `not_determined`
   que a config traz, o `VIRAL_KEEP_TIERS` expande para os cinco tiers do CheckV
   e este braço ficaria sempre verdadeiro — o portão inteiro viraria no-op.
3. tem **≥ 5000 bp** (`viral_min_contig`).

Racional, com literatura: é a mesma regra do **MVP** (Coclet, Camargo & Roux,
2024, *mSystems* 9:e00888-24, https://doi.org/10.1128/msystems.00888-24),
pipeline do grupo que escreveu o geNomad, o CheckV e o IMG/VR —
*"selects low-quality genomes larger than 5 kb or complete, high-, or
medium-quality and larger than 1 kb"*. O 5 kb é o corte do IMG/VR e do Earth's
Virome Protocol (Roux et al., 2021, *NAR* 49:D764,
https://doi.org/10.1093/nar/gkaa946). O 3000 da tabela antiga não tinha
respaldo na literatura — era número interno.

**VirSorter2 e geNomad continuam vendo tudo desde `min_contig` (agora 1000
bp)** — o filtro NÃO roda antes da detecção viral. Custo em compute conhecido
e aceito: baixar `min_contig` para 1000 multiplica a contagem de contigs, o
mapeamento fica bem mais pesado.

O que mudou, por peça:

| onde | hoje | antes |
|---|---|---|
| MEGAHIT / COBRA / QUAST | `--min-contig-len {MIN_CONTIG}`, **1000** | 3000 |
| `merge_contigs` | removido em (d) — peso morto, MEGAHIT já filtrava | filtrava 3000 de novo |
| detecção viral (VS2/geNomad) | sem piso, vê tudo desde 1000 bp | sem piso, via 3000 bp |
| `viral_nonredundant` (per-sample, pós-vRhyme) | **portão composto** (bin OU qualidade CheckV OU `viral_min_contig=5000`) | nada — todo não-binado passava |
| `coassembly_viral_nonredundant` (grupo short-read, pós-vRhyme) | **mesmo portão composto de 3 braços** — paridade com o per-sample | regra não existia; grupo nunca chegava ao catálogo via bins |
| grupo long-read (sem vRhyme de grupo) | segue vindo do `coassembly_viral_trimmed` pré-binning, sem portão | idem |
| binners | clampam sozinhos (vRhyme 2000, MetaBAT2 1500) | inalterado |

Implementação: `scripts/viral_length_gate.py` (função pura, testada por
unittest fora do Snakemake) é compartilhada por `rule viral_nonredundant`
(`rules/viral_binning.smk`) e `rule coassembly_viral_nonredundant`
(`rules/coassembly.smk`). A armadilha corrigida junto: `coassembly_vrhyme`
passava `-l {MIN_CONTIG}` **sem** o clamp de 2000 que a regra por amostra já
tinha, e engolia falha com `|| true` + `touch done.txt` incondicional — com
`MIN_CONTIG=1000` isso faria a trilha viral de grupo produzir zero bins em
silêncio. Corrigido: mesmo clamp `max(MIN_CONTIG, 2000)` e status real
(`ok`/`failed: vRhyme exit $RC`) em vez do `touch` cego.

**Follow-up 2026-08-18 (mesma sessão): duas pendências fechadas.**

1. **Cadeia de clustering interna do grupo unificada.** `coassembly_skani_votu`
   / `coassembly_skani_cluster` / `coassembly_viral_votu_reps` ainda liam
   `coassembly_viral_trimmed` (pré-vRhyme, pré-portão) mesmo depois de
   `_catalog_sources()` (`rules/votu_catalog.smk`) já ter sido movido para o
   fasta pós-portão. Reapontadas para `_coas_skani_source`, uma variável
   módulo definida uma vez em `rules/coassembly.smk`: para grupos short-read
   é `coassembly_viral_nonredundant.output.fasta` (o mesmo que o catálogo
   global usa); para long-read (sem vRhyme de grupo) continua
   `coassembly_viral_trimmed.output.fasta`, comportamento preservado. Exigiu
   mover o bloco `coassembly_vrhyme` / `coassembly_checkv_vrhyme` /
   `coassembly_viral_nonredundant` para ANTES do bloco skani no arquivo —
   `rules.coassembly_viral_nonredundant` só existe depois de definido, e os
   dois blocos tinham condições de guarda (`if`) diferentes (`not LONG_READS`
   vs sempre), então reabrir o `if COASSEMBLY_ENABLED and COASSEMBLY_VIRAL:`
   depois do bloco `not LONG_READS` foi necessário para não prender
   skani_votu/cluster/viral_votu_reps/viral_taxonomy dentro do guard errado.
   `coassembly_organize_outputs` também copiava, sob o nome
   `viral_nonredundant.fasta` em `final/`, o arquivo ERRADO
   (`{group}_viral_trimmed.fasta`, pré-portão) para grupos short-read —
   corrigido para copiar `{group}_viral_nonredundant.fasta` (pós-portão)
   nesse caso, mantendo o trimmed só para long-read.

2. **FASTA + TSV de descarte.** O que o portão composto derruba deixou de
   simplesmente desaparecer. `rule viral_nonredundant` e
   `rule coassembly_viral_nonredundant` agora também escrevem
   `{sample|group}_viral_discarded.fasta` (header = `contig_id` puro, SEM
   motivo codificado nele — decisão explícita do usuário: um header anotado
   quebra qualquer join de ferramenta downstream que rode sobre o conjunto
   descartado) e `{sample|group}_viral_discarded.tsv` (uma linha por
   sequência: `contig_id, length, checkv_quality, checkv_completeness,
   in_vrhyme_bin, source_id` — schema em
   `scripts/viral_length_gate.py::DISCARD_TSV_COLUMNS`, idêntico nas duas
   trilhas). `checkv_quality`/`checkv_completeness` ficam vazios (não `0`)
   quando o CheckV nunca pontuou o contig — distinguir "nunca avaliado" de
   "avaliado como zero" é o próprio motivo de existir do TSV: um
   "Not-determined" curto é exatamente a cara de um vírus genuinamente novo.
   Entregue em `final/viral/viral_discarded.{fasta,tsv}` (per-sample e por
   grupo), ao lado de `viral_nonredundant.fasta`.

### (f) `bacphlip_votu` — repor lifestyle — FEITO

O BACPHLIP **já está instalado** (`envs/env_reads_classify.yaml`, via pip) e o
script `scripts/reads_classify/bacphlip_lifestyle.py` já existe — mas hoje recebe
genomas de referência do sylph, não contigs montados.

Regra nova: entra o FASTA viral filtrado por qualidade CheckV, sai lifestyle;
religa a coluna do `make_votu_table`.

**Parâmetro `bacphlip_min_quality`, default `HQ+Complete`.** Reusar o vocabulário
que já existe em `VIRAL_KEEP_TIERS` (`rules/votu_catalog.smk:202` — o tier `mq`
já está definido como "Complete/HQ/MQ or completeness >= 50%").

**Por que o default é conservador:** o BACPHLIP decide pela presença de domínios
associados à lisogenia (integrase, recombinase, repressor). Em MQ falta até
metade do genoma; se o pedaço ausente carregava a integrase, um fago lisogênico
é chamado de virulento. **O erro não é simétrico — a fragmentação enviesa na
direção de "virulento"** e infla artificialmente o percentual de líticos. Rodar
em MQ é permitido, mas é escolha consciente do usuário. Registrar em coluna qual
tier foi usado.

### (g) `eggnog_viral` — repor AMG como *putative* — FEITO

**Nada a instalar.** A vapor já tem as duas peças:
- `rule prodigal_viral` (`taxonomy.smk:76`) já produz os ORFs virais (`.faa`)
- eggNOG-mapper 2.1.15 já instalado e `eggnog_db` já baixado (`config.yaml:179`),
  com KEGG, CAZy e GO

AMG, na definição operacional do VIBRANT e do DRAM-v, é **gene metabólico KEGG
em genoma viral**. Rodar eggNOG sobre o `.faa` viral e filtrar KOs metabólicos
entrega o mesmo tipo de candidato com zero DB novo.

**Nomear como "putative AMGs" no relatório e nas colunas, nunca como "AMGs"
cravado** — decisão explícita do usuário. Justificativa: existe literatura de
2025 alertando que chamada de AMG sofre muita má-anotação e exige inspeção
manual do contexto genômico. Isso vale igual para DRAM-v e VIBRANT: não estamos
trocando método curado por bruto, estamos trocando um não-curado por outro mais
barato.

### (h) Princípio: computar no representante, herdar no membro

> **Executar ANTES de (d).** As duas mexem nos mesmos consumidores; fazer
> (h) primeiro reduz o conjunto que (d) precisa reapontar — ver
> "Ordem de execução".

**Esta é a mudança estruturalmente mais importante do roadmap.** Não é um
conserto pontual — é o modelo que o resto das regras deveria seguir.

#### O princípio

Um vOTU é um agrupamento a 95% ANI / 85% AF — nível de espécie pelo padrão ICTV.
Todo membro de um vOTU é, por construção, a mesma entidade biológica que o
representante. Logo:

> **Toda análise que depende só da sequência roda UMA vez, no representante
> global. O membro herda o resultado por join na provenance.**

Só permanece por amostra o que depende de **dados daquela amostra** — não da
sequência.

| natureza da análise | onde roda | exemplos |
|---|---|---|
| depende só da sequência | **uma vez, no representante global** | predição de ORF, anotação funcional, taxonomia, lifestyle, AMG, qualidade |
| depende de reads/bins da amostra | **por amostra** | recrutamento, cobertura, abundância, presença/ausência, predição de hospedeiro |

A tabela por amostra deixa de ser um recálculo e passa a ser uma **visão**: junta
a anotação global (via `provenance.tsv`) com a abundância local.

#### Por que isso é possível hoje

A parte difícil já está pronta. O `votu_catalog_pool` prefixa cada contig com sua
origem (*"Contig IDs are only unique within an assembly, so they are prefixed
with their source"*) e emite `provenance.tsv`. **Pool global com rastreio de
origem já existe** — é exatamente o que essa arquitetura precisa.

As regras `bacphlip_votu` e `eggnog_viral` (itens f e g) **já nasceram assim** e
servem de molde: global, no representante, com o membro herdando no
`make_votu_table.py`.

#### O ganho

Com 32 amostras, cada regra migrada deixa de rodar 32 vezes e passa a rodar 1.
Para `pharokka`/`phold` isso é a diferença entre horas e minutos. E elimina uma
classe inteira de inconsistência: hoje nada garante que duas amostras anotem o
mesmo representante do mesmo jeito.

#### Estado verificado (2026-08-18)

Sete regras por amostra consomem `votu_catalog_reps`. Verificado input a input.
**Cuidado com o método:** referências `rules.X.output` são por amostra sem conter
a string `{sample}` — uma checagem por texto literal dá falso positivo. Foi assim
que a primeira contagem desta análise saiu errada.

O `CLAUDE.md` registra que o catálogo global de vOTU *"Replaces the former
per-sample clustering"*. A parte difícil está pronta: o `votu_catalog_pool`
prefixa cada contig com sua origem (*"Contig IDs are only unique within an
assembly, so they are prefixed with their source"*) e emite `provenance.tsv` —
ou seja, **pool global com rastreio de origem já existe**.

Mas nem toda regra a jusante foi movida junto. Sete regras por amostra consomem
`votu_catalog_reps`. Verificado input a input (cuidado: referências
`rules.X.output` são por amostra sem conter a string `{sample}` — uma checagem
por texto literal dá falso positivo):

| regra | inputs | veredito |
|---|---|---|
| `prodigal_viral` | **só o FASTA global** | **N× idêntico — desperdício confirmado** |
| `pharokka` | global + CheckV **por amostra** | **BUG confirmado — não anota nada** |
| `viral_taxonomy` | + geNomad/MMseqs por amostra | **BUG confirmado — geNomad sai do merge** |
| `genome_map_virus` | + CheckV por amostra | **BUG confirmado — zero mapas** |
| `genome_map_phage` | + phold/pharokka por amostra | quebrado em cascata pelo `pharokka` |
| `phist` | + `bins_dir` por amostra | legítimo |
| `make_votu_table` | dados por amostra | **correto** — já namespaceia (ver abaixo) |

**`prodigal_viral`**: entrada exclusivamente global, saída por amostra. Com 32
amostras são 32 execuções produzindo bytes idênticos. Deve virar regra global.

#### O bug de namespace de ID (analisado em 2026-08-18)

O `pharokka` **não estava só suspeito — está quebrado**, e não sozinho.

**Causa raiz.** O commit `8ef8bb4` (2026-08-13, *"remove a cadeia per-amostra e
reponta os 9 consumidores no catálogo global"*) trocou o input dos consumidores
de `viral_votu_reps` (por amostra) para `votu_catalog_reps` (global). Os dois
produzem **espaços de nome de ID diferentes**:

| fonte | ID |
|---|---|
| `viral_votu_reps` (antigo, por amostra) | `MEGAHIT_k141_10006` |
| `votu_catalog_reps` (global) | `S1\|MEGAHIT_k141_10006` |

O `build_pool` prefixa por origem de propósito. Mas os *joins por ID* dos
consumidores não foram movidos junto: continuam cruzando o catálogo global com
tabelas por amostra que têm ID nu. Nenhum deles falha — todos devolvem conjunto
vazio, que passa por resultado biológico.

**Três regras afetadas, verificadas uma a uma:**

1. **`pharokka`** (`annotation.smk:82-107`) monta `hq_set` do CheckV da amostra
   (IDs nus) e filtra o FASTA global (namespaced). Zero match → `hq_fa` vazio →
   cai no `if not hq_set or getsize == 0` → grava *"No HQ phages found —
   skipping"* e toca outputs vazios. Sai como *"essa amostra não tem fago de alta
   completude"*. Em cascata: `phold` recebe GBK vazio, `genome_map_phage`
   também.
2. **`genome_map_virus`** (`genome_map_universal.py:714-733`), independente do
   pharokka: `comp_map` vem do CheckV (nu), `seqs` do FASTA global (namespaced),
   então `comp_map.get(gid, 0.0)` sempre devolve 0.0 e nenhum genoma passa do
   `min_comp`. Zero mapas.
3. **`viral_taxonomy`** (`taxonomy.smk:376-434`) itera os contigs do FASTA global
   (namespaced) e busca `genomad_tax.get(contig)`, cujas chaves são o `seq_name`
   do geNomad **daquela amostra** (nu). O geNomad roda sobre o `rep_seq.fasta` da
   amostra, então nunca casa — **o geNomad sai do merge de taxonomia**. O MMseqs
   escapa porque roda sobre as proteínas do `prodigal_viral`, que já é global.

**Reprodução** (cópia literal da lógica do `pharokka`):

```
hq_set (do CheckV da amostra) : ['MEGAHIT_k141_10006', 'MEGAHIT_k141_777']
headers do catálogo global    : ['S1|MEGAHIT_k141_10006', 'S2|...', 'S1|...']
sequências extraídas          : 0
```

**Por que ninguém viu.** O `results/` da corrida atual **não tem
`votu_catalog/`** — os outputs de pharokka em disco são da arquitetura anterior
(headers `k141_2550`, sem prefixo). O código atual nunca rodou até o fim; a
regressão está latente desde 13/08.

**A convenção era conhecida e só não foi propagada:** `votu_catalog.smk:79` monta
a chave como `f"{source_id}|{contig}"`, e `make_votu_table.py:309` compara com
`f"{sample}|{mem}" == rep` enquanto busca o CheckV com o ID nu. Esses dois estão
certos e servem de molde.

**O conserto não é adicionar o prefixo no join.** Mesmo com o namespace certo, a
semântica continua errada: cada amostra filtraria o catálogo global pelo CheckV
dela, e um representante originário da S2 seria aceito ou descartado pela tabela
da S1. O conserto é o próprio (h) — as regras viram globais e usam a completude
do próprio catálogo (`_load_catalog_completeness()`, que já namespaceia certo). O
namespace se resolve como consequência.

**Isso reclassifica o (h): não é só otimização de custo, é correção de
resultado.**

As duas regras adicionadas em (f) e (g) já nasceram globais e não aumentam esta
dívida.

#### Fechamento (2026-08-18): `defensefinder_viral` e `dbapis_viral` também globalizadas

Auditoria adicional encontrou o mesmo padrão em duas regras que já liam a
saída global do `votu_prodigal` (movida na primeira metade do (h)) mas
continuavam sendo regras **por amostra**, com gêmeas por grupo em
`rules/coassembly.smk`: `defensefinder_viral` (`rules/defense_amr.smk:181`)
e `dbapis_viral` (`:294`).

Duas consequências, iguais em espécie às do `pharokka`/`viral_taxonomy`/
`genome_map_virus` documentadas acima:

1. **Desperdício N+G**: rodavam sobre entrada byte-idêntica em toda
   amostra/grupo, o mesmo problema diagnosticado para `prodigal_viral`, só
   que um nível abaixo.
2. **Mudança de significado silenciosa nos grupos**: como cada regra de
   grupo também lia o `.faa` global, `{group}/viral/defensefinder/
   viral_defense_systems.tsv` continha os sistemas de defesa de **todo o
   catálogo**, não do grupo — e `rules/coassembly.smk` copiava esse
   arquivo para `final/viral/defense_amr/` como se fosse do grupo. Mesma
   classe de bug do namespace de ID, mas na camada de finalize em vez de
   join.

Movidas para `rules/votu_catalog.smk` como `votu_defensefinder_viral` e
`votu_dbapis_viral`, saída sob `{OUTDIR}/votu_catalog/{defensefinder,dbapis}/`,
lógica interna preservada integralmente (mesmo `--antidefensefinder`, cache
de modelos, auto-download do dbAPIS, split defense/antidefense). As regras
por amostra e as gêmeas de grupo foram apagadas; os consumidores
(`Snakefile`, `rules/finalize.smk`, `rules/coassembly.smk`, `rules/report.smk`,
`scripts/report/data_loaders.py`, `scripts/report/renderer.py`,
`scripts/report/components/hostdefense.js`) foram repontados para o caminho
global. `rules/defense_amr.smk:rule defensefinder` (procariótico, por bin) e
suas gêmeas de grupo **não foram tocadas** — bins não são dedupados num
catálogo global, então essa regra continua legitimamente por amostra.

Verificado com dry-run nos quatro perfis (`config_dagtest`/PE,
`config_se`/SE, `config_lr`/LR-ONT, `config_lr_hifi`/LR-HiFi): delta exato
de −2 jobs por amostra, −2 por grupo, +2 globais em todos os quatro, sem
mudança de contagem em nenhuma outra regra (`defensefinder` procariótico
incluído).

---

## Achados que não devem ser redescobertos

### O "consenso" viral hoje é uma união

`viral_consensus_mode: "hybrid"` com `min_viral_tools: 2`, mas a lógica
(`viral_detection.smk:252-287`) faz `len(tools) >= 2 **OR** n in high_conf`, e:
- geNomad roda com `--min-score 0.7` contra um portão `score_genomad_min: 0.5`
  → **toda** chamada do geNomad passa
- VIBRANT entrava no `high_conf` se `type` fosse lytic/lysogenic → **toda**
  chamada passava

Ou seja o consenso é **≈ união** das ferramentas, e o `min_viral_tools` é
praticamente decorativo nesse modo. Isso explica o excesso de FP. Depois de (a),
revisar se `hybrid` ainda faz sentido ou se vira `score` puro.

**RESOLVIDO em 2026-08-19.** Duas conclusões, e a segunda é um bug.

1. **O `hybrid` foi mantido.** Ele não perdeu sentido por sobrarem duas
   ferramentas — perde sentido com `min_viral_tools: 1`, porque aí
   `len(t) >= 1 OR high_conf` é literalmente idêntico a `count`, e ambos são a
   união pura (todo contig em `tool_hits` tem ≥ 1 ferramenta). Com
   `min_viral_tools: 2` os três modos continuam distintos e ordenados:
   `count` = interseção, `score` = qualquer um acima do corte,
   `hybrid` = `count ∪ score`. São duas linhas de código que sustentam um
   meio-termo real ("as duas concordam, OU uma está muito confiante").

   **Mas com os cortes de hoje, trocar `score` por `hybrid`+2 é um no-op.** O
   geNomad roda com `--min-score 0.7` (`viral_detection.smk:87`) e o portão do
   consenso é `score_genomad_min: 0.7` — todo `virus_summary` já vem ≥ 0.7,
   então **toda** chamada do geNomad passa por construção e
   `geNomad ⊆ high_conf`. Daí:

       count(2)  = VS2 ∩ geNomad ⊆ geNomad ⊆ score
       hybrid(2) = count(2) ∪ score = score

   Idênticos. O VS2 é o oposto: roda a 0.5 contra um portão de 0.7, então ali o
   portão filtra de verdade. **A alavanca não é o modo, é `score_genomad_min`** —
   para o geNomad ser seletivo no consenso ele teria de subir acima de 0.7
   (0.8, espelhando o vão 0.5→0.7 do VS2). Mantido em 0.7 de propósito: vírus
   novos de ambiente pontuam baixo, e o portão composto do item (e) + CheckV já
   controlam FP a jusante.

2. **O modo `score` era geNomad sozinho** (commit `c87e13c`). O bloco
   `high_conf` guardava os nomes crus dos TSV de score, enquanto as chaves de
   `tool_hits` já vinham normalizadas. Como o `seqname` do VirSorter2 carrega
   `||full`/`||partial` — documentado em `results-reference/file_schemas.json`
   — o teste `n in high_conf` **nunca** casava para o VS2. Proviroses do
   geNomad (`contig|provirus_X_Y`) caíam pelo mesmo motivo. Consequência: o
   braço de score do `hybrid` nunca acrescentou um contig sequer, e o modo
   `score` filtrava só pelo geNomad. Mesma família dos bugs de namespace deste
   roadmap: duas normalizações escritas inline em um lugar só e não repetidas
   no outro. Viraram `_norm_vs2`/`_norm_genomad`, definidas uma vez.

   Nota para reexecuções: com a correção, `score` e `hybrid` passam a manter
   **mais** contigs do que mantinham antes (os de score alto só do VS2). Um
   resultado viral anterior nesses modos não é comparável ao novo.

### Descartado: strobealign como mapeador — por ora

Avaliado a fundo em 2026-08-18. Ver `docs/ANALISE_TOOLS_VOMIX_METAFUN.md` e as
referências abaixo.

- Nem o vOMIX nem o metaFun usam strobealign (o vOMIX usa **minimap2 + CoverM**).
- Existe paper de metagenômica — **AEMB** (Pan et al. 2025) — mas ele é sobre
  viabilizar binning **multi-against-multi** (N² alinhamentos), não sobre trocar
  o alinhador.
- A vantagem só aparece com **≥30 amostras**. Nos dados atuais: 32 amostras em 7
  grupos de co-assembly (~4-5 por grupo) — **nenhum dos dois modos da vapor é o
  modo em que o AEMB ganha**.
- No oceano (n=109) o AEMB deu **−11,3% bins HQ** vs Bowtie2(concat). Efeito
  dependente de ambiente.
- **Limitação estrutural:** Methods 4.2.5 do paper diz textualmente que o AEMB
  *"cannot give the per-base coverage"*. O MetaBAT2 da vapor depende do
  `jgi_summarize_bam_contig_depths`, que usa profundidade **e variância** por
  base. O AEMB não pode alimentar o MetaBAT2.
- O BAM continua necessário de todo jeito (vRhyme, CoverM, COBRA, MetaBAT2).

**Reabrir só se** a vapor for reestruturada para binning multi-against-multi.

### Descartado: `Binette(MetaBAT2 + SemiBin2)` está correto

O metaFun mediu: MetaBAT2 sozinho dá menos MAGs e **pior ARI**; o refinamento do
DAS Tool melhora acurácia de bp **mas reduz ARI** no marinho e na rizosfera (no
gut humano, o oposto). A lição **não** é "menos binners" — é que refinamento não
é de graça e o efeito depende do ambiente. A escolha atual da vapor está alinhada
com a evidência. **Não mexer.**

### Descartado: PPanGGOLiN

Exige vários MAGs da mesma espécie e responde pergunta de genômica comparativa.
O eixo do metaFun é invertido ao nosso (ele é comparativa com metagenômica de
apoio; a vapor é viroma com procarioto de apoio).

### Assimetria do vOMIX

Ele prega "uma ferramenta bem escolhida" no viral e usa **quatro binners**
(VAMB + MetaBAT2 + MaxBin2 + CONCOCT → DAS Tool) no procariótico, sem benchmark
que sustente. Não importar o lado procarioto deles por autoridade.

### As 6 ferramentas virais do vOMIX são só benchmark

Gating verificado no `Snakefile` deles: `viral-benchmark` não está no
`viral-end-to-end` (o módulo default). Padrão que vale copiar — mover VirSorter2
para um módulo opcional de benchmark em vez do caminho default.

---

## Pendências menores

- 28 `touch {output.done}` cosméticos restantes (regras que falham alto; o
  relatório vê "existe" em vez de "ok"). Baixa prioridade.
- Apagar as branches `backup/master-pre-rewrite` e
  `backup/unit-wildcard-pre-rewrite` — **só depois** de confirmar que o pull deu
  certo em todas as máquinas. São a rede de segurança da reescrita de autoria.
- Apagar a branch `refactor/unit-wildcard` (já mergeada).
- `config_amazon_12-08-26.yaml` (untracked) tem chaves órfãs `vcontact3_db` e
  `vcontact3_ver`; ganhará `vibrant_base` órfã depois de (a). Inofensivas — o
  Snakemake ignora chaves desconhecidas.
- Se algum dia entrar strobealign ou inStrain: a imagem BioContainers está bem
  atrás do bioconda (strobealign 0.9.0 vs 0.17.0), então pinar hoje faria o
  `scripts/check_env_container_sync.py` acusar divergência. Resolver com imagem
  custom no GHCR (a infra já existe: `ghcr.io/lymf/vapor-genome-map` com
  `custom: true`, mais `.github/workflows/build-containers.yml`).

---

## Ideias ainda não avaliadas

Nenhuma foi analisada a fundo — não tratar como aprovadas.

### Pendências abertas pelo (h) — TODAS CONSERTADAS em 2026-08-18

Quatro defeitos de **empacotamento/disponibilidade/observabilidade**, nenhum
produz resultado biológico errado — mas o (4) esconderia um que produzisse.
Dois são anteriores ao (h).

1. ~~**`final/viral/defense_amr/` virou diretório vazio**~~ — **FEITO**. Entrou a
   regra `finalize_votu_catalog` (molde: `finalize_reads_classify`), que copia
   os outputs globais para `final/votu_catalog/` e loga item a item o que foi
   copiado e o que estava ausente — para "sumiu" e "não rodou" não se
   confundirem. Os dois `mkdir` órfãos saíram. Descrição do problema original:
   O `finalize.smk:176` e o `coassembly.smk:1755` continuam criando o
   diretório, mas com o defensefinder/dbAPIS virais globalizados **nenhuma
   regra escreve dentro dele** — verificado por grep, não há um único writer.
   Conserto: uma regra `finalize_votu_catalog` copiando
   `votu_catalog/{defensefinder,dbapis}/*.tsv` (e, pela mesma lógica, pharokka,
   phold e os genome maps) para `final/votu_catalog/`, mais remover o `mkdir`
   órfão dos dois organize. Precedente pronto: `finalize_reads_classify`
   (`finalize.smk:295-307`). **É decisão de layout do `final/`** — por isso não
   foi feito junto.

2. ~~**Corrida entre o relatório e a anotação**~~ — **FEITO**. O `report.smk`
   ganhou arestas para `votu_pharokka`, `votu_phold` e os dois
   `votu_genome_map_*`, declaradas só como restrição de ordem (os loaders
   seguem resolvendo os próprios caminhos). Descrição do problema original:
   `rules/report.smk` não tem aresta `input:` para `votu_pharokka`,
   `votu_phold` nem `votu_genome_map_*`, e o `generate_report` não os alcança
   transitivamente. Com `--cores` alto, o `report.html` pode ser escrito antes
   do pharokka existir: gráfico PHROGS vazio e nenhum genome map, **sem erro no
   log**. O `load_phrogs`/`load_genome_maps` hardcodam o caminho global sem
   aresta. Conserto: arestas gated por `VOTU_CATALOG_ENABLED`.

3. ~~**`votu_catalog_enabled: false` ficou auto-contraditório**~~ — **FEITO, pela
   raiz: a flag foi removida.** Ela nasceu em 13/08 quando o catálogo global
   era opcional e o clustering por amostra ainda existia como alternativa;
   esse clustering foi removido, então a flag guardava uma alternativa
   inexistente. E nunca desligou nada de verdade — gateava só os alvos do
   `rule all`, nunca as definições de regra, que os demais consumidores
   puxavam de volta. O catálogo passa a ser estágio obrigatório. Descrição do
   problema original: Os alvos de
   pharokka/phold/genome maps foram para dentro do bloco
   `if VOTU_CATALOG_ENABLED:` (`Snakefile:444`), quando os antigos alvos por
   amostra eram incondicionais — com a flag desligada o relatório perde PHROGS
   e todos os genome maps virais. Ao mesmo tempo os dois inputs novos do
   `report.smk` são incondicionais e arrastam a cadeia do catálogo de volta.
   A flag já era meio decorativa; o (h) alargou a inconsistência.

4. ~~**Nenhuma das regras globais novas era rastreada no relatório**~~ — **FEITO**
   (2026-08-19). As dez regras que o (h) moveu para o catálogo global ficaram
   fora de `STATUS_TRACKED_GLOBAL_TOOLS` (`scripts/report/data_loaders.py`):
   pharokka, phold, os dois `genome_map`, defensefinder e dbAPIS virais, mais
   prodigal e os três de taxonomia. Se qualquer uma quebrasse, o relatório não
   mostrava **uma linha sequer** — a aba apenas ficava vazia, que é exatamente
   a aparência de "não havia nada a mostrar".

   O conserto tem duas metades, e só a segunda é óbvia:

   - **Seis regras nem escreviam status de verdade** — `votu_pharokka`,
     `votu_genome_map_virus`, `votu_defensefinder_viral`,
     `votu_mmseqs_taxonomy` e seu gêmeo `_custom` faziam
     `Path(done).touch()`/`touch {output.done}`. Registrá-las no dicionário
     sem antes corrigir isso teria produzido o erro **oposto**: `done.txt`
     vazio é lido como `unknown`, e `tool_failed()` trata `unknown` como
     lacuna — toda execução bem-sucedida passaria a aparecer como falha.
     Todas passaram a usar `write_status()`.
   - Duas engoliam a falha ativamente: `votu_defensefinder_viral` tinha
     `|| echo WARNING` e `votu_genome_map_virus` tinha `|| true`, ambos
     seguidos de `touch` cego. Um DefenseFinder quebrado escrevia tabelas
     vazias e era lido como "nenhum sistema anti-defesa neste viroma" — a
     mesma classe de erro do AMRFinderPlus com disco cheio que originou a
     convenção do `write_status`. Agora capturam a exceção e gravam
     `failed: <erro>`, sem derrubar o DAG (o gêmeo `votu_genome_map_phage` já
     fazia isso via `|| RC=$?`).

   Guarda automática: `test_every_global_done_file_is_status_tracked` extrai
   por regex todo `f"{CATALOG_DIR}/…done….txt"` do `.smk` e exige que esteja
   no dicionário. Hoje são 14 de 14. A próxima regra migrada não consegue
   repetir o esquecimento sem quebrar a suíte.

### Sétima manifestação — os bins do vRhyme nunca foram encontrados

**Verificada contra os dados da Amazônia em 2026-08-18** (7 grupos de
co-assembly, 2870 headers de bin). Anterior a todo este roadmap.

O layout real do vRhyme não é o que o código assumia:

| | real | o que o código procurava |
|---|---|---|
| diretório | `vrhyme/vRhyme_best_bins_fasta/` | `vrhyme/` |
| arquivo | `vRhyme_bin_100.fasta` | `vRhyme_best_bins.*.fasta` |
| header | `>vRhyme_100__k141_246312` | `k141_246312` |

O glob casava **zero** arquivos, em 100% dos casos. O padrão parecia plausível
numa leitura rápida porque `vRhyme_best_bins.19.*` **existe** naquele
diretório — só que como `.membership.tsv` e `.summary.tsv`, nunca `.fasta`.

Cinco pontos usavam o glob errado: `viral_nonredundant` (loop e contador),
o par de grupo em `coassembly.smk`, o `split_viral_fastas.py` (PHIST) e o
`finalize.smk`. O `checkv_vrhyme` e o `data_loaders.py` acertavam — o que
tornava a divergência ainda menos visível.

**Consequências, todas silenciosas:**
- nenhum vMAG entrava no conjunto não-redundante;
- o contador de bins do log reportava 0, então "não há bins" e "não achei os
  bins" eram indistinguíveis;
- o dedup do `split_viral_fastas.py` nunca disparava (a correção de namespace
  do commit `a317707` estava certa mas era insuficiente — o glob não achava
  nada de qualquer jeito);
- e, depois do item (e), o braço "está num bin" do portão nunca disparava.

**Impacto medido no grupo TAP** (1429 sequências), simulando o portão do (e)
com dado real:

| | mantidos | descartados |
|---|---|---|
| com o glob quebrado | 378 | **1051** |
| com o glob corrigido | 942 | 487 |

Ou seja: 860 contigs binados seriam descartados por falta de evidência que
existia e não estava sendo lida. O item (e) não causou isso — mas teria sido
quem transformaria o bug latente em perda de dado.

**Conserto:** `scripts/vrhyme_bins.py`, com `bin_fastas()`,
`contig_from_bin_header()` e `binned_contigs()`. Um lugar só, para não haver
uma sexta cópia do glob. Verificado: 860 de 860 contigs binados do TAP passam
a casar (era 0).

### Hospedeiro no track de reads: duas fontes — FEITO (parte) 2026-08-19

Fechado o item que ficara aberto na caça: o `sylph-tax merge` roda com
`--column relative_abundance`, então a coluna `Virus_host (if viral)` **nunca
chegava** ao `collapse_by_host`. Entrou `rule reads_host_map`
(`scripts/reads_classify/build_host_map.py`), que lê a fonte original — os
`.sylphmpa` por amostra.

Medido nos 32 arquivos da Amazônia: 12.716 linhas, **66 com hospedeiro
não-nulo, 24 clados únicos**. Depois de resolver por gênero: **8 táxons virais
com hospedeiro real** (`Acinetobacter` 5, `UBA3064` 3) de 1407. Ou seja, a
anotação do banco cobre ~0,6% — é exatamente por isso que o PHIST é necessário
para o resto, e por isso as duas fontes ficam em colunas separadas em vez de
uma só.

O `viral_abundance_by_host.tsv` passa a ter `host_source` (`db`, `phist`,
`db,phist`, `none`) e ganha o sidecar `viral_host_assignments.tsv`, com
`clade_name | host_db | host_phist | host_genus | host_source` — a proveniência
de cada número agregado, mesma filosofia do sidecar de descarte do item (e).

A agregação é por **gênero só**, não por `(gênero, fonte)`: agrupar pelos dois
partiria o mesmo gênero em duas linhas quando parte dos seus vírus vem do banco
e parte do PHIST, e no gráfico isso apareceria como dois táxons homônimos.

Precedência: **banco primeiro, PHIST depois**. O banco é curado e descreve o
genoma de referência que o sylph de fato detectou; o PHIST entra onde o banco
cala. O contrário deixaria uma predição por k-mer sobrescrever uma atribuição
publicada.

**Ainda falta: PHIST no track de reads.** `_RC_PHIST_MAP` /
`reads_classify_phist_map` já existe e alimenta a coluna; nenhuma regra o
produz. Duas coisas travam, e a segunda é uma armadilha:

1. **Formato.** O PHIST aqui é `kmer-db build` + `kmer-db new2all` + `phist`.
   Exige **um FASTA por genoma** dos dois lados (o kmer-db dá uma linha de
   resultado por ARQUIVO, não por sequência), mais dois arquivos de lista com
   um caminho por linha. É por isso que existe o `split_viral_fastas.py`.
2. **Os resultados de PHIST que já existem não servem.** O `rule phist` roda
   sobre contigs MONTADOS (`k141_...`) contra os MAGs da amostra. Os clados do
   sylph são genomas de REFERÊNCIA do IMG/VR (`t__IMGVR_UViG_...`) — entidades
   diferentes, sem chave em comum. Reaproveitar a tabela existente produziria
   junção vazia (mais um caso da família). É preciso um PHIST novo, sobre os
   genomas detectados pelo sylph, extraídos do `reads_classify_genome_fasta`
   pelo `Contig_name` — a mesma extração que o `bacphlip_lifestyle.py` faz, e
   cujo bug de ID foi corrigido hoje.

Decisão pendente do usuário: quem são os hospedeiros candidatos — os MAGs da
própria amostragem (liga o perfil por reads às bactérias montadas) ou os
genomas GTDB que o sylph detectou (independe de montagem)?

### Caça a bugs de tratamento de dados — 2026-08-19

Varredura pedida depois da oitava manifestação, com o mesmo método: **verificar
cada suposição de formato contra os arquivos reais**, nunca contra o comentário
do código. Dados usados: `~/global/results` (run completa, com `report.html`) e
a run da Amazônia em `/media/nas1/.../amazon/results`.

**Cinco defeitos confirmados, todos medidos.**

1. **`build_host_collapse` inventava o gênero do hospedeiro.** Derivava de
   `host.split('_')[0]`, apoiado no comentário "primeiro token do nome do
   arquivo (ex. `Bacteroides_fragilis.fa` → `Bacteroides`)". Isso vale para
   genomas de *referência* nomeados por espécie, não para os MAGs desta
   pipeline. Prova no `report.html` real: todo `Host` é `binette_binN`, logo o
   "gênero" de **todos** os hospedeiros era a constante `"binette"` — um único
   feixe no gráfico. No co-assembly os bins são números e sairiam como gêneros
   `1`, `136`. O dado certo estava ao lado, no mesmo relatório: o mesmo
   hospedeiro aparece em `HOST_DEFENSE_LINKS` com `Host_genus: "Neisseria"`,
   resolvido pelo GTDB-Tk. Passa a usar essa fonte.

2. **`load_phist` não removia `.fna`.** Os bins do VAMB (o track co-assembly usa
   `bin_ext=".fna"`) chamam-se `1.fna`, então `Host` ficava `1.fna` e a junção
   com o `user_genome` do GTDB-Tk (`1`) falhava **para o track inteiro** —
   taxonomia de hospedeiro vazia em todo grupo, sem erro. O track por amostra
   escapou por acaso: bins do Binette são `.fa`. De quebra, `.fa` era removido
   antes de `.fasta`, então `x.fasta` virava `xsta`.

3. **BACPHLIP nunca funcionaria com banco IMG/VR.** Os IDs saíam de
   `basename(Genome_file).split(".")[0]`, o que serve para o GTDB (um arquivo
   por genoma) mas não para o IMG/VR, onde **todo vírus vem do mesmo arquivo**:
   nas 111 linhas de sylph da Amazônia, `Genome_file` é `imgvr_reps.fna` em
   todas as virais. O conjunto detectado colapsava para um id inexistente,
   `imgvr_reps`, zero sequências eram extraídas e o script abortava. Quem
   distingue os vírus é `Contig_name`. Dormente na config atual
   (`reads_classify_genome_fasta` vazio), mas é mina para quem ligar.

4. **O colapso viral por hospedeiro do sylph produzia tabela vazia.** O filtro
   exigia `d__Viruses` **e** `s__`. A taxonomia do IMG/VR no sylph-tax começa no
   **realm** (`r__Duplodnaviria`) e as linhagens **pulam espécie** — a folha é
   `t__IMGVR_UViG_...`. Ou seja, as duas condições falhavam. Medido: 1482 clados
   virais na tabela mesclada, **zero** casando, e o `viral_abundance_by_host.tsv`
   em disco tem exatamente **uma linha, o cabeçalho**. Depois da correção, 1407
   folhas virais. A seleção passou a ser "folha da hierarquia" em vez de um rank
   fixo, o que funciona nas duas taxonomias sem saber qual banco gerou a tabela.

5. **`enrich_taxonomy_with_checkv` perderia todo profago** — consequência da
   correção da oitava manifestação: com o aparo do CheckV funcionando, o
   `Genome` da taxonomia pode ser `k141_9_1` enquanto o `quality_summary` só
   conhece `k141_9`. Resolvido com o mesmo `checkv_provirus.resolve_original_id`.

**Verificados e CORRETOS** (registrado para não se re-investigar):

- `prok_binning.smk:106` usa `seq.split("|")[0]` para provírus do geNomad — e o
  geNomad **de fato** usa `contig|provirus_START_END`. Confirmado no único caso
  do dataset: `k141_251078|provirus_255_24262`. CheckV e geNomad têm convenções
  diferentes; foi aplicar a do geNomad ao arquivo do CheckV que gerou a oitava
  manifestação.
- `votu_abundance` (`abundance.smk`) despreza o prefixo do catálogo
  explicitamente antes de casar com a tabela do CoverM. Correto.
- As 4 categorias de `_PHROGS_HALLMARK_CATEGORIES` existem literalmente na saída
  real do pharokka.
- `rgi_card`, `galah_derep` e `gtdbtk` escrevem status de verdade — os
  `done.txt` vazios em disco são de uma run anterior à convenção.
- `viral_depth` casa `contigName` do depth com `contig_id` do CheckV; ambos nus.

**Aberto, precisa de decisão de layout:** o `sylph-tax merge` não propaga a
coluna `Virus_host (if viral)` dos `.sylphmpa` por amostra para a tabela
mesclada, então `has_host` é sempre False e o agrupamento por hospedeiro nunca
roda de fato — escreve-se a saída em nível de taxon. Os `.sylphmpa` por amostra
**têm** a coluna (107 linhas com valor na Amazônia), mas todos os valores são
`UNKNOWN;UNKNOWN;...` — o IMG/VR não atribuiu hospedeiro a esses vírus. Levar o
host até o `collapse_by_host` exige mudar o merge; com o IMG/VR sozinho o ganho
seria nulo, com o UHGV não.

### Oitava manifestação — o aparo de provírus do CheckV nunca acontecia

**Medida contra os dados da Amazônia em 2026-08-19** (7 grupos de co-assembly).
Anterior a todo este roadmap. Consertada nos commits desta data.

O código lia `proviruses.fna` e derivava o contig de origem com
`hdr_id.rsplit('|', 1)[0]`, apoiado num comentário que afirmava que o header é
`"orig_id|start_end"`. **Não é.** O formato real:

    >k141_219139_1 1-13933/18998
    >k141_97527_1 1-2446/3389

É `{contig}_{n}`, e o intervalo está na *descrição*, depois do espaço — que todo
leitor descarta com `.split()[0]`. Sem nenhum `|` no header, o `rsplit` era um
no-op: devolvia o header inteiro, a chave nunca batia com o contig do consenso,
e o efeito foi duplo — **o fragmento aparado era descartado E a sequência
original, com o DNA de hospedeiro flanqueando o profago, era emitida como
"viral"**.

| | |
|---|---|
| provírus nos 7 grupos | 118 |
| resolvidos antes | **0** |
| resolvidos depois | 118 |
| pb de hospedeiro carregados como "viral" | **216.458** |

Consequências a jusante: comprimento inflado no portão do item (e), ANI errada
no clustering de vOTU, e genes de hospedeiro anotados como virais por
pharokka/phold. Este último importa em particular para o eixo
defesa/anti-defesa: genes de defesa **do hospedeiro** dentro do flanco podiam
ser contados como anti-defesa **do fago**.

Três cópias da mesma linha errada (`viral_binning.smk`, `coassembly.smk`, e o
`lookup_key` do portão de grupo) — foi ter três cópias que permitiu a
divergência sobreviver. Agora há um módulo só, `scripts/checkv_provirus.py`,
com 8 testes. Ele **não adivinha delimitador**: testa candidatos contra o
conjunto de ids que se sabe existirem, e devolve um flag `resolvido` que o
chamador conta. Se o CheckV mudar o formato de novo, as regras abortam com
mensagem explícita em vez de emitir flanco de hospedeiro em silêncio.

**Uma armadilha no conserto, que quase virou um bug novo:** como o aparo nunca
acontecia, o catálogo global **nunca via headers com sufixo**. Consertar só
essa camada faria os ids `k141_219139_1` chegarem ao pool, onde
`_load_catalog_completeness()` os chavearia contra um `quality_summary.tsv` que
só conhece `k141_219139` — todo profago cairia para completude 0.0, sumindo do
nível `mq` e do pharokka. Ou seja, o conserto teria movido o descasamento de
lugar em vez de removê-lo. Por isso `_load_catalog_completeness()` ganhou uma
segunda passada sobre o `provenance.tsv`, que copia a evidência do contig de
origem para o id do fragmento.

### Quinta manifestação do bug de namespace — `split_viral_fastas.py` — CONSERTADA

**Anterior ao (h)** (vem do `8ef8bb4`), e é o mesmo padrão. O roadmap marcou o
`phist` como "legítimo" e não pegou.

O `phist` passa ao script o `mq_fasta` **global** (headers namespaced) junto com
os bins de `{sample}/bins/vrhyme/`, cujos headers são do conjunto viral **por
amostra** (nus). O script monta `binned_contigs` a partir dos bins (nus) e testa
`if seq_name not in binned_contigs` contra nomes namespaced — **o teste nunca
dispara**. Resultado: todo contig já binado é escrito de novo como FASTA
individual, e o PHIST emite linha duplicada para a mesma sequência — exatamente
a redundância que o docstring do script diz eliminar. A jusante, o
`build_host_collapse` conta `n_viruses` em dobro e soma `total_rpkm` duas vezes
por gênero.

Efeito colateral relacionado: as linhas de bin do PHIST trazem
`Virus = vRhyme_best_bins.3`, que não casa com nada nos joins de antidefesa e
abundância — essas linhas mostram `—` na matriz para sempre. Inerente a tratar
bin como genoma, não é o mesmo bug.

**Conserto (2026-08-18):** o script passou a receber o `source_id` como
argumento e a registrar no `binned_contigs` **as duas formas** do nome de cada
contig binado — a nua e a namespaced. Foi essa a direção escolhida, e não
tirar o prefixo dos nomes do catálogo: ID nu só é único **dentro** de uma
fonte, então strippar deixaria um contig de outra amostra colidir com os bins
desta. Ganhou também um WARNING explícito para o caso "bins carregados mas
nenhum contig casou como já binado", que é a assinatura exata do
descasamento.

- **`eggnog_viral` roda um Prodigal próprio** sobre o mesmo `mq_fasta`
  (`votu_catalog.smk`), duplicando o `votu_prodigal` que o (h) tornou global.
  Achado pela auditoria do (h) em 2026-08-18. Pequeno e no mesmo espírito do
  (h): trocar pelo `.faa` global e conferir se o eggNOG precisa de algum flag
  de predição que o `votu_prodigal` não usa (prodigal-gv vs prodigal padrão —
  **verificar antes de assumir que são intercambiáveis**).

- **Catálogo global de MAGs procarióticos** — o análogo do (h) do lado
  procariótico, e o maior item em aberto. Verificado em 2026-08-18: o
  `galah_derep` é **por amostra** (desreplica os bins da própria amostra) e
  alimenta **só o GTDB-Tk** (`prok_binning.smk:711`). Não existe nada
  equivalente ao `votu_catalog`. Consequência: MAGs da mesma espécie
  recuperados em amostras diferentes são anotados repetidamente por bakta,
  eggNOG e pelas ferramentas de AMR, e nada garante que dois deles recebam a
  mesma anotação.

  A biologia que justifica o (h) vale igual — 95% ANI é nível de espécie para
  procarioto também. Mas o custo é outro: no viral o catálogo global já
  existia e a metade difícil estava pronta; aqui seria preciso **construir** um
  galah cruzando amostras, com provenance e herança próprios.

  **Diferença qualitativa que importa para priorizar: aqui nada está saindo
  errado.** É redundância de compute e risco de inconsistência entre anotações
  do mesmo organismo — não resultado zerado passando por resultado biológico,
  que era o caso do `pharokka`. Não tem a urgência que o (h) tinha.

  Não confundir com as regras de AMR e anotação procariótica em si
  (`amrfinderplus`, `rgi_card`, `deeparg`, `abricate`, `defensefinder`,
  `bakta`, `eggnog_prok`): elas consomem `prok_bin_proteins`, que é derivado
  dos bins **daquela** amostra — conteúdo genuinamente diferente por amostra.
  Rodar N vezes ali não é desperdício, é o trabalho. E não há ali a classe de
  bug do (h): o consumo é por amostra e o dado é por amostra, sem espaço para
  ID nu bater contra namespaced.

- **Retirar os genome maps** (`genome_map_phage`, `genome_map_virus`,
  `genome_map_prok`, `scripts/genome_map*.py`, painel do relatório) — intenção
  declarada pelo usuário em 2026-08-18. Enquanto isso não acontece, fica uma
  dívida conhecida e aceita: com os mapas globalizados em (h), o
  `load_genome_maps` devolve a mesma lista sob toda chave de amostra e o
  `renderer.py` serializa com `json.dumps`, que não deduplica — os SVGs inline
  passam a ser replicados N vezes no `report.html`. Com 32 amostras isso são
  dezenas de MB de HTML. **Não consertar**: some junto com os mapas.

- **dbCAN** (CAZymes) — lacuna real, o metaFun tem e a vapor não. Disponível:
  bioconda 5.2.9, quay.io `5.2.9--pyhdfd78af_0`.
- **hostile** para remoção de hospedeiro (bioconda 2.0.2) — substituiria
  bwa-mem2 + filtro samtools.
- **PhaBOX2** (PhaTYP, CHERRY) — bioconda 2.1.13.
- **CheckV-PyHMMER** — só compensa acima de ~10⁵ contigs.
- **inStrain** — nível de cepa.

---

## Como verificar cada passo

Invariante do DAG: o conjunto ordenado de arquivos de output do dry-run.

```bash
conda activate snakemake
SP=/tmp/.../scratchpad   # helpers: snap.sh, ruleedit.py, config_dagtest.yaml
snakemake -n -p --configfile $SP/config_dagtest.yaml > $SP/dag_X.txt
bash $SP/snap.sh $SP/dag_X.txt $SP/outputs_X.txt
diff $SP/outputs_ANTERIOR.txt $SP/outputs_X.txt
```

Baseline em 2026-08-18 (pós status, pré-VIBRANT): **1616 jobs, 2625 outputs**.

**Regra de ouro:** mudança que só mexe em corpo de `shell:` deve dar diff
**vazio**. Mudança que remove ferramenta deve dar um delta **explicável item a
item** — nunca aceitar "mudou mais ou menos o esperado".

### Lição registrada sobre método

Já houve um erro caro neste projeto: uma comparação automatizada que colapsava
diffs de corpo executável **não viu** um `*.fa` vs `*.fna` no `gunc`, o que teria
feito a trilha de grupo achar zero bins **em silêncio**. Ver §0 de
`docs/AUDITORIA_COASSEMBLY_PARES.md`.

> Contagem de diff normalizado serve para **priorizar**, nunca para **decidir**.

Vale igual para trabalho delegado a subagente: a spec pode ser mecânica, a
revisão não.

---

## Referências

Pan, S., Tolstoganov, I., Sahlin, K., Martin, M., Zhao, X.-M., & Coelho, L. P. (2025). AEMB: efficient abundance estimation for metagenomic binning. *bioRxiv*. https://doi.org/10.1101/2025.07.30.667338

Sahlin, K. (2022). Strobealign: flexible seed size enables ultra-fast and accurate read alignment. *Genome Biology*, 23, 260. https://doi.org/10.1186/s13059-022-02831-7

Aroney, S. T. N., Newell, R. J. P., Nissen, J. N., Camargo, A. P., Tyson, G. W., & Woodcroft, B. J. (2025). CoverM: read alignment statistics for metagenomics. *Bioinformatics*, 41(4), btaf147. https://doi.org/10.1093/bioinformatics/btaf147

Pratama, A. A., Bolduc, B., Zayed, A. A., Zhong, Z.-P., Guo, J., Vik, D. R., Gazitúa, M. C., Wainaina, J. M., Roux, S., & Sullivan, M. B. (2021). Expanding standards in viromics: in silico evaluation of dsDNA viral genome identification, classification, and auxiliary metabolic gene curation. *PeerJ*, 9, e11447. https://doi.org/10.7717/peerj.11447

Shaffer, M., Borton, M. A., McGivern, B. B., Zayed, A. A., La Rosa, S. L., Solden, L. M., Liu, P., Narrowe, A. B., Rodríguez-Ramos, J., Bolduc, B., Gazitúa, M. C., Daly, R. A., Smith, G. J., Vik, D. R., Pope, P. B., Sullivan, M. B., Roux, S., & Wrighton, K. C. (2020). DRAM for distilling microbial metabolism to automate the curation of microbiome function. *Nucleic Acids Research*, 48(16), 8883–8900. https://doi.org/10.1093/nar/gkaa621

Hockenberry, A. J., & Wilke, C. O. (2021). BACPHLIP: predicting bacteriophage lifestyle from conserved protein domains. *PeerJ*, 9, e11396. https://doi.org/10.7717/peerj.11396

### Documentos irmãos neste repo

- `docs/ANALISE_TOOLS_VOMIX_METAFUN.md` — análise ferramenta a ferramenta das duas pipelines de referência
- `docs/BENCHMARK_VOMIX_METAFUN.md` — arquitetura e performance
- `docs/VAPOR_TOOLS_MAP.md` — inventário por estágio das ferramentas da vapor
- `docs/AUDITORIA_COASSEMBLY_PARES.md` — laudo par a par da unificação amostra/grupo
