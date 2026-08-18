#!/usr/bin/env python3
"""Generate an editable SVG of the VAPOR pipeline (swimlane), scientifically corrected."""
import html

W, H = 2040, 1500
EL = []  # svg fragments

# palette
C = {
    'ink':'#0f172a','sub':'#475569','line':'#e6ebf0','lane':'#f8fafc',
    'teal':'#0d9488','tealbg':'#ccfbf1','tealink':'#134e4a',
    'green':'#16a34a','greenbg':'#dcfce7','greenink':'#14532d',
    'purple':'#7c3aed','purplebg':'#ede9fe','purpleink':'#4c1d95',
    'amber':'#d97706','amberbg':'#fef3c7','amberink':'#78350f',
    'indigo':'#4f46e5','indigobg':'#eef0ff','indigoink':'#312e81',
    'slate':'#64748b','iobg':'#f1f5f9','arr':'#94a3b8',
}
CLS = {
    'io':(C['iobg'],'#94a3b8',C['ink']),
    'viral':(C['tealbg'],C['teal'],C['tealink']),
    'prok':(C['greenbg'],C['green'],C['greenink']),
    'reads':(C['purplebg'],C['purple'],C['purpleink']),
    'coas':('#ffffff',C['indigo'],C['indigoink']),
    'hub':(C['amberbg'],C['amber'],C['amberink']),
    'plain':('#ffffff','#cbd5e1',C['ink']),
}

def esc(s): return html.escape(str(s), quote=True)

def box(x,y,w,h,lines,cls='io',dashed=False,bold0=True,fs=13):
    fill,stroke,tc = CLS[cls]
    sw = 2 if cls=='hub' else 1.4
    da = ' stroke-dasharray="5 4"' if dashed else ''
    EL.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="9" ry="9" '
              f'fill="{fill}" stroke="{stroke}" stroke-width="{sw}"{da}/>')
    if isinstance(lines,str): lines=[lines]
    n=len(lines); lh=15
    ty = y + h/2 - (n-1)*lh/2 + 4
    for i,ln in enumerate(lines):
        weight = 'bold' if (i==0 and bold0) else 'normal'
        size = fs if i==0 else fs-1.5
        col = tc if i==0 else C['sub'] if cls in ('io','plain') else tc
        EL.append(f'<text x="{x+w/2}" y="{ty+i*lh}" text-anchor="middle" '
                  f'font-size="{size}" font-weight="{weight}" fill="{col}" '
                  f'font-family="Helvetica, Arial, sans-serif">{esc(ln)}</text>')
    return {'x':x,'y':y,'w':w,'h':h,'cx':x+w/2,'cy':y+h/2,'r':x+w,'b':y+h}

def diamond(x,y,w,h,lines,stroke=None):
    stroke = stroke or C['slate']
    cx,cy=x+w/2,y+h/2
    EL.append(f'<polygon points="{cx},{y} {x+w},{cy} {cx},{y+h} {x},{cy}" '
              f'fill="#ffffff" stroke="{stroke}" stroke-width="1.6"/>')
    if isinstance(lines,str): lines=[lines]
    n=len(lines); ty=cy-(n-1)*7+4
    for i,ln in enumerate(lines):
        EL.append(f'<text x="{cx}" y="{ty+i*14}" text-anchor="middle" font-size="12" '
                  f'font-weight="600" fill="{C["ink"]}" font-family="Helvetica, Arial, sans-serif">{esc(ln)}</text>')
    return {'x':x,'y':y,'w':w,'h':h,'cx':cx,'cy':cy,'r':x+w,'b':y+h}

def arrow(x1,y1,x2,y2,color=None,dashed=False):
    color = color or C['arr']
    da = ' stroke-dasharray="5 4"' if dashed else ''
    EL.append(f'<path d="M {x1} {y1} L {x2} {y2}" fill="none" stroke="{color}" '
              f'stroke-width="1.6" marker-end="url(#ah)"{da}/>')

def elbow(x1,y1,x2,y2,color=None,dashed=False):
    """orthogonal: horizontal then vertical"""
    color=color or C['arr']
    da=' stroke-dasharray="5 4"' if dashed else ''
    midx=x2
    EL.append(f'<path d="M {x1} {y1} L {midx} {y1} L {x2} {y2}" fill="none" '
              f'stroke="{color}" stroke-width="1.6" marker-end="url(#ah)"{da}/>')

def a_rl(a,b,color=None,dashed=False): arrow(a['r'],a['cy'],b['x'],b['cy'],color,dashed)  # right->left
def a_tb(a,b,color=None,dashed=False): arrow(a['cx'],a['b'],b['cx'],b['y'],color,dashed)  # top->bottom
def lanelabel(y,tn,td,h):
    EL.append(f'<text x="18" y="{y+22}" font-size="12" font-weight="700" fill="{C["teal"]}" '
              f'font-family="Helvetica, Arial, sans-serif">{esc(tn)}</text>')
    yy=y+38
    for ln in td:
        EL.append(f'<text x="18" y="{yy}" font-size="10.5" font-weight="600" fill="{C["sub"]}" '
                  f'font-family="Helvetica, Arial, sans-serif">{esc(ln)}</text>')
        yy+=13
    EL.append(f'<line x1="8" y1="{y}" x2="{W-455}" y2="{y}" stroke="{C["line"]}" stroke-width="1"/>')

LX=150   # left of main flow
COASX=1600  # co-assembly column left

# ---- READS TRACK (separate top lane) ----
ry=70
EL.append(f'<rect x="{LX-4}" y="{ry-6}" width="{1440}" height="118" rx="12" fill="{C["purplebg"]}" '
          f'fill-opacity="0.35" stroke="{C["purple"]}" stroke-width="1.4" stroke-dasharray="6 4"/>')
EL.append(f'<text x="18" y="{ry+16}" font-size="12" font-weight="700" fill="{C["purple"]}" '
          f'font-family="Helvetica, Arial, sans-serif">READS TRACK</text>')
for i,ln in enumerate(['assembly-free','(independent of','the assembly pipe)']):
    EL.append(f'<text x="18" y="{ry+34+i*13}" font-size="10.5" font-weight="600" fill="{C["sub"]}" '
              f'font-family="Helvetica, Arial, sans-serif">{esc(ln)}</text>')
r1=box(LX+8, ry+18, 120,54,['raw reads'],'reads')
r2=box(LX+168, ry+18, 150,54,['sylph','sketch → profile'],'reads')
r3=box(LX+358, ry+18, 210,54,['sylph databases','IMG/VR · UHGV · GTDB · custom'],'reads')
r4=box(LX+608, ry+18, 130,54,['sylph-tax','taxonomy'],'reads')
r5=box(LX+778, ry+18, 300,54,['merged outputs','relative/sequence abundance · OTU','viral-abundance-by-host · BACPHLIP'],'reads')
for a,b in [(r1,r2),(r2,r3),(r3,r4),(r4,r5)]: a_rl(a,b,C['purple'])
rrep=box(LX+1118,ry+18,120,54,['→ report'],'reads')
a_rl(r5,rrep,C['purple'])

# ================= MAIN ASSEMBLY-BASED PIPELINE =================
EL.append(f'<text x="18" y="{ry+150}" font-size="11" font-weight="700" fill="{C["sub"]}" '
          f'font-family="Helvetica, Arial, sans-serif">ASSEMBLY-BASED PIPELINE</text>')

# ---- TIER 0 ----
y0=ry+165
lanelabel(y0,'TIER 0',['Input &','track select'],120)
b_in=box(LX+8,y0+22,150,64,['Sequencing input','Illumina PE/SE','ONT / PacBio HiFi'],'io')
d_rt=diamond(LX+185,y0+26,120,58,['Read','type?'])
b_sr=box(LX+335,y0+16,120,32,['Short reads'],'plain',fs=12)
b_lr=box(LX+335,y0+62,120,32,['Long reads'],'plain',fs=12)
d_tr=diamond(LX+485,y0+26,118,58,['tracks','(config)'],stroke=C['teal'])
b_tkV=box(LX+632,y0+18,150,30,['VIRAL — assembly'],'viral',fs=11)
b_tkP=box(LX+632,y0+54,150,30,['PROK — assembly'],'prok',fs=11)
a_rl(b_in,d_rt); a_rl(d_rt,b_sr); a_rl(d_rt,b_lr)
arrow(b_sr['r'],b_sr['cy'],d_tr['x'],d_tr['cy']-8); arrow(b_lr['r'],b_lr['cy'],d_tr['x'],d_tr['cy']+8)
a_rl(d_tr,b_tkV); a_rl(d_tr,b_tkP)
# reads selection routes UP to the separate reads lane
arrow(d_tr['cx'],d_tr['y'],d_tr['cx'],ry+118,C['purple'],dashed=True)
EL.append(f'<text x="{d_tr["cx"]+6}" y="{ry+140}" font-size="10" fill="{C["purple"]}" '
          f'font-family="Helvetica, Arial, sans-serif">reads track (separate)</text>')

# ---- TIER 1 ----
y1=y0+110
lanelabel(y1,'TIER 1',['Quality','control'],110)
b_fp=box(LX+8,y1+18,160,50,['fastp (SR)','Q≥20 · length · adapters'],'io')
b_hr1=box(LX+200,y1+18,150,50,['host removal','bwa-mem2 (opt)'],'io',dashed=True)
a_rl(b_fp,b_hr1)
b_np=box(LX+470,y1+18,110,50,['NanoPlot'],'io')
b_pc=box(LX+600,y1+18,130,50,['Porechop_ABI'],'io')
b_ft=box(LX+750,y1+18,110,50,['Filtlong'],'io')
b_hr2=box(LX+880,y1+18,150,50,['host removal','bwa-mem2 (opt)'],'io',dashed=True)
for a,b in [(b_np,b_pc),(b_pc,b_ft),(b_ft,b_hr2)]: a_rl(a,b)
a_tb(b_tkV,b_fp)   # viral/prok tracks → QC
# ---- TIER 2 ----
y2=y1+90
lanelabel(y2,'TIER 2',['Assembly &','dedup (hub)'],150)
sr_asm=box(LX+8,y2+14,300,44,['SR assembler','MEGAHIT'],'io')
lr_asm=box(LX+340,y2+14,360,44,['LR assemblers','metaFlye · hifiasm-meta · metaMDBG · Medaka(opt)'],'io')
a_tb(b_hr1,sr_asm); a_tb(b_ft,lr_asm)
mrg=box(LX+8,y2+78,120,44,['merge'],'io')
lf=box(LX+150,y2+78,140,44,['length filter','≥ MIN_CONTIG'],'io')
qa=box(LX+312,y2+78,120,44,['QUAST QC'],'io',dashed=True)
dd=box(LX+454,y2+78,170,44,['MMseqs2 dedup','95% identity'],'io')
hub=box(LX+648,y2+70,300,60,['rep_seq.fasta','CENTRAL HUB (deduplicated reference)','feeds all downstream modules'],'hub')
a_tb(sr_asm,mrg)
for a,b in [(mrg,lf),(lf,qa),(qa,dd),(dd,hub)]: a_rl(a,b)

# ---- TIER 3 ----
y3=y2+150
lanelabel(y3,'TIER 3',['Detection &','binning'],250)
# viral column
EL.append(f'<text x="{LX+8}" y="{y3+16}" font-size="11" font-weight="700" fill="{C["teal"]}" font-family="Helvetica, Arial, sans-serif">VIRAL TRACK</text>')
v_det=box(LX+8,y3+24,340,50,['detection','VirSorter2 · GeNomad · VIBRANT'],'viral')
v_con=box(LX+8,y3+86,340,64,['consensus  (mode: count / score / hybrid)','count: ≥ N tools  ·  score: VS2 or GeNomad ≥ 0.5','hybrid (default): N tools OR score'],'viral',fs=12)
v_cv=box(LX+8,y3+162,215,50,['CheckV → provirus trim','→ viral_nonredundant.fasta'],'viral')
v_vo=box(LX+240,y3+162,220,50,['skani vOTU','95% ANI + 85% AF → reps + table'],'viral')
v_gate=diamond(LX+240,y3+224,120,54,['viral_min_','quality gate'],stroke=C['teal'])
v_sub=box(LX+380,y3+228,215,46,['annotation subset','reps ≥ tier → Tier 4 tax/host/annot'],'viral',fs=11)
v_vr=box(LX+8,y3+228,215,46,['vRhyme → vMAGs','CheckV on vMAGs'],'viral')
a_tb(v_det,v_con); a_tb(v_con,v_cv); a_rl(v_cv,v_vo); a_tb(v_vo,v_gate); a_rl(v_gate,v_sub)
arrow(v_cv['cx'],v_cv['b'],v_vr['cx'],v_vr['y'],C['teal'])
a_tb(hub,v_det)
EL.append(f'<text x="{LX+8}" y="{y3+292}" font-size="10" font-style="italic" fill="{C["sub"]}" '
          f'font-family="Helvetica, Arial, sans-serif">vConTACT3 keeps its own HQ+/≥10 kb gate</text>')
# prok column
px=LX+640
EL.append(f'<text x="{px}" y="{y3+16}" font-size="11" font-weight="700" fill="{C["green"]}" font-family="Helvetica, Arial, sans-serif">PROKARYOTIC TRACK</text>')
p_map=box(px,y3+24,330,50,['read mapping + depth','bwa-mem2 (SR) / minimap2 (LR)'],'prok')
p_flt=box(px,y3+86,330,50,['filter_viral_for_prok','remove free-living viral (keep provirus)'],'prok')
p_ld=diamond(px,y3+148,140,54,['low_depth_','mode?'],stroke=C['green'])
p_bin=box(px+160,y3+140,300,44,['binning MetaBAT2 + SemiBin2','→ Binette → final_bins'],'prok',fs=11)
p_ps=box(px+160,y3+192,300,34,['pseudo-genome (whole contig set)'],'prok',fs=11)
p_qc=box(px,y3+236,470,46,['CheckM2 MIMAG (HQ ≥90/≤5 · MQ ≥50/≤10) · GUNC chimera','galah derep 95% ANI · GTDB-Tk bac120/ar53'],'prok',fs=11)
a_tb(p_map,p_flt); a_tb(p_flt,p_ld)
arrow(p_ld['r'],p_ld['cy']-8,p_bin['x'],p_bin['cy']); arrow(p_ld['r'],p_ld['cy']+8,p_ps['x'],p_ps['cy'])
arrow(p_bin['cx'],p_bin['b'],p_qc['cx'],p_qc['y'],C['green'])
a_tb(hub,p_map)

# ---- TIER 4 ----
y4=y3+320
lanelabel(y4,'TIER 4',['Taxonomy &','annotation'],150)
t_vt=box(LX+8,y4+16,300,64,['Viral taxonomy — deepest-rank consensus','vConTACT3 · MMseqs2/INPHARED','MMseqs2/custom · GeNomad → merged'],'viral',fs=11)
t_va=box(LX+320,y4+16,220,64,['Viral annotation','Pharokka (PHROGS) · Phold','genome maps'],'viral',fs=11)
t_vd=box(LX+552,y4+16,220,64,['Viral defense','DefenseFinder · dbAPIS','defense islands'],'viral',fs=11)
t_pt=box(LX+8,y4+92,300,58,['Prok taxonomy','Prodigal/bin · GTDB-Tk','MMseqs2 custom LCA → merged'],'prok',fs=11)
t_pa=box(LX+320,y4+92,220,58,['Prok annotation','Bakta (gate) · eggNOG','KEGG-Decoder'],'prok',fs=11)
t_amr=box(LX+552,y4+92,300,58,['AMR & mobile elements','AMRFinderPlus · RGI · DeepARG · argNorm → consensus','ABRicate → VFDB · PlasmidFinder'],'prok',fs=10.5)
t_pd=box(LX+864,y4+92,220,58,['Prok defense','DefenseFinder · AntiDefenseFinder','defense islands'],'prok',fs=11)

# ---- TIER 4b ----
y4b=y4+165
lanelabel(y4b,'TIER 4b',['Integration','viral ↔ prok'],70)
i_ph=box(LX+8,y4b+14,360,44,['PHIST host prediction','k-mer similarity · vOTUs vs prok MAGs'],'io')
i_ar=box(LX+400,y4b+14,420,44,['host ↔ defense arms-race cross-links','link viral anti-defense to host defense systems'],'io')
EL.append(f'<path d="M {i_ph["r"]} {i_ph["cy"]} L {i_ar["x"]} {i_ar["cy"]}" fill="none" stroke="{C["arr"]}" stroke-width="1.6" marker-start="url(#ah)" marker-end="url(#ah)"/>')

# ---- TIER 5 ----
y5=y4b+80
lanelabel(y5,'TIER 5',['Abundance &','diversity'],80)
d_cov=box(LX+8,y5+16,300,50,['CoverM','RPKM / TPM / mean / covered_fraction'],'io')
d_ab=box(LX+320,y5+16,200,50,['vOTU + bin','abundance tables'],'io')
d_al=box(LX+540,y5+16,230,50,['Alpha (per sample)','Shannon · Simpson · Chao1 · richness'],'io')
d_be=box(LX+788,y5+16,220,50,['Beta (between samples)','Bray–Curtis → PCoA'],'io')
d_pr=box(LX+1026,y5+16,210,50,['Procrustes','viral vs prok congruence'],'io')
a_rl(d_cov,d_ab)

# ---- TIER 6 ----
y6=y5+90
lanelabel(y6,'TIER 6',['Outputs &','report'],70)
r_fin=box(LX+8,y6+14,200,46,['finalize & organize','per-sample & per-group'],'io')
r_rep=box(LX+228,y6+14,1010,46,['Interactive HTML report (ECharts + D3)  ·  + MultiQC',
  'Overview · Read QC · Assembly · Viral · Prokaryotic · Taxonomy · Host & Defense · Bins · Abundance · Annotation · Co-assembly · Diversity · About'],'io',fs=11)
a_rl(r_fin,r_rep)

# ================= CO-ASSEMBLY COLUMN =================
cx=COASX
EL.append(f'<rect x="{cx-8}" y="{y0-8}" width="{430}" height="{y5+80-y0}" rx="14" fill="{C["indigobg"]}" '
          f'fill-opacity="0.5" stroke="{C["indigo"]}" stroke-width="1.6"/>')
EL.append(f'<text x="{cx+205}" y="{y0+16}" text-anchor="middle" font-size="12" font-weight="700" fill="{C["indigo"]}" '
          f'font-family="Helvetica, Arial, sans-serif">CO-ASSEMBLY BRANCH (opt-in)</text>')
EL.append(f'<text x="{cx+205}" y="{y0+32}" text-anchor="middle" font-size="10.5" fill="{C["sub"]}" '
          f'font-family="Helvetica, Arial, sans-serif">per sample group via metadata</text>')
cw=400
c1=box(cx,y0+44,cw,44,['Co-assembly (per group)','MEGAHIT (SR) / metaFlye (LR)'],'coas')
c_note=box(cx,y0+98,cw,26,['from viral-filtered contigs'],'coas',bold0=False,fs=11)
c2=box(cx,y0+132,190,58,['Co-binning (SR)','VAMB differential','coverage → group MAGs'],'coas',fs=11)
c3=box(cx+210,y0+132,190,58,['Multi-split','VAMB on concatenated','per-sample assemblies'],'coas',fs=11)
c4=box(cx,y0+200,cw,44,['Group MAG QC','CheckM2 · GTDB-Tk · GUNC · galah derep'],'coas',fs=11)
c5=box(cx,y0+254,cw,58,['Viral consumer (per group)','full viral pipeline: detection→vOTU→vRhyme','taxonomy · annotation · host · defense'],'coas',fs=11)
c6=box(cx,y0+322,cw,44,['Functional layer (per group)','AMR · defense · annotation · PHIST'],'coas',fs=11)
c7=box(cx,y0+376,cw,36,['Per-group final/ outputs'],'coas')
a_tb(c1,c_note)
arrow(c_note['cx']-90,c_note['b'],c2['cx'],c2['y'],C['indigo']); arrow(c_note['cx']+90,c_note['b'],c3['cx'],c3['y'],C['indigo'])
arrow(c2['cx'],c2['b'],c4['cx'],c4['y'],C['indigo']); arrow(c3['cx'],c3['b'],c4['cx'],c4['y'],C['indigo'])
for a,b in [(c4,c5),(c5,c6),(c6,c7)]: a_tb(a,b)
EL.append(f'<text x="{cx+205}" y="{y0+430}" text-anchor="middle" font-size="10" font-style="italic" fill="{C["sub"]}" '
          f'font-family="Helvetica, Arial, sans-serif">runs in parallel to per-sample — no global dereplication</text>')

# ================= TITLE + LEGEND =================
EL.insert(0,f'<text x="{W/2}" y="40" text-anchor="middle" font-size="21" font-weight="700" fill="{C["ink"]}" '
          f'font-family="Helvetica, Arial, sans-serif">VAPOR — modular metagenomic &amp; virome recovery: '
          f'selectable reads / viral / prokaryotic tracks with per-sample and co-assembly modes</text>')

ly=y6+90
EL.append(f'<line x1="8" y1="{ly-14}" x2="{W-40}" y2="{ly-14}" stroke="{C["line"]}"/>')
EL.append(f'<text x="18" y="{ly+6}" font-size="12" font-weight="700" fill="{C["sub"]}" font-family="Helvetica, Arial, sans-serif">LEGEND</text>')
leg=[('Reads track (assembly-free)',C['purplebg'],C['purple'],False),
     ('Viral track',C['tealbg'],C['teal'],False),
     ('Prokaryotic track',C['greenbg'],C['green'],False),
     ('Central hub',C['amberbg'],C['amber'],False),
     ('Co-assembly branch (opt-in)','#ffffff',C['indigo'],False),
     ('Optional module','#ffffff','#94a3b8',True),
     ('Standard / shared step',C['iobg'],'#94a3b8',False)]
lx=110
for label,fill,stroke,dash in leg:
    da=' stroke-dasharray="4 3"' if dash else ''
    EL.append(f'<rect x="{lx}" y="{ly-8}" width="18" height="15" rx="4" fill="{fill}" stroke="{stroke}" stroke-width="1.4"{da}/>')
    EL.append(f'<text x="{lx+24}" y="{ly+5}" font-size="11.5" fill="{C["sub"]}" font-family="Helvetica, Arial, sans-serif">{esc(label)}</text>')
    lx += 34 + len(label)*6.6
# decision glyph
EL.append(f'<polygon points="{lx+8},{ly-9} {lx+18},{ly-1} {lx+8},{ly+7} {lx-2},{ly-1}" fill="#fff" stroke="{C["slate"]}" stroke-width="1.4"/>')
EL.append(f'<text x="{lx+26}" y="{ly+5}" font-size="11.5" fill="{C["sub"]}" font-family="Helvetica, Arial, sans-serif">Decision</text>')

svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" '
       f'font-family="Helvetica, Arial, sans-serif">'
       f'<defs><marker id="ah" markerWidth="9" markerHeight="9" refX="7" refY="3" orient="auto" markerUnits="strokeWidth">'
       f'<path d="M0,0 L7,3 L0,6 Z" fill="{C["arr"]}"/></marker></defs>'
       f'<rect x="0" y="0" width="{W}" height="{H}" fill="#ffffff"/>'
       + ''.join(EL) + '</svg>')
open('/tmp/claude-1000/-home-lucas-metagen-pipe-final-claude/06e92905-01ab-42a0-b0a3-f2b3b41aff0a/scratchpad/vapor_figure.svg','w').write(svg)
print("wrote vapor_figure.svg", len(svg), "bytes")
