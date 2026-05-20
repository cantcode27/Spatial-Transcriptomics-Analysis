# Spatial Transcriptomics Analysis — 10x Genomics Visium & Xenium

Four notebooks covering spatial transcriptomics from first principles to advanced spatial statistics. Built on top of official Scanpy and Squidpy tutorials, extended with additional visualizations and biological annotation throughout.

The common thread across all four: spatial context is not just a pretty overlay on top of gene expression — it's a biological dimension in its own right that constrains, validates, and often surprises what you'd conclude from expression data alone.

---

## Background

Standard scRNA-seq tells you what genes a cell expresses. What it can't tell you is where that cell lives in the tissue, what its neighbors are, or how its gene expression relates to visible anatomy. Spatial transcriptomics closes that gap by pairing the gene expression measurement with a physical coordinate.

Two platforms are covered here, and they represent fundamentally different tradeoffs:

**Visium** places tissue on a slide printed with spatially barcoded capture spots (55 µm diameter). RNA from the tissue diffuses into the spots and is sequenced. You get the whole transcriptome, but each spot contains multiple cells — so you're measuring the average expression of whatever cell types happen to fall over that spot. Near-single-cell resolution, not true single-cell.

**Xenium** works differently: fluorescent probes are hybridized directly to RNA molecules inside the intact tissue, and each transcript is detected and localized in situ. True single-cell resolution (each transcript assigned to a segmented cell), but limited to a targeted panel of ~480 genes rather than the full transcriptome.

Neither platform is strictly better. They're complementary. Visium gives broad transcriptome coverage with some spatial averaging; Xenium gives precise single-cell spatial mapping with targeted gene coverage.

---

## Repository Layout

```
spatial_transcriptomics_analysis/
│
├── 01_scanpy_basic/            # Visium: QC → clustering → spatial visualization
│   ├── scanpy_visium_basic.ipynb
│   └── results/
│
├── 02_squidpy_visium_fluo/     # Visium: fluorescence image analysis & feature extraction
│   ├── squidpy_visium_fluo.ipynb
│   └── results/
│
├── 03_squidpy_visium_hne/      # Visium: spatial statistics on H&E data
│   ├── squidpy_visium_hne.ipynb
│   └── results/
│
└── 04_squidpy_xenium/          # Xenium: true single-cell spatial analysis
    ├── tutorial_xenium.ipynb
    └── results/
```

---

## Notebooks

| # | Dataset | Key Methods |
|---|---------|-------------|
| 1 | Human Lymph Node (Visium) | QC, normalization, Leiden clustering, spatial overlay, marker genes |
| 2 | Mouse Brain — fluorescence (Visium) | Nucleus segmentation, image feature extraction, cross-modal clustering |
| 3 | Mouse Brain — H&E (Visium) | Neighborhood enrichment, co-occurrence, Moran's I |
| 4 | Human Lung Cancer (Xenium) | Centrality scores, co-occurrence, neighborhood enrichment, Moran's I |

---

## Notebook 1 — Basic Spatial Analysis with Scanpy

**Dataset:** Human Lymph Node, 10x Genomics Visium  
**Notebook:** `01_scanpy_basic/scanpy_visium_basic.ipynb`

The starting point. This notebook runs a standard single-cell analysis pipeline on a Visium dataset, but with the spatial coordinate system carried along the whole way so that every result — clusters, QC metrics, individual gene expression — can be mapped back onto the tissue image.

Dataset loaded via `sc.datasets.visium_sge("V1_Human_Lymph_Node")` which returns an AnnData object containing:
- `adata.X` — count matrix (spots × genes)
- `adata.obsm['spatial']` — (x, y) coordinates of each spot on the tissue
- `adata.uns['spatial']` — the H&E tissue image and scale factors

Starting dimensions: **4,035 spots × 36,601 genes**

---

### QC

Three metrics calculated per spot: number of genes detected (`n_genes_by_counts`), total UMI count (`total_counts`), and mitochondrial read percentage (`pct_counts_mt`). Mitochondrial genes are flagged by the `MT-` prefix.

#### Violin Plots
![Violin Plots](01_scanpy_basic/results/violin_plots.png)

Each violin shows the distribution of a QC metric across all spots. The width at any value indicates how many spots have that value. What to look for: spots with very few genes are probably empty tissue or damaged; spots with very high mitochondrial percentage are likely dying cells where cytoplasmic RNA has degraded, leaving mitochondrial RNA disproportionately behind.

#### Count & Gene Histograms
![QC Histograms](01_scanpy_basic/results/genes_by_counts_Histogram.png)

Four panels arranged 2×2. The full-range panels show the overall count/gene distribution across all spots. The zoomed panels crop to spots below 10,000 counts and 4,000 genes respectively — necessary because a handful of high-count outlier spots compress the x-axis in the full-range view, hiding the low-end structure where filtering decisions actually need to be made.

**Filtering applied:**

| Parameter | Threshold | Rationale |
|-----------|-----------|-----------|
| Min total counts | 5,000 | Below this → likely empty or damaged tissue |
| Max total counts | 35,000 | Above this → likely doublets |
| Max MT% | 20% | High MT% → stressed/dying cells |
| Min cells per gene | 10 | Removes genes with negligible coverage |

Results: 44 spots removed for low counts, 130 removed for high counts, 16,916 genes filtered. **3,861 spots × 19,685 genes** remain.

---

### Normalization & Feature Selection

Library-size normalization (`normalize_total`) scales each spot to 10,000 total counts, correcting for sequencing depth differences. Log1p transformation (`log1p`) compresses the dynamic range. Highly variable gene selection (Seurat method, top 2,000 genes) retains only the genes that vary most across spots — the ones that carry signal for distinguishing cell populations.

---

### Dimensionality Reduction & Clustering

PCA on the 2,000 HVGs → neighborhood graph → UMAP for visualization → Leiden clustering for community detection. Standard single-cell workflow; the key difference here is that spots rather than cells are the unit of analysis.

#### UMAP Plots
![UMAP Plots](01_scanpy_basic/results/UMAP_plots.png)

Three panels: total counts, genes by counts, and Leiden clusters. The first two panels are a sanity check — if sequencing depth or gene detection were driving the cluster separation, high-count spots would form their own island in UMAP space. Here both metrics are distributed evenly across the clusters, confirming that the Leiden communities reflect transcriptional biology rather than technical variation.

---

### Spatial Visualization

#### QC Metrics on Tissue
![Spatial Overlay — QC](01_scanpy_basic/results/spatial_overlay.png)

Total counts and gene detection mapped back onto the tissue image. Even coverage across the section confirms the slide preparation was uniform — no technical gradients confounding downstream results.

#### Leiden Clusters on Tissue
![Spatial Overlay — Clusters](01_scanpy_basic/results/spatial_overlay2.png)

The key result of notebook 1. Leiden clustering was computed purely in gene expression space — no coordinates involved. Yet the resulting clusters are spatially coherent: spots in the same cluster tend to occupy the same anatomical region of the lymph node. Spatial coherence that emerges from a coordinate-blind analysis is biological validation that the clusters are real.

#### Zoomed View — Clusters 5 and 9
![Spatial Overlay — Zoom](01_scanpy_basic/results/spatial_overlay3_alpha.png)

Upper-right region of the tissue, `alpha=0.5` transparency on the colored spots so the underlying H&E morphology shows through. Cluster 9 spots are concentrated in small, discrete, circular foci — morphologically consistent with germinal centers.

---

### Marker Gene Analysis

`sc.tl.rank_genes_groups` runs a t-test per cluster against all others, ranking genes by differential expression. Top 10 markers for cluster 9 visualized in a heatmap.

#### Marker Gene Heatmap
![Heatmap](01_scanpy_basic/results/heatmap.png)

Each column is a gene, each row is a cluster. Deep color in the cluster 9 row with low expression elsewhere indicates genes specific to that cluster. Specificity here validates that Leiden found biologically distinct populations, not arbitrary groupings.

#### CR2 Expression on Tissue
![Spatial — CR2](01_scanpy_basic/results/spatial_overlay4_CR2.png)

CR2 (Complement Receptor 2 / CD21) is a surface marker of mature B cells and follicular dendritic cells. Its spatial expression pattern mirrors cluster 9 exactly — confirming that cluster 9 corresponds to germinal center B cells. Germinal centers are where B cells undergo affinity maturation in lymph nodes; their location in discrete foci matches the morphology seen in the zoomed spatial overlay. This gene-to-anatomy match is the strongest validation the pipeline can produce.

#### COL1A2 and SYPL1
![Spatial — COL1A2/SYPL1](01_scanpy_basic/results/spatial_overlay5_COL1A2_SYPL1.png)

Two additional genes at `alpha=0.7` transparency. COL1A2 (Collagen Type I Alpha 2) marks stromal fibroblasts — expected in the fibrous capsule and trabeculae surrounding the lymph node, which is what the spatial plot shows. SYPL1 (Synaptophysin-like 1) marks a distinct, non-overlapping region, illustrating that different cell populations occupy spatially segregated niches.

---

## Notebook 2 — Fluorescence Image Analysis with Squidpy

**Dataset:** Mouse Brain Coronal Section — Visium + fluorescence (pre-processed crop)  
**Notebook:** `02_squidpy_visium_fluo/squidpy_visium_fluo.ipynb`

Where notebook 1 treated the tissue image as a backdrop for displaying gene expression results, this notebook treats the image itself as data. The question being answered: does the fluorescence image, analyzed independently of gene expression, recover the same spatial structure that gene expression clustering finds?

Two linked objects: `adata` (gene expression + pre-annotated clusters) and `img` (an `ImageContainer` holding the fluorescence tissue image). The crop covers hippocampal and cortical regions of the mouse brain.

---

### Spatial Cluster Reference
![Cluster Annotation](02_squidpy_visium_fluo/results/1_cluster_annotation_spatialoverlay.png)

Pre-annotated gene expression clusters visualized as the reference. Major anatomical regions are clearly separated — hippocampal layers, cortex, fiber tracts, lateral ventricles — each assigned distinct clusters based on transcriptional profile.

---

### Fluorescence Channels
![Three Fluorescence Channels](02_squidpy_visium_fluo/results/2_three_fluoresence_channels.png)

Three channels in the image, each targeting a different biological component:

**Channel 0 — DAPI:** DNA-binding dye that stains all cell nuclei. Present in every cell; its intensity and density reflect cell packing across the tissue. Used as the channel for nucleus segmentation.

**Channel 1 — anti-NeuN:** Antibody against Neuronal Nuclei (NeuN), a protein expressed in the nuclei and perikarya of most post-mitotic neurons. High NeuN signal → neuronal tissue. Cortex and hippocampus should be strongly NeuN-positive.

**Channel 2 — anti-GFAP:** Antibody against Glial Fibrillary Acidic Protein, expressed in astrocytes. High GFAP signal → glial-rich areas such as fiber tracts and the subventricular zone around the lateral ventricles.

Each channel is independently informative about a different cell class. Together they offer a protein-level spatial cell-type map that is entirely orthogonal to gene expression — making cross-validation meaningful.

---

### Nucleus Segmentation
![Watershed Segmentation](02_squidpy_visium_fluo/results/3_segmented_watershed.png)

**Pre-processing:** DAPI channel smoothed with `sq.im.process(method="smooth")` to reduce noise before segmentation.

**Segmentation:** `sq.im.segment(method="watershed")` applied to the smoothed DAPI channel. The watershed algorithm treats the image as a topographic surface — brighter pixels are peaks — and fills from local intensity minima, separating touching nuclei along their shared intensity ridges.

Left panel: raw DAPI signal in a 500×500 px region, nuclei visible as bright ellipsoidal objects. Right panel: same region post-segmentation, each nucleus labeled with a distinct integer. Qualitative check: each visible nucleus should become exactly one labeled object with no merging across adjacent nuclei and no fragmentation within a single nucleus.

---

### Segmentation Features
`sq.im.calculate_image_features(features="segmentation")` extracts per-Visium-spot measurements from the segmentation:

- Nucleus count per spot (`segmentation_label`) — cell density proxy
- Mean DAPI intensity within segmented nuclei per spot
- Mean NeuN intensity within segmented nuclei per spot
- Mean GFAP intensity within segmented nuclei per spot

![Segmentation Feature Plots](02_squidpy_visium_fluo/results/4_features_segmentation_plots.png)

Four spatial panels:

**Cell count per spot (top-left):** Cell density is spatially heterogeneous. The hippocampal pyramidal layer shows substantially higher cell counts per spot than surrounding tissue — a level of sub-regional detail that gene expression clusters, which group the entire hippocampus into a single cluster, cannot capture. Image-derived cell counting adds sub-cluster spatial resolution.

**Gene clusters (top-right):** Reference annotation for comparison.

**NeuN intensity (bottom-left):** Cortex_1 and Cortex_3 show substantially higher anti-NeuN signal than other clusters — consistent with their transcriptional identity as neuronal populations. Two independent modalities (gene expression and protein fluorescence) agree.

**GFAP intensity (bottom-right):** Fiber_tracts and lateral ventricles show elevated anti-GFAP signal — consistent with the known glial enrichment of white matter and the subventricular zone. Again, gene expression and fluorescence agree.

---

### Multi-Scale Image Feature Extraction & Clustering

Three feature types extracted to characterize image texture beneath each spot:

| Feature type | What it captures |
|---|---|
| Summary | Per-channel pixel intensity statistics (mean, std, percentiles) |
| Histogram | Full intensity distribution across bins — spread, not just mean |
| Texture | GLCM statistics: homogeneity, contrast, entropy — spatial patterning in the image |

Extracted at three scales: spot-boundary masked at full resolution, full-resolution with neighborhood context, and downsampled to 25% resolution. All features concatenated → PCA → neighbor graph → Leiden, producing image-derived cluster annotations independent of gene expression.

![Feature Clustering](02_squidpy_visium_fluo/results/5_features_summary_histogram_texture_cluster.png)

**Summary clusters:** Recover major anatomical regions correctly. Broadly separate hippocampus from cortex and fiber tracts. Tend to subdivide the hippocampus into sub-regions that gene clusters merge into one — reflecting real heterogeneity in cell density and staining intensity within the hippocampal field.

**Histogram clusters:** Using the full intensity distribution rather than mean values makes cortical layering more distinct — intensity distributions vary across cortical depths even when mean intensities are similar.

**Texture clusters:** GLCM features detect microstructural patterns invisible to intensity-based methods. The hippocampal pyramidal layer (densely packed, regularly arranged neurons) vs. surrounding strata becomes clearly visible — sub-structure absent from gene cluster annotations entirely.

**Conclusion:** Image-derived and gene-derived clusterings agree on major regions, validating both. Image features expose finer sub-structure within regions that appear transcriptionally homogeneous. The two modalities are complementary.

---

### Quantitative Cross-Modal Validation — Violin Plots
![Violin Plots](02_squidpy_visium_fluo/results/violinplots.png)

Four metrics (DAPI, NeuN, GFAP intensity, cell count per spot) plotted as violins per gene expression cluster. This converts the spatial maps' visual impression into a quantitative comparison. Neuron-annotated clusters should be NeuN-high; glial-annotated clusters should be GFAP-high; cell-dense anatomical regions should have higher per-spot cell counts. These violins confirm all three expectations — gene expression clusters and fluorescence image features carry consistent biological information about the tissue.

---

## Notebook 3 — Spatial Graph Statistics with Squidpy (H&E)

**Dataset:** Mouse Brain Coronal Section — Visium + H&E stain (full section)  
**Notebook:** `03_squidpy_visium_hne/squidpy_visium_hne.ipynb`

The previous two notebooks established how to cluster spots and visualize gene expression spatially. This notebook asks formal spatial biology questions: which cell types are statistically enriched as neighbors? At what distance do clusters co-occur? Which genes are non-randomly distributed?

H&E staining (Haematoxylin and Eosin) stains nuclei blue and cytoplasm pink — revealing tissue morphology without molecular specificity. The gene expression data provides the molecular layer on top.

---

### Cluster Reference + Image Feature Clusters
![Spatial Cluster Annotation](03_squidpy_visium_hne/results/1_spatial_cluster_annotation.png)

Pre-annotated clusters cover major mouse brain structures: Hippocampus, Pyramidal_layer, Pyramidal_layer_dentate_gyrus, Cortex (multiple sub-regions), Fiber_tracts, Lateral_ventricles.

![Gene vs Image Clusters](03_squidpy_visium_hne/results/2_spatial_cluster_and_features_cluster.png)

Summary image features extracted at scales 1.0 and 2.0 → PCA → Leiden, producing image-derived clusters (left) compared to gene expression clusters (right).

**Agreement:** Fiber_tracts are well-recovered by image features — myelinated axon bundles have a distinctive pale, uniform H&E appearance that makes them morphologically distinct. The hippocampal region is broadly recapitulated.

**Divergence:** The cortex is the informative disagreement. Gene expression reveals the layered cortical structure — different transcriptional identities at different cortical depths reflecting the well-established six-layer laminar architecture. Image features instead group the cortex by broader spatial region, because H&E staining at Visium resolution doesn't clearly differentiate adjacent cortical layers morphologically. Gene expression captures cell identity; image features capture morphological character. These are not the same thing.

---

### Neighborhood Enrichment

Which cluster pairs are spatially adjacent more or less than expected by chance?

`sq.gr.spatial_neighbors` builds a connectivity graph connecting each spot to its physical neighbors on the tissue. `sq.gr.nhood_enrichment` runs 1,000 permutations — shuffling cluster labels each time — and compares the observed co-adjacency count for each cluster pair against this null distribution. Positive score = more adjacent than chance; negative = spatially segregated.

![Neighbourhood Enrichment](03_squidpy_visium_hne/results/3_neighbourhood_enrichment.png)

Warm colors = frequent neighbors; cool colors = spatially separated.

The strongest positive enrichment is between **Pyramidal_layer**, **Pyramidal_layer_dentate_gyrus**, and **Hippocampus**. Anatomically expected: the pyramidal layer is the principal cell body layer of the hippocampus, and the dentate gyrus is a hippocampal subfield. Their physical proximity is reflected directly in the enrichment scores — the analysis recovers known neuroanatomy from spatial coordinates and cluster labels alone.

Negative scores between distant anatomical regions (cortex ↔ hippocampus) confirm the graph is meaningful: structures that are anatomically separated are correctly identified as non-neighbors.

---

### Co-occurrence Analysis

Co-occurrence extends neighbor analysis across a range of distances:

$$\text{score} = \frac{p(\text{exp} \mid \text{cond})}{p(\text{exp})}$$

Score > 1: the two clusters co-occur more than chance at that radius. Score = 1: random. Score < 1: spatial avoidance. Run at increasing radii to reveal both *whether* two clusters are associated and *at what spatial scale*.

![Co-occurrence Scores](03_squidpy_visium_hne/results/4_co-occurence_score_plot.png)

Hippocampus as the anchor cluster. **Pyramidal_layer** shows the highest co-occurrence at short distances, peaking at small radii and decaying toward 1 as distance increases — Pyramidal_layer spots are predominantly immediately surrounding Hippocampus spots, consistent with the pyramidal cell layer forming the inner boundary of the hippocampal region.

**Pyramidal_layer_dentate_gyrus** shows a similar but shifted profile, reflecting the dentate gyrus's more distal position relative to the core hippocampal field.

Cortex and Fiber_tracts converge near 1 at all distances — spatially independent of hippocampal spots at every scale.

The distance-resolved view adds information that neighborhood enrichment alone doesn't provide. Two clusters can show similar enrichment scores but completely different distance decay profiles, revealing whether their spatial association is tight and boundary-restricted or spread across a broader zone.

---

### Spatially Variable Genes — Moran's I

Moran's I measures spatial autocorrelation: the degree to which nearby spots have similar expression values.

- **+1** — perfectly clustered: neighbors have very similar expression
- **0** — spatially random
- **-1** — perfectly dispersed: neighbors have opposite expression

`sq.gr.spatial_autocorr(mode="moran")` run on the top 1,000 highly variable genes. Genes with the highest I scores show the strongest, most consistent spatial clustering in their expression patterns.

![Spatially Variable Genes](03_squidpy_visium_hne/results/5_expression_levels.png)

Three top-scoring genes visualized spatially alongside the cluster annotation:

**Olfm1** (Olfactomedin-1): A secreted glycoprotein involved in neural development and synaptic signaling. Expression enriched in hippocampal and cortical regions, consistent with its role in neuronal maturation and the dense synaptic architecture of these areas. High Moran's I reflects tight spatial confinement to neuronal regions.

**Plp1** (Proteolipid Protein 1): The most abundant protein in myelin, produced by oligodendrocytes. Expression spatially confined to fiber tracts — the white matter regions of densely myelinated axons. A classical oligodendrocyte marker mapping precisely to the Fiber_tracts cluster. High Moran's I is expected here: myelinated regions are anatomically discrete, so Plp1-expressing spots are tightly clustered.

**Itpka** (Inositol 1,4,5-trisphosphate 3-kinase A): Involved in calcium signaling and actin dynamics, highly expressed in neurons, particularly pyramidal layers. Spatial pattern tracks the Pyramidal_layer clusters.

All three genes map onto biologically coherent structures. Moran's I successfully identifies genes whose spatial distribution is shaped by tissue architecture — and each recovered gene corroborates the cluster annotations derived from unbiased transcriptional clustering.

---

## Notebook 4 — Xenium Single-Cell Spatial Analysis

**Dataset:** Human Lung Cancer — 10x Genomics Xenium  
**Notebook:** `04_squidpy_xenium/tutorial_xenium.ipynb`

The only notebook using Xenium rather than Visium — and the methodological differences matter. Where Visium spots average expression over multiple cells, Xenium assigns each detected transcript to a segmented individual cell. The result is true single-cell spatial resolution, at the cost of a targeted gene panel (480 genes) rather than whole transcriptome.

Data loaded via `spatialdata-io`'s Xenium reader into a `SpatialData` object, then converted to Zarr format for efficient access. The AnnData table pulled from `sdata.tables` has dimensions:

```
AnnData: 161,000 cells × 480 genes
obs: cell_id, transcript_counts, control_probe_counts,
     cell_area, nucleus_area, nucleus_count
obsm: spatial   ← cell coordinates already set by the Xenium reader
```

Negative control probe contamination: 0.005% (negative DNA probes), 0.0025% (negative decoding codewords) — very low background signal, confirming high assay quality.

---

### QC & Preprocessing

QC metrics calculated per cell: total transcripts, unique transcripts, cell area, nucleus-to-cell area ratio. These four distributions (visualized as histograms) inform filtering thresholds — cells with too few transcripts are likely empty or fragment-sized segmentation artifacts; cells with unusual area ratios may have segmentation errors.

Filtering with `sc.pp.filter_cells` and `sc.pp.filter_genes`, thresholds set from the QC distribution plots. Normalization → log1p → PCA → neighborhood graph → UMAP → Leiden clustering follows the same standard single-cell pipeline.

At 161,000 cells the neighborhood graph computation is expensive — the co-occurrence step alone took approximately 17 minutes with `n_splits=39` set automatically to prevent creation of a 79,670 × 79,670 distance matrix.

---

### Spatial Statistics

With 161,000 individually segmented cells and their coordinates, the spatial statistics here operate at single-cell resolution rather than the spot-averaged resolution of Visium.

#### Centrality Scores

`sq.gr.centrality_scores` computes three graph-theoretic measures per cluster using a Delaunay triangulation-based spatial neighborhood graph (`coord_type="generic"`, `delaunay=True`):

- **Closeness centrality** — how close the cluster's cells are to all other cells in the graph
- **Degree centrality** — fraction of cells from other clusters directly connected to this cluster's cells
- **Clustering coefficient** — degree to which cells of this cluster form tight local neighborhoods with each other

These scores characterize each cell population's spatial relationship to the rest of the tissue — whether it forms isolated islands, is spatially integrated with neighboring types, or occupies a central position in the tissue graph.

#### Co-occurrence (Cluster 3 anchor)

Co-occurrence probabilities computed across increasing radii with `sq.gr.co_occurrence`. At single-cell resolution, the distance axis resolves associations at sub-cellular to tissue-level scales — a more granular picture of spatial organization than the Visium co-occurrence analysis in Notebook 3.

#### Neighborhood Enrichment

`sq.gr.nhood_enrichment` with 1,000 permutations identifies which Leiden clusters are statistically enriched as spatial neighbors. At single-cell resolution, this operates on individual cell-cell relationships rather than spot-level neighborhoods — capturing tight intercellular spatial associations invisible at Visium's 55 µm spot size.

#### Moran's I (Top Results)

Spatial autocorrelation computed on 100 genes across the 161,000-cell spatial graph (5 min 49 sec runtime). Top Moran's I genes:

| Gene | Moran's I | Biological role |
|------|-----------|-----------------|
| AREG | 0.696 | Amphiregulin — EGFR ligand, expressed by cancer epithelial cells |
| MET | 0.683 | Hepatocyte growth factor receptor — frequently amplified in lung cancer |
| ANXA1 | 0.667 | Annexin A1 — anti-inflammatory; expressed in myeloid and epithelial cells |
| EPCAM | 0.633 | Epithelial Cell Adhesion Molecule — pan-epithelial marker |
| DMBT1 | 0.588 | Deleted in Malignant Brain Tumors 1 — tumor suppressor |

All five top-scoring genes are biologically coherent for lung cancer tissue: AREG and MET are established oncogenes in lung adenocarcinoma; EPCAM marks carcinoma cells that would be spatially clustered by definition; ANXA1 marks spatially organized myeloid infiltrates in the tumor microenvironment. High Moran's I for these genes reflects the spatial organization of tumor epithelium and immune infiltrate within the lung cancer section.

---

## Summary of Methods Across Notebooks

| Analysis | NB1 | NB2 | NB3 | NB4 |
|----------|-----|-----|-----|-----|
| QC + filtering | ✓ | ✓ | ✓ | ✓ |
| Normalization + HVGs | ✓ | ✓ | ✓ | ✓ |
| Leiden clustering | ✓ | ✓ | ✓ | ✓ |
| UMAP | ✓ | — | — | ✓ |
| Spatial cluster overlay | ✓ | ✓ | ✓ | ✓ |
| Marker genes | ✓ | — | — | — |
| Image feature extraction | — | ✓ | ✓ | — |
| Nucleus segmentation | — | ✓ | — | — |
| Centrality scores | — | — | — | ✓ |
| Neighborhood enrichment | — | — | ✓ | ✓ |
| Co-occurrence | — | — | ✓ | ✓ |
| Moran's I | — | — | ✓ | ✓ |

---

## Dependencies

```bash
# Notebooks 1–3
pip install scanpy squidpy igraph leidenalg scikit-image dask

# Notebook 4
pip install spatialdata spatialdata-io spatialdata-plot
```

All notebooks developed and tested on Google Colab.

---

## References

- Wolf FA et al. (2018) SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology*. https://doi.org/10.1186/s13059-017-1382-0
- Palla G et al. (2022) Squidpy: a scalable framework for spatial omics analysis. *Nature Methods*. https://doi.org/10.1038/s41592-021-01358-2
- Traag VA, Waltman L, van Eck NJ (2019) From Louvain to Leiden: guaranteeing well-connected communities. *Scientific Reports*. https://doi.org/10.1038/s41598-019-41695-z
- McInnes L, Healy J, Melville J (2018) UMAP: Uniform Manifold Approximation and Projection. *arXiv*. https://arxiv.org/abs/1802.03426
- Moran PAP (1950) Notes on continuous stochastic phenomena. *Biometrika*. https://doi.org/10.2307/2332142
- Scanpy spatial tutorial: https://scanpy-tutorials.readthedocs.io/en/latest/spatial/basic-analysis.html
- Squidpy Visium fluorescence tutorial: https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_visium_fluo.html
- Squidpy Visium H&E tutorial: https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_visium_hne.html
- Squidpy Xenium tutorial: https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_xenium.html
