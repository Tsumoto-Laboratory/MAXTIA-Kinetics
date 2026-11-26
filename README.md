# MAXTIA-Kinetics
A Simple, Fast, and Comprehensive Data Processing Module for MAXTIA

# MAXTIA Kinetics 3.6.3

This English web workbench compares independent binding experiments on a shared WT protein and binding partner. The single-batch view uses the same structure color scale as the batch reports. The separate 3.4.1 release can be used independently.

## Start locally

Extract the entire archive. Start the included launcher with `python3 /path/to/start_local.py` (Windows: `py -3 path/to/start_local.py`). It finds its own dist folder and opens an available loopback port.

Alternatively, run the following command **in the folder containing this README and the dist folder**. On Windows replace python3 with py -3.

    python3 -m http.server 8002 --bind 127.0.0.1 --directory dist

Open http://127.0.0.1:8002/ in your browser. A 404 at the root usually means the server is looking for dist in the wrong folder. You can also pass its absolute path to --directory. Stop the server with Ctrl+C. All analysis assets are bundled and study data are processed in browser memory. Save the study ZIP before leaving the page.

## Run a study

1. Enter the shared binding partner and, optionally, the shared WT sequence and AlphaFold structure files. The selected chain applies to the whole study. If WT is inferred, unambiguous C-labeled WT sequence records must agree.
2. Add each experiment with a unique name, one Carterra LSA workbook, its NGS CSV, its barcode block, and optionally a v3.4.1 QC policy. A combined NGS file can be reused with a different block. Repeated wells or identical barcode blocks in different batches are supported. Exact duplicate LSA exports are rejected.
3. Describe the experimental condition. Use the same label for batches with equivalent conditions. Different declared conditions or explicit experiment metadata remain separate. Every unnamed condition is unique to its batch. Only one protein, binding partner, and analyte set per workbook are supported per study.
4. Analyze the study. Inspect variants at study, condition and batch level. The built-in example illustrates overlapping variants, variants unique to one run, different barcode blocks and batch-specific WT values.
5. Save the study ZIP. integrated_wells.csv keeps batch ID, normalized well, source NGS well and original row numbers. variant_by_batch.csv contains counts, QC, raw variant and WT medians, and eligible fold values. variant_by_condition.csv and variant_overview.csv cover the union of variants with reliable genotypes. study_manifest.json records the batch IDs, condition groups, original input hashes, effective plate mappings, and WT and QC summaries. Numbered batch folders provide full single-batch outputs and effective policies. A shared structure additionally enables condition_structure_mapping.csv and per-condition SVGs for mapped single substitutions.

## Interpreting fold changes

Calculate each variant median and WT median from eligible measurements in the same batch. The batch variant and WT must both pass aggregate QC for its ratio to contribute. Each qualifying batch contributes once per variant and metric. A condition center is the back-transformed median of log2 batch folds. The descriptive study center applies the same rule to available condition centers. Ranges always span the contributing batches. One quantitative batch stays marked as such. Conditions with different variant coverage can yield centers based on different subsets of experiments.

KD folds below 0.80 and above 1.25 in different batches flag opposing effects. Active and Inactive evidence for a variant is retained; no artificial KD is assigned to Inactive wells. These discordances withhold the combined center, leaving batch measurements, ranges and counts visible. A missing variant is unobserved. A missing or QC-limited WT prevents normalized estimates in the affected batch. Apparent agreement or disagreement is descriptive. The study does not estimate a causal batch correction, a condition-independent affinity, or statistical significance.

## Variant ROIs and kinetics

The ranking tables show up to five eligible single substitutions on each side of the WT KD thresholds. Each batch uses its own WT. Median koff (s⁻¹), kon (M⁻¹ s⁻¹), Rmax (RU), and KD use the same active measurements that passed kinetic QC.

In variant_summary.csv and variant_by_batch.csv, observed_rois lists distinct Carterra ROI identifiers for genotype-confirmed observations, including Inactive, Complex, and QC-limited measurements. median_rois lists the identifiers of the measurements used to calculate the medians. An empty field means no eligible observations for that list. An unreported identifier indicates a missing ROI in the LSA export. integrated_records.csv preserves each original well, ROI, and QC decision.

The study report presents these ranks and raw medians separately for each batch. variant_by_condition.csv and variant_overview.csv include observed_rois_by_batch and median_rois_by_batch as JSON objects keyed by batch name, so a repeated ROI number can be traced to its experiment.

## Structure colors

For each residue, use the signed binding ΔΔG of the eligible single substitution with the largest absolute ΔΔG. Ties use ascending substitution label. The per-batch raw value is RT ln(KD variant / KD WT) in kcal/mol, using the configured temperature (default 298.15 K). Each condition uses the median of its eligible batch ΔΔG values per substitution. Conflicting evidence for the same substitution follows the study QC rules.

Only the color coordinate is normalized. Negative energies are divided by the magnitude of the most negative mapped residue energy; positive energies by the largest positive mapped residue energy. The most negative value is blue, zero is white, and the most positive value is red. All batch, condition and integrated study maps share endpoints. The numerical scale does not set a white interval around zero. The KD-fold thresholds in the analytical report retain their existing purpose.

If no eligible substitution supplies an energy, review is purple, inactive-only is yellow, and observed positions without a passing effect are magenta. Gray identifies positions without attributable observations. Uncertain NGS calls cannot establish a tested residue. The summary for a residue observed only in a multi-substitution clone is magenta because the independent effect is unresolved.

The Analysis report lists the count and WT residue labels for No data, Review, Inactive only, and No PASS data. The integrated map and each condition have their own lists in a multi-batch study. Positions follow the input WT sequence, starting at 1.

Use residue_colors.csv (or residue_colors_C1.csv, etc.) to find the selected substitution, raw ΔΔG, signed color coordinate and hex value. The SVG and ChimeraX scripts use these same values. Open the CXC in a fresh ChimeraX session with the matching structure_model.cif beside it. The ZIP includes a CXC for each condition. structure_color_scale.json records the raw endpoints.

## Hosting

Serve the contents of dist/ from an HTTPS static server. Use application/wasm for the bundled WebAssembly asset. The worker and its dependencies are same-origin. No server upload or LLM is used. Pyodide 0.27.7 and Python license notices are in dist/runtime/.

## Integrated report and structures

The study report ranks the five strongest eligible single substitutions in each direction, using the condition-balanced KD fold (<0.80 for improvements and >1.25 for losses). It shows KD, kon and koff folds, raw ΔΔG, batch and condition support, and observed ROIs by batch. A variant observed in only one batch can appear. Check its coverage before interpreting it as replicated. The same selections are in study_ranked_variants.csv. The original batch reports contain the raw medians and ROI provenance.

With a shared structure, study_structure_mapping.csv joins individual single substitutions to exact WT positions. residue_colors_study.csv contains one row per WT position. Open **structure_study.html** in a browser after extracting the study ZIP to revisit the integrated 3D map, switch to a condition or batch, and inspect the variants. All coordinates, measurements, and viewer code are embedded. The single-run equivalent is **structure_explorer.html**. For ChimeraX, run **structure_study.cxc** with structure_model.cif in the working directory. A PyMOL script (structure_study.pml) and a static SVG are also included.

The interactive explorer opens on All batches. Drag to rotate, scroll or pinch to zoom, or choose Atoms and bonds for side-chain detail. Hover to preview substitutions. Click or select a WT residue to pin its measurements across maps. Resume hover releases the selection. Each card shows variant/WT ratios for KD, kon and koff. A lower KD indicates stronger binding, a larger kon faster association, and a larger koff faster dissociation. Each measurement in a multi-substitution clone belongs to the entire clone. The map assigns individual residue colors from independent single substitutions. Exact WT-to-structure matches are selectable; model identifiers in the inspector preserve the original author numbering and insertion codes.

Study spatial candidates use fixed study-wide defaults: at least one eligible single substitution with study KD fold ≥2.0, or exclusively inactive-only single substitutions at that position. Members with Cα atoms at most 8 Å apart are connected, including via intermediate candidates. The report distinguishes multi-position regions from singletons and counts quantitative loss and qualitative inactive-only evidence separately. These descriptive clusters cannot establish direct contacts. study_clusters.csv lists the components; study_cluster_members.csv gives each residue and its evidence. Per-batch cluster files use that batch’s policy.

Browser rendering uses the bundled 3Dmol.js 2.5.5. Notices are in dist/vendor/ and in each portable HTML map.

## Structure display across browsers

Safari, the in-app browser and downloaded HTML maps use the same molecular scene, depth buffer, lighting and starting camera. The ribbon uses the cartoon surface generated by the bundled 3Dmol.js library from the selected chain. Atoms and bonds uses the same parsed atoms and bonds with shaded spheres and sticks. Canvas 2D displays the calculated image locally. Drag to rotate, scroll or pinch to zoom, and use Reset view to restore the original orientation. Hover or click the visible surface to inspect its WT residue. Lighting shades the surface; the residue base colors correspond to the exported color table.
