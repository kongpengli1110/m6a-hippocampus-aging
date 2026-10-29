# Third-party data, attribution and licence boundaries

The base LICENSE is retained for inherited code. It does not relicense any data, annotation, third-party Python distribution or original publication. No new authorship verification was performed.

## Allen Institute / SEA-AD

This research handoff contains limited processed metadata, feature/index/target-count tables and derived gene-background lookups needed to audit the analysis. It does not contain the h5ad, cell-level matrix or a complete clinical/cognitive workbook. Donor IDs are public resource identifiers; do not attempt re-identification.

- Exact dataset: SEAAD HIP RNAseq final-nuclei 2026-06-22; public URL/object version/hash are in `input_manifest.json`.
- Resource publication: Gabitto MI et al. *Integrated multimodal cell atlas of Alzheimer's disease*. Nature Neuroscience 27, 2366–2383 (2024). DOI: https://doi.org/10.1038/s41593-024-01774-5. This resource citation does not imply the 2026 HIP dataset was analysed in that 2024 paper.
- Data/download page: https://brain-map.org/consortia/sea-ad/our-data
- Terms: https://alleninstitute.org/legal/terms-of-use/
- Citation policy: https://alleninstitute.org/citation-policy/

The official data page and Terms of Use were read on 2026-10-28. The terms permit use/copy/distribution/derivatives for research or other noncommercial purposes, subject to the full terms, citation obligations, third-party conditions and the Freedom to Innovate section; commercial redistribution requires written permission. **Do not treat these data as MIT or unrestricted commercial-use data.** This archive is a local noncommercial research-code handoff, not an assertion that every possible redistribution use is licensed.

## HGNC

The bundled exact complete-set snapshot is a required nomenclature lookup, not expression data. Credit **HUGO Gene Nomenclature Committee at the University of Cambridge**, https://www.genenames.org/. Its official policy (read 2026-10-28) states that no restrictions are imposed on access or use and requests attribution: https://www.genenames.org/about. Preserve snapshot identity rather than silently replacing it with today's nomenclature. Derived human crosswalks carry the same source attribution.

## GEO / Affymetrix GPL570

The small original GEO annotation and derived probe/sample lookup tables are provided for research provenance. GSE48350 matrix is not bundled. Cite GEO accession https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48350 and GPL570; the normal-aging cohort lineage is Berchtold et al., PNAS (2008), DOI: https://doi.org/10.1073/pnas.0806883105. Submitted third-party GEO material is not made MIT by this archive. NCBI's policy distinguishes US-government-created public-domain material from third-party content: https://www.ncbi.nlm.nih.gov/home/about/policies/.

## GENCODE

Full GTF is omitted; exact v32 URL, bytes and hash are supplied. Bundled source-gene exon-union lengths and matching diagnostics were derived with the included original technical code. Attribute GENCODE human release 32 and retain its applicable resource terms: https://www.gencodegenes.org/human/release_32.html. No new annotation licence is asserted here.

## Python and inherited mouse inputs

Third-party dependencies are installed by the recipient from their publishers and retain their own licences. No vendored dependencies are included. Base mouse/candidate attribution is inherited from the original package and has not been independently validated by the human-code supplement.
