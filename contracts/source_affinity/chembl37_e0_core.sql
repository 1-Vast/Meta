SELECT
    act.activity_id,
    ass.assay_id,
    ass.chembl_id AS assay_chembl_id,
    COALESCE(doc.chembl_id, '') AS document_chembl_id,
    COALESCE(doc.doi, '') AS document_doi,
    doc.pubmed_id AS document_pmid,
    doc.year AS document_year,
    parent.molregno AS parent_molregno,
    parent.chembl_id AS molecule_chembl_id,
    struct.canonical_smiles,
    struct.standard_inchi_key,
    td.chembl_id AS target_chembl_id,
    td.target_type,
    tc.component_id AS target_component_id,
    cs.accession AS target_accession,
    COALESCE(vs.sequence, cs.sequence) AS protein_sequence,
    ass.variant_id AS target_variant_id,
    vs.mutation AS target_variant_mutation,
    act.standard_type AS endpoint_family,
    act.standard_relation,
    act.standard_value,
    act.standard_units,
    act.type AS published_type,
    act.relation AS published_relation,
    act.value AS published_value,
    act.units AS published_units,
    act.pchembl_value AS pchembl_value_reported,
    act.bao_endpoint,
    act.activity_comment,
    act.standard_flag,
    act.data_validity_comment,
    act.potential_duplicate,
    act.src_id,
    ass.assay_type,
    ass.confidence_score AS assay_confidence,
    ass.description AS assay_description,
    ass.assay_organism,
    ass.relationship_type,
    ass.bao_format,
    ass.cell_id,
    ass.tissue_id,
    ass.assay_subcellular_fraction AS subcellular_fraction,
    ass.variant_id
FROM activities AS act
JOIN assays AS ass ON ass.assay_id = act.assay_id
JOIN target_dictionary AS td ON td.tid = ass.tid
JOIN target_components AS tc ON tc.tid = td.tid
JOIN (
    SELECT tid
    FROM target_components
    GROUP BY tid
    HAVING COUNT(*) = 1
) AS single_component ON single_component.tid = td.tid
JOIN component_sequences AS cs ON cs.component_id = tc.component_id
LEFT JOIN variant_sequences AS vs ON vs.variant_id = ass.variant_id
LEFT JOIN molecule_hierarchy AS hierarchy ON hierarchy.molregno = act.molregno
JOIN molecule_dictionary AS parent
    ON parent.molregno = COALESCE(hierarchy.parent_molregno, act.molregno)
JOIN compound_structures AS struct ON struct.molregno = parent.molregno
LEFT JOIN docs AS doc ON doc.doc_id = COALESCE(act.doc_id, ass.doc_id)
WHERE ass.assay_type = 'B'
  AND ass.confidence_score = 9
  AND td.target_type = 'SINGLE PROTEIN'
  AND act.standard_type IN ('Ki', 'Kd')
  AND act.standard_relation = '='
  AND act.standard_value IS NOT NULL
  AND act.standard_value > 0
  AND act.standard_units IN ('M', 'mM', 'uM', 'um', 'µM', 'nM', 'pM', 'fM')
  AND act.standard_flag = 1
  AND act.data_validity_comment IS NULL
  AND COALESCE(act.potential_duplicate, 0) = 0
  AND COALESCE(vs.sequence, cs.sequence) IS NOT NULL
  AND struct.canonical_smiles IS NOT NULL
  AND struct.standard_inchi_key IS NOT NULL
  AND doc.chembl_id IS NOT NULL
ORDER BY act.activity_id;
