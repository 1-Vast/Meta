SELECT
    act.activity_id,
    ass.assay_id,
    ass.chembl_id AS assay_chembl_id,
    doc.chembl_id AS document_chembl_id,
    act.standard_type AS endpoint_family,
    ass.assay_organism,
    ass.bao_format,
    ass.cell_id,
    ass.tissue_id,
    ass.assay_subcellular_fraction AS subcellular_fraction,
    ass.relationship_type,
    ass.variant_id,
    ass.description AS assay_description
FROM selected_activity AS selected
JOIN activities AS act ON act.activity_id = selected.activity_id
JOIN assays AS ass ON ass.assay_id = act.assay_id
LEFT JOIN docs AS doc ON doc.doc_id = COALESCE(act.doc_id, ass.doc_id)
ORDER BY act.activity_id;
