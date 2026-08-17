-- Run with: paperclip sql -s chembl "$(< queries/photo_clenbuterol_chembl.sql)"
-- PubMed 40492834 / DOI 10.1021/acs.jmedchem.5c00792
SELECT
  md.chembl_id,
  cr.compound_key,
  cs.canonical_smiles,
  cs.standard_inchi_key,
  a.activity_id,
  a.assay_id,
  ass.chembl_id AS assay_chembl_id,
  ass.description AS assay_description,
  ass.assay_type,
  ass.assay_organism,
  ass.assay_cell_type,
  a.standard_type,
  a.standard_relation,
  a.standard_value,
  a.standard_units,
  a.pchembl_value,
  a.activity_comment
FROM chembl.activities AS a
JOIN chembl.assays AS ass ON ass.assay_id = a.assay_id
JOIN chembl.docs AS d ON d.doc_id = a.doc_id
JOIN chembl.compound_records AS cr ON cr.record_id = a.record_id
JOIN chembl.molecule_dictionary AS md ON md.molregno = a.molregno
JOIN chembl.compound_structures AS cs ON cs.molregno = a.molregno
WHERE d.pubmed_id = 40492834
ORDER BY md.chembl_id, a.assay_id, a.standard_type, a.activity_id;
