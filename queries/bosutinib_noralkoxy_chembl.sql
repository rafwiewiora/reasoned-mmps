-- Paperclip: paperclip sql -s chembl "$(< queries/bosutinib_noralkoxy_chembl.sql)"
-- PMC10726475 / DOI 10.1021/acsmedchemlett.3c00479 / ChEMBL doc_id 128032
SELECT
  md.chembl_id,
  cr.compound_key,
  cs.canonical_smiles,
  cs.standard_inchi_key,
  a.activity_id,
  a.assay_id,
  ass.chembl_id AS assay_chembl_id,
  ass.description,
  a.standard_type,
  a.standard_relation,
  a.standard_value,
  a.standard_units,
  a.pchembl_value
FROM chembl.activities AS a
JOIN chembl.assays AS ass ON ass.assay_id = a.assay_id
JOIN chembl.compound_records AS cr ON cr.record_id = a.record_id
JOIN chembl.molecule_dictionary AS md ON md.molregno = a.molregno
JOIN chembl.compound_structures AS cs ON cs.molregno = a.molregno
WHERE a.doc_id = 128032
  AND md.chembl_id IN ('CHEMBL288441', 'CHEMBL5435819')
ORDER BY a.assay_id, md.chembl_id, a.activity_id;
