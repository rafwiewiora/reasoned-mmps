-- Paperclip: paperclip sql -s chembl "$(< queries/biib129_chembl.sql)"
-- PMC11129193 / DOI 10.1021/acs.jmedchem.4c00220
-- ChEMBL document CHEMBL5500379 / raw doc_id 128824
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
WHERE a.doc_id = 128824
ORDER BY a.assay_id, cr.compound_key, a.activity_id;
