# A1 SPKOP strict-firewall correction

Date: 2026-07-26. This correction was frozen after the completion audit detected that the first A1
run held out KLIFS families and Bemis–Murcko scaffolds but did not explicitly merge cross-scaffold
high-Tanimoto neighbours or use the project's full-sequence homology components. The first result is
retained but is not the authoritative strict-dual-cold result.

The estimator, seed, eight-protein/eight-ligand neighbourhood, label transformation, arms, metric
and pass/fail gates are unchanged. Only the split components change:

- target components: the existing project-wide full-sequence 4-mer-containment union-find,
  threshold 0.40 against the shorter sequence;
- ligand components: a union-find edge for canonical parent identity, equal Bemis–Murcko scaffold,
  or Morgan radius-2 Tanimoto >=0.50.

The correction is deliberately conservative: every high-similarity edge is transitively closed
before five-fold assignment. Assay/document/publication isolation is not identifiable inside the
single KIRHub release. Therefore this remains a within-source mechanism probe, never a
confirmation of cross-assay DTA generalization. Sharing the source makes the task easier; failure
under this condition is conservative.

Literature basis for the frozen standard components:

- ESM-2 representation: Lin et al., *Science* 2023, DOI 10.1126/science.ade2574.
- ECFP/Morgan representation: Rogers and Hahn, *JCIM* 2010, DOI 10.1021/ci100050t.
- Molecular frameworks: Bemis and Murcko, *J. Med. Chem.* 1996,
  DOI 10.1021/jm9602928.
- Tanimoto fingerprint similarity: Bajusz et al., *J. Cheminformatics* 2015,
  DOI 10.1186/s13321-015-0069-3.
- Alignment-free k-mer comparison: Zielezinski et al., *Genome Biology* 2017,
  DOI 10.1186/s13059-017-1319-7.
- Local weighted-average estimator: Nadaraya, *Theory Probab. Appl.* 1964,
  DOI 10.1137/1109020.

This correction consumes no new candidate slot and is not an A1 mechanism revision. It is one
additional low-cost evaluation of the same candidate, still within the round maximum of four.

