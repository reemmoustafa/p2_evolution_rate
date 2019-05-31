from Bio.Phylo.PAML import codeml
cml = codeml.Codeml()
cml.alignment = "Alg_NucSeq_PY_project_2_input_SRSF1.phy"
cml.tree = "project2_tree.newick"
cml.out_file = "results_final.out"
cml.working_dir = r"C:\Users\Reem\bin\final_PAML"
cml.set_options(CodonFreq=2, model=0, seqtype=1, NSsites='0', runmode=0, cleandata=0)
results = cml.run()
ns_sites = results.get("NSsites")
m0 = ns_sites.get(0)
m0_params = m0.get("parameters")
print(m0_params.get("omega"))



