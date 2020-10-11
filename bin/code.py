
import os
import sys
import csv


def get_lab_panels_from_file(filepath):
    panels = {}
    with open(filepath) as fh:
        reader = csv.reader(fh, delimiter='\t')
        for line in reader:
            panel, _, gene = line
            if panel.startswith("GEL_"):
                continue
            panels.setdefault(panel, set()).add(gene)
    return panels


def get_panelapp_panels_from_dir(directory):
    panels = {}
    for filename in os.listdir(directory):
        if not filename.endswith(".tsv"):
            continue
        filepath = os.path.join(directory, filename)
        panel = get_panelapp_panel_from_file(filepath) 
        panels.update(panel)   
    return panels


def get_panelapp_panel_from_file(filepath):
    panel_genes = {}
    with open(filepath) as fh:
        reader = csv.reader(fh, delimiter='\t')
        for line in reader:
            _, panel, gene = line
            panel_genes.setdefault(panel, set()).add(gene)
    return panel_genes


def compare_panels(lab_panels, panelapp_panels):
    
    pairwise_scores = {}

    for lab_panel_name, lab_panel_genes in sorted(lab_panels.items()):
        for pa_panel_name, pa_panel_genes in sorted(panelapp_panels.items()):
            
            lab_in_pa, pa_in_lab = pairwise_comparison(lab_panel_genes, pa_panel_genes)
            
            pairwise_scores.setdefault(lab_panel_name, {}) 
            pairwise_scores[lab_panel_name][pa_panel_name] = {"lab_in_pa":lab_in_pa,
                                                               "pa_in_lab":pa_in_lab,
                                                               "score":(lab_in_pa * pa_in_lab)}
    return pairwise_scores


def pairwise_comparison(lab_genes, pa_genes):
    shared = lab_genes.intersection(pa_genes)
    num_shared = len(shared)
    frac_of_lab_in_pa = float( num_shared / len(lab_genes) )

    missing_from_pa = lab_genes.difference(pa_genes)

    missing_from_lab = pa_genes.difference(lab_genes)
    frac_of_pa_in_lab = float( num_shared / len(pa_genes) )
    
    return (frac_of_lab_in_pa, frac_of_pa_in_lab)


def output_results(pairwise_scores):
    for lab_panel, score_data in pairwise_scores.items():
        print("\n")
        print(lab_panel)
        scores = [(x, score_data[x]["score"]) for x in score_data.keys()]
        scores = sorted(scores, key=lambda tup: tup[1], reverse=True)
        for pa_panel, pa_score in (scores[:10]):
            lab_in_pa = "%.4f" % score_data[pa_panel]["lab_in_pa"]
            pa_in_lab = "%.4f" % score_data[pa_panel]["pa_in_lab"] 
            print("\t".join(["%.4f" % pa_score, lab_in_pa, pa_in_lab, pa_panel]))



panelapp_panels = get_panelapp_panels_from_dir(sys.argv[1])
lab_panels = get_lab_panels_from_file(sys.argv[2])qq
pairwise_scores = compare_panels(lab_panels, panelapp_panels)
output_results(pairwise_scores)