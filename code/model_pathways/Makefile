
; Analysis workflow

DATA_DIR="../data"
LEAR_CODE_DIR="./learn_GWAS_paths"



; Load all necessary MSIGDB data.

$[DATA_DIR]/remaining_ids_symbols.txt, $[DATA_DIR]/remaining_ids.txt, $[DATA_DIR]/entrez_id_to_symbol.pkl <- $[PASCAL_DIR]/resources/genesets/msigdb/msigBIOCARTA_KEGG_REACTOME.gmt $[PASCAL_DIR]/resources/genesets/msigdb/msigdb.v4.0.entrez.gmt
    cd $[LEARN_CODE_DIR]; ./load_msigdb.py;


; Work through the pathway data and store it, create a visual picture of q-value distribution.

$[DATA_PASCAL_OUT]/path_scores.pkl $[DATA_PASCAL_OUT]/qvals_collated.txt $[DATA_PASCAL_OUT]/qvals_collated.pdf $[DATA_PASCAL_OUT]/qvals_persublist.txt $[DATA_PASCAL_OUT]/qvals_persublist.pdf $[DATA_PASCAL_OUT]/qvals_perpathway.txt $[DATA_PASCAL_OUT]/qvals_perpathway.pdf <-
    cd $[LEARN_CODE_DIR]; ./load_PASCAL_pathway_scores.py; R CMD BATCH ./hist_q_vals.R;


; Load the networks

