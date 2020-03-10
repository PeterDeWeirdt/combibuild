#' return all combinations of list, as a tibble
#'
#' @param ar array of any type
#' @param m number of elements to choose
#' @param self whether to include the same element twice
#' @param reverse whether to flip the construct
#' @return tibble with element1 and element2
get_combos <- function(ar, self = T) {
  ar_combos <- combn(ar, 2)
  combo_mat <- t(ar_combos)
  colnames(combo_mat) <- c('el1', 'el2')
  combo_tibble <- tibble::as_tibble(combo_mat)
  if (self) {
    self_tib <- tibble::tibble(el1 = ar, el2 = ar)
    combo_tibble <- bind_rows(combo_tibble, self_tib)
  }
  return(combo_tibble)
}

design_gene_combos <- function(all_genes,
                               all_by_all_gene, row_genes, col_genes,
                               gene_pairs, ref_genes, dual_orientation) {
  aba_combos <- NULL
  row_col_combos <- NULL
  ref_combos <- NULL
  if (all_by_all_gene) {
    aba_combos <- get_combos(all_genes) %>%
      rename(gene_x = el1, gene_y = el2)
  }
  if (!is.null(row_genes) & !is.null(col_genes)) {
    row_col_combos <- tidyr::crossing(col_genes, row_genes) %>%
      rename(gene_x = col_genes, gene_y = row_genes)
  }
  if (!is.null(gene_pairs)) {
    colnames(gene_pairs) <- c('gene_x', 'gene_y')
  }
  if (!is.null(ref_genes)) {
    query_genes <- all_genes[!(all_genes %in% ref_genes)]
    ref_combos <- tidyr::crossing(query_genes, ref_genes) %>%
      rename(gene_x = query_genes, gene_y = ref_genes)
  }
  gene_combos <- bind_rows(aba_combos, row_col_combos, gene_pairs, ref_combos) %>%
    distinct()
  if (dual_orientation) {
    rev_combos <- gene_combos %>%
      mutate(temp_x = gene_y, temp_y = gene_x) %>%
      select(-gene_x, -gene_y) %>%
      rename(gene_x = temp_x, gene_y = temp_y)
    gene_combos <- bind_rows(gene_combos, rev_combos) %>%
      distinct()
  }
  return(gene_combos)
}

design_guide_combos <- function(gene_combos, guide_tibble, gene_col, guide_pairing) {
  guide_combos <- gene_combos %>%
    left_join(guide_tibble, by = c('gene_x' = gene_col)) %>%
    left_join(guide_tibble, by = c('gene_y' = gene_col),
              suffix = c('_x', '_y'))
  if (guide_pairing == 'all') {
    # Don't do anything
  } else if (guide_pairing == 'rank') {
    guide_combos <- guide_combos %>%
      filter((rank_x == rank_y) | (is.na(rank_x) & !is.na(rank_y)) |
               (is.na(rank_y) & !is.na(rank_x)))
  } else {
    stop('guide_pairing argument not recognized. Options: all, rank')
  }
  guide_combos <- ungroup(guide_combos) %>%
    filter((guide_x != guide_y) | is.na(rank_x) | is.na(rank_y)) %>%
    select(-c(rank_x, rank_y)) %>%
    mutate(guide_x = tidyr::replace_na(guide_x, ''),
           guide_y = tidyr::replace_na(guide_y, ''))
  return(guide_combos)
}

#' design a combinatorial library
#'
#' @importFrom magrittr %>%
#' @import dplyr
#' @param design_file from GPP's CRISPick
#' @param all_by_all_gene design the simplest combinatorial library, pairing all
#' genes in the design file with all other genes
#' @param row_genes genes to be paired with col_genes
#' @param col_genes genes to be paired with row_genes
#' @param ref_genes all genes get paired with these
#' @param gene_pairs programmed gene pairs
#' @param guide_pairing either 'all' or 'rank'
#' @param dual_orientation pair guides in both directions
#' @param gene_col column name with genes from CRISPick
#' @param guide_col column name with guides from CRISPick
#' @param guide_rank column name with guide ranks from CRISPick
#' @return combinatorial designs with columns gene1, gene2, guide1, guide2
#' @export
design_combo_lib <- function(design_file,
                             all_by_all_gene = F, row_genes = NULL,
                             col_genes = NULL, ref_genes = NULL,
                             gene_pairs = NULL,
                             guide_pairing = 'all', dual_orientation = F,
                             gene_col = 'Target Gene Symbol',
                             guide_col = 'sgRNA Sequence',
                             guide_rank = 'Pick Order') {
  minimal_designs <- design_file %>% select(gene_col, guide_col, guide_rank)
  all_genes <- unique(minimal_designs[[gene_col]])

  # Genes
  gene_combos <- design_gene_combos(all_genes,
                                    all_by_all_gene, row_genes, col_genes,
                                    gene_pairs, ref_genes, dual_orientation)

  # Guides
  guide_tibble <- minimal_designs %>%
    rename(guide = guide_col, rank = guide_rank)
  guide_combos <- design_guide_combos(gene_combos, guide_tibble, gene_col, guide_pairing)
  return(guide_combos)
}
