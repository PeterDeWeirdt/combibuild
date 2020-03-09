context('design')
library(combibuild)

test_that('get_combos', {
  expect_equal(nrow(get_combos(c('A', 'B', 'C'), self = T)),
               choose(3,2) + 3)
})

test_that('design gene combos', {
  gene_col <- 'Target Gene Symbol'
  guide_col <- 'sgRNA Sequence'
  guide_rank <- 'Pick Order'
  minimal_designs <- example_designs %>% select(gene_col, guide_col, guide_rank)
  all_genes <- unique(minimal_designs[[gene_col]])

  expect_equal(choose(3,2) + 3,
               nrow(design_gene_combos(all_genes, T, NULL, NULL, NULL, NULL, F)))
  expect_equal(2*2,
               nrow(design_gene_combos(all_genes, F, c('EEF2', 'BCL2L1'),
                                       c('EEF2', 'MCL1'), NULL, NULL, F)))
  expect_equal(2,
               nrow(design_gene_combos(all_genes, F, NULL, NULL,
                                       tibble(gene1 = c('EEF2', 'BCL2L1'),
                                              gene2 = c('MCL1', 'EEF2')), NULL, F)))
  expect_equal(3, nrow(design_gene_combos(all_genes, F, NULL, NULL, NULL,
                                          c('EEF2'), F)))
  expect_equal(3*3,
               nrow(design_gene_combos(all_genes, T, NULL, NULL, NULL, NULL, T)))

})

test_that('design guide combos', {
  gene_col <- 'Target Gene Symbol'
  guide_col <- 'sgRNA Sequence'
  guide_rank <- 'Pick Order'
  minimal_designs <- example_designs %>% select(gene_col, guide_col, guide_rank)
  gene_combos <- tibble(gene_x = 'BCL2L1', gene_y = 'MCL1')
  guide_tibble <- minimal_designs %>%
    rename(guide = guide_col, rank = guide_rank)
  expect_equal(5^2, nrow(design_guide_combos(gene_combos, guide_tibble, gene_col, 'all')))
  expect_equal(5, nrow(design_guide_combos(gene_combos, guide_tibble, gene_col, 'rank')))
})

test_that('library design', {
  expect_equal(2*choose(3,2)*5^2 + 3*(5^2-5), nrow(design_combo_lib(example_designs,
                                                      all_by_all_gene = T,
                                                      dual_orientation = T)))
  expect_equal(2*5 + 3*5 + 5, nrow(design_combo_lib(example_designs,
                                        row_genes = c('EEF2'),
                                        col_genes = c('BCL2L1', 'MCL1'),
                                        ref_genes = c(''),
                                        gene_pairs = tibble(gene1 = c('MCL1'),
                                                            gene2 = c('BCL2L1')),
                                        guide_pairing = 'rank')))
})
