#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=6) {
  stop("Two arguments must be supplied (input file).n", call.=FALSE)
}

mincl <- as.integer(args[4])
minout <- as.integer(args[5])
maxout <- as.integer(args[6])

#install.packages("tidyverse")
library("tidyverse");

d <- readr::read_tsv(
              args[1],
              col_types = cols(
                `Sequence length` = col_integer(),
                `Collection date` = col_character(),
                `Variant` = col_character()
              )
            ) %>%
  dplyr::rename(pangolin_lineage=`Pango lineage`);

special_lineages <- c("B.1.1.7", "B.1.351", "P.1", "A.23.1", "B.1.525")

compute_nsel <- function(pangolin_lineage, cnt) {
  ifelse(
    cnt <= minout,
    cnt,
    ifelse(
      pangolin_lineage %in% special_lineages,
      pmin(5*maxout, pmax(5*minout, as.integer(cnt * 0.01))),
      pmin(maxout, pmax(minout, as.integer(cnt * 0.01)))
    )
  )
}


lineage_counts <- d %>% group_by(pangolin_lineage) %>% summarize(cnt = n())

selected_lineages <- ( lineage_counts
  %>% filter(cnt > mincl | pangolin_lineage %in% special_lineages)
  %>% mutate(nsel = compute_nsel(pangolin_lineage, cnt))
                                        # Get 1% of the sequences for each group
                                        # but at least minout and at most maxout
                                        # (special lineages have 5 times more seqs)
)
selected_lineages %>% summarize(sum(nsel))

set.seed(423) # Allow to replicate the rnd sampling
selected <- ( d
  %>% inner_join(selected_lineages)
  %>% group_by(pangolin_lineage)
  %>% mutate(pos = sample(1:n())) # Randomize order
  %>% filter(pos <= nsel)         # Pick the first nsel for each group
  %>% select(-pos, -nsel, -cnt)
  %>% ungroup()
)
selected

readr::write_tsv(selected, file=args[2])
write_lines(str_replace_all(selected$`Virus name`, fixed(" "), "_"), file=args[3])

## Execute the following command in the shell to extract the sequences
# ./faSomeRecords ./mmsa_2021-03-02/mmsa_2021-03-02.fa.gz pangolin_selected_gisaid_epi_isl.txt pangolin_selected.fa
