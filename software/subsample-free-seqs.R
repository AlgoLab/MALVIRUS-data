#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=5) {
  stop("Two arguments must be supplied (input file).n", call.=FALSE)
}

minout <- as.integer(args[4])
maxout <- as.integer(args[5])

#install.packages("tidyverse")
library("tidyverse");

d <- readr::read_csv(
              args[1],
              col_types = cols(
                lineage = col_factor(),
                status = col_factor(),
                note = col_character(),
                conflict = col_double(),
                ambiguity_score = col_double()
              )
            );
d
summary(d)

d <- ( d
  %>% filter(status == "passed_qc")
)
d
summary(d)

special_lineages <- c(
  "A.23.1",
  "AV.1",
  "B.1.1.7",
  "B.1.1.318",
  "B.1.324.1",
  "B.1.351",
  "B.1.427",
  "B.1.429",
  "B.1.525",
  "B.1.526",
  "B.1.617.1",
  "B.1.617.2",
  "B.1.617.3",
  "P.1",
  "P.2",
  "P.3"
)

compute_nsel <- function(lineage, cnt) {
  pmin(cnt,
  ifelse(
    cnt <= minout,
    cnt,
    ifelse(
      lineage %in% special_lineages,
      pmin(5*maxout, pmax(5*minout, as.integer(cnt * 0.01))),
      pmin(maxout, pmax(minout, as.integer(cnt * 0.01)))
    )
  ))
}


lineage_counts <- d %>% group_by(lineage) %>% summarize(cnt = n())
lineage_counts %>% arrange(desc(cnt))
lineage_counts %>% filter(lineage %in% special_lineages)

selected_lineages <- ( lineage_counts
  %>% mutate(nsel = compute_nsel(lineage, cnt))
                                        # Get 1% of the sequences for each group
                                        # but at least minout and at most maxout
                                        # (special lineages have 5 times more seqs)
)
selected_lineages %>% summarize(sum(nsel))

set.seed(423) # Allow to replicate the rnd sampling
selected <- ( d
  %>% inner_join(selected_lineages)
  %>% group_by(lineage)
  %>% mutate(pos = sample(1:n())) # Randomize order
  %>% filter(pos <= nsel)         # Pick the first nsel for each group
  %>% select(-pos, -nsel, -cnt)
  %>% ungroup()
)
selected

readr::write_csv(selected, file=args[2], na = "")
write_lines(str_replace_all(selected$taxon, fixed(" "), "_"), file=args[3])

selected <- selected %>% group_by(lineage) %>% summarize(cnt = n())
selected %>% arrange(desc(cnt))
selected %>% filter(lineage %in% special_lineages)

## Execute the following command in the shell to extract the sequences
# ./faSomeRecords ./mmsa_2021-03-02/mmsa_2021-03-02.fa.gz pangolin_selected_gisaid_epi_isl.txt pangolin_selected.fa
