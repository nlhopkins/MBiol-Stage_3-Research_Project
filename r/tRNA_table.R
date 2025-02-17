library(tidyverse)

load("environments/clusters.RData")

clusters <- c(
    "ALOXE",
    "",
    "",
    "",
    "Ebersole",
    "",
    "",
    "",
    "HES7",
    "Per1",
    "",
    "TMEM107",
    "",
    "",
    "",
    "Arg-CCG",
    "",
    "",
    "",
    "Glu-TTC",
    "",
    "",
    "",
    "",
    "",
    "",
    "iMET",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "Met",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "SeC",
    ""
)


trnas <- c(
    "tRNA-Lys-TTT-3-5",
    "tRNA-Gln-CTG-1-5",
    "tRNA-Leu-TAG-1-1",
    "tRNA-Arg-TCT-2-1",
    "tRNA-Glu-CTC-1-5",
    "tRNA-Gly-TCC-2-5",
    "tRNA-Asp-GTC-2-5",
    "tRNA-Leu-CAG-1-5",
    "tRNA-Gly-GCC-2-6",
    "tRNA-Ser-CGA-1-1",
    "tRNA-Thr-AGT-5-1",
    "tRNA-Trp-CCA-3-3",
    "tRNA-Ser-GCT-4-3",
    "tRNA-Thr-AGT-1-1",
    "tRNA-Ile-AAT-5-5",
    "tRNA-Arg-CCG-1-1",
    "tRNA-Arg-CCG-1-2",
    "tRNA-Arg-CCG-1-3",
    "tRNA-Arg-CCG-2-1",
    "tRNA-Glu-TTC-1-1",
    "tRNA-Glu-TTC-1-2",
    "tRNA-Glu-TTC-2-1",
    "tRNA-Glu-TTC-2-2",
    "tRNA-Glu-TTC-3-1",
    "tRNA-Glu-TTC-4-1",
    "tRNA-Glu-TTC-4-2",
    "tRNA-iMet-CAT-1-1",
    "tRNA-iMet-CAT-1-2",
    "tRNA-iMet-CAT-1-3",
    "tRNA-iMet-CAT-1-4",
    "tRNA-iMet-CAT-1-5",
    "tRNA-iMet-CAT-1-6",
    "tRNA-iMet-CAT-1-7",
    "tRNA-iMet-CAT-1-8",
    "tRNA-iMet-CAT-2-1",
    "tRNA-Met-CAT-1-1",
    "tRNA-Met-CAT-2-1",
    "tRNA-Met-CAT-3-1",
    "tRNA-Met-CAT-3-2",
    "tRNA-Met-CAT-4-1",
    "tRNA-Met-CAT-4-2",
    "tRNA-Met-CAT-4-3",
    "tRNA-Met-CAT-5-1",
    "tRNA-Met-CAT-6-1",
    "tRNA-Met-CAT-7-1",
    "tRNA-SeC-TCA-1-1",
    "tRNA-SeC-TCA-2-1"
)


func <- c(
    "Insulator[@raab2011; @sizer2022]",
    "",
    "",
    "",
    "Insulator[@Ebersole2011; @sizer2022]",
    "",
    "",
    "",
    "Proximate to ALOXE",
    "Per1",
    "",
    "Insulator[@raab2011; @sizer2022]",
    "",
    "",
    "",
    "Cancer Progression[@Goodarzi2016]",
    "",
    "",
    "",
    "Cancer Progression[@Goodarzi2016]",
    "",
    "",
    "",
    "",
    "Translation Initiation",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "iMet Control",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "REDOX[@Sangha2022]",
    ""
)



df <- list(aloxe,
           ebersole,
           hes,
           per1,
           tmem,
           argccg,
           gluttc,
           imet,
           met,
           sec) %>%
    reduce(full_join) %>%
    mutate(ed_activity = ifelse(ed_activity == "1", "Active", "Inactive")) %>%
    mutate(dox_activity = ifelse(dox_activity == "1", "Active", "Inactive")) %>%
    mutate_if(is.numeric, ~ round(., 2)) %>%
    rename(
        "-Dox Activity" = ed_activity,
        "+Dox Activity" = dox_activity,
        "FOXA1 Fold Change" = fox_fold,
        "H3K27ac Fold Change" = h3_fold
    ) %>% select(c("group", "name","func", "+Dox Activity", "-Dox Activity", "FOXA1 Fold Change","H3K27ac Fold Change"))

df <- df[, c("group", "name", "func", "-Dox Activity", "+Dox Activity", "FOXA1 Fold Change", "H3K27ac Fold Change")]

group <-
    df %>% select("group")

func <-
    df %>% select("func")

dox <-
    df %>% select("+Dox Activity")

ed <-
    df %>% select("-Dox Activity")
fox <-
    df %>% select("FOXA1 Fold Change")

h3 <-
    df %>% select("H3K27ac Fold Change")

save.image(file = "environments/trna.RData")
