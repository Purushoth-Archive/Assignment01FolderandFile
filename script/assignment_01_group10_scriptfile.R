# ==============================================================================
# Title: BINF*6210 Assignment #2 - Collaborative Script Revision
# Primary Author: Yichen Shen (GitHub: yichenhellenshen1201)
#   - Role: Wrote the original Assignment #1 R script, including data import,
#     filtering, summarization, binding, and initial exploratory visualizations.
# Secondary Contributor: Purushoth Puvanendran (GitHub: Purushoth-Archive)
#   - Role: Revised comment formatting to organize outline display, 
# ==============================================================================

##***************************
## Assignment 01 - Group 10 
## 2026-11-02
##***************************

## Packages used -------
library(tidyverse)
conflicted::conflict_prefer("filter", "dplyr")
library(viridis)
# + scale_color/fill_viridis(discrete = T/F)
theme_set(theme_light())


rm(list = ls())


## Import TE annotation datasets -----

## _Import T. patagoniensis dataset
df_patagoniensis <- read_csv("../data/GCA_001447655.1.csv")
df_patagoniensis 

## _Import T. nelsoni dataset
df_nelsoni <- read_csv(file = "../data/trichinella_nelsoni_te.csv") 
df_nelsoni

## _Import T. spiralis dataset
df_spiralis = read_csv("../data/trichinella_spiralis_te.csv")
df_spiralis 


## Clean dataset columns -----
df_patagoniensis =
  df_patagoniensis |> 
  mutate(
    qry_left = as.numeric(str_replace_all(qry_left, "\\(|\\)", "")),
    in_repeat_end = as.numeric(str_replace_all(in_repeat_end, "\\(|\\)", "")),
    no_bp_in_complement = as.numeric(str_replace_all(no_bp_in_complement, "\\(|\\)", ""))
  )

df_nelsoni =
  df_nelsoni |> 
  mutate(
    qry_left = as.numeric(str_replace_all(qry_left, "\\(|\\)", "")),
    in_repeat_end = as.numeric(str_replace_all(in_repeat_end, "\\(|\\)", "")),
    no_bp_in_complement = as.numeric(str_replace_all(no_bp_in_complement, "\\(|\\)", ""))
  )

df_spiralis =
  df_spiralis |> 
  mutate(
    qry_left = as.numeric(str_replace_all(qry_left, "\\(|\\)", "")),
    in_repeat_end = as.numeric(str_replace_all(in_repeat_end, "\\(|\\)", ""))
  )


## Define non-TE classes to exclude -----
ls_non_te = c("simple_repeat", "Simple_repeat",
              "Low_complexity", "Satellite",
              "rRNA", "snRNA", "ARTEFACT",
              "Unknown")
ls_non_te


## Sumarize data -----

## _Summarize T. spiralis TE annotations -----
df_spiralis_summary =
  df_spiralis |>
  filter(matching_class %notin% ls_non_te) |>
  group_by(matching_class) |>
  summarise(
    annotation_count = n(),
    total_length = sum(qry_width, na.rm = TRUE)
  ) |>
  mutate(species = "T. spiralis") |> 
  arrange(desc(annotation_count))


df_spiralis_summary


## _Summarize T. patagoniensis TE annotations -----
df_patagoniensis_summary =
  df_patagoniensis |>
  filter(matching_class %notin% ls_non_te) |>
  group_by(matching_class) |>
  summarise(
    annotation_count = n(),
    total_length = sum(qry_width, na.rm = TRUE)
  ) |>
  mutate(species = "T. patagoniensis") |> 
  arrange(desc(annotation_count))
df_patagoniensis_summary

## _Summarize T. nelsoni TE annotations -----
df_nelsoni_summary =
  df_nelsoni |>
  filter(matching_class %notin% ls_non_te) |>
  group_by(matching_class) |>
  summarise(
    annotation_count = n(),
    total_length = sum(qry_width, na.rm = TRUE)
  ) |>
  mutate(species = "T. nelsoni") |> 
  arrange(desc(annotation_count))
df_nelsoni_summary


## COMBINING -----

## _Combine species summaries -----
df_comparison =
  bind_rows(df_spiralis_summary,
            df_patagoniensis_summary, 
            df_nelsoni_summary)
df_comparison


## _Combine raw TE annotation datasets -----
df_TE_combined =
  bind_rows(
    df_spiralis |> mutate(species = "T. spiralis"),
    df_patagoniensis |> mutate(species = "T. patagoniensis"),
    df_nelsoni |> mutate(species = "T. nelsoni")
  )
df_TE_combined


## RANKING TE CLASSES -----

## _Rank TE classes by annotation count -----
df_comparison |>
  group_by(species) |>
  arrange(desc(annotation_count), .by_group = TRUE)


## _Rank TE classes by total annotated length -----
df_comparison |>
  group_by(species) |>
  arrange(desc(total_length), .by_group = TRUE)


# _Top 3 TE classes by total length for each species
df_comparison |> 
  group_by(species) |> 
  arrange(desc(total_length), .by_group = TRUE) |> 
  slice_head(n = 3)


## PLOTTING -----

## _Plot distribution of TE annotation lengths -----
df_TE_combined |>
  filter(matching_class %notin% ls_non_te) |>
  ggplot(aes(x = qry_width)) +
  geom_histogram() +
  facet_grid(~ species) +
  labs(
    title = "Distribution of TE Annotation Lengths",
    x = "Annotated TE length (bp)",
    y = "Frequency"
  )


## _Plot annotation counts against total annotated lengths -----
df_comparison |> 
  ggplot(aes(x = annotation_count, y = total_length)) +
  geom_point() +
  geom_smooth() +
  facet_grid(~ species) +
  labs(title = "Comparison of TE Abundance & Length",
       x = "annotated count",
       y = "annotated sequence length(bp)")


## _Plot average TE length by class -----
df_comparison |>
  ggplot(aes(x = matching_class,
             y = total_length / annotation_count)) +
  geom_col() +
  coord_flip() +
  facet_grid(~ species) +
  labs(
    title = "Average TE Length by Class",
    x = "TE class",
    y = "Average annotated length (bp)"
  )

