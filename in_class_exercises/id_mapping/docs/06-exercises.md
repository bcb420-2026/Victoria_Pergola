# Exercises

These are intended to be done **after** completing the worked examples.


``` r
library(dplyr)
library(readr)
library(tibble)
library(knitr)
library(stringr)
library(biomaRt)

source("./fetch_geo_supp.R")

strip_ensembl_version <- function(x) sub("\\..*$", "", x)

safe_read <- function(file) {
	df <- tryCatch(
		readr::read_tsv(file, show_col_types = FALSE),
		error = function(e) NULL
	)

	if (is.null(df)) {
		return(readr::read_table(file, show_col_types = FALSE))
	}

	probs <- readr::problems(df)
	if (nrow(probs) > 0) {
		return(readr::read_table(file, show_col_types = FALSE))
	}

	df
}

find_data_file <- function(gse) {
	fetch_geo_supp(gse = gse)
	path <- file.path("data", gse)
	files <- list.files(
		path,
		pattern = "\\.txt\\.gz$|\\.tsv\\.gz$|\\.csv\\.gz$",
		full.names = TRUE,
		recursive = TRUE
	)

	if (length(files) == 0) {
		stop("No supported supplementary expression file found for ", gse)
	}

	files[1]
}

map_to_hgnc <- function(ids, mart, cache_key) {
	cleaned <- unique(strip_ensembl_version(ids))
	cleaned <- cleaned[grepl("^ENSG", cleaned)]

	if (length(cleaned) == 0) {
		return(tibble(ensembl_gene_id = character(), hgnc_symbol = character()))
	}

	dir.create("cache", showWarnings = FALSE, recursive = TRUE)
	cache_file <- file.path("cache", paste0("hgnc_map_", cache_key, ".rds"))

	live_map <- tryCatch(
		getBM(
			attributes = c("ensembl_gene_id", "hgnc_symbol"),
			filters = "ensembl_gene_id",
			values = cleaned,
			mart = mart
		),
		error = function(e) NULL
	)

	if (!is.null(live_map)) {
		saveRDS(live_map, cache_file)
		return(live_map)
	}

	if (file.exists(cache_file)) {
		message("BioMart unavailable; using cached HGNC mapping from ", cache_file)
		return(readRDS(cache_file))
	}

	message("BioMart unavailable and no HGNC cache found; using empty HGNC fallback.")
	tibble::tibble(
		ensembl_gene_id = as.character(cleaned),
		hgnc_symbol = NA_character_
	)
}

connect_ensembl <- function() {
	mirrors <- c("www", "useast", "uswest", "asia")

	for (m in mirrors) {
		mart <- tryCatch(
			biomaRt::useEnsembl(
				biomart = "genes",
				dataset = "hsapiens_gene_ensembl",
				mirror = m
			),
			error = function(e) NULL
		)

		if (!is.null(mart)) {
			return(mart)
		}
	}

	stop("Could not connect to Ensembl using available mirrors.")
}

safe_getBM <- function(attributes, filters, values, mart) {
	tryCatch(
		getBM(
			attributes = attributes,
			filters = filters,
			values = values,
			mart = mart
		),
		error = function(e) {
			message("BioMart query unavailable in this run; using fallback with empty HGNC symbols.")
			tibble::tibble(
				ensembl_gene_id = as.character(values),
				hgnc_symbol = NA_character_
			)
		}
	)
}

ensembl <- connect_ensembl()
```

## Exercise 1 - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119732

Using **GSE119732**, confirm whether the ID column contains Ensembl IDs with version suffixes.


``` r
gse1_file <- find_data_file("GSE119732")
```

```
## Setting options('download.file.method.GEOquery'='auto')
```

```
## Setting options('GEOquery.inmemory.gpl'=FALSE)
```

```
## Using locally cached version of supplementary file(s) GSE119732 found here:
## data/GSE119732/GSE119732_count_table_RNA_seq.txt.gz
```

``` r
gse1 <- safe_read(gse1_file)

gse1_id_col <- names(gse1)[1]
gse1_ids <- as.character(gse1[[1]])

gse1_first20 <- tibble(raw_id = head(gse1_ids, 20))
gse1_dot_count <- sum(grepl("\\.", gse1_first20$raw_id))

gse1_with_stripped <- gse1 %>%
	mutate(ensembl_gene_id = strip_ensembl_version(.data[[gse1_id_col]]))

gse1_map <- map_to_hgnc(gse1_ids, ensembl, cache_key = "GSE119732")
```

```
## BioMart unavailable; using cached HGNC mapping from cache/hgnc_map_GSE119732.rds
```

``` r
gse1_mapped <- gse1_with_stripped %>%
	left_join(gse1_map, by = "ensembl_gene_id")
```

```
## Warning in left_join(., gse1_map, by = "ensembl_gene_id"): Detected an unexpected many-to-many relationship between `x` and `y`.
## ℹ Row 29117 of `x` matches multiple rows in `y`.
## ℹ Row 53969 of `y` matches multiple rows in `x`.
## ℹ If a many-to-many relationship is expected, set `relationship =
##   "many-to-many"` to silence this warning.
```

1. First 20 IDs:


``` r
knitr::kable(gse1_first20)
```



|raw_id            |
|:-----------------|
|ENSG00000223972.5 |
|ENSG00000227232.5 |
|ENSG00000278267.1 |
|ENSG00000243485.4 |
|ENSG00000237613.2 |
|ENSG00000268020.3 |
|ENSG00000240361.1 |
|ENSG00000186092.4 |
|ENSG00000238009.6 |
|ENSG00000239945.1 |
|ENSG00000233750.3 |
|ENSG00000268903.1 |
|ENSG00000269981.1 |
|ENSG00000239906.1 |
|ENSG00000241860.6 |
|ENSG00000222623.1 |
|ENSG00000241599.1 |
|ENSG00000279928.1 |
|ENSG00000279457.3 |
|ENSG00000273874.1 |

2. Number of first 20 IDs containing a `.`: **20**

3. New column with versions stripped (preview):


``` r
knitr::kable(head(gse1_with_stripped[, c(gse1_id_col, "ensembl_gene_id")], 10))
```



|gene_id           |ensembl_gene_id |
|:-----------------|:---------------|
|ENSG00000223972.5 |ENSG00000223972 |
|ENSG00000227232.5 |ENSG00000227232 |
|ENSG00000278267.1 |ENSG00000278267 |
|ENSG00000243485.4 |ENSG00000243485 |
|ENSG00000237613.2 |ENSG00000237613 |
|ENSG00000268020.3 |ENSG00000268020 |
|ENSG00000240361.1 |ENSG00000240361 |
|ENSG00000186092.4 |ENSG00000186092 |
|ENSG00000238009.6 |ENSG00000238009 |
|ENSG00000239945.1 |ENSG00000239945 |

4. Mapping to HGNC symbols (preview):


``` r
gse1_mapped_nonempty <- gse1_mapped %>%
	filter(!is.na(hgnc_symbol), hgnc_symbol != "") %>%
	distinct(ensembl_gene_id, .keep_all = TRUE)

if (nrow(gse1_mapped_nonempty) > 0) {
	knitr::kable(
		head(gse1_mapped_nonempty[, c(gse1_id_col, "ensembl_gene_id", "hgnc_symbol")], 10),
		caption = "GSE119732: first 10 successfully mapped IDs"
	)
} else {
	knitr::kable(
		head(gse1_mapped[, c(gse1_id_col, "ensembl_gene_id", "hgnc_symbol")], 10),
		caption = "GSE119732: no HGNC mappings returned in this run"
	)
}
```



Table: (\#tab:ex1_map_preview)GSE119732: first 10 successfully mapped IDs

|gene_id           |ensembl_gene_id |hgnc_symbol |
|:-----------------|:---------------|:-----------|
|ENSG00000223972.5 |ENSG00000223972 |DDX11L1     |
|ENSG00000227232.5 |ENSG00000227232 |WASH7P      |
|ENSG00000278267.1 |ENSG00000278267 |MIR6859-1   |
|ENSG00000243485.4 |ENSG00000243485 |MIR1302-2HG |
|ENSG00000237613.2 |ENSG00000237613 |FAM138A     |
|ENSG00000268020.3 |ENSG00000268020 |OR4G4P      |
|ENSG00000240361.1 |ENSG00000240361 |OR4G11P     |
|ENSG00000186092.4 |ENSG00000186092 |OR4F5       |
|ENSG00000233750.3 |ENSG00000233750 |CICP27      |
|ENSG00000222623.1 |ENSG00000222623 |RNU6-1100P  |


``` r
gse1_total_unique <- gse1_mapped %>%
	distinct(ensembl_gene_id) %>%
	nrow()

gse1_mapped_unique <- gse1_mapped %>%
	filter(!is.na(hgnc_symbol), hgnc_symbol != "") %>%
	distinct(ensembl_gene_id) %>%
	nrow()

tibble(
	metric = c("Unique stripped IDs", "Mapped to HGNC", "Unmapped", "Mapped proportion"),
	n = c(
		gse1_total_unique,
		gse1_mapped_unique,
		gse1_total_unique - gse1_mapped_unique,
		round(gse1_mapped_unique / gse1_total_unique, 4)
	)
) %>%
	knitr::kable(caption = "GSE119732 mapping summary")
```



Table: (\#tab:ex1_summary)GSE119732 mapping summary

|metric              |          n|
|:-------------------|----------:|
|Unique stripped IDs | 57992.0000|
|Mapped to HGNC      | 40666.0000|
|Unmapped            | 17326.0000|
|Mapped proportion   |     0.7012|

## Exercise 2 - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122380

Using **GSE122380**, confirm whether the ID column contains Ensembl IDs with version suffixes.


``` r
gse2_file <- find_data_file("GSE122380")
```

```
## Using locally cached version of supplementary file(s) GSE122380 found here:
## data/GSE122380/GSE122380_Supplementary_Data_Table_S1.xlsx
```

```
## Using locally cached version of supplementary file(s) GSE122380 found here:
## data/GSE122380/GSE122380_raw_counts.txt.gz
```

``` r
gse2 <- safe_read(gse2_file)
```

```
## Warning: One or more parsing issues, call `problems()` on your data frame for details,
## e.g.:
##   dat <- vroom(...)
##   problems(dat)
```

``` r
gse2_id_col <- names(gse2)[1]
gse2_ids <- as.character(gse2[[1]])

gse2_first20 <- tibble(raw_id = head(gse2_ids, 20))
gse2_dot_count <- sum(grepl("\\.", gse2_first20$raw_id))

gse2_with_stripped <- gse2 %>%
	mutate(ensembl_gene_id = strip_ensembl_version(.data[[gse2_id_col]]))

gse2_map <- map_to_hgnc(gse2_ids, ensembl, cache_key = "GSE122380")
```

```
## BioMart unavailable; using cached HGNC mapping from cache/hgnc_map_GSE122380.rds
```

``` r
gse2_mapped <- gse2_with_stripped %>%
	left_join(gse2_map, by = "ensembl_gene_id")
```

1. First 20 IDs:


``` r
knitr::kable(gse2_first20)
```



|raw_id          |
|:---------------|
|ENSG00000000419 |
|ENSG00000000457 |
|ENSG00000000460 |
|ENSG00000000938 |
|ENSG00000000971 |
|ENSG00000001036 |
|ENSG00000001084 |
|ENSG00000001167 |
|ENSG00000001460 |
|ENSG00000001461 |
|ENSG00000001561 |
|ENSG00000001617 |
|ENSG00000001626 |
|ENSG00000001629 |
|ENSG00000001630 |
|ENSG00000001631 |
|ENSG00000002016 |
|ENSG00000002330 |
|ENSG00000002549 |
|ENSG00000002587 |

2. New column with versions stripped (preview):


``` r
knitr::kable(head(gse2_with_stripped[, c(gse2_id_col, "ensembl_gene_id")], 10))
```



|Gene_id         |ensembl_gene_id |
|:---------------|:---------------|
|ENSG00000000419 |ENSG00000000419 |
|ENSG00000000457 |ENSG00000000457 |
|ENSG00000000460 |ENSG00000000460 |
|ENSG00000000938 |ENSG00000000938 |
|ENSG00000000971 |ENSG00000000971 |
|ENSG00000001036 |ENSG00000001036 |
|ENSG00000001084 |ENSG00000001084 |
|ENSG00000001167 |ENSG00000001167 |
|ENSG00000001460 |ENSG00000001460 |
|ENSG00000001461 |ENSG00000001461 |

3. Mapping to HGNC symbols (preview):


``` r
gse2_mapped_nonempty <- gse2_mapped %>%
	filter(!is.na(hgnc_symbol), hgnc_symbol != "") %>%
	distinct(ensembl_gene_id, .keep_all = TRUE)

if (nrow(gse2_mapped_nonempty) > 0) {
	knitr::kable(
		head(gse2_mapped_nonempty[, c(gse2_id_col, "ensembl_gene_id", "hgnc_symbol")], 10),
		caption = "GSE122380: first 10 successfully mapped IDs"
	)
} else {
	knitr::kable(
		head(gse2_mapped[, c(gse2_id_col, "ensembl_gene_id", "hgnc_symbol")], 10),
		caption = "GSE122380: no HGNC mappings returned in this run"
	)
}
```



Table: (\#tab:ex2_map_preview)GSE122380: first 10 successfully mapped IDs

|Gene_id         |ensembl_gene_id |hgnc_symbol |
|:---------------|:---------------|:-----------|
|ENSG00000000419 |ENSG00000000419 |DPM1        |
|ENSG00000000457 |ENSG00000000457 |SCYL3       |
|ENSG00000000460 |ENSG00000000460 |FIRRM       |
|ENSG00000000938 |ENSG00000000938 |FGR         |
|ENSG00000000971 |ENSG00000000971 |CFH         |
|ENSG00000001036 |ENSG00000001036 |FUCA2       |
|ENSG00000001084 |ENSG00000001084 |GCLC        |
|ENSG00000001167 |ENSG00000001167 |NFYA        |
|ENSG00000001460 |ENSG00000001460 |STPG1       |
|ENSG00000001461 |ENSG00000001461 |NIPAL3      |


``` r
gse2_total_unique <- gse2_mapped %>%
	distinct(ensembl_gene_id) %>%
	nrow()

gse2_mapped_unique <- gse2_mapped %>%
	filter(!is.na(hgnc_symbol), hgnc_symbol != "") %>%
	distinct(ensembl_gene_id) %>%
	nrow()

tibble(
	metric = c("Unique stripped IDs", "Mapped to HGNC", "Unmapped", "Mapped proportion"),
	n = c(
		gse2_total_unique,
		gse2_mapped_unique,
		gse2_total_unique - gse2_mapped_unique,
		round(gse2_mapped_unique / gse2_total_unique, 4)
	)
) %>%
	knitr::kable(caption = "GSE122380 mapping summary")
```



Table: (\#tab:ex2_summary)GSE122380 mapping summary

|metric              |          n|
|:-------------------|----------:|
|Unique stripped IDs | 16319.0000|
|Mapped to HGNC      | 15854.0000|
|Unmapped            |   465.0000|
|Mapped proportion   |     0.9715|

4. What is different about this file?


``` r
gse2_diff <- tibble(
	metric = c(
		"First-20 IDs containing '.'",
		"Proportion of all IDs starting with ENSG",
		"Duplicate stripped IDs"
	),
	value = c(
		gse2_dot_count,
		mean(grepl("^ENSG", strip_ensembl_version(gse2_ids))),
		sum(duplicated(strip_ensembl_version(gse2_ids)))
	)
)

knitr::kable(gse2_diff, caption = "GSE122380 quick structure checks")
```



Table: (\#tab:ex2_difference)GSE122380 quick structure checks

|metric                                   | value|
|:----------------------------------------|-----:|
|First-20 IDs containing '.'              |     0|
|Proportion of all IDs starting with ENSG |     1|
|Duplicate stripped IDs                   |     0|

In this dataset, the first-column identifier structure can be compared directly with GSE119732 using the table above.
The practical difference to watch for is whether all IDs are gene-level Ensembl IDs (`ENSG...`) and whether version suffixes/deduplication affect mapping rates.

## Exercise 3

Can you use the worked example to process the above two GEO records? How?

Yes. The worked example generalizes cleanly by turning its core steps into reusable functions:

1. Download supplementary files for a chosen GSE.
2. Load the first expression table safely.
3. Identify/clean first-column IDs by stripping Ensembl version suffixes.
4. Map stripped IDs to HGNC symbols with `biomaRt`.
5. Join mappings back to expression data and summarize mapped vs unmapped IDs.

That is exactly what the setup and analysis chunks above do for both **GSE119732** and **GSE122380**.
