#################  ADDITIONAL SOURCED SCRIPTS #################################
source("./tissue_annotations.R")
source("./cell_annotations.R")


#################  JSON PARSING AND ANNOTATION READING ########################

# Extracts annotation data from MetaAnnotator JSON object
extract_annotation <- function(json) {
  tibble(
    BioSample = map_chr(json, "biosampleAccession"),
    attributes = map(json, "attributes")) %>%
    mutate(
      attributeName = map(
        attributes, function(x) pluck(x, 1,"attributeValue")
      ),
      attributeName_annotated = map(
        attributes, function(x) pluck(x, 1,"attributeValueTermLabel")
      )
    ) %>%
    select(attributeName, attributeName_annotated) %>%
    mutate_all(as.character)
}

# Converts extracted JSON metadata into a datafram
create_annotationDF_fromDir <- function(jsonPath){
  tibble(file = list.files(jsonPath)) %>%
    mutate(filepath = paste0(jsonPath, file)) %>%
    rowwise() %>%
    mutate(json = list(fromJSON(file = filepath))) %>%
    ungroup() %>%
    mutate(json = map_depth(json, 1, extract_annotation)) %>%
    pull(json) %>%
    reduce(rbind)
}

# Recommends a disease annotation based on (1) MONDO, then (2) MEDDRA
disease_recommender <- function(mondo, meddra){
  ifelse(!is.na(mondo), mondo,
         ifelse(!is.na(meddra), meddra, 
                NA)
  )
}

#################  TISSUE RECOMMENDATION FUNCTIONS ############################

# Tissue recommender
tissue_recommender <- function(bto, uberon, clo, mondo, meddra){
  if (bto %in% tissue.annotations){
    bto
  } else if (uberon %in% tissue.annotations){
    uberon
  } else if (clo %in% tissue.annotations){
    clo
  } else if (mondo %in% tissue.annotations){
    mondo
  } else if (meddra %in% tissue.annotations){
    meddra
  } 
  else
    NA
}

#################  CELL RECOMMENDATION FUNCTIONS ##############################

# Cell recommender
cell_recommender <- function(bto, uberon, clo, mondo, meddra){
  if (bto %in% cell.annotations){
    bto
  } else if (uberon %in% cell.annotations){
    uberon
  } else if (clo %in% cell.annotations){
    clo
  } else if (mondo %in% cell.annotations){
    mondo
  } else if (meddra %in% cell.annotations){
    meddra
  } 
  else
    NA
}

make_immune_lineage_key <-  function() {
  # Primary patterns
  tcell <- c("t-cell", "tcell", " t cell", "^t cell", "t-lymphocyte", "t lymphocyte",
             "^cd4$", " th[12]", "^th[12]", "cd4+ ", "cd4t", "cd8t", "helper", "cd4\\+",
             "cd8[ _psa\\+]", "cd8$", "cd3\\+", "cd3 p", "treg", "regulatory[ _]t",
             "regulatory cell", "teff", "thymocyte", "mait", "cd4-[pt]",
             "[( ^]th0[ )]", "^th0$", "temra", "tfh") %>% paste(.,collapse = "|")
  bcell <- c("b-cell", "bcell", " b cell", "^b cell", "b lympho", "b memory", 
             "plasmablast", "plasma cell", "cd19\\+", "cd19 ", 
             "cd19-b") %>% paste(.,collapse = "|")
  ilc <- c("innate lymph", "^ilc", " ilc") %>% paste(.,collapse = "|")
  nk <- c("nk cell", "nkcell", "natural killer", "cd56[db\\+]") %>% paste(.,collapse = "|")
  dc <- c("dendritic.cell", "-cdc", " cdc", "plasmacytoid", "^pdc", 
          "[ _(]pdc", "myeloid dc", "+ dc", "predc", "langerhans cells", 
          "modc", "cd11c\\+", "dentritic") %>% paste(.,collapse = "|")
  macs <- c("cd11b[ \\+]", "macrop") %>% paste(.,collapse = "|")
  mono <- c("monocyte") %>% paste(.,collapse = "|")
  pmn <- c("^pmn", "neutroph") %>% paste(.,collapse = "|")
  baso <- c("basophils", "^basophil$") %>% paste(.,collapse = "|")
  eos <- c("eosinophil ", "eosinophils", "eosinophil$") %>% paste(.,collapse = "|")
  mast <- c("mast cell")
  
  #Key assembly
  immune_key <- as_tibble(matrix(c(
    "t-lymphocyte", tcell,
    "b-lymphocyte", bcell,
    "innate lymphoid cell", ilc,
    "natural killer cell", nk,
    "dendritic cell", dc,
    "macrophage", macs,
    "monocyte", mono,
    "neutrophil", pmn,
    "basophil", baso,
    "eosinophil", eos,
    "mast cell", mast
  ), ncol = 2, byrow = T, dimnames = list(NULL ,c("id", "pattern"))))
  
  immune_key
}

immune_lineage_match <- function(input, key){
  key %>%
    mutate(match = map2_chr(pattern, id, 
                            function(pattern, id) ifelse(grepl(pattern, input), id, NA))) %>%
    filter(!is.na(match)) %>%
    pull(match) %>%
    sort() %>% 
    paste(., collapse = "_") %>%
    ifelse(.=="", NA, .)
}

make_immune_subset_key <- function(){
  #T-cell subtypes
  t_treg <- c("treg", "regulatory[ _]t", "t regu") %>% paste(.,collapse = "|")
  t_cd4 <- c("cd4[ \\+t_]", "cd4$", "cd4\\-[pt]", "help") %>% paste(.,collapse = "|")
  t_cd8 <- c("cd8[^n\\-]", "cd8$", "cd8\\-[pt]", "cytotoxic") %>% paste(.,collapse = "|")
  t_tfh <- c("follicular", "tfh") %>% paste(.,collapse = "|")
  th0 <- c("th0", "0 cell") %>% paste(.,collapse = "|")
  th1 <- c("th1$", "th1[a-z _\\-]", "1 cell") %>% paste(.,collapse = "|")
  th2 <- c("th2", "2 cell") %>% paste(.,collapse = "|")
  th17 <- c("th17", "17 cell")  %>% paste(.,collapse = "|")
  t_eff <- c("effector", "tem$", "tem ", "teff") %>% paste(.,collapse = "|")
  t_central <- c("central", "tcm") %>% paste(.,collapse = "|")
  t_mem <- c("memory",  "tem$", "tem ", "tcm", "cd45ro[ \\+]") %>% paste(.,collapse = "|")
  t_naive <- c("na[iÃ¯]ve", "cd45ra[ \\+]") %>% paste(.,collapse = "|")
  t_mait <- c("mait", "invariant") %>% paste(.,collapse = "|")
  t_inkt <- c("inkt")
  t_gd <- c("gamma[ \\-]*delta($| )", "gd t", "Î³Î´") %>% paste(.,collapse = "|")
  t_temra <- c("temra")
  t_prog <- c("prec","prog", "cd34[ \\+]") %>% paste(.,collapse = "|")
  
  tcell_key <- as_tibble(matrix(c(
    "cd4+ regulatory", t_treg,
    "cd4+", t_cd4,
    "cd8+", t_cd8,
    "follicular", t_tfh,
    "th0 cd4+", th0,
    "th1 cd4+", th1,
    "th2 cd4+", th2,
    "th17 cd4+", th17,
    "effector", t_eff,
    "central", t_central,
    "memory", t_mem,
    "naive", t_naive,
    "mait", t_mait,
    "inkt", t_inkt,
    "gamma-delta", t_gd,
    "temra memory", t_temra,
    "progenitor", t_prog
  ), ncol = 2, byrow = T, dimnames = list(NULL ,c("id", "pattern")))) %>% 
    mutate(id = paste0(id, " t-lymphocyte")) %>% 
    mutate(group = "t-lymphocyte") 
  
  #B-cell subtypes
  b_mem <- c("memory",  "cd27\\+") %>% paste(.,collapse = "|")
  b_plasma <- c("plasma ")
  b_plasmablast <- c("plasmablast")
  b_naive <- c("na[iÃ¯]ve", "nave") %>% paste(.,collapse = "|")
  b_germ <- c("germ")
  b_fol <- c("follic")
  b_marg <- c("marginal")
  b_igg <- c("igg[ \\+]")
  b_igd <- c("igd[ \\+]")
  b_igm <- c("igm[h \\+]")
  b_iga <- c("iga[ \\+]")
  b_prog <- c("prec","prog", "cd34\\+", "cd34 (^n|p)") %>% paste(.,collapse = "|")
  
  bcell_key <- as_tibble(matrix(c(
    "naive", b_naive,
    "memory", b_mem,
    "plasma cell", b_plasma,
    "plasmablast", b_plasmablast,
    "germinal center", b_germ,
    "follicular", b_fol,
    "marginal zone", b_marg,
    "igg+", b_igg,
    "igd+", b_igd,
    "igm+", b_igm,
    "iga+", b_iga,
    "progenitor", b_prog
  ), ncol = 2, byrow = T, dimnames = list(NULL ,c("id", "pattern")))) %>% 
    mutate(id = paste0(id, " b-lymphocyte"))%>% 
    mutate(group = "b-lymphocyte") 
  
  #ILC subtypes
  ilc_1 <- c("group 1", "ilc1", " 1$") %>% paste(.,collapse = "|")
  ilc_2 <- c("group 2", "ilc2") %>% paste(.,collapse = "|")
  ilc_3 <- c("group 2", "ilc3", " 3$") %>% paste(.,collapse = "|")
  
  ilc_key <- as_tibble(matrix(c(
    "group 1", ilc_1,
    "group 2", ilc_2,
    "group 3", ilc_3
  ), ncol = 2, byrow = T, dimnames = list(NULL ,c("id", "pattern")))) %>% 
    mutate(id = paste0(id, " innate lymphoid cell"))%>% 
    mutate(group = "innate lymphoid cell") 
  
  #DC subtypes
  dc_pdc <- c("plasmacytoid", "pdc") %>% paste(.,collapse = "|")
  dc_prog <- c("pre[dc-]","prog") %>% paste(.,collapse = "|")
  dc_mono <- c("monocyte", "modc") %>% paste(.,collapse = "|")
  dc_cdc <- c("conv", "cdc", "classic") %>% paste(.,collapse = "|")
  dc_cd141 <- c("cd141")
  dc_cd1c <- c("cd1c")
  dc_bdca1 <- c("bdca1")
  dc_lang <- c("langerhans cells", "lang", "epiderm") %>% paste(.,collapse = "|")
  dc_my <- c("myeloid")
  
  dc_key <- as_tibble(matrix(c(
    "plasmacytoid", dc_pdc,
    "progenitor", dc_prog,
    "monocyte-derived", dc_mono,
    "langerhans", dc_lang,
    "conventional", dc_cdc,
    "cd141+ conventional", dc_cd141,
    "cd1c+ conventional", dc_cd1c,
    "bdca1+ conventional", dc_bdca1,
    "myeloid", dc_my
  ), ncol = 2, byrow = T, dimnames = list(NULL ,c("id", "pattern")))) %>% 
    mutate(id = paste0(id, " dendritic cell")) %>% 
    mutate(group = "dendritic cell") 
  
  #NK subtypes
  nk_56b <- c("cd56high", "cd56\\+", "cd56bright") %>% paste(.,collapse = "|")
  nk_56d <- c("cd56dim", "cd56\\-") %>% paste(.,collapse = "|")
  nk_16p <- c("cd16\\+") %>% paste(.,collapse = "|")
  nk_16n <- c("cd16\\-") %>% paste(.,collapse = "|")
  nk_inv <- c("invariant")
  
  nk_key <- as_tibble(matrix(c(
    "cd56bright", nk_56b,
    "cd56dim", nk_56d,
    "cd16+", nk_16p,
    "cd16-", nk_16n,
    "invarient", nk_inv
  ), ncol = 2, byrow = T, dimnames = list(NULL ,c("id", "pattern")))) %>% 
    mutate(id = paste0(id, " natural killer cell"))%>% 
    mutate(group = "natural killer cell") 
  
  #monocyte subtypes
  mono_classic <- c("^class", " class") %>% paste(.,collapse = "|")
  mono_nonclassic <- c("nonclass", "non class", "non\\-class") %>% paste(.,collapse = "|")
  mono_cd14p <- c("cd14\\+", "cd14[ p]") %>% paste(.,collapse = "|")
  mono_cd14d <- c("cd14\\-[c ]", "cd14dim") %>% paste(.,collapse = "|")
  mono_cd16p <- c("cd16\\+", "cd16\\-p", "cd16p") %>% paste(.,collapse = "|")
  mono_cd16n <- c("cd16\\-[ n]", "cd16n") %>% paste(.,collapse = "|")
  
  mono_key <- as_tibble(matrix(c(
    "classical", mono_classic,
    "non-classical", mono_nonclassic,
    "cd14+", mono_cd14p,
    "cd14dim", mono_cd14d,
    "cd16+", mono_cd16p,
    "cd16-", mono_cd16n
  ), ncol = 2, byrow = T, dimnames = list(NULL ,c("id", "pattern")))) %>% 
    mutate(id = paste0(id, " monocyte"))%>% 
    mutate(group = "monocyte")
  
  #macrophage subtypes
  mac_m1 <- c("m1", "classic") %>% paste(.,collapse = "|")
  mac_m2 <- c("m2", "alternative") %>% paste(.,collapse = "|")
  mac_m0 <- c("m0") %>% paste(.,collapse = "|")
  mac_prog <- c("prec","prog") %>% paste(.,collapse = "|")
  
  mac_key <- as_tibble(matrix(c(
    "m1", mac_m1,
    "m2", mac_m2,
    "m0", mac_m0,
    "progenitor", mac_prog
  ), ncol = 2, byrow = T, dimnames = list(NULL ,c("id", "pattern")))) %>% 
    mutate(id = paste0(id, " macrophage"))%>% 
    mutate(group = "macrophage") 
  
  #Compile
  subset_key <- rbind(tcell_key, bcell_key, ilc_key, dc_key, mono_key, nk_key, mac_key)
  subset_key
}

make_immune_name_order <- function(){
  t_order <- list(
    c("progenitor",  
      "temra", "naive", "effector", "central", "memory", 
      "th0", "th1", "th2", "th17", "follicular", "regulatory", 
      "cd4+", "cd8+", 
      "mait", "inkt", "gamma-delta",
      "t-lymphocyte"))
  
  b_order <- list(
    c("progenitor",  
      "igd+", "igm+", "igg+", "iga+",
      "naive", "memory", 
      "follicular", "marginal", "zone", "germinal", "center",
      "plasma", "cell", "plasmablast", 
      "b-lymphocyte"))
  
  ilc_order <- list(c("group", "1", "2", "3", "innate", "lymphoid", "cell"))
  
  dc_order <- list(
    c("progenitor",  
      "cd1c+", "cd141+", "bdca1+", 
      "plasmacytoid", "monocyte-derived", "langerhans", "conventional", "myeloid",
      "dendritic", "cell"))
  
  nk_order <- list(
    c("cd56bright", "cd56dim", "cd16+", "cd16-",
      "invariant",
      "natural", "killer", "cell"))
  
  mono_order <- list(
    c("cd14dim", "cd14+", "cd16+", "cd16-",
      "classical", "non-classical",
      "monocyte"))
  
  mac_order <- list(
    c("progenitor",  
      "m0", "m1", "m2",
      "macrophage"))
  
  order_key <- tibble(
    "t-lymphocyte" = t_order,
    "b-lymphocyte" = b_order,
    "innate lymphoid cell" = ilc_order,
    "dendritic cell" = dc_order,
    "natural killer cell" = nk_order,
    "monocyte" = mono_order,
    "macrophage" = mac_order,
    index = 1) %>% 
    pivot_longer(-index, names_to = "id", values_to = "order") %>% 
    select(-index)
  order_key
}

immune_subset_match <- function(input_value, group_value, pattern_key, order_key){
  
  order <- order_key %>% 
    filter(id == group_value) %>% 
    pull(order) %>% unlist()
  
  pattern_match <- pattern_key %>%
    mutate(match = map2_chr(pattern, id, 
                            function(pattern, id) ifelse(grepl(pattern, input_value), id, NA))) %>%
    filter(!is.na(match) & group == group_value) %>%
    pull(match) 
  
  pattern_match <- unique(
    unlist(str_split(pattern_match, " "))
  ) 
  
  sort(factor(pattern_match, levels = order, ordered=TRUE)) %>%
    paste(., collapse = " ")
}  

stem_pattern <- c(
  "stem[ -]cell",
  "embryonic stem",
  "ietic stem",
  "chymal stem",
  "hsc(s|$)",
  "(^|[ h_-])esc",
  "stem-like",
  "ipsc",
  "es[ -]cell",
  "[bhi]msc"
) %>% paste(.,collapse = "|")

make_stem_subset_key <-  function() {
  # Primary patterns
  esc <- c("(^|[h _-])esc", "embryonic", "es[ -]cell") %>% paste(.,collapse = "|")
  msc <- c("mesenchym", "[bh]msc") %>% paste(.,collapse = "|")
  ipsc <- c("ipsc", "induced") %>% paste(.,collapse = "|")
  hsc <- c("hsc", "hema", "haem", "peripheral blood") %>% paste(.,collapse = "|")
  csc <- c("cancer", "cml", "leukem") %>% paste(.,collapse = "|")
  cbsc <- c("cord", "jelly") %>% paste(.,collapse = "|")
  musc_sc <- c("musc") %>% paste(.,collapse = "|")
  card_sc <- c("cardiac") %>% paste(.,collapse = "|")
  gbm_sc <- c("gbm", "glio") %>% paste(.,collapse = "|")
  bsc <- c("biliary") %>% paste(.,collapse = "|")
  nsc <- c("neural", "neuronal") %>% paste(.,collapse = "|")
  epsc <- c("epidermal") %>% paste(.,collapse = "|")
  hpsc <- c("hepatic")%>% paste(.,collapse = "|")
  prsc <- c("prosta") %>% paste(.,collapse = "|")
  
  #Key assembly
  stem_key <- as_tibble(matrix(c(
    "embryonic stem cell", esc,
    "mesenchymal stem cell", msc,
    "induced pluripotent stem cell", ipsc,
    "hematopoietic stem cell", hsc,
    "cancer stem cell", csc,
    "cord blood stem cell", cbsc,
    "muscle stem cell", musc_sc,
    "cardiac stem cell", card_sc,
    "glioblastoma stem cell", gbm_sc,
    "biliary stem cell", bsc,
    "neuronal stem cell", nsc,
    "epidermal stem cell", epsc,
    "hepatic stem cell", hpsc,
    "prostate stem cell", prsc
  ), ncol = 2, byrow = T, dimnames = list(NULL ,c("id", "pattern"))))
  
  stem_key
}

stem_subset_match <- function(input, key){
  key %>%
    mutate(match = map2_chr(pattern, id, 
                            function(pattern, id) ifelse(grepl(pattern, input), id, NA))) %>%
    filter(!is.na(match)) %>%
    pull(match) %>%
    sort() %>% 
    paste(., collapse = "_") %>%
    ifelse(.=="", NA, .)
}

#################  SEX RECOMMENDATION FUNCTIONS  ##############################

make_sex_key <- function(){
  female <- c("^female","[ _]female", "^f$") %>% paste(., collapse = "|")
  male <- c("^male","[ _]male[ _,;]", "[_ \\-]male$", "^m$") %>% paste(., collapse = "|")
  sex_key <- tibble("id" = c("male", "female"), "pattern" = c(male, female))
  sex_key
}

sex_match <- function(input, key){
  key %>%
    mutate(match = map2_chr(pattern, id,
                            function(pattern, id) ifelse(grepl(pattern, input), id, NA))) %>%
    filter(!is.na(match)) %>%
    pull(match) %>%
    sort() %>%
    paste(., collapse = "_") %>%
    ifelse(.=="", NA, .)
}

#################  AGE RECOMMENDATION FUNCTIONS  ##############################

# Age unit finder
unit.finder <- function(x){
  age.categories <- 
    c("adolescent", "adult", "adult (unknown)", "adult liver", "birth", 
      "blastocyst", "blastocyst; icm and polar trophectoderm", "blastocyst; mural trophectoderm", 
      "child", "children", "days", "early embryo", "early pregnancy", 
      "elderly", "emboryonic", "embroyonic kidney", "embryo", "embryo stage", 
      "embryonic", "embryos", "fetal", "fetal liver", "fetal lung", 
      "fetal progenitor organoid", "fetal stage", "fetus", "first trimeste", 
      "foetus", "from third trimester placenta after delivery", "full term", 
      "high", "hours", "human adult", "human embryo", "infant", "juvenile", 
      "juvenile stage", "late", "late embryo", "late embryo stage", 
      "late embryonic stage", "late-stage", "long-lived", "long-lived subjects", 
      "low", "mature", "mid", "mid-age", "mid-stage", "middle adult", 
      "middle-aged", "mixed", "mixed sample (non-demultiplexed)", "morula", 
      "neonatal", "neonatal, reprogrammed, and differentiated", "neonatal,reprogrammed", 
      "neonate", "new born", "newborn", "old", "old adult", "older", 
      "oocyte", "pediatric", "pool", "pooled", "post-partum", "postimplantation embryo", 
      "preimplantation-stage embryo", "prenatal", "school age child", 
      "term delivery", "term placenta", "three monthes", "toddler", 
      "trophoblast", "twelve weeks before pregnancy", "young", "young adult", 
      "young control", "zygote")

  unit_key <- as_tibble(matrix(c(
    "year", "year|yr|[0-9 ]y$|y/o| yo|y//.o//| y.|[0-9]+y",
    "month","month|[0-9 ]+m$|mo$|mos" ,
    "week", "week|wk|[0-9 ]w$|gw",
    "day", "day|[0-9]d$",
    "pooled", "[0-9] and [0-9]",
    "unspecified", " to |[0-9]{1,2}\\-[0-9]|^[0-9]{1,2}$|^[0-9]{1,2}\\.[0-9]"), 
    ncol = 2, byrow = T, dimnames = list(NULL ,c("unit", "pattern"))))
  
  if (x %in% age.categories){
    "category"
  } else {
    #Finding unit using matching patters in unit_key
    unit_key %>% 
      mutate(match = 
               map2_chr(pattern, unit, 
                        function(pattern, unit) ifelse(grepl(pattern, x), unit, NA))
      ) %>%
      # Collapsing unit_key matches into a string
      filter(!is.na(match)) %>% 
      pull(match) %>% 
      sort(.) %>% 
      paste(., collapse = "_") %>% 
      ifelse(.=="", NA, .) %>% 
      # Removing _unspecified from values with any other unit match
      ifelse(
        ((str_count(.,"_") >= 1) & grepl("_unspecified", .)),
        str_remove(., "_unspecified"), .
      )
  }
}

# Time cleaning function
time.cleaner <- function(unit, value){
  value <- value %>% 
    gsub(" to ", "-", .) %>%  # convert "to" to "-"
    gsub("less|under", "<", .) %>%  # convert "less" to "<"
    gsub("older", ">", .) # convert "older" to ">"
  age.value <- value %>% 
    str_remove_all(. , "[a-zA-Z ]+") %>%  # delete characters and spaces
    str_remove_all(. , "^[+~:\\-]") %>% #remove all beginning "+", "-", ":", and "~"
    str_remove_all(. , "\\(\\)") %>% #remove all "()"
    str_remove_all(. , "[\\-]+$|\\/$") %>% # remove all trailing "-" and "/"
    str_remove_all(. , "^,") %>% # remove all starting ","
    ifelse( grepl(">", .), # rename all values with "older than" meaning 
            paste0( str_remove_all(. , "[^0-9]"), " and older"), 
            . ) %>% 
    ifelse( grepl("<", .), # rename all values with "younger than" meaning
            paste0( str_remove_all(. , "[^0-9]"), " and younger"),
            . ) %>% 
    # Only keep above changes for the following units
    ifelse(unit %in% c("year", "month", "week", "day"),
           . ,
           value ) %>% 
    # Converting day_week values to weeks only
    ifelse( unit == "day_week",
            value %>% 
              map_dbl(. , 
                      function(y) unlist(str_split(y, "[ a-z\\-\\,]+")) %>%
                        .[. != ""] %>% as.numeric() %>%
                        (function(x) round(x[1] + (x[2]/7), 2))),
            . ) %>% 
    # Converting year_month values to years only
    ifelse( unit == "month_year",
            value %>% 
              map_dbl(. , 
                      function(y) unlist(str_split(y, "[ a-z\\-\\,]+")) %>%
                        .[. != ""] %>% as.numeric() %>%
                        (function(x) round(x[1] + (x[2]/12), 2))),
            . )
  # Relabeling compound units
  unit <- unit %>% 
    ifelse(. == "day_week", "week", .) %>% 
    ifelse(. == "month_year", "year", .)
  
  # Return age value
  paste(value, unit, age.value, sep = "_____")
}


#################  SINGLE CELL IDENTIFICATION RESOURCES ######################
sc_pattern <- c("single[ -_]cell", "sc[_-]?rna", "single[ -_]nucle[ui]") %>% paste0(.,collapse = "|")
not_sc_pattern <- c("^no$", "^none$", "^false$") %>% paste0(.,collapse = "|")


#################  FINAL CLEANING FUNCTIONS ###################################

# Used for final age cleaning to return a pattern if it is found in a combo value
# that results from pivoting (e.g. year_____category --> year)
return_if_match_in_string <- function(string, pattern){
  result <- unlist(
    map(
      strsplit(string, "_____"), 
      function(x) pattern %in% x))
  return(ifelse(result, pattern, string))
}

# test <- c("category_____year", "month")
# return_if_match_in_string(test, "year")

# Patterns to seach final columns for cell line like names
cell_line_pattern <- c("cell[ _]line", "(^|[ _-])line([ _-]|\\.ORIG$)") %>% paste(., collapse = "|")



#################  ARCHIVE FUNCTIONS ##########################################

# make_immune_key <-  function() {
#   # Primary patterns
#   tcell <- c("t-cell", "tcell", " t cell", "^t cell", "t-lymphocyte", "t lymphocyte",
#              "^cd4$", " th[12]", "^th[12]", "cd4+ ", "cd4t", "cd8t", "helper", "cd4\\+",
#              "cd8[ _psa\\+]", "cd8$", "cd3\\+", "cd3 p", "treg", "regulatory[ _]t",
#              "regulatory cell", "teff", "thymocyte", "mait", "cd4-[pt]",
#              "[( ^]th0[ )]", "^th0$", "temra", "tfh") %>% paste(.,collapse = "|")
#   bcell <- c("b-cell", "bcell", " b cell", "^b cell", "b lympho", "b memory", 
#              "plasmablast", "plasma cell", "cd19\\+", "cd19 ", 
#              "cd19-b") %>% paste(.,collapse = "|")
#   ilc <- c("innate lymph", "^ilc", " ilc") %>% paste(.,collapse = "|")
#   nk <- c("nk cell", "nkcell", "natural killer", "cd56[db\\+]") %>% paste(.,collapse = "|")
#   dc <- c("dendritic.cell", "-cdc", " cdc", "plasmacytoid", "^pdc", 
#           "[ _(]pdc", "myeloid dc", "+ dc", "predc", "langerhans cells", 
#           "modc", "cd11c\\+", "dentritic") %>% paste(.,collapse = "|")
#   macs <- c("cd11b[ \\+]", "macrop") %>% paste(.,collapse = "|")
#   mono <- c("monocyte") %>% paste(.,collapse = "|")
#   pmn <- c("^pmn", "neutroph") %>% paste(.,collapse = "|")
#   baso <- c("basophils", "^basophil$") %>% paste(.,collapse = "|")
#   eos <- c("eosinophil ", "eosinophils", "eosinophil$") %>% paste(.,collapse = "|")
#   mast <- c("mast cell")
#   
#   #Key assembly
#   immune_key <- as_tibble(matrix(c(
#     "t-lymphocyte", tcell,
#     "b-lymphocyte", bcell,
#     "innate lymphoid cell", ilc,
#     "natural killer cell", nk,
#     "dendritic cell", dc,
#     "macrophage", macs,
#     "monocyte", mono,
#     "neutrophil", pmn,
#     "basophil", baso,
#     "eosinophil", eos,
#     "mast cell", mast
#   ), ncol = 2, byrow = T, dimnames = list(NULL ,c("id", "pattern"))))
#   
#   immune_key
# }
# 
# 
# 
# make_immune_subset_key <- function(){
#   
#   # immune_key <- make_immune_key()
#   
#   #T-cell subtypes
#   t_treg <- c("treg", "regulatory[ _]t", "t regu") %>% paste(.,collapse = "|")
#   t_cd4 <- c("^cd4$", " th[12]", "^th[12]", "cd4+ ", "cd4t") %>% paste(.,collapse = "|")
#   t_cd8 <- c("cd8t") %>% paste(.,collapse = "|")
#   t_tfh <- c("follicular") %>% paste(.,collapse = "|")
#   
#   #DC subtypes
#   dc_pdc <- c("plasmacytoid", "^pdc", "[ _(]pdc") %>% paste(.,collapse = "|")
#   dc_prog <- c("pre[d-]","prog") %>% paste(.,collapse = "|")
#   dc_mono <- c("monocyte[ -]") %>% paste(.,collapse = "|")
#   dc_cdc <- c("conv", "[( ]cdc") %>% paste(.,collapse = "|")
#   dc_lang <- c("langerhans cells") %>% paste(.,collapse = "|")
#   
#   #B-cell subtypes
#   b_plasma <- c("plasmablast", "plasma cell") %>% paste(.,collapse = "|")
#   
#   tcell_key <- as_tibble(matrix(c(
#     "regulatory t-lymphocyte", t_treg,
#     "cd4+ t-lymphocyte", t_cd4,
#     "cd8+ t-lymphocyte", t_cd8,
#     "follicular t-lymphocyte", t_tfh
#   ), ncol = 2, byrow = T, dimnames = list(NULL ,c("subtype", "pattern")))) %>% mutate(id = "t-lymphocyte")
#   
#   dc_key <- as_tibble(matrix(c(
#     "plasmacytoid dendritic cell", dc_pdc,
#     "dendritic cell precursor", dc_prog,
#     "monocyte-derived dendritic cell", dc_mono,
#     "conventional denritic cell", dc_cdc,
#     "langerhans cell", dc_lang
#   ), ncol = 2, byrow = T, dimnames = list(NULL ,c("subtype", "pattern")))) %>% mutate(id = "dendritic cell")
#   
#   bcell_key <- as_tibble(matrix(c(
#     "plasma cell", b_plasma
#   ), ncol = 2, byrow = T, dimnames = list(NULL ,c("subtype", "pattern")))) %>% mutate(id = "b-lymphocyte")
#   
#   key <- immune_key %>% left_join(bind_rows(
#     tcell_key,
#     dc_key
#   ), by = "id")
# }
# 
# 
# 
# 
# 
# pattern_match <- function(input, key){
#   key %>%
#     mutate(match = map2_chr(pattern, id, 
#                             function(pattern, id) ifelse(grepl(pattern, input), id, NA))) %>%
#     filter(!is.na(match)) %>%
#     pull(match) %>%
#     sort() %>% 
#     paste(., collapse = "_") %>%
#     ifelse(.=="", NA, .)
# }

# pattern_match <-  function(x) {
#   # Primary patterns
#   tcell <- c("t-cell", "tcell", " t cell", "^t cell", "t-lymphocyte", "t lymphocyte",
#              "^cd4$", " th[12]", "^th[12]", "cd4+ ", "cd4t", "cd8t", "helper", "cd4\\+",
#              "cd8[ _psa\\+]", "cd8$", "cd3\\+", "cd3 p", "treg", "regulatory[ _]t",
#              "regulatory cell", "teff") %>% paste(.,collapse = "|")
#   bcell <- c("b-cell", "bcell", " b cell", "^b cell", "b lympho", "b memory", 
#              "plasmablast", "plasma cell", "cd19\\+", "cd19 l", "cd19-b") %>% paste(.,collapse = "|")
#   ilc <- c("innate lymphoid", "^ilc", " ilc") %>% paste(.,collapse = "|")
#   nk <- c("nk cell", "nkcell", "natural killer", "cd56[db\\+]") %>% paste(.,collapse = "|")
#   dc <- c("dendritic.cell", "-cdc", " cdc", "plasmacytoid", "^pdc", "[ _(]pdc", "myeloid dc", 
#           "+ dc", "predc", "langerhans cells") %>% paste(.,collapse = "|")
#   macs <- c("cd11b[ \\+]", "macrop") %>% paste(.,collapse = "|")
#   pmn <- c("^pmn", "neutroph") %>% paste(.,collapse = "|")
#   baso <- c("^basophil[s]?")
#   eos <- c("eosinophil ", "eosinophils", "eosinophil$") %>% paste(.,collapse = "|")
#   mast <- c("mast cell")
# 
#   #T-cell subtypes
#   t_treg <- c("treg", "regulatory[ _]t", "t regu") %>% paste(.,collapse = "|")
#   t_cd4 <- c("^cd4$", " th[12]", "^th[12]", "cd4+ ", "cd4t") %>% paste(.,collapse = "|")
#   t_cd8 <- c("cd8t") %>% paste(.,collapse = "|")
#   t_tfh <- c("follicular") %>% paste(.,collapse = "|")
# 
#   #DC subtypes
#   dc_pdc <- c("plasmacytoid", "^pdc", "[ _(]pdc") %>% paste(.,collapse = "|")
#   dc_prog <- c("pre[d-]","prog") %>% paste(.,collapse = "|")
#   dc_mono <- c("monocyte[ -]") %>% paste(.,collapse = "|")
#   dc_cdc <- c("conv", "[( ]cdc") %>% paste(.,collapse = "|")
#   dc_lang <- c("langerhans cells") %>% paste(.,collapse = "|")
# 
#   #B-cell subtypes
#   b_plasma <- c("plasmablast", "plasma cell") %>% paste(.,collapse = "|")
# 
#   #Key assembly
#   immune_key <- as_tibble(matrix(c(
#     "t-lymphocyte", tcell,
#     "b-lymphocyte", bcell,
#     "innate lymphoid cell", ilc,
#     "natural killer cell", nk,
#     "dendritic cell", dc,
#     "macrophage", macs,
#     "neutrophil", pmn,
#     "basophil", baso,
#     "eosinophil", eos,
#     "mast cell", mast
#   ), ncol = 2, byrow = T, dimnames = list(NULL ,c("cell", "pattern"))))
#   
#   tcell_key <- as_tibble(matrix(c(
#     "regulatory t-lymphocyte", t_treg,
#     "cd4+ t-lymphocyte", t_cd4,
#     "cd8+ t-lymphocyte", t_cd8,
#     "follicular t-lymphocyte", t_tfh
#   ), ncol = 2, byrow = T, dimnames = list(NULL ,c("subtype", "pattern")))) %>% mutate(cell = "t-lymphocyte")
# 
#   dc_key <- as_tibble(matrix(c(
#     "plasmacytoid dendritic cell", dc_pdc,
#     "dendritic cell precursor", dc_prog,
#     "monocyte-derived dendritic cell", dc_mono,
#     "conventional denritic cell", dc_cdc,
#     "langerhans cell", dc_lang
#   ), ncol = 2, byrow = T, dimnames = list(NULL ,c("subtype", "pattern")))) %>% mutate(cell = "dendritic cell")
# 
#   bcell_key <- as_tibble(matrix(c(
#     "plasma cell", b_plasma
#   ), ncol = 2, byrow = T, dimnames = list(NULL ,c("subtype", "pattern")))) %>% mutate(cell = "b-lymphocyte")
# 
#   key <- immune_key %>% left_join(bind_rows(
#     tcell_key,
#     dc_key
#   ), by = "cell")
#   
#   
#   first_pass <- immune_key$pattern %>% paste(., collapse = "|")
#   
#   tibble(value = x) %>% 
#     mutate(match = grepl(first_pass, x)) %>% # first pass filter for speed
#     mutate(recommended = ifelse(
#       match == TRUE,
#       # Matches an input value to pattern in immune_key and pulls the key value out
#       map_chr(value, function(y) {
#       immune_key %>%
#         mutate(match = map2_chr(pattern, cell, function(pattern, cell) ifelse(grepl(pattern, y), cell, NA))) %>%
#         filter(!is.na(match)) %>%
#         pull(match) %>%
#         paste(., collapse = "_") %>%
#         ifelse(.=="", NA, .)
#       }),
#       NA
#     )) %>% 
#     pull(recommended)
#   
# }


# immune_key %>%
#   mutate(match = map2_chr(pattern, cell, function(pattern, cell) ifelse(grepl(pattern, x), cell, NA))) %>%
#   filter(!is.na(match)) %>%
#   pull(match) %>%
#   paste(., collapse = "_") %>%
#   ifelse(.=="", NA, .),




# first_pass <- "t.?cell|lympho|cd4|th[12]|cd8|helper|b.?cell|memory|plasma|ilc|nk|natural killer|dendritic|dc|langerhans"
# first_pass <- c(tcell, bcell, ilc, nk, dc, macs, pmn, baso, eos, mast) %>% paste(., collapse = "|")

# # Old immune cell recommender
# immuneCell_finder <- function(value){
#   tcell_pattern <- c("t-cell", "tcell", " t cell", "^t cell", "t-lymphocyte", "t lymphocyte", 
#                      "^cd4$", " th[12]", "^th[12]", "cd4+ ", "cd4t", "cd8t")
#   tcell_pattern <- paste(tcell_pattern, collapse = "|")
#   
#   bcell_pattern <- c("b-cell", "bcell", " b cell", "^b cell", "b lympho", "b memory")
#   bcell_pattern <- paste(bcell_pattern, collapse = "|")
#   
#   plasma_pattern <- c("plasmablast", "plasma cell")
#   plasma_pattern <- paste(plasma_pattern, collapse = "|")
#   
#   ilc_pattern <- c("innate lymphoid", "^ilc", " ilc")
#   ilc_pattern <- paste(ilc_pattern, collapse = "|")
#   
#   nk_pattern <- c("nk cell", "nkcell", "natural killer")
#   nk_pattern <- paste(nk_pattern, collapse = "|")
#   
#   dc_pattern <- c("dendritic.cell", "-cdc", " cdc", "plasmacytoid", "^pdc", "[ _(]pdc", "myeloid dc", "+ dc", "predc")
#   dc_pattern <- paste(dc_pattern, collapse = "|")
#   
#   
#   if (grepl(tcell_pattern, value)){
#     "t-lymphocyte"
#   } else if (grepl(bcell_pattern, value)){
#     "b-lymphocyte"
#   } else if (grepl(plasma_pattern, value)){
#     "plasma cell"
#   } else if (grepl(ilc_pattern, value)){
#     "innate lymphoid cell"
#   } else if (grepl(nk_pattern, value)){
#     "natural killer cell"
#   } else if (grepl(dc_pattern, value)){
#     "dendritic cell"
#   } else if (grepl("macrophage", value)){
#     "macrophage"
#   } else if (grepl("monocyte", value)){
#     "monocyte"
#   } else if (grepl("granulocyte", value)){
#     "granulocyte"
#   } else if (grepl("neutrophil", value)){
#     "neutrophil"
#   } else if (grepl("mast cell", value)){
#     "mast cell"
#   } else {
#     NA
#   }
# }

# pattern_match <-  function(value) {
#   # Patterns
#   tcell <- c("t-cell", "tcell", " t cell", "^t cell", "t-lymphocyte", "t lymphocyte", 
#              "^cd4$", " th[12]", "^th[12]", "cd4+ ", "cd4t", "cd8t", "helper") %>% paste(.,collapse = "|")
#   bcell <- c("b-cell", "bcell", " b cell", "^b cell", "b lympho", "b memory", "plasmablast", "plasma cell") %>% paste(.,collapse = "|")
#   ilc <- c("innate lymphoid", "^ilc", " ilc") %>% paste(.,collapse = "|")
#   nk <- c("nk cell", "nkcell", "natural killer") %>% paste(.,collapse = "|")
#   # dc <- c("dendritic.cell", "-cdc", " cdc", "plasmacytoid", "^pdc", "[ _(]pdc", "myeloid dc", "+ dc", "predc", "langerhans cells") %>% paste(.,collapse = "|")
#   dc <- c("dendritic.cell", "-cdc", " cdc", "plasmacytoid", "^pdc", "[ _(]pdc", "myeloid dc",  "predc", "langerhans cells") %>% paste(.,collapse = "|")
#   
#   
#   #T-cell subtypes
#   t_treg <- c("treg", "regulatory[ _]t", "t regu") %>% paste(.,collapse = "|")
#   t_cd4 <- c("^cd4$", " th[12]", "^th[12]", "cd4+ ", "cd4t") %>% paste(.,collapse = "|")
#   t_cd8 <- c("cd8t") %>% paste(.,collapse = "|")
#   t_tfh <- c("follicular") %>% paste(.,collapse = "|")
#   
#   #DC subtypes
#   dc_pdc <- c("plasmacytoid", "^pdc", "[ _(]pdc") %>% paste(.,collapse = "|")
#   dc_prog <- c("pre[d-]","prog") %>% paste(.,collapse = "|")
#   dc_mono <- c("monocyte[ -]") %>% paste(.,collapse = "|")
#   dc_cdc <- c("conv", "[( ]cdc") %>% paste(.,collapse = "|")
#   dc_lang <- c("langerhans cells") %>% paste(.,collapse = "|")
#   
#   #B-cell subtypes
#   b_plasma <- c("plasmablast", "plasma cell") %>% paste(.,collapse = "|")
#   
#   
#   if (grepl(tcell, value)){
#     "t-lymphocyte"
#   } else if (grepl(bcell, value)){
#     "b-lymphocyte"
#   } else if (grepl(ilc, value)){
#     "innate lymphoid cell"
#   } else if (grepl(nk, value)){
#     "natural killer cell"
#   } else if (grepl(dc, value)){
#     "dendritic cell"
#   } else if (grepl("macrophage", value)){
#     "macrophage"
#   } else if (grepl("monocyte", value)){
#     "monocyte"
#   } else if (grepl("granulocyte", value)){
#     "granulocyte"
#   } else if (grepl("neutrophil", value)){
#     "neutrophil"
#   } else if (grepl("mast cell", value)){
#     "mast cell"
#   } else {
#     NA
#   }
# }
# 
# 
