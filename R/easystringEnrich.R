#' easystringEnrich: Protein Enrichment Analysis Function
#'
#' @description
#' Performs enrichment analysis on the input protein data. This function integrates protein aliases, protein information,
#' and protein enrichment terms data (string database v12.0), calculates p-values and FDR using the hypergeometric test, and generates pathway
#' enrichment plots.
#'
#' @param dt A data table (or data frame) containing protein ID information.
#' @param proteincol.id A string specifying the column name that stores protein IDs. Default is "Accessions".
#' @param sep.pattern A string specifying the delimiter used to separate protein IDs. Default is ";".
#' @param savepath A string specifying the directory where the results and plots will be saved. Default is "enrichoutputfile".
#' @param backgrounds Optional background protein data table; if NULL, all proteins in the background will be used.
#' @param cutoff Numeric value; proteins with counts below this threshold will be filtered out. Default is 2.
#' @param protein_aliases_path A string providing the local file path for the protein aliases file
#'   (e.g., "9606.protein.aliases.v12.0.txt.gz").
#' @param protein_Info_path A string providing the local file path for the protein info file
#'   (e.g., "9606.protein.info.v12.0.txt.gz").
#' @param protein_enrichment_terms_path A string providing the local file path for the protein enrichment terms file
#'   (e.g., "9606.protein.enrichment.terms.v12.0.txt.gz").
#'
#' @return A list containing:
#' \describe{
#'   \item{enrichment}{The enrichment analysis results as a tibble.}
#'   \item{gene_pathway_matrix}{The gene-pathway mapping matrix as a tibble.}
#'   \item{enrichplots}{A named list of ggplot objects, each corresponding to a pathway category.}
#' }
#'
#' @importFrom data.table fread uniqueN setnames fwrite as.data.table
#' @importFrom dplyr mutate select filter arrange slice_head
#' @importFrom tidyr separate
#' @importFrom purrr map
#' @importFrom ggplot2 ggplot geom_rect geom_segment geom_point scale_size labs theme_bw
#' @importFrom ggplot2 scale_fill_gradient scale_color_gradient element_text unit
#' @importFrom stringr str_wrap
#' @importFrom tibble as_tibble
#' @importFrom stats phyper p.adjust
#'
#' @export
#'
easystringEnrich=function(
  dt,
  proteincol.id="Accessions",
  sep.pattern=";",
  savepath="enrichoutputfile",
  backgrounds=NULL, # backgrounds=readr::read_csv("enrich_demo_background.csv")
  cutoff=2,
  protein_aliases_path,# protein_aliases_path="9606.protein.aliases.v12.0.txt.gz"
  protein_Info_path,# protein_Info_path="9606.protein.info.v12.0.txt.gz"
  protein_enrichment_terms_path # protein_enrichment_terms_path="9606.protein.enrichment.terms.v12.0.txt.gz"

){

  if(missing(protein_aliases_path) || missing(protein_Info_path) || missing(protein_enrichment_terms_path)){
    stop("\nYou must provide the local file paths for `protein.aliases`, `protein.Info`, and `protein.link`--zcp.")
  }

  protein.aliases=data.table::fread(protein_aliases_path)
  protein.Info=data.table::fread(protein_Info_path)
  protein.enrichment=data.table::fread(protein_enrichment_terms_path)


  .create_dir =function(path) {
    if (!dir.exists(path)) {
      dir.create(path, recursive = TRUE)
    }
  }
  .create_dir(savepath)


  map_id=function(input_hits,lib){
    mappedid=(data.table::as.data.table(lib))[alias%in%input_hits]
    mappedid=mappedid[, .SD[1], by = alias]
    missing_hits = length(setdiff(input_hits, lib$alias))
    message(sprintf("%.2f%% of IDs mapping failed.--zcp", missing_hits*100/length(input_hits) ))
    return(mappedid)
  }

  #差异蛋白id信息整合



  dt= dt |>#dplyr::mutate(!!rlang::sym(proteincol.id)=purrr::map_chr(stringr::str_split(proteincol.id, sep.pattern), ~ .x[1]))
    tidyr::separate( col = !!rlang::sym(proteincol.id), into = c(proteincol.id, "delcol"), sep =sep.pattern) |>
    dplyr::select(!delcol)
  mppping_id=map_id(input_hits=unique(dt[[proteincol.id]]),lib=protein.aliases)
  protein_info_unique = protein.Info[!duplicated(`#string_protein_id`)]
  mppping_id= merge(mppping_id, protein_info_unique, by = "#string_protein_id", all.x = TRUE)
  hits=unique(mppping_id$`#string_protein_id`)
 #背景蛋白 如果有就进行Id整合 没有就默认全部的背景
  if (is.null(backgrounds)) {
    background=protein.enrichment
  }else{
    bkhits=backgrounds |> tidyr::separate( col = !!rlang::sym(proteincol.id), into = c(proteincol.id, "delcol"), sep =sep.pattern) |>
      dplyr::select(!delcol)
    #背景蛋白的stringid
    mppping_id_bk=map_id(input_hits=unique(bkhits[[proteincol.id]]),lib=protein.aliases)
    #按照背景蛋白的stringid将Protein.enrichment通路中的数据找出来
    background=protein.enrichment[(`#string_protein_id`)%in%(mppping_id_bk[["#string_protein_id"]])]
  }




  #统计背景数据库中全部蛋白的数量
  total_protein=background[, .(unique_id_count = data.table::uniqueN(`#string_protein_id`)), by = category]

  #统计背景数据库中每个pathway蛋白的数量:number_of_genes_in_background
  proteincounts.per.pahtway.all=background[, .(number_of_genes_in_background = .N), by = .(category, term)]


  #统计差异蛋白在每个通路中的数量
  diff.protein.database.clean=background[`#string_protein_id` %in% hits]
  proteincounts.per.pahtway.diff=diff.protein.database.clean[, .(number_of_genes = .N), by = .(category, term)]

  #数据整合合并
  proteincounts.merged = data.table::as.data.table(proteincounts.per.pahtway.all)[
    data.table::as.data.table(proteincounts.per.pahtway.diff),
    on = c("category", "term"),
    nomatch = 0L,  #只保留common的Pathway值
    allow.cartesian = FALSE #避免重复Pathway
  ]

  #增加每个通路背景蛋白的数量
  proteincounts.merged=merge(
    proteincounts.merged,
    total_protein,
    by = "category",
    all.x = TRUE
  )
  #改名
  data.table::setnames(proteincounts.merged, "unique_id_count", "background.proteincounts")

  #增加差异蛋白数量
  proteincounts.merged$diff.proteincounts=length(unique(hits))
  #超几何分布
  result=data.table::as.data.table(proteincounts.merged)[, p_value :=  stats::phyper(
    q = number_of_genes - 1, #每条通路中的差异蛋白数量-1
    m = number_of_genes_in_background, #通路中的背景蛋白数
    n = background.proteincounts - number_of_genes_in_background, #非通路背景蛋白数
    k = diff.proteincounts, #差异蛋白总数
    lower.tail = FALSE
  )]
  result=data.table::as.data.table(result)[, FDR := p.adjust(p_value, method = "BH")]
  result=data.table::as.data.table(result)[, Enrich_factor := round((number_of_genes / number_of_genes_in_background) /
                                                                      (diff.proteincounts / background.proteincounts),3)]



  #计算 Strength：取 observed/expected 的对数（以 10 为底）
  result[, Strength := log10(Enrich_factor)]
  result=result[Strength > 0.2  ]
  result=result[FDR <=0.05  ]
  w1 <- 1
  w2 <- 2
  result[, Signal := (w1 + w2) / ((w1 / Strength) + (w2 / (-log10(FDR))))]


  #先给差异蛋白信息富集结果增加额外的注释信息

  protein_info_subset = protein.Info[, c("#string_protein_id", "preferred_name")]
  diff.protein.database.clean_string = merge(
    diff.protein.database.clean,
    protein_info_subset,
    by = "#string_protein_id",
    all.x = TRUE
  )

  #补充gene信息
  genes=data.table::as.data.table(diff.protein.database.clean_string)[, .(preferredNames = paste(preferred_name, collapse = "/")), by = c("category","term")]

  #补充description
  description=unique(background[, setdiff(names(background), "#string_protein_id"), with = FALSE], by = c("category", "term", "description"))



  #补充蛋白信息
  #未完待续
  Final.results.geneInfo= data.table::as.data.table(result)[
    data.table::as.data.table(genes),
    on = c("category","term"),
    nomatch = 0L,  #只保留common的Pathway值
    allow.cartesian = FALSE #避免重复Pathway
  ] |>    dplyr::filter(number_of_genes>=cutoff)

  Final.results= data.table::as.data.table(Final.results.geneInfo)[
    data.table::as.data.table(description),
    on = c("category","term"),
    nomatch = 0L,  #只保留common的Pathway值
    allow.cartesian = FALSE #避免重复Pathway
  ]


  Final.results |> data.table::fwrite(paste0(savepath,"/Enrichment.tsv"),sep = "\t")
  diff.protein.database.clean_string |> dplyr::select(!`#string_protein_id`)|> data.table::fwrite(paste0(savepath,"/gene_pathway_matrix.tsv"),sep = "\t")


  categories =unique(Final.results$category)


  plots = purrr::map(categories, function(cat_val) {

    cat("plotting:", cat_val, "...--zcp\n")


    dt.plot=Final.results |> dplyr::filter(category==cat_val)  |>
      dplyr::slice_max(order_by = number_of_genes, n = 10)   |>
      dplyr::mutate(description = stringr::str_wrap(description, width = 30))

    p =
      ggplot2::ggplot() +
      ggplot2::geom_rect(data = dt.plot, ggplot2::aes(xmin = -1,
                                                      xmax = as.numeric(max(.data[["number_of_genes"]]))*1.1,
                                                      ymin = as.numeric(reorder(description, number_of_genes)) - 0.5, #背景色
                                                      ymax = as.numeric(reorder(description, number_of_genes)) + 0.5,
                                                      fill = FDR),
                         alpha = 0.3,
                         color="white",
                         linewidth=1,
                         show.legend = F,
                         inherit.aes = FALSE
      ) +
      ggplot2::geom_segment(data = dt.plot,
                            ggplot2::aes(x = 0,
                                         xend = number_of_genes,
                                         color= FDR,
                                         y = reorder(description,number_of_genes),
                                         yend = reorder(description,number_of_genes)),
                            linewidth = 1,
                            show.legend = F
      ) +
      ggplot2::geom_point(data=dt.plot,
                          ggplot2::aes(x = number_of_genes,
                                       y = reorder(description, number_of_genes),
                                       #color = fdr,
                                       fill= FDR,
                                       size = number_of_genes),
                          shape=21,
                          color="#2FA4C2"
      ) +
      ggplot2::scale_size(range = c(2, 4)) +
      ggplot2::labs(title = cat_val,
                    subtitle = "Top 10 term by counts",
                    x = "counts",
                    y = NULL,
                    fill = "FDR",
                    color="FDR",
                    size = "Protein counts") +
      ggplot2::theme_bw() +
      ggplot2::scale_fill_gradient(high = "#D6EFB3", low = "#2FA4C2") +
      ggplot2:: scale_color_gradient(high = "#D6EFB3", low = "#2FA4C2") +

      ggplot2::theme(axis.text.y =  ggplot2::element_text(size = 8),
                     plot.title =  ggplot2::element_text(hjust = 0.5,size = 10),
                     plot.subtitle =  ggplot2::element_text(hjust = 0.5,size = 7),
                     panel.grid.major = ggplot2:: element_blank(),#去掉主网格线
                     panel.grid.minor = ggplot2:: element_blank(),#去掉次网格线
                     legend.position = "right",
                     #缩小图例
                     legend.key.size =ggplot2:: unit(0.4, "cm"),#图例符号尺寸
                     legend.text =  ggplot2::element_text(size = 8),#图例文字大小
                     legend.title = ggplot2:: element_text(size = 8), #图例标题大小
                     legend.spacing.x =ggplot2:: unit(0.2, "cm"), #图例项之间的水平间距
                     legend.spacing.y =ggplot2:: unit(0.2, "cm") #图例项之间的垂直间距
      )



    file_name = gsub("\\s|\\(|\\)", "", cat_val)
    save_zcp(Fig = p, FigName = file_name, outputfile = savepath, widths = 6, heights = 4)
    return(p)
  })

  names(plots)=gsub("\\s|\\(|\\)", "", unique(Final.results$category))




  return(list(
    enrichment=tibble::as_tibble(Final.results),
    gene_pathway_matrix=tibble::as_tibble(
      diff.protein.database.clean_string |>
        dplyr::select(!`#string_protein_id`)),
    enrichplots=plots
              ))




}
