#' Generate Protein-Protein Interaction (PPI) Network Report
#'
#' @description
#' This function processes input protein data, maps protein IDs using the provided local files,
#' generates PPI networks (both at protein and gene levels), saves intermediate files and network plots,
#' and returns a list containing the mapping result and network data.
#'
#' @param dt A data.frame containing protein identifiers. The column specified by \code{proteincol.id} is used.
#' @param proteincol.id A character string indicating the column name in \code{dt} that contains the protein IDs.
#' Default is "Accessions".
#' @param sep.pattern A regular expression used to separate multiple protein IDs in the column defined by \code{proteincol.id}.
#' Default is ";\\s*".
#' @param score_threshold An integer specifying the minimum combined score for filtering PPI interactions. Default is 400.
#' @param savepath A character string specifying the directory where the output files (report and plots) will be saved.
#' Default is "report".
#' @param protein_aliases_path A character string specifying the local file path for the protein aliases file.
#' This parameter is required.
#' @param protein_Info_path A character string specifying the local file path for the protein information file.
#' This parameter is required.
#' @param protein_link_path A character string specifying the local file path for the protein link file.
#' This parameter is required.
#'
#' @return A list with the following components:
#' \item{id_map}{A data.table containing the mapping result (with protein alias and additional information).}
#' \item{ppi_protein_data}{A data.table of protein-level PPI interactions with columns renamed to \code{From} and \code{To}.}
#' \item{ppi_gene_data}{A data.table of gene-level PPI interactions with columns renamed to \code{From} and \code{To}.}
#' \item{plt_ppi_protein}{A ggraph object representing the protein-level PPI network plot.}
#' \item{plt_ppi_gene}{A ggraph object representing the gene-level PPI network plot.}
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Checks whether \code{protein_aliases_path}, \code{protein_Info_path} and \code{protein_link_path} are provided. 
#'         If any of these are missing, it throws an error with an appropriate message.
#'   \item Reads the local files using \code{data.table::fread} to obtain protein aliases, protein information, and PPI links.
#'   \item Processes the input data (\code{dt}) by separating multiple protein IDs based on the provided delimiter.
#'   \item Maps protein IDs to obtain corresponding aliases and merges with protein information (keeping only the first occurrence for duplicates).
#'   \item Filters the PPI network based on the score threshold and replaces string IDs with protein aliases or gene names.
#'   \item Rescales the combined scores and generates network plots using \code{igraph} and \code{ggraph}.
#'   \item Saves the intermediate mapping table and interaction data, as well as the generated plots, to the specified output directory.
#' }
#'
#' @examples
#' \dontrun{
#' demo <- data.frame(
#'   Accessions = c("P01042", "P07911", "P02760", "Q08380", "O14686")
#' )
#'
#' result <- easyPPI(
#'   dt = demo,
#'   proteincol.id = "Accessions",
#'   sep.pattern = ";\\s*",
#'   score_threshold = 400,
#'   savepath = "report",
#'   protein_aliases_path = "9606.protein.aliases.v12.0.txt.gz",
#'   protein_Info_path = "9606.protein.info.v12.0.txt.gz",
#'   protein_link_path = "9606.protein.links.v12.0.txt.gz"
#' )
#' }
#'
#' @importFrom data.table fread as.data.table fwrite
#' @importFrom tidyr separate
#' @importFrom dplyr select rename any_of
#' @importFrom scales rescale
#' @importFrom igraph graph_from_data_frame
#' @importFrom ggraph ggraph geom_edge_link scale_edge_color_gradient geom_node_point geom_node_text
#' @importFrom ggplot2 aes theme_void ggtitle theme element_text
#' @importFrom rlang sym
#' @importFrom readr read_csv
#'
#' @export
#' 
easyPPI=function(dt,
        proteincol.id="Accessions",
        sep.pattern=";\\s*",
        score_threshold=400,
        savepath="report",
        protein_aliases_path,# protein_aliases_path="9606.protein.aliases.v12.0.txt.gz"
        protein_Info_path,# protein_Info_path="9606.protein.info.v12.0.txt.gz"
        protein_link_path# protein_link_path="9606.protein.links.v12.0.txt.gz"
        ){

  if(missing(protein_aliases_path) || missing(protein_Info_path) || missing(protein_link_path)){
    stop("\nYou must provide the local file paths for `protein.aliases`, `protein.Info`, and `protein.link`--zcp.")
  }
  
  protein.aliases=data.table::fread(protein_aliases_path)  
  protein.Info=data.table::fread(protein_Info_path)
  protein.link=data.table::fread(protein_link_path)

  
  

#proteincol.id="Accessions"#p1
 #p2:分隔符处理 定义分隔符符号
 #sep.pattern=";\\s*"#p2 拆分蛋白列

#p3 score_threshold :PPI互作网络卡值 默认400
 #score_threshold=400
 
#p4 图片和数据保存路径
 #savepath="report"
 
 
 
 #路径生成
 .create_dir =function(path) {
   if (!dir.exists(path)) {
     dir.create(path, recursive = TRUE)
   }
 }
 .create_dir(savepath)
 
 
dt=dt |>
  tidyr::separate( col = !!rlang::sym(proteincol.id), into = c(proteincol.id, "delcol"), sep =sep.pattern) |>
  dplyr::select(!delcol)

#功能一:id map 获取backgroup和hits的gene name -uniprot id-  string id
map_id=function(input_hits,lib){
  mappedid=(data.table::as.data.table(lib))[alias%in%input_hits]
  mappedid=mappedid[, .SD[1], by = alias]
  missing_hits = length(setdiff(input_hits, lib$alias)) 
  message(sprintf("%.2f%% of IDs mapping failed.--zcp", missing_hits*100/length(input_hits) ))
  return(mappedid)
}
mppping_id=map_id(input_hits=unique(dt[[proteincol.id]]),lib=protein.aliases)
 #合并基因信息
protein_info_unique = protein.Info[!duplicated(`#string_protein_id`)]

mppping_id= merge(mppping_id, protein_info_unique, by = "#string_protein_id", all.x = TRUE)# 返回1 id映射结果----

mppping_id |> dplyr::select(!dplyr::any_of(c("#string_protein_id","source"))) |> 
  data.table::fwrite(paste0(savepath,"/mppping_id.tsv"),sep = "\t")

#功能二:PPI网络
#linklib=protein.link
#hit_id=mppping_id$`#string_protein_id`
#score=score_threshold
PPI_interaction=function(hit_id,linklib,score){
Interaction=linklib[
  ((linklib$protein1) %in% hit_id)&( (linklib$protein2) %in% hit_id)
                     ]
Interaction=Interaction[combined_score >= score]
return(Interaction)
}

ppi_netdat.pro=PPI_interaction(hit_id =mppping_id$`#string_protein_id`,linklib=protein.link,score=score_threshold )

temp_mapping = setNames(mppping_id$alias, mppping_id$`#string_protein_id`)
#把stringid变成uniprotid
ppi_netdat.pro[, protein1 := temp_mapping[protein1]]
ppi_netdat.pro[, protein2 := temp_mapping[protein2]]
ppi_netdat.pro[, scaled_score := scales::rescale(combined_score, to = c(0.5, 3))]# 返回值2 ppi_pro数据----

ppi_netdat.pro |> dplyr::rename(From=protein1,To=protein2) |> 
data.table::fwrite(paste0(savepath,"/PPI_protein_interaction.tsv"),sep = "\t")


ppi_network_pro <- igraph::graph_from_data_frame(d = ppi_netdat.pro, directed = T)
plt.ppi.protein=
  ggraph::ggraph(ppi_network_pro, layout = "fr") +  # fr circle tree linear matrix dendrogram reingold.tilford
  # ggraph::geom_edge_link(ggplot2::aes(width=scaled_score/3,color = scaled_score),alpha=0.3) + 
  ggraph::geom_edge_link(ggplot2::aes(alpha=scaled_score/3,color = scaled_score),width=0.4,linetype="solid") + 
  ggraph::scale_edge_color_gradient(
    low = "#B3D0EB",  # 低值颜色
    high ="#00ABF0",        #高值颜色
    guide = "colourbar"
  )+
  ggraph::geom_node_point(
    shape = 21,     
    fill =   "#25848E",   # 填充颜色（节点内部）
    color = "white" ,      # 边框颜色
    size = 2,             # 节点大小
    stroke = 0.3            # 边框宽度
  ) +            # 节点样式
  ggraph::geom_node_text(ggplot2::aes(label = name), 
                         check_overlap = TRUE,vjust = 1.5, size = 2) +  # 节点标签
  ggplot2::theme_void() +                                             # 去掉背景和坐标轴
  ggplot2::ggtitle("PPI Network(protein)")+
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, size = 10, face = "bold")  # 标题居中
  )# 返回值 3 ppi.protein图----

save_zcp(Fig = plt.ppi.protein,FigName = "PPI_network_protein",outputfile =savepath,widths = 4,heights =4,ppt = TRUE)

#绘制基因的ppi

ppi_netdat_gene=PPI_interaction(hit_id =mppping_id$`#string_protein_id`,linklib=protein.link,score=score_threshold )

temp_mapping = setNames(mppping_id$preferred_name, mppping_id$`#string_protein_id`)
#把stringid变成gene name
ppi_netdat_gene[, protein1 := temp_mapping[protein1]]
ppi_netdat_gene[, protein2 := temp_mapping[protein2]]
ppi_netdat_gene[, scaled_score := scales::rescale(combined_score, to = c(0.5, 3))]# 返回值4 ppi_gene数据----

ppi_netdat_gene |> dplyr::rename(From=protein1,To=protein2) |> 
  data.table::fwrite(paste0(savepath,"/PPI_gene_interaction.tsv"),sep = "\t")

ppi_network_gene <- igraph::graph_from_data_frame(d = ppi_netdat_gene, directed = T)
plt.ppi.gene=
  ggraph::ggraph(ppi_network_gene, layout = "fr") +  # fr circle tree linear matrix dendrogram reingold.tilford
  # ggraph::geom_edge_link(ggplot2::aes(width=scaled_score/3,color = scaled_score),alpha=0.3) + 
  ggraph::geom_edge_link(ggplot2::aes(alpha=scaled_score/3,color = scaled_score),width=0.4,linetype="solid") + 
  ggraph::scale_edge_color_gradient(
    low = "#B3D0EB", 
    high ="#00ABF0",      
    guide = "colourbar"
  )+
  ggraph::geom_node_point(
    shape = 21,     
    fill =   "#25848E",  
    color = "white" ,     
    size = 2,             
    stroke = 0.3          
  ) +           
  ggraph::geom_node_text(ggplot2::aes(label = name), 
                         check_overlap = TRUE,vjust = 1.5, size = 2) +  
  ggplot2::theme_void() +                                             
  ggplot2::ggtitle("PPI Network(gene)")+
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, size = 10, face = "bold")  
  )# 返回值 5 ppi.gene图----
save_zcp(Fig = plt.ppi.gene,FigName = "PPI_network_gene",outputfile =savepath,widths = 4,heights =4,ppt = TRUE)

return(list(
  id_map=mppping_id |> dplyr::select(!dplyr::any_of(c("#string_protein_id","source"))) ,
  ppi_protein_data=ppi_netdat.pro |> dplyr::rename(From=protein1,To=protein2),
  ppi_gene_data=ppi_netdat_gene |> dplyr::rename(From=protein1,To=protein2),
  plt_ppi_gene=plt.ppi.gene,
  plt_ppi_protein=plt.ppi.protein
))

}

