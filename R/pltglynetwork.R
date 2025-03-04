#' Plot Glycan Network
#'
#' This function processes glycoprotein data to extract glycan types and glycosylation site information, then visualizes the glycan network using ggplot2 and ggnetwork. The network comprises nodes representing glycan compositions and proteins, with edges indicating their associations.
#'
#' @param dt A data.frame containing the glycoprotein data.
#' @param peptide.col Character. The column name in \code{dt} representing the peptide sequence. Default is "Annotated Sequence".
#' @param glycancomposition.col Character. The column name in \code{dt} representing the glycan composition. Default is "Glycan composition".
#' @param protein.col Character. The column name in \code{dt} representing the protein identifier. Default is "Master Protein Accessions".
#' @param outdir Character. Directory path to save output files. Default is "outputfile_test".
#' @param SiteperProteinsCutoff Numeric. A cutoff value for the number of glycosylation sites per protein. Default is 4.
#' @param edge_cex Numeric. The line width for the network edges in the plot. Default is 0.5.
#' @param edge_alpha Numeric, range from 0 to 1.The transparency for the network edges in the plot. Default is 0.3.
#' @param node_alpha Numeric, range from 0 to 1.The transparency for the network nodes in the plot. Default is 0.5
#' @param node_size Numeric. The size of the network nodes. Default is 3.
#' @param Title Character. The title of the network plot. Default is "test".
#'
#' @return A ggplot object representing the glycan network plot.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Selects and renames the required columns.
#'   \item Computes the number of glycosylation sites per protein.
#'   \item Processes the glycan composition by replacing specific substrings and extracting numeric values.
#'   \item Classifies glycan types based on computed counts.
#'   \item Constructs nodes and edges for the glycan network and computes a custom layout.
#'   \item Saves intermediate data files and outputs the final network plot.
#' }
#'
#' @examples
#' \dontrun{
#'   # Example usage:
#'   library(ZCPstates)
#'   data("Glycan_network_demo", package = "ZCPstates")
#'   plt <- pltglynetwork(Glycan_network_demo,
#'                        peptide.col = "Annotated Sequence",
#'                        glycancomposition.col = "Glycan composition",
#'                        protein.col = "Master Protein Accessions",
#'                        Title = "Glycan Network")
#'   print(plt)
#' }
#'
#' @export
#'
#' @import ggplot2
#' @import ggnetwork
#' @import igraph
#' @importFrom data.table fwrite
#' @importFrom dplyr select group_by mutate case_when summarise n arrange desc
#' @importFrom stringr str_replace str_extract
#' @importFrom tidyr extract
#' @importFrom rlang sym
#' @importFrom ggsci pal_npg
#' @importFrom useful cart2pol pol2cart
#' @importFrom vctrs vec_c

pltglynetwork=function(dt,
                       peptide.col="Annotated Sequence",
                       glycancomposition.col="Glycan composition",
                       protein.col="Master Protein Accessions",
                       outdir="outputfile_test",
                       SiteperProteinsCutoff=4,
                       edge_alpha=0.3,
                       node_alpha=0.5,
                       edge_cex=0.5,
                       node_size=3,
                       Title="test"
                       ){

  #使用sym() 和!!来动态引用列名
  peptide = rlang::sym(peptide.col)
  glycancomposition = rlang::sym(glycancomposition.col)
  protein=rlang::sym(protein.col)

  #按照顺序提取列名
  dat=dt |> dplyr::select(!!peptide,!!glycancomposition,!!protein)
  #更改列名
  colnames(dat)=c("peptide","glycancomposition","protein")


  #每个蛋白的糖基化修饰位点个数
  dat_sites_per_peotein <- dat |>
    dplyr::group_by(protein)  |>
    dplyr:: mutate(Glyco_site=length(unique(peptide)))#有n种肽段 Glyco_site就等于n



  #糖型分类
  #第一步---
  #把HexNAc替换为A
  tst <-stringr::str_replace( dat_sites_per_peotein$glycancomposition,"HexNAc","A")
  #把Hex替换为B
  tst1 <- stringr::str_replace(tst,"Hex","B")
  #把NeuAc替换成C
  tst2 <- stringr::str_replace(tst1,"NeuAc","C")
  #把Fuc替换成F
  tst3 <- stringr::str_replace(tst2,"Fuc","F")
  dat_sites_per_peotein$tst3 <- tst3
  #提取出每个糖型后的数字，从字符问题转化为数据问题
  #HexNAc
  HexNAc <- stringr::str_extract(dat_sites_per_peotein$tst3,"(A\\(\\d{1,3}\\))") |>
    as.data.frame()  |>  rlang::set_names("HexNAc")
  HexNAc <- tidyr::extract(HexNAc,HexNAc,"A(","(\\d{1,3})") |>  as.data.frame() |>
    rlang::set_names("HexNAc")#HexNAc个数
  #Hex
  Hex <- stringr::str_extract(dat_sites_per_peotein$tst3,"(B\\(\\d{1,3}\\))") |>
    as.data.frame()  |> rlang::set_names("Hex")
  Hex <- tidyr::extract(Hex,Hex,"B(","(\\d{1,3})")  |>  as.data.frame()  |>
    rlang::set_names("Hex")#Hex个数
  #NeuAc
  NeuAc <- stringr::str_extract(dat_sites_per_peotein$tst3,"(C\\(\\d{1,3}\\))") |>
    as.data.frame()  |>   rlang::set_names("NeuAc")
  NeuAc <- tidyr::extract(NeuAc,NeuAc,"C(","(\\d{1,3})") |>  as.data.frame()  |>  rlang::set_names("NeuAc")#NeuAc个数
  #Fuc
  Fuc <- stringr::str_extract(dat_sites_per_peotein$tst3,"(F\\(\\d{1,3}\\))")  |>
    as.data.frame()  |> rlang:: set_names("Fuc")
  Fuc <- tidyr::extract(Fuc,Fuc,"F(","(\\d{1,3})")  |>  as.data.frame()  |>  rlang::set_names("Fuc") #Fuc个数
  #合并再把NA替换成0
  dat_temp <- cbind(HexNAc,Hex,NeuAc,Fuc)  |>  as.data.frame()
  dat_temp[is.na(dat_temp)] <- 0
  dat_temp <- as.data.frame(lapply(dat_temp,as.numeric))#必须把提取出来的数字转化成Numeric

  dat_all <- dat_sites_per_peotein  |> cbind(dat_temp) |>  as.data.frame()


  if(!dir.exists(outdir)){dir.create(outdir)}

  data.table::fwrite(dat_all,paste0(outdir,"/glycan_type_count.tsv"),sep = "\t")


  #第二步----标记糖型
  dat_sites_gly = dat_all  |>
    dplyr::mutate(`Glycan Type`=dplyr::case_when(
      HexNAc<=2&Hex==0&NeuAc==0~'Paucimannose',
      HexNAc<=2&Hex<=3&NeuAc==0~'Paucimannose',
      HexNAc<3&Hex>3&Fuc==0&NeuAc==0~'High Mannose',
      HexNAc>=3&Hex>=3&Fuc==0&NeuAc==0~'Complex/Hybrid',
      NeuAc!=0~'Sialyated',
      TRUE~'Fucosylated'
    ))
  dat$Glyco_site = dat_sites_per_peotein$Glyco_site
  dat$`Glycan_Type` = dat_sites_gly$`Glycan Type`

  data.table::fwrite(dat,paste0(outdir,"/GlycanType_GlycanSites.tsv"),sep = "\t")


  dat =unique(dat)#去除完全重复的行

    dat_glycosite= dat[!duplicated(dat[c('protein', 'Glyco_site')]), ]  |>
      dplyr::select(protein,SiteperProt_pre=Glyco_site)


    dat_glycosite = dat_glycosite  |>
      dplyr::mutate(SiteperProt=ifelse(SiteperProt_pre>SiteperProteinsCutoff, paste0(">",SiteperProteinsCutoff),SiteperProt_pre ))

    dat_proteinsitecount=dat_glycosite  |>
      dplyr::group_by(SiteperProt)  |>
      dplyr::summarise(
        Prot_Num = dplyr::n()
      )

    dat_proteinsitecount$SiteperProt= factor(dat_proteinsitecount$SiteperProt, levels = c(as.character(1:SiteperProteinsCutoff), paste0(">",SiteperProteinsCutoff)))

    glycan.node =dat |> dplyr::select( glycancomposition, Glycan_Type)
    protein.node = dplyr::arrange(dat_glycosite, dplyr::desc(SiteperProt))  |>
      dplyr::select(protein)
    protein.node =protein.node |>dplyr:: mutate(Glycan_Type="protein")
    glycan.node = glycan.node[!duplicated(glycan.node[c('glycancomposition')]), ]
    glycan.order=c("Fucosylated", "High Mannose","Complex/Hybrid", "Paucimannose","Sialyated")
    glycan.node =glycan.node  |>
      dplyr::mutate(Glycan_Type = factor(Glycan_Type, levels = glycan.order)) |>
      dplyr::arrange(Glycan_Type)
    glycan.node = dplyr::arrange(glycan.node, Glycan_Type) # 糖型排序
    protein.node = protein.node[!duplicated(protein.node[c('protein')]), ]
    colnames(glycan.node)[1] <- "node.name"
    colnames(protein.node)[1] <- "node.name"
    node = vctrs::vec_c(glycan.node, protein.node)  |>
      dplyr::select(node.name)

    gplinks = dplyr::select(dat, glycancomposition, protein, Glycan_Type)
    gplinks = gplinks[!duplicated(gplinks[c('glycancomposition','protein')]), ]

    colnames(gplinks)= c("from","to", "Glycan_Type")

    g_glycopro_pairs =igraph::graph_from_data_frame(d=gplinks, directed = FALSE,vertices = node )
    g_glycopro_pairs_2 = igraph::graph_from_data_frame(d=gplinks, directed = FALSE,vertices = node )
    g_glycan =igraph::make_empty_graph()+igraph::vertices(glycan.node$node.name)
    g_protein = igraph::make_empty_graph()+igraph::vertices(protein.node$node.name)

    l_glycan0 <- igraph::layout_in_circle(g_glycan)
    l_glycan_polar <- useful::cart2pol(l_glycan0[, 1], l_glycan0[, 2])[, 1:2]
    theta <- l_glycan_polar$theta
    theta <- ifelse((theta>=0) & (theta<pi/2), theta*0.8 ,
                    ifelse((theta>=pi/2) & (theta<pi), pi-(pi-theta)*0.8,
                           ifelse((theta>=pi) & (theta<pi*3/2), pi+(theta-pi)*0.8, 2*pi-(2*pi-theta)*0.8)))
    l_glycan_polar$theta <- theta
    l_glycan <- as.matrix(useful::pol2cart(l_glycan_polar$r, l_glycan_polar$theta)[, 1:2])

    l_glycan[, 1] <- l_glycan[, 1]+sign(l_glycan[, 1])*0.1

    prot_len <- length(protein.node$node.name)
    l_protein <- matrix(0, prot_len, 2)
    l_protein[, 2] <- seq(-1, 1,length.out=prot_len);
    l_coordinate <- as.matrix(vctrs::vec_c(l_glycan, l_protein)[, 1:2])

    myColors =ggsci::pal_npg("nrc", alpha =0.6)(9)
    names(myColors) <- c("Fucosylated", "High Mannose","Complex/Hybrid", "Paucimannose","Sialyated")

    #bar_color1 <- c( "#8dd3c7", "#62C2B2","#389081","#2A6C61","#1E4E46")
    bar_color2 <- c("#CEC8E2","#a59aca","#8B7CBA","#6755A1","#44386A")

    fillColors <- bar_color2
    names(fillColors) <-  c(as.character(1:SiteperProteinsCutoff), paste0(">",SiteperProteinsCutoff))
    colScale <- ggplot2::scale_colour_manual(name = "Glycan Type",values = myColors)
    colFill <- ggplot2::scale_fill_manual(name = "Site/Protein", values=fillColors)
    net_test <- ggnetwork::ggnetwork(g_glycopro_pairs,layout = l_coordinate)  |>  dplyr::filter(!is.na(Glycan_Type))
    write.csv(net_test,file = paste0(outdir,"/",Title,"_glycan_type_count.csv"),row.names = F)



    fig <- ggplot2::ggplot(net_test) +
      ggnetwork::geom_edges(linewidth= edge_cex, ggplot2::aes(x, y, xend = xend, yend = yend, colour=Glycan_Type), alpha=edge_alpha) +
      ggnetwork::geom_nodes(ggplot2::aes(x, y, colour=Glycan_Type),size = node_size,alpha=node_alpha) +
      ggplot2::geom_bar(data=dat_proteinsitecount,stat = "identity",
               width = 0.03,ggplot2::aes(x=0.5, y = Prot_Num/sum(Prot_Num),
                                fill=as.factor(SiteperProt)))+
      colScale+
      colFill+
      ggplot2::coord_fixed(ratio = 1/1.1)+
      ggplot2::ggtitle(Title)+
      ggnetwork::theme_blank()+
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    print(fig)

  return(plt_network=fig)

}
