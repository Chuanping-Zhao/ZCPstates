#' Easy PCA with Grouping and Scaling
#'
#' Performs a principal component analysis (PCA) on a given data frame with options for data scaling
#' and multiple grouping variables. This function returns the PCA model, the transformed PCA data with
#' group annotations, and ggplot2 visualizations for easy interpretation of PCA results.
#'
#' @param df A data frame containing the input data. It should include a column for sample identifiers,
#'   one or more grouping variable columns, and numeric feature columns for PCA.
#' @param sample.id A string specifying the column name in `df` that contains the sample identifier for each observation.
#' @param mutigoup.id A character vector of column names in `df` corresponding to grouping variables.
#'   These variables will be used to label or color points in the PCA plots (e.g., different colors for each group).
#'   If more than one grouping variable is provided, multiple visualizations may be generated.
#' @param scaling A string specifying the scaling or transformation method to apply to the numeric data before PCA.
#'   Options are:
#'   \itemize{
#'     \item `"pareto"`: Pareto scaling (mean-center each variable and divide by the square root of its standard deviation).
#'     \item `"uv"`: Unit variance scaling (mean-center each variable and divide by its standard deviation, i.e., standardize to variance = 1).
#'     \item `"logtrans"`: Log transformation of each numeric value (e.g., base-10 logarithm), useful for reducing skewness.
#'     \item `"noop"`: No scaling or transformation (use the raw data as is).
#'   }
#' @param procomp.center Logical; if `TRUE`, center the data (subtract the mean of each variable) before performing PCA.
#'   This is passed to the `center` argument of [stats::prcomp()] and is generally recommended unless the data is already centered.
#' @param procomp.scale Logical; if `TRUE`, scale the data to unit variance (divide by the standard deviation of each variable) before PCA.
#'   This is passed to the `scale.` argument of [stats::prcomp()]. Using `procomp.center = TRUE` and `procomp.scale = TRUE` together is equivalent to standardizing the data.
#' @param pointsize Numeric value indicating the size of points in the PCA scatter plot(s).
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{\code{pca_model}}{The PCA model object returned by \code{\link[stats]{prcomp}}. This includes the rotation (loadings for each principal component), standard deviations of principal components, center and scale used, and the PCA scores (accessible via `pca_model$x`).}
#'   \item{\code{pca_data}}{A data frame (or tibble) of PCA-transformed data. Each row corresponds to a sample from the original dataset and contains its principal component scores (e.g., PC1, PC2, ...), along with the original sample ID and grouping variable values for reference.}
#'   \item{\code{plots}}{A named list of \code{\link[ggplot2]{ggplot}} objects for visualizing the PCA results. For example, this may include a scatter plot of PC1 vs PC2 with points colored by one of the grouping variables. If multiple grouping variables are provided in `mutigoup.id`, the list may contain a separate PCA plot for each grouping variable (with list names corresponding to the grouping variable names).}
#' }
#'
#' @importFrom stats prcomp
#' @importFrom ggplot2 ggplot aes geom_point labs theme
#' @importFrom dplyr select mutate n
#' @importFrom purrr map
#' @importFrom tidyr unnest
#'
#' @examples
#' #demo
#'
#' #Load example dataset and prepare data
#' data(iris)
#' dt=iris |>  dplyr::mutate(group2 = sample(c("A", "B"), dplyr::n(), replace = TRUE)) |> dplyr::rename(group1=Species) |> dplyr::mutate(sample=paste0(group1,"_",dplyr::row_number()))
#' dt$group1 <- as.character(dt$group1)
#' dt$group2 <- as.character(dt$group2)
#'
#' colnames(dt)[1:4] <- paste0("Feature", 1:4)
#'
#' #pca
#' result <- easypca(
#' df = dt,
#' sample.id = "sample",
#' mutigoup.id = c("group1","group2"),
#' scaling = "uv",
#' procomp.center = TRUE,
#' procomp.scale = TRUE,
#' pointsize = 1.5
#' )
#'
#' #Inspect the PCA model and data
#' print(result$pca_model)    #PCA model details (prcomp output)
#' head(result$pca_data)      #First few rows of PCA scores with grouping info
#'
#' #Visualize PCA plots
#' print(result$plots$group1)        #Print the muti-group PCA analysis
#' print(result$plots$group2)
#'
#' @export
#'
easypca=function(df,
             sample.id="File ID",
             #group.id="gender",
             mutigoup.id=c("gender","age_range"),#如果mutigroup=T 就需要指定
             scaling="pareto", # pareto uv logtrans noop无任何操作
             procomp.center=TRUE,
             procomp.scale=FALSE,
             pointsize=1.5
){

  #定义scaling
  if (is.null(scaling)) {
    cal =function(x){(x-mean(x, na.rm = T))/sd(x, na.rm = T)}#uv
  }else if(scaling=="pareto"){
    cal=function(x){(x-mean(x,na.rm=T))/sqrt(sd(x, na.rm = T))}#pareto
  }else if(scaling=="uv"){
    cal =function(x){(x-mean(x, na.rm = T))/sd(x, na.rm = T)}#uv
  }else if(scaling=="logtrans"){
    cal = function(x) { ifelse(is.na(x), NA,log2(pmax(x, 0) + 1))}
  }else if(scaling=="noop"){
    cal=function(x){ifelse(is.na(x), NA, x)}
  }



  samplecol <- sym(sample.id)


    #如果矩阵中存在多种分组方案 需要指出这些分组方案的列名 否则抛出错误
   # if(length(unique(mutigoup.id))<2){stop("The column names provided for the `mutigoup.id` are fewer than 2. Please set `mutigoup.id=FALSE`.--zcp")}
    #然后判断提供的分组列名信息是否存在
    if(!all(mutigoup.id %in% colnames(df))){stop("The column names provided for the `mutigoup.id` are not in `colnames(df)`. Please check it.--zcp")}
    #去掉样本id和分组id计算pca
    plot.dat <- df  |>
      dplyr::select(!dplyr::all_of(c(as.character(samplecol),mutigoup.id)))  |>
      dplyr:: mutate( dplyr::across( dplyr::everything(), as.numeric))

    plot.dat=as.data.frame(apply(plot.dat, 2, cal))

    plot.dat[is.na(plot.dat)]=0

    #提取group信息 提取sample信息
    groupInfo=df |> dplyr::select(dplyr::all_of(c(as.character(samplecol),mutigoup.id)))

    group =df  |> dplyr::select(dplyr::all_of(mutigoup.id))
    plot.dat = plot.dat[, apply(plot.dat, 2, sd, na.rm = TRUE) > 0]#保留sd不为0的数
    pca=prcomp(plot.dat, center =procomp.center, scale. =procomp.scale)
    pca.data = data.frame(pca$x) |> cbind(groupInfo) |> dplyr::select(dplyr::any_of(colnames(group)),dplyr::everything())#return
    pca.variance= pca$sdev^2 / sum(pca$sdev^2)


    generate_pca_plots = function(pca_data, group_vars, variance, sample_col, pointsize) {
      color_palettes= c(
        "#00ABF0", "#39B243", "#25848E", "#34618DFF", "#CB3E71FF",
        "#440154FF", "#E54C5E", "#EE822F", "#B5379A", "#CC2020",
        "#0505E7", "#A50021", "palevioletred3", "#ceca7c", "#c59fc9",
        "#84b59f", "cornflowerblue", "salmon3","#5698c4", "#d88c9a",
        "#3b374c", "#44598e", "#64a0c0", "#7ec4b7", "#deebcd","#073f82",
        "#1b71b4", "#58a4cf", "#a2cbe3", "#f2f9fe","#492952", "#82677e",
        "#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7"
        , "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"
      )

      purrr::map(group_vars, function(group) {#group="gender"
        #计算椭圆数据
        #  ellipse_data <- pca_data |>
        #   dplyr::group_by(!!rlang::sym(group)) |>
        #   dplyr::summarize(ggplot2::StatEllipse$compute_group(
        #   data = data.frame(x = PC1, y = PC2),
        #   scales = list(),
        #   level = 0.95
        #  ))

        #内置函数计算置信区间
        .compute_ellipse_points = function(df, x_col, y_col, group_col, level = 0.95) {
          ellipse_data <- df  |>
            dplyr::group_by(!!rlang::sym(group_col)) |>
            dplyr::summarise(ellipse = list(ggplot2::StatEllipse$compute_group(
              data = tibble(x = !!rlang::sym(x_col), y = !!rlang::sym(y_col)),
              scales = NULL,
              level = level
            )), .groups = "drop")  |>
            tidyr::unnest(cols = c(ellipse))

          return(ellipse_data)
        }
        #计算
        ellipse_data = .compute_ellipse_points(
          df = pca_data,
          x_col = "PC1",
          y_col = "PC2",
          group_col = group
        )

        #计算 x 和 y 轴范围
        x_min=  min(ellipse_data$x, na.rm = TRUE)
        x_max =  max(ellipse_data$x, na.rm = TRUE)
        y_min=  min(ellipse_data$y, na.rm = TRUE)
        y_max =  max(ellipse_data$y, na.rm = TRUE)

        #计算 buffer，避免边界贴近图像
        buffer_x= (x_max - x_min) * 0.01
        buffer_y=  (y_max - y_min) * 0.01

        ggplot2::ggplot(data = pca_data, mapping = ggplot2::aes(
          x = PC1, y = PC2,
          fill =  !!rlang::sym(group),
          color =  !!rlang::sym(group)
        )) +
          ggplot2::stat_ellipse(linewidth = 0.5, level = 0.95, geom = "polygon", alpha = 0.1,color=NA) +
          ggplot2::geom_point(size = pointsize) +
          ggplot2::theme_bw() +
          ggplot2::scale_color_manual(values = color_palettes) +
          ggplot2::scale_fill_manual(values = color_palettes) +
          ggplot2::labs(
            x = paste0("PC1: ", signif(variance[1] * 100, 3), "%"),
            y = paste0("PC2: ", signif(variance[2] * 100, 3), "%"),
            title = "PCA Analysis"
          ) +
          ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5),
            plot.background = ggplot2::element_rect(fill = "white", color = NA),
            panel.background = ggplot2::element_rect(fill = "white", color = NA),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid = ggplot2::element_blank(),
            panel.ontop = FALSE
          ) +
          ggplot2::geom_text(ggplot2::aes(label = !!rlang::sym(sample_col)),# color = !!rlang::sym(mutigoup.id[1])),
                             check_overlap = T,
                             show.legend = FALSE, vjust = -0.5, size = 2) +
          ggplot2::xlim(x_min - buffer_x, x_max + buffer_x) +
          ggplot2::ylim(y_min - buffer_y, y_max + buffer_y)


      }) |> stats::setNames(group_vars)
    }


    plot_list <- generate_pca_plots(
      pca_data = pca.data,
      group_vars = mutigoup.id,
      variance = pca$sdev^2 / sum(pca$sdev^2),
      sample_col = sample.id,
      pointsize = pointsize
    )

    return(list(
      pca_model = pca,
      pca_data = pca.data,
      plots = plot_list
    ))



}

