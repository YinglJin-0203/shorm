#' Function to visualize SHARP test results
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 annotate
#' @importFrom ggplot2 theme_minimal
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 coord_cartesian
#' @importFrom ggplot2 .pt
#' @importFrom scales trans_new
#'
#' @param pinc,pdec,pconc,pconv SHARP test p-values corresponding to increasing, decreasing, concave and convex
#' @param alpha significance threshold. alpha=NULL corresponding to no significance threshold
#' @param scale how axis should be scaled for better visual presentation. Choose between "none", "log" and "linear". See Details for the transformation fomula.
#' @param label labels of the point to visualize. label=NULL means no labeling.
#' @param size_point,size_label the size of points and label. size_lable is used when label is not NULL.
#' @param ... additional parameters passed to ggplot
#' @return The plotted ggplot object
#'
#' @details The shape scatter plot represents the four p-values (corresponding to four shape types) of SHARP test on a 2D space.
#' it can not only visualize the conclusions of the test, but also the confidence level of the test.
#' The plot surface is divided into two parts along two perpendicular directions, corresponding to two opposite shapes.
#' The coordinates of points are derived from p-values such that:
#' * Inconclusive shapes are placed in the center area within the significance thresholds (p > 0.05)
#' * Significant shapes are placed in the edge area out of the significance thresholds (p < 0.05)
#' * The closer points are to the significance thresholds, the more ambiguous shape conclusion is
#'
#' Axis can be scaled for better visual presentation. When scale = "log", a log-based transformation is implemented as follows:
#' \deqn{f(x) = 1-\frac{log(1+a(1-|x|))}{log(1+a)}} where \eqn{a=\frac{1-2\alpha}{\alpha^2}}.
#'
#' When scale = "linear", a piece-wise linear transformation is implemented:
#' \deqn{f(x) = 1-0.5\frac{(1-|x|)}{\alpha}I(|x|>1-\alpha)+\frac{0.5|x|}{1-\alpha}I(|x| \leq \alpha)}
#' In both cases, the significant and insignificant regions are also evenly divided in each test.

#'
#' @export
#' @examples
#' # Simulate dose-response data
#' y <- 2*sqrt(x)+rnorm(48)
#' curve$rep <- rep(1:3, each = 16)
#'
#' # Fixed-model based test
#' sharpt <- SHARPtest(curve, xName = "x", yName = "y", niter = 100)
#'
#' # Plot the graph
#' SharpScatter(sharpt[1], sharpt[2], sharpt[3], sharpt[4], niter = 100)

SharpScatter <- function(pinc, pdec, pconc, pconv, alpha = 0.05, scale = "log",
                         label = NULL, size_point=2, size_label = 5,...){
  # coordinates for points
  ycoord = 0*(pinc==pdec)+(1-pinc)*(pinc<pdec)+(pdec-1)*(pinc>pdec)
  xcoord = 0*(pconv==pconc)+(pconc<pconv)*(1-pconc)+(pconc>pconv)*(pconv-1)

  # scatter plot
  sharp_scatter <-ggplot(data.frame(xcoord, ycoord), aes(x=xcoord, y=ycoord,...))+
    geom_point(size = size_point)+
    geom_hline(yintercept = 0, linewidth = 1) +  # center y-axis
    geom_vline(xintercept = 0, linewidth = 1) +  # center x-axis
    labs(x=" ", y=" ")+
    coord_cartesian(clip = "off")+
    #theme_minimal()+
    theme(panel.grid = element_blank(),
          axis.ticks = element_blank())

  # visualization threshold
  if(!is.null(alpha)){
    sharp_scatter <- sharp_scatter+
      geom_hline(yintercept = c(1-alpha, -1+alpha), linetype = "dashed", linewidth = 0.5)+
      geom_vline(xintercept = c(1-alpha, -1+alpha), linetype = "dashed", linewidth = 0.5)
  }

  # scale significance threshold
  if(scale == "linear"){
    # piecewise linear transformation
    NewScale <- function(x){
      newx <- ifelse(abs(x)>(1-alpha), 1-((1-abs(x))/alpha)*0.5, abs(x)*0.5/(1-alpha))
      newx = sign(x)*newx
      return(newx)
    }
    InverseNewScale <- function(y){
      origy <- ifelse(abs(y) > 0.5, 1 - ((1 - abs(y)) * alpha / 0.5), abs(y) * (1-alpha) / 0.5)
      origy = sign(y) * origy
      return(origy)
    }
    NewTrans <- trans_new(
      name = "NewTrans",
      transform = NewScale,
      inverse = InverseNewScale
    )
    sharp_scatter <- sharp_scatter+
      scale_x_continuous(trans = NewTrans, limits = c(-1, 1),
                         breaks = c(alpha-1, 1-alpha), labels = c("p=0.05", "p=0.05"))+
      scale_y_continuous(trans = NewTrans, limits = c(-1, 1),
                         breaks = c(alpha-1, 1-alpha), labels = c("p=0.05", "p=0.05"))
  } else if(scale == "log") {
    # log-base transformation
    a <- (1 - 2 * alpha) / alpha^2
    NewScale <- function(x) {
      sign(x) * (1 - log(1 + a * (1 - abs(x))) / log(1 + a))
    }
    InverseNewScale <- function(y) {
      sign(y) * ((1 + a - (1 + a)^(1 - abs(y))) / a)
    }
    NewTrans <- trans_new(
      name = "NewTrans",
      transform = NewScale,
      inverse = InverseNewScale
    )
    sharp_scatter <- sharp_scatter+
      scale_x_continuous(trans = NewTrans, limits = c(-1, 1),
                         breaks = c(alpha-1, 1-alpha), labels = c("p=0.05", "p=0.05"))+
      scale_y_continuous(trans = NewTrans, limits = c(-1, 1),
                         breaks = c(alpha-1, 1-alpha), labels = c("p=0.05", "p=0.05"))

  } else {
    # no transformation
    sharp_scatter <- sharp_scatter+
      scale_x_continuous(breaks = c(alpha-1, 1-alpha), labels = c("p=0.05", "p=0.05"))+
      scale_y_continuous(breaks = c(alpha-1, 1-alpha), labels = c("p=0.05", "p=0.05"))

  }


  # point label
  if(!is.null(label)){
    sharp_scatter <- sharp_scatter+geom_text(aes(label=label), size = size_label)
  }

  # add axis label
  sharp_scatter <- sharp_scatter+
    annotate("text", x = c(-1, 1), y = c(0,0), label = c("Convex", "Concave"),
             vjust = 1.1, size = 15/.pt)+
    annotate("text", y = c(-1, 1), x = c(0,0), label = c("Decrease", "Increase"),
             hjust = 1.1, size = 15/.pt)

 return(sharp_scatter)

}





