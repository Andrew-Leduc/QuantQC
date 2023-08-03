
#' test1.
#'
#' This function takes two numeric inputs and returns their sum.
#'
#' @param x A numeric value.
#' @param y A numeric value.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' add_numbers(2, 3)
#' @export

PlotSlideLayout_celltype <- function(nPOP_obj){

  ggplot(nPOP_obj@cellenONE.meta) +
    geom_point(aes(x = dropXPos,y = dropYPos,color = sample)) +
    geom_text(aes(x = pickupXPos_numb,y = pickupYPos_numb,label = injectWell,size = 5),hjust= .5, vjust=-.6) +
    facet_wrap(~field,ncol = 4)+
    scale_y_reverse()

}


#' test1.
#'
#' This function takes two numeric inputs and returns their sum.
#'
#' @param x A numeric value.
#' @param y A numeric value.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' add_numbers(2, 3)
#' @export
PlotSlideLayout_label <- function(nPOP_obj){

  # print the mTRAQ labels overlayed on the positions of the slide
  ggplot(nPOP_obj@cellenONE.meta, aes(x = dropXPos,y = dropYPos,color = label)) +
     geom_point() +scale_y_reverse()+ facet_wrap(~field,ncol = 4)

}


