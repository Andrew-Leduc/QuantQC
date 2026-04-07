### Report Generation

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
Gen_QQC_report_DDA <- function(data_path, linker_path, plex_exp,output_path = "output_report.html", isolation, CV_thresh = .4) {

    # Set up parameters for the report
    file_paths <- list(data_path = data_path, linker_path = linker_path, plex = plex_exp, isolation = isolation, CV_thresh = CV_thresh)


    # Generate the report using knitr
    rmarkdown::render(input = system.file("rmarkdown/QuantQC_DDA.Rmd", package = "QuantQC"),
                      output_format = "html_document",
                      output_file = output_path,
                      params = file_paths)

}

#' Generate a QC HTML report for a DIA (plexDIA / mTRAQ) experiment.
#'
#' Runs the full QuantQC DIA preprocessing and QC pipeline and renders the
#' results as a self-contained HTML report.
#'
#' @param data_path Path to the DIA-NN \code{report.tsv} file.
#' @param linker_path Path to the linker CSV mapping raw file names to wells/plates.
#' @param isolation Named list of paths to CellenONE \code{.xls} isolation files,
#'   one entry per cell type.
#' @param prot_vis_umap Optional UniProt accession to visualise on the UMAP.
#'   If \code{NULL} the UMAP section is skipped. Default \code{NULL}.
#' @param output_path Full path (including filename) for the output HTML report.
#'   Default \code{"output_report.html"}.
#' @param ChQ Channel Q-value threshold applied during peptide matrix construction.
#'   Default \code{0.1}.
#' @param plex_exp Integer plex number (2 or 3). Default \code{3}.
#' @param carrier_used Logical; whether a carrier channel was used. Default \code{FALSE}.
#' @return Invisibly returns the path to the rendered HTML file.
#' @export
Gen_QQC_report_DIA <- function(data_path, linker_path, isolation, prot_vis_umap = NULL, output_path = "output_report.html", ChQ = .1,
                               plex_exp = 3,carrier_used = F) {

  # Set up parameters for the report
  param_vals <- list(data_path = data_path, linker_path = linker_path, isolation = isolation,
                     prot_vis_umap = prot_vis_umap, ChQ = ChQ,plex_exp = plex_exp,
                     carrier_used = carrier_used)


  # Generate the report using knitr
  rmarkdown::render(input = system.file("rmarkdown/QuantQC_DIA.Rmd", package = "QuantQC"),
                    output_format = "html_document",
                    output_file = output_path,
                    params = param_vals)

}


#' Generate a QC HTML report for a Miceotopes (2-plex mTRAQ, no carrier) experiment.
#'
#' Runs the full Miceotopes preprocessing and QC pipeline — including
#' heavy/light peptide trimming and LC batch correction — and renders the
#' results as a self-contained HTML report. The report contains the same
#' sections as \code{Gen_QQC_report_DIA} but uses
#' \code{Miceotope_cellXpeptide} for matrix construction and omits the
#' carrier-channel statistics section.
#'
#' @param data_path Path to the DIA-NN \code{report.tsv} file.
#' @param linker_path Path to the linker CSV mapping raw file names to wells/plates.
#' @param isolation Named list of paths to CellenONE \code{.xls} isolation files,
#'   one entry per cell type.
#' @param prot_vis_umap Optional UniProt accession to visualise on the UMAP.
#'   If \code{NULL} the UMAP section is skipped. Default \code{NULL}.
#' @param output_path Full path (including filename) for the output HTML report.
#'   Default \code{"output_report.html"}.
#' @param ChQ Channel Q-value threshold. Default \code{1} (no filtering).
#' @param plex_exp Integer plex number. Default \code{2}.
#' @param t Minimum number of observations required to retain a heavy/light
#'   peptide pair in \code{Miceotope_cellXpeptide}. Default \code{5}.
#' @return Invisibly returns the path to the rendered HTML file.
#' @export
Gen_QQC_report_Miceotopes <- function(data_path, linker_path, isolation,
                                      prot_vis_umap = NULL,
                                      output_path = "output_report.html",
                                      ChQ = 1, plex_exp = 2, t = 5) {

  param_vals <- list(
    data_path     = data_path,
    linker_path   = linker_path,
    isolation     = isolation,
    prot_vis_umap = prot_vis_umap,
    ChQ           = ChQ,
    plex_exp      = plex_exp,
    t             = t
  )

  rmarkdown::render(
    input         = system.file("rmarkdown/QuantQC_Miceotopes.Rmd", package = "QuantQC"),
    output_format = "html_document",
    output_file   = output_path,
    params        = param_vals
  )

}


