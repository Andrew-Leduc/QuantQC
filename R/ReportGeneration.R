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
Gen_QQC_report_from_list <- function(data_path, linker_path, isolation1, isolation2 = NULL) {
  # Load necessary libraries

  # Set up parameters for the report
  file_paths <- list(data_path = data_path, linker_path = linker_path,isolation1 = isolation1, isolation2 = isolation2)

  # Generate the report using knitr
  rmarkdown::render(input = system.file("rmarkdown/QuantQC_DDA.Rmd", package = "QuantQC"),
                    output_format = "html_document",
                    output_file = "output_report.html",
                    params = file_paths)
}

