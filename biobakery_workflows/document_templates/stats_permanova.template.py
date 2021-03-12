#' # One against all: Permanova test

#' <% if vars["beta_diversity_plots"]["univariate"]: print("## Univariate") %>

#' <% utilities.show_all_permanova(vars["permanova_plots"]) %>

#' <% utilities.show_all_variate_plots("univariate",vars["beta_diversity_plots"]) %>

#' <% if vars["beta_diversity_plots"]["multivariate"] and pdf_format: print("\clearpage") %>

#' <% if vars["beta_diversity_plots"]["multivariate"]: print("## Multivariate") %>

#' <% if vars["beta_diversity_plots"]["multivariate"]: print("For the multivariate model the following covariate equation was provided: 'bray ~ "+vars["covariate_equation"]+"' .") %>

#' <% utilities.show_all_variate_plots("multivariate",vars["beta_diversity_plots"]) %>

#' <% if pdf_format: print("\clearpage") %>
