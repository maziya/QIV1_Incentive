library(ggvenn)
library(patchwork)

# df is the responder groups based on adjFC
strains <- c("HongKong", "Victoria", "Phuket", "Washington")

# ---- Function to build sets for a response ----
get_sets_response <- function(df, response) {
  df_logic <- df
  for(strain in strains) {
    df_logic[[strain]] <- df[[strain]] == response
  }
  
  sets <- lapply(strains, function(strain) {
    df_logic$SubjectID[df_logic[[strain]]]
  })
  names(sets) <- strains
  return(sets)
}

# ---- Build sets ----
sets_LR <- get_sets_response(df, "LR")
sets_MR <- get_sets_response(df, "MR")
sets_HR <- get_sets_response(df, "HR")

# ---- Plot Venn diagrams ----
plot_LR <- ggvenn(sets_LR,
                  fill_color = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF"),
                  stroke_size = 0.8,
                  set_name_size = 4,
                  show_percentage = FALSE) + ggtitle("Low Responders (LR)")

plot_MR <- ggvenn(sets_MR,
                  fill_color = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF"),
                  stroke_size = 0.8,
                  set_name_size = 4,
                  show_percentage = FALSE) + ggtitle("Medium Responders (MR)")

plot_HR <- ggvenn(sets_HR,
                  fill_color = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF"),
                  stroke_size = 0.8,
                  set_name_size = 4,
                  show_percentage = FALSE) + ggtitle("High Responders (HR)")


