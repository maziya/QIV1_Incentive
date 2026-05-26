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


library(ggvenn)
library(gridExtra)
library(dplyr)

strains  <- c("HongKong", "Victoria", "Phuket", "Washington")
colors_4 <- c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")
responses <- c("LR", "MR", "HR")
vax_groups <- c("Covaxin", "Covishield", "No_covid_vaccine")


get_sets_response <- function(df, response) {
  sets <- lapply(strains, function(strain) {
    col <- paste0("Responder_", strain)
    df$SubjectID[df[[col]] == response & !is.na(df[[col]])]
  })
  names(sets) <- strains
  return(sets)
}


plot_list <- list()

for (vax in vax_groups) {
  df_sub <- QIV1_meta_filt_cond_V1 %>% filter(covid_vaccinated == vax)
  
  for (resp in responses) {
    sets <- get_sets_response(df_sub, resp)
    p <- ggvenn(sets,
                fill_color      = colors_4,
                stroke_size     = 0.8,
                set_name_size   = 5,
                text_size       = 5,
                show_percentage = FALSE) +
      ggtitle(paste0(vax, " — ", resp)) +
      theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5))
    plot_list[[paste0(vax, "_", resp)]] <- p
  }
}


g <- arrangeGrob(
  plot_list[["Covaxin_LR"]],          plot_list[["Covaxin_MR"]],          plot_list[["Covaxin_HR"]],
  plot_list[["Covishield_LR"]],       plot_list[["Covishield_MR"]],       plot_list[["Covishield_HR"]],
  plot_list[["No_covid_vaccine_LR"]], plot_list[["No_covid_vaccine_MR"]], plot_list[["No_covid_vaccine_HR"]],
  ncol = 3,
  top  = textGrob(
    "Responder overlap across strains by COVID-19 vaccination status",
    gp = gpar(fontsize = 18, fontface = "bold")
  )
)
ggsave("venn_responders_by_vaccination.png", plot = g,
        width = 18, height = 16, units = "in", dpi = 300, bg = "white")
