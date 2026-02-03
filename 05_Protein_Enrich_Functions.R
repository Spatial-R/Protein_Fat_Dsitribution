# Helper function for reordering within facets
reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun)
}

scale_y_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  ggplot2::scale_y_discrete(labels = function(x) gsub(reg, "", x), ...)
}

#' Set SCI-quality theme for all plots
set_sci_theme <- function(base_size = 10) {
  theme_set(
    theme_classic(base_size = base_size) +
      theme(
        # Text elements
        text = element_text(family = "serif", color = "black"),
        axis.text = element_text(color = "black", size = base_size * 0.9),
        axis.title = element_text(color = "black", size = base_size, face = "bold"),
        plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = base_size, hjust = 0.5),
        legend.title = element_text(size = base_size * 0.9, face = "bold"),
        legend.text = element_text(size = base_size * 0.8),
        strip.text = element_text(size = base_size * 0.9, face = "bold"),
        
        # Axis lines
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        
        # Panel
        panel.background = element_rect(fill = "white", color = NA),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        
        # Margins
        #plot.margin = margin(10, 10, 10, 10, "pt"),
        
        # Legends
        legend.position = "right",
        legend.background = element_rect(fill = "white", color = "grey80"),
        legend.key = element_rect(fill = "white", color = NA),
        legend.key.size = unit(0.4, "cm"),
        legend.spacing = unit(0.2, "cm")
      )
  )
}

set_sci_theme()


convert_protein_ids <- function(protein_names, species = org.Hs.eg.db) {
  message("Converting protein names to Entrez IDs...")
  
  id_mapping <- bitr(protein_names,
                     fromType = "SYMBOL",
                     toType = c("ENTREZID", "UNIPROT", "ENSEMBL"),
                     OrgDb = species,
                     drop = FALSE)
  
  original_count <- length(protein_names)
  mapped_count <- length(unique(id_mapping$SYMBOL))
  
  message(paste("Successfully mapped:", mapped_count, "/", original_count, 
                "proteins (", round(mapped_count/original_count*100, 1), "%)"))
  
  # Report unmapped proteins
  unmapped <- setdiff(protein_names, id_mapping$SYMBOL)
  if (length(unmapped) > 0) {
    message(paste("Unmapped proteins:", length(unmapped)))
    if (length(unmapped) <= 10) {
      message(paste("  ", paste(unmapped, collapse = ", ")))
    } else {
      message(paste("  ", paste(head(unmapped, 10), collapse = ", "), "..."))
    }
  }
  
  return(id_mapping)
}


run_go_enrichment <- function(entrez_ids, species = org.Hs.eg.db, 
                              p_cutoff = 0.05, q_cutoff = 0.2) {
  message("Running GO enrichment analysis...")
  
  go_results <- list()
  
  # Biological Process
  go_results$BP <- enrichGO(gene = entrez_ids,
                            OrgDb = species,
                            ont = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff = p_cutoff,
                            qvalueCutoff = q_cutoff,
                            readable = TRUE)
  
  # Cellular Component
  go_results$CC <- enrichGO(gene = entrez_ids,
                            OrgDb = species,
                            ont = "CC",
                            pAdjustMethod = "BH",
                            pvalueCutoff = p_cutoff,
                            qvalueCutoff = q_cutoff,
                            readable = TRUE)
  
  # Molecular Function
  go_results$MF <- enrichGO(gene = entrez_ids,
                            OrgDb = species,
                            ont = "MF",
                            pAdjustMethod = "BH",
                            pvalueCutoff = p_cutoff,
                            qvalueCutoff = q_cutoff,
                            readable = TRUE)
  
  # Report significant results
  message(paste("Significant GO terms:"))
  message(paste("  BP:", sum(go_results$BP$p.adjust < p_cutoff)))
  message(paste("  CC:", sum(go_results$CC$p.adjust < p_cutoff)))
  message(paste("  MF:", sum(go_results$MF$p.adjust < p_cutoff)))
  
  return(go_results)
}


run_kegg_enrichment <- function(entrez_ids, organism = "hsa", 
                                p_cutoff = 0.05, q_cutoff = 0.2) {
  message("Running KEGG pathway enrichment analysis...")
  
  kegg_result <- enrichKEGG(gene = entrez_ids,
                            organism = organism,
                            pvalueCutoff = p_cutoff,
                            qvalueCutoff = q_cutoff)
  
  message(paste("Significant KEGG pathways:", 
                ifelse(!is.null(kegg_result), sum(kegg_result$p.adjust < p_cutoff), 0)))
  
  return(kegg_result)
}


prepare_go_data <- function(go_result, ontology = "BP", top_n = 15) {
  
  if (is.null(go_result) || nrow(go_result) == 0) {
    return(NULL)
  }
  
  # Convert to data frame and filter significant terms
  go_df <- as.data.frame(go_result) %>%
    filter(p.adjust < 0.05) %>%
    arrange(p.adjust)
  
  if (nrow(go_df) == 0) {
    return(NULL)
  }
  
  # Take top N terms
  go_df <- go_df %>%
    head(min(top_n, nrow(go_df)))
  
  # Calculate -log10(p.adjust)
  go_df$log10_padj <- -log10(go_df$p.adjust)
  
  # Calculate gene ratio as numeric
  go_df <- go_df %>%
    rowwise() %>%
    mutate(
      GeneRatio_num = as.numeric(strsplit(GeneRatio, "/")[[1]][1]) / 
        as.numeric(strsplit(GeneRatio, "/")[[1]][2]),
      BgRatio_num = as.numeric(strsplit(BgRatio, "/")[[1]][1]) / 
        as.numeric(strsplit(BgRatio, "/")[[1]][2]),
      Enrichment = GeneRatio_num / BgRatio_num
    ) %>%
    ungroup()
  
  # Add ontology label
  go_df$Ontology <- case_when(
    ontology == "BP" ~ "Biological Process",
    ontology == "CC" ~ "Cellular Component", 
    ontology == "MF" ~ "Molecular Function"
  )
  
  # Shorten long descriptions
  go_df$ShortDescription <- sapply(go_df$Description, function(x) {
    words <- strsplit(x, " ")[[1]]
    if (length(words) > 6) {
      paste(paste(words[1:6], collapse = " "), "...", sep = "")
    } else {
      x
    }
  })
  
  # Reorder by significance
  go_df <- go_df %>%
    mutate(ShortDescription = fct_reorder(ShortDescription, p.adjust))
  
  return(go_df)
}


create_go_bubble_plot <- function(go_data, ontology_color = "#1f78b4", 
                                  show_labels = TRUE) {
  
  if (is.null(go_data) || nrow(go_data) == 0) {
    return(NULL)
  }
  
  # Create the plot
  p <- ggplot(go_data, 
              aes(x = GeneRatio_num, 
                  y = reorder(ShortDescription, GeneRatio_num),
                  size = Count,
                  color = log10_padj)) +
    geom_point(alpha = 0.9) +
    scale_color_gradientn(
      colors = c("lightblue", ontology_color, "darkblue"),
      name = expression(-log[10](p[adj])),
      limits = c(0, max(go_data$log10_padj) * 1.1),
      guide = guide_colorbar(order = 1)
    ) +
    scale_size_continuous(
      range = c(1, 5),
      name = "Gene Count",
      breaks = pretty_breaks(n = 4),
      guide = guide_legend(order = 2)
    ) +
    labs(
      x = "Gene Ratio",
      y = "",
      title = paste("GO", go_data$Ontology[1], "Enrichment")
    ) +
    theme(
      axis.text.y = element_text(size = 10, lineheight = 0.9),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      legend.position = "right",
      legend.box = "vertical",
      legend.spacing = unit(0.2, "cm")
    )
  
  # Add value labels if requested
  if (show_labels) {
    p <- p + 
      geom_text(aes(label = round(log10_padj, 1)),
                color = "black", size = 3, vjust = -1.5)
  }
  return(p)
}


create_kegg_dot_plot <- function(kegg_result, top_n = 15) {
  
  if (is.null(kegg_result) || nrow(kegg_result) == 0) {
    return(NULL)
  }
  
  # Prepare data
  kegg_df <- as.data.frame(kegg_result) %>%
    filter(p.adjust < 0.05) %>%
    arrange(p.adjust) %>%
    head(top_n)
  
  if (nrow(kegg_df) == 0) {
    return(NULL)
  }
  
  # Calculate gene ratio as numeric
  kegg_df <- kegg_df %>%
    rowwise() %>%
    mutate(
      GeneRatio_num = as.numeric(strsplit(GeneRatio, "/")[[1]][1]) / 
        as.numeric(strsplit(GeneRatio, "/")[[1]][2]),
      log10_padj = -log10(p.adjust)
    ) %>%
    ungroup()
  
  # Shorten pathway names
  kegg_df$ShortDescription <- sapply(kegg_df$Description, function(x) {
    if (nchar(x) > 50) {
      paste0(substr(x, 1, 47), "...")
    } else {
      x
    }
  })
  
  # Reorder by gene ratio
  kegg_df <- kegg_df %>%
    mutate(ShortDescription = fct_reorder(ShortDescription, GeneRatio_num))
  
  # Create plot
  p <- ggplot(kegg_df, 
              aes(x = GeneRatio_num, 
                  y = reorder(ShortDescription, GeneRatio_num),
                  size = Count,
                  color = log10_padj)) +
    geom_point(alpha = 0.9) +
    scale_color_viridis_c(
      option = "plasma",
      name = expression(-log[10](p[adj])),
      guide = guide_colorbar(order = 1)
    ) +
    scale_size_continuous(
      range = c(1, 5),
      name = "Gene Count",
      breaks = pretty_breaks(n = 4),
      guide = guide_legend(order = 2)
    ) +
    labs(
      x = "Gene Ratio",
      y = ""#,
      #title = "KEGG Pathway Enrichment",
      #subtitle = paste("Top", nrow(kegg_df), "significant pathways")
    ) +
    theme(
      axis.text.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      legend.position = "right"
    )
  
  return(p)
}


create_multi_panel_go_figure <- function(go_results, 
                                         output_file = NULL,
                                         width = 14, 
                                         height = 12, 
                                         dpi = 300) {
  
  # Prepare data for each ontology
  bp_data <- prepare_go_data(go_results$BP, "BP", 10)
  cc_data <- prepare_go_data(go_results$CC, "CC", 10)
  mf_data <- prepare_go_data(go_results$MF, "MF", 10)
  
  all_data <- bind_rows(
    if (!is.null(bp_data)) bp_data %>% mutate(panel = "BP"),
    if (!is.null(cc_data)) cc_data %>% mutate(panel = "CC"),
    if (!is.null(mf_data)) mf_data %>% mutate(panel = "MF")
  )
  
  if (nrow(all_data) == 0) {
    warning("No data to plot")
    return(NULL)
  }
  scale_y_reordered <- function(..., sep = "___") {
    reg <- paste0(sep, ".+$")
    ggplot2::scale_y_discrete(labels = function(x) gsub(reg, "", x), ...)
  }
  
  combined_plot <- ggplot(all_data, 
         aes(x = GeneRatio_num, 
             y = reorder_within(ShortDescription, GeneRatio_num, panel),
             size = Count,
             color = log10_padj)) +
    geom_point(alpha = 0.85, shape = 16) +
    facet_grid(panel ~ ., scales = "free_y", space = "free_y") +
    scale_color_viridis_c(
      option = "plasma",
      name = expression(-log[10](p[adj])),
      limits = c(0, max(all_data$log10_padj, na.rm = TRUE) * 1.1),
      guide = guide_colorbar(
        barwidth = unit(0.5, "cm"),
        barheight = unit(2, "cm")
      )
    ) +
    scale_size_continuous(
      range = c(1, 5),
      name = "Gene Count",
      breaks = pretty_breaks(n = 4),
      limits = range(all_data$Count, na.rm = TRUE)
    ) +
    scale_y_reordered() +
    labs(
      x = "Gene Ratio",
      y = "GO Terms"#,
      #title = "GO Enrichment Analysis with Standardized Scales",
      #subtitle = "Shared color and size scales for direct comparison across ontologies"
    ) +
    theme_pubr(base_family = "serif") +
    theme(
      strip.text = element_text(face = "bold", size = 12),
      strip.background = element_rect(fill = "grey90", color = "black"),
      panel.spacing = unit(1, "lines"),
      legend.position = "right",
      axis.text.y = element_text(size = 9, lineheight = 0.9))

  # Save if output file specified
  if (!is.null(output_file)) {
    ggsave(output_file, combined_plot,
           width = width, height = height,
           dpi = dpi, bg = "white")
    message("Multi-panel GO figure saved to: ", output_file)
  }
  return(combined_plot)
}


create_comprehensive_enrichment_figure <- function(go_results, kegg_result,
                                                   output_dir = "SCI_Figures",
                                                   study_title = "Proteomics Study") {
  
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 1. Multi-panel GO figure
  go_figure <- create_multi_panel_go_figure(
    go_results,output_file = file.path(output_dir, "Figure_GO_Enrichment.tiff")
  )
  
  # 2. KEGG pathway figure
  if (!is.null(kegg_result) && nrow(kegg_result) > 0) {
    kegg_plot <- create_kegg_dot_plot(kegg_result, top_n = 15)
    
    if (!is.null(kegg_plot)) {
      ggsave(file.path(output_dir, "Figure_KEGG_Pathways.tiff"),
             kegg_plot, width = 12, height = 10, dpi = 300, bg = "white")
    }
  }
  
  # 3. GO network visualization
  if (!is.null(go_results$BP) && nrow(go_results$BP) >= 5) {
    network_plot <- tryCatch({
      cnetplot(go_results$BP,
               showCategory = 10,
               colorEdge = TRUE,
               circular = FALSE,
               size_category = 1.5,
               node_label_size = "Count",
               #node_label = "category",  # Display only category labels to reduce overlap
               cex_label_category = 0.9, # Category label size
               cex_label_gene = 0.7,     # Gene label size
               color_category = "firebrick",  # Category node color
               color_item = "steelblue") +    # Gene node color
        # scale_color_gradientn(
        #   name = "Adjusted p-value",
        #   colors = c("#2E86AB", "#A23B72", "#F18F01"),  # Custom color gradient
        #   breaks = c(0.01, 0.05, 0.1),
        #   labels = c("1e-02", "5e-02", "1e-01")) +
        
        theme_minimal(base_size = 12) +
        theme(
          plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
          legend.position = "right",
          legend.title = element_text(face = "bold"),
          legend.text = element_text(size = 9),
          panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          plot.background = element_rect(fill = "white", color = NA))
    }, error = function(e) {
      message("Network plot failed: ", e$message)
      return(NULL)
    })
    
    if (!is.null(network_plot)) {
      ggsave(file.path(output_dir, "Figure_GO_Network.tiff"),
             network_plot, width = 16, height = 12, dpi = 300, bg = "white")
    }
  }
  
  # 4. Create summary dashboard
  summ_fig <- create_enrichment_dashboard(go_results, kegg_result, output_dir)
  
  return(list(
    go_figure = go_figure,
    summ_fig = summ_fig,
    kegg_plot = if(exists("kegg_plot")) kegg_plot else NULL,
    network_plot = if(exists("network_plot")) network_plot else NULL
  ))
}


create_enrichment_dashboard <- function(go_results, kegg_result, output_dir) {
  
  # Prepare data for summary plots
  summary_data <- data.frame()
  
  # Add GO summary
  for (ontology in c("BP", "CC", "MF")) {
    if (!is.null(go_results[[ontology]]) && nrow(go_results[[ontology]]) > 0) {
      go_df <- as.data.frame(go_results[[ontology]])
      sig_count <- sum(go_df$p.adjust < 0.05)
      summary_data <- rbind(summary_data,
                            data.frame(Category = paste("GO", ontology),
                                       Count = sig_count,
                                       Type = "GO"))
    }
  }
  
  # Add KEGG summary
  if (!is.null(kegg_result) && nrow(kegg_result) > 0) {
    kegg_df <- as.data.frame(kegg_result)
    sig_kegg <- sum(kegg_df$p.adjust < 0.05)
    summary_data <- rbind(summary_data,
                          data.frame(Category = "KEGG Pathways",
                                     Count = sig_kegg,
                                     Type = "Pathway"))
  }
  
  # Create summary bar plot
  if (nrow(summary_data) > 0) {
    summary_plot <- ggplot(summary_data,
                           aes(x = reorder(Category, Count),
                               y = Count,
                               fill = Type)) +
      geom_bar(stat = "identity", width = 0.7) +
      geom_text(aes(label = Count), 
                hjust = -0.3, size = 4, fontface = "bold") +
      scale_fill_manual(values = c("GO" = "#1f78b4", "Pathway" = "#33a02c")) +
      labs(x = "", y = "Number of Significant Terms",
           title = "Enrichment Analysis Summary",
           subtitle = "Count of significantly enriched terms and pathways") +
      coord_flip() +
      theme(
        legend.position = "none",
        axis.text = element_text(size = 11, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold")
      )
    
    ggsave(file.path(output_dir, "Figure_Enrichment_Summary.tiff"),
           summary_plot, width = 10, height = 6, dpi = 600, bg = "white")
  }
  return(summary_plot)
}


export_enrichment_results <- function(id_mapping, 
                                      go_results, 
                                      kegg_result = NULL, 
                                      output_dir = "results") {
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  message(paste("Exporting results to:", output_dir))
  
  # 1. Save ID mapping
  write.csv(id_mapping, 
            file.path(output_dir, "id_mapping.csv"), 
            row.names = FALSE,
            quote = FALSE)
  
  # 2. Save detailed GO results
  for (ontology in c("BP", "CC", "MF")) {
    if (!is.null(go_results[[ontology]]) && nrow(go_results[[ontology]]) > 0) {
      go_df <- as.data.frame(go_results[[ontology]]) %>%
        mutate(across(where(is.numeric), round, 4))
      
      readr::write_csv(go_df,
                file.path(output_dir, paste0("GO_", ontology, "_detailed.csv")))
    }
  }
  
  # 3. Save KEGG results
  if (!is.null(kegg_result) && nrow(kegg_result) > 0) {
    kegg_df <- as.data.frame(kegg_result) %>%
      mutate(across(where(is.numeric), round, 4))
    
    readr::write_csv(kegg_df,
              file.path(output_dir, "KEGG_pathways_detailed.csv"))
  }
  
  # 4. Create summary table
  summary_table <- data.frame()
  
  for (ontology in c("BP", "CC", "MF")) {
    if (!is.null(go_results[[ontology]]) && nrow(go_results[[ontology]]) > 0) {
      top_terms <- as.data.frame(go_results[[ontology]]) %>%
        filter(p.adjust < 0.05) %>%
        arrange(p.adjust) %>%
        head(5) %>%
        dplyr::select(Description, GeneRatio, p.adjust, Count) %>%
        mutate(Ontology = paste("GO", ontology),
               p.adjust = format(p.adjust, scientific = TRUE, digits = 3))
      
      summary_table <- bind_rows(summary_table, top_terms)
    }
  }
  
  if (!is.null(kegg_result) && nrow(kegg_result) > 0) {
    kegg_top <- as.data.frame(kegg_result) %>%
      filter(p.adjust < 0.05) %>%
      arrange(p.adjust) %>%
      head(5) %>%
      dplyr::select(Description, GeneRatio, p.adjust, Count) %>%
      mutate(Ontology = "KEGG Pathways",
             p.adjust = format(p.adjust, scientific = TRUE, digits = 3))
    
    summary_table <- bind_rows(summary_table, kegg_top)
  }
  
  if (nrow(summary_table) > 0) {
    write.csv(summary_table,
              file.path(output_dir, "enrichment_summary_top_terms.csv"),
              row.names = FALSE,
              quote = FALSE)
  }
}


run_complete_enrichment_analysis <- function(protein_names, 
                                             organism = "hsa", 
                                             output_dir = "SCI_Enrichment_Results",
                                             p_cutoff = 0.05,
                                             q_cutoff = 0.2,
                                             create_plots = TRUE) {
  
  # Step 1: Convert protein names to Entrez IDs
  id_mapping <- convert_protein_ids(protein_names)
  entrez_ids <- unique(id_mapping$ENTREZID[!is.na(id_mapping$ENTREZID)])
  
  if (length(entrez_ids) < 5) {
    warning("Insufficient proteins mapped (<5). Enrichment analysis may not be reliable.")
    return(NULL)
  }
  
  message(paste("Proceeding with", length(entrez_ids), "mapped proteins..."))
  
  # Step 2: Run enrichment analyses
  go_results <- run_go_enrichment(entrez_ids, 
                                  p_cutoff = p_cutoff, 
                                  q_cutoff = q_cutoff)
  
  kegg_result <- run_kegg_enrichment(entrez_ids, 
                                     organism = organism,
                                     p_cutoff = p_cutoff,
                                     q_cutoff = q_cutoff)
  
  # Step 3: Create SCI-quality figures
  if (create_plots) {
    message("\nCreating SCI-quality visualization figures...")
    
    figures <- create_comprehensive_enrichment_figure(
      go_results = go_results,
      kegg_result = kegg_result,
      output_dir = output_dir,
      study_title = "Proteomics Enrichment Analysis"
    )
    
    message("SCI-quality figures created successfully.")
  }
  
  # Step 4: Export results
  export_enrichment_results(id_mapping, go_results, kegg_result, output_dir)
  
  # Step 5: Summary

  message("ANALYSIS COMPLETE")
  message(strrep("=", 50))
  
  cat("\nSUMMARY STATISTICS:\n")
  cat(paste("Input proteins:", length(protein_names), "\n"))
  cat(paste("Mapped proteins:", length(entrez_ids), "\n"))
  cat(paste("Mapping rate:", round(length(entrez_ids)/length(protein_names)*100, 1), "%\n"))
  
  if (!is.null(go_results$BP)) {
    cat(paste("Significant GO BP terms:", sum(go_results$BP$p.adjust < 0.05), "\n"))
  }
  if (!is.null(go_results$CC)) {
    cat(paste("Significant GO CC terms:", sum(go_results$CC$p.adjust < 0.05), "\n"))
  }
  if (!is.null(go_results$MF)) {
    cat(paste("Significant GO MF terms:", sum(go_results$MF$p.adjust < 0.05), "\n"))
  }
  if (!is.null(kegg_result)) {
    cat(paste("Significant KEGG pathways:", sum(kegg_result$p.adjust < 0.05), "\n"))
  }
  
  cat(paste("\nResults saved in:", normalizePath(output_dir), "\n"))
  
  # Return comprehensive results
  return(list(
    id_mapping = id_mapping,
    go_results = go_results,
    kegg_result = kegg_result,
    entrez_ids = entrez_ids,
    figures = figures
  ))
}
