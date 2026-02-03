# Load necessary R packages
library(STRINGdb)
library(clusterProfiler)
library(org.Hs.eg.db)
library(igraph)
library(dplyr)
library(biomaRt)
library(ggplot2)
library(tidygraph)
library(ggraph)



validate_causal_mechanisms <- function(kocmi_network, ppi_network) {
  
  # 提取显著的因果边
  significant_causal_edges <- kocmi_network %>%
    filter(p_adj < 0.05 & abs(cs) > 0.1)
  
  # 在PPI网络中寻找直接物理证据
  physical_evidence <- ppi_network %>%
    filter(protein_A %in% significant_causal_edges$regulator &
             protein_B %in% significant_causal_edges$target)
  
  # 分类结果
  validation_categories <- list(
    strongly_supported = physical_evidence,  # 有直接物理证据
    potentially_indirect = setdiff(significant_causal_edges, physical_evidence)  # 可能通过中介
  )
  
  return(validation_categories)
}




# Improved protein mapping function
map_via_uniprot <- function(protein_names, species = "hsapiens") {
  
  tryCatch({
    ensembl <- useMart("ensembl", dataset = paste0(species, "_gene_ensembl"))
    
    uniprot_mapping <- getBM(
      attributes = c("hgnc_symbol", "uniprotswissprot", "uniprotsptrembl"),
      filters = "hgnc_symbol", 
      values = protein_names,
      mart = ensembl
    )
    
    # Handle empty values
    uniprot_mapping$uniprotswissprot[uniprot_mapping$uniprotswissprot == ""] <- NA
    uniprot_mapping$uniprotsptrembl[uniprot_mapping$uniprotsptrembl == ""] <- NA
    
    uniprot_mapping_clean <- uniprot_mapping %>%
      mutate(final_uniprot = ifelse(!is.na(uniprotswissprot), 
                                    uniprotswissprot, uniprotsptrembl)) %>%
      group_by(hgnc_symbol) %>%
      dplyr::slice(1) %>%  
      ungroup() %>%
      filter(!is.na(final_uniprot))  # Remove proteins that cannot be mapped
    
    cat("Successfully mapped", nrow(uniprot_mapping_clean), "proteins to UniProt\n")
    
    return(data.frame(uniprot_mapping_clean))
    
  }, error = function(e) {
    warning("UniProt mapping failed: ", e$message)
    return(NULL)
  })
}

# Improved enrichment analysis function
perform_enrichment_analysis <- function(mapped_proteins, 
                                        pvalue_cutoff = 0.05,
                                        qvalue_cutoff = 0.2,
                                        min_genes = 3) {
  
  # Get gene symbols
  gene_symbols <- na.omit(unique(mapped_proteins$hgnc_symbol))
  
  if(length(gene_symbols) < min_genes) {
    warning("Too few genes (", length(gene_symbols), "), skipping enrichment analysis")
    return(NULL)
  }
  
  # Convert to ENTREZ ID
  entrez_mapping <- tryCatch({
    clusterProfiler::bitr(gene_symbols, 
                          fromType = "SYMBOL", 
                          toType = "ENTREZID", 
                          OrgDb = org.Hs.eg.db)
  }, error = function(e) {
    warning("Gene symbol conversion failed: ", e$message)
    return(data.frame(SYMBOL = character(), ENTREZID = character()))
  })
  
  if(nrow(entrez_mapping) < min_genes) {
    warning("Too few successfully converted genes (", nrow(entrez_mapping), "), skipping enrichment analysis")
    return(NULL)
  }
  
  entrez_ids <- entrez_mapping$ENTREZID
  
  # Perform enrichment analysis
  enrichment_results <- list()
  
  # GO enrichment analysis
  tryCatch({
    enrichment_results$go_bp <- enrichGO(
      gene = entrez_ids,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = pvalue_cutoff,
      qvalueCutoff = qvalue_cutoff,
      minGSSize = min_genes,
      maxGSSize = 500
    )
    cat("GO enrichment analysis completed\n")
  }, error = function(e) {
    warning("GO enrichment analysis failed: ", e$message)
  })
  
  # KEGG enrichment analysis
  tryCatch({
    enrichment_results$kegg <- enrichKEGG(
      gene = entrez_ids,
      organism = 'hsa',
      pAdjustMethod = "BH",
      pvalueCutoff = pvalue_cutoff,
      qvalueCutoff = qvalue_cutoff,
      minGSSize = min_genes,
      maxGSSize = 500
    )
    cat("KEGG enrichment analysis completed\n")
  }, error = function(e) {
    warning("KEGG enrichment analysis failed: ", e$message)
  })
  
  return(list(
    proteins = gene_symbols,
    entrez_ids = entrez_ids,
    enrichment = enrichment_results
  ))
}

# Main PPI network analysis function - without module analysis
analyze_ppi_network_simple <- function(causal_proteins, 
                                       species = 9606, 
                                       score_threshold = 0.004,  # Lower threshold to get more connections
                                       timeout = 600) {
  
  # Set timeout
  options(timeout = timeout)
  
  cat("Starting analysis of PPI network containing", length(causal_proteins), "proteins...\n")
  
  # Initialize STRING database
  string_db <- STRINGdb$new(
    version = "11.5",
    species = species,
    score_threshold = score_threshold * 1000,  # Use more lenient threshold
    input_directory = ""
  )
  
  # Protein mapping
  uniprot_result <- map_via_uniprot(causal_proteins)
  
  if(is.null(uniprot_result) || nrow(uniprot_result) == 0) {
    stop("Protein mapping failed: no proteins successfully mapped")
  }
  
  cat("Successfully mapped", nrow(uniprot_result), "proteins to UniProt\n")
  
  # STRING ID mapping
  tryCatch({
    string_mapped <- string_db$map(
      uniprot_result, 
      "final_uniprot", 
      removeUnmappedRows = TRUE
    )
    cat("Successfully mapped", nrow(string_mapped), "proteins to STRING ID\n")
  }, error = function(e) {
    stop("STRING mapping failed: ", e$message)
  })
  
  # Build PPI network
  tryCatch({
    net_igraph <- string_db$get_subnetwork(string_mapped$STRING_id)
    cat("Constructed network contains", vcount(net_igraph), "nodes and", ecount(net_igraph), "edges\n")
    
    # If network has no edges, create minimal connected network for visualization
    if(ecount(net_igraph) == 0) {
      warning("No edge connections in network, creating minimal connected network for visualization")
      # Create a star network for visualization
      net_igraph <- make_star(vcount(net_igraph), mode = "undirected")
      V(net_igraph)$name <- string_mapped$STRING_id
    }
  }, error = function(e) {
    stop("Network construction failed: ", e$message)
  })
  
  # Overall enrichment analysis
  enrichment_results <- perform_enrichment_analysis(string_mapped, min_genes = 3)
  
  # Network statistics
  network_stats <- list(
    nodes = vcount(net_igraph),
    edges = ecount(net_igraph),
    density = edge_density(net_igraph),
    average_degree = mean(degree(net_igraph)),
    components = components(net_igraph)$no
  )
  
  # Return results
  return(list(
    string_db = string_db,
    network_igraph = net_igraph,
    enrichment_results = enrichment_results,
    mapped_proteins = string_mapped,
    network_stats = network_stats
  ))
}


# Identify fat-related pathways
identify_fat_pathways <- function(ppi_result) {
  
  fat_related_pathways <- list()
  
  tryCatch({
    enrichment <- ppi_result$enrichment_results
    
    if(is.null(enrichment)) {
      message("No enrichment analysis results")
      return(fat_related_pathways)
    }
    
    # Fat-related pathway keywords
    fat_keywords <- c("lipid", "fat", "adipose", "fatty", "obesity", 
                      "adipocytokine", "ppar", "insulin resistance", 
                      "metabolism", "leptin", "adiponectin", "cholesterol",
                      "triglyceride", "lipoprotein", "sterol")
    
    fat_pattern <- paste(fat_keywords, collapse = "|")
    
    # Check KEGG enrichment results
    if(!is.null(enrichment$enrichment$kegg)) {
      kegg_data <- enrichment$enrichment$kegg
      if(class(kegg_data)[1] == "enrichResult" && nrow(kegg_data@result) > 0) {
        kegg_terms <- kegg_data@result
        fat_terms <- kegg_terms[grepl(fat_pattern, kegg_terms$Description, ignore.case = TRUE), ]
        
        if(nrow(fat_terms) > 0) {
          fat_related_pathways[["kegg"]] <- list(
            pathways = fat_terms,
            fat_pathway_count = nrow(fat_terms),
            total_pathway_count = nrow(kegg_terms)
          )
          cat("✅ Found", nrow(fat_terms), "fat-related KEGG pathways\n")
        }
      }
    }
    
    # Check GO enrichment results
    if(!is.null(enrichment$enrichment$go_bp)) {
      go_data <- enrichment$enrichment$go_bp
      if(class(go_data)[1] == "enrichResult" && nrow(go_data@result) > 0) {
        go_terms <- go_data@result
        fat_terms <- go_terms[grepl(fat_pattern, go_terms$Description, ignore.case = TRUE), ]
        
        if(nrow(fat_terms) > 0) {
          fat_related_pathways[["go_bp"]] <- list(
            pathways = fat_terms,
            fat_pathway_count = nrow(fat_terms),
            total_pathway_count = nrow(go_terms)
          )
          cat("✅ Found", nrow(fat_terms), "fat-related GO biological processes\n")
        }
      }
    }
    
    return(fat_related_pathways)
    
  }, error = function(e) {
    message("Error identifying fat-related pathways: ", e$message)
    return(fat_related_pathways)
  })
}


# Generate enrichment analysis summary - simplified version
generate_enrichment_summary <- function(ppi_results) {
  
  enrichment <- ppi_results$enrichment_results
  
  if(is.null(enrichment)) {
    return(NULL)
  }
  
  summary_list <- list()
  
  # GO BP summary
  if(!is.null(enrichment$enrichment$go_bp) && nrow(enrichment$enrichment$go_bp@result) > 0) {
    go_summary <- enrichment$enrichment$go_bp@result %>%
      dplyr::arrange(p.adjust) %>%
      head(10) %>%  # Increase display count
      dplyr::mutate(type = "GO_BP") %>%
      dplyr::select(type, ID, Description, p.adjust, Count, geneID)
    
    summary_list[["GO_BP"]] <- go_summary
  }
  
  # KEGG summary
  if(!is.null(enrichment$enrichment$kegg) && nrow(enrichment$enrichment$kegg@result) > 0) {
    kegg_summary <- enrichment$enrichment$kegg@result %>%
      dplyr::arrange(p.adjust) %>%
      head(10) %>%  # Increase display count
      dplyr::mutate(type = "KEGG") %>%
      dplyr::select(type, ID, Description, p.adjust, Count, geneID)
    
    summary_list[["KEGG"]] <- kegg_summary
  }
  
  if(length(summary_list) > 0) {
    return(bind_rows(summary_list, .id = "database"))
  } else {
    return(NULL)
  }
}


ppi_network_vis <- function(ppi_results, plot_title = "PPI Network Visualization") {
  
  net <- ppi_results$network_igraph
  
  # Convert igraph object to tidygraph object for ggplot2 processing
  net_tidy <- as_tbl_graph(net)
  
  safe_cut <- function(x, breaks_probs = c(0, 0.33, 0.66, 1),
                       labels = c("Low", "Medium", "High")) {
    # Calculate quantiles
    breaks <- quantile(x, probs = breaks_probs, na.rm = TRUE)
    if (any(duplicated(breaks))) {
      # If there are duplicate quantiles, use equally spaced breakpoints
      breaks <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = length(breaks_probs))
    }
    cut(x, breaks = breaks, labels = labels, include.lowest = TRUE)
  }
  
  # 检查可用的边属性
  edge_attrs <- edge_attr_names(net)
  cat("Available edge attributes:", paste(edge_attrs, collapse = ", "), "\n")
  
  # 确定使用哪个属性作为边权重
  edge_weight_attr <- NULL
  possible_weight_attrs <- c("weight", "combined_score", "score", "interaction_score")
  
  for (attr in possible_weight_attrs) {
    if (attr %in% edge_attrs) {
      edge_weight_attr <- attr
      break
    }
  }
  
  # 如果有边权重属性，先添加到tidygraph对象中
  if (!is.null(edge_weight_attr)) {
    cat("Using", edge_weight_attr, "for edge weights\n")
    
    # 直接从igraph对象获取边权重
    edge_weights <- edge_attr(net, edge_weight_attr)
    
    # 将边权重添加到tidygraph对象
    net_tidy <- net_tidy %>%
      activate(edges) %>%
      mutate(edge_weight = edge_weights)
  }
  
  # Add node attributes
  net_tidy <- net_tidy %>%
    activate(nodes) %>%
    mutate(
      # Add gene symbols
      gene_symbol = ppi_results$mapped_proteins$hgnc_symbol[match(name, ppi_results$mapped_proteins$STRING_id)],
      # Calculate degree centrality
      degree = centrality_degree(),
      # Calculate betweenness centrality
      betweenness = centrality_betweenness(),
      # Set node size based on degree
      node_size = sqrt(degree) * 2 + 3,
      # Set node color based on betweenness centrality
      node_color = safe_cut(betweenness)
    )
  
  # Use ggraph to plot the network
  p <- ggraph(net_tidy, layout = 'fr')
  
  if (!is.null(edge_weight_attr)) {
    p <- p + 
      geom_edge_link(
        aes(width = edge_weight),
        color = "gray60",
        show.legend = TRUE
      ) +
      scale_edge_width_continuous(
        name = "Interaction\nStrength",
        range = c(0.5, 3)
      ) +
      scale_edge_alpha_continuous(
        name = "Interaction\nStrength",
        guide = "none"  # Hide alpha legend to avoid duplication
      )
  } else {
    cat("No edge weight attribute found. Using fixed edge width.\n")
    p <- p + 
      geom_edge_link(
        aes(alpha = after_stat(index)),
        color = "gray60",
        width = 0.5,
        show.legend = FALSE
      )
  }
  
  # Continue with the rest of the plot
  p <- p +
    # Draw nodes
    geom_node_point(
      aes(size = node_size, fill = node_color),
      shape = 21,
      color = "white",
      stroke = 0.5
    ) +
    # Add node labels
    geom_node_text(
      aes(label = ifelse(is.na(gene_symbol), name, gene_symbol)),
      size = 3,
      repel = TRUE,
      max.overlaps = 20,
      color = "black",
      fontface = "bold"
    ) +
    # Set color scale
    scale_fill_manual(
      values = c("Low" = "#4E79A7", "Medium" = "#F28E2B", "High" = "#E15759"),
      name = "Betweenness\nCentrality"
    ) +
    # Set size scale
    scale_size_continuous(
      name = "Degree",
      range = c(3, 8),
      guide = "none"  # Usually only show one in legend
    ) +
    # Theme and labels
    labs(
      title = plot_title,
      subtitle = paste("Network with", vcount(net), "nodes and", ecount(net), "edges")
    ) +
    theme_void(base_family = "serif") +
    theme(
      plot.title = element_text(
        size = 16,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 10)
      ),
      plot.subtitle = element_text(
        size = 12,
        hjust = 0.5,
        color = "gray50",
        margin = margin(b = 15)
      ),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 9),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  return(p)
}


convert_ppi_to_network_format <- function(ppi_results, kocmi_network) {
  
  # Extract network information from PPI results and convert to standard format
  ppi_igraph <- ppi_results$network_igraph
  mapped_proteins <- ppi_results$mapped_proteins
  
  # Check the column structure of mapped_proteins to determine how to get original gene names
  cat("Checking mapped_proteins columns:", paste(colnames(mapped_proteins), collapse = ", "), "\n")
  
  # Directly use hgnc_symbol as gene identifier
  if ("hgnc_symbol" %in% colnames(mapped_proteins)) {
    gene_col <- "hgnc_symbol"
    cat("Using 'hgnc_symbol' as gene identifier\n")
  } else {
    # If hgnc_symbol column is unexpectedly missing, use other available columns
    for (col in c("gene_symbol", "gene_name", "symbol", "final_uniprot")) {
      if (col %in% colnames(mapped_proteins)) {
        gene_col <- col
        cat("Using column '", gene_col, "' as gene identifier (hgnc_symbol not found)\n")
        break
      }
    }
  }  
  
  # Convert igraph network to edge list data frame
  if (vcount(ppi_igraph) > 0 && ecount(ppi_igraph) > 0) {
    # Get edge list
    edge_list <- as_data_frame(ppi_igraph, what = "edges")
    
    # Convert STRING IDs back to gene names
    edge_list <- edge_list %>%
      left_join(
        mapped_proteins %>% 
          dplyr::select(STRING_id, gene_name = !!sym(gene_col)),
        by = c("from" = "STRING_id")
      ) %>%
      dplyr::rename(protein_A = gene_name) %>%
      left_join(
        mapped_proteins %>% 
          dplyr::select(STRING_id, gene_name = !!sym(gene_col)),
        by = c("to" = "STRING_id")
      ) %>%
      dplyr::rename(protein_B = gene_name) %>%
      filter(!is.na(protein_A) & !is.na(protein_B)) %>%
      mutate(
        interaction_score = case_when(
          !is.na(combined_score) ~ combined_score / 1000,
          TRUE ~ 0.5),  # Default score
        interaction_type = ifelse(interaction_score > 0.3, "high_confidence", "Other")
      ) %>%
      dplyr::select(protein_A, protein_B, interaction_score, interaction_type)
  } else {
    # If network is empty, create empty data frame
    edge_list <- data.frame(
      protein_A = character(),
      protein_B = character(), 
      interaction_score = numeric(),
      interaction_type = character(),
      stringsAsFactors = FALSE
    )
    cat("Empty PPI network - no edges found\n")
  }
  
  # Standardize PPI network format - handle undirected nature
  ppi_standardized <- edge_list %>%
    standardize_protein_names() %>%
    # For PPI network, create edges in both directions to ensure matching with KOCMI
    bind_rows(
      edge_list %>%
        standardize_protein_names() %>%
        dplyr::rename(protein_A = protein_B, protein_B = protein_A)  # Add reverse edges
    ) %>%
    distinct(protein_A, protein_B, .keep_all = TRUE) %>%  # Remove duplicates
    mutate(
      edge_id = paste(protein_A, protein_B, sep = "|"),
      is_directed = FALSE  # Mark as undirected
    ) 
  
  ppi_proteins <- unique(c(ppi_standardized$protein_A, ppi_standardized$protein_B))
  cat("PPI network contains", length(ppi_proteins), "unique proteins\n")
  
  # Standardize KOCMI network format - maintain directionality
  kocmi_standardized <- kocmi_network %>%
    filter(!is.na(regulator) & !is.na(target)) %>%
    standardize_protein_names() %>%
    # Only keep protein pairs that exist in the PPI network
    filter(regulator %in% ppi_proteins & target %in% ppi_proteins) %>%
    mutate(
      edge_id = paste(regulator, target, sep = "|"),
      direction = ifelse(cs > 0, "activation", "inhibition"),
      is_directed = TRUE  # Mark as directed
    ) 
  
  # Create merged network - handle directional differences
  merged_network <- create_merged_network(kocmi_standardized, ppi_standardized)
  
  # Print network statistics
  cat("Converted PPI network:", nrow(ppi_standardized), "edges (bidirectional)\n")
  cat("KOCMI network:", nrow(kocmi_standardized), "directed edges\n")
  cat("Merged network:", nrow(merged_network$edges), "edges\n")
  cat("Overlapping edges:", sum(merged_network$edges$has_ppi_support & merged_network$edges$has_kocmi_support), "\n")
  
  return(list(
    kocmi = kocmi_standardized,
    ppi = ppi_standardized,
    merged = merged_network,
    ppi_raw = ppi_results,  # Keep original PPI results for reference
    conversion_info = list(
      gene_column_used = gene_col,
      mapped_proteins_sample = head(mapped_proteins, 5)
    )
  ))
}

# Helper function: Create merged network, handle directionality
create_merged_network <- function(kocmi_net, ppi_net) {
  
  # Create node list
  kocmi_nodes <- unique(c(kocmi_net$regulator, kocmi_net$target))
  ppi_nodes <- unique(c(ppi_net$protein_A, ppi_net$protein_B))
  all_nodes <- unique(c(kocmi_nodes, ppi_nodes))
  
  # Create merged edge list
  merged_edges <- bind_rows(
    # KOCMI edges
    kocmi_net %>%
      dplyr::select(from = regulator, to = target, 
                    kocmi_score = cs, kocmi_pvalue = p_adj, 
                    edge_id, is_directed) %>%
      mutate(
        has_kocmi_support = TRUE,
        has_ppi_support = edge_id %in% ppi_net$edge_id | 
          paste(to, from, sep = "|") %in% ppi_net$edge_id
      ),
    
    # PPI edges (only add edges not present in KOCMI)
    ppi_net %>%
      dplyr::select(from = protein_A, to = protein_B, 
                    ppi_score = interaction_score, 
                    ppi_type = interaction_type, edge_id, is_directed) %>%
      mutate(
        has_ppi_support = TRUE,
        has_kocmi_support = edge_id %in% kocmi_net$edge_id
      ) %>%
      dplyr::filter(!has_kocmi_support)  # Only add edges not present in KOCMI
  ) %>%
    # For edges that exist in both, merge information
    group_by(from, to) %>%
    summarise(
      has_kocmi_support = any(has_kocmi_support, na.rm = TRUE),
      has_ppi_support = any(has_ppi_support, na.rm = TRUE),
      kocmi_score = dplyr::first(na.omit(kocmi_score)),
      kocmi_pvalue = dplyr::first(na.omit(kocmi_pvalue)),
      ppi_score = dplyr::first(na.omit(ppi_score)),
      ppi_type = dplyr::first(na.omit(ppi_type)),
      is_directed = any(is_directed, na.rm = TRUE),  # If any network is directed, mark as directed
      .groups = "drop"
    ) %>%
    mutate(
      edge_id = paste(from, to, sep = "|"),
      # Calculate combined score
      combined_score = case_when(
        has_kocmi_support & has_ppi_support ~ (kocmi_score + ppi_score) / 2,
        has_kocmi_support ~ kocmi_score,
        has_ppi_support ~ ppi_score,
        TRUE ~ 0
      ),
      confidence_level = case_when(
        has_kocmi_support & has_ppi_support ~ "high",
        has_kocmi_support | has_ppi_support ~ "medium",
        TRUE ~ "low"
      )
    )
  
  # Create node attributes
  node_attributes <- data.frame(
    node_id = all_nodes,
    in_kocmi = all_nodes %in% kocmi_nodes,
    in_ppi = all_nodes %in% ppi_nodes,
    degree_kocmi = sapply(all_nodes, function(x) 
      sum(kocmi_net$regulator == x) + sum(kocmi_net$target == x)),
    degree_ppi = sapply(all_nodes, function(x) 
      sum(ppi_net$protein_A == x) + sum(ppi_net$protein_B == x)),
    stringsAsFactors = FALSE
  )
  
  return(list(
    edges = merged_edges,
    nodes = node_attributes,
    network_type = "mixed_directed_undirected"
  ))
}

# Helper function: Standardize protein names (if not already defined)
standardize_protein_names <- function(df) {
  # This should contain your protein name standardization logic
  # For example: convert to uppercase, remove version numbers, etc.
  if ("protein_A" %in% colnames(df)) {
    df <- df %>%
      mutate(
        protein_A = toupper(protein_A),
        protein_B = toupper(protein_B)
      )
  }
  if ("regulator" %in% colnames(df) & "target" %in% colnames(df)) {
    df <- df %>%
      mutate(
        regulator = toupper(regulator),
        target = toupper(target)
      )
  }
  return(df)
}


validate_causal_edges_with_ppi <- function(networks, p_threshold = 0.05, cs_threshold = 0.1) {
  
  # Filter significant causal edges
  significant_causal_edges <- networks$kocmi %>%
    filter(p_adj < p_threshold & abs(cs) > cs_threshold) %>%
    dplyr::select(regulator, target, cs, p_adj, direction, edge_id)
  
  # Find direct physical evidence in PPI network
  physical_evidence <- networks$ppi %>%
    filter(interaction_type == "high_confidence") %>%
    mutate(
      edge_id_forward = paste(protein_A, protein_B, sep = "|"),
      edge_id_reverse = paste(protein_B, protein_A, sep = "|")
    )
  
  # Match causal edges with physical interactions
  validated_edges <- significant_causal_edges %>%
    left_join(
      physical_evidence %>%
        filter(edge_id_forward %in% significant_causal_edges$edge_id | 
                 edge_id_reverse %in% significant_causal_edges$edge_id) %>%
        mutate(
          matched_edge = ifelse(edge_id_forward %in% significant_causal_edges$edge_id, 
                                edge_id_forward, edge_id_reverse)
        ) %>%
        dplyr::select(matched_edge, interaction_score, interaction_type),
      by = c("edge_id" = "matched_edge")
    ) %>%
    mutate(
      has_physical_evidence = !is.na(interaction_score),
      validation_status = case_when(
        has_physical_evidence & interaction_type == "high_confidence" ~ "strongly_supported",
        has_physical_evidence ~ "moderately_supported", 
        TRUE ~ "no_physical_evidence"
      )
    )
  
  # Statistical validation results
  validation_summary <- validated_edges %>%
    group_by(validation_status) %>%
    summarise(
      count = n(),
      avg_causal_strength = mean(abs(cs)),
      proportion = round(n() / nrow(validated_edges) * 100, 1)
    )
  
  cat("=== CAUSAL EDGE VALIDATION WITH PPI ===\n")
  print(validation_summary)
  
  return(list(
    validated_edges = validated_edges,
    summary = validation_summary
  ))
}


analyze_indirect_path_mechanisms <- function(networks, indirect_paths) {
  
  # Helper function: Check physical interactions
  check_physical_interaction <- function(ppi_net, protein1, protein2) {
    
    protein1_upper <- toupper(protein1)
    protein2_upper <- toupper(protein2)
    
    interaction <- ppi_net %>%
      filter(
        (protein_A == protein1_upper & protein_B == protein2_upper) |
          (protein_A == protein2_upper & protein_B == protein1_upper)
      )
    
    if (nrow(interaction) > 0) {
      return(list(
        exists = TRUE,
        interaction_score = interaction$interaction_score[1],
        interaction_type = interaction$interaction_type[1]
      ))
    } else {
      return(list(exists = FALSE, interaction_score = 0, interaction_type = "none"))
    }
  }
  
  # Helper function: Calculate mechanism confidence score - FIXED VERSION
  calculate_mechanism_confidence <- function(physical_links) {
    total_possible_links <- 0
    confirmed_links <- 0
    total_score <- 0
    
    # Check source to first mediator
    if (!is.null(physical_links$source_first_mediator) && 
        !is.null(physical_links$source_first_mediator$exists)) {
      total_possible_links <- total_possible_links + 1
      if (physical_links$source_first_mediator$exists) {
        confirmed_links <- confirmed_links + 1
        total_score <- total_score + physical_links$source_first_mediator$interaction_score
      }
    }
    
    # Check between mediators
    if (!is.null(physical_links$mediator_links)) {
      for (link_name in names(physical_links$mediator_links)) {
        link <- physical_links$mediator_links[[link_name]]
        if (!is.null(link$exists)) {
          total_possible_links <- total_possible_links + 1
          if (link$exists) {
            confirmed_links <- confirmed_links + 1
            total_score <- total_score + link$interaction_score
          }
        }
      }
    }
    
    # Check last mediator to target
    if (!is.null(physical_links$last_mediator_target) && 
        !is.null(physical_links$last_mediator_target$exists)) {
      total_possible_links <- total_possible_links + 1
      if (physical_links$last_mediator_target$exists) {
        confirmed_links <- confirmed_links + 1
        total_score <- total_score + physical_links$last_mediator_target$interaction_score
      }
    }
    
    if (total_possible_links == 0) return(0)
    
    # Calculate weighted confidence score
    connection_ratio <- confirmed_links / total_possible_links
    avg_score <- ifelse(confirmed_links > 0, total_score / confirmed_links, 0)
    
    # Return weighted score (connection ratio * average interaction score)
    return(connection_ratio * avg_score)
  }
  
  # Helper function: Interpret confidence score
  interpret_mechanism_confidence <- function(score) {
    if (score >= 0.7) return("High confidence: Strong physical evidence supports the causal path")
    if (score >= 0.4) return("Medium confidence: Partial physical evidence exists")
    if (score >= 0.1) return("Low confidence: Limited physical evidence")
    return("Very low confidence: Little to no physical evidence")
  }
  
  # Main analysis logic
  path_analyses <- list()
  
  for (i in 1:nrow(indirect_paths)) {
    path <- indirect_paths[i,]
    
    # Extract proteins from the path
    mediators <- unlist(path$mediators)
    path_proteins <- c(path$source, mediators, path$target)
    
    # Check physical interactions between adjacent proteins in the path
    physical_links <- list()
    
    # Check source protein with first mediator
    if (length(path$mediators) > 0) {
      first_mediator <- mediators[1]
      physical_links$source_first_mediator <- check_physical_interaction(
        networks$ppi, path$source, first_mediator
      )
    }
    
    # Check interactions between mediators
    mediator_links <- list()
    if (length(mediators) > 1) {
      for (j in 1:(length(mediators) - 1)) {
        mediator1 <- mediators[j]
        mediator2 <- mediators[j + 1]
        link_name <- paste("mediator", j, j+1, sep = "_")
        mediator_links[[link_name]] <- check_physical_interaction(
          networks$ppi, mediator1, mediator2)
      }
    }
    
    physical_links$mediator_links <- mediator_links
    
    # Check last mediator with target
    if (length(mediators) > 0) {
      last_mediator <- mediators[length(mediators)]
      physical_links$last_mediator_target <- check_physical_interaction(
        networks$ppi, last_mediator, path$target
      )
    }
    
    # Calculate mechanism confidence score
    confidence_score <- calculate_mechanism_confidence(physical_links)
    
    path_analyses[[i]] <- list(
      path_index = i,
      path_description = path$full_path,
      path_strength = path$path_strength,
      mediators = mediators,
      physical_links = physical_links,
      mechanism_confidence = confidence_score,
      interpretation = interpret_mechanism_confidence(confidence_score)
    )
  }
  
  # Sort by confidence score
  path_analyses <- path_analyses[order(-sapply(path_analyses, function(x) x$mechanism_confidence))]
  
  return(path_analyses)
}

visualize_indirect_path_analysis <- function(path_analyses, top_n = 10) {
  
  # 1. 提取数据用于可视化
  analysis_df <- do.call(rbind, lapply(path_analyses, function(x) {
    data.frame(
      path_index = x$path_index,
      path_description = x$path_description,
      path_strength = x$path_strength,
      mechanism_confidence = x$mechanism_confidence,
      interpretation = x$interpretation,
      n_mediators = length(strsplit(x$path_description, " -> ")[[1]]) - 2,
      stringsAsFactors = FALSE
    )
  }))
  
  # 2. 置信度分布图
  p1 <- ggplot(analysis_df, aes(x = mechanism_confidence, fill = interpretation)) +
    geom_histogram(bins = 20, alpha = 0.7) +
    scale_fill_manual(values = c("low" = "#E41A1C", "medium" = "#377EB8", "high" = "#4DAF4A")) +
    labs(
      title = "Distribution of Mechanism Confidence Scores",
      x = "Mechanism Confidence Score",
      y = "Count",
      fill = "Confidence Level"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # 3. 路径强度 vs 置信度散点图
  p2 <- ggplot(analysis_df, aes(x = path_strength, y = mechanism_confidence, 
                                color = interpretation, size = n_mediators)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("low" = "#E41A1C", "medium" = "#377EB8", "high" = "#4DAF4A")) +
    labs(
      title = "Path Strength vs Mechanism Confidence",
      x = "Path Strength",
      y = "Mechanism Confidence",
      color = "Confidence Level",
      size = "Number of Mediators"
    ) +
    theme_minimal()
  
  # 4. 前N条路径的置信度条形图
  top_paths <- head(analysis_df[order(-analysis_df$mechanism_confidence), ], top_n)
  
  p3 <- ggplot(top_paths, aes(x = reorder(path_description, mechanism_confidence), 
                              y = mechanism_confidence, fill = interpretation)) +
    geom_col() +
    scale_fill_manual(values = c("low" = "#E41A1C", "medium" = "#377EB8", "high" = "#4DAF4A")) +
    labs(
      title = paste("Top", top_n, "Paths by Mechanism Confidence"),
      x = "Path",
      y = "Mechanism Confidence Score"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    coord_flip()
  
  # 5. 单个高置信度路径的详细网络图
  if (length(path_analyses) > 0) {
    best_path <- path_analyses[[1]]  # 置信度最高的路径
    p4 <- visualize_single_path_network(best_path)
  } else {
    p4 <- ggplot() + 
      annotate("text", x = 1, y = 1, label = "No paths to visualize") +
      theme_void()
  }
  
  # 6. 组合所有图形
  combined_plot <- (p1 | p2) / (p3 | p4) +
    plot_annotation(title = "Indirect Path Mechanism Analysis Results",
                    theme = theme(plot.title = element_text(size = 16, face = "bold")))
  
  print(combined_plot)
  
  return(list(
    summary_plot = combined_plot,
    confidence_distribution = p1,
    strength_vs_confidence = p2,
    top_paths_barplot = p3,
    best_path_network = p4,
    analysis_data = analysis_df
  ))
}

# 可视化单个路径的网络图
visualize_single_path_network <- function(path_analysis) {
  
  # 提取路径信息
  path_elements <- strsplit(path_analysis$path_description, " -> ")[[1]]
  source_node <- path_elements[1]
  mediators <- path_elements[2:(length(path_elements)-1)]
  target_node <- path_elements[length(path_elements)]
  
  # 创建节点数据
  nodes <- data.frame(
    name = c(source_node, mediators, target_node),
    type = c("source", rep("mediator", length(mediators)), "target"),
    stringsAsFactors = FALSE
  )
  
  # 创建边数据（因果边）
  edges <- data.frame(
    from = c(source_node, mediators),
    to = c(mediators, target_node),
    type = "causal",
    stringsAsFactors = FALSE
  )
  
  # 添加物理相互作用边
  physical_edges <- data.frame(
    from = character(),
    to = character(),
    type = character(),
    stringsAsFactors = FALSE
  )
  
  # 检查源与第一个中介的物理相互作用
  if (!is.null(path_analysis$physical_links$source_first_mediator) &&
      path_analysis$physical_links$source_first_mediator) {
    physical_edges <- rbind(physical_edges, 
                            data.frame(from = source_node, to = mediators[1], type = "physical"))
  }
  
  # 检查中介之间的物理相互作用
  if (length(mediators) > 1 && !is.null(path_analysis$physical_links$mediator_links)) {
    for (j in 1:(length(mediators) - 1)) {
      link_name <- paste("mediator", j, j+1, sep = "_")
      if (!is.null(path_analysis$physical_links$mediator_links[[link_name]]) &&
          path_analysis$physical_links$mediator_links[[link_name]]) {
        physical_edges <- rbind(physical_edges,
                                data.frame(from = mediators[j], to = mediators[j+1], type = "physical"))
      }
    }
  }
  
  # 创建igraph对象
  g <- graph_from_data_frame(edges, directed = TRUE, vertices = nodes)
  
  # 添加物理相互作用边（无向）
  if (nrow(physical_edges) > 0) {
    for (i in 1:nrow(physical_edges)) {
      g <- add_edges(g, c(physical_edges$from[i], physical_edges$to[i]), 
                     type = "physical")
    }
  }
  
  # 设置节点和边的可视化属性
  V(g)$color <- ifelse(V(g)$type == "source", "#E41A1C", 
                       ifelse(V(g)$type == "target", "#4DAF4A", "#377EB8"))
  V(g)$size <- ifelse(V(g)$type == "source", 10, 
                      ifelse(V(g)$type == "target", 10, 8))
  
  E(g)$color <- ifelse(E(g)$type == "causal", "black", "red")
  E(g)$linetype <- ifelse(E(g)$type == "causal", "solid", "dashed")
  E(g)$arrow.size <- ifelse(E(g)$type == "causal", 0.5, 0)
  
  # 绘制网络图
  p <- ggraph(g, layout = "linear") +
    geom_edge_arc(aes(color = color, linetype = linetype),
                  arrow = arrow(length = unit(2, 'mm'), type = "closed"),
                  strength = 0.2) +
    geom_node_point(aes(color = color, size = size)) +
    geom_node_text(aes(label = name), repel = TRUE, size = 3) +
    scale_edge_color_identity() +
    scale_color_identity() +
    scale_size_identity() +
    labs(title = paste("Best Path:", path_analysis$path_description),
         subtitle = paste("Confidence:", round(path_analysis$mechanism_confidence, 3))) +
    theme_void() +
    theme(plot.title = element_text(size = 12, face = "bold"),
          plot.subtitle = element_text(size = 10))
  
  return(p)
}

# 创建物理相互作用检查函数（如果尚未定义）
check_physical_interaction <- function(ppi_network, protein1, protein2) {
  # 检查两个蛋白质之间是否存在物理相互作用
  forward_id <- paste(protein1, protein2, sep = "|")
  reverse_id <- paste(protein2, protein1, sep = "|")
  
  interaction_exists <- any(ppi_network$edge_id == forward_id | 
                              ppi_network$edge_id == reverse_id)
  
  return(interaction_exists)
}

# 计算机制置信度函数（如果尚未定义）
calculate_mechanism_confidence <- function(physical_links) {
  total_links <- 0
  confirmed_links <- 0
  
  # 计算所有可能的连接
  if (!is.null(physical_links$source_first_mediator)) {
    total_links <- total_links + 1
    if (physical_links$source_first_mediator) confirmed_links <- confirmed_links + 1
  }
  
  if (!is.null(physical_links$mediator_links)) {
    total_links <- total_links + length(physical_links$mediator_links)
    confirmed_links <- confirmed_links + sum(unlist(physical_links$mediator_links))
  }
  
  if (total_links == 0) return(0)
  
  return(confirmed_links / total_links)
}

# 解释置信度函数（如果尚未定义）
interpret_mechanism_confidence <- function(confidence_score) {
  if (confidence_score >= 0.7) return("high")
  if (confidence_score >= 0.4) return("medium")
  return("low")
}


networks_target <- convert_ppi_to_network_format(ppi_results = ppi_result,
                                                 kocmi_network = res_kocmi)
validation_results <- validate_causal_edges_with_ppi(networks = networks_target)
plot_validation_results(validation_results)
View(validation_results$validated_edges)

resm <- analyze_indirect_path_mechanisms(networks = networks_target,
                                         indirect_paths = path_res$indirect_paths[1:10,])
visualization_results <- visualize_indirect_path_analysis(resm, top_n = 10)
