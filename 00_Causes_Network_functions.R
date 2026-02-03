analyze_pathways <- function(weightNet, outcome_var, 
                             p_threshold = 0.05, 
                             cs_threshold = 0.1,
                             max_path_length = 3,
                             path_strength_threshold = 0.05) {
  library(igraph)
  library(dplyr)
  
  # Filter significant edges
  sig_edges <- weightNet %>% 
    dplyr::filter(p_adj < p_threshold & cs >= cs_threshold)
  
  if(nrow(sig_edges) == 0) {
    warning("No significant edges found for visualization.")
    return(NULL)
  }
  
  # Create graph
  g <- graph_from_data_frame(sig_edges[, c("regulator", "target")], directed = TRUE)
  E(g)$weight <- sig_edges$cs
  E(g)$p_adj <- sig_edges$p_adj
  
  # Check if outcome variable exists
  if(!outcome_var %in% V(g)$name) {
    warning(paste("Outcome variable", outcome_var, "not found in the network."))
    return(NULL)
  }
  
  # Find direct regulators
  direct_regulators <- neighbors(g, outcome_var, mode = "in")$name
  direct_edges <- sig_edges %>% 
    dplyr::filter(target == outcome_var) %>%
    arrange(desc(cs))
  
  # Find indirect paths
  indirect_paths <- find_all_indirect_paths(
    weightNet,
    target_nodes = outcome_var,
    max_path_length = max_path_length,
    min_strength = path_strength_threshold,
    p_threshold = p_threshold,
    cs_threshold = cs_threshold
  )
  
  # Select top indirect paths
  if(nrow(indirect_paths) > 0) {
    top_indirect_paths <- indirect_paths[order(-indirect_paths$path_strength), ]
  } else {
    top_indirect_paths <- data.frame()
  }
  
  # Prepare nodes for visualization
  all_nodes <- unique(c(
    outcome_var,
    direct_edges$regulator,
    unlist(top_indirect_paths$source),
    unlist(top_indirect_paths$mediators)
  ))
  
  # Create subgraph
  sub_g <- induced_subgraph(g, intersect(all_nodes, V(g)$name))
  
  # Calculate node properties
  node_degree <- degree(sub_g, mode = "in")
  betweenness_centrality <- betweenness(sub_g)
  
  # Create node categories
  node_categories <- data.frame(
    id = V(sub_g)$name,
    category = "other",
    stringsAsFactors = FALSE
  )
  
  # Mark direct regulators
  node_categories$category[node_categories$id %in% direct_edges$regulator] <- "direct"
  
  # Mark indirect sources
  indirect_sources <- setdiff(unique(top_indirect_paths$source), direct_edges$regulator)
  node_categories$category[node_categories$id %in% indirect_sources] <- "indirect_source"
  
  # Mark mediators
  all_mediators <- unique(unlist(top_indirect_paths$mediators))
  mediators_not_direct <- setdiff(all_mediators, direct_edges$regulator)
  node_categories$category[node_categories$id %in% mediators_not_direct] <- "mediator"
  
  # Mark outcome
  node_categories$category[node_categories$id == outcome_var] <- "outcome"
  
  V(sub_g)$category <- node_categories$category
  V(sub_g)$betweenness <- betweenness_centrality
  
  # Create node data frame for ggraph
  nodes_df <- data.frame(
    name = V(sub_g)$name,
    category = node_categories$category,
    degree = node_degree,
    betweenness = betweenness_centrality,
    stringsAsFactors = FALSE
  )
  
  # 修复：正确创建边数据框
  edges_df <- data.frame(
    from = character(),
    to = character(),
    weight = numeric(),
    p_adj = numeric(),
    stringsAsFactors = FALSE
  )
  
  # 从子图中提取边信息
  if (ecount(sub_g) > 0) {
    # 获取边的端点
    edge_list <- get.edgelist(sub_g)
    
    # 创建边数据框
    edges_df <- data.frame(
      from = edge_list[, 1],
      to = edge_list[, 2],
      weight = E(sub_g)$weight,
      p_adj = E(sub_g)$p_adj,
      stringsAsFactors = FALSE
    )
    
    # 从原始sig_edges中补充信息
    for(i in 1:nrow(edges_df)) {
      edge_data <- sig_edges %>% 
        filter(regulator == edges_df$from[i] & target == edges_df$to[i])
      if(nrow(edge_data) > 0) {
        edges_df$weight[i] <- edge_data$cs[1]
        edges_df$p_adj[i] <- edge_data$p_adj[1]
      }
    }
  }
  
  # Create summary statistics
  summary_stats <- list(
    outcome = outcome_var,
    total_nodes = length(V(sub_g)),
    direct_regulators = length(direct_regulators),
    indirect_paths_found = nrow(indirect_paths),
    top_direct_edges = direct_edges,
    top_indirect_paths = top_indirect_paths
  )
  
  # Return all analysis results
  return(list(
    graph = sub_g,
    nodes_df = nodes_df,
    edges_df = edges_df,
    direct_edges = direct_edges,
    indirect_paths = top_indirect_paths,
    summary = summary_stats,
    sig_edges = sig_edges,
    node_categories = node_categories
  ))
}

visualize_pathways <- function(analysis_results, 
                               top_direct = 10,
                               top_indirect = 15,
                               plot_theme = "classic") {
  library(ggplot2)
  library(ggraph)
  library(patchwork)
  library(gridExtra)
  library(dplyr)
  library(tidyr)
  library(igraph)
  library(ggrepel)
  
  category_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#999999")
  category_shapes <- c(18, 15, 17, 19, 19)
  
  # Check if analysis results are valid
  if(is.null(analysis_results)) {
    stop("No analysis results provided for visualization.")
  }
  
  # Extract components from analysis results
  sub_g <- analysis_results$graph
  
  category_config <- list(
    outcome = list(color = "#E41A1C", shape = 18, label = "Outcome"),
    direct = list(color = "#377EB8", shape = 15, label = "Direct Regulator"),
    indirect_source = list(color = "#4DAF4A", shape = 17, label = "Indirect Source"),
    mediator = list(color = "#984EA3", shape = 19, label = "Mediator"),
    other = list(color = "#999999", shape = 19, label = "Other")
  )
  
  # 提取颜色和形状向量
  category_colors <- sapply(category_config, function(x) x$color)
  category_shapes <- sapply(category_config, function(x) x$shape)
  category_labels <- sapply(category_config, function(x) x$label)
  
  # 安全地提取数据
  nodes_df <- analysis_results$nodes_df
  edges_df <- analysis_results$edges_df
  
  # 确保我们不超过实际可用的行数
  if(!is.null(analysis_results$direct_edges) && nrow(analysis_results$direct_edges) > 0) {
    n_direct <- min(top_direct, nrow(analysis_results$direct_edges))
    direct_edges <- arrange(analysis_results$direct_edges, desc(cs))[1:n_direct,]
  } else {
    direct_edges <- data.frame()
  }
  
  if(!is.null(analysis_results$indirect_paths) && nrow(analysis_results$indirect_paths) > 0) {
    n_indirect <- min(top_indirect, nrow(analysis_results$indirect_paths))
    top_indirect_paths <- arrange(analysis_results$indirect_paths, desc(path_strength))[1:n_indirect,]
  } else {
    top_indirect_paths <- data.frame()
  }
  
  summary_stats <- analysis_results$summary
  outcome_var <- ifelse(!is.null(summary_stats$outcome), summary_stats$outcome, "Outcome")
  
  # 从分析结果中提取边权重信息
  # 创建权重查找表
  if(!is.null(edges_df) && nrow(edges_df) > 0) {
    edge_weights <- edges_df %>%
      dplyr::select(from, to, weight) %>%
      rename(cs = weight)
  } else {
    # 从graph对象中提取
    if(!is.null(sub_g)) {
      edge_weights <- as.data.frame(igraph::as_edgelist(sub_g))
      colnames(edge_weights) <- c("from", "to")
      if("weight" %in% edge_attr_names(sub_g)) {
        edge_weights$cs <- E(sub_g)$weight
      } else {
        edge_weights$cs <- 1  # 默认值
      }
      edge_weights$weight <- edge_weights$cs
    } else {
      edge_weights <- data.frame()
    }
  }
  
  # Plot 1: Static network visualization
  network_plot <- ggraph(sub_g, layout = "fr") +
    geom_edge_fan(
      aes(width = weight, alpha = weight),
      arrow = arrow(length = unit(2, 'mm'), type = "closed"),
      end_cap = circle(3, 'mm'),
      color = "gray50",
      show.legend = TRUE) +
    geom_node_point(
      aes(color = category, shape = category, size = betweenness),
      alpha = 0.8) +
    geom_node_text(
      aes(label = name),
      repel = TRUE,
      size = 3,
      max.overlaps = 20) +
    scale_color_manual(
      name = "Node Type",
      values = category_colors, labels = category_labels) +
    scale_shape_manual(
      name = "Node Type", 
      values = category_shapes, labels = category_labels) +
    scale_size_continuous(name = "Betweenness Centrality") +
    scale_edge_width_continuous(name = "Causal Strength", range = c(0.5, 2)) +
    scale_edge_alpha_continuous(name = "Causal Strength", range = c(0.3, 1)) +
    theme_void() +
    labs(
      title = paste("Causal Network for", outcome_var),
      subtitle = paste("Direct regulators:", 
                       ifelse(is.null(summary_stats$direct_regulators), 0, summary_stats$direct_regulators), 
                       "| Indirect paths:", 
                       ifelse(is.null(summary_stats$indirect_paths_found), 0, summary_stats$indirect_paths_found))) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "bottom") +
    guides(
      color = guide_legend(
        nrow = 2,
        byrow = TRUE,
        title.position = "top",
        title.hjust = 0.5),
      shape = guide_legend(
        nrow = 2,
        byrow = TRUE,
        title.position = "top",
        title.hjust = 0.5),
      size = guide_legend(
        nrow = 1,
        title.position = "top",
        title.hjust = 0.5),
      edge_width = guide_legend(
        nrow = 1,
        title.position = "top",
        title.hjust = 0.5),
      edge_alpha = guide_legend(
        nrow = 1,
        title.position = "top",
        title.hjust = 0.5))
  
  # Plot 2: Direct effects bar plot
  if(nrow(direct_edges) > 0) {
    direct_effects_plot <- ggplot(direct_edges, aes(x = reorder(regulator, cs), y = cs)) +
      geom_col(aes(fill = cs), width = 0.7) +
      #geom_text(aes(label = format(cs, scientific = TRUE, digits = 2)), hjust = -0.1, size = 3) +
      scale_fill_gradient(low = "lightblue", high = "darkblue", 
                          name = "Causal Strength",guide = "none") +
      coord_flip() +
      labs(
        title = "Direct Effects on Outcome",
        x = "Regulator",
        y = "Causal Strength (CS)") +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold"),
        axis.text.y = element_text(size = 9))
  } else {
    direct_effects_plot <- ggplot() + 
      annotate("text", x = 1, y = 1, label = "No direct effects found") +
      theme_void()
  }
  
  # Plot 3: Indirect pathways strength
  if(nrow(top_indirect_paths) > 0) {
    # Prepare data for indirect pathways plot
    indirect_plot_data <- top_indirect_paths %>%
      arrange(path_strength) %>%
      mutate(
        path_label = substr(full_path, 1, 50),  # Truncate long paths
        path_label = factor(path_label, levels = path_label))
    
    indirect_pathways_plot <- ggplot(indirect_plot_data, 
                                     aes(x = path_label, y = path_strength)) +
      geom_col(aes(fill = path_strength), width = 0.7) +
      #geom_text(aes(label = path_length), hjust = -0.2, size = 3) +
      scale_fill_gradient(low = "lightblue", high = "darkblue", 
                          name = "Path Strength",guide = "none") +
      coord_flip() +
      labs(
        title = "Top Indirect Pathways",
        x = "Path (Source → ... → Outcome)",
        y = "Path Strength") +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold"),
        axis.text.y = element_text(size = 8))
  } else {
    indirect_pathways_plot <- ggplot() + 
      annotate("text", x = 1, y = 1, label = "No indirect pathways found") +
      theme_void()
  }
  
  # Plot 4: Pathway length distribution
  if(nrow(top_indirect_paths) > 0) {
    path_length_dist <- ggplot(top_indirect_paths, aes(x = factor(path_length))) +
      geom_bar(aes(fill = factor(path_length)), alpha = 0.7) +
      geom_text(stat = 'count', aes(label = ..count..), vjust = -0.5) +
      scale_fill_brewer(palette = "Set2", name = "Path Length") +
      labs(
        title = "Distribution of Pathway Lengths",
        x = "Path Length",
        y = "Count") +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold"))
  } else {
    path_length_dist <- ggplot() + 
      annotate("text", x = 1, y = 1, label = "No pathway length data") +
      theme_void()
  }
  
  # 5. Mediator Importance Analysis
  mediator_analysis_result <- NULL
  if(nrow(top_indirect_paths) > 0) {
    # Calculate importance of each mediator node
    mediator_importance <- top_indirect_paths %>%
      tidyr::unnest(mediators) %>%
      group_by(mediators) %>%
      summarise(
        n_paths = n(),
        avg_strength = mean(path_strength),
        total_strength = sum(path_strength)
      ) %>%
      arrange(desc(total_strength)) 
    
    if(nrow(mediator_importance) > 0) {
      mediator_importance_plot <- ggplot(mediator_importance, 
                                         aes(x = n_paths, y = avg_strength, color = total_strength)) +
        geom_point(alpha = 0.7) +
        ggrepel::geom_text_repel(aes(label = mediators), size = 4, max.overlaps = 20) +
        scale_color_gradient(low = "#3498DB", high = "#E74C3C", 
                             name = "Total Path Strength",guide = "none") +
        labs(
          title = "Mediator Node Importance Analysis",
          x = "Number of Involved Paths",
          y = "Average Strength") +
        theme_minimal() +
        theme(
          plot.title = element_text(face = "bold", size = 14),
          legend.position = "right")
      
      mediator_analysis_result <- list(
        plot = mediator_importance_plot,
        mediator_stats = mediator_importance)
    }
  }
  
  hierarchy_plot <- NULL
  if(nrow(top_indirect_paths) > 0 || nrow(direct_edges) > 0) {
    # 创建层次数据
    hierarchy_data <- data.frame(
      node = character(),
      level = numeric(),
      type = character(),
      stringsAsFactors = FALSE
    )
    
    hierarchy_data <- rbind(hierarchy_data, 
                            data.frame(node = outcome_var, level = 0, type = "outcome"))

    if(nrow(direct_edges) > 0) {
      hierarchy_data <- rbind(hierarchy_data,
                              data.frame(node = direct_edges$regulator, level = 1, type = "direct"))
    }

    if(nrow(top_indirect_paths) > 0) {
      all_mediators <- unique(unlist(top_indirect_paths$mediators))
      all_sources <- unique(top_indirect_paths$source)
      
      hierarchy_data <- rbind(hierarchy_data,
                              data.frame(node = all_mediators, level = 1, type = "mediator"))
      hierarchy_data <- rbind(hierarchy_data,
                              data.frame(node = all_sources, level = 2, type = "source"))
    }
    
    hierarchy_data <- distinct(hierarchy_data, node, .keep_all = TRUE)
    
    node_importance <- hierarchy_data %>%
      left_join(
        top_indirect_paths %>%
          tidyr::unnest(mediators) %>%
          group_by(mediators) %>%
          summarise(n_paths = n()) %>%
          dplyr::rename(node = mediators) %>%
          bind_rows(
            if(nrow(direct_edges) > 0) {
              direct_edges %>%
                dplyr::select(node = regulator) %>%
                mutate(n_paths = 1)
            } else {
              data.frame(node = character(), n_paths = numeric())
            }
          ) %>%
          bind_rows(
            top_indirect_paths %>%
              group_by(source) %>%
              summarise(n_paths = n()) %>%
              dplyr::rename(node = source)
          ) %>%
          group_by(node) %>%
          summarise(n_paths = sum(n_paths)),
        by = "node"
      ) %>%
      mutate(
        n_paths = ifelse(is.na(n_paths), 1, n_paths),
        size = case_when(
          type == "outcome" ~ 8,
          type == "direct" ~ 5 + n_paths * 0.5,
          type == "mediator" ~ 4 + n_paths * 0.5,
          type == "source" ~ 3 + n_paths * 0.3,
          TRUE ~ 4
        )
      )
    
    connections <- data.frame(
      from = character(),
      to = character(),
      cs = numeric(),
      stringsAsFactors = FALSE
    )
  
    if(nrow(direct_edges) > 0) {
      connections <- rbind(connections,
                           data.frame(from = direct_edges$regulator, 
                                      to = outcome_var, 
                                      cs = direct_edges$cs))
    }
    
    if(nrow(top_indirect_paths) > 0 && nrow(edge_weights) > 0) {
      for(i in 1:nrow(top_indirect_paths)) {
        path <- top_indirect_paths[i, ]
        mediators <- unlist(path$mediators)
        source <- path$source
        
        for(mediator in mediators) {
          # 查找中介节点到结果节点的CS值
          mediator_to_outcome_cs <- edge_weights %>%
            filter(from == mediator & to == outcome_var) %>%
            pull(cs)
          
          if(length(mediator_to_outcome_cs) > 0) {
            cs_value <- mean(mediator_to_outcome_cs)  # 如果有多个，取平均值
          } else {
            cs_value <- 0.1  # 默认值
          }
          
          connections <- rbind(connections,
                               data.frame(from = mediator, 
                                          to = outcome_var, 
                                          cs = cs_value))
        }
        
        # 为源节点到每个中介节点的边添加CS值
        for(mediator in mediators) {
          # 查找源节点到中介节点的CS值
          source_to_mediator_cs <- edge_weights %>%
            filter(from == source & to == mediator) %>%
            pull(cs)
          
          if(length(source_to_mediator_cs) > 0) {
            cs_value <- mean(source_to_mediator_cs)  # 如果有多个，取平均值
          } else {
            cs_value <- 0.1  # 默认值
          }
          
          connections <- rbind(connections,
                               data.frame(from = source, 
                                          to = mediator, 
                                          cs = cs_value))
        }
      }
    }
    
    # 移除重复连接，对重复边的CS值取平均值
    connections <- connections %>%
      group_by(from, to) %>%
      summarise(cs = mean(cs, na.rm = TRUE), .groups = "drop")
    
    # 创建图对象
    if(nrow(connections) > 0) {
      graph <- graph_from_data_frame(connections, vertices = node_importance)
      
      # 创建层次布局
      create_hierarchy_layout <- function(node_importance, outcome_var) {
        n_nodes <- nrow(node_importance)
        layout_matrix <- matrix(0, nrow = n_nodes, ncol = 2)
        
        levels <- sort(unique(node_importance$level))
        max_level <- max(levels)
        
        for(i in 1:n_nodes) {
          node <- node_importance$node[i]
          level <- node_importance$level[i]
          type <- node_importance$type[i]
          
          level_nodes <- which(node_importance$level == level)
          level_node_count <- length(level_nodes)
          
          if(level_node_count == 1) {
            x_pos <- 0
          } else {
            node_idx_in_level <- which(level_nodes == i)
            x_pos <- -1 + 2 * ((node_idx_in_level - 1) / (level_node_count - 1))
          }
          
          y_pos <- max_level - level
          
          layout_matrix[i, 1] <- x_pos
          layout_matrix[i, 2] <- y_pos
        }
        
        outcome_idx <- which(node_importance$node == outcome_var)
        if(length(outcome_idx) > 0) {
          layout_matrix[outcome_idx, 1] <- 0
          layout_matrix[outcome_idx, 2] <- max_level + 1.5
        }
        return(layout_matrix)
      }
      
      layout_matrix <- create_hierarchy_layout(node_importance, outcome_var)
      
      node_importance_with_coords <- node_importance %>%
        mutate(
          x = layout_matrix[, 1],
          y = layout_matrix[, 2]
        )
      
      # 创建带坐标的图
      graph_with_coords <- graph_from_data_frame(connections, 
                                                 vertices = node_importance_with_coords)
      
      # 创建hierarchy plot，线条颜色根据CS值设置
      hierarchy_plot <- ggraph(graph_with_coords, layout = "manual", 
                               x = layout_matrix[,1], y = layout_matrix[,2]) +
        geom_edge_link(
          arrow = arrow(length = unit(2, 'mm'), type = "closed"),
          end_cap = circle(3, 'mm'),
          aes(edge_width = cs, alpha = cs),  # 只用宽度和透明度
          color = "#2C3E50",  # 固定颜色 - 深蓝色
          show.legend = TRUE
        ) +
        geom_node_point(
          aes(color = type, size = size),
          alpha = 0.8
        ) +
        geom_node_text(
          aes(label = name),
          repel = TRUE,
          size = 3,
          max.overlaps = 200, 
          min.segment.length = 0.2,
          box.padding = 0.6,
          point.padding = 0.6,
          segment.size = 0.4,
          segment.alpha = 0.5,
          segment.color = "grey40",
          force = 1.5,
          force_pull = 0.8
        ) +
        # 边宽度梯度（CS值）
        scale_edge_width_continuous(
          name = "Causal Strength (CS)",
          range = c(0.5, 3),  # 线条宽度范围，CS值越大线条越粗
          breaks = c(0.1, 0.5, 1.0, 1.5, 2.0),
          # breaks = seq(min(connections$cs, na.rm = TRUE), 
          #             max(connections$cs, na.rm = TRUE), 
          #            length.out = 4),
          guide = guide_legend(title.position = "top", title.hjust = 0.5)
        ) +
        scale_edge_alpha_continuous(
          name = "Causal Strength (CS)",
          range = c(0.4, 1),
          guide = "none"  # 透明度与宽度一致，不需要单独的图例
        ) +
        # 节点颜色
        scale_color_manual(
          values = c(
            "outcome" = "#E74C3C",  # 红色
            "direct" = "#3498DB",   # 蓝色
            "mediator" = "#F39C12", # 橙色
            "source" = "#27AE60"    # 绿色
          ),
          name = "Node Type",
          guide = guide_legend(title.position = "top", title.hjust = 0.5)
        ) +
        scale_size_continuous(
          name = "Node Importance",
          #range = c(3, 10),
          #breaks = c(3, 6, 9),
          #labels = c("Low", "Medium", "High"),
          range = c(3, 10),   # 这里我们可以保持固定的点的大小范围，或者根据实际数据调整
          breaks = pretty(node_importance$size, n = 5), 
          guide = "none"
          #guide = guide_legend(title.position = "top", title.hjust = 0.5)
        ) +
        theme_void() +
        labs(
          title = "Influence Path Hierarchy Structure"#,
          #  subtitle = "Line width indicates causal strength (thicker = stronger); Node size indicates importance"
        ) +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 9),
          legend.position = "bottom",
          legend.box = "vertical",
          legend.direction = "horizontal",
          legend.spacing = unit(0.1, "cm"),
          legend.margin = margin(0, 0, 5, 0),
          plot.margin = margin(20, 20, 20, 20)
        ) +
        coord_cartesian(clip = "off")
    }
  }
  
  # 7. Arrange all plots
  plot_list <- list(network_plot = network_plot)
  
  if(nrow(direct_edges) > 0) {
    plot_list$direct_effects_plot <- direct_effects_plot
  }
  
  if(nrow(top_indirect_paths) > 0) {
    plot_list$indirect_pathways_plot <- indirect_pathways_plot
    plot_list$path_length_distribution <- path_length_dist
  }
  
  if(!is.null(mediator_analysis_result)) {
    plot_list$mediator_importance_plot <- mediator_analysis_result$plot
    plot_list$mediator_stats <- mediator_analysis_result$mediator_stats
  }
  
  if(!is.null(hierarchy_plot)) {
    plot_list$hierarchy_plot <- hierarchy_plot
  }
  
  if(length(plot_list) >= 4) {
    basic_plots <- plot_list[names(plot_list) %in% 
                               c("network_plot", "direct_effects_plot","hierarchy_plot",
                                 "indirect_pathways_plot", "mediator_importance_plot")]
    
    if(length(basic_plots) >= 3) {
      # combined_plot <- basic_plots$hierarchy_plot + 
      #   (basic_plots$direct_effects_plot / plot_list$mediator_importance_plot) +
      #   plot_layout(ncol = 2, widths = c(2, 1)) +
      #   plot_annotation(
      #     title = paste("Comprehensive Pathway Analysis for", outcome_var),
      #     theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
      #   )
      combined_plot <- plot_grid(basic_plots$hierarchy_plot,
                plot_grid(basic_plots$direct_effects_plot,
                          plot_list$mediator_importance_plot,ncol = 2,labels = c("B","C")),
                rel_heights = c(0.7,0.3),
                ncol = 1,labels = c("A",""))
      plot_list$combined_plot <- combined_plot
    }
  }
  
  # 8. Print summary
  cat("=== PATHWAY VISUALIZATION SUMMARY ===\n")
  cat("Outcome variable:", outcome_var, "\n")
  cat("Total nodes in visualization:", ifelse(is.null(summary_stats$total_nodes), "N/A", summary_stats$total_nodes), "\n")
  cat("Direct regulators:", ifelse(is.null(summary_stats$direct_regulators), "0", summary_stats$direct_regulators), "\n")
  cat("Indirect paths found:", ifelse(is.null(summary_stats$indirect_paths_found), "0", summary_stats$indirect_paths_found), "\n")
  
  if(nrow(direct_edges) > 0) {
    cat("Top direct regulators:", paste(direct_edges$regulator, collapse = ", "), "\n")
  }
  
  if(nrow(top_indirect_paths) > 0) {
    cat("\nTop indirect pathways:\n")
    for(i in 1:min(5, nrow(top_indirect_paths))) {
      cat(i, ")", top_indirect_paths$full_path[i], 
          "(strength:", round(top_indirect_paths$path_strength[i], 3), ")\n")
    }
  }
  
  if(!is.null(mediator_analysis_result) && nrow(mediator_analysis_result$mediator_stats) > 0) {
    cat("\nTop mediators by total path strength:\n")
    top_mediators <- head(mediator_analysis_result$mediator_stats, 5)
    for(i in 1:nrow(top_mediators)) {
      cat(i, ")", top_mediators$mediators[i], 
          "(total strength:", round(top_mediators$total_strength[i], 3), 
          ", paths:", top_mediators$n_paths[i], ")\n")
    }
  }
  
  return(plot_list)
}


find_all_indirect_paths <- function(weightNet, source_nodes = NULL, target_nodes = NULL, 
                                    max_path_length = 4, min_strength = 0.1, 
                                    p_threshold = 0.05, cs_threshold = 0.1) {
  library(igraph)
  library(dplyr)
  
  # Keep only significant edges
  sig_edges <- weightNet %>% 
    filter(p_adj < p_threshold & cs >= cs_threshold) %>%
    dplyr::select(regulator, target, weight = cs)
  
  # Check if there are enough significant edges
  if(nrow(sig_edges) == 0) {
    warning("No significant edges found with the given thresholds.")
    return(data.frame())
  }
  
  # Create directed graph
  g <- graph_from_data_frame(sig_edges, directed = TRUE)
  E(g)$weight <- sig_edges$weight
  
  # If source and target nodes are not specified, use all nodes
  if(is.null(source_nodes)) source_nodes <- V(g)$name
  if(is.null(target_nodes)) target_nodes <- V(g)$name
  
  all_paths <- list()
  path_count <- 0
  
  # Find paths for each source-target node pair
  for(source in source_nodes) {
    # Check if source node exists in graph
    if(!source %in% V(g)$name) next
    
    for(target in target_nodes) {
      if(source == target) next
      
      # Check if target node exists in graph
      if(!target %in% V(g)$name) next
      
      # Find all simple paths
      paths <- tryCatch({
        all_simple_paths(g, from = source, to = target, 
                         mode = "out", cutoff = max_path_length)
      }, error = function(e) {
        # If no paths found, return empty list
        list()
      })
      
      for(path in paths) {
        path_length <- length(path) - 1  # Path length = number of edges
        
        # Only keep paths with length >= 2 (indirect paths)
        if(path_length >= 2) {
          path_names <- names(path)
          
          # Calculate path strength (using geometric mean)
          edge_weights <- sapply(1:path_length, function(i) {
            edge_weight <- E(g)[path_names[i] %->% path_names[i+1]]$weight
            if(length(edge_weight) == 0) 1 else edge_weight
          })
          
          path_strength <- exp(mean(log(edge_weights)))
          
          # Only keep paths with sufficient strength
          if(path_strength >= min_strength) {
            path_count <- path_count + 1
            
            all_paths[[path_count]] <- data.frame(
              source = path_names[1],
              target = path_names[length(path_names)],
              mediators = I(list(path_names[2:(length(path_names)-1)])),
              path_length = path_length,
              path_strength = path_strength,
              full_path = paste(path_names, collapse = " -> "),
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
  }
  
  # Convert to data frame and sort
  if(length(all_paths) > 0) {
    paths_df <- do.call(rbind, all_paths)
    paths_df <- paths_df[order(-paths_df$path_strength), ]
    rownames(paths_df) <- NULL
    return(paths_df)
  } else {
    return(data.frame())
  }
}


create_detailed_pathway_tables <- function(analysis_results,
                                           top_direct = 10,
                                           top_indirect = 15) {
  library(knitr)
  library(kableExtra)
  library(DT)
  library(dplyr)
  
  if(is.null(analysis_results)) {
    warning("No analysis results provided.")
    return(NULL)
  }
  
  # 从analyze_pathways结果中提取数据
  direct_edges <- arrange(analysis_results$direct_edges,desc(cs))[1:top_direct,]
  indirect_paths <- arrange(analysis_results$indirect_paths,desc(path_strength))[1:top_indirect,]
  
  summary_stats <- analysis_results$summary
  
  tables_list <- list()
  
  # 1. 创建直接效应表格
  if(!is.null(direct_edges) && nrow(direct_edges) > 0) {
    direct_table <- direct_edges %>%
      mutate(
        cs = round(cs, 3),
        p_adj = format(p_adj, scientific = TRUE, digits = 2),
        significance = ifelse(p_adj < 0.001, "***", 
                              ifelse(p_adj < 0.01, "**",
                                     ifelse(p_adj < 0.05, "*", "")))
      ) %>%
      dplyr::select(regulator, target, cs, p_adj, significance) %>%
      dplyr::rename(
        Regulator = regulator,
        Target = target,
        `Causal Strength` = cs,
        `Adjusted P-value` = p_adj,
        Significance = significance
      )
    
    # 打印直接效应表格
    cat("\n=== DIRECT EFFECTS ===\n")
    print(kable(direct_table, format = "simple", row.names = FALSE))
    
    # 创建交互式表格
    direct_dt <- datatable(
      direct_table,
      options = list(
        pageLength = 20,
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel')
      ),
      caption = "Direct Causal Effects on Outcome Variable",
      rownames = FALSE
    ) %>%
      formatStyle(
        'Causal Strength',
        background = styleColorBar(range(direct_table$`Causal Strength`, na.rm = TRUE), 'lightblue'),
        backgroundSize = '98% 88%',
        backgroundRepeat = 'no-repeat',
        backgroundPosition = 'center'
      )
    
    tables_list$direct_effects <- direct_table
    tables_list$direct_dt <- direct_dt
  } else {
    direct_table <- data.frame(Message = "No significant direct effects found")
    cat("\n=== DIRECT EFFECTS ===\n")
    cat("No significant direct effects found\n")
    tables_list$direct_effects <- direct_table
    tables_list$direct_dt <- NULL
  }
  
  # 2. 创建间接通路表格
  if(!is.null(indirect_paths) && nrow(indirect_paths) > 0) {
    # 处理中介节点列表，将其转换为字符串
    indirect_table <- indirect_paths
    
    # 如果mediators是列表，转换为逗号分隔的字符串
    if("mediators" %in% names(indirect_table) && is.list(indirect_table$mediators)) {
      indirect_table$mediators_str <- sapply(indirect_table$mediators, function(x) {
        if(length(x) > 0) paste(x, collapse = ", ") else "None"
      })
    } else {
      indirect_table$mediators_str <- "Not available"
    }
    
    indirect_table <- indirect_table %>%
      mutate(
        path_strength = round(path_strength, 3),
        mediators_count = sapply(mediators, length)
      ) %>%
      dplyr::select(source, target, path_length, mediators_count, path_strength, full_path, mediators_str) %>%
      dplyr::rename(
        Source = source,
        Target = target,
        `Path Length` = path_length,
        `Mediators Count` = mediators_count,
        `Path Strength` = path_strength,
        `Full Path` = full_path,
        `Mediators` = mediators_str
      )
    
    # 打印间接通路表格
    cat("\n=== TOP INDIRECT PATHWAYS ===\n")
    print(kable(indirect_table %>% dplyr::select(-Mediators), format = "simple", row.names = FALSE))
    
    # 创建交互式表格（包含中介节点信息）
    indirect_dt <- datatable(
      indirect_table,
      options = list(
        pageLength = 20,
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel'),
        columnDefs = list(
          list(targets = 6, width = '200px')  # 为Mediators列设置固定宽度
        )
      ),
      caption = "Top Indirect Causal Pathways to Outcome Variable",
      rownames = FALSE
    ) %>%
      formatStyle(
        'Path Strength',
        background = styleColorBar(range(indirect_table$`Path Strength`, na.rm = TRUE), 'lightgreen'),
        backgroundSize = '98% 88%',
        backgroundRepeat = 'no-repeat',
        backgroundPosition = 'center'
      )
    
    tables_list$indirect_pathways <- indirect_table
    tables_list$indirect_dt <- indirect_dt
  } else {
    indirect_table <- data.frame(Message = "No significant indirect pathways found")
    cat("\n=== INDIRECT PATHWAYS ===\n")
    cat("No significant indirect pathways found\n")
    tables_list$indirect_pathways <- indirect_table
    tables_list$indirect_dt <- NULL
  }
  
  # 3. 创建节点统计表格
  if(!is.null(analysis_results$nodes_df) && nrow(analysis_results$nodes_df) > 0) {
    node_stats <- analysis_results$nodes_df %>%
      mutate(
        degree = as.numeric(degree),
        betweenness = round(as.numeric(betweenness), 2)
      ) %>%
      dplyr::select(name, category, degree, betweenness) %>%
      dplyr::rename(
        Node = name,
        Category = category,
        `In-Degree` = degree,
        `Betweenness Centrality` = betweenness
      ) %>%
      arrange(desc(`Betweenness Centrality`))
    
    cat("\n=== NODE STATISTICS (Top 10 by Betweenness) ===\n")
    print(kable(head(node_stats, 10), format = "simple", row.names = FALSE))
    
    node_dt <- datatable(
      node_stats,
      options = list(
        pageLength = 20,
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel')
      ),
      caption = "Node Statistics in Causal Network",
      rownames = FALSE
    ) %>%
      formatStyle(
        'Betweenness Centrality',
        background = styleColorBar(range(node_stats$`Betweenness Centrality`, na.rm = TRUE), 'lightcoral'),
        backgroundSize = '98% 88%',
        backgroundRepeat = 'no-repeat',
        backgroundPosition = 'center'
      )
    
    tables_list$node_statistics <- node_stats
    tables_list$node_dt <- node_dt
  }
  
  # 4. 创建网络统计表格
  if(!is.null(summary_stats)) {
    network_stats <- data.frame(
      Metric = c(
        "Outcome Variable", 
        "Total Nodes", 
        "Direct Regulators", 
        "Indirect Pathways Found",
        "Average Direct Causal Strength",
        "Average Indirect Path Strength"
      ),
      Value = c(
        summary_stats$outcome,
        summary_stats$total_nodes,
        summary_stats$direct_regulators,
        summary_stats$indirect_paths_found,
        ifelse(!is.null(direct_edges) && nrow(direct_edges) > 0, 
               round(mean(direct_edges$cs, na.rm = TRUE), 3), "N/A"),
        ifelse(!is.null(indirect_paths) && nrow(indirect_paths) > 0, 
               round(mean(indirect_paths$path_strength, na.rm = TRUE), 3), "N/A")
      )
    )
    
    cat("\n=== NETWORK STATISTICS ===\n")
    print(kable(network_stats, format = "simple", row.names = FALSE))
    
    tables_list$network_statistics <- network_stats
  }
  
  # 5. 创建中介节点重要性表格（如果存在）
  if(!is.null(indirect_paths) && nrow(indirect_paths) > 0 && 
     "mediators" %in% names(indirect_paths)) {
    
    # 计算中介节点重要性
    mediator_importance <- indirect_paths %>%
      tidyr::unnest(mediators) %>%
      group_by(mediators) %>%
      summarise(
        n_paths = n(),
        avg_strength = mean(path_strength, na.rm = TRUE),
        total_strength = sum(path_strength, na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      mutate(
        avg_strength = round(avg_strength, 3),
        total_strength = round(total_strength, 3)
      ) %>%
      arrange(desc(total_strength)) %>%
      dplyr::rename(
        Mediator = mediators,
        `Number of Paths` = n_paths,
        `Average Path Strength` = avg_strength,
        `Total Path Strength` = total_strength
      )
    
    if(nrow(mediator_importance) > 0) {
      cat("\n=== MEDIATOR IMPORTANCE ===\n")
      print(kable(head(mediator_importance, 10), format = "simple", row.names = FALSE))
      
      mediator_dt <- datatable(
        mediator_importance,
        options = list(
          pageLength = 20,
          dom = 'Bfrtip',
          buttons = c('copy', 'csv', 'excel')
        ),
        caption = "Mediator Node Importance in Indirect Pathways",
        rownames = FALSE
      ) %>%
        formatStyle(
          'Total Path Strength',
          background = styleColorBar(range(mediator_importance$`Total Path Strength`, na.rm = TRUE), 'lightgoldenrodyellow'),
          backgroundSize = '98% 88%',
          backgroundRepeat = 'no-repeat',
          backgroundPosition = 'center'
        )
      
      tables_list$mediator_importance <- mediator_importance
      tables_list$mediator_dt <- mediator_dt
    }
  }
  return(tables_list)
}


save_pathway_analysis <- function(visualization_result, analysis_results = NULL, filename_prefix = "pathway_analysis") {
  if(is.null(visualization_result)) {
    warning("No visualization results provided.")
    return(NULL)
  }
  
  if(!dir.exists(dirname(filename_prefix))) {
    dir.create(dirname(filename_prefix), recursive = TRUE)
  }
  
  if(!is.null(visualization_result$combined_plot)) {
    tryCatch({
      ggsave(
        paste0(filename_prefix, "_combined_plot.png"),
        visualization_result$combined_plot,
        width = 14, height = 8, dpi = 300
      )
      cat("Saved:", paste0(filename_prefix, "_combined_plot.png"), "\n")
    }, error = function(e) {
      warning("Failed to save combined plot: ", e$message)
    })
  }
  
  # 保存网络图
  if(!is.null(visualization_result$network_plot)) {
    tryCatch({
      ggsave(
        paste0(filename_prefix, "_network.png"),
        visualization_result$network_plot,
        width = 12, height = 10, dpi = 300
      )
      cat("Saved:", paste0(filename_prefix, "_network.png"), "\n")
    }, error = function(e) {
      warning("Failed to save network plot: ", e$message)
    })
  }
  
  # 保存直接效应图
  if(!is.null(visualization_result$direct_effects_plot)) {
    tryCatch({
      ggsave(
        paste0(filename_prefix, "_direct_effects.png"),
        visualization_result$direct_effects_plot,
        width = 8, height = 6, dpi = 300
      )
      cat("Saved:", paste0(filename_prefix, "_direct_effects.png"), "\n")
    }, error = function(e) {
      warning("Failed to save direct effects plot: ", e$message)
    })
  }
  
  # 保存间接通路图
  if(!is.null(visualization_result$indirect_pathways_plot)) {
    tryCatch({
      ggsave(
        paste0(filename_prefix, "_indirect_pathways.png"),
        visualization_result$indirect_pathways_plot,
        width = 8, height = 6, dpi = 300
      )
      cat("Saved:", paste0(filename_prefix, "_indirect_pathways.png"), "\n")
    }, error = function(e) {
      warning("Failed to save indirect pathways plot: ", e$message)
    })
  }
  
  # 保存路径长度分布图
  if(!is.null(visualization_result$path_length_distribution)) {
    tryCatch({
      ggsave(
        paste0(filename_prefix, "_path_length_distribution.png"),
        visualization_result$path_length_distribution,
        width = 8, height = 6, dpi = 300
      )
      cat("Saved:", paste0(filename_prefix, "_path_length_distribution.png"), "\n")
    }, error = function(e) {
      warning("Failed to save path length distribution plot: ", e$message)
    })
  }
  
  # 保存中介重要性图
  if(!is.null(visualization_result$mediator_importance_plot)) {
    tryCatch({
      ggsave(
        paste0(filename_prefix, "_mediator_importance.png"),
        visualization_result$mediator_importance_plot,
        width = 10, height = 8, dpi = 300
      )
      cat("Saved:", paste0(filename_prefix, "_mediator_importance.png"), "\n")
    }, error = function(e) {
      warning("Failed to save mediator importance plot: ", e$message)
    })
  }
  
  # 保存层次结构图
  if(!is.null(visualization_result$hierarchy_plot)) {
    tryCatch({
      ggsave(
        paste0(filename_prefix, "_hierarchy.png"),
        visualization_result$hierarchy_plot,
        width = 12, height = 10, dpi = 300
      )
      cat("Saved:", paste0(filename_prefix, "_hierarchy.png"), "\n")
    }, error = function(e) {
      warning("Failed to save hierarchy plot: ", e$message)
    })
  }
  
  # 保存PDF版本（高质量）
  if(!is.null(visualization_result$combined_plot)) {
    tryCatch({
      ggsave(
        paste0(filename_prefix, "_combined_plot.pdf"),
        visualization_result$combined_plot,
        width = 16, height = 10
      )
      cat("Saved:", paste0(filename_prefix, "_combined_plot.pdf"), "\n")
    }, error = function(e) {
      warning("Failed to save combined plot PDF: ", e$message)
    })
  }
  
  # 保存数据表格
  tables <- NULL
  if(!is.null(analysis_results)) {
    tables <- create_detailed_pathway_tables(analysis_results)
  } else {
    # 如果没有提供analysis_results，尝试从visualization_result中获取
    if(!is.null(visualization_result$summary)) {
      tables <- create_detailed_pathway_tables(visualization_result)
    }
  }
  
  if(!is.null(tables)) {
    # 保存直接效应表格
    if(!is.null(tables$direct_effects) && 
       nrow(tables$direct_effects) > 0 && 
       !"Message" %in% names(tables$direct_effects)) {
      tryCatch({
        write.csv(tables$direct_effects, 
                  paste0(filename_prefix, "_direct_effects.csv"), 
                  row.names = FALSE)
        cat("Saved:", paste0(filename_prefix, "_direct_effects.csv"), "\n")
      }, error = function(e) {
        warning("Failed to save direct effects table: ", e$message)
      })
    }
    
    # 保存间接通路表格
    if(!is.null(tables$indirect_pathways) && 
       nrow(tables$indirect_pathways) > 0 && 
       !"Message" %in% names(tables$indirect_pathways)) {
      tryCatch({
        write.csv(tables$indirect_pathways, 
                  paste0(filename_prefix, "_indirect_pathways.csv"), 
                  row.names = FALSE)
        cat("Saved:", paste0(filename_prefix, "_indirect_pathways.csv"), "\n")
      }, error = function(e) {
        warning("Failed to save indirect pathways table: ", e$message)
      })
    }
    
    # 保存节点统计表格
    if(!is.null(tables$node_statistics) && nrow(tables$node_statistics) > 0) {
      tryCatch({
        write.csv(tables$node_statistics, 
                  paste0(filename_prefix, "_node_statistics.csv"), 
                  row.names = FALSE)
        cat("Saved:", paste0(filename_prefix, "_node_statistics.csv"), "\n")
      }, error = function(e) {
        warning("Failed to save node statistics table: ", e$message)
      })
    }
    
    # 保存中介重要性表格
    if(!is.null(tables$mediator_importance) && nrow(tables$mediator_importance) > 0) {
      tryCatch({
        write.csv(tables$mediator_importance, 
                  paste0(filename_prefix, "_mediator_importance.csv"), 
                  row.names = FALSE)
        cat("Saved:", paste0(filename_prefix, "_mediator_importance.csv"), "\n")
      }, error = function(e) {
        warning("Failed to save mediator importance table: ", e$message)
      })
    }
    
    # 保存网络统计表格
    if(!is.null(tables$network_statistics) && nrow(tables$network_statistics) > 0) {
      tryCatch({
        write.csv(tables$network_statistics, 
                  paste0(filename_prefix, "_network_statistics.csv"), 
                  row.names = FALSE)
        cat("Saved:", paste0(filename_prefix, "_network_statistics.csv"), "\n")
      }, error = function(e) {
        warning("Failed to save network statistics table: ", e$message)
      })
    }
  }
  
  # 保存R数据文件（包含完整结果）
  tryCatch({
    saveRDS(visualization_result, file = paste0(filename_prefix, "_visualization_results.rds"))
    cat("Saved:", paste0(filename_prefix, "_visualization_results.rds"), "\n")
    
    if(!is.null(analysis_results)) {
      saveRDS(analysis_results, file = paste0(filename_prefix, "_analysis_results.rds"))
      cat("Saved:", paste0(filename_prefix, "_analysis_results.rds"), "\n")
    }
  }, error = function(e) {
    warning("Failed to save R data files: ", e$message)
  })
  
  # 创建保存摘要
  saved_files <- list.files(path = dirname(filename_prefix), 
                            pattern = paste0(basename(filename_prefix), ".*"), 
                            full.names = TRUE)
  
  cat("\n=== SAVE SUMMARY ===\n")
  cat("All pathway analysis results saved with prefix:", filename_prefix, "\n")
  cat("Total files saved:", length(saved_files), "\n")
  cat("Files:\n")
  for(file in saved_files) {
    cat(" -", basename(file), "\n")
  }
  
  return(list(
    tables = tables,
    saved_files = saved_files
  ))
}



get_network_statistics <- function(weightNet, p_adj_threshold = 0.05, cs_threshold = 0.1) {
  
  sig_edges <- weightNet %>%
    filter(p_adj < p_adj_threshold & abs(cs) > cs_threshold)
  
  if (nrow(sig_edges) == 0) {
    return(list(edges = 0, nodes = 0, message = "No significant network"))
  }
  
  g <- graph_from_data_frame(sig_edges[, c("regulator", "target")], directed = TRUE)
  
  stats <- list(
    n_edges = ecount(g),
    n_nodes = vcount(g),
    density = graph.density(g),
    transitivity = transitivity(g),
    diameter = diameter(g),
    avg_path_length = average.path.length(g),
    positive_edges = sum(sig_edges$t1 > 0),
    negative_edges = sum(sig_edges$t1 < 0),
    avg_causal_strength = mean(abs(sig_edges$cs)),
    hub_nodes = names(sort(degree(g, mode = "out"), decreasing = TRUE)[1:5]))
  
  return(stats)
}



plot_driver_ranking <- function(weightNet, top_n = 20) {
  # 计算每个调节因子的总因果强度
  driver_strength <- aggregate(cs ~ regulator, data = weightNet, 
                               function(x) sum(abs(x)))
  
  # 计算每个调节因子的显著连接数
  driver_connections <- aggregate(p_adj ~ regulator, data = weightNet,
                                  function(x) sum(x < 0.05))
  
  # 合并数据
  driver_stats <- merge(driver_strength, driver_connections, by = "regulator")
  colnames(driver_stats) <- c("regulator", "total_cs", "sig_connections")
  
  # 综合评分（可根据需要调整权重）
  driver_stats$score <- driver_stats$total_cs * log(driver_stats$sig_connections + 1)
  
  # 取top N
  top_drivers <- head(driver_stats[order(-driver_stats$score), ], top_n)
  
  # 绘制条形图
  ggplot(top_drivers, aes(x = reorder(regulator, score), y = score)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
    coord_flip() +
    labs(title = "Top Key Driver Factors",
         x = "Regulator", 
         y = "Driver Score (Causal Strength × Connections)") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 10))
}



visualize_whole_network <- function(weightNet, 
                                       p_adj_threshold = 0.05, 
                                       cs_threshold = 0.1,
                                       node_size_range = c(2, 8), 
                                       edge_width_range = c(0.1, 1),
                                       layout_algorithm = "stress",
                                       color_palette = "Set1") {
  
  # 过滤显著边
  sig_edges <- weightNet %>%
    filter(p_adj < p_adj_threshold & abs(cs) > cs_threshold) %>%
    mutate(
      direction = ifelse(t1 > 0, "Positive", "Negative"),
      edge_alpha = scales::rescale(abs(cs), to = c(0.1, 0.8)),  # 透明度基于因果强度
      edge_width = scales::rescale(abs(cs), to = edge_width_range)
    )
  
  if (nrow(sig_edges) == 0) {
    warning("No significant edges found with the current thresholds.")
    return(NULL)
  }
  
  # 创建图对象
  g <- graph_from_data_frame(
    d = sig_edges[, c("regulator", "target", "cs", "t1", "direction", "edge_alpha", "edge_width")],
    directed = TRUE
  )
  
  # 计算节点中心性指标
  V(g)$out_degree <- degree(g, mode = "out")
  V(g)$in_degree <- degree(g, mode = "in")
  V(g)$betweenness <- betweenness(g, directed = TRUE)
  
  # 设置节点属性
  V(g)$size <- scales::rescale(V(g)$out_degree, to = node_size_range)
  V(g)$label <- V(g)$name
  V(g)$community <- as.factor(cluster_walktrap(g, steps = 4)$membership)
  
  # 设置边的属性
  E(g)$color <- ifelse(E(g)$t1 > 0, "#E41A1C", "#377EB8")  # 红色正相关，蓝色负相关
  E(g)$alpha <- E(g)$edge_alpha
  E(g)$width <- E(g)$edge_width
  
  # 创建ggraph绘图
  p <- ggraph(g, layout = layout_algorithm) +
    
    # 绘制边
    geom_edge_fan(
      aes(edge_color = color, edge_width = width, edge_alpha = alpha),
      arrow = arrow(length = unit(2, "mm"), type = "closed"),
      end_cap = circle(3, "mm"),
      start_cap = circle(3, "mm")
    ) +
    
    # 绘制节点
    geom_node_point(
      aes(size = size, fill = community),
      shape = 21,
      color = "white",
      stroke = 0.5
    ) +
    
    # 节点标签
    geom_node_text(
      aes(label = label),
      size = 3,
      repel = TRUE,
      box.padding = 0.5,
      point.padding = 0.2,
      max.overlaps = 20
    ) +
    
    # 标度和主题
    scale_edge_color_identity() +
    scale_edge_alpha_identity() +
    scale_edge_width_identity() +
    scale_size_identity() +
    scale_fill_brewer(palette = color_palette) +
    
    # 主题设置
    theme_void() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      legend.position = "bottom",
      plot.margin = unit(c(1, 1, 1, 1), "cm")
    ) +
    
    # 标题和图例
    labs(
      title = "Causal Protein Interaction Network",
      subtitle = paste("KOCMI Analysis |", nrow(sig_edges), "significant edges"),
      fill = "Functional Modules",
      size = "Regulatory Activity (Out-degree)",
      edge_width = "Causal Strength",
      edge_color = "Causal Direction"
    )
  
  return(p)
}


visualize_outcome_network <- function(weightNet, outcome_var, top_n = 20, show_legend = TRUE) {
  library(igraph)
  library(ggraph)
  library(ggplot2)
  library(tidygraph)
  library(dplyr)
  
  # Extract edges related to the outcome variable
  outcome_related <- weightNet[
    weightNet$target == outcome_var | weightNet$regulator == outcome_var, 
  ]
  
  # Select top_n strongest relationships
  outcome_related <- outcome_related[order(-outcome_related$cs), ]
  outcome_related <- head(outcome_related, top_n)
  
  # Create graph object
  g <- graph_from_data_frame(outcome_related[, c("regulator", "target")], 
                             directed = TRUE)
  
  # Add edge attributes
  E(g)$causal_strength <- outcome_related$cs
  E(g)$p_adjusted <- outcome_related$p_adj
  E(g)$significance <- ifelse(outcome_related$p_adj < 0.01, "p < 0.01",
                              ifelse(outcome_related$p_adj < 0.05, "p < 0.05", "p ≥ 0.05"))
  
  # Add node attributes
  V(g)$node_type <- ifelse(V(g)$name == outcome_var, "outcome", "regulator")
  V(g)$degree <- degree(g)
  V(g)$betweenness <- betweenness(g)
  
  # Convert to tidygraph for better integration with ggraph
  tg <- as_tbl_graph(g)
  
  # Create the visualization
  network_plot <- ggraph(tg, layout = "stress") +
    # Edges with varying width and color based on causal strength and significance
    geom_edge_fan(
      aes(edge_width = causal_strength, 
          edge_alpha = causal_strength,
          color = significance),
      arrow = arrow(length = unit(2, "mm"), type = "closed"),
      end_cap = circle(3, "mm"),
      show.legend = show_legend  # 根据参数控制边图例显示
    ) +
    # Nodes with varying size and color
    geom_node_point(
      aes(color = node_type,
          fill = node_type),
      shape = 21, size = 3,
      stroke = 1.5,
      show.legend = show_legend  # 根据参数控制节点图例显示
    ) +
    # Node labels
    geom_node_text(
      aes(label = name,
          color = node_type),
      size = 3,
      repel = TRUE,
      fontface = "bold",
      show.legend = FALSE  # 文本标签永远不显示图例
    ) +
    # Scales and colors
    scale_edge_width_continuous(
      name = "Causal Strength: ",
      range = c(0.5, 2),
      guide = if(show_legend) guide_legend(order = 1) else "none"
    ) +
    scale_edge_alpha_continuous(
      name = "Causal Strength: ",
      guide = if(show_legend) guide_legend(order = 1) else "none"
    ) +
    scale_edge_color_manual(
      name = "Significance: ",
      values = c("p < 0.01" = "#E31A1C", 
                 "p < 0.05" = "#FDBF6F",
                 "p ≥ 0.05" = "#B2B2B2"),
      guide = if(show_legend) guide_legend(order = 2) else "none"
    ) +
    scale_color_manual(
      name = "Node Type: ",
      values = c("outcome" = "#E31A1C", 
                 "regulator" = "#1F78B4"),
      labels = c("outcome" = "Outcome Variable", 
                 "regulator" = "Regulator Protein"),
      guide = if(show_legend) guide_legend(order = 3) else "none"
    ) +
    scale_fill_manual(
      name = "Node Type: ",
      values = c("outcome" = "#FB9A99", 
                 "regulator" = "#A6CEE3"),
      labels = c("outcome" = "Outcome Variable", 
                 "regulator" = "Regulator Protein"),
      guide = if(show_legend) guide_legend(order = 3) else "none"
    ) +
    # 主题设置 - 根据参数控制图例位置
    theme_void() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      plot.caption = element_text(size = 9, hjust = 0.5, color = "gray50"),
      legend.position = if(show_legend) "bottom" else "none",
      #plot.margin = unit(c(1, 1, 1, 1), "cm"),
      legend.box = if(show_legend) "vertical" else NULL,
      legend.text = if(show_legend) element_text(size = 8) else element_blank(),
      legend.title = if(show_legend) element_text(size = 9, face = "bold") else element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  # 如果显示图例，添加图例覆盖设置
  if(show_legend) {
    network_plot <- network_plot +
      guides(
        edge_width = guide_legend(override.aes = list(edge_color = "black")),
        edge_alpha = guide_legend(override.aes = list(edge_color = "black")),
        size = guide_legend(override.aes = list(color = "black"))
      )
  }
  return(network_plot)
}

# Additional function to create a companion centrality bar plot
create_centrality_companion_plot <- function(weightNet, outcome_var, top_n = 20) {
  library(ggplot2)
  library(dplyr)
  
  # Extract and prepare data
  outcome_related <- weightNet[
    weightNet$target == outcome_var | weightNet$regulator == outcome_var, 
  ]
  outcome_related <- outcome_related[order(-outcome_related$cs), ]
  outcome_related <- head(outcome_related, top_n)
  
  # Create bar plot for causal strengths
  bar_plot <- outcome_related %>%
    mutate(
      relationship = paste(regulator, "→", target),
      relationship = factor(relationship, levels = relationship[order(cs)])
    ) %>%
    ggplot(aes(x = cs, y = reorder(relationship, cs), fill = cs)) +
    geom_col(width = 0.7) +
    #geom_text(aes(label = formatC(cs, format = "f", digits = 2)), hjust = -0.1, size = 3) +
    scale_fill_gradient(
      low = "#DEEBF7", 
      high = "#2171B5",
      name = "Causal\nStrength"
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(
      x = "Causal Strength (cs)",
      y = "Protein Relationships"#,
      #title = "Causal Strength of Protein-Outcome Relationships",
      #subtitle = paste("Top", top_n, "relationships for", outcome_var)
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.text.y = element_text(size = 9),
      # plot.margin = unit(c(0.5, -0.5, 0.5, 0.5), "cm"),
      panel.grid.major.y = element_blank(),
      legend.position = "none"
    )
  
  return(bar_plot)
}

result_summary <- function(weightNet) {
  cat("=== KOCMI 网络推断结果统计 ===\n")
  cat("总边数:", nrow(weightNet), "\n")
  cat("显著边数 (p_adj < 0.05):", sum(weightNet$p_adj < 0.05), "\n")
  cat("强因果关系 (CS ≥ 0.3):", sum(weightNet$cs >= 0.3,na.rm = T), "\n")
  cat("中等因果关系 (0.2 ≤ CS < 0.3):", sum(weightNet$cs >= 0.2 & weightNet$cs < 0.3,na.rm = T), "\n")
  cat("弱因果关系 (0.1 ≤ CS < 0.2):", sum(weightNet$cs >= 0.1 & weightNet$cs < 0.2,na.rm = T), "\n")
  
  # 显著性水平分布
  sig_levels <- c(
    "p_adj < 0.001" = sum(weightNet$p_adj < 0.001),
    "0.001 ≤ p_adj < 0.01" = sum(weightNet$p_adj >= 0.001 & weightNet$p_adj < 0.01),
    "0.01 ≤ p_adj < 0.05" = sum(weightNet$p_adj >= 0.01 & weightNet$p_adj < 0.05),
    "0.05 ≤ p_adj < 0.1" = sum(weightNet$p_adj >= 0.05 & weightNet$p_adj < 0.1)
  )
  cat("\n--- 显著性水平分布 ---\n")
  print(sig_levels)
}


main_analysis_vis <- function(weightNet, 
                                 outcome_var, 
                                 top_n = 20,
                                 labels_sub = "",
                                 show_legend = TRUE) {
  
  # Create network visualization
  network_viz <- visualize_outcome_network(weightNet, outcome_var, top_n,show_legend)
  
  # Create companion bar plot
  bar_viz <- create_centrality_companion_plot(weightNet, outcome_var, top_n)
  
  combined_plot <- plot_grid(network_viz,bar_viz,ncol = 2,labels = labels_sub)
  
  return(combined_plot)
}


causal_network_vis_enhanced <- function(weightNet, p_adj_threshold = 0.05,
                                         cs_threshold = 0.1,
                                         layout_type = "stress",  # "stress", "fr", "kk", "dh"
                                         node_color_by = "community",  # "community", "degree", "betweenness"
                                         show_edge_labels = FALSE,
                                         highlight_hubs = TRUE,
                                         max_nodes = 30) {
  
  # 过滤和准备数据
  sig_edges <- weightNet %>%
    filter(p_adj < p_adj_threshold & abs(cs) > cs_threshold)
  
  # 如果节点太多，选择最重要的子网络
  all_nodes <- unique(c(sig_edges$regulator, sig_edges$target))
  if (length(all_nodes) > max_nodes) {
    node_importance <- table(c(sig_edges$regulator, sig_edges$target))
    important_nodes <- names(sort(node_importance, decreasing = TRUE)[1:max_nodes])
    sig_edges <- sig_edges %>%
      filter(regulator %in% important_nodes & target %in% important_nodes)
    warning(paste("Network too large. Showing top", max_nodes, "most connected nodes."))
  }
  
  if (nrow(sig_edges) == 0) {
    stop("No significant edges found with the current thresholds.")
  }
  
  # 创建图对象
  g <- graph_from_data_frame(sig_edges[, c("regulator", "target", "cs", "t1")], 
                             directed = TRUE)
  
  # 计算网络指标
  V(g)$out_degree <- degree(g, mode = "out")
  V(g)$in_degree <- degree(g, mode = "in")
  V(g)$total_degree <- degree(g, mode = "total")
  V(g)$betweenness <- betweenness(g, directed = TRUE)
  V(g)$eigen_centrality <- eigen_centrality(g, directed = TRUE)$vector
  
  # 社区检测
  communities <- cluster_louvain(as.undirected(g))
  V(g)$community <- as.factor(communities$membership)
  
  # 设置节点颜色
  if (node_color_by == "community") {
    V(g)$color_var <- V(g)$community
    color_title <- "Functional Modules"
  } else if (node_color_by == "degree") {
    V(g)$color_var <- cut(V(g)$total_degree, breaks = 5, labels = FALSE)
    color_title <- "Connectivity Level"
  } else if (node_color_by == "betweenness") {
    V(g)$color_var <- cut(V(g)$betweenness, breaks = 5, labels = FALSE)
    color_title <- "Betweenness Centrality"
  }
  
  # 设置节点大小（基于出度，表示调控活性）
  V(g)$size <- scales::rescale(V(g)$out_degree, to = c(4, 12))
  
  # 设置边的属性
  E(g)$direction <- ifelse(E(g)$t1 > 0, "Positive", "Negative")
  E(g)$width <- scales::rescale(abs(E(g)$cs), to = c(0.8, 3))
  E(g)$alpha <- scales::rescale(abs(E(g)$cs), to = c(0.4, 1))
  
  # 创建绘图
  p <- ggraph(g, layout = layout_type) +
    
    # 边图层
    geom_edge_link(
      aes(edge_width = width, edge_alpha = alpha, 
          color = direction),
      arrow = arrow(length = unit(3, "mm"), type = "closed"),
      end_cap = circle(4, "mm"),
      start_cap = circle(4, "mm")
    ) +
    
    # 可选：边标签
    {if (show_edge_labels) {
      geom_edge_label(
        aes(label = sprintf("%.2f", cs)),
        size = 2.5,
        label.padding = unit(0.1, "lines"),
        alpha = 0.7
      )
    }} +
    
    # 节点图层
    geom_node_point(
      aes(size = size, fill = color_var),
      shape = 21,
      color = "white",
      stroke = 0.8
    ) +
    
    # 突出显示枢纽节点
    {if (highlight_hubs) {
      hub_nodes <- V(g)[out_degree > quantile(out_degree, 0.8)]
      geom_node_point(
        data = function(nodes) { nodes %>% filter(name %in% hub_nodes$name) },
        aes(size = size * 1.2),
        shape = 21,
        color = "gold",
        stroke = 1.5,
        fill = NA
      )
    }} +
    
    # 节点标签
    geom_node_text(
      aes(label = name, size = size),
      repel = TRUE,
      box.padding = 0.8,
      point.padding = 0.3,
      max.overlaps = 15,
      fontface = "bold"
    ) +
    
    # 标度设置
    scale_edge_width_continuous(
      range = range(E(g)$width),
      name = "Causal Strength",
      guide = guide_legend(override.aes = list(edge_alpha = 0.8))
    ) +
    scale_edge_alpha_identity() +
    scale_edge_color_manual(
      values = c(Positive = "#D62728", Negative = "#1F77B4"),
      name = "Causal Direction"
    ) +
    scale_size_identity() +
    scale_fill_viridis_d(
      name = color_title,
      option = "plasma",
      begin = 0.2,
      end = 0.8
    ) +
    
    # 主题和布局
    theme_void() +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
      legend.position = "bottom",
      legend.box = "vertical",
      #legend.margin = margin(t = 10),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    ) +
    
    # 标题和图例
    labs(
      title = "Causal Protein Interaction Network",
      subtitle = sprintf("KOCMI Analysis | %d nodes, %d edges | p < %.3f, |cs| > %.2f", 
                         vcount(g), ecount(g), p_adj_threshold, cs_threshold),
      caption = "Node size: Regulatory activity (out-degree)\nEdge width: Causal strength magnitude"
    )
  return(p)
}

