library(shiny)
library(ggplot2)
library(purrr)
library(RANN2)
library(colorspace)
library(grid)
source("data.R")

closest.points = function(emb) {
    wann = WANN(emb %>% as.matrix)
    function(x, y, n, max=Inf) {
        r = wann$query(cbind(x, y), n, 0)
        as.numeric(r$nn.idx)[as.numeric(r$nn.dists) < max]
    }
}

closest.points.cite = closest.points(cite_umap_emb)

colour.scheme = function(mode) {
    switch(mode, "gene"="red", "protein"="blue", "predicted_spatial"="purple")
}

plot.umap = function(title, color, mode) {
    high = ifelse(max(color)==0, "grey", colour.scheme(mode))
    
    if (mode == "predicted_spatial") {
        ggplot(cite_umap_emb) +
            geom_point(aes(x=V2, y=V3, color=color), size=0.3, show.legend=FALSE) +
            geom_point(data = cite_umap_emb[which(color != 0),],
                       aes(x=V2, y=V3, color=color[which(color != 0)]), size=3, show.legend=FALSE) +
            scale_color_gradient(low = "grey", high = high) +
            ggtitle(title) +
            theme_void() + #coord_fixed()
            coord_fixed(ratio=539/759) +
            theme(plot.margin = unit(c(0,0,10,0),'lines'))
    } else {
      ggplot(cite_umap_emb) +
        geom_point(data = cite_umap_emb[which(color == 0),],
                   aes(x=V2, y=V3, color=color[which(color == 0)]), size=0.3, show.legend=FALSE) +
        geom_point(data = cite_umap_emb[which(color != 0),],
                   aes(x=V2, y=V3, color=color[which(color != 0)]), size=0.3, show.legend=FALSE) +
        scale_color_gradient(low = "grey", high = high) +
        ggtitle(title) +
        theme_void() + #coord_fixed()
        coord_fixed(ratio=539/759) +
        theme(plot.margin = unit(c(0,0,10,0),'lines'))
    }
}

plot.pilot.umap = function(title, color, mode) {
    high = ifelse(max(color)==0, "grey", colour.scheme(mode))
    ggplot(tumap, aes(V2,V3, color=color)) + geom_point(show.legend = FALSE, size=0.4) + theme_void() + scale_color_gradient(low = "grey", high = high) + ggtitle(title) + coord_fixed(ratio=1425/2331) +
      theme(plot.margin = unit(c(0,0,14,0),'lines'))
}
plot.spatial = function(spatial, guide.image, value, title, mode) {
    plot_data = as.data.frame(cbind(x=spatial$x, y=spatial$y, expression=value))
    
    if (mode == "protein") {
      ggplot(plot_data, aes(x=spatial$x,y=spatial$y,color=expression)) +
        geom_point(size=0.3) +
        scale_color_gradient(low="white", high=colour.scheme(mode)) +
        guides(color = FALSE) +
        ylim(max(spatial$y), 0) +
        ggtitle(title) +
        theme_void() + coord_fixed() +
        annotation_custom(rasterGrob(guide.image, interpolate=TRUE, x = unit(0.5, "npc"), y = unit(0.5, "npc"),
                                     height=unit(0.92, "npc"), width=unit(0.92, "npc")),
                                     xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)
    }
    else { # duplication
      high = ifelse(max(value)==0, "grey94", colour.scheme(mode))
      sub = plot_data[which(value != 0),]
      ggplot(plot_data) +
          geom_point(aes(x=x,y=y,color=expression), size=0.3, show.legend=FALSE) +
          geom_point(data = sub,
                     aes(x=x,y=y, color=expression), size=0.5, show.legend=FALSE) +
          scale_color_gradient(low="grey94", high=high) +
          guides(color = FALSE) +
          ylim(max(spatial$y), 0) +
          ggtitle(title) +
          theme_void() + coord_fixed() +
          annotation_custom(rasterGrob(guide.image, interpolate=TRUE, x = unit(0.5, "npc"),
                                       y = unit(0.5, "npc"), height=unit(0.92, "npc"), width=unit(0.92, "npc")),
                                       xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)
    }
}

cite.plot = function(mode, input) {
    renderPlot({
        if (mode$mode == "protein") {
            if (input$protein %in% colnames(cite_protein)) {
                plot.umap(paste(input$protein, "protein"), cite_protein[, input$protein], "protein")
            } else if (input$protein %in% colnames(codex_exp_clean_BALBc1)) {
                plot.umap(paste(input$protein, "protein (unavailable)"), 0, "protein")
            }
        } else if (mode$mode == "gene") {
            if (input$gene %in% rownames(cite_gene)) {
                #plot.umap(paste(input$gene, "gene"), cite_gene[input$gene, ], "gene")
              plot.umap(paste(input$gene, "gene"), log(1 + 5000 * cite_gene[input$gene, ]), "gene")
            }
        } else { # predicted_spatial
            selected <- rep_len(0, ncol(cite_gene))
            selected[mode$cite.cells[1]] <- 1
            plot.umap("Selected cell", selected, "predicted_spatial")
        }
    }, width=500, height=500)
}

cite_pilot.plot = function(mode, input) {
    renderPlot({
        if (mode$mode == "gene") {
            if (input$gene %in% rownames(cite_gene)) {
              expression = log2(1+10000*tcounts[input$gene,row.names(tumap)]/colSums(tcounts[,row.names(tumap)]))
              plot.pilot.umap(paste(input$gene, "gene"), expression, "gene")
            }
        } else {
            expression = rep_len(0, nrow(tumap))
            plot.pilot.umap("Protein and spatial data not available", expression, mode$mode)
        }
    }, width=500, height=500)
}

codex.plot = function(mode, input, spatial, guide.image, codex_exp_clean, cite_nn) {
    renderPlot({
      if (mode$mode == "protein") {
        if (input$protein %in% colnames(codex_exp_clean)) {
          plot.spatial(spatial, guide.image, codex_exp_clean[, input$protein], paste(input$protein, "protein"), "protein")
        }
      } else if (mode$mode == "gene") {
        if (input$gene %in% rownames(cite_gene)) {
          codex_gene_expr <- log(1 + 1000 * (cite_nn %*% cite_gene[input$gene, ]))
          plot.spatial(spatial, guide.image, codex_gene_expr %>% as.numeric, paste(input$gene, "predicted gene expression"), "gene")
        }
      } else { # predicted_spatial
        # if no cells within 'max' distance of mouse, nothing is plotted
        if (length(mode$cite.cells) > 1) {
          codex.cells = cite_nn[, mode$cite.cells] %>% rowSums
          plot.spatial(spatial, guide.image, codex.cells, "Predicted CODEX coordinates", "predicted_spatial")
        } else if (length(mode$cite.cells) == 1) {
          codex.cells = cite_nn[, mode$cite.cells]
          plot.spatial(spatial, guide.image, codex.cells, "Predicted CODEX coordinates", "predicted_spatial")
        }
      }
    })
}

server = function(input, output, session) {
    mode <- reactiveValues(mode="protein", click.x=0, click.y=0)
    
    genes = rownames(cite_gene) %>% sort
    proteins = colnames(codex_exp_clean_BALBc1) %>% sort
    
    updateSelectizeInput(session = session, inputId = 'gene', choices = c("", genes), selected = "Mmp12", server = TRUE)
    updateSelectizeInput(session = session, inputId = 'protein', choices = c("", proteins), selected="", server = TRUE)
 
    output$cite = cite.plot(mode, input)
    output$cite_pilot = cite_pilot.plot(mode, input)

    output$codex1 = codex.plot(mode, input, spatial1, guide.image1, codex_exp_clean_BALBc1, cite_nn1)
    output$codex2 = codex.plot(mode, input, spatial2, guide.image2, codex_exp_clean_BALBc2, cite_nn2)
    output$codex3 = codex.plot(mode, input, spatial3, guide.image3, codex_exp_clean_BALBc3, cite_nn3)

    observeEvent(input$gene, ignoreInit = T, {
        if ("gene" != isolate(mode$mode) && "" != isolate(input$gene)) {
            mode$mode <- "gene"
            updateSelectizeInput(session = session, inputId = "protein", choices = proteins, selected="", server = TRUE)
        }
    })

    observeEvent(input$protein, ignoreInit = T, {
        if ("protein" != isolate(mode$mode) && "" != isolate(input$protein)) {
            mode$mode <- "protein"
            updateSelectizeInput(session = session, inputId = "gene", choices = genes, selected="", server = TRUE)
        }
    })
    
    observeEvent(input$cite_click, ignoreInit = T, {
      mode$mode <- "predicted_spatial"
      mode$click.x = input$cite_click$x
      mode$click.y = input$cite_click$y
      mode$cite.cells = closest.points.cite(mode$click.x, mode$click.y, n=1, max=5)
      #mode$cite.cells = closest.points.cite(mode$click.x, mode$click.y, n=10, max=5)
      updateSelectizeInput(session = session, inputId = "gene", choices = genes, selected="", server = TRUE)
      updateSelectizeInput(session = session, inputId = "protein", choices = proteins, selected="", server = TRUE)
    })
}
