library(shiny)

build.ui = function() {
  fluidPage(
    titlePanel('A Spatially-Resolved Single-Cell Multi-Omics Atlas of the Murine Spleen'),
    sidebarLayout(
      sidebarPanel(
        selectizeInput(
          inputId = 'gene',
          label = 'Gene',
          choices = NULL), 
        
        selectizeInput(
          inputId = 'protein',
          label = 'Protein',
          choices = NULL)
      ),
      
      mainPanel(
        fluidRow(style="font-size: medium;",
                 column(7, wellPanel(
                          tags$p("We profiled 7,097 cells from the murine spleen of 9-month BALBc mice with CITE-seq using a 30-antibody panel. The CITE-seq data was mapped with STvEA to three matching murine splenic sections profiled with CODEX in",
                                 tags$i("Goltsev, Samusik, et al., Cell 174, 968-981 (2018)."),
                                 "In addition we profiled 6,399 cells from the murine spleen of 2-month C57BL/6 mice using single-cell RNA-seq (10x)."),
                           tags$ul(
                           tags$li("Select a gene in the left panel to color the mRNA UMAP representations and the 
                                   splenic sections respectively with the measured and inferred expression levels of the gene."), tags$br(),
                           tags$li("Select a protein in the left panel to color the mRNA UMAP representation and the splenic
                                   sections of BALBc mice with the expression levels of the protein."), tags$br(),
                           tags$li("Click on a region in the BALBc mRNA UMAP representation to see the location of cells predicted
                                   to have a similar transcriptional profile in the three tissue sections.")), tags$br(),
                           tags$p("More information about STvEA can be found in:"),
                           tags$ul(tags$li("K. W. Govek*, E. C. Troisi*, Z. Miao, R. G. Aubin, S. Woodhouse, and P. G. Camara. ",
                                        tags$i("Single-Cell Transcriptomic Analysis of mIHC Images via Antigen Mapping. "),
                                        tags$b("Science Advances"),
                                        " 7 (2021) 10. ",
                                        a("DOI: 10.1126/sciadv.abc5464", href="https://doi.org/10.1126/sciadv.abc5464"),
                                        ". [*authors contributed equally]."))))))),
    
    fluidRow(column(6, h3('BALBc CITE-seq'),
                       img(src="clusters.svg", alt = "UMAP clusters", style="width:75%;")),
             column(6, h3('C57BL/6 scRNA-seq'),
                       img(src="C57BL.png", alt = "UMAP clusters", style="width:75%;"))),

    fluidRow(column(6, h3('BALBc CITE-seq mRNA UMAP'),
                       plotOutput('cite', click="cite_click")),
             column(6, h3('C57BL/6 scRNA-seq mRNA UMAP'),
                       plotOutput('cite_pilot'))),

    fluidRow(column(12, h3('BALBc CODEX spatial plots'))),
    fluidRow(
      column(3, plotOutput('codex1')),
      column(3, plotOutput('codex2')),
      column(3, plotOutput('codex3'))
    )
  )
}
