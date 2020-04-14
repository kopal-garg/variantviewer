library(shiny)
library(DT)
library(shinyFiles)
library(readr)
library(shinydashboard)
library(dplyr)
library(sparkline)
library(formattable)
library(DBI)
library(RMySQL)
library(shinyjs)
library(shinydashboardPlus)
library(shinyFeedback)




hpo_terms <- read_tsv('DB_COMPILED_HPO_PHENOTYPE_GENE')
colnames(hpo_terms)<-c("hp_id", "phen", "genes")
omim <- read_tsv('DB_OMIM_GENE_DISEASE')
omim <- subset(omim, select = c(GENE, DISEASE))
phenotype_list <- sort(cbind(hpo_terms$phen, omim$DISEASE))
phenotype_list<-as.data.frame(phenotype_list)
phenotype_list <- phenotype_list %>% filter(!grepl('Diamond-Blackfan', phenotype_list))
phenotype_list <- phenotype_list %>% filter(!grepl('Diamond-blackfan', phenotype_list))
phenotype_list <- phenotype_list %>% filter(!grepl('Diamond blackfan', phenotype_list))
phenotype_list <- phenotype_list %>% filter(!grepl('Diamond Blackfan', phenotype_list))
phenotype_list <- phenotype_list %>% filter(!grepl('Bone marrow failure', phenotype_list))
phenotype_list <- phenotype_list %>% filter(!grepl('Shwachman-Diamond syndrome', phenotype_list))
dba <- c()
dba$phenotype_list <- c('Diamond-Blackfan anemia', 'Bone marrow failure', 'Shwachman-Diamond syndrome')
dba <- as.data.frame(dba)
phenotype_list <- rbind(phenotype_list, dba)
colnames(phenotype_list) <- 'Select Phenotype'
gene_list <- read.table('gene_list.tsv', header = FALSE)
colnames(gene_list) <- 'Select Gene'
impact_list <- c('HIGH', 'MODERATE', 'LOW', 'MODIFIER')
ui <- dashboardPage(
  dashboardHeader(title = HTML(paste(''))
                  
                  
                  
                  
                  ),
  dashboardSidebar(width = 300,collapsed = FALSE,
                   
                   shinyjs::useShinyjs(),
                   fileInput('dataset', 'Choose TSV or TXT File',accept=c('text/tsv', 'text/tab-separated- values,text/plain', '.tsv', '.tsv', '.txt', 'TSV', '.TSV')),
                   selectInput(inputId = 'phenotype', choices = phenotype_list, label = 'Select Phenotype Keywords', multiple = TRUE, selectize = TRUE),
                   selectInput(inputId = 'gene',choices = gene_list, label ='Select Gene', multiple = TRUE, selectize = TRUE),
                   sliderInput(inputId = 'AF', label = 'AF', min = 0, max= 0.01, value = 0.01, width= 300),
                   sliderInput(inputId = 'AC', label = 'AC', min = 0, max = 150, value = 150, width = 300),
                   div(id="range_var_css",
                   sliderInput(inputId = 'SIFT', label = 'SIFT', min = 0, max = 1, value = c(0,1), step = 0.001,width = 300)),
                   div(id="range_var_css1",
                   sliderInput(inputId = 'PolyPhen', label = 'PolyPhen', min  = 0, max = 1, value = c(0,1), step =0.001,width = 300)),
                   div(id="range_var_css2",
                    sliderInput(inputId = 'spliceai', label ='SpliceAI', min = 0, max =1, value = c(0,1), width = 300)),
                   selectInput(inputId = 'impact_type', label = 'IMPACT', choices = impact_list, multiple = TRUE, selected = c('HIGH', 'LOW', 'MODERATE', 'MODIFIER')),
                   valueBoxOutput('num', width = 300),
                   actionButton(inputId = 'reset',label = 'Reset',width = 200)
),
  

  dashboardBody(
    
    tabBox(id= 'tables', width = 100,
           tabPanel('Autosomal De Novo',
                    DT::dataTableOutput(outputId = 'table_output_denovo')),
           tabPanel('Autosomal Recessive',
                    DT::dataTableOutput(outputId = 'table_output_recessive')),
           tabPanel('X-Linked De Novo',
                    DT::dataTableOutput(outputId = 'table_output_x_denovo')),
           tabPanel('X-Linked Recessive',
                    DT::dataTableOutput(outputId = 'table_output_x_recessive')),
           tabPanel('Compound Heterozygous',
                    DT::dataTableOutput(outputId = 'table_output_comp_het')),
           tabPanel('Based on Known Phenotype',
                    DT::dataTableOutput(outputId = 'table_output_known_genes')),
           tabPanel('Based on Gene(s) of Interest',
                    DT::dataTableOutput(outputId = 'table_based_on_input_genes'))
    )
  ) )
server <- function(input, output, session) {  
  
  #reset
  
  observeEvent(input$reset,{
    updateSliderInput(session,inputId = 'SIFT', value=c(0,1))
    updateSliderInput(session,inputId = 'PolyPhen', value=c(0,1))
    updateSliderInput(session, inputId = 'AF', value=0.01)
    updateSliderInput(session, inputId = 'AC', value = 150)
    updateSliderInput(session, inputId = 'spliceai', value = c(0,1))
    
  })

  #gradiant color for slider input

  observeEvent(input$SIFT, {
    if (input$SIFT[2] <= .4) {
      runjs(paste0('$("#range_var_css .irs-bar").css("background-color"," red")'))
    }

    if (input$SIFT[2] > .4 & input$SIFT[2] <= 1) {
      runjs(paste0('$("#range_var_css .irs-bar").css("background-color"," green")'))
    }
  })
  
  
  observeEvent(input$PolyPhen, {
    if (input$PolyPhen[1] >= 0 && input$PolyPhen[2]>=0.15) {
      runjs(paste0('$("#range_var_css1 .irs-bar").css("background-color"," green")'))
    }
    if (input$PolyPhen[1] >= .15) {
      runjs(paste0('$("#range_var_css1 .irs-bar").css("background-color"," orange")'))
    }
    
    if (input$PolyPhen[1] > .85) {
      runjs(paste0('$("#range_var_css1 .irs-bar").css("background-color"," red")'))
    }
  })
  
  observeEvent(input$spliceai, {
    if (input$spliceai[1] >= 0 && input$spliceai[2]>=0.2) {
      runjs(paste0('$("#range_var_css2 .irs-bar").css("background-color"," yellow")'))
    }
    if (input$spliceai[1] >= .2) {
      runjs(paste0('$("#range_var_css2 .irs-bar").css("background-color"," orange")'))
    }
    
    if (input$spliceai[1] > .8) {
      runjs(paste0('$("#range_var_css2 .irs-bar").css("background-color","red")'))
    }
  })

    
  
  
  

  output$num <- renderValueBox({
    
    a <- rbind(mydata_denovo(), mydata_recessive(), mydata_x_denovo(), mydata_x_recessive(), mydata_comp_het())
    a_dim <- dim(a)[1]
    valueBox(
      a_dim,
      'Total',
      color = 'light-blue'
    )
  })
  

  library(vcfR)

  #De Novo    
  mydata_denovo <- reactive({
    inFile <- input$dataset
    if (is.null(inFile))
      return(NULL)
    max_AF <- input$AF
    max_AC <- input$AC
    max_PolyPhen <- input$PolyPhen[2]
    min_PolyPhen <- input$PolyPhen[1]
    max_SpliceAI <- input$spliceai[2]
    min_SpliceAI <- input$spliceai[1]
    max_SIFT <- input$SIFT[2]
    min_SIFT <- input$SIFT[1]
    impact_type <- input$impact_type
    
    tbl <- read_tsv(inFile$datapath)
    colnames(tbl)[colnames(tbl)=="SIFT Score"]<-"SIFT"
    colnames(tbl)[colnames(tbl)=="PolyPhen Score"]<-"PolyPhen"
    colnames(tbl)[colnames(tbl)=="Interpretation"]<-"gnomad/ExAC"
    colnames(tbl)[colnames(tbl)=="Variant"]<-"Coordinate"
    colnames(tbl)[colnames(tbl)=="SYMBOL_link"]<-"Symbol"
    colnames(tbl)[colnames(tbl)=="Function"]<-"Type"
    colnames(tbl)[colnames(tbl)=="gnomad_popmax_af"]<-"AF"
    colnames(tbl)[colnames(tbl)=="gene_description_2"]<-"Gene Description"
    colnames(tbl)[colnames(tbl)=="pli"]<-"pLI"
    colnames(tbl)[colnames(tbl)=="oe"]<-"o/e"
    denovo<-tbl%>%filter(grepl('^denovo$',mode))
    denovo<-subset(denovo, select = c(Coordinate, Symbol, AF, AC,`REF:ALT`, Type,IMPACT, `HGVSc/p`,SIFT,PolyPhen,spliceai, `gnomad/ExAC`,pLI, `o/e`,`Gene Description` ))
    denovo$SIFT <- as.numeric(denovo$SIFT)
    denovo$PolyPhen <- as.numeric(denovo$PolyPhen)
    denovo$spliceai <- as.numeric(denovo$spliceai)
    denovo$pLI <- as.numeric(denovo$pLI)
    denovo$`o/e` <- as.numeric(denovo$`o/e`)
    denovo$AC <- as.numeric(denovo$AC)
    denovo <- denovo %>% filter((AC <= max_AC))
    denovo <- denovo %>% filter((AF <= max_AF))
    
    #SIFT
    na_values_SIFT <- denovo %>% filter(is.na(SIFT))
    SIFT <- denovo %>% filter(SIFT >= min_SIFT) %>% filter(SIFT <= max_SIFT)
    SIFT <- rbind(SIFT, na_values_SIFT)
    denovo <- intersect(SIFT, denovo)
    
    #PolyPhen
    na_values_PolyPhen <- denovo %>% filter(is.na(PolyPhen))
    PolyPhen <- denovo %>% filter(PolyPhen >= min_PolyPhen) %>% filter(PolyPhen <= max_PolyPhen)
    PolyPhen <- rbind(PolyPhen, na_values_PolyPhen)
    denovo <- intersect(PolyPhen, denovo)
    
    #SpliceAI
    na_values_spliceai <- denovo %>% filter(is.na(spliceai))
    spliceai <- denovo %>% filter(spliceai >= min_SpliceAI) %>% filter(spliceai <= max_SpliceAI)
    spliceai <- rbind(spliceai, na_values_spliceai)
    denovo <- intersect(spliceai, denovo)
    
    #Impact Type
    impact_type <- input$impact_type
    impact <- c()
    for (i in 1:length(impact_type)){
    a <- denovo %>% filter(grepl(impact_type[i], IMPACT))
    impact <- rbind(impact, a)
    i = i+1
    }
    denovo <- impact
    return(denovo)
  })
  
  
  color.picker.SIFT <- function(z){
    if(is.na(z)){return("white")}
    else if(z <= 0.05){return("red")}
    else if( z > 0.05){return("lightgreen")}
    else {return("black")}
  }
  
  color.picker.PolyPhen <- function(z){
    if(is.na(z)){return("white")}
    else if(z <= 0.15){return("lightgreen")}
    else if( z > 0.15 && z < 0.85){return("orange")}
    else {return("red")}
  }
  color.picker.Type <- function(z){
    if(is.na(z)){return("white")}
    else if(z == "HIGH"){return("red")}
    else if(z == "MODERATE"){return("orange")}
    else if( z == 'MODIFIER' || z == "LOW"){return("lightblue")}
    else {return("white")}
  }
  
  
  output$table_output_denovo <- renderDataTable({
    tryCatch({
    as.datatable(formattable(mydata_denovo(), list(
      SIFT = formatter("span",
                       style = x ~ style(
                         display = "block",
                         "border-radius" = "4px",
                         "padding-right" = "4px",
                         "background-color" = sapply(x,color.picker.SIFT)
                       )),
      PolyPhen = formatter("span",
                           style = x ~ style(
                             display = "block",
                             "border-radius" = "4px",
                             "padding-right" = "4px",
                             "background-color" = sapply(x,color.picker.PolyPhen))),
      IMPACT = formatter('span',
                         style = x ~ style(
                           display = "block",
                           "border-radius" = "4px",
                           "background-color" = sapply(x,color.picker.Type)))
    )), caption =htmltools::tags$caption(
      style = 'caption-side: bottom; text-align: left;', htmltools::em('SIFT: 0.0 to 0.05 (deleterious) to 0.05 to 1.0 (tolerated);', 'PolyPhen-2: 0.0 to 0.15 (tolerated) to 0.15 to 0.85 (possibly damaging) to 0.85 to 1.0 (damaging);', 'SpliceAI: 0.2 (high recall/likely pathogenic), 0.5 (recommended/pathogenic), and 0.8 (high precision/pathogenic) cutoffs;','o/e: low oe values are indicative of strong intolerance;','pLI >= 0.9 are an extremely intolerant set of transcripts')),
    extensions = c('Buttons','ColReorder','RowReorder','FixedColumns'),escape = FALSE,options = list(rowReorder = TRUE,colReorder = TRUE,paging  = FALSE,
                                                                                                                    fixedColumns = TRUE,autowidth=TRUE, searchHighlight=TRUE,dom = 'Bfrtip',scrollX = TRUE,
                                                                                                                    buttons = 
                                                                                                                      list(
                                                                                                                        extend = 'collection',
                                                                                                                        buttons = c('csv', 'excel', 'pdf'),
                                                                                                                        text = 'Download'),I('colvis'))
    )
      
    
    }, error=function(e){})
  })
  
  
  
  #Autosomal Recessive 
  mydata_recessive <- reactive({
    inFile <- input$dataset
    a <- c()
    if (is.null(inFile))
      return(a)
    max_AF <- input$AF
    max_AC <- input$AC
    max_PolyPhen <- input$PolyPhen[2]
    min_PolyPhen <- input$PolyPhen[1]
    max_SpliceAI <- input$spliceai[2]
    min_SpliceAI <- input$spliceai[1]
    max_SIFT <- input$SIFT[2]
    min_SIFT <- input$SIFT[1]
    tbl <- read_tsv(inFile$datapath)
    colnames(tbl)[colnames(tbl)=="SIFT Score"]<-"SIFT"
    colnames(tbl)[colnames(tbl)=="PolyPhen Score"]<-"PolyPhen"
    colnames(tbl)[colnames(tbl)=="Interpretation"]<-"gnomad/ExAC"
    colnames(tbl)[colnames(tbl)=="Variant"]<-"Coordinate"
    colnames(tbl)[colnames(tbl)=="SYMBOL_link"]<-"Symbol"
    colnames(tbl)[colnames(tbl)=="Function"]<-"Type"
    colnames(tbl)[colnames(tbl)=="gnomad_popmax_af"]<-"AF"
    colnames(tbl)[colnames(tbl)=="gene_description_2"]<-"Gene Description"
    colnames(tbl)[colnames(tbl)=="pli"]<-"pLI"
    colnames(tbl)[colnames(tbl)=="oe"]<-"o/e"
    recessive<-tbl%>%filter(grepl('^recessive$',mode))
    recessive<-subset(recessive, select = c(Coordinate, Symbol, AF, AC,`REF:ALT`, Type,IMPACT, `HGVSc/p`,SIFT,PolyPhen,spliceai, `gnomad/ExAC`,pLI, `o/e`,`Gene Description` ))
    
    recessive$SIFT <- as.numeric(recessive$SIFT)
    recessive$PolyPhen <- as.numeric(recessive$PolyPhen)
    recessive$spliceai <- as.numeric(recessive$spliceai)
    recessive$pLI <- as.numeric(recessive$pLI)
    recessive$`o/e` <- as.numeric(recessive$`o/e`)
    recessive <- recessive %>% filter((AC <= max_AC))
    recessive <- recessive %>% filter((AF <= max_AF))
    
    #SIFT
    na_values_SIFT <- recessive %>% filter(is.na(SIFT))
    SIFT <- recessive %>% filter(SIFT >= min_SIFT) %>% filter(SIFT <= max_SIFT)
    SIFT <- rbind(SIFT, na_values_SIFT)
    recessive <- intersect(SIFT, recessive)
    
    #PolyPhen
    na_values_PolyPhen <- recessive %>% filter(is.na(PolyPhen))
    PolyPhen <- recessive %>% filter(PolyPhen >= min_PolyPhen) %>% filter(PolyPhen <= max_PolyPhen)
    PolyPhen <- rbind(PolyPhen, na_values_PolyPhen)
    recessive <- intersect(PolyPhen, recessive)
    
    #SpliceAI
    na_values_spliceai <- recessive %>% filter(is.na(spliceai))
    spliceai <- recessive %>% filter(spliceai >= min_SpliceAI) %>% filter(spliceai <= max_SpliceAI)
    spliceai <- rbind(spliceai, na_values_spliceai)
    recessive <- intersect(spliceai, recessive)
    
    #Impact Type
    impact_type <- input$impact_type
    impact <- c()
    for (i in 1:length(impact_type)){
      a <- recessive %>% filter(grepl(impact_type[i], IMPACT))
      impact <- rbind(impact, a)
      i = i+1
    }
    recessive <- impact
    return(recessive)
  })
  
  output$table_output_recessive <- renderDataTable({
    tryCatch({
    as.datatable(formattable(mydata_recessive(), list(
      SIFT = formatter("span",
                       style = x ~ style(
                         display = "block",
                         "border-radius" = "4px",
                         "padding-right" = "4px",
                         "background-color" = sapply(x,color.picker.SIFT)
                       )),
      PolyPhen = formatter("span",
                           style = x ~ style(
                             display = "block",
                             "border-radius" = "4px",
                             "padding-right" = "4px",
                             "background-color" = sapply(x,color.picker.PolyPhen))),
      IMPACT = formatter('span',
                         style = x ~ style(
                           display = "block",
                           "border-radius" = "4px",
                           "background-color" = sapply(x,color.picker.Type)))
    )), caption =htmltools::tags$caption(
      style = 'caption-side: bottom; text-align: left;', htmltools::em('SIFT: 0.0 to 0.05 (deleterious) to 0.05 to 1.0 (tolerated);', 'PolyPhen-2: 0.0 to 0.15 (tolerated) to 0.15 to 0.85 (possibly damaging) to 0.85 to 1.0 (damaging);', 'SpliceAI: 0.2 (high recall/likely pathogenic), 0.5 (recommended/pathogenic), and 0.8 (high precision/pathogenic) cutoffs;','o/e: low oe values are indicative of strong intolerance;','pLI >= 0.9 are an extremely intolerant set of transcripts')),
    extensions = c('Buttons','ColReorder','RowReorder','FixedColumns'),escape = FALSE,options = list(rowReorder = TRUE,colReorder = TRUE,paging  = FALSE,
                                                                                                                    fixedColumns = TRUE,autowidth=TRUE, searchHighlight=TRUE,dom = 'Bfrtip',scrollX = TRUE,
                                                                                                                    buttons = 
                                                                                                                      list(
                                                                                                                        extend = 'collection',
                                                                                                                        buttons = c('csv', 'excel', 'pdf'),
                                                                                                                        text = 'Download'),I('colvis'))
    )
    
    
    }, error=function(e){})
    
  })
  # X-Linked De Novo
  mydata_x_denovo <- reactive({
    inFile <- input$dataset
    if (is.null(inFile))
      return(NULL)
    max_AF <- input$AF
    max_AC <- input$AC
    max_PolyPhen <- input$PolyPhen[2]
    min_PolyPhen <- input$PolyPhen[1]
    max_SpliceAI <- input$spliceai[2]
    min_SpliceAI <- input$spliceai[1]
    max_SIFT <- input$SIFT[2]
    min_SIFT <- input$SIFT[1]
    tbl <- read_tsv(inFile$datapath)
    colnames(tbl)[colnames(tbl)=="SIFT Score"]<-"SIFT"
    colnames(tbl)[colnames(tbl)=="PolyPhen Score"]<-"PolyPhen"
    colnames(tbl)[colnames(tbl)=="Interpretation"]<-"gnomad/ExAC"
    colnames(tbl)[colnames(tbl)=="Variant"]<-"Coordinate"
    colnames(tbl)[colnames(tbl)=="SYMBOL_link"]<-"Symbol"
    colnames(tbl)[colnames(tbl)=="Function"]<-"Type"
    colnames(tbl)[colnames(tbl)=="gnomad_popmax_af"]<-"AF"
    colnames(tbl)[colnames(tbl)=="gene_description_2"]<-"Gene Description"
    colnames(tbl)[colnames(tbl)=="pli"]<-"pLI"
    colnames(tbl)[colnames(tbl)=="oe"]<-"o/e"
    x_denovo<-tbl%>%filter(grepl('^x_denovo$',mode))
    x_denovo<-subset(x_denovo, select = c(Coordinate, Symbol, AF, AC,`REF:ALT`, Type, IMPACT,`HGVSc/p`,SIFT,PolyPhen,spliceai, `gnomad/ExAC`,pLI, `o/e`,`Gene Description` ))
    x_denovo$SIFT <- as.numeric(x_denovo$SIFT)
    x_denovo$PolyPhen <- as.numeric(x_denovo$PolyPhen)
    x_denovo$spliceai <- as.numeric(x_denovo$spliceai)
    x_denovo$pLI <- as.numeric(x_denovo$pLI)
    x_denovo$`o/e` <- as.numeric(x_denovo$`o/e`)
    x_denovo <- x_denovo %>% filter((AC <= max_AC))
    x_denovo <- x_denovo %>% filter((AF <= max_AF))
    
    #SIFT
    na_values_SIFT <- x_denovo %>% filter(is.na(SIFT))
    SIFT <- x_denovo %>% filter(SIFT >= min_SIFT) %>% filter(SIFT <= max_SIFT)
    SIFT <- rbind(SIFT, na_values_SIFT)
    x_denovo <- intersect(SIFT, x_denovo)
    
    #PolyPhen
    na_values_PolyPhen <- x_denovo %>% filter(is.na(PolyPhen))
    PolyPhen <- x_denovo %>% filter(PolyPhen >= min_PolyPhen) %>% filter(PolyPhen <= max_PolyPhen)
    PolyPhen <- rbind(PolyPhen, na_values_PolyPhen)
    x_denovo <- intersect(PolyPhen, x_denovo)
    
    #SpliceAI
    na_values_spliceai <- x_denovo %>% filter(is.na(spliceai))
    spliceai <- x_denovo %>% filter(spliceai >= min_SpliceAI) %>% filter(spliceai <= max_SpliceAI)
    spliceai <- rbind(spliceai, na_values_spliceai)
    x_denovo <- intersect(spliceai, x_denovo)
    
    #Impact Type
    impact_type <- input$impact_type
    impact <- c()
    for (i in 1:length(impact_type)){
      a <- x_denovo %>% filter(grepl(impact_type[i], IMPACT))
      impact <- rbind(impact, a)
      i = i+1
    }
    x_denovo <- impact
    return(x_denovo)
  })
  
  output$table_output_x_denovo <- renderDataTable({
    tryCatch({
    as.datatable(formattable(mydata_x_denovo(), list(
      SIFT = formatter("span",
                       style = x ~ style(
                         display = "block",
                         "border-radius" = "4px",
                         "padding-right" = "4px",
                         "background-color" = sapply(x,color.picker.SIFT)
                       )),
      PolyPhen = formatter("span",
                           style = x ~ style(
                             display = "block",
                             "border-radius" = "4px",
                             "padding-right" = "4px",
                             "background-color" = sapply(x,color.picker.PolyPhen))),
      IMPACT = formatter('span',
                         style = x ~ style(
                           display = "block",
                           "border-radius" = "4px",
                           "background-color" = sapply(x,color.picker.Type)))
    )), caption =htmltools::tags$caption(
      style = 'caption-side: bottom; text-align: left;', htmltools::em('SIFT: 0.0 to 0.05 (deleterious) to 0.05 to 1.0 (tolerated);', 'PolyPhen-2: 0.0 to 0.15 (tolerated) to 0.15 to 0.85 (possibly damaging) to 0.85 to 1.0 (damaging);', 'SpliceAI: 0.2 (high recall/likely pathogenic), 0.5 (recommended/pathogenic), and 0.8 (high precision/pathogenic) cutoffs;','o/e: low oe values are indicative of strong intolerance;','pLI >= 0.9 are an extremely intolerant set of transcripts')),
    extensions = c('Buttons','ColReorder','RowReorder','FixedColumns'),escape = FALSE,options = list(rowReorder = TRUE,colReorder = TRUE,paging  = FALSE,
                                                                                                                    fixedColumns = TRUE,autowidth=TRUE, searchHighlight=TRUE,dom = 'Bfrtip',scrollX = TRUE,
                                                                                                                    buttons = 
                                                                                                                      list(
                                                                                                                        extend = 'collection',
                                                                                                                        buttons = c('csv', 'excel', 'pdf'),
                                                                                                                        text = 'Download'),I('colvis'))
    )
    
    
    }, error=function(e){})
    
  })
  
  #X-Linked Recessive 
  mydata_x_recessive <- reactive({
    inFile <- input$dataset
    if (is.null(inFile))
      return(NULL)
    max_AF <- input$AF
    max_AC <- input$AC
    max_PolyPhen <- input$PolyPhen[2]
    min_PolyPhen <- input$PolyPhen[1]
    max_SpliceAI <- input$spliceai[2]
    min_SpliceAI <- input$spliceai[1]
    max_SIFT <- input$SIFT[2]
    min_SIFT <- input$SIFT[1]
    tbl <- read_tsv(inFile$datapath)
    colnames(tbl)[colnames(tbl)=="SIFT Score"]<-"SIFT"
    colnames(tbl)[colnames(tbl)=="PolyPhen Score"]<-"PolyPhen"
    colnames(tbl)[colnames(tbl)=="Interpretation"]<-"gnomad/ExAC"
    colnames(tbl)[colnames(tbl)=="Variant"]<-"Coordinate"
    colnames(tbl)[colnames(tbl)=="SYMBOL_link"]<-"Symbol"
    colnames(tbl)[colnames(tbl)=="Function"]<-"Type"
    colnames(tbl)[colnames(tbl)=="gnomad_popmax_af"]<-"AF"
    colnames(tbl)[colnames(tbl)=="gene_description_2"]<-"Gene Description"
    colnames(tbl)[colnames(tbl)=="pli"]<-"pLI"
    colnames(tbl)[colnames(tbl)=="oe"]<-"o/e"
    x_recessive<-tbl%>%filter(grepl("^x_recessive$", mode))
    x_recessive<-subset(x_recessive, select = c(Coordinate, Symbol, AF, AC,`REF:ALT`, Type,IMPACT, `HGVSc/p`,SIFT,PolyPhen,spliceai, `gnomad/ExAC`,pLI, `o/e`,`Gene Description` ))
    
    x_recessive$SIFT <- as.numeric(x_recessive$SIFT)
    x_recessive$PolyPhen <- as.numeric(x_recessive$PolyPhen)
    x_recessive$spliceai <- as.numeric(x_recessive$spliceai)
    x_recessive$pLI <- as.numeric(x_recessive$pLI)
    x_recessive$`o/e` <- as.numeric(x_recessive$`o/e`)
    x_recessive <- x_recessive %>% filter((AC <= max_AC))
    x_recessive <- x_recessive %>% filter((AF <= max_AF))
   
    #SIFT
    na_values_SIFT <- x_recessive %>% filter(is.na(SIFT))
    SIFT <- x_recessive %>% filter(SIFT >= min_SIFT) %>% filter(SIFT <= max_SIFT)
    SIFT <- rbind(SIFT, na_values_SIFT)
    x_recessive <- intersect(SIFT, x_recessive)
    
    #PolyPhen
    na_values_PolyPhen <- x_recessive %>% filter(is.na(PolyPhen))
    PolyPhen <- x_recessive %>% filter(PolyPhen >= min_PolyPhen) %>% filter(PolyPhen <= max_PolyPhen)
    PolyPhen <- rbind(PolyPhen, na_values_PolyPhen)
    x_recessive <- intersect(PolyPhen, x_recessive)
    
    #SpliceAI
    na_values_spliceai <- x_recessive %>% filter(is.na(spliceai))
    spliceai <- x_recessive %>% filter(spliceai >= min_SpliceAI) %>% filter(spliceai <= max_SpliceAI)
    spliceai <- rbind(spliceai, na_values_spliceai)
    x_recessive <- intersect(spliceai, x_recessive)
    
    #Impact Type
    impact_type <- input$impact_type
    impact <- c()
    for (i in 1:length(impact_type)){
      a <- x_recessive %>% filter(grepl(impact_type[i], IMPACT))
      impact <- rbind(impact, a)
      i = i+1
    }
    x_recessive <- impact
    return(x_recessive)
  })
  
  output$table_output_x_recessive <- renderDataTable({
    tryCatch({
    as.datatable(formattable(mydata_x_recessive(), list(
      SIFT = formatter("span",
                       style = x ~ style(
                         display = "block",
                         "border-radius" = "4px",
                         "padding-right" = "4px",
                         "background-color" = sapply(x,color.picker.SIFT)
                       )),
      PolyPhen = formatter("span",
                           style = x ~ style(
                             display = "block",
                             "border-radius" = "4px",
                             "padding-right" = "4px",
                             "background-color" = sapply(x,color.picker.PolyPhen))),
      IMPACT = formatter('span',
                         style = x ~ style(
                           display = "block",
                           "border-radius" = "4px",
                           "background-color" = sapply(x,color.picker.Type)))
    )), caption =htmltools::tags$caption(
      style = 'caption-side: bottom; text-align: left;', htmltools::em('SIFT: 0.0 to 0.05 (deleterious) to 0.05 to 1.0 (tolerated);', 'PolyPhen-2: 0.0 to 0.15 (tolerated) to 0.15 to 0.85 (possibly damaging) to 0.85 to 1.0 (damaging);', 'SpliceAI: 0.2 (high recall/likely pathogenic), 0.5 (recommended/pathogenic), and 0.8 (high precision/pathogenic) cutoffs;','o/e: low oe values are indicative of strong intolerance;','pLI >= 0.9 are an extremely intolerant set of transcripts')),
    extensions = c('Buttons','ColReorder','RowReorder','FixedColumns'),escape = FALSE,options = list(rowReorder = TRUE,colReorder = TRUE,paging  = FALSE,
                                                                                                                    fixedColumns = TRUE,autowidth=TRUE, searchHighlight=TRUE,dom = 'Bfrtip',scrollX = TRUE,
                                                                                                                    buttons = 
                                                                                                                      list(
                                                                                                                        extend = 'collection',
                                                                                                                        buttons = c('csv', 'excel', 'pdf'),
                                                                                                                        text = 'Download'),I('colvis'))
    )
    
    }, error=function(e){})
    
    
  })
  
  #Compound Het 
  mydata_comp_het <- reactive({
    inFile <- input$dataset
    if (is.null(inFile))
      return(NULL)
    max_AF <- input$AF
    max_AC <- input$AC
    max_PolyPhen <- input$PolyPhen[2]
    min_PolyPhen <- input$PolyPhen[1]
    max_SpliceAI <- input$spliceai[2]
    min_SpliceAI <- input$spliceai[1]
    max_SIFT <- input$SIFT[2]
    min_SIFT <- input$SIFT[1]
    tbl <- read_tsv(inFile$datapath)
    colnames(tbl)[colnames(tbl)=="SIFT Score"]<-"SIFT"
    colnames(tbl)[colnames(tbl)=="PolyPhen Score"]<-"PolyPhen"
    colnames(tbl)[colnames(tbl)=="Interpretation"]<-"gnomad/ExAC"
    colnames(tbl)[colnames(tbl)=="Variant"]<-"Coordinate"
    colnames(tbl)[colnames(tbl)=="SYMBOL_link"]<-"Symbol"
    colnames(tbl)[colnames(tbl)=="Function"]<-"Type"
    colnames(tbl)[colnames(tbl)=="gnomad_popmax_af"]<-"AF"
    colnames(tbl)[colnames(tbl)=="gene_description_2"]<-"Gene Description"
    colnames(tbl)[colnames(tbl)=="pli"]<-"pLI"
    colnames(tbl)[colnames(tbl)=="oe"]<-"o/e"
    comp_het<-tbl%>%filter(grepl("^comp_het$", mode))
    comp_het<-subset(comp_het, select = c(Coordinate, Symbol, AF, AC,`REF:ALT`, Type,IMPACT, `HGVSc/p`,SIFT,PolyPhen,spliceai, `gnomad/ExAC`,pLI, `o/e`,`Gene Description` ))
    
    
    
    comp_het$SIFT <- as.numeric(comp_het$SIFT)
    comp_het$PolyPhen <- as.numeric(comp_het$PolyPhen)
    comp_het$spliceai <- as.numeric(comp_het$spliceai)
    comp_het$pLI <- as.numeric(comp_het$pLI)
    comp_het$`o/e` <- as.numeric(comp_het$`o/e`)
    comp_het <- comp_het %>% filter((AC <= max_AC))
    comp_het <- comp_het %>% filter((AF <= max_AF))
    
    #SIFT
    na_values_SIFT <- comp_het %>% filter(is.na(SIFT))
    SIFT <- comp_het %>% filter(SIFT >= min_SIFT) %>% filter(SIFT <= max_SIFT)
    SIFT <- rbind(SIFT, na_values_SIFT)
    comp_het <- intersect(SIFT, comp_het)
    
    #PolyPhen
    na_values_PolyPhen <- comp_het %>% filter(is.na(PolyPhen))
    PolyPhen <- comp_het %>% filter(PolyPhen >= min_PolyPhen) %>% filter(PolyPhen <= max_PolyPhen)
    PolyPhen <- rbind(PolyPhen, na_values_PolyPhen)
    comp_het <- intersect(PolyPhen, comp_het)
    
    #SpliceAI
    na_values_spliceai <- comp_het %>% filter(is.na(spliceai))
    spliceai <- comp_het %>% filter(spliceai >= min_SpliceAI) %>% filter(spliceai <= max_SpliceAI)
    spliceai <- rbind(spliceai, na_values_spliceai)
    comp_het <- intersect(spliceai, comp_het)
    
    #Impact Type
    impact_type <- input$impact_type
    impact <- c()
    for (i in 1:length(impact_type)){
      a <- comp_het %>% filter(grepl(impact_type[i], IMPACT))
      impact <- rbind(impact, a)
      i = i+1
    }
    comp_het <- impact
    return(comp_het)
  })
  #Associated Phenotype
  mydata_associated_genes <- reactive({
    
    pheno <- input$phenotype

    max_AF <- input$AF
    max_AC <- input$AC
    max_PolyPhen <- input$PolyPhen[2]
    min_PolyPhen <- input$PolyPhen[1]
    max_SpliceAI <- input$spliceai[2]
    min_SpliceAI <- input$spliceai[1]
    max_SIFT <- input$SIFT[2]
    min_SIFT <- input$SIFT[1]
    hpo_terms <- read_tsv('DB_COMPILED_HPO_PHENOTYPE_GENE')
    colnames(hpo_terms)<-c("hp_id", "phen", "genes")
    hpo_terms <- subset(hpo_terms, select = c('genes', 'phen'))
    omim <- read_tsv('DB_OMIM_GENE_DISEASE')
    omim <- subset(omim, select = c(GENE, DISEASE))
    omim$GENE<-gsub(x=omim$GENE, pattern = ', ', replacement = ',')
    colnames(omim)<-c("genes", "phen")
    hpo_terms <- unique(rbind(hpo_terms, omim))
    hpo_subset <- c()

    
        
    for (i in 1:length(pheno)){
      if (pheno[i] == 'Diamond-Blackfan anemia'){
        
        a <- hpo_terms %>% filter(grepl('Diamond-Blackfan',phen))
        hpo_subset <- rbind(hpo_subset, a)
        i = i+1 
        
      }
      else if (pheno[i] == 'Bone marrow failure'){
        a <- hpo_terms %>% filter(grepl('Bone marrow failure', phen))
        hpo_subset <- rbind(hpo_subset, a)
        i = i+1 
      }
      else if (pheno[i] == 'Shwachman-Diamond syndrome'){
        a <- hpo_terms %>% filter(grepl('Shwachman-Diamond syndrome', phen))
        hpo_subset <- rbind(hpo_subset, a)
      }
      
      else {
        
        a <- hpo_terms %>% filter(grepl(pheno[i], phen, ignore.case =TRUE, fixed = TRUE))
        hpo_subset <- rbind(hpo_subset, a)
        i = i+1
        
        
      }

    }
    known_gene_list  <- unique(hpo_subset$genes)
    known_gene_list <- strsplit(known_gene_list, split=",") %>% unlist()
    inFile <- input$dataset
    if (is.null(inFile))
      return(NULL)
    df <- read_tsv(inFile$datapath)
    candidate_gene_list <- unique(df$SYMBOL_orig)
    common_genes <- intersect(candidate_gene_list, known_gene_list)
    df_new <- c()
    if (length(common_genes!=0)){
      for (i in 1:length(common_genes)){
        a <- df %>% filter(grepl(common_genes[i], SYMBOL_orig))
        df_new <- rbind(df_new,a)
        i = i+1 }} 
    tbl <- df_new
    colnames(tbl)[colnames(tbl)=="SIFT Score"]<-"SIFT"
    colnames(tbl)[colnames(tbl)=="PolyPhen Score"]<-"PolyPhen"
    colnames(tbl)[colnames(tbl)=="Interpretation"]<-"gnomad/ExAC"
    colnames(tbl)[colnames(tbl)=="Variant"]<-"Coordinate"
    colnames(tbl)[colnames(tbl)=="SYMBOL_link"]<-"Symbol"
    colnames(tbl)[colnames(tbl)=="Function"]<-"Type"
    colnames(tbl)[colnames(tbl)=="gnomad_popmax_af"]<-"AF"
    colnames(tbl)[colnames(tbl)=="gene_description_2"]<-"Gene Description"
    colnames(tbl)[colnames(tbl)=="pli"]<-"pLI"
    colnames(tbl)[colnames(tbl)=="oe"]<-"o/e"
    tbl$mode <- gsub(x=tbl$mode, pattern='^denovo$',replacement = 'De Novo')
    tbl$mode <- gsub(x=tbl$mode, pattern='^recessive$', replacement = 'Autosomal Recessive')
    tbl$mode <- gsub(x=tbl$mode, pattern='^x_denovo$', replacement = 'X-Linked Denovo')
    tbl$mode <- gsub(x=tbl$mode, pattern='^x_recessive$', replacement = 'X-Linked Recessive')
    tbl$mode <- gsub(x=tbl$mode, pattern='^comp_het$', replacement = 'Compound Heterozygote')
    
    colnames(tbl)[colnames(tbl)=='mode']<-'Mode'
    tbl<-subset(tbl, select = c(Mode,Coordinate, Symbol, AF, AC,`REF:ALT`, Type,IMPACT, `HGVSc/p`,SIFT,PolyPhen,spliceai, `gnomad/ExAC`,pLI, `o/e`,`Gene Description` ))
    tbl$SIFT <- as.numeric(tbl$SIFT)
    tbl$PolyPhen <- as.numeric(tbl$PolyPhen)
    tbl$spliceai <- as.numeric(tbl$spliceai)
    tbl$pLI <- as.numeric(tbl$pLI)
    tbl$`o/e` <- as.numeric(tbl$`o/e`)
    tbl <- tbl %>% filter((AC <= max_AC))
    tbl <- tbl %>% filter((AF <= max_AF))
    
    #SIFT
    na_values_SIFT <- tbl %>% filter(is.na(SIFT))
    SIFT <- tbl %>% filter(SIFT >= min_SIFT) %>% filter(SIFT <= max_SIFT)
    SIFT <- rbind(SIFT, na_values_SIFT)
    tbl <- intersect(SIFT, tbl)
    
    #PolyPhen
    na_values_PolyPhen <- tbl %>% filter(is.na(PolyPhen))
    PolyPhen <- tbl %>% filter(PolyPhen >= min_PolyPhen) %>% filter(PolyPhen <= max_PolyPhen)
    PolyPhen <- rbind(PolyPhen, na_values_PolyPhen)
    tbl <- intersect(PolyPhen, tbl)
    
    #SpliceAI
    na_values_spliceai <- tbl %>% filter(is.na(spliceai))
    spliceai <- tbl %>% filter(spliceai >= min_SpliceAI) %>% filter(spliceai <= max_SpliceAI)
    spliceai <- rbind(spliceai, na_values_spliceai)
    tbl <- intersect(spliceai, tbl)
    
    #Impact Type
    impact_type <- input$impact_type
    impact <- c()
    for (i in 1:length(impact_type)){
      a <- tbl %>% filter(grepl(impact_type[i], IMPACT))
      impact <- rbind(impact, a)
      i = i+1
    }
    tbl <- impact
    return(tbl)
    
    
    return(tbl)
  })
  
  #Associated Gene
  known_genes <- reactive({
    gene <- input$gene
    max_AF <- input$AF
    max_AC <- input$AC
    max_PolyPhen <- input$PolyPhen[2]
    min_PolyPhen <- input$PolyPhen[1]
    max_SpliceAI <- input$spliceai[2]
    min_SpliceAI <- input$spliceai[1]
    max_SIFT <- input$SIFT[2]
    min_SIFT <- input$SIFT[1]
    data <- input$dataset
    data <- read_tsv(data$datapath)
    gene_subset <- c()
    gene_list <- read_tsv('gene_list.tsv')
    colnames(gene_list) <-""
    for (i in 1:length(gene)){
      a <- data %>% filter(grepl(gene[i], SYMBOL_orig))
      gene_subset <- rbind(a, gene_subset)
      i = i+1
    }
    tbl <- gene_subset
    colnames(tbl)[colnames(tbl)=="SIFT Score"]<-"SIFT"
    colnames(tbl)[colnames(tbl)=="PolyPhen Score"]<-"PolyPhen"
    colnames(tbl)[colnames(tbl)=="Interpretation"]<-"gnomad/ExAC"
    colnames(tbl)[colnames(tbl)=="Variant"]<-"Coordinate"
    colnames(tbl)[colnames(tbl)=="SYMBOL_link"]<-"Symbol"
    colnames(tbl)[colnames(tbl)=="Function"]<-"Type"
    colnames(tbl)[colnames(tbl)=="gnomad_popmax_af"]<-"AF"
    colnames(tbl)[colnames(tbl)=="gene_description_2"]<-"Gene Description"
    colnames(tbl)[colnames(tbl)=="pli"]<-"pLI"
    colnames(tbl)[colnames(tbl)=="oe"]<-"o/e"
    tbl$mode <- gsub(x=tbl$mode, pattern='^denovo$',replacement = 'De Novo')
    tbl$mode <- gsub(x=tbl$mode, pattern='^recessive$', replacement = 'Autosomal Recessive')
    tbl$mode <- gsub(x=tbl$mode, pattern='^x_denovo$', replacement = 'X-Linked Denovo')
    tbl$mode <- gsub(x=tbl$mode, pattern='^x_recessive$', replacement = 'X-Linked Recessive')
    tbl$mode <- gsub(x=tbl$mode, pattern='^comp_het$', replacement = 'Compound Heterozygote')
    
    colnames(tbl)[colnames(tbl)=='mode']<-'Mode'
    tbl<-subset(tbl, select = c(Mode,Coordinate, Symbol, AF, AC,`REF:ALT`, Type,IMPACT, `HGVSc/p`,SIFT,PolyPhen,spliceai, `gnomad/ExAC`,pLI, `o/e`,`Gene Description` ))
    tbl$SIFT <- as.numeric(tbl$SIFT)
    tbl$PolyPhen <- as.numeric(tbl$PolyPhen)
    tbl$spliceai <- as.numeric(tbl$spliceai)
    tbl$pLI <- as.numeric(tbl$pLI)
    tbl$`o/e` <- as.numeric(tbl$`o/e`)
    tbl <- tbl %>% filter((AC <= max_AC))
    tbl <- tbl %>% filter((AF <= max_AF))
    
    #SIFT
    na_values_SIFT <- tbl %>% filter(is.na(SIFT))
    SIFT <- tbl %>% filter(SIFT >= min_SIFT) %>% filter(SIFT <= max_SIFT)
    SIFT <- rbind(SIFT, na_values_SIFT)
    tbl <- intersect(SIFT, tbl)
    
    #PolyPhen
    na_values_PolyPhen <- tbl %>% filter(is.na(PolyPhen))
    PolyPhen <- tbl %>% filter(PolyPhen >= min_PolyPhen) %>% filter(PolyPhen <= max_PolyPhen)
    PolyPhen <- rbind(PolyPhen, na_values_PolyPhen)
    tbl <- intersect(PolyPhen, tbl)
    
    #SpliceAI
    na_values_spliceai <- tbl %>% filter(is.na(spliceai))
    spliceai <- tbl %>% filter(spliceai >= min_SpliceAI) %>% filter(spliceai <= max_SpliceAI)
    spliceai <- rbind(spliceai, na_values_spliceai)
    tbl <- intersect(spliceai, tbl)
    
    #Impact Type
    impact_type <- input$impact_type
    impact <- c()
    for (i in 1:length(impact_type)){
      a <- tbl %>% filter(grepl(impact_type[i], IMPACT))
      impact <- rbind(impact, a)
      i = i+1
    }
    tbl <- impact
    return(tbl)
  })
  
  
  
  output$table_output_comp_het <- renderDataTable({
    
    tryCatch({
    as.datatable(formattable(mydata_comp_het(), list(
      SIFT = formatter("span",
                       style = x ~ style(
                         display = "block",
                         "border-radius" = "4px",
                         "padding-right" = "4px",
                         "background-color" = sapply(x,color.picker.SIFT)
                       )),
      PolyPhen = formatter("span",
                           style = x ~ style(
                             display = "block",
                             "border-radius" = "4px",
                             "padding-right" = "4px",
                             "background-color" = sapply(x,color.picker.PolyPhen))),
      IMPACT = formatter('span',
                         style = x ~ style(
                           display = "block",
                           "border-radius" = "4px",
                           "background-color" = sapply(x,color.picker.Type)))
    )), caption =htmltools::tags$caption(
      style = 'caption-side: bottom; text-align: left;', htmltools::em('SIFT: 0.0 to 0.05 (deleterious) to 0.05 to 1.0 (tolerated);', 'PolyPhen-2: 0.0 to 0.15 (tolerated) to 0.15 to 0.85 (possibly damaging) to 0.85 to 1.0 (damaging);', 'SpliceAI: 0.2 (high recall/likely pathogenic), 0.5 (recommended/pathogenic), and 0.8 (high precision/pathogenic) cutoffs;','o/e: low oe values are indicative of strong intolerance;','pLI >= 0.9 are an extremely intolerant set of transcripts')),
    extensions = c('Buttons','ColReorder','RowReorder','FixedColumns'),escape = FALSE,options = list(rowReorder = TRUE,colReorder = TRUE,paging  = FALSE,
                                                                                                                    fixedColumns = TRUE,autowidth=TRUE, searchHighlight=TRUE,dom = 'Bfrtip',scrollX = TRUE,
                                                                                                                    buttons = 
                                                                                                                      list(
                                                                                                                        extend = 'collection',
                                                                                                                        buttons = c('csv', 'excel', 'pdf'),
                                                                                                                        text = 'Download'),I('colvis'))
    )
    
    
    }, error=function(e){})
    
  })
  
  output$table_output_known_genes <- renderDataTable({
    
    tryCatch({
    as.datatable(formattable(mydata_associated_genes(), list(
      SIFT = formatter("span",
                       style = x ~ style(
                         display = "block",
                         "border-radius" = "4px",
                         "padding-right" = "4px",
                         "background-color" = sapply(x,color.picker.SIFT)
                       )),
      PolyPhen = formatter("span",
                           style = x ~ style(
                             display = "block",
                             "border-radius" = "4px",
                             "padding-right" = "4px",
                             "background-color" = sapply(x,color.picker.PolyPhen))),
      IMPACT = formatter('span',
                         style = x ~ style(
                           display = "block",
                           "border-radius" = "4px",
                           "background-color" = sapply(x,color.picker.Type)))
    )), caption =htmltools::tags$caption(
      style = 'caption-side: bottom; text-align: left;', htmltools::em('SIFT: 0.0 to 0.05 (deleterious) to 0.05 to 1.0 (tolerated);', 'PolyPhen-2: 0.0 to 0.15 (tolerated) to 0.15 to 0.85 (possibly damaging) to 0.85 to 1.0 (damaging);', 'SpliceAI: 0.2 (high recall/likely pathogenic), 0.5 (recommended/pathogenic), and 0.8 (high precision/pathogenic) cutoffs;','o/e: low oe values are indicative of strong intolerance;','pLI >= 0.9 are an extremely intolerant set of transcripts')),
    extensions = c('Buttons','ColReorder','RowReorder','FixedColumns'),escape = FALSE,options = list(rowReorder = TRUE,colReorder = TRUE,paging  = FALSE,
                                                                                                                    fixedColumns = TRUE,autowidth=TRUE, searchHighlight=TRUE,dom = 'Bfrtip',scrollX = TRUE,
                                                                                                                    buttons = 
                                                                                                                      list(
                                                                                                                        extend = 'collection',
                                                                                                                        buttons = c('csv', 'excel', 'pdf'),
                                                                                                                        text = 'Download'),I('colvis'))
    )
    
    }, error=function(e){})
    
    
  })
  
  output$table_based_on_input_genes <- renderDataTable({
  tryCatch({
    as.datatable(formattable(known_genes(), list(
      SIFT = formatter("span",
                       style = x ~ style(
                         display = "block",
                         "border-radius" = "4px",
                         "padding-right" = "4px",
                         "background-color" = sapply(x,color.picker.SIFT)
                       )),
      PolyPhen = formatter("span",
                           style = x ~ style(
                             display = "block",
                             "border-radius" = "4px",
                             "padding-right" = "4px",
                             "background-color" = sapply(x,color.picker.PolyPhen))),
      IMPACT = formatter('span',
                         style = x ~ style(
                           display = "block",
                           "border-radius" = "4px",
                           "background-color" = sapply(x,color.picker.Type)))
    )), caption =htmltools::tags$caption(
      style = 'caption-side: bottom; text-align: left;', htmltools::em('SIFT: 0.0 to 0.05 (deleterious) to 0.05 to 1.0 (tolerated);', 'PolyPhen-2: 0.0 to 0.15 (tolerated) to 0.15 to 0.85 (possibly damaging) to 0.85 to 1.0 (damaging);', 'SpliceAI: 0.2 (high recall/likely pathogenic), 0.5 (recommended/pathogenic), and 0.8 (high precision/pathogenic) cutoffs;','o/e: low oe values are indicative of strong intolerance;','pLI >= 0.9 are an extremely intolerant set of transcripts')),
    extensions = c('Buttons','ColReorder','RowReorder','FixedColumns'),escape = FALSE,options = list(rowReorder = TRUE,colReorder = TRUE,paging  = FALSE,
                                                                                                     fixedColumns = TRUE,autowidth=TRUE, searchHighlight=TRUE,dom = 'Bfrtip',scrollX = TRUE,
                                                                                                     buttons = 
                                                                                                       list(
                                                                                                         extend = 'collection',
                                                                                                         buttons = c('csv', 'excel', 'pdf'),
                                                                                                         text = 'Download'),I('colvis'))
    )
    
    
  }, error=function(e){})

  })
  
  
  
  
  
}
shinyApp(ui, server)
