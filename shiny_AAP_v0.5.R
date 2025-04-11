#---------------------------------------------------------------------
# 作者：Xu Runlai
# PS: 该shiny app仅简单打包先前使用的统计分析和作图代码，未进行整合和去冗余，请谨慎使用！如有意见和建议请不吝反馈。
#---------------------------------------------------------------------

library(shiny)
library(ggplot2)
library(RColorBrewer)
library(viridis)

# ---- load file
ui_upload <- sidebarLayout(
  sidebarPanel(
    fileInput(inputId = "file", label="Input your data", buttonLabel = "Upload...")),
  mainPanel(
    h3("Raw data"),
    DT::DTOutput("preview1")
  )
)
# # ----melt
# ui_melt <- sidebarLayout(
#   sidebarPanel(
#     checkboxInput("melt", "melt?")
#   ),
#   mainPanel(
#     h3("Melt data"),
#     DT::DTOutput("preview4")
# ))
# ---- choose hsd method & download
ui_sig <- sidebarLayout(
  sidebarPanel(
    selectInput(inputId = "state", "Choose hsd methods", c('duncan','lsd', 'tukey')),
    downloadButton(outputId = "download", label = "Download ANOVA results!", class = "btn-block")
  ),
  mainPanel(
    h3("One-way ANOVA results"),
    DT::DTOutput("preview2")
  )
)


# ----two way anova
# ui_twoway <- sidebarLayout(
#   selectInput(inputId = "fac1", "Choose factor1", c(colnames(raw()))),
#   selectInput(inputId = "fac2", "Choose factor2", c(colnames(raw()))),
#   mainPanel(
#     h3("Two-way ANOVA results"),
#     DT::DTOutput("preview3")
#   )
# )
# ui_twoway <- fluidPage(
#   h3("Two-way ANOVA results"),
#   DT::DTOutput("preview3")
#   )

# ----change plot and group type
ui_plot1 <- fluidPage(
  titlePanel("Quick visualization"),
  fluidRow(
    column(6,selectInput("dataset", "Plot type", choices = c('box','violin','bar','halves',
                                                             'line','line_noer','heatmap',
                                                             'bar_stack','bar_stackp'))),
    column(6,selectInput("gtype", "Group type", choices = c('none','facet','group'))),
    fluidRow(column(6),column(6))
  ))
ui_plot2 <- fluidPage(
  fluidRow(
    plotOutput("plot")),
  fluidRow(column(6),column(6))
)

# ---- change x/y axis, siglab position, and output width and height
ui_axis <- fluidPage(
  fluidRow(
    column(6,textInput("x", label="x-axis",value = "x-axis")),
    column(6,textInput("y", label="y-axis",value = "y-axis"))),
  fluidRow(
    column(6,textInput("palette", label="Palette",value = "")),
    column(6,selectInput("theme", "Theme", choices = c('classic','bw','bw_blank')))
  ),
  # selectInput("palettetype", "Palette type", choices = c('color','fill')),
  fluidRow(
    column(6,selectInput("er", "Error bar", choices = c('se','sd'))),
    column(6,selectInput("jitter", "Jitter", choices = c('no','yes')))
  ),
  fluidRow(
    column(6,selectInput("sigmethod", "Sigmethod", choices = c('none','anova'))),
    column(6,sliderInput("siglab", "Sig label position", min = 0, max = 1, value = 1))
  ),
  fluidRow(
    column(6,sliderInput("width", "Output width", min = 4, max = 20, value = 5)),
    column(6,sliderInput("height", "Output height", min = 4, max = 20, value = 5)
    )))

# ----
ui_downloadp <- fluidRow(
  column(width = 12, downloadButton(outputId = "downloadp", label = "Download plot!", class = "btn-block")))

# ----
ui <- fluidPage(
  ui_upload,
  # ui_melt,
  ui_sig,
  # ui_twoway,
  ui_plot1,ui_plot2,ui_axis,
  ui_downloadp
)

# ----
server <- function(input, output, session) {
  # Upload ---------------------------------------------------------
  raw <- reactive({
    req(input$file)
    vroom::vroom(input$file$datapath)
  })
  output$preview1 <- DT::renderDT(raw(), options = list(pageLength = 5))
  # melt-----------------------------------------------------------
  melt <- reactive({
    out <- raw()
    if (length(colnames(out)) > 2 & !'facet' %in% colnames(out)) {
      out <- reshape2::melt(out, id.vars = c(colnames(out)[colnames(out)=='group']),
                            measure.vars = c(colnames(out)[colnames(out)!='group']),
                            variable.name = c('facet'),
                            value.name = 'value')
    }
    out
  })
  # output$preview4 <- DT::renderDT(melt(),options = list(pageLength = 5))
  # sig -----------------------------------------------------------
  tidied <- reactive({
    DVana <- function(data, method = 'duncan'){
      library(agricolae)
      if (method %in% c('duncan', 'lsd', 'tukey')) {
        if ('facet' %in% colnames(data)) {
          facetgroup <- data.frame(unique(data$facet))
          colnames(facetgroup) <- "facetgroup"
          
          resulte <-  data.frame(group="",mean="",sd="",se="",label="",max="",facet = "")[-1,]
          for (i in c(1:length(unique(data$facet)))){
            facetgroup1 <- facetgroup[i,]
            datafacet <- data[data$facet == facetgroup1,]
            oneway <- aov(value~group, data = datafacet)
            anova(oneway)
            if (method == 'duncan') {
              result <- duncan.test(oneway, 'group')
            } else if (method == 'lsd') {
              result <- LSD.test(oneway, 'group')
            } else if (method == 'tukey') {
              result <- HSD.test(oneway, 'group')
            }
            mar<-result$groups
            rownamemar<-row.names(mar)
            newmar<-data.frame(rownamemar,mar[,1],mar$groups)
            sort<-newmar[order(newmar$rownamemar),]
            rowname<-row.names(result$means)
            mean<-result$means[,1]
            sd<-result$means[,2]
            se <- sd/sqrt(result$means[,3])
            max <- result$means[,ifelse(method == 'lsd',8,6)]
            marker<-sort$mar.groups
            facet <- rep(facetgroup1,length(rowname))
            result2<-data.frame(rowname,mean,sd,se,marker,max,facet)
            colnames(result2) <- c('group','mean','sd','se','label','max','facet')
            resulte <- rbind(resulte,result2)
          }
          return(resulte)
        } else{
          oneway <- aov(value ~ group, data = data)
          anova(oneway)
          if (method == 'duncan') {
            result <- duncan.test(oneway, 'group')
          } else if (method == 'lsd') {
            result <- LSD.test(oneway, 'group')
          } else if (method == 'tukey') {
            result <- HSD.test(oneway, 'group')
          }
          mar<-result$groups
          rownamemar<-row.names(mar)
          newmar<-data.frame(rownamemar,mar[,1],mar$groups)
          sort<-newmar[order(newmar$rownamemar),]
          rowname<-row.names(result$means)
          mean<-result$means[,1]
          sd<-result$means[,2]
          se <- sd/sqrt(result$means[,3])
          max <- result$means[,ifelse(method == 'lsd',8,6)]
          marker<-sort$mar.groups
          result2<-data.frame(rowname,mean,sd,se,marker,max)
          colnames(result2) <- c('group','mean','sd','se','label','max')
          return(result2)
        }
      } else if(method %in% c('t.test', 'wilcox.test', "kruskal.test")) {
        result_purb <- compare_means(value ~ group, data = data, method = method)
        return(result_purb)
      }
    }
    out1 <- DVana(melt(), method = input$state)
  })
  output$preview2 <- DT::renderDT(tidied(), options = list(pageLength = 5))
  
  # two-way anova--------------------------------------------------
  # twoway <- reactive({
  #   data <- melt()
  #   data2 <- data[!colnames(data) %in% c('group','value')]
  #   twoways <- aov(data$value~data2[,1]+data2[,2])
  #   summary(twoways)
  # })
  # output$preview3 <- DT::renderDT(twoway(), options = list(pageLength = 5))
  
  # Download -------------------------------------------------------
  output$download <- downloadHandler(
    filename = function() {
      paste0(tools::file_path_sans_ext(input$file$name),'_anova_result', ".csv") 
    },
    content = function(file) {
      # vroom::vroom_write(tidied(), file) # .TSV则用
      write.csv(tidied(), file,row.names = F)
    }
  )
  # Plot -------------------------------------------------------
  pldata <- reactive({
    DVisual <- function(data, 
                        type = 'box', group_type = 'none', theme = 'bw', sigmethod = F, sigdata, er='se',
                        x_axis='x_axis',y_axis='y_axis'){
      palette <- unlist(ifelse(input$palette != "", strsplit(input$palette, split=","), input$palette))

      oa <- ifelse(input$jitter=='yes',0,1)
      ## facet and none----
      if (group_type %in% c('none','facet')){
        if (type == 'box') {
          p1 <- ggplot(data = data,aes(x=group,y = value))+
            geom_boxplot(aes(fill=group),color="black",width=0.5,size = 0.7,alpha= 1,position = position_dodge(1),notch = F,
                         outlier.colour = "black", outlier.size=1.5, outlier.alpha = oa)+
            stat_summary(fun.y='mean',geom='point',shape=4,size=3,color='black')
          # scale_y_continuous(limits = c(min(na.omit(data$value))*0.5, max(na.omit(data$value))*1.5))
        } else if (type == 'violin') {
          p1 <- ggplot(data = data,mapping = aes(x=group,y = value))+
            geom_violin(mapping = aes(fill = group),alpha=0.5,width=0.8,size=0.6,position=position_dodge(width=0.8), trim = T)+
            geom_boxplot(mapping = aes(fill = group),alpha= 1, width=0.1,size = 0.6,position = position_dodge(width=0.8),
                         outlier.size=1.5, outlier.alpha = oa, outlier.colour = "black")+
            stat_summary(fun.y='mean',geom='point',shape=4,size=3,color='black')
          # scale_y_continuous(limits = c(2000,9000))+
          # facet_wrap(~type,scales = "free_y",2,2)+
        } else if (type == 'halves') {
          library(gghalves)
          library(ggdist)
          p1 <- ggplot(data = data, mapping = aes(x = group, y = value, color = group))+
            geom_half_point(mapping = aes(color = group), side = "l",
                            range_scale =0.5, alpha =0.5)+
            geom_boxplot(aes(color=group),alpha= 0.5, outlier.size=1,outlier.alpha = oa,outlier.color = "black",
                         width=0.1,size = 0.8,position = position_dodge(1))+
            stat_halfeye(mapping = aes(fill=group),width = 0.25, justification = -0.65, .width = 0,point_colour = NA)+ # 中位数
            stat_summary(fun.y='mean',geom='point',shape=4,size=3,color="black") # 标记均值
          # scale_y_continuous(expand = c(0.05,0.05)
          #                    # ,limits = c(22,35)
          # )+
          # facet_grid(sample~.)
        } else if (type == 'bar') {
          # p1 <- ggplot(data = data,aes(x=group,y = value))+
          #   stat_summary(color = 'black', geom='errorbar',width=0.2, size = 0.5,
          #                fun.min=function(x){mean(x)-sd(x)},fun.max = function(x){mean(x)+sd(x)})+
          #   stat_summary(aes(fill=group), geom='bar',fun=mean,width = 0.6, alpha= 1,size = 0.5)+
          #   scale_y_continuous(expand=c(0,0),limits = c(0,max(na.omit(data$value))*1.1))
          ifelse(er == 'se',
                 p1 <- ggplot(sigdata, aes(x = group, y = mean, fill = group))+
                   geom_errorbar(aes(x=group,ymax = mean+se, ymin = mean-se),# errorbar:sd or se
                                 position = position_dodge(width = 0.8),width = 0.4,size = 0.5)+
                   geom_bar(stat = "identity",position = "dodge",color="transparent", width = 0.8, alpha= 1,size = 0.5)+
                   # scale_fill_manual(values = palette)+
                   # scale_color_manual(values = palette)+
                   scale_y_continuous(expand=c(0,0),limits = c(0,max(na.omit(sigdata$max))*1.1)),
                 p1 <- ggplot(sigdata, aes(x = group, y = mean, fill = group))+
                   geom_errorbar(aes(x=group,ymax = mean+sd, ymin = mean-sd),# errorbar:sd or se
                                 position = position_dodge(width = 0.8),width = 0.4,size = 0.5)+
                   geom_bar(stat = "identity",position = "dodge",color="transparent", width = 0.8, alpha= 1,size = 0.5)+
                   # scale_fill_manual(values = palette)+
                   # scale_color_manual(values = palette)+
                   scale_y_continuous(expand=c(0,0),limits = c(0,max(na.omit(sigdata$max))*1.1))
          )
        } else if (type == 'line'){
          ifelse(er == 'se',
                 p1 <- ggplot(sigdata, aes(x = group, y = mean, color = group, group = group))+
                   geom_errorbar(aes(x=group,ymax = mean+se, ymin = mean-se),# errorbar:sd or se
                                 width = 0.25,size = 0.5)+
                   geom_point(aes(x=group,color=group),size=2)+
                   geom_line(aes(x=group,colour = group))+
                   # scale_fill_manual(values = palette)+
                   # scale_color_manual(values = palette)+
                   scale_y_continuous(expand=c(0,0),limits = c(0,max(na.omit(sigdata$max))*1.1)),
                 p1 <- ggplot(sigdata, aes(x = group, y = mean, color = group, group = group))+
                   geom_errorbar(aes(x=group,ymax = mean+sd, ymin = mean-sd),# errorbar:sd or se
                                 width = 0.25,size = 0.5)+
                   geom_point(aes(x=group,color=group),size=2)+
                   geom_line(aes(x=group,colour = group))+
                   # scale_fill_manual(values = palette)+
                   # scale_color_manual(values = palette)+
                   scale_y_continuous(expand=c(0,0),limits = c(0,max(na.omit(sigdata$max))*1.1))
                 
          )
        } else if (type == 'line_noer'){
          p1 <- ggplot(sigdata, aes(x = group, y = mean, color = group, group = group))+
            geom_point(aes(x=group,color=group),size=2)+
            geom_line(aes(x=group,colour = group))+
            scale_y_continuous(expand=c(0,0),limits = c(0,max(na.omit(sigdata$max))*1.1))
        } else if (type == 'heatmap'){
          p1 <- ggplot(sigdata, aes(x = group, y = facet, fill = mean)) +
            geom_tile() +
            geom_text(mapping = aes(label=paste(mean,'',label)))+
            scale_fill_viridis(discrete=F, option=input$palette)+
            # scale_fill_gradient(low = "red", high = "yellow", na.value = NA) +
            # theme(axis.text.x = element_text(angle = 45, hjust = 1))
            labs(x=x_axis, y=y_axis)+
            theme_minimal()+
            theme(legend.position="right")+
            theme(text = element_text(color = "black",size=16))
          return(p1)
        }
        # facet
        if (group_type == 'none') {
          p2 <- p1
        } else if (group_type == 'facet') {
          p2 <- p1+
            facet_wrap(~facet,scales = "free_y")
        }
        # jitter
        if (input$jitter == 'yes') {
          p3 <- p2+
            geom_jitter(aes(color=group),position=position_jitter(0.2),size=2, alpha=0.5)
        } else if (input$jitter == 'no') {
          p3 <- p2
        }
        # theme
        if (input$theme == 'bw_blank') {
          ifelse(palette != '',
                 p4 <- p3+
                   scale_fill_manual(values = c(palette))+
                   scale_color_manual(values = c(palette))+
                   labs(x=x_axis, y=y_axis)+
                   theme_bw()+
                   theme(panel.grid=element_blank())+
                   theme(legend.position="right")+
                   theme(text = element_text(color = "black",size=16)),
                 p4 <- p3+
                   # scale_fill_manual(values = palette)+
                   # scale_color_manual(values = palette)+
                   labs(x=x_axis, y=y_axis)+
                   theme_bw()+
                   theme(panel.grid=element_blank())+
                   theme(legend.position="right")+
                   theme(text = element_text(color = "black",size=16)))
        } else if (input$theme == 'classic') {
          ifelse(palette != '',
                 p4 <- p3+
                   scale_fill_manual(values = c(palette))+
                   scale_color_manual(values = c(palette))+
                   labs(x=x_axis, y=y_axis)+
                   theme_classic()+
                   theme(legend.position="right")+
                   theme(text = element_text(color = "black",size=16)),
                 p4 <- p3+
                   # scale_fill_manual(values = palette)+
                   # scale_color_manual(values = palette)+
                   labs(x=x_axis, y=y_axis)+
                   theme_classic()+
                   theme(legend.position="right")+
                   theme(text = element_text(color = "black",size=16)))
          
        } else if (input$theme == 'bw') {
          ifelse(palette != '',
                 p4 <- p3+
                   scale_fill_manual(values = c(palette))+
                   scale_color_manual(values = c(palette))+
                   labs(x=x_axis, y=y_axis)+
                   theme_bw()+
                   # theme(panel.grid=element_blank())+
                   theme(legend.position="right")+
                   theme(text = element_text(color = "black",size=16)),
                 p4 <- p3+
                   # scale_fill_manual(values = palette)+
                   # scale_color_manual(values = palette)+
                   labs(x=x_axis, y=y_axis)+
                   theme_bw()+
                   # theme(panel.grid=element_blank())+
                   theme(legend.position="right")+
                   theme(text = element_text(color = "black",size=16)))
        }
        
        # significant
        if (sigmethod == 'none') {
          p5 <- p4
        } else if (sigmethod == 'anova'){
          p5 <- p4+
            geom_text(data=sigdata,aes(x = group, y = max*input$siglab, label = label), vjust = -0.5, color = 'black', size = 6)
        } else if (sigmethod %in% c('t.test', 'wilcox.test', "kruskal.test")) {
          group = levels(factor(data$group))
          data$group = factor(data$group, levels = group)
          comp = combn(group,2)
          my_comparisons=list()
          for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
          p5 <- p4+
            stat_compare_means(comparisons = my_comparisons, aes(label = after_stat(p.signif)),
                               method = sigmethod, size = 5, vjust = 1)
        }
        return(p5)
      }
      # group ---------------------------------
      if (group_type == 'group') {
        if (type == 'box') {
          p1 <- ggplot(data = data,aes(x=facet,y = value,fill=group))+
            geom_boxplot(aes(fill=group),width=0.5,size = 0.7,alpha= 1,notch = F,position = position_dodge(width=0.8),
                         outlier.colour = "black", outlier.size=1.5, outlier.alpha = oa)+
            stat_summary(fun.y='mean',geom='point',shape=4,size=3,color='black',position = position_dodge(0.8))
          # geom_point(data = sigdata,
          #            aes(x = facet, y = mean),
          #            position = position_dodge(0.8), shape=5,size=1.5,color='black')
          # scale_y_continuous(limits = c(min(na.omit(data$value))*0.5, max(na.omit(data$value))*1.5))
        } else if (type == 'violin') {
          p1 <- ggplot(data = data,mapping = aes(x=facet,y = value,fill=group))+
            geom_violin(mapping = aes(fill = group),alpha=0.5,width=0.8,size=0.6,position=position_dodge(width=0.8), trim = T)+
            geom_boxplot(mapping = aes(fill = group),alpha= 1, width=0.1,size = 0.6,position = position_dodge(width=0.8),
                         outlier.size=1.5, outlier.alpha = oa, outlier.colour = "black")+
            stat_summary(fun.y='mean',geom='point',shape=4,size=3,color='black',position = position_dodge(0.8))
          # scale_y_continuous(limits = c(2000,9000))+
          # facet_wrap(~type,scales = "free_y",2,2)+
        } else if (type == 'halves') {
          library(gghalves)
          library(ggdist)
          p1 <- ggplot(data = data, mapping = aes(x = facet,y = value, color = group))+
            geom_half_point(mapping = aes(color = group), side = "l",
                            range_scale =0.5, alpha =0.5)+
            geom_boxplot(aes(color=group),alpha= 0.5, outlier.size=1,outlier.alpha = oa,outlier.color = "black",
                         width=0.1,size = 0.8,position = position_dodge(1))+
            stat_halfeye(mapping = aes(fill=group),width = 0.25, justification = -0.65, .width = 0,point_colour = NA)+ # 中位数
            stat_summary(fun.y='mean',geom='point',shape=4,size=3,color="black",position = position_dodge(0.8)) # 标记均值
          # scale_y_continuous(expand = c(0.05,0.05)
          #                    # ,limits = c(22,35)
          # )+
          # facet_grid(sample~.)
        } else if (type == 'bar') {
          # p1 <- ggplot(data = data,aes(x=group,y = value))+
          #   stat_summary(color = 'black', geom='errorbar',width=0.2, size = 0.5,
          #                fun.min=function(x){mean(x)-sd(x)},fun.max = function(x){mean(x)+sd(x)})+
          #   stat_summary(aes(fill=group), geom='bar',fun=mean,width = 0.6, alpha= 1,size = 0.5)+
          #   scale_y_continuous(expand=c(0,0),limits = c(0,max(na.omit(data$value))*1.1))
          ifelse(er == 'se',
                 p1 <- ggplot(sigdata, aes(x = facet, y = mean, fill = group))+
                   geom_errorbar(aes(x=facet,ymax = mean+se, ymin = mean-se),# errorbar:sd or se
                                 position = position_dodge(width = 0.8),width = 0.4,size = 0.5)+
                   geom_bar(stat = "identity",position = "dodge",color="transparent", width = 0.8, alpha= 1,size = 0.5)+
                   # scale_fill_manual(values = palette)+
                   # scale_color_manual(values = palette)+
                   scale_y_continuous(expand=c(0,0),limits = c(0,max(na.omit(sigdata$max))*1.1)),
                 p1 <- ggplot(sigdata, aes(x = facet, y = mean, fill = group))+
                   geom_errorbar(aes(x=facet,ymax = mean+sd, ymin = mean-sd),# errorbar:sd or se
                                 position = position_dodge(width = 0.8),width = 0.4,size = 0.5)+
                   geom_bar(stat = "identity",position = "dodge",color="transparent", width = 0.8, alpha= 1,size = 0.5)+
                   # scale_fill_manual(values = palette)+
                   # scale_color_manual(values = palette)+
                   scale_y_continuous(expand=c(0,0),limits = c(0,max(na.omit(sigdata$max))*1.1))
          )
        } else if (type == 'line'){
          ifelse(er == 'se',
                 p1 <- ggplot(sigdata, aes(x = facet, y = mean, color = group, group = group))+
                   geom_errorbar(aes(x=facet,ymax = mean+se, ymin = mean-se),# errorbar:sd or se
                                 width = 0.25,size = 0.5)+
                   geom_point(aes(x=facet,color=group),size=2)+
                   geom_line(aes(x=facet,colour = group))+
                   scale_y_continuous(expand=c(0,0),limits = c(0,max(na.omit(sigdata$max))*1.1)),
                 p1 <- ggplot(sigdata, aes(x = facet, y = mean, color = group, group = group))+
                   geom_errorbar(aes(x=facet,ymax = mean+sd, ymin = mean-sd),# errorbar:sd or se
                                 width = 0.25,size = 0.5)+
                   geom_point(aes(x=facet,color=group),size=2)+
                   geom_line(aes(x=facet,colour = group))+
                   scale_y_continuous(expand=c(0,0),limits = c(0,max(na.omit(sigdata$max))*1.1))
          )
        } else if (type == 'line_noer'){
          p1 <- ggplot(sigdata, aes(x = facet, y = mean, color = group, group = group))+
            # geom_errorbar(aes(x=facet,ymax = mean+se, ymin = mean-se),# errorbar:sd or se
            #               width = 0.25,size = 0.5)+
            geom_point(aes(x=facet,color=group),size=2)+
            geom_line(aes(x=facet,colour = group))+
            scale_y_continuous(expand=c(0,0),limits = c(0,max(na.omit(sigdata$max))*1.1))
        } else if (type == 'bar_stackp'){
          p1 <- ggplot(sigdata, aes(x = facet, y = mean, fill = group))+
            geom_bar(stat = "identity",position = "fill",color="transparent", width = 0.8, alpha= 1,size = 0.5)+
            scale_y_continuous(expand=c(0,0))
        } else if (type == 'bar_stack'){
          p1 <- ggplot(sigdata, aes(x = facet, y = mean, fill = group))+
            geom_bar(stat = "identity",position = "stack",color="transparent", width = 0.8, alpha= 1,size = 0.5)+
            scale_y_continuous(expand=c(0,0))
        }
        
        p2 <- p1
        # jitter
        if (input$jitter == 'yes') {
          p3 <- p2+
            geom_jitter(aes(color=group),position = position_jitterdodge(0.3),size=2, alpha=0.5)
        } else if (input$jitter == 'no') {
          p3 <- p2
        }
        # theme
        if (input$theme == 'bw_blank') {
          ifelse(palette != '',
                 p4 <- p3+
                   scale_fill_manual(values = c(palette))+
                   scale_color_manual(values = c(palette))+
                   labs(x=x_axis, y=y_axis)+
                   theme_bw()+
                   theme(panel.grid=element_blank())+
                   theme(axis.text.x = element_text(angle =45, hjust = 0.5,vjust=0.5))+
                   theme(legend.position="right")+
                   theme(text = element_text(color = "black",size=16)),
                 p4 <- p3+
                   # scale_fill_manual(values = palette)+
                   # scale_color_manual(values = palette)+
                   labs(x=x_axis, y=y_axis)+
                   theme_bw()+
                   theme(panel.grid=element_blank())+
                   theme(axis.text.x = element_text(angle =45, hjust = 0.5,vjust=0.5))+
                   theme(legend.position="right")+
                   theme(text = element_text(color = "black",size=16)))
        } else if (input$theme == 'classic') {
          ifelse(palette != '',
                 p4 <- p3+
                   scale_fill_manual(values = c(palette))+
                   scale_color_manual(values = c(palette))+
                   labs(x=x_axis, y=y_axis)+
                   theme_classic()+
                   theme(axis.text.x = element_text(angle =45, hjust = 0.5,vjust=0.5))+
                   theme(legend.position="right")+
                   theme(text = element_text(color = "black",size=16)),
                 p4 <- p3+
                   # scale_fill_manual(values = palette)+
                   # scale_color_manual(values = palette)+
                   labs(x=x_axis, y=y_axis)+
                   theme_classic()+
                   theme(axis.text.x = element_text(angle =45, hjust = 0.5,vjust=0.5))+
                   theme(legend.position="right")+
                   theme(text = element_text(color = "black",size=16)))
          
        } else if (input$theme == 'bw') {
          ifelse(palette != '',
                 p4 <- p3+
                   scale_fill_manual(values = c(palette))+
                   scale_color_manual(values = c(palette))+
                   labs(x=x_axis, y=y_axis)+
                   theme_bw()+
                   theme(axis.text.x = element_text(angle =45, hjust = 0.5,vjust=0.5))+
                   # theme(panel.grid=element_blank())+
                   theme(legend.position="right")+
                   theme(text = element_text(color = "black",size=16)),
                 p4 <- p3+
                   # scale_fill_manual(values = palette)+
                   # scale_color_manual(values = palette)+
                   labs(x=x_axis, y=y_axis)+
                   theme_bw()+
                   theme(axis.text.x = element_text(angle =45, hjust = 0.5,vjust=0.5))+
                   # theme(panel.grid=element_blank())+
                   theme(legend.position="right")+
                   theme(text = element_text(color = "black",size=16)))
        }
        
        # significant
        if (sigmethod == 'none') {
          p5 <- p4
        } else if (sigmethod == 'anova'){
          p5 <- p4+
            geom_text(data=sigdata,aes(x = facet, y = max *input$siglab, label = label), vjust = -0.5, color = 'black', size = 6,
                      position = position_dodge(0.8))
        } else if (sigmethod %in% c('t.test', 'wilcox.test', "kruskal.test")) {
          group = levels(factor(data$group))
          data$group = factor(data$group, levels = group)
          comp = combn(group,2)
          my_comparisons=list()
          for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
          p5 <- p4+
            stat_compare_means(comparisons = my_comparisons, aes(label = after_stat(p.signif)),
                               method = sigmethod, size = 5, vjust = 1)
        }
        return(p5)
      }
    }
    pldata <- DVisual(data = melt(), type = input$dataset, group_type = input$gtype, theme = 'classic', 
                      sigmethod = input$sigmethod, sigdata = tidied(),x_axis=input$x,y_axis=input$y, er=input$er)
  })
  output$plot <- renderPlot({
    plot(pldata())
  }, res = 96) # margin
  
  # DOWNLOAD PLOT----
  output$downloadp <- downloadHandler(
    filename = function() {
      paste0(tools::file_path_sans_ext(input$file$name), ".pdf")
    },
    content = function(file) {
      # vroom::vroom_write(tidied(), file) # .TSV则用
      ggsave(file,pldata(), width = input$width, height = input$height, units = 'in',dpi = 600)
    }
  )
  
}
shinyApp(ui, server)