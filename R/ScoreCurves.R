ScoreCurve <- function(data, colors, title, stats, colormatch, alt.heights){
  require(ggplot2)
  require(ggpubr)
  require(dplyr)
  require(tidyr)
  dft = function(df) {
    tdf = t(df)
    colnames(tdf) = rownames(df)
    rownames(tdf) = colnames(df)
    return(data.frame(tdf, check.names=FALSE, stringsAsFactors=FALSE))
  }
  EAEdata <- na.omit(data)
  colnames(EAEdata) <- gsub('X', '', colnames(EAEdata))
  EAEdata$TGS <- paste(EAEdata$Treatment, EAEdata$Sex)
  EAElong <- EAEdata %>% gather('Day', 'Score', -'Treatment', -"Sex", -"TGS")
  EAElong$Day <- as.numeric(EAElong$Day)
  EAElong$Score <- as.numeric(EAElong$Score)
  if(missing(stats)){
    resA <- compare_means(Score ~ TGS ,EAElong, method = 'kruskal.test', group.by = c("Day"), paired = F)
    resT <- compare_means(Score ~ TGS, EAElong, method = 'wilcox.test', group.by = c("Day"), paired = F)
  } else if (stats == 'parametric'){
    resA <- compare_means(Score ~ TGS ,EAElong, method = 'anova', group.by = c("Day"), paired = F)
    resT <- compare_means(Score ~ TGS, EAElong, method = 't.test', group.by = c("Day"), paired = F)
  } else {
    resA <- compare_means(Score ~ TGS ,EAElong, method = 'kruskal.test', group.by = c("Day"), paired = F)
    resT <- compare_means(Score ~ TGS, EAElong, method = 'wilcox.test', group.by = c("Day"), paired = F)
  }


  write.csv(resT, file = "Two-Group-Comparisons.csv")
  write.csv(resA, file = "All-Group-Comparisons.csv")

  Mean <- EAElong %>%
    group_by(Day, TGS) %>%
    summarise_all(funs(mean))

  MaxScores <- aggregate(Score ~ Day, data = Mean, max)
  resA <- merge(resA, MaxScores, by=c('Day', 'Day'))
  colnames(resA)[colnames(resA)=="Score"] <- "max"

  resA <- resA[resA$p.adj != 'NaN', ]
  if (missing(alt.heights)){
    resA$p.adj.star1 <- ifelse(resA$p < 0.05, '*', "")
    resA$p.adj.star1.height <- resA$max + 0.5
    resA$p.adj.star2 <- ifelse(resA$p < 0.01, '*', "")
    resA$p.adj.star2.height <- resA$max + 0.65
    resA$p.adj.star3 <- ifelse(resA$p < 0.001, '*', "")
    resA$p.adj.star3.height <- resA$max + 0.80

  } else {
    resA$p.adj.star1 <- ifelse(resA$p < 0.05, '*', "")
    resA$p.adj.star1.height <- resA$max + alt.heights[1]
    resA$p.adj.star2 <- ifelse(resA$p < 0.01, '*', "")
    resA$p.adj.star2.height <- resA$max +  alt.heights[2]
    resA$p.adj.star3 <- ifelse(resA$p < 0.001, '*', "")
    resA$p.adj.star3.height <- resA$max +  alt.heights[3]
  }

  resA$Day <- as.numeric(resA$Day)
  EAElong$Day <- as.numeric(EAElong$Day)
  unique2 <- unique(EAElong$TGS)

  ggob = ggline(EAElong,
                y = "Score",
                x = "Day", group = "TGS", add = "mean_se", width = 5,
                color = "TGS")
  #sel = (0.01 < resA$p.adj & resA$p.adj < 0.05)
  #ggob = ggob + geom_signif(data=anno_df, aes(xmin = group1, xmax = group2, annotations = p.adj, y_position = y_pos), manual= TRUE)
  if(missing(colormatch)){
    ggob = ggob + color_palette(colors)
  } else if (missing(colors)){
    ggob = ggob + scale_color_manual(
      values = colormatch)
  }
  ggob = ggob + ylab('Disease Score')
  ggob = ggob + ggtitle(title)
  # ggob = ggob + labs(fill = "Treatment Group and Sex")
  ggob = ggob + annotate('text', x= resA$Day - min(EAElong$Day) + 1, y=resA$p.adj.star1.height, label=resA$p.adj.star1, size=6)
  ggob = ggob + annotate('text', x= resA$Day - min(EAElong$Day) + 1, y=resA$p.adj.star2.height, label=resA$p.adj.star2, size=6)
  ggob = ggob + annotate('text', x= resA$Day - min(EAElong$Day) + 1, y=resA$p.adj.star3.height, label=resA$p.adj.star3, size=6)
  # ggob = ggob + annotate('text', x= resA$Day, y=resA$p.adj.star1.height, label=resA$p.adj.star1, size=6)
  # ggob = ggob + annotate('text', x= resA$Day, y=resA$p.adj.star2.height, label=resA$p.adj.star2, size=6)
  # ggob = ggob + annotate('text', x= resA$Day, y=resA$p.adj.star3.height, label=resA$p.adj.star3, size=6)
  ggob = ggob + theme(
    axis.text.x.bottom = element_text(size=8)
  )
  #pdf(paste0(title, ".pdf"), width=8, height=4, res = 300)
  print(ggob)
  ggsave(file= "ScoreCurve.eps", plot = last_plot(), h=4, w=8, dpi=320, units = c('in'), device = "eps")
  dev.off()
  ##########################################################################################################################
  unique <- split(resT, list(resT$group1, resT$group2))
  unique <- unique[sapply(unique, function(x) dim(x)[1]) > 0]
  i=1
  while (i<1+(length(unique))){
    loaded_df <- as.data.frame(unique[i])
    fullname <- colnames(loaded_df)[2]
    name <- gsub('..y.', '', fullname)
    name <- gsub('\\.', '_', name)

    group1 <- unique(loaded_df[3])
    group1 <- as.list(group1[1])
    group2 <- unique(loaded_df[4])
    group2 <- as.list(group2[1])

    Mean1 <- Mean[Mean$TGS == group1 | Mean$TGS == group2,]

    MaxScores1 <- aggregate(Score ~ Day, data = Mean1, max)
    colnames(loaded_df)[1] <- 'Day'
    colnames(loaded_df)[5] <- 'p'
    resTloaded <- merge(loaded_df, MaxScores1, by.x = 'Day', by.y = 'Day')
    if (missing(alt.heights)){
      resTloaded$p.adj.star1 <- ifelse(resTloaded$p < 0.05, '*', "")
      resTloaded$p.adj.star1.height <- resTloaded$Score + 0.45
      resTloaded$p.adj.star2 <- ifelse(resTloaded$p < 0.01, '*', "")
      resTloaded$p.adj.star2.height <- resTloaded$Score + 0.6
      resTloaded$p.adj.star3 <- ifelse(resTloaded$p < 0.001, '*', "")
      resTloaded$p.adj.star3.height <- resTloaded$Score + 0.75
    } else {
      resA$p.adj.star1 <- ifelse(resA$p < 0.05, '*', "")
      resA$p.adj.star1.height <- resA$max + alt.heights[1]
      resA$p.adj.star2 <- ifelse(resA$p < 0.01, '*', "")
      resA$p.adj.star2.height <- resA$max +  alt.heights[2]
      resA$p.adj.star3 <- ifelse(resA$p < 0.001, '*', "")
      resA$p.adj.star3.height <- resA$max +  alt.heights[3]
    }



    EAElongloaded <- EAElong[EAElong$TGS == group1 | EAElong$TGS == group2,]

    ggob = ggline(EAElongloaded,
                  y = "Score",
                  x = "Day", group = "TGS", add = "mean_se", width = 5,
                  color = "TGS")
    #ggob = ggob + geom_signif(data=anno_df, aes(xmin = group1, xmax = group2, annotations = p.adj, y_position = y_pos), manual= TRUE)
    if(missing(colormatch)){
      ggob = ggob + color_palette(colors)
    } else if (missing(colors)){
      ggob = ggob + scale_color_manual(
        values = colormatch)
    }
    ggob = ggob + ylab('Disease Score')
    ggob = ggob + ggtitle(title)
    # ggob = ggob + labs(fill = "Treatment Group and Sex")
    ggob = ggob + annotate('text', x= resTloaded$Day - min(EAElong$Day) + 1, y=resTloaded$p.adj.star1.height, label=resTloaded$p.adj.star1, size=6)
    ggob = ggob + annotate('text', x= resTloaded$Day - min(EAElong$Day) + 1, y=resTloaded$p.adj.star2.height, label=resTloaded$p.adj.star2, size=6)
    ggob = ggob + annotate('text', x= resTloaded$Day - min(EAElong$Day) + 1, y=resTloaded$p.adj.star3.height, label=resTloaded$p.adj.star3, size=6)
    # ggob = ggob + annotate('text', x= resA$Day, y=resA$p.adj.star1.height, label=resA$p.adj.star1, size=6)
    # ggob = ggob + annotate('text', x= resA$Day, y=resA$p.adj.star2.height, label=resA$p.adj.star2, size=6)
    # ggob = ggob + annotate('text', x= resA$Day, y=resA$p.adj.star3.height, label=resA$p.adj.star3, size=6)
    ggob = ggob + theme(
      axis.text.x.bottom = element_text(size=8))
    print(ggob)
    ggsave(paste0(name, "ScoreCurve.eps"), plot = last_plot(), h=4, w=8, dpi=320, units = c('in'), device = "eps")
    dev.off()
    i = i+1



  }


}
