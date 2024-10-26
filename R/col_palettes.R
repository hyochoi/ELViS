
generate_palettes<-FALSE
if(!generate_palettes){
  choi_discrete_palettes <-
    list(
      Rushmore1 = c("#E1BD6D", "#EABE94", "#0B775E", "#35274A", "#F2300F")
      ,Royal1 = c("#899DA4", "#C93312", "#FAEFD1", "#DC863B")
      ,Royal2 = c("#9A8822", "#F5CDB4", "#F8AFA8", "#FDDDA0", "#74A089")
      ,Zissou1 = c("#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00")
      ,Darjeeling1 = c("#FF0000", "#00A08A", "#F2AD00", "#F98400", "#5BBCD6")
      ,Darjeeling2 = c("#ECCBAE", "#046C9A", "#D69C4E", "#ABDDDE", "#000000")
      ,Chevalier1 = c("#446455", "#FDD262", "#D3DDDC", "#C7B19C")
      ,FantasticFox1 = c("#DD8D29", "#E2D200", "#46ACC8", "#E58601", "#B40F20")
      ,Moonrise2 = c("#798E87", "#C27D38", "#CCC591", "#29211F")
      ,Moonrise3 = c("#85D4E3", "#F4B5BD", "#9C964A", "#CDC08C", "#FAD77B")
      ,Cavalcanti1 = c("#D8B70A", "#02401B", "#A2A475", "#81A88D", "#972D15")
      ,GrandBudapest1 = c("#F1BB7B", "#FD6467", "#5B1A18", "#D67236")
      ,GrandBudapest2 = c("#E6A0C4", "#C6CDF7", "#D8A499", "#7294D4")
      ,aurora = c("#BF616A", "#D08770", "#EBCB8B", "#A3BE8C", "#B48EAD")
      ,Accent = c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17", "#666666")
      ,Dark2 = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")
      ,Pastel1 = c("#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC", "#F2F2F2")
      ,Pastel2 = c("#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9", "#FFF2AE", "#F1E2CC", "#CCCCCC")
      ,Set1 = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")
      ,Set2 = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3")
      ,Set3 = c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")
    )
  choi_paired_palettes <-
    list(
      blue = c("#A6CEE3", "#1F78B4")
      ,green = c("#B2DF8A", "#33A02C")
      ,red = c("#FB9A99", "#E31A1C")
      ,orange = c("#FDBF6F", "#FF7F00")
      ,purple = c("#CAB2D6", "#6A3D9A")
      ,brown = c("#FFFF99", "#B15928")
    )

  # choi_ordered_palettes %>% lapply(paste,collapse='", "') %>% lapply(\(x) paste0('"',x,'"')) %>% {tmp=.;names(tmp) %>% sapply(\(x) paste0(",",x," = c(",tmp[[x]],")\n"))} %>% cat
  choi_ordered_palettes <-
    list(
      lumina = c("#EDDAEB", "#AD8CAE", "#4F93B8", "#306489", "#222B4C")
      ,silver_mine = c("#4B644B", "#647D4B", "#E1E1E1", "#7D96AF", "#647D96")
      ,lake_superior = c("#7D4B19", "#C89664", "#C87d4B", "#4B647D", "#324B64", "#19324B")
      ,victory_bonds = c("#AF1900", "#C83200", "#E19600", "#193264", "#001964")
      ,halifax_harbor = c("#E1C8AF", "#C8AF96", "#AF967D", "#967D7D", "#644B64", "#4B324b")
      ,moose_pond = c("#4B3232", "#7D4B32", "#966432", "#AF7D32", "#E19632", "#E1AF4B", "#C8C896", "#4B4B4B")
      ,algoma_forest = c("#4B4B4B", "#967D4B", "#AFAF7D", "#C89632", "#647D64", "#96AFAF", "#7D96AF")
      ,red_mountain = c("#7D3232", "#7D4B4B", "#7D6464", "#AF967D", "#FAC87D", "#E1AF64", "#C8964B", "#32324B")
      ,afternoon_prarie = c("#486090", "#6078A8", "#7890A8", "#90A8C0", "#F0D8C0", "#D6BBCF", "#A8C0C0", "#C0D8D8", "#A8A890")
    )

  choi_continuous_palettes <-
    list(
      RdYlBu = c("#A50026", "#FFFFBF", "#313695")
      ,RdGy = c("#67001F", "#FFFFFF", "#1A1A1A")
      ,PRGn = c("#40004B", "#F7F7F7", "#00441B")
      ,PiYG = c("#8E0152", "#F7F7F7", "#276419")
      ,BrBG = c("#543005", "#F5F5F5", "#003C30")
      ,YlGB = c("yellow", "grey", "black")
      ,BlGn = c("#1a1334", "#26294a", "#01545a", "#017351", "#03c383", "#aad962")
      ,PuYl = c("#110141", "#710162", "#a12a5e", "#ed0345", "#ef6a32", "#fbbf45")
      ,BuGn = c("#3e71a8", "#577f9f", "#698e96", "#779d8d", "#84ad83", "#8fbd77", "#99cd6b", "#a2dd5c", "#aaee49", "#b2ff2e")
      ,rainbow = c("#D12600", "#DB6A00", "#B2FF2E", "#00AD00", "#9CCADE", "#005B94", "#1E2085", "#610052", "#953272")
    )



  candicol1 <-
    c("#FFFFCC", "#FFF2AE", "#FFFFE5", "#FFFFCC", "#FFEDA0", "#FFF5F0", "#FFF7F3", "#FFF7EC", "#FEE8C8", "#FFF5EB", "#FEE6CE", "aliceblue"
    ) # candidate colors for exonic regions

  candicol2 <-
    c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F"
    ) # candidiate colors for regions with shape changes

  candicol3 <-
    c("#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9", "#FFF2AE", "#F1E2CC", "#CCCCCC"
    ) # candidiate colors for regions with shape changes

  candicol <- c(candicol2,candicol3);
  exon.col <- candicol1[9]

  #  piratepal(palette = "info2") %>% paste(collapse='", "') %>% paste0('"',.,'"')  %>% cat
  col_yarrr_info2  <- c("#006A40FF", "#F08892FF", "#75B41EFF", "#95828DFF", "#708C98FF", "#8AB8CFFF", "#007E7FFF", "#358359FF", "#8BA1BCFF", "#5A5895FF", "#F2990CFF", "#5A5895FF", "#E5BA3AFF", "#D86C4FFF")

}else{
  pkgs_undetected <- c()
  for(pkg in c("nord","wesanderson","RColorBrewer")){
    if (!requireNamespace(pkg, quietly = TRUE)) {
      pkgs_undetected <- c(pkgs_undetected,pkg)
    }
  }

  if(length(pkgs_undetected)>0){
    pkgs_undetected_str <- paste0("'",paste(pkgs_undetected,collapse="', '"),"'")
    stop(glue("The '{pkgs_undetected_str}' package is required but not installed. Please install it."))
  }


  if(!requireNamespace("nord"))
    choi_discrete_palettes <- c(wes_palettes[c("Rushmore1","Royal1","Royal2","Zissou1",
                                               "Darjeeling1","Darjeeling2",
                                               "Chevalier1","FantasticFox1","Moonrise2",
                                               "Moonrise3","Cavalcanti1",
                                               "GrandBudapest1","GrandBudapest2")],
                                nord_palettes["aurora"])
  choi_discrete_palettes$Accent <- RColorBrewer::brewer.pal(8,"Accent")
  choi_discrete_palettes$Dark2 <- RColorBrewer::brewer.pal(8,"Dark2")
  choi_discrete_palettes$Pastel1 <- RColorBrewer::brewer.pal(9,"Pastel1")
  choi_discrete_palettes$Pastel2 <- RColorBrewer::brewer.pal(8,"Pastel2")
  choi_discrete_palettes$Set1 <- RColorBrewer::brewer.pal(9,"Set1")
  choi_discrete_palettes$Set2 <- RColorBrewer::brewer.pal(8,"Set2")
  choi_discrete_palettes$Set3 <- RColorBrewer::brewer.pal(12,"Set3")

  choi_paired_palettes <- as.list(NULL)
  choi_paired_palettes$blue <- RColorBrewer::brewer.pal(12,"Paired")[c(1,2)]
  choi_paired_palettes$green <- RColorBrewer::brewer.pal(12,"Paired")[c(3,4)]
  choi_paired_palettes$red <- RColorBrewer::brewer.pal(12,"Paired")[c(5,6)]
  choi_paired_palettes$orange <- RColorBrewer::brewer.pal(12,"Paired")[c(7,8)]
  choi_paired_palettes$purple <- RColorBrewer::brewer.pal(12,"Paired")[c(9,10)]
  choi_paired_palettes$brown <- RColorBrewer::brewer.pal(12,"Paired")[c(11,12)]

  choi_ordered_palettes <- nord_palettes[c(5,7,8,9,10,11,12,14,16)]
  # lumina, silver_mine, lake_superior, victory_bonds, halifax_harbor,
  # moose_pond, algoma_forest, red_mountain, afternnon_prairie
  choi_continuous_palettes <- as.list(NULL)
  choi_continuous_palettes$RdYlBu <- RColorBrewer::brewer.pal(11,"RdYlBu")[c(1,6,11)]
  choi_continuous_palettes$RdGy <- RColorBrewer::brewer.pal(11,"RdGy")[c(1,6,11)]
  choi_continuous_palettes$PRGn <- RColorBrewer::brewer.pal(11,"PRGn")[c(1,6,11)]
  choi_continuous_palettes$PiYG <- RColorBrewer::brewer.pal(11,"PiYG")[c(1,6,11)]
  choi_continuous_palettes$BrBG <- RColorBrewer::brewer.pal(11,"BrBG")[c(1,6,11)]
  choi_continuous_palettes$YlGB <- c("yellow","grey","black")
  choi_continuous_palettes$BlGn <- c("#1a1334", "#26294a", "#01545a", "#017351", "#03c383",
                                     "#aad962")
  choi_continuous_palettes$PuYl <- c("#110141", "#710162", "#a12a5e", "#ed0345", "#ef6a32",
                                     "#fbbf45")
  choi_continuous_palettes$BuGn <- c("#3e71a8", "#577f9f", "#698e96", "#779d8d", "#84ad83",
                                     "#8fbd77", "#99cd6b", "#a2dd5c", "#aaee49", "#b2ff2e")
  choi_continuous_palettes$rainbow <- c(rosso_corsa = "#D12600", spanish_orange = "#DB6A00",
                                        green_yellow = "#B2FF2E", green = "#00AD00", pale_cerulean = "#9CCADE",
                                        sea_blue = "#005B94", st_patricks_blue = "#1E2085", tyrian_purple = "#610052",
                                        amaranth_deep_purple = "#953272")






  candicol1 <- c(RColorBrewer::brewer.pal(9,"Pastel1")[6], # candidate colors for exonic regions
                 RColorBrewer::brewer.pal(8,"Pastel2")[6],
                 RColorBrewer::brewer.pal(9,"YlOrBr")[1],
                 RColorBrewer::brewer.pal(9,"YlOrRd")[1],
                 RColorBrewer::brewer.pal(9,"YlOrRd")[2],
                 RColorBrewer::brewer.pal(9,"Reds")[1],
                 RColorBrewer::brewer.pal(9,"RdPu")[1],
                 RColorBrewer::brewer.pal(9,"OrRd")[1],
                 RColorBrewer::brewer.pal(9,"OrRd")[2],
                 RColorBrewer::brewer.pal(9,"Oranges")[1],
                 RColorBrewer::brewer.pal(9,"Oranges")[2],
                 "aliceblue");
  candicol2 <- RColorBrewer::brewer.pal(12,"Set3") # candidiate colors for regions with shape changes
  candicol3 <- RColorBrewer::brewer.pal(8,"Pastel2") # candidiate colors for regions with shape changes
  candicol <- c(candicol2,candicol3);
  exon.col <- candicol1[9]


  col_yarrr_info2 <- piratepal(palette = "info2")



}
