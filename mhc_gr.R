library(beanplot)

library(dplyr)

library(ggplot2)

library(RColorBrewer)

read_mhc <- function(folder)
{
  result = NULL
  models <- list.files(folder)
  for (m in models)
  {
    files <- list.files(paste0(folder, "/", m))
    for (f in files)
    {
      name <- strsplit(f, "\\.")[[1]][[1]]
      table <- read.csv(paste0(folder, "/", m, "/", f), sep = "\t", stringsAsFactors = FALSE)
      table$method <- m
      table$sample <- name
      if (is.null(result))
      {
        result <- table
      }
      else
       result <- rbind(result, table) 
    }
  }
  return(result)
}

fisher <- function(rand, data, threshold=100)
{
  r0 = sum(rand>=threshold)
  r1 = sum(rand<threshold)
  d0 = sum(data>=threshold)
  d1 = sum(data<threshold)
  
  ft = fisher.test(data.frame(c(r0,r1), c(d0,d1)), alternative = "g")
  return(ft$p.value)
}

filt <- function(vec, th = 10000)
{
  na.omit(ifelse(vec > th, yes = NA, no = vec))
}

filt1 <- function(vec, th = 20000)
{
  na.omit(ifelse(vec > th, yes = NA, no = vec))
}

##################################### I #############################################


indir <- "h2"
refset <- "randMB"
outdir <- paste0(indir, "_out")


mhci <- read_mhc(indir)
dir.create(outdir)

methods <- unique(mhci$method)
samples <- unique(mhci$sample)
alleles <- unique(mhci$model)
samples_all <- samples
samples <- samples[!startsWith(samples, "rand")]


options(digits=2)

#mhci$ic50 <- ifelse(mhci$ic50 > 5000, yes = NA, no = mhci$ic50)


report <- data.frame(row.names = samples)
ralleles <- data.frame(row.names = samples)
relalleles <- data.frame(row.names = samples)


for (m in methods)
{
  col = list()
  acol = list()
  
  for (s in samples)
  {
    print(paste(m, s))
    pdf(file = paste0(outdir, "/", s, "_", m, ".pdf"))
    sc <- sapply(alleles, function(x){
      mb <- na.omit(mhci[mhci$method == m & mhci$model == x & mhci$sample == "randMB", "ic50"])
      dt <- na.omit(mhci[mhci$method == m & mhci$model == x & mhci$sample == s, "ic50"])
      #wt <- wilcox.test(dt, mb, alternative = "less")
      ftpv <- fisher(mb,dt)
      return(ftpv)
    }) 
    
    alleles <- alleles[order(sc)]
    col[s] <- alleles[1]  
    acol <- paste(unique(na.omit(sapply(strsplit(mhci[mhci$method == m & mhci$ic50 < 50 & mhci$sample == s, "model"], "/"),
                                        "[", 1))), collapse = ",")
    for (a in alleles)
    {
      randhs <- na.omit(mhci[mhci$method == m & mhci$model == a & mhci$sample == "randHS", "ic50"])
      randnc <- na.omit(mhci[mhci$method == m & mhci$model == a & mhci$sample == "randNC", "ic50"])
      randmb <- na.omit(mhci[mhci$method == m & mhci$model == a & mhci$sample == "randMB", "ic50"])
      
      data <- na.omit(mhci[mhci$method == m & mhci$model == a & mhci$sample == s, "ic50"])
      if (length(data)>0)
      {
        wt <- wilcox.test(data, randmb, alternative = "less")
    
        ks <- ks.test(data, randmb, alternative = "less", exact = NULL) 
      
        ftpv = fisher(randmb, data)
      
        beanplot(filt1(data), filt1(randhs), filt1(randmb), filt1(randnc), names = c(s, "HS", "MB", "NC"), log = "")
        title(paste0(a, ", Wt=", round(wt$p.value, 3), ", KS=", round(ks$p.value,3), ", FT =", round(ftpv,3)))
      }
    }
    dev.off()
  }
  
  report[m] <- col
  ralleles[m] <- acol
}
write.table(report, file = paste0(outdir, "/mhc.txt"), sep = "\t", row.names = TRUE)
write.table(ralleles, file = paste0(outdir, "/mhc_a.tsv"), sep = "\t", row.names = TRUE)


mhci <- na.omit(mhci)
#mhci$score <- 2*50/(50+mhci$ic50)
mhci$score <- 2*1./(1+exp(mhci$ic50/100))
#mhci <- mhci[!startsWith(mhci$sample, "rand"),]
for (m in methods)
{
  s <- mhci %>% subset(method==m) %>% group_by(model, sample) %>% summarise(score = sum(score))
  pl <- lapply(samples, function(cs)
  {
    ss <- s[s$sample == cs,c("model", "score")]
    p <- ggplot(ss) + geom_bar(aes(x=model, y = score), stat="identity")+
      ggtitle(cs)+
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 4))
  })
  names(pl) <- samples
  pdf(paste0(outdir, "/all_", m, ".pdf"))
  invisible(lapply(pl, print))
  dev.off()
}

#### delta  ####
colourCount = length(unique(sapply(strsplit(mhci$model, '/'), "[", 1)))
mycolors <- rep(brewer.pal(name="Dark2", n = 8), colourCount %/% 8 +1)
relalleles <- list()
for (m in methods)
{
  s <- mhci %>% subset(method==m) %>% group_by(model, sample) %>% summarise(score = sum(score))
  if (m=="smm")
  {
    s <- s[s$model!="HLA-B*53:01/8",]  ## dirty
  }
  
  pivot <- s[s$sample == refset ,c("model", "score")]
  pivot <- rename(pivot, pivot = score)
  
  relNC <- s[s$sample == "randNC",c("model", "score")]
  relNC <- merge(relNC, pivot)
  relNC$rel <- relNC$score / relNC$pivot
  
  relHS <- s[s$sample == "randHS",c("model", "score")]
  relHS <- merge(relHS, pivot)
  relHS$rel <- relHS$score / relHS$pivot
  
  relMB <- s[s$sample == "randMB",c("model", "score")]
  relMB <- merge(relMB, pivot)
  relMB$rel <- relMB$score / relMB$pivot
  
  if (refset == "randHS"){
    std <- sd(relMB$rel)
  }else if (refset == "randMB"){
    std <- sd(relHS$rel)
  }else if (refset == "randNC"){
    std <- sd(relMB$rel)
  }else
    stop("Bad refset")
  
  s <- merge(s, pivot, by = "model", all.x = TRUE)
   
   rcol <- list()
    for (cs in samples)
   {
     ss <- s[s$sample == cs,c("model", "score", "pivot")]
     ss$rel <- ss$score/ss$pivot
     ss$rel <- ss$rel/mean(ss$rel)
     rcol[cs] <- paste(collapse = ";", ss[ss$rel > 1+3*std, "model"])
   }
   relalleles[[m]] <- rcol
   
  
  pl <- lapply(samples_all, function(cs)
  {
    ss <- s[s$sample == cs,c("model", "score", "pivot")]
    ss$rel <- ss$score/ss$pivot
    ss$rel <- ss$rel/mean(ss$rel)
    ss$al <- factor(sapply(strsplit(ss$model, '/'), "[", 1))
    p <- ggplot(ss) + geom_bar(aes(x=model, y = rel, fill = al), stat="identity")+
      geom_hline(yintercept=Filter(function(x) x>0, 
                                   c(mean(ss$rel)-std, mean(ss$rel+std))),
                 size = 1)+
      geom_hline(yintercept=Filter(function(x) x>0, 
                                   c(mean(ss$rel)-2*std, mean(ss$rel+2*std))),
                 size = 0.5)+
      geom_hline(yintercept=Filter(function(x) x>0, 
                                   c(mean(ss$rel)-3*std, mean(ss$rel+3*std))),
                 size = 0.25)+
      ggtitle(cs)+
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 4))+
      scale_fill_manual(values = mycolors)+
      guides(fill=FALSE)
  })
  names(pl) <- samples
  pdf(paste0(outdir, "/all_rel_", m, ".pdf"))
  invisible(lapply(pl, print))
  dev.off()
}
tt <- as.data.frame(sapply(relalleles, unlist))
write.table(tt, file = paste0(outdir, "/mhc_3sigma.tsv"), sep = "\t", row.names = TRUE)

####################################### II ###############################################
# mhcii <- read_mhc("h2")
# 
# 
# methods <- unique(mhcii$method)
# samples <- unique(mhcii$sample)
# alleles <- unique(mhcii$model)
# samples_all <- samples
# samples <- samples[!startsWith(samples, "rand")]
# options(digits=2)
# 
# #mhcii$ic50 <- ifelse(mhcii$ic50 > 20000, yes = NA, no = mhcii$ic50)
# 
# 
# report <- data.frame(row.names = samples)
# ralleles <- data.frame(row.names = samples)
# 
# for (m in methods)
# {
#   col = list()
#   acol = list()
#   for (s in samples)
#   {
#     print(paste(m, s))
#     pdf(file = paste0("hii/", s, "_", m, ".pdf"))
#     sc <- sapply(alleles, function(x){
#       mb <- na.omit(mhcii[mhcii$method == m & mhcii$model == x & mhcii$sample == "randHS", "ic50"])
#       dt <- na.omit(mhcii[mhcii$method == m & mhcii$model == x & mhcii$sample == s, "ic50"])
#       #wt <- wilcox.test(dt, mb, alternative = "less")
#       ftpv <- fisher(mb,dt)
#       return(ftpv)
#     }) 
#     
#     alleles <- alleles[order(sc)]
#     col[s] <- alleles[1]  
#     acol <- paste(unique(na.omit(sapply(strsplit(mhcii[mhcii$method == m & mhcii$ic50 < 50 & mhcii$sample == s, "model"], "/"),
#                                 "[", 1))), collapse = ",")
#     for (a in alleles)
#     {
#       randhs <- na.omit(mhcii[mhcii$method == m & mhcii$model == a & mhcii$sample == "randHS", "ic50"])
#       randnc <- na.omit(mhcii[mhcii$method == m & mhcii$model == a & mhcii$sample == "randNC", "ic50"])
#       randmb <- na.omit(mhcii[mhcii$method == m & mhcii$model == a & mhcii$sample == "randMB", "ic50"])
#       
#       data <- na.omit(mhcii[mhcii$method == m & mhcii$model == a & mhcii$sample == s, "ic50"])
#       if (length(data)>0)
#       {
#         wt <- wilcox.test(data, randmb, alternative = "less")
#         
#         ks <- ks.test(data, randmb, alternative = "less", exact = NULL) 
#         
#         ftpv = fisher(randmb, data)
#         
#         
#         beanplot(filt(data), filt(randhs), filt(randmb), filt(randnc), names = c(s, "HS", "MB", "NC"), log = "")
#         title(paste0(a, ", Wt=", round(wt$p.value, 3), ", KS=", round(ks$p.value,3), ", FT =", round(ftpv,3)))
#       }
#     }
#     dev.off()
#   }
#   
#   report[m] <- col
#   ralleles[m] <- acol
# }
# 
# mhcii <- na.omit(mhcii)
# #mhcii$score <- 2*50/(50+mhcii$ic50)
# #mhcii <- mhcii[!startsWith(mhcii$sample, "rand"),]
# mhcii$score <- 2*1./(1+exp(mhcii$ic50/100))
# for (m in methods)
# {
#   s <- mhcii %>% subset(method==m) %>% group_by(model, sample) %>% summarise(score = sum(score))
#   pl <- lapply(samples_all, function(cs)
#   {
#     ss <- s[s$sample == cs,c("model", "score")]
#     p <- ggplot(ss) + geom_bar(aes(x=model, y = score), stat="identity")+
#       ggtitle(cs)+
#       theme(axis.text.x = element_text(angle = 90, hjust = 1))
#   })
#   names(pl) <- samples
#   pdf(paste0("hii/all_", m, ".pdf"))
#   invisible(lapply(pl, print))
#   dev.off()
# }
# 
# #### delta  ####
# colourCount = length(unique(sapply(strsplit(mhcii$model, '/'), "[", 1)))
# mycolors <- rep(brewer.pal(name="Dark2", n = 8), colourCount %/% 8 +1)
# for (m in methods)
# {
#   s <- mhcii %>% subset(method==m) %>% group_by(model, sample) %>% summarise(score = sum(score))
#   pivot <- s[s$sample == "randHS",c("model", "score")]
#   pivot <- rename(pivot, pivot = score)
#   
#   relNC <- s[s$sample == "randNC",c("model", "score")]
#   relNC <- merge(relNC, pivot)
#   relNC$rel <- relNC$score / relNC$pivot
#   
#   relMB <- s[s$sample == "randMB",c("model", "score")]
#   relMB <- merge(relMB, pivot)
#   relMB$rel <- relMB$score / relMB$pivot
#   
#   std <- mean(c(sd(relMB$rel), sd(relNC$rel)))
#   
#   s <- merge(s, pivot, by = "model", all.x = TRUE)
#   pl <- lapply(samples_all, function(cs)
#   {
#     ss <- s[s$sample == cs,c("model", "score", "pivot")]
#     ss$rel <- ss$score/ss$pivot
#     ss$rel <- ss$rel/mean(ss$rel)
#     ss$al <- factor(sapply(strsplit(ss$model, '/'), "[", 1))
#     p <- ggplot(ss) + geom_bar(aes(x=model, y = rel, fill = al), stat="identity")+
#       geom_hline(yintercept=Filter(function(x) x>0, 
#                                    c(mean(ss$rel)-std, mean(ss$rel+std))),
#                  size = 1)+
#       geom_hline(yintercept=Filter(function(x) x>0, 
#                                    c(mean(ss$rel)-2*std, mean(ss$rel+2*std))),
#                  size = 0.5)+
#       geom_hline(yintercept=Filter(function(x) x>0, 
#                                    c(mean(ss$rel)-3*std, mean(ss$rel+3*std))),
#                  size = 0.25)+
#       ggtitle(cs)+
#       theme(axis.text.x = element_text(angle = 90, hjust = 1))+
#       scale_fill_manual(values = mycolors)+
#       guides(fill = FALSE)
#   })
#   names(pl) <- samples
#   pdf(paste0("hii/all_rel_", m, ".pdf"))
#   invisible(lapply(pl, print))
#   dev.off()
# }
# 
# #s_ <- mhcii %>% group_by(model, sample) %>% summarise(score = mean(score))
# #ss_ <- s_[s_$sample == "PMI-1",c("model", "score")]
# #ggplot(ss_) + geom_bar(aes(x=model, y = score), stat="identity")+theme(axis.text.x = element_text(angle = 90, hjust = 1))
# 
# 
# write.table(report, file = "hii/mhcii.txt", sep = "\t", row.names = TRUE)
# write.table(ralleles, file = "hii/mhcii_a.tsv", sep = "\t", row.names = TRUE)

