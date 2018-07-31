library(readxl)
library(shiny)
library(tidyverse)
library(readxl)
library(stringr)
library(knitr)
library(kableExtra)

function(input, output) {
  p <- reactive({(2 * input$num1 + input$num2)/(2 * (input$num1 + input$num2 + input$num3))})
  q <- reactive({(2 * input$num3 + input$num2)/(2 * (input$num1 + input$num2 + input$num3))})
  alleles_freq <- reactive({
    req(input$num1, input$num2, input$num3)
    af <- cbind(Frequency_of_allele_A = p(),
                Frequency_of_allele_B = q(),
                sum = p() + q())
  })
  
  HW <- reactive({
    req(input$num1, input$num2, input$num3)
    p2 <- p()^2
    pq <- 2*p()*q()
    q2 <- q()^2
    N_obs <- rbind(input$num1, input$num2, input$num3)
    N_exp <- rbind(p2, pq, q2) * (input$num1 + input$num2 + input$num3)
    result <- data.frame(Genotype = c("Common homozygotes", "Heterozygotes", "Rare homozygotes"),
                         Observed_N = N_obs,
                         Expected_N = N_exp,
                         Observed_freq = N_obs / (input$num1 + input$num2 + input$num3),
                         Expected_freq = N_exp / (input$num1 + input$num2 + input$num3))
    final_result <- rbind(result, data.frame(Genotype="Total",t(colSums(result[,-1]))))
    list(df1 = N_obs, df2 = N_exp, df3 = final_result)
  })
                  
  output$pq <- renderTable(digits = 3, width = "200%", {
    alleles_freq()
  }, caption = "Table 1: Frequency of alleles",
  caption.placement = getOption("xtable.caption.placement", "top"), 
  caption.width = getOption("xtable.caption.width", NULL))
  
  output$table1 <- renderTable(digits = 3, width = "200%", {
    HW()[[3]]
  }, caption = "Table 2: Analysis of the number and frequency of genotypes",
  caption.placement = getOption("xtable.caption.placement", "top"), 
  caption.width = getOption("xtable.caption.width", NULL))
  
  output$table2 <- renderTable(digits = 3, width = "200%", {
    chi_squared_test <- sum(((HW()[[1]] - HW()[[2]]) ^ 2 / HW()[[2]]))
    p_value <- pchisq(chi_squared_test, 1, lower.tail=FALSE) 
    statistic <- data.frame(chi_square = chi_squared_test,
                            df = 1,
                            p_value = p_value)
  }, caption = "Table 3: Chi - squared test for deviation",
  caption.placement = getOption("xtable.caption.placement", "top"), 
  caption.width = getOption("xtable.caption.width", NULL))
  
  N_ccx <- reactive({data.frame(Genotype = c("Common homozygotes", "Heterozygotes", "Rare homozygotes"),
                                Obs_c = rbind(input$num4, input$num5, input$num6),
                                Obs_t = rbind(input$num7, input$num8, input$num9))})
  N_ccy <- reactive({data.frame(Genotype = c(paste("Common homozygotes", "Heterozygotes", sep = " + "), "Rare homozygotes"),
                                Obs_c = rbind((input$num4 + input$num5), input$num6),
                                Obs_t = rbind((input$num7 + input$num8), input$num9))})
  N_ccz <- reactive({data.frame(Genotype = c("Common homozygotes", paste("Heterozygotes", "Rare homozygotes", sep = " + ")),
                                Obs_c = rbind(input$num4, (input$num5 + input$num6)),
                                Obs_t = rbind(input$num7, (input$num8 + input$num9)))})
  N_cc <- reactive({
    req(input$num4, input$num5, input$num6, input$num7, input$num8, input$num9)
    N_cc <- N_ccx()
    N_cc <- cbind(N_cc, Sum_cols = rowSums(N_cc[,-1]))
    N_cc <- rbind(N_cc, c(NA, sum(N_cc$Obs_c), sum(N_cc$Obs_t), sum(N_cc$Sum_cols)))
    N_exp <- data.frame(Genotype = c("Common homozygotes", "Heterozygotes", "Rare homozygotes"),
                       Exp_c = round((N_cc$Obs_c[4] * N_cc$Sum_cols[-4])/N_cc$Sum_cols[4],3),
                       Exp_t = round((N_cc$Obs_t[4] * N_cc$Sum_cols[-4])/N_cc$Sum_cols[4],3))
    N_exp <- cbind(N_exp, Sum = rowSums(N_exp[,-1]))
    N_exp <- rbind(N_exp, c(NA, sum(N_exp$Exp_c), sum(N_exp$Exp_t), sum(N_exp$Sum)))
    Expected_controle <- paste(round((N_exp$Exp_c / N_exp$Exp_c[4]) * 100, 3), "%", sep = "")
    Expected_treatment <- paste(round((N_exp$Exp_t / N_exp$Exp_t[4]) * 100, 3), "%", sep = "")
    N_cc <- cbind(N_cc, Exp_c = Expected_controle)
    N_cc <- cbind(N_cc, Exp_t = Expected_treatment)
    N_all = N_cc[, colnames(N_cc[c(1:2,5,3,6,4)])]
    rownames(N_all$Genotype <- c("Common homozygotes", "Heterozygotes", "Rare homozygotes", "Sum_rows"))
    list(df1 = N_cc, df2 = N_exp, df3 = N_all)
  })
  N_cc2 <- reactive({
    req(input$num4, input$num5, input$num6, input$num7, input$num8, input$num9)
    N_cc2 <- N_ccy()
    N_cc2 <- cbind(N_cc2, Sum_cols = rowSums(N_cc2[,-1]))
    N_cc2 <- rbind(N_cc2,c(NA, sum(N_cc2$Obs_c), sum(N_cc2$Obs_t), sum(N_cc2$Sum_cols)))
    N_exp2 <- data.frame(Genotyp = c(paste("Common homozygotes", "Heterozygotes", sep = " + "), "Rare homozygotes"),
                         Exp_c = round((N_cc2$Obs_c[3] * N_cc2$Sum_cols[-3])/N_cc2$Sum_cols[3],3),
                         Exp_t = round((N_cc2$Obs_t[3] * N_cc2$Sum_cols[-3])/N_cc2$Sum_cols[3],3))
    N_exp2 <- cbind(N_exp2, Sum = rowSums(N_exp2[,-1]))
    N_exp2 <- rbind(N_exp2,c(NA, sum(N_exp2$Exp_c), sum(N_exp2$Exp_t), sum(N_exp2$Sum)))
    Expected_controle2 <- paste(round((N_exp2$Exp_c/N_exp2$Exp_c[3])*100,3), "%", sep = "")
    Expected_treatment2 <- paste(round((N_exp2$Exp_t/N_exp2$Exp_t[3])*100,3), "%", sep = "")
    N_cc2 <- cbind(N_cc2, Exp_c = Expected_controle2)
    N_cc2 <- cbind(N_cc2, Exp_t = Expected_treatment2)
    N_all2 = N_cc2[, colnames(N_cc2)[c(1:2,5,3,6,4)]]
    rownames(N_all2$Genotype <- c(paste("Common homozygotes", "Heterozygotes", sep = " + "), "Rare homozygotes", "Sum_rows"))
    list(df1 = N_cc2, df2 = N_exp2, df3 = N_all2)
  })
  N_cc3 <- reactive({
    req(input$num4, input$num5, input$num6, input$num7, input$num8, input$num9)
    N_cc3 <- N_ccz()
    N_cc3 <- cbind(N_cc3, Sum_cols = rowSums(N_cc3[,-1]))
    N_cc3 <- rbind(N_cc3,c(NA, sum(N_cc3$Obs_c), sum(N_cc3$Obs_t), sum(N_cc3$Sum_cols)))
    N_exp3 <- data.frame(Genotyp = c("Common homozygotes", paste("Heterozygotes", "Rare homozygotes", sep = " + ")),
                         Exp_c = round((N_cc3$Obs_c[3] * N_cc3$Sum_cols[-3])/N_cc3$Sum_cols[3],3),
                         Exp_t = round((N_cc3$Obs_t[3] * N_cc3$Sum_cols[-3])/N_cc3$Sum_cols[3],3))
    N_exp3 <- cbind(N_exp3, Sum = rowSums(N_exp3[,-1]))
    N_exp3 <- rbind(N_exp3,c(NA, sum(N_exp3$Exp_c), sum(N_exp3$Exp_t), sum(N_exp3$Sum)))
    Expected_controle3 <- paste(round((N_exp3$Exp_c/N_exp3$Exp_c[3])*100,3), "%", sep = "")
    Expected_treatment3 <- paste(round((N_exp3$Exp_t/N_exp3$Exp_t[3])*100,3), "%", sep = "")
    N_cc3 <- cbind(N_cc3, Exp_c = Expected_controle3)
    N_cc3 <- cbind(N_cc3, Exp_t = Expected_treatment3)
    N_all3 = N_cc3[, colnames(N_cc3)[c(1:2,5,3,6,4)]]
    rownames(N_all3$Genotype <- c("Common homozygotes", paste("Heterozygotes", "Rare homozygotes", sep = " + "), "Sum_rows"))
    list(df1 = N_cc3, df2 = N_exp3, df3 = N_all3)
  })
  
  output$table3 <- renderTable(width = "200%",{
    N_cc()[[3]]
  })
  
  output$table4 <- renderTable(width = "200%",{
    N_cc2()[[3]]
  })
  
  output$table5 <- renderTable(width = "200%",{
    N_cc3()[[3]]
  })
  
  output$table6 <- renderTable(digits = 4, width = "200%", {
    chi_test <- sum(((N_cc()[[1]][1:3,2:3] - N_cc()[[2]][1:3,2:3]) ^ 2) / N_cc()[[2]][1:3,2:3])
    df = 2
    p_val <- pchisq(chi_test, df, lower.tail = FALSE)
    
    chi_test2 <- sum(((N_cc2()[[1]][1:2,2:3] - N_cc2()[[2]][1:2,2:3]) ^ 2) / N_cc2()[[2]][1:2,2:3])
    df2 = 1
    p_val2 <- pchisq(chi_test2, df2, lower.tail = FALSE) 
    or2 = (N_cc2()[[1]]$Obs_t[1] / N_cc2()[[1]]$Obs_t[2]) / (N_cc2()[[1]]$Obs_c[1]/N_cc2()[[1]]$Obs_c[2])
    
    chi_test3 <- sum(((N_cc3()[[1]][1:2,2:3] - N_cc3()[[2]][1:2,2:3]) ^ 2) / N_cc3()[[2]][1:2,2:3])
    df3 = 1
    p_val3 <- pchisq(chi_test3, df3, lower.tail = FALSE) 
    or3 = (N_cc3()[[1]]$Obs_t[1] / N_cc3()[[1]]$Obs_t[2]) / (N_cc3()[[1]]$Obs_c[1] / N_cc3()[[1]]$Obs_c[2])
    
    chi_test_vec <- c(chi_test, chi_test2, chi_test3)
    df_vec <- c(df, df2, df3)
    p_value_vec <- c(p_val, p_val2, p_val3)
    or_vec <- c('-', or2, or3)
    stat <- data.frame(chi_square = chi_test_vec,
                       df = df_vec,
                       p_value = p_value_vec,
                       OR = or_vec)
    rownames(stat) <- c("Tabela 1", "Tabela 2", "Tabela 3")
    stat
  })
  
  file_input <- reactive({
    req(input$file1)
    tryCatch(
      {
        df <- read_excel(input$file1$datapath)
        dataset <- df[,input$colnum, drop=FALSE]
        dataset <- na.omit(dataset)
        colnames(dataset) <- c("genotype")
      },
      error = function(e) {
        stop(safeError(e))
      }
    )
    list(df1 = df, df2 = dataset)
  })
  
  HW_auto <- reactive({
    req(input$gen1, input$gen2, input$gen3)
    dataset <- file_input()[[2]]
    freq <- as.data.frame(table(dataset))
    freq_sort <- rbind(freq[freq$dataset == input$gen1,],
                       freq[freq$dataset == input$gen3,])
    freq_sort <- freq_sort[order(freq_sort$Freq),]
    N_obs <- rbind(freq_sort,
                   freq[freq$dataset == input$gen2,])
    rownames(N_obs) <- NULL
    p <- (2 * N_obs$Freq[2] + N_obs$Freq[3])/(nrow(dataset)*2)
    q <- (2 * N_obs$Freq[1] + N_obs$Freq[3])/(nrow(dataset)*2)
    alleles_freq <- cbind(Frequency_of_allele_A = p,
                          Frequency_of_allele_B = q,
                          sum = p + q)
    p2 <- p^2
    pq <- 2*p*q
    q2 <- q^2
    N_exp <- (rbind(q2, p2, pq) * nrow(dataset))
    result <- data.frame(Genotype = N_obs$dataset,
                         N_observed = N_obs$Freq,
                         N_expected = N_exp,
                         Observed_freq = N_obs$Freq/nrow(dataset),
                         Expected_freq = N_exp/nrow(dataset))
    final_result <- rbind(result, data.frame(Genotype="Total",t(colSums(result[,-1]))))
    list(df1 = alleles_freq, df2 = final_result)
  })
  
  output$file_head <- renderTable(width = "200%", {
    head(file_input()[[1]],5)
  }, caption = "Data from the loaded file has the following form:",
  caption.placement = getOption("xtable.caption.placement", "top"), 
  caption.width = getOption("xtable.caption.width", NULL))
  
  output$pq2 <- renderTable(digits = 3, width = "200%", {
    HW_auto()[[1]]
  }, caption = "Table 1: Frequency of alleles",
  caption.placement = getOption("xtable.caption.placement", "top"), 
  caption.width = getOption("xtable.caption.width", NULL))
  
  output$result <- renderTable(digits = 3, width = "200%", {
    HW_auto()[[2]]
  }, caption = "Table 2: Analysis of the number and frequency of genotypes",
  caption.placement = getOption("xtable.caption.placement", "top"), 
  caption.width = getOption("xtable.caption.width", NULL))
  
  output$chi <- renderTable(digits = 3, width = "200%", {
    chi_squared_test <- sum(((HW_auto()[[2]]$N_observed - HW_auto()[[2]]$N_expected) ^ 2 / HW_auto()[[2]]$N_expected))
    p_value <- pchisq(chi_squared_test, 1, lower.tail=FALSE) 
    statistic <- data.frame(chi_square = chi_squared_test,
                            df = 1,
                            p_value = p_value)
  }, caption = "Table 3: Chi - squared test for deviation",
  caption.placement = getOption("xtable.caption.placement", "top"), 
  caption.width = getOption("xtable.caption.width", NULL))
  
  
  file_input2 <- reactive({
    req(input$fileLD)
    tryCatch(
      {
        df <- read_excel(input$fileLD$datapath)
        dataset <- na.omit(df)
        colnames(dataset) <- c("SNP1","SNP2")
        dataset
        },
      error = function(e) {
        stop(safeError(e))
        }
      )
    })
  
  LD <- reactive({
    dataset <- file_input2()
    freq1 <- as.data.frame(table(dataset$SNP1))
    freq2 <- as.data.frame(table(dataset$SNP2))
    levels_gen1 = as.vector(freq1$Var1)
    levels_gen2 = as.vector(freq2$Var1)
    freq1_sort <- rbind(freq1[freq1$Var1 == levels_gen1[1],],
                        freq1[freq1$Var1 == levels_gen1[3],])
    freq1_sort <- rbind(freq1_sort[order(freq1_sort$Freq),],
                        freq1[freq1$Var1 == levels_gen1[2],])
    freq2_sort <- rbind(freq2[freq2$Var1 == levels_gen2[1],],
                        freq2[freq2$Var1 == levels_gen2[3],])
    freq2_sort <- rbind(freq2_sort[order(freq2_sort$Freq),],
                        freq2[freq2$Var1 == levels_gen2[2],])
    genotypes <- cbind(freq1_sort, freq2_sort)
    colnames(genotypes) <- c("Genotype_1", "N1", "Genotype_2", "N2")
    rownames(genotypes) <- c(1:3)
    split1 <- str_split_fixed(dataset$SNP1, "/", 2)
    split2 <- str_split_fixed(dataset$SNP2, "/", 2)
    split_all <- data.frame(cbind(split1, split2))
    colnames(split_all) <- c("SNP1.1", "SNP1.2", "SNP2.1", "SNP2.2")
    haplo <- data.frame(cbind(split_all, h1 = paste(split_all$SNP1.1, split_all$SNP2.1, sep = ""),
                              h2 = paste(split_all$SNP1.1, split_all$SNP2.2, sep = ""),
                              h3 = paste(split_all$SNP1.2, split_all$SNP2.1, sep = ""),
                              h4 = paste(split_all$SNP1.2, split_all$SNP2.2, sep = "")))
    binding <- unlist(haplo[c("h1","h2","h3","h4")], use.names = FALSE)
    obs = sort(table(binding))
    obsvec = as.vector(obs)
    freq = obsvec/length(binding)
    pA <- (2 * freq1_sort$Freq[2] + freq1_sort$Freq[3])/(nrow(dataset) * 2)
    qa <- (2 * freq1_sort$Freq[1] + freq1_sort$Freq[3])/(nrow(dataset) * 2)
    pB <- (2 * freq2_sort$Freq[2] + freq2_sort$Freq[3])/(nrow(dataset) * 2)
    qb <- (2 * freq2_sort$Freq[1] + freq2_sort$Freq[3])/(nrow(dataset) * 2)
    pAB <- (pA * pB)  
    pAb <- (pA * qb) 
    paB <- (qa * pB) 
    pab <- (qa * qb) 
    alleles_freq <- data.frame(Allele = rbind(levels(split_all[,1])[2],levels(split_all[,1])[1],"Alleles_freq"),
                               a = rbind(round(pAB,5), round(paB,5), round(pB,5)),
                               b = rbind(round(pAb,5), round(pab,5), round(qb,5)),
                               SUM = rbind(round(pA,5), round(qa,5), pA + qa))
    colnames(alleles_freq) <- c("Allele", levels(split_all[,3])[2], levels(split_all[,3])[1], "Alleles_freq")
    haplotypes <- data.frame(Haplotype = rbind(names(obs)[4],
                                             names(obs)[3],
                                             names(obs)[2],
                                             names(obs)[1]),
                            Exp_N = rbind(pAB * (nrow(dataset) * 4),
                                               pAb * (nrow(dataset) * 4),
                                               paB * (nrow(dataset) * 4),
                                               pab * (nrow(dataset) * 4)),
                            Exp_freq = rbind(pAB,
                                          pAb,
                                          paB,
                                          pab),
                            Obs_N = rbind(obsvec[4],
                                                obsvec[3],
                                                obsvec[2],
                                                obsvec[1]),
                            Obs_freq = rbind(obsvec[4]/sum(obsvec),
                                          obsvec[3]/sum(obsvec),
                                          obsvec[2]/sum(obsvec),
                                          obsvec[1]/sum(obsvec)))
    D_value = haplotypes$Obs_freq - haplotypes$Exp_freq
    haplotypes <- cbind(haplotypes, D_value)
    haplotypes <- rbind(haplotypes, data.frame(Haplotype="Sum",round(t(colSums(haplotypes[,-1])),3)))
    rownames(haplotypes) <- NULL
    D = haplotypes$Obs_freq[1]*haplotypes$Obs_freq[4] - haplotypes$Obs_freq[2]*haplotypes$Obs_freq[3]
    alleles_freq_D<- data.frame(Allele = rbind(levels(split_all[,1])[2],levels(split_all[,1])[1],"SUM"),
                                a = rbind(pAB + D, paB - D, pB),
                                b = rbind(pAb - D, pab + D, qb),
                                SUM = rbind(pA, qa, pA + qa))
    colnames(alleles_freq_D) <- c("Allele", levels(split_all[,3])[2], levels(split_all[,3])[1], "SUM")
    rownames(alleles_freq_D) <- NULL
    Dmax = 0
    if (D >= 0){
      a <- pA * qb
      b <- qa * pB
      if (a < b){
        Dmax = a
      }
      else {
        Dmax = b
      }
    } else if(D < 0){
      a <- -(pA * pB)
      b <- -(qa * qb)
      if (a > b){
        Dmax = a
      }
      else {
        Dmax = b
      }
    }
    Dprim <- D / Dmax
    r2 <- D ^ 2 / (pA * qa * pB * qb)
    chi_square <- r2 * (nrow(dataset) * 4)
    p_val <- pchisq(chi_square, 1, lower.tail=FALSE)
    ld <- cbind(D = D, D_prim = Dprim, r_squared = r2, 
                chi_square = chi_square, p_value = p_val)
    nk1 = levels(split_all[,1])[1] #C
    nk2 = levels(split_all[,1])[2] #T
    nk3 = levels(split_all[,3])[1] #A
    nk4 = levels(split_all[,3])[2] #G
    n1 <- nrow(subset(split_all, !(SNP1.1==nk1 | SNP1.2==nk1 | SNP2.1==nk3 | SNP2.2==nk3))) #TTGG
    n2 <- nrow(subset(split_all, !(SNP1.1==nk1 | SNP1.2==nk1 | SNP2.1==nk4 | SNP2.2==nk4))) #TTAA
    n3 <- nrow(subset(split_all, !(SNP1.1==nk2 | SNP1.2==nk2 | SNP2.1==nk3 | SNP2.2==nk3))) #CCGG
    n4 <- nrow(subset(split_all, !(SNP1.1==nk2 | SNP1.2==nk2 | SNP2.1==nk4 | SNP2.2==nk4))) #CCAA
    n5 <- nrow(subset(split_all, !(SNP1.1==nk1 | SNP1.2==nk2 | SNP2.1==nk3 | SNP2.2==nk4))) #TCGA
    n6 <- nrow(subset(split_all, !(SNP1.1==nk1 | SNP1.2==nk2 | SNP2.1==nk3 | SNP2.2==nk3))) #TCGG
    n7 <- nrow(subset(split_all, !(SNP1.1==nk1 | SNP1.2==nk2 | SNP2.1==nk4 | SNP2.2==nk4))) #TCAA
    n8 <- nrow(subset(split_all, !(SNP1.1==nk1 | SNP1.2==nk1 | SNP2.1==nk3 | SNP2.2==nk4))) #TTGA
    n9 <- nrow(subset(split_all, !(SNP1.1==nk2 | SNP1.2==nk2 | SNP2.1==nk3 | SNP2.2==nk4))) #CCGA
    genotypes2 <- data.frame(Genotype = rbind(levels(genotypes$Genotype_1)[3], levels(genotypes$Genotype_1)[2],
                                            levels(genotypes$Genotype_1)[1]),
                            a = rbind(n2, n7, n4),
                            b = rbind(n8, n5, n9),
                            c = rbind(n1, n6, n3))
    colnames(genotypes2) <- c("Genotype", levels(genotypes$Genotype_2)[3], levels(genotypes$Genotype_2)[2],
                             levels(genotypes$Genotype_2)[1])
    genotypes2 <- cbind(genotypes2, SUM = rowSums(genotypes2[,-1]))
    genotypes2 <- rbind(genotypes2,c(NA, colSums(genotypes2[,-1])))
    rownames(genotypes2) <- NULL
    gen2_exp <- data.frame(Genotype = rbind(levels(genotypes$Genotype_1)[3], levels(genotypes$Genotype_1)[2],
                                            levels(genotypes$Genotype_1)[1]),
                           a = cbind((genotypes2[4,2]*genotypes2$SUM[-4])/genotypes2$SUM[4]),
                           b = cbind((genotypes2[4,3]*genotypes2$SUM[-4])/genotypes2$SUM[4]),
                           c = cbind((genotypes2[4,4]*genotypes2$SUM[-4])/genotypes2$SUM[4]))
    colnames(gen2_exp) <- c("Genotype", levels(genotypes$Genotype_2)[3], levels(genotypes$Genotype_2)[2],
                                 levels(genotypes$Genotype_2)[1])
    gen2_exp <- cbind(gen2_exp, SUM = rowSums(gen2_exp[,-1]))
    gen2_exp <- rbind(gen2_exp,c(NA, colSums(gen2_exp[,-1])))
    chi_squaredld = sum(((genotypes2[1:3,2:4] - gen2_exp[1:3,2:4]) ^ 2) / gen2_exp[1:3,2:4])
    df = 4
    p_value2 = pchisq(chi_squaredld,df,lower.tail=FALSE) 
    stat <- cbind(chi_square = chi_squaredld,
                  df = df,
                  p_value = p_value2)
    list(df1 = genotypes, df2 = alleles_freq, df3 = haplotypes, df4 = alleles_freq_D, 
         df5 = ld, df6 = genotypes2, df7 = gen2_exp, df8 = stat)
  })
  
  output$table_LD1 <- renderTable(digits = 3, width = "200%", {
    LD()[[1]]
  }, caption = "Table 1: Number of individual genotypes",
  caption.placement = getOption("xtable.caption.placement", "top"), 
  caption.width = getOption("xtable.caption.width", NULL))
  
  output$table_LD2 <- renderTable(digits = 3, width = "200%", {
    LD()[[2]]
  }, caption = "Table 2: The expected frequency of haplotypes + alleles frequency",
  caption.placement = getOption("xtable.caption.placement", "top"), 
  caption.width = getOption("xtable.caption.width", NULL))
  
  output$table_LD3 <- renderTable(digits = 3, width = "200%", {
    LD()[[3]]
  }, caption = "Table 3: The expected and observed frequency (and the number) of haplotypes",
  caption.placement = getOption("xtable.caption.placement", "top"), 
  caption.width = getOption("xtable.caption.width", NULL))
  
  output$table_LD4 <- renderTable(digits = 3, width = "200%", {
    LD()[[4]]
  }, caption = "Table 4: The expected frequency of haplotypes in terms of linkage disequilibrium 
                          + alleles frequency",
  caption.placement = getOption("xtable.caption.placement", "top"), 
  caption.width = getOption("xtable.caption.width", NULL))
  
  output$table_LD5 <- renderTable(digits = 3, width = "200%", {
    LD()[[5]]
  }, caption = "Table 5: Basic statistics of linkage disequilibrium",
  caption.placement = getOption("xtable.caption.placement", "top"), 
  caption.width = getOption("xtable.caption.width", NULL))
  
  output$table_LD6 <- renderTable(width = "200%", {
    LD()[[6]]
  }, caption = "Table 6: Combined (observed) number of genotypes",
  caption.placement = getOption("xtable.caption.placement", "top"), 
  caption.width = getOption("xtable.caption.width", NULL))
  
  output$table_LD7 <- renderTable(digits = 3, width = "200%", {
    LD()[[7]]
  }, caption = "Table 7: Combined (expected) number of genotypes",
  caption.placement = getOption("xtable.caption.placement", "top"), 
  caption.width = getOption("xtable.caption.width", NULL))
  
  output$table_LD8 <- renderTable(digits = 3, width = "200%", {
    LD()[[8]]
  }, caption = "Table 8: Chi - square test for linkage disequilibrium between two loci",
  caption.placement = getOption("xtable.caption.placement", "top"), 
  caption.width = getOption("xtable.caption.width", NULL))
  
  
  library(rmarkdown)
  output$report <- downloadHandler(
    filename = "report.html",
    
    content = function(file) {
      src <- normalizePath('report.Rmd')
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      file.copy(src, 'report.Rmd', overwrite = TRUE)
      
      out <- render("report.Rmd", "html_document" )
      file.rename(out, file)
    })
}