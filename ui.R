library(shinythemes)
library(shiny)


tagList(
  navbarPage(theme = shinytheme("yeti"), 
             strong("Population Genetics Calculator"),
             navbarMenu(strong("Hardy - Weinberg Equilibrium"),
                      tabPanel("Enter data manually",
                               sidebarPanel(
                                p(strong("Chi-sq Hardy-Weinberg Equilibrium Test Calculator"), "for biallelic markers (like SNPs, indels etc.)"),
                                p("You only need to", span(strong("enter the observed numbers for each genotype"), style = "color:SteelBlue"),
                                  "and the result will update automatically."),
                                numericInput("num1", "Common homozygotes", value = NA),
                                numericInput("num2", "Heterozygotes", value = NA),
                                numericInput("num3", "Rare Homozygotes", value = NA),
                                p("Under the results there is a brief explanation of their meaning.")),
                               mainPanel(
                                 h4(strong("RESULTS PANEL"), align = "left"),
                                 hr(),
                                 tableOutput("pq"),
                                 tableOutput("table1"),
                                 tableOutput("table2"),
                                 hr(),
                                 p(strong("*** Basic information about how to interpret the results"), align = "right"),
                                 p("Testing deviation from the Hardy - Weinberg principle is generally performed using", 
                                   strong("Pearson's chi-squared test."), 
                                   "The following hypotheses are tested: ", strong("the null hypothesis"), 
                                   "is that the population is in Hardy Weinberg proportions and", 
                                   strong("the alternative hypothesis"), 
                                   "is that the population is not in Hardy Weinberg proportions."),
                                 p("There is 1 degree of freedom (df = number of genotypes variants - number of alleles variants). 
                            Assuming the level of significance at the level of 5%, we consider two options:"),
                                 p(strong("p - value > 0.05"), "means there is no reason to reject the null hypothesis."),
                                 p(strong("p - value < 0.05"), "means we reject the null hypothesis in favor of the alternative hypothesis.")
                               )),       
                      tabPanel("Load the file",
                               sidebarPanel(
                                 p(strong("Chi-sq Hardy-Weinberg Equilibrium Test Calculator"), "for biallelic markers (like SNPs, indels etc.)"),
                                 p("You only need to", span(strong("load a .xls or .xlsx file"), style = "color:SteelBlue"),
                                   "and then", span(strong("select the number of the column"), style = "color:SteelBlue"), 
                                   "whose content you want to analyze. In the last step,", 
                                   span(strong("select from the list the genotype variants"), style = "color:SteelBlue"), 
                                   "that are in the selected column."),
                                 p(code("The file must contain genotypes separated by a separator '/', e.g. G/A")),
                                 fileInput("file1", "Load the file",
                                           multiple = FALSE,
                                           accept = c(".xls", ".xlsx")),
                                 numericInput("colnum", "Choose the column number", value = 1),
                                 h6("Select the appropriate genotype variants:"),
                                 selectInput("gen1", "Homozygote 1",
                                              c(Choose = "",
                                                "AA" = "A/A",
                                                "GG" = "G/G",
                                                "CC" = "C/C",
                                                "TT" = "T/T"),
                                              selectize=TRUE),
                                 selectInput("gen2", "Homozygote 2",
                                              c(Choose = "",
                                                "AA" = "A/A",
                                                "GG" = "G/G",
                                                "CC" = "C/C",
                                                "TT" = "T/T"),
                                             selectize = TRUE),
                                 selectInput("gen3", "Heterozygote",
                                              c(Choose = "",
                                                "AG" = "A/G",
                                                "GA" = "G/A",
                                                "CT" = "C/T",
                                                "TC" = "T/C"),
                                             selectize = TRUE)),
                               mainPanel(
                                 h4(strong("RESULTS PANEL"), align = "left"),
                                 hr(),
                                 tableOutput("file_head"),
                                 tableOutput("pq2"),
                                 tableOutput("result"),
                                 tableOutput("chi"),
                                 hr(),
                                 p(strong("*** Basic information about how to interpret the results"), align = "right"),
                                 p("Testing deviation from the Hardy - Weinberg principle is generally performed using", 
                                   strong("Pearson's chi-squared test."), 
                                   "The following hypotheses are tested: ", strong("the null hypothesis"), 
                                   "is that the population is in Hardy Weinberg proportions and", 
                                   strong("the alternative hypothesis"), 
                                   "is that the population is not in Hardy Weinberg proportions."), 
                                p("There is 1 degree of freedom (df = number of genotypes variants - number of alleles variants). 
                                  Assuming the level of significance at the level of 5%, we consider two options:"),
                                p(strong("p - value > 0.05"), "means there is no reason to reject the null hypothesis."),
                                p(strong("p - value < 0.05"), "means we reject the null hypothesis in favor of the alternative hypothesis.")
                                ))),
             tabPanel(strong("Case - Control Studies"),
                      sidebarPanel(
                        p(strong("Case - Control Studies Calculator"), "is a tool that allows you to get the expected number of genotypes in the 
                          control and treatment(test) group on the basis of given observed numbers."),
                        p("You only need to", span(strong("enter the observed numbers for each genotype in the control and treatment group"), 
                                                   style = "color:SteelBlue"), "and the result will update automatically."),
                        hr(),
                        h6("For the control group:"),
                        numericInput("num4", "Common homozygote", value = NA),
                        numericInput("num5", "Heterozygote", value = NA),
                        numericInput("num6", "Rare Homozygote", value = NA),
                        hr(),
                        h6("For the treatment group:"),
                        numericInput("num7", "Common homozygote", value = NA),
                        numericInput("num8", "Heterozygote", value = NA),
                        numericInput("num9", "Rare Homozygote", value = NA)
                      ),
                      mainPanel(
                        h4(strong("RESULTS PANEL"), align = "left"),
                        hr(),
                        tableOutput("table3"),
                        tableOutput("table4"),
                        tableOutput("table5"),
                        tableOutput("table6"),
                        hr(),
                        p(strong("*** Explanation of the abbreviations included in the results"), align = "right"),
                        p(code("Obs_c"), "- observed counts in the control group", code("Obs_t"), "- observed counts in the treatment group",
                          code("Exp_c"), "- expected counts in the control group", code("Exp_t"), "- expected counts in the treatment group",
                          code("chi_square"), "- the Chi Square statistic", code("df"), "- degrees of freedom", code("OR"), "- odds ratio")
                      )),
             tabPanel(strong("Linkage Disequilibrium"),
                      sidebarPanel(
                        p(strong("Linkage Disequilibrium Calculator"), "between any pair of molecular loci (including SNPs)"),
                        p("You only need to", span(strong("load a .xls or .xlsx file."), style = "color:SteelBlue"), 
                          "Note: The file must contain 2 columns with genotypes separated by a separator '/', e.g. G/A"),
                        fileInput("fileLD", "Load the file",
                                  multiple = FALSE,
                                  accept = c(".xls", ".xlsx")),
                        p("You can", span(strong("save the received results"), style = "color:SteelBlue"), 
                          "in the form of HTML report by clicking on the button below"),
                        downloadButton("report", "Generate report")),
                      mainPanel(
                        h4(strong("RESULTS PANEL"), align = "left"),
                        hr(),
                        tableOutput("table_LD1"),
                        tableOutput("table_LD2"),
                        tableOutput("table_LD3"),
                        tableOutput("table_LD4"),
                        tableOutput("table_LD5"),
                        tableOutput("table_LD6"),
                        tableOutput("table_LD7"),
                        tableOutput("table_LD8"),
                        hr(),
                        p(strong("*** Explanation of the abbreviations included in the results"), align = "right"),
                        p(code("Alleles_freq"), "- alleles frequency (corresponds to p and q)", 
                          code("Exp_N"), "- expected number of haplotypes",
                          code("Exp_freq"), "- expected frequency of haplotypes",
                          code("Obs_N"), "- observed number of haplotypes",
                          code("Obs_freq"), "- observed frequency of haplotypes",
                          code("D_value"), "- the coefficient of linkage disequilibrium",
                          code("D_prim"), "- theoretical maximum difference between the observed and expected allele frequencies",
                          code("r_squared"), "- the correlation coefficient between pairs of loci",
                          code("chi_square"), "- the Chi Square statistic", code("df"), "- degrees of freedom")
                        )
             ))
)