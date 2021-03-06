shinyUI(fluidPage(
  titlePanel("driftR: Population Genetic Simulations in R"),
  sidebarLayout(
    sidebarPanel(
      numericInput("p","Starting allele frequency",value=0.5,min=0,max=1),
      sliderInput("Uab","Mutation Rate",value=0,min=0,max=0.1),
      sliderInput("Waa","Fitness of genotype AA",value=1,min=0,max=1),
      sliderInput("Wab","Fitness of genotype AB",value=1,min=0,max=1),
      sliderInput("Wbb","Fitness of genotype BB",value=1,min=0,max=1),
      sliderInput("m","Migration Rate",0,min=0,max=0.35),
      numericInput("n","Population Size",100,min=1,max=1e5),
      numericInput("nPop","Number of Populations",2,min=1,max=100),
      numericInput("gen","Number of Generations",100,min=1,max=5000),
      checkboxGroupInput(inputId="plotStats",label="plot:",choices=c("p","Hs","Ht","Fis","Fst"),inline=T,selected="p"),
      div(style="margin-left:40%;",checkboxInput("drift","Drift?",value = T)),
      actionButton("go","go",width="100%"),
      div(helpText("driftR simulates allele and genotype frequencies for a single biallelic locus in biological populations. 
                    Core functions adapted from the Java program popG (http://evolution.gs.washington.edu/popgen/popg.html). 
                    Full code available on github: https://github.com/cjbattey/popSim"),style="font-size:75%")
      ),
    
    mainPanel(
      plotOutput("plot"),
      textOutput("nLost"),
      helpText("Click the button below to run 100, 100-generation simulations of 2 populations using the current 
               parameters."),
      actionButton("runSim","Run Replicate Simulations"),
      tableOutput("meanTable"),
      tableOutput("varTable"),
      div(tableOutput("sumTable"), style = "font-size: 75%; width: 75%;")
    )
    )
  
))