#nbp server.R

shinyServer(function(input,output,session){
  
  source("./dev.R")
  
  plot.data <- eventReactive(input$go,{
    validate(
      need(input$gen<=5000,"Please select < 5000 generations."),
      need(input$nPop<=100,"Please select < 100 populations"),
      need(input$n<=1e6,"Please select n < 1e6"),
      need(input$plotStats!="","Select a variable to plot.")
      )
    runPopSim(gen=input$gen,p=input$p,Waa=input$Waa,Wab=input$Wab,Wbb=input$Wbb,n=input$n,nPop=input$nPop,m=input$m,stats=input$plotStats)
  })
  
  output$plot <- renderPlot({
    plotSingleRun(plot.data(),nPop=input$nPop,gen=input$gen)
  })
  
  sumTable <- eventReactive(input$runSim,{
    sumTable <- data.frame(matrix(ncol=3*input$nPop+8))
    withProgress(message="simulating populations...",value=0,{
      for(i in 1:100){
        df <- runPopSim.noMelt(gen=input$gen,p=input$p,Waa=input$Waa,Wab=input$Wab,Wbb=input$Wbb,n=input$n,nPop=2,m=input$m)
        names(sumTable) <- names(df)
        sumTable[i,] <- df[nrow(df),]
        incProgress(1/100)
      }
    })
    sumTable
  })
 
  tableData <- reactive({
    tbl <- colMeans(sumTable(),na.rm=T)
    tbl <- tbl[c("Fis","Hs","Ht","Fst")]
    tbl
  })
  
  output$table <- renderTable(t(tableData()),colnames = T,digits=4,caption = "Mean state at final generation:",
                              caption.placement = getOption("xtable.caption.placement", "top"),
                              caption.width = getOption("xtable.caption.width", NULL))
  output$sumTable <- renderTable(sumTable())
})
