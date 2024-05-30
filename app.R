#--------------------------------#
# SHINY APP TO USE SCHISTO MODEL #
#--------------------------------#

# load packages
pacman::p_load(shiny, plotly, ggplot2, dplyr, bslib, thematic, SchistoTransmissionModel)
# library(shiny)
# library(plotly)
# library(ggplot2)
# library(dplyr)
# library(mvtnorm)
# library(RcppNumerical)
# library(bslib)     #bs_theme()
# library(thematic)  #thematic_shiny()
# library(tidyverse)
# 
# # set working directory
# setwd("~/GitHub/shiny-schistomod")
# 
# ## load model from package (so compilation not needed)
# library(SchistoTransmissionModel)  #RunTransmissionModel()

# Define the user interface ----------------------------------------------------
ui <- fluidPage(
  
  # set app theme
  theme = bs_theme(bootswatch = "darkly"),

  # create title
  titlePanel("Schistosomiasis MDA simulator"),
  
  fluidRow(
    
    column(
      h4(style = "padding:5px;"), #add vertical offset
      
      width = 4,
      # Input sliders
      # Choose R0
      sliderInput("R0", "Basic reproduction number (R0)", min = 2, max = 10, step = 1, value = 3),
      # Choose immunity strength
      sliderInput("immune_strength", "Strength of acquired immunity", min = 0, max = 1, step = 0.1, value = 0),
      # Choose no. treatments
      sliderInput("numtx", "Number of MDA rounds", min = 1, max = 20, step = 1, value = 10), 
      # Choose frequency of MDA
      sliderInput("freqtx", "Years between MDA rounds", min = 0.25, max = 2, step = 0.25, value = 1)
    
    ), 
    
    # Output plot
    column(
      
      width = 8,
      plotlyOutput("plot")
      
    )
    
  ), 
  
  fluidRow(
    
    column(
      
      width = 4,
      # Choose SAC treatment coverage
      sliderInput("coverage_sac", "MDA coverage in SAC (%)", min = 5, max = 100, step = 5, value = 75),
      # Choose adult treatment coverage
      sliderInput("coverage_adult", "MDA coverage in adults (%)", min = 0, max = 100, step = 5, value = 20)
    ), 
    
    column(
      
      width = 3, 
      # offset = 1, #add horizontal offset
      # Choose model outputs
      checkboxGroupInput("choose_outputs", "Model output(s)",
                         choiceNames = list("EPG", "Prevalence", "Worm burden", "TAL1-IgE"),
                         choiceValues = list("epg", "prev", "W", "ige"), 
                         select = c("epg", "prev", "W", "ige"))
      
    ),
    
    column(
      
      width = 4, 
      # Choose time format
      radioButtons(
        inputId = "years_post_last_tx",
        label = "Number of years to simulate after MDA",
        choices = seq(5,20,5)
      )
      
    ),
    
    # Footer to define acronyms
    hr(),
    print(
      "Acronyms: epg, eggs per gram of faeces; 
      MDA, mass drug administration; 
      SAC, school-age children (5-15 years old);
      smTAL1, Schistosoma mansoni tegument allergen-like protein 1"
      )
  )
  
)

# Define the server ------------------------------------------------------------
server <- function(input, output, session) {
  
  # ensures ggplot2 automatically matches the app theme
  thematic_shiny()
  
  # set parameters outside of reactive environment
  params <- c(
    R0 = 3,
    R0_weight = 0.6, 
    NhNs = 0.6, 
    kW = 0.4, 
    decay_immunity = 7.5,
    protection_immunity = 0, 
    epg_constant = 5.81, 
    ige_constant = 0.5
  )
  
  # Reactive function to calculate and update the plot based on input changes
  reactive_plot <- reactive({
    
    # set reactive parameters
    params["R0"] <- input$R0
    params["protection_immunity"] <- input$immune_strength
    
    # get treatment parameters from input
    tx_pars <- c(
      input_tx_times= 0,    #toggle for user inputting of treatment times (1) or calculation from input parameters (0)
      toggle_tx     = 1,    #toggle MDA treatment (1) or not (0)
      start_tx      = 20,   #model time to start treatment
      n_tx          = input$numtx,  #no. treatments (if input_tx_times==0)
      freq_tx       = input$freqtx, #no. years between each treatment event (if input_tx_times==0)
      sac_coverage  = input$coverage_sac/100,  #treatment coverage in SAC
      efficacy      = 0.95, #praziquantel efficacy
      min_tx_age    = 2,    #minimum treatment age 
      max_tx_age    = 51,   #maximum treatment age 
      cov_weight    = input$coverage_adult / input$coverage_sac
    )
    
    # get treatment times
    tx_times <- seq(tx_pars["start_tx"], tx_pars["start_tx"] + 
                      (tx_pars["n_tx"]-1)*tx_pars["freq_tx"], tx_pars["freq_tx"])
    
    # simulate the model using the input values
    sim <- RunTransmissionModel(
      theta    = params,
      tx_pars  = tx_pars, 
      runtime  = tx_pars["start_tx"] + (max(tx_times) - tx_pars["start_tx"]) + as.numeric(input$years_post_last_tx), 
      stepsize = 1/8, 
      user_tx_times = NA,
      user_cov_weight = NA,
      time_extract_states = NA,
      init_female_states = matrix(NA),
      init_male_states = matrix(NA)
    )
    
    # calculate indices to extract time-dependent model outputs
    indices <- sim$time > (tx_pars["start_tx"]-1)
    
    # store the output as a dataframe in long format
    dat <- data.frame(
      "age"  = sim$age             [indices],
      "time" = sim$time            [indices] - tx_pars["start_tx"], #present as "years since first treatment"
      "W"    = sim$worm_burden     [indices],
      "epg"  = sim$epg             [indices],
      "prev" = (sim$prevalence*100)[indices], 
      "ige"  = sim$IgE             [indices]
    ) %>% tidyr::pivot_longer(cols=c("W", "epg", "prev", "ige"))
    
    return(dat)
  })
  
  # Plot the output reactively
  output$plot <- renderPlotly({
    
    # load in the reactive data 
    dat <- reactive_plot() %>%
      # filter data by user-selected model outputs
      filter(name %in% input$choose_outputs) %>%
      # rename variables for plotting
      mutate(
        name = case_when(
          name == "W"    ~ "Mean worm burden", 
          name == "epg"  ~ "Mean infection intensity (epg)",
          name == "prev" ~ "Mean prevalence (%)",
          name == "ige"  ~ "Mean TAL1-IgE optical densitity"
        )
      ) 
    
    # make the plot 
    ggplotly(
      ggplot(dat,
             aes(
               x = time,
               y = value,
               group = 1,  #use dummy group variable so geom_line works properly
               text = paste0("Years since first treatment : ", time, "\n", 
                             name, ": ", round(value, 3))
             )
      ) +
        geom_line() +
        labs(x = "Years since first treatment", y="") +
        facet_wrap(~name, scales = "free_y", ncol = 2) +
        theme(
          axis.text        = element_text(size = 14),
          axis.title       = element_text(size = 14, face = "bold"),
          strip.text       = element_text(size = 14),
          strip.background = element_blank(),
          panel.spacing    = unit(2, "lines"),  #space between facets
          plot.margin      = grid::unit(c(0, 1, 0, 0), "cm") #space on right
        ),
      dynamicTicks = TRUE,
      tooltip = "text"
    )
    

  })
  
}

# Run the app ------------------------------------------------------------------
shinyApp(
  ui = ui, server = server
  , options = list(launch.browser = TRUE)
)
