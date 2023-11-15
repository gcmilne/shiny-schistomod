#--------------------------------#
# SHINY APP TO USE SCHISTO MODEL #
#--------------------------------#

# load packages
library(shiny)
library(plotly)
library(ggplot2)
library(dplyr)
library(mvtnorm)
library(RcppNumerical)
library(bslib)     #bs_theme()
library(thematic)  #thematic_shiny()

# set working directory
setwd("~/GitHub/schisto-immunity-v2")

## load model from package (so compilation not needed)
library(SchistoTransmissionModel)  # RunTransmissionModel() to run transmission model

# Define the user interface ----------------------------------------------------
ui <- fluidPage(
  
  # change width of numeric input boxes
  tags$head(
    tags$style(HTML("
  input[type=\"number\"] {
    width: 100px;
  }

"))
  ),
  
  # set app theme
  theme = bs_theme(bootswatch = "darkly"),

  # create title
  titlePanel("Schistosomiasis MDA simulator"),
  
  fluidRow(
    
    column(
      width = 4,
      # Input sliders
      # Choose R0
      sliderInput("R0", "Basic reproduction number (R0)", min = 2, max = 10, step = 1, value = 3),
      # Choose immunity strength
      sliderInput("immune_strength", "Strength of acquired immunity", min = 0, max = 1, step = 0.1, value = 0),
      # Choose no. treatments
      sliderInput("numtx", "Total no. MDA rounds", min = 1, max = 20, step = 1, value = 10), 
      # Choose frequency of MDA
      sliderInput("freqtx", "No. years between MDA rounds", min = 0.25, max = 2, step = 0.25, value = 1),
      # Choose SAC treatment coverage
      sliderInput("coverage_sac", "MDA coverage in SAC (%)", min = 5, max = 100, step = 5, value = 75),
      # Choose adult treatment coverage
      sliderInput("coverage_adult", "MDA coverage in adults (%)", min = 0, max = 100, step = 5, value = 20),
      # Choose treatment ages
      sliderInput("txages", "Treatment ages (years)", min = 0, max = 50, step = 2, value = c(2, 35))
    ), 
    
    # Output plot
    column(
      width = 8,
      plotlyOutput("plot")
    )
    
  ), 
  
  fluidRow(
    # Choose model outputs
    checkboxGroupInput("choose_outputs", "Model output(s)",
                       choiceNames = list("EPG", "Prevalence", "Worm burden", "TAL1-IgE"),
                       choiceValues = list("epg", "prev", "W", "ige"), 
                       select = c("epg", "prev", "W", "ige")),
    # Display dynamics x years post-treatment
    numericInput(
      "years_post_last_tx",
      "No. years to simulate post-treatment",
      value = 5,
      min = 0,
      max = 20,
      step = NA
    )
  )
  
)

# Define the server ------------------------------------------------------------
server <- function(input, output, session) {
  
  # ensures ggplot2 automatically matches the app theme
  thematic_shiny()
  
  # set parameters outside of reactive environment
  params <- log(c(
    R0 = 3,
    R0_weight = 0.6, 
    NhNs = 0.6, 
    kW = 0.4, 
    decay_immunity = 7.5,
    protection_immunity = 0, 
    epg_constant = 5.81, 
    ige_constant = 0.5
  ))
  
  # Reactive function to calculate and update the plot based on input changes
  reactive_plot <- reactive({
    
    # set reactive parameters
    params["R0"] <- log(input$R0)
    params["protection_immunity"] <- log(input$immune_strength)
    
    # get treatment parameters from input
    tx_pars <- c(toggle_tx     = 1,
                 toggle_covdat = 0,   
                 start_tx      = 20,   
                 n_tx          = input$numtx,
                 freq_tx       = input$freqtx,
                 coverage      = input$coverage_sac/100, 
                 efficacy      = 0.95, 
                 min_tx_age    = min(input$txages),   
                 max_tx_age    = max(input$txages), 
                 cov_weight    = input$coverage_adult / input$coverage_sac
    )
    
    # get treatment times
    tx_times <- seq(tx_pars["start_tx"], tx_pars["start_tx"] + 
                      (tx_pars["n_tx"]-1)*tx_pars["freq_tx"], tx_pars["freq_tx"])
    
    # simulate the model using the input values
    sim <- RunTransmissionModel(
      theta    = params, 
      runtime  = tx_pars["start_tx"] + (max(tx_times) - tx_pars["start_tx"]) + input$years_post_last_tx, 
      stepsize = 1/8, 
      alltimes = 1, 
      tx_pars  = tx_pars, 
      tx_times = tx_times, 
      coverage_data = NA
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
          name == "epg"  ~ "Mean infection intensity (EPG)",
          name == "prev" ~ "Mean prevalence (%)",
          name == "ige"  ~ "Mean TAL1-IgE optical densitity"
        )
      ) 
    
    # create a new data frame with rounded values for hover-over display
    # dat_rounded <- dat %>%
    #   mutate(value_rounded = round(value, 1)) %>%
    #   mutate(text = paste(name, ": ", value_rounded))
    
    # make the plot
    ggplotly(
      # ggplot(dat_rounded, aes(x = time, y = value_rounded)) +
      ggplot(dat, aes(x = time, y = value)) +
        geom_line() + 
        labs(x = "Years since first treatment", y="") + 
        facet_wrap(~name, scales = "free_y", ncol = 2) + 
        theme( axis.text = element_text( size = 14 ),
               axis.text.x = element_text( size = 14 ),
               axis.title = element_text( size = 14, face = "bold" ),
               strip.background = element_blank(),
               strip.text = element_text(size = 14), 
               panel.spacing = unit(2, "lines"))#,
      # hoverinfo = 'text',  # Set hoverinfo to 'text' to use only the custom text defined in dat_rounded
      # text = ~text  # Specify the text to be displayed in the hover-over info
    )

  })
  
}

# Run the app ------------------------------------------------------------------
shinyApp(ui = ui, server = server)
# run_with_themer(shinyApp(ui = ui, server = server))
