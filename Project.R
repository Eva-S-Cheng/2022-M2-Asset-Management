# Libraries
library(dplyr)
library(shiny)
library(shinydashboard)
library(quadprog)
library(Hmisc)
library(DT)
library(reshape2)  
#conflict_prefer("filter", "dplyr")

## Loading the dataset
stock_price <- read.csv(file = 'Price_per-stock.csv', header=T)
names(stock_price)[1] = 'Date'
stock_price$Date=as.Date(stock_price$Date,"%d/%m/%Y")
industries <- read.csv(file = "Stocks_CAC.csv", header=T)
names(industries)[1] = 'Ticker'

## Important variables 
industry_list <- industries %>%
  count(Industry)

## Creating the functions 
# Calculate the returns
Compute_return <- function(Price){
  n = length(Price)
  ret = (Price[2:n]/Price[1:(n-1)]) - 1
  ret = c(0, ret)
  return (ret)
}

# Compound Annual Growth Rate
Compute_CAGR <- function(Price){
  n = length(Price)
  # In one year, we assume there are 252 trading days
  ret = (Price[n]/Price[1])^(252/(n-1))-1
  return (ret)
}

# Abnormal returns compared to the CACT
Compute_ab_return <- function(r, CACT=returns[,ncol(returns)-1]){
  n = length(r)
  # In one year, we assume there are 252 trading days
  ret = r[1:n] - CACT[1:n]
  return (ret)
}

# Get the datasets by industry
Get_data_industry <- function(industry_title, list_industries, dataset){
  tempo <-list_industries %>%
    filter(Industry %in% industry_title)%>%
    select(Ticker)
  ret <- dataset %>% select(c((tempo[,1])))
}

# Return the minimum
minimum <- function(a, b) {
  if (a > b) {
    return (b)
  }
  else {
    return (a)
  }
}

# Sortino : It is the risk adjusted return depending on how much return the user expects
Compute_sortino <- function(adjusted_return) {
  n = length(adjusted_return)
  mean_adjusted_ret = mean(adjusted_return[1:n])
  sqrt_mean_sq_adjusted_ret_zero = sqrt(mean(minimum(adjusted_return[1:n],0)^2))
  ret <- mean_adjusted_ret/sqrt_mean_sq_adjusted_ret_zero
  return (ret)
}

# Sharpe ratio : Average retrun depending on the standard deviation
Compute_sharpe <- function(returns) {
  n = length(returns)
  average_return = mean(Compute_ab_return(returns[1:n]))
  std = sd(Compute_ab_return(returns[1:n]))
  ret <- average_return/std
  return (ret)
}

# Annual covariance (used for Markowitz)
Compute_sigma <- function(returns) {
  return (cov(returns)*252)
}

# Return according to volatility efficient frontier, matrice as input
Get_data_fictive_portfolio <- function(selected_stock) {
  selected_returns <- apply(selected_stock, 2, FUN = Compute_return)
  # Number of assets selected
  num_asset = dim(selected_returns)[2]
  # Lists that contain the returns and the volatilities
  list_Ret = c()
  list_Vol = c()
  # Compute the CAGRS
  cagrs_mat = apply(selected_stock, 2, FUN = Compute_CAGR)
  # Covariance matrix
  sigma_mat = Compute_sigma(selected_returns)

  # Try with 10000 values
  for(i in 1:10000){
    # Allocation / weights are randomly given 
    weights <- runif(num_asset)
    weights <- as.matrix(weights/sum(weights))
    
    # Volatility and return 
    ret_portfolio = cagrs_mat %*% weights
    vol_portfolio = sqrt(t(weights) %*% (sigma_mat %*% weights))
    
    list_Ret = c(list_Ret, as.numeric(ret_portfolio))
    list_Vol = c(list_Vol, as.numeric(vol_portfolio))
  }
  return (data.frame(list_Ret, list_Vol))
}

# Using quadratic programming to get the efficient frontier 
Get_Markowitz_frontier <- function(selected_stock) {
  # sequence that will be used for the efficient frontier, 1000 values
  gammas = seq(from = 0.001, to = 0.999, by = 0.001)
  selected_returns <- apply(selected_stock, 2, FUN = Compute_return)
  # Number of assets selected
  num_asset = dim(selected_returns)[2]
  # Lists that contain the returns and the volatilities
  list_Ret_QP = c()
  list_Vol_QP = c()
  list_Sharp_QP = c()
  # Compute the CAGRS, matrix used for the function that will be minimized
  cagrs_mat = apply(selected_stock, 2, FUN = Compute_CAGR)
  # Covariance matrix will be Dmat, symetric matrix
  sigma_mat = Compute_sigma(selected_returns)
  # Matrix reprezenting and dependent the number of assets, because we want to change the weights of each asset
  Amat = cbind(rep(1, num_asset), diag(num_asset))
  # Constraint 1
  meq = 1
  # Default value, 0
  bvec = c(1, rep(0, num_asset))
  
  # Vector containing the CAGRS of each asset (constraint)
  for(gamma in gammas){
    # Vector containing the CAGRS of each asset (constraint)
    dvec = cagrs_mat*gamma
    qp.out = solve.QP(Dmat=sigma_mat, dvec=dvec, Amat=Amat, bvec=bvec, meq = meq)
    # Solution weights
    weights =qp.out$solution
    # Get the return and volatility
    ret_portfolio = cagrs_mat %*% weights
    vol_portfolio = sqrt(t(weights) %*% (sigma_mat %*% weights))
    sharpe = ret_portfolio/vol_portfolio
    # Adding to list of optimized return to volatility
    list_Ret_QP = c(list_Ret_QP, as.numeric(ret_portfolio))
    list_Vol_QP = c(list_Vol_QP, as.numeric(vol_portfolio))
    # Used to pick the right gamma, actually, 
    list_Sharp_QP = c(list_Sharp_QP, as.numeric(sharpe))
  }
  # return a data frame witht the volatility and associated return
  return (data.frame(list_Ret_QP, list_Vol_QP, list_Sharp_QP,gammas))
}

# Markowtiz Allocation, gamma, the stocks, rebalancing window and regulation factor
Markowitz_optimisation <- function(gamma, selected_stock, rebalancing, rho){
  # Calculate the returns 
  returns_matrix <-  apply(selected_stock, 2, FUN = Compute_return)
  # Number of observations / number of stock we want to invest in 
  num_Observation = dim(returns_matrix)[1]
  num_Asset = dim(returns_matrix)[2]
  # Matrix that will contain the weights by Markowitz optimisation
  weights = matrix(nrow = num_Observation, ncol = num_Asset)
  # Equally weighted at the beginning
  weights_before_op = rep(1.0/num_Asset, num_Asset)
  # Rebalancing index
  rebalancing_index_day = seq(253, num_Observation - rebalancing, by = rebalancing)
  
  # Optimization with rebalacing the portfolio every period
  for(indx in rebalancing_index_day){
    # Values based on 1 year
    values_temp = selected_stock[(indx-252):indx, ]
    return_temp = returns_matrix[(indx-252):indx, ]
    # Calculate CAGR and sigma, we should base the investment period on at least 1 year
    CAGRs_temp = apply(values_temp, 2, FUN = Compute_CAGR)
    sigma_temp = cov(return_temp) * 252
    # Quadratic programming
    Dmat = sigma_temp + 2*rho*diag(num_Asset)
    dvec = (gamma*CAGRs_temp + 2*rho*weights_before_op)
    Amat = cbind(rep(1, num_Asset), diag(num_Asset))
    bvec = c(1, rep(0, num_Asset))
    meq = 1
    # get te solution
    qp.out = solve.QP(Dmat=Dmat, dvec=dvec, Amat=Amat, bvec=bvec, meq = meq)
    weights[(indx+1):(indx+rebalancing), ] = matrix(rep(qp.out$solution, rebalancing), ncol = num_Asset,  byrow = TRUE)
    weights_before_op = qp.out$solution
  }
  weights = weights[rowSums(is.na(weights)) != ncol(weights), ]
  return (data.frame(rebalancing_index_day, weights))
}

# We assume equally weighted
Get_data_equally_weighted_portfolio <- function(data_portfolio){
  return (data.frame(rowSums(data_portfolio)))
}

# Calculate the cumulative returns of a specific portfolio
Get_cum_ret <- function(return_port) {
  return(cumprod(1 + return_port)-1)
}

# Selection of the stocks according to certain criteria
Select_best_stocks <- function(number_wanted, stock_price) {
  # Name of the selected stocks 
  selected = c()
  number_of_selection = minimum(number_wanted, ncol(stock_price))
  if(number_of_selection == ncol(stock_price)) {
    return(colnames(stock_price))
  }
  # Select the best stocks by cummulative returns
  cummulative_returns <- Get_cum_ret(as.data.frame(apply(stock_price, 2, FUN = Compute_return)))    
  # Sorting by cummulative returns
  cummulative_returns <- cummulative_returns[,order(cummulative_returns[nrow(cummulative_returns),],decreasing = TRUE)]
  remaining = c(colnames(cummulative_returns))
  # The more diversifiedn the better
  correlation_matrix = as.data.frame(cor(cummulative_returns))
  # Add the first stock, assumed it's the best one
  selected = c(selected, colnames(cummulative_returns[1]))
  remaining = remaining[!remaining %in% selected]
  index_last_added = 1
  # Afterwards, we try to take the ones with the highest cummulated returns with a low correlation to the others
  # Diversification and high return strategy
  if (number_of_selection > 1) {
    for (i in 2:number_of_selection) {
      for (j in index_last_added+1:ncol(cummulative_returns)) {
        average_corr = 0
        for (k in 1:length(selected)) {
          average_corr = average_corr + correlation_matrix[rownames(correlation_matrix) == selected[k], colnames(correlation_matrix) == colnames(cummulative_returns[j])]
        }
        average_corr = average_corr/length(selected)
        # Correlation lower than 0.4 then add to portfolio
        if (average_corr < 0.4) {
          selected = c(selected, colnames(cummulative_returns[j]))
          remaining = remaining[!remaining %in% selected]
          index_last_added = j
          # Out of the loop
          {break}
        }
      }
    }
  }
  # In case the number of stock selected by the loop are not sufficent, add the stocks with te best cummulative return
  if (length(selected) < number_of_selection) {
    # Add remaining
    to_add = number_of_selection - length(selected)
    for (i in 1:to_add) {
      selected = c(selected, remaining[i])
    }
    remaining = remaining[!remaining %in% selected]
  }
  return(selected)
}

# Average expected excess return 
Compute_avg_excess_return <- function(Compute_ab_return) {
  return(colMeans(Compute_ab_return))
}

# Maximum Drawdown / Maximum Drawdown Duration
Compute_MDD_portfolio <- function(portfolio) {
  # Beginning values
  mdd = 0
  beg = 0
  end = 0
  mdd_temp = 0
  for (i in 2:nrow(portfolio)) {
    for (j in 2:i) {
      # Calculate the drawdown
      mdd_temp = (portfolio[i,1] - portfolio[j,1])/portfolio[j,1]
      if (mdd_temp < mdd) {
        # Drawdown higher
        mdd <- mdd_temp
        beg = j
        end = i
      }
    }
  }
  mdd = -mdd
  return (c(mdd, beg, end))
}

# Value at risk with gaussian distribution assuming the expected return is 0
Compute_gaussian_VaR_95 <- function(std){
  return(1.65*std)
}

Compute_gaussian_VaR_99 <- function(std){
  return(2.33*std)
}

# To get the data from a specific industry we can use Get_data_industry("INSERT INDUSTRY", industries, data)
# For several industrieswe use Get_data_industry(c(ind1, ind2...), industries, data)
# test_ind = Get_data_industry(c("Technology", "Energy"), industries, stock_price)

## Separting the values by industries
# Dates 
dates <- stock_price %>% select(c(Date))

# All except Index
ind_all <- stock_price %>% select(c(-.CACT, -.FCHI))

# Indexes
ind_index <- stock_price %>% select(c(.CACT, .FCHI))

# Results
returns = as.data.frame(apply(stock_price[,-1], 2, FUN = Compute_return))
CAGRs = as.data.frame(t(apply(stock_price[,-1], 2, FUN = Compute_CAGR)))
ab_returns = as.data.frame(apply(returns, 2, FUN = Compute_ab_return)) 
cum_ret = as.data.frame(apply(returns, 2, FUN = Get_cum_ret))
# TEST : equally_weighted_tech =  Get_data_equally_weighted_portfolio(Get_data_industry("Technology",industries,stock_price))
# TEST : sortino = as.data.frame(t(apply(returns - 0.001, 2, FUN = Compute_sortino)))
# TEST : sharpe = as.data.frame(t(apply(returns, 2, FUN = Compute_sharpe)))
# TEST : sigmas = Compute_sigma(returns)
# TEST : portfolios = Get_data_fictive_portfolio(Get_data_industry("Technology",industries,stock_price))
# TEST : portfolio_sol = Get_Markowitz_frontier(Get_data_industry("Technology",industries,stock_price))
# position = which.min(portfolio_sol[,2] < 0.285)
# acceptable = portfolio_sol[1:position,]
# maximum_sharp_index = which.max(acceptable[,3])
# associated_gamma = portfolio_sol[maximum_sharp_index,4]
# TEST : weight_markowitz = Markowitz_optimisation(associated_gamma,Get_data_industry("Technology",industries,stock_price),10,0.02)
# test = cbind(weight_markowitz, dates[weight_markowitz[,1],1])
# colnames(test) <- c('Rebalancing', colnames(Get_data_industry("Technology",industries,stock_price)),'date')
# TEST : t1 = Select_best_stocks(5, stock_price[,-1])
# TEST : MDD_CAC = Compute_MDD_portfolio(as.data.frame(stock_price[,ncol(stock_price)]))
# Maximum_duration = MDD_CAC[3] - MDD_CAC[2] 


# Layouts
# User interface
header <- dashboardHeader(title = "Investment Advisory")
# Side bar with basically the inputs and the different tabs
sidebar <- dashboardSidebar(
          # Tabs where the user can navigate
          sidebarMenu(
            menuItem("About the project", tabName = "about", icon = icon("question")),
            menuItem("Overview", tabName = "overv", icon = icon("search")),
            menuItem("Risk measures and returns", tabName = "results", icon = icon("table")),
            menuItem("Dashboard Markowitz", tabName = "dashboard", icon = icon("dashboard")),
            # Inputs
            menuItem("Inputs", icon = icon("input-cursor"),
              checkboxGroupInput(inputId="industry_select", label = "Select industries", industry_list$Industry,
                                 selected = industry_list$Industry),
              numericInput("acc_return", "Minimum acceptable return ", 0, min = -1, max = 1, step=0.001), 
              numericInput("acc_volatility", "Maximum acceptable volatility for Markowitz ", 0.28, min = 0, max = 1, step=0.001), 
              numericInput("stock_bn", "How many stocks would you like to invest in ?", 3, min = 2, max = 10, step=1),
              sliderInput("rebalance_days", "Rebalancing frequency for Markowitz", 2, min = 0, max = 300, step=1),
              dateRangeInput(inputId="date_select", label="Investment period (at least 1 year)", start = "2018-06-12", end = "2022-05-05", min = "2018-06-12", max = "2022-05-05", format = "yyyy-mm-dd", language="en"),
              numericInput("rho_in", "Regulation coefficient (rho) :", 0.01, min = 0, max = 100, step=0.001)
            )
          )
      )

body <- dashboardBody(
  mainPanel(width = 12,
    # Information about the project
    tabItems(
      tabItem(tabName = "about", tabBox(
      title = "Project and Author",
      tabPanel(title = "Purpose", HTML("<H3>Financial Project</H3>  
              <p>The <a href = https://evasyutyicheng.shinyapps.io/Project/>project</a> has been created as part of the Financial Visualization course of the Data Management Master taught by Professor Alexandre Garel. The purpose of this application is to provide more informations about the French Stock Market (CAC all).
               </p><p>The application calculate risk measures such as the Gaussian Value at Risk, the volatility, the abnormal returns. But also some to tell whether the investment is worth it thanks to Sharpe ratio and Sortino ratio. </p>
               <p>There is also a part dedicatied to the theory of Markowitz and the efficent portfolio allocation regarding returns/risks. The tools allow the user to input some variables, especially selecting an industry or choosing the number of stock they want. </p><p>For the good usage of the tool, it 
               is important to select an investment period which is above 1 year. At the begining, we will use an equally weighted portfolio, this is not a rigourous way to represent an industry since it depends on its cap. Calculations may take time especially for the risk part, please be patient. </p>")),
      tabPanel(title = "About me", HTML("<H3>Eva Cheng - DMF 2021/2022</H3><p>Student in the DMF program, studied Computer Science and IT Engineering with Master majoring in Market Finance.</p><p>I have chosen to relate the project to finance and particularly with my other major (market finance) to do something really meaningful for my career. </p><p>Visit my
                                      <a href = https://github.com/Syutyi>GitHub page</a> and add me on LinkedIn : <a href = https://www.linkedin.com/in/eva-syutyi-cheng/>LinkedIn</a>.  You can also contact me by <a href=mailto:syutyi.cheng@gmail.com?subject=Contact>Email</a>.</p>")),
      tabPanel(title = "Return measures and Sharpe", HTML("<p>Average return and Average excess return : These measures indicate how much return the investor is expected to get every day and how much return in excess compared to the CAC all tradable (the market).</p>
      <p>CAGR (Compound annual growth rate) : It is the yearly or annual growth rate estimated for this period of investment (longer than 1 year usually, this is why I asked the user to input a period longer than 1 year).</p>
      <p>Total cummulated return : It is the cummulated return of the last day of investment, considering the period of investment chosen by the user.</p>
      <p>Sharpe ratio : A measure of risk adjusted return, It tells about how much excess return the investor is going get while investing in riskier assets</p>")),
      tabPanel(title = "Sortino ratio and risk measures", HTML("<p>Sortino ratio (from the Economics Times) : It measures measures the performance of the investment relative to the downward deviation. Unlike Sharpe, it doesn't take into account the total volatility in the investment but it considers a user-specified target.</p>
      <p>Maximum drawdown : It represents the biggest decline from a peak in historical data, and the duration tells how many days the decline has lasted.</p>
      <p>Value at risk : It measures the greatest risk of loss of an investment with a certain probability level. For example a VaR99 that is 6% can be interpreted as, the investor as probability of 99% to have a loss that is at most equal to 6%.</p>
      <p>Volatility : It measures the degree of variation of the return (in positive but also in negative), it can be interpreted as the amplitude of variation.</p>"))
      ),
      tabBox(
        title = "How was it built ?",
        tabPanel("Data", HTML("<H3>Resources</H3><p>The data got collected from Refinitiv EIKON taking in account the CAC all tradable, the closing price are only the one which have been kept and before starting the project, the dataset has been cleaned and the NULL values have been removed. Therefore, we only have data starting 2018 
                 since some stocks or index (CAC40 or CACT) did not exist before the first day contained in the dataset.</p> 
                 <p>The industries have also been collected from Refinitiv EIKON, these kind of data are hard to find on Kaggle or in other website that do not relate to finance. </p>
                 The dataset has been collected from and EIKON and cleaned. It is available <a href =  https://github.com/Syutyi/2022-M2-Asset-Management/blob/main/Price_per-stock.csv> here</a>. Similarly to the stock price data, the dataset is available <a href =  https://github.com/Syutyi/2022-M2-Asset-Management/blob/main/Stocks_CAC.csv>here</a>. 
                  </p><p>The class material was provided by Professor Garel and it covers how to do an application in R and data visualization in R language. Some references are also available on shiny.rstudio.com, the website has been used a lot to complete the knowledge brought by the class.</p>
                              <p>The project is published online following the link <a href = https://evasyutyicheng.shinyapps.io/Project/>here</a>, but also on my gitHub page, repository named 2022-M2-Asset_Management_Visualization.")),
        tabPanel("Skills & Techniques", HTML("<H3>Capabilities</H3><p>The project has been done with R programming and the course of Professor Garel covers most of the technical specificities of the R language. </p><p>The financial notions are the one that are mainly used in quantitative market finance. I did not really use simulation methods such as Monte
                 Carlo except when I tried to simulate the CAPM of Markowitz. </p><p>The risk measures calculated are the most commonly used in finance : Gaussian Value at Risk, volatility, Maximum drawdown which is the maximum loss during the period and its duration (also called MDD).
                 </p><p>DISCLAIMER : the project is not about prediction on the CAC but more about visualization of the past perfomance and risk measures that may be used in the future.</p>
                  </p><p>The data is only based on stock prices from June 2018 to May 2022 therefore, the situation may have changed since the collection of the data and EIKON do not provide all the stock prices of the CAC all tradable, we only have 60 stocks.</p>
                                             <p>The project is published online following the link <a href = https://evasyutyicheng.shinyapps.io/Project/>here</a>."),), 
        # Hide the double scroll bar
        tags$head(tags$style(
          HTML('.wrapper {height: auto !important; position:relative; overflow-x:hidden; overflow-y:hidden}'),
          tags$head(tags$style(HTML('.skin-blue .left-side, .skin-blue .wrapper {background-color: #ecf0f5;}')))
        ))
      ),
      tabsetPanel(
        id = 'dataset',
        tabPanel("Industries", DT::dataTableOutput("table_ind"), width = "100%"),
        tabPanel("Stock prices", DT::dataTableOutput("table_stock"), width = "100%")
      )
    ),
    tabItem(tabName = "overv", 
            # Get the input of the user
            valueBoxOutput("in_acc_return"), 
            valueBoxOutput("in_acc_volatility"),
            valueBoxOutput("select_ind"),
            valueBoxOutput("start_date"),
            valueBoxOutput("in_stocks"), 
            valueBoxOutput("end_date"),
            box(title = "Cummulative return of selected stocks", width = 8, solidHeader = TRUE, status = "primary",
                plotOutput("cum_plot_2", click = "plot_click")
            ),
            box(
              title = "Explanation of the over view", width = 4, background = "light-blue",
              HTML("<p>The idea of this tab is to give an overview of the user's input and the portfolio generated by the machine based on the user's choice except for the volatility that has not been taken in account for the selection unfortunately because the idea was to diversify the portfolio as much as possible.</p> 
                 <p>From what we can observe, we can say that the generated portfolio's cumulated return is performing better than the market, the industry and even the indexes related to these.</p>
                 <p>The reason is that I am considering the stocks with the highest cummulated returns but also the ones that are not really correlated in other to diversify the portfolio as much as possible. Diversification is important since we do not want to put all the eggs in the same basket.
                  </p><p>Of course, the overview here is nothing about prediction and as a disclaimer, I want to say that this is not the only portfolio one should invest in since it is based on past performance and well-performing stocks can because stocks with bad performance. Implementing some machine learning may be useful 
                   in that case. </p>
                   <p>The portfolio is equally weighted for this example, later in the part with Markowitz, I will compare equally weighted to optimized portfolio."))
            ,
            tabsetPanel(
              id = 'test1',
              tabPanel("Portfolio return", DT::dataTableOutput("table1"), width = "100%")
            )),
    tabItem(tabName = "results",
            # Risk measures and performance
            valueBoxOutput("avg_ret"), 
            valueBoxOutput("excess_ret"),
            valueBoxOutput("cagr"), 
            valueBoxOutput("cum_ret_port"),
            valueBoxOutput("vol_port"),
            valueBoxOutput("sortino_port"),
            valueBoxOutput("sharpe_port"),
            valueBoxOutput("mdd_port"),
            valueBoxOutput("mdd_dur_port"),
            valueBoxOutput("var_95"),
            valueBoxOutput("var_99"),
            valueBoxOutput("best_perf"),
            box(title = "Cummulative return of selected stocks", width = 7, solidHeader = TRUE, status = "primary",
                plotOutput("cum_plot", click = "plot_click")
            ),
            box(title = "Histogram of the returns", width = 5, solidHeader = TRUE, status = "primary",
                plotOutput("dist_return")
            )
           ),
    tabItem(tabName = "dashboard",
            box(title = "Efficient frontier between return and volatility", width = 6, solidHeader = TRUE, status = "primary",
                plotOutput("mark_2", click = "plot_click"),
                verbatimTextOutput("markowitz_coords_2")
            ),
            box(title = "Markowitz optimized portfolio", width = 6, solidHeader = TRUE, status = "primary",
                plotOutput("data_opti",click = "plot_click"),
                verbatimTextOutput("markowitz_op")
            ),
            box(title = "Cumulative returns of the two portfolio", width = 12, solidHeader = TRUE, status = "primary",
                plotOutput("cum_ret_mark")
            ),
            tabsetPanel(
              id = 'test1',
              tabPanel("Portfolio allocation", DT::dataTableOutput("table_test"), width = "100%")
            ))
  )
)
)


ui <- dashboardPage(
  skin = "blue",
  header, sidebar,body
)


# Calculation zone and server
server <- function(input, output, session) {
  # Output the user inputs
  output$table_ind = DT::renderDataTable(
    industries,
    options = list(scrollX = TRUE)
  )
  output$table_stock = DT::renderDataTable(
    stock_price,
    options = list(scrollX = TRUE)
  )
  output$in_acc_return <- renderValueBox({
    valueBox(
      paste0(input$acc_return*100, "%"), "Min Return", icon = icon("repeat"),
      color = "blue"
    )
  })
  output$in_acc_volatility <- renderValueBox({
    valueBox(
      paste0(input$acc_volatility*100, "%"), "Max Volatility", icon = icon("warning"),
      color = "blue", width = 10
    )
  })
  output$in_stocks <- renderValueBox({
    valueBox(
      input$stock_bn, "Number of assets", icon = icon("piggy-bank"),
      color = "blue"
    )
  })
  output$select_ind <- renderValueBox({
    valueBox(
      length(input$industry_select), "Industries selected", icon = icon("calendar"),
      color = "blue"
    )
  })
  output$start_date <- renderValueBox({
    valueBox(
      "Start", input$date_select[1], icon = icon("flag"),
      color = "green"
    )
  })
  output$end_date <- renderValueBox({
    valueBox(
      "End", input$date_select[2], icon = icon("stop"),
      color = "maroon"
    )
  })
  
  # Calculation zone for the returns 
  # Select the stocks of the selected industry
  user_industry_list = reactive({
    if (is.null(input$industry_select)){
      return(NULL)
    }
    data <- Get_data_industry(input$industry_select, industries, stock_price)
    return (data)})
  # Adding the dates
  stock = reactive(cbind(dates, user_industry_list()))
  # Eliminating the rows that are not in the date range
  stock_in_range = reactive(stock()[stock()$Date >= input$date_select[1] & stock()$Date <= input$date_select[2],])
  # Selection of the stocks 
  stock_selected = reactive(Select_best_stocks(input$stock_bn, stock_in_range()[,-1])) 
  # filtering the dataset with the selected stock
  final_stock_data <- reactive(stock_in_range()[,names(stock_in_range()) %in% stock_selected()])
  #dataset_industry_user = Get_data_industry(input$industry_select, industries, stock_price)
  portfolio_user <- reactive({data = Get_data_equally_weighted_portfolio(final_stock_data())
                            names(data)[1] = "Price_portfolio"
                            return(data)})
  return_portfolio <- reactive({data = as.data.frame(Compute_return(portfolio_user()[,1]))
                        names(data)[1] = "Return_portfolio"
                        return(data)})
  user_cumulative_return <- reactive({data = as.data.frame(Get_cum_ret(return_portfolio()[,1]))
                            names(data)[1] = "Cum_ret_portfolio"
                            return(data)})
  # Measures
  portfolio_avg_excess_return = reactive(mean(Compute_ab_return(return_portfolio()[,1])))
  total_cum_ret <- reactive(user_cumulative_return()[nrow(user_cumulative_return()),1])
  CAGR <-  reactive(Compute_CAGR(portfolio_user()[,1]))
  volatility_port <- reactive(sd(return_portfolio()[,1]))
  max_dd <- reactive(Compute_MDD_portfolio(as.data.frame(portfolio_user()[,1])))
  duration_mdd = reactive(max_dd()[3] - max_dd()[2])
  sharpe_port = reactive(Compute_sharpe(return_portfolio()[,1]))
  sortino_port = reactive(Compute_sortino(return_portfolio()[,1] - input$acc_return))
  avg_ret = reactive(mean(return_portfolio()[,1]))
  VaR_1 = reactive(Compute_gaussian_VaR_95(volatility_port()[1]))
  VaR_2 = reactive(Compute_gaussian_VaR_99(volatility_port()[1]))
  
  # Value box for the calculations
  output$avg_ret <- renderValueBox({
    color <- 'blue'
    if(avg_ret()[1] < 0) color <- 'red'
    if(avg_ret()[1] > 0.05) color <- 'green'
    valueBox(
      paste0(round(avg_ret()[1]*100,3), "%"), "Average return", icon = icon("repeat"),
      color = color, width = 10
    )
  })
  output$excess_ret <- renderValueBox({
    color <- 'blue'
    if(portfolio_avg_excess_return()[1] < 0) color <- 'red'
    if(portfolio_avg_excess_return()[1] > 0.03) color <- 'green'
    valueBox(
      paste0(round(portfolio_avg_excess_return()[1]*100,3), "%"), "Excess return", icon = icon("plus"),
      color = color, width = 10
    )
  })
  output$cagr <- renderValueBox({
    color <- 'blue'
    if(CAGR()[1] < 0) color <- 'red'
    if(CAGR()[1] > 0.1) color <- 'green'
    valueBox(
      paste0(round(CAGR()[1]*100,3), "%"), "CAGR", icon = icon("arrow-up"),
      color = color, width = 10
    )
  })
  output$cum_ret_port <- renderValueBox({
    color <- 'blue'
    if(total_cum_ret()[1] < 0) color <- 'red'
    if(total_cum_ret()[1] > 0.5) color <- 'green'
    valueBox(
      paste0(round(total_cum_ret()[1]*100,2), "%"), "Cum return", icon = icon("piggy-bank"),
      color = color, width = 10
    )
  })
  output$vol_port <- renderValueBox({
    color <- 'blue'
    if(volatility_port()[1] > 0.3) color <- 'red'
    if(volatility_port()[1] < 0.05) color <- 'green'
    valueBox(
      paste0(round(volatility_port()[1]*100,3), "%"), "Volatility", icon = icon("dove"),
      color = color, width = 10
    )
  })
  output$sortino_port <- renderValueBox({
    color <- 'blue'
    if(sortino_port()[1] < 0.2) color <- 'red'
    if(sortino_port()[1] < 0) color <- 'black'
    if(sortino_port()[1] > 0.5) color <- 'green'
    valueBox(
      paste0(round(sortino_port()[1]*100,3), "%"), "Sortino ratio", icon = icon("exclamation"),
      color = color, width = 10
    )
  })
  output$sharpe_port <- renderValueBox({
    color <- 'blue'
    if(sharpe_port()[1] < 0.2) color <- 'red'
    if(sharpe_port()[1] > 0.5) color <- 'green'
    valueBox(
      paste0(round(sharpe_port()[1]*100,3), "%"), "Sharpe ratio", icon = icon("exclamation"),
      color = color, width = 10
    )
  })
  output$mdd_port <- renderValueBox({
    color <- 'blue'
    if(max_dd()[1] > 0.4) color <- 'red'
    if(max_dd()[1] < 0.1) color <- 'green'
    valueBox(
      paste0(round(max_dd()[1]*100,2), "%"), "Max DD", icon = icon("arrow-down"),
      color = color, width = 10
    )
  })
  output$mdd_dur_port <- renderValueBox({
    valueBox(
      duration_mdd(), "DD duration (days)", icon = icon("calendar"),
      color = "blue", width = 10
    )
  })
  output$var_95 <- renderValueBox({
    color <- 'blue'
    if(VaR_1() > 0.06) color <- 'red'
    if(VaR_1() < 0.02) color <- 'green'
    valueBox(
      paste0(round(VaR_1()*100,3), "%"), "Value at risk 95%", icon = icon("exclamation"),
      color = color, width = 10
    )
  })
  output$var_99 <- renderValueBox({
    color <- 'blue'
    if(VaR_2() > 0.06) color <- 'red'
    if(VaR_2() < 0.02) color <- 'green'
    valueBox(
      paste0(round(VaR_2()*100,3), "%"), "Value at risk 99%", icon = icon("exclamation"),
      color = color, width = 10
    )
  })
  output$best_perf <- renderValueBox({
    valueBox(
      stock_selected()[1], "Best performing", icon = icon("star"),
      color = "blue", width = 10
    )
  })
  output$selection <- renderDataTable(
    final_stock_data(),
    options = list(scroller = TRUE, scrollX = TRUE),
    fillContainer = TRUE
  )
  
  # Cumulative of the returns
  
  data_portfolio_stock = reactive((as.data.frame(apply(as.data.frame(apply(final_stock_data(), 2, FUN = Compute_return)), 2, FUN = Get_cum_ret))))
  data_port_and_stock = reactive({data = cbind(stock_in_range()[,1],data_portfolio_stock(),user_cumulative_return())
                                  names(data)[1] = "Date"
                                  return (data)})
  data_long <- reactive(melt(data_port_and_stock(), id = "Date"))
  output$cum_plot <- renderPlot({
    ggplot(data_long(),            
           aes(x = Date,
               y = value,
               color = variable)) +  geom_line()
  })
  # Create the histogram for the returns
  output$dist_return <- renderPlot({
    hist(return_portfolio())
  })
  
  # Create several dataset to represent the overview : Cumulative return of the industry, the portfolio, the CACall, CAC40, the equally weighted cac
  # Industry
  equally_weighted_porfolio_industry = reactive(Get_data_equally_weighted_portfolio(stock_in_range()[,-1]))
  return_industry <- reactive({data = as.data.frame(Compute_return(equally_weighted_porfolio_industry()[,1]))
  names(data)[1] = "Return_industry"
  return(data)})
  industry_cumulative_return <- reactive({data = as.data.frame(Get_cum_ret(return_industry()[,1]))
  names(data)[1] = "Cum_ret_industry"
  return(data)})
  # Portfolio already done
  # CACall
  cac_all_price = stock_price[,c(1,(ncol(stock_price)-1))]
  cac_all_price_range = reactive(cac_all_price[cac_all_price$Date >= input$date_select[1] & cac_all_price$Date <= input$date_select[2],])
  cac_all_retrun <- reactive({data = as.data.frame(Compute_return(cac_all_price_range()[,2]))
  names(data)[1] = "CAC_return"
  return(data)})
  cac_all_cumulative_return <- reactive({data = as.data.frame(Get_cum_ret(cac_all_retrun()[,1]))
  names(data)[1] = "CAC_cum_ret"
  return(data)})
  # CAC 40
  cac_40_price = stock_price[,c(1,(ncol(stock_price)))]
  cac_40_price_range = reactive(cac_40_price[cac_40_price$Date >= input$date_select[1] & cac_40_price$Date <= input$date_select[2],])
  cac_40_retrun <- reactive({data = as.data.frame(Compute_return(cac_40_price_range()[,2]))
  names(data)[1] = "CAC40_return"
  return(data)})
  cac_40_cumulative_return <- reactive({data = as.data.frame(Get_cum_ret(cac_40_retrun()[,1]))
  names(data)[1] = "CAC40_cum_ret"
  return(data)})
  # Equally weighted CAC 
  all_data = stock_price[,c(-(ncol(stock_price)-1), -ncol(stock_price))]
  all_price_range = reactive(all_data[all_data$Date >= input$date_select[1] & all_data$Date <= input$date_select[2],])
  equally_weighted_all = reactive(Get_data_equally_weighted_portfolio(all_price_range()[,-1]))
  all_return <- reactive({data = as.data.frame(Compute_return(equally_weighted_all()[,1]))
                names(data)[1] = "Equal_return"
                return(data)})
  cac_all_cum_ret <- reactive({data = as.data.frame(Get_cum_ret(all_return()[,1]))
  names(data)[1] = "Equal_cac_cum_ret"
  return(data)})
  all_cum = reactive({data = cbind(stock_in_range()[,1], user_cumulative_return(), industry_cumulative_return(), cac_all_cumulative_return(), cac_40_cumulative_return(),
                            cac_all_cum_ret())
                      names(data)[1] = "Date"
                      return (data)})
  output$table1 = DT::renderDataTable(
    datatable(data_port_and_stock()),
    options = list(scrollX = TRUE),
    fillContainer = TRUE
  )
  # Create a plot
  data_overview <- reactive(melt(all_cum(), id = "Date"))
  output$cum_plot_2 <- renderPlot({
    ggplot(data_overview(),            
           aes(x = Date,
               y = value,
               color = variable)) +  geom_line()
  })
  
  # Markowitz part
  # Simulation, Monte Carlo
  markowitz_fictive_port <- reactive(Get_data_fictive_portfolio(final_stock_data()))
  # Efficient frontier
  markowitz_frontier <- reactive(Get_Markowitz_frontier(final_stock_data()))
  # Get the gamma for the optimization depending on the level of risk 
  position = reactive(which.min(markowitz_frontier()[,2] < input$acc_volatility))
  acceptable = reactive(markowitz_frontier()[1:position(),])
  # Get the maximum sharpe ratio based on this level of risk
  maximum_sharp_index = reactive(which.max(acceptable()[,3]))
  associated_gamma = reactive(markowitz_frontier()[maximum_sharp_index(),4])
  # Optimize
  optimized = reactive(Markowitz_optimisation(associated_gamma(),final_stock_data(), input$rebalance_days, input$rho_in))
  optimized_dates = reactive({dates[c(optimized()[,1]),]})
  optimized_data = reactive({data = cbind(optimized_dates(), optimized())
                colnames(data)[1] = "Date"
                return(data)})
  # Assume we invest in 100 stocks
  volume_opti = reactive({data = cbind(optimized_data()[,1], optimized_data()[,3:ncol(optimized_data())])
                          colnames(data)[1] = "Date"
                          colnames(data)[2:ncol(data)] = colnames(final_stock_data())
                          return(data)})
  volume_not_opti = reactive(1/input$stock_bn)
  price_at_date =reactive(final_stock_data()[c(optimized()[,1]),])
  new_portfolio = reactive(price_at_date()*round(volume_opti()[,-1]*100))
  new_portfolio_with_dates = reactive({data = cbind(optimized_dates(), new_portfolio())
                            colnames(data)[1] = "Date"
                            return (data)})
  other_portfolio = reactive({data = cbind(optimized_dates(),price_at_date()*round(volume_not_opti()*100))
                      colnames(data)[1] = "Date"
                      return (data)})
  # Calculate the cummulated return of the portfolio
  # do the sum
  price_new_pf = reactive(Get_data_equally_weighted_portfolio(new_portfolio_with_dates()[,-1]))
  price_old_pf = reactive(Get_data_equally_weighted_portfolio(other_portfolio()[,-1]))
  # Returns 
  new_pf_return <- reactive({data = as.data.frame(Compute_return(price_new_pf()[,1]))
                                                  names(data)[1] = "new_pf_ret"
                                                  return(data)})
  old_pf_return <- reactive({data = as.data.frame(Compute_return(price_old_pf()[,1]))
                                                  names(data)[1] = "old_pf_ret"
                                                  return(data)})
  # Cumulated returns
  new_pf_cum_ret <- reactive({data = as.data.frame(Get_cum_ret(new_pf_return()[,1]))
  names(data)[1] = "new_pf_cum_ret"
  return(data)})
  old_pf_cum_ret <- reactive({data = as.data.frame(Get_cum_ret(old_pf_return()[,1]))
  names(data)[1] = "old_pf_cum_ret"
  return(data)})
  
  # Markowitz plot
  output$mark_1 <- renderPlot({
    plot(markowitz_fictive_port()$list_Vol, markowitz_fictive_port()$list_Ret, type = 'p', pch= 16, cex = 0.7, xlab = 'Volatility', ylab = 'Return', xlim = c(0.1, 0.4), ylim = c(-0.05, 0.3))
  })
  output$markowitz_coords_1 <- renderText({
    paste0("Volatility = ", input$plot_click$x, "\nReturn = ", input$plot_click$y)
  })
  output$mark_2 <- renderPlot({
    plot(markowitz_fictive_port()$list_Vol, markowitz_fictive_port()$list_Ret, type = 'p', pch= 16, cex = 0.7, xlab = 'Volatility', ylab = 'Return', xlim = c(0.18, 0.34), ylim = c(-0.05, 0.3))
    lines(markowitz_frontier()$list_Vol_QP, markowitz_frontier()$list_Ret_QP, col = 'red', lwd=4)
  })
  output$markowitz_coords_2 <- renderText({
    paste0("Volatility = ", input$plot_click$x, "\nReturn = ", input$plot_click$y)
  })
  output$markowitz_op <- renderText({
    paste0("Date = ", input$plot_click$x, "\nInvestment = ", input$plot_click$y)
  })
  
  # Cummulative returns comparison
  output$cum_ret_mark <- renderPlot({
    plot(new_pf_cum_ret()[,1], type = 'l', xlab = 'Investment days', ylab = 'Cumulated returns', pch= 16)
    lines(old_pf_cum_ret()[,1], col = 'red', lwd=4, pch= 16)
    legend("bottomleft", legend=c("Markowitz optimized", "Old method"),col=c("black", "red"), lty=1:2, cex=0.8)
  })
  
  # Portfolio allocation
  data_opti <- reactive(melt(volume_opti(), id = "Date"))
  output$data_opti <- renderPlot({
    ggplot(data = data_opti(),            
           aes(x = Date,
               y = value,
               fill = variable)) +  geom_area(position="fill") + ylim(0,1)
  })
  
  output$table_test = DT::renderDataTable(
    datatable(volume_opti()[,-1]),
    options = list(scroller = TRUE, scrollX = TRUE),
    fillContainer = TRUE
  )
}
shinyApp(ui, server)
