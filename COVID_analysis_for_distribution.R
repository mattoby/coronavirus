
####################
#### File names ####
####################

warning("

pull the github repo from Johns Hopkins COVID-19: 
        
https://github.com/CSSEGISandData/COVID-19 and redo filenames for where you have the repo
        
everything else should run fine then.")

fname_COVID_confirmed <- 'GitHub/COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Confirmed.csv'
fname_COVID_deaths    <- 'GitHub/COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Deaths.csv'
fname_COVID_recovered <- 'GitHub/COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Recovered.csv'


fname_Hubei_confirmed_actual <- 'case_numbers_perday_Hubei_03142020.txt' 
fname_Wuhan_confirmed        <- 'case_numbers_perday_Wuhan_03152020.txt' 

#################
#### Imports ####
#################

#suppressPackageStartupMessages
library(devtools)
library(tidyverse)
library(zoo)
library(lubridate)
library(reshape2)
library(arm)
library(caret)
library(ordinal) 
library(car)
library(scales)
library(RColorBrewer)
library(assertthat)

set.seed(123456)

###################
#### functions ####
###################

convert_COVID_colname_to_date <- function(d){
  d1 <- substring(d, 2, nchar(d))
  d2 <- strsplit(d1, '\\.')[[1]]
  d3 <- as.Date(ISOdate(as.numeric(paste0("20", d2[3])), as.numeric(d2[1]), as.numeric(d2[2])))
  return(d3)
}

read_COVID_data <- function(fname, cases_col_name=NULL){
  dfc <- read.delim(file=fname         , header = TRUE, sep = "," , stringsAsFactors = FALSE) 
  current_columns <- colnames(dfc)
  
  ### deal with reshaping data so it's in EAV format, with city, date, and count
  row_ids <- current_columns[1:4]
  if (any(row_ids != c("Province.State", "Country.Region", "Lat", "Long"))){ 
    stop('check columns in input file')
  }
  dfc2 <- melt(dfc, id=row_ids) 
  
  ### deal with dates and base column names:
  dates_col <- as.character(dfc2$variable)
  dates_col <- (lapply(dates_col, convert_COVID_colname_to_date) %>% melt %>% dplyr::select(value))$value
  dfc2$dates <- dates_col
  dfc2$cases <- dfc2$value
  
  # clean up:
  dfc_out <- dfc2[c('Province.State','Country.Region','Lat','Long','dates', 'cases')]
  colnames(dfc_out) <- c('province','country','Lat','Long','dates', 'cases')
  dfc_out$Location <- paste0(dfc_out$province, "_", dfc_out$country)
  
  # rename the 'cases' column if desired (based on 2nd input to function):
  if (!is.null(cases_col_name)){
    names(dfc_out)[names(dfc_out) == 'cases'] <- cases_col_name
  }
  
  return(dfc_out)
}

# this reads in the data from https://jamanetwork.com/journals/jama/fullarticle/2762130
# and reposted in chart 7 of 
# https://medium.com/@tomaspueyo/coronavirus-act-today-or-people-will-die-f4d3d9cd99ca
read_Hubei_data <- function(fname){
  dfh <- read.delim(file=fname, header = TRUE, sep = "\t" , stringsAsFactors = FALSE) 
  dfh$dates <- as.Date(dfh$date, format="%d-%h-%y")
  dfh$date <- NULL
  dfh$log_confirmed <- log(dfh$confirmed)
  dfh$log_actual    <- log(dfh$actual   )
  return(dfh)
}


# this reads in the data from the WHO-China report Figure 2 first panel
# wuhan confirmed cases total outbreak (aside from first few days)
read_Wuhan_data <- function(fname){
  dfw <- read.delim(file=fname, header = TRUE, sep = "\t" , stringsAsFactors = FALSE) 
  dfw$dates <- as.Date(dfw$date, format="%d-%h-%y")
  dfw$date <- NULL
  dfw$confirmed_daily <- dfw$confirmed_wuhan
  dfw$confirmed_wuhan <- NULL

  return(dfw)
}
# organize the covid dfs for plotting 
# (eav style with attribute being confirmed, deaths, recovered)
combine_covid_dfs_for_plotting <- function(dfc, dfd, dfr){
  
  cat('Note: cases_col_names should be: confirmed, deaths, and recovered (we do a non-rigorous check)')
  stopifnot('confirmed' %in% names(dfc))
  stopifnot('deaths'    %in% names(dfd))
  stopifnot('recovered' %in% names(dfr))
  
  dfcplot <- dfc
  dfcplot$attribute <- 'confirmed'
  dfcplot$value <- dfcplot$confirmed
  dfcplot$confirmed <- NULL
  
  dfdplot <- dfd
  dfdplot$attribute <- 'deaths'
  dfdplot$value <- dfdplot$deaths
  dfdplot$deaths <- NULL
  
  dfrplot <- dfr
  dfrplot$attribute <- 'recovered'
  dfrplot$value <- dfrplot$recovered
  dfrplot$recovered <- NULL
  
  dfplot <- rbind(dfcplot, dfdplot, dfrplot)
  
  return(dfplot)
}

# validate that Location works as a primary key still 
# (same # unique rows as all key cols)
# this could potentially break in future with data uploads - 
# not sure what Hopkins data team is using as primary key
stopifnot_primary_key_valid <- function(df, 
                            primary_key_col='Location', 
                            all_key_cols=c('province', 
                                           'country', 
                                           'Lat', 
                                           'Long', 
                                           'Location'),
                            print_result=TRUE
                            ){
  num_unique__Location      <- nrow(unique(df %>% dplyr::select_(primary_key_col)))
  num_unique__total_records <- nrow(unique(df %>% dplyr::select_(.dots=as.list(all_key_cols))))
 
  stopifnot(num_unique__Location == num_unique__total_records)
  
  if (print_result){
    cat(paste0('Okay to use ', 
               primary_key_col, 
               ' as the primary key for ', 
               deparse(substitute(df)), 
               '\n'))
  }
}


remove_redundant_timepoints_per_Location <- function(df, remove_col){
  all_Locations <- unique(df$Location)
  
  df_new <- data.frame()
  for (i in seq(1, length(all_Locations))){
    loc = all_Locations[i]
    
    # grab Location and order by date
    df_curr <- df %>% 
      dplyr::filter(Location==loc)
    df_curr <- df_curr[order(df_curr$dates),]
    
    # stop if repeated dates, which shouldn't happen. this should be fixed 
    if (length(df_curr$dates)!=length(unique(df_curr$dates))){
      stop(paste0('you have repeated dates in the dataframe for this Location:\n\n', 
                  loc, 
                  '\n\nThis should be fixed upon data import.'))
    }
    
    # remove values from remove_col if they have the same value as the previous timepoint
    inds_to_keep = which(diff(df_curr[[remove_col]])!=0) + 1
    df_curr <- df_curr[inds_to_keep, ]
    
    # put df_curr back into df_new (the output)
    df_new <- rbind(df_new, df_curr)
    
  }
  
  return(df_new)
}



# Get params from a log-linear model fit
## assumes that your linear regression formula is: 
#
# log_Y ~ t
#
# (giving a log-linear fit to Y. note, t is the date)
# 
# this gives a model of the form: 
#
# log_Y = b0 + b1*t
#
# equiv to the form: 
#
# Y = exp(r * (t - t0 ) ) 
# 
# how interpret results:
#
# to get t0 (start point): t0 = -b0/b1
# to get r (growth rate) : r  =  b1
# to get doubling time   : tD =  log(2)/r
#
# This function assumes that the x-axis is called 'dates'
get_params_from_loglinear_model <- function(model){
  
  b0 = model$coefficients['(Intercept)']
  b1 = model$coefficients['dates']
  t0 = -(b0 / b1)
  r  = b1
  tD = log(2) / r
  
  names(t0) <- NULL
  names(r ) <- NULL
  names(tD) <- NULL
  

  # alternate forms of t0 (note, none of these lose information)
  t0_as_date <- as.Date(t0)
  t0_as_datetime <- as_datetime(t0_as_date)
  
  return(list(t0=t0,
              r=r,
              tD=tD,
              t0_as_date=t0_as_date,
              t0_as_datetime=t0_as_datetime))
}


##########################
#### global variables ####
##########################

shutdown_date_Wuhan = as.Date('2020-01-23')
shutdown_date_other = as.Date('2020-01-24')
last_day_cases_Wuhan = as.Date('2020-02-20') # from WHO report
last_day_cases_other = as.Date('2020-02-14') # from WHO report
italy_school_close_date = as.Date('2020-03-04')
italy_lockdown_date = as.Date('2020-03-09')

#################################################################
#################################################################

###########################
#### import COVID data ####
###########################


# data load:
dfc <- read_COVID_data(fname_COVID_confirmed, cases_col_name = 'confirmed')
dfd <- read_COVID_data(fname_COVID_deaths   , cases_col_name = 'deaths'   )
dfr <- read_COVID_data(fname_COVID_recovered, cases_col_name = 'recovered')
dfh <- read_Hubei_data(fname_Hubei_confirmed_actual)
dfw <- read_Wuhan_data(fname_Wuhan_confirmed       )

stopifnot_primary_key_valid(dfc)
stopifnot_primary_key_valid(dfd)
stopifnot_primary_key_valid(dfr)

#### combine data on confirmed, deaths, and recovered into one dataframe
dfcurr <- dfc
dfcurr <- merge(x=dfcurr, y=dfd, by=c('province', 'country', 'Location', 'Lat', 'Long', 'dates'), all.x=TRUE, all.y=TRUE)
dfcurr <- merge(x=dfcurr, y=dfr, by=c('province', 'country', 'Location', 'Lat', 'Long', 'dates'), all.x=TRUE, all.y=TRUE)
dfcurr$log_confirmed = log(dfcurr$confirmed)
dfcurr$log_deaths    = log(dfcurr$deaths   )
dfcurr$log_recovered = log(dfcurr$recovered)


stopifnot_primary_key_valid(dfcurr)

#### combine into dataframe with true entity-attribute-value for plotting
dfplot <- combine_covid_dfs_for_plotting(dfc, dfd, dfr)

stopifnot_primary_key_valid(dfplot)

#### important warning
warning("
 ####
 From now on, we will hardcode everything to use the 'Location' column as primary key! 
 ####
 
 That is, as primary key for locations (city/provinces), the original rows in the COVID files).
 
 Note, this is because writing code using standard evaluation is a headache. thanks R. 
 so, not specifying the column name of the primary key upfront leads to much more complex code.
 to simplify code for analysis, we assume that the primary key is Location, and if 
 that stops working on new data updates, we will fix the import functions to make 
 Location still unique on rows. 

")


#######################
#### summary stats ####
#######################

current_values <- dfcurr %>% filter(dates==max(dfcurr$dates)) %>% summarize(curr_confirmed=sum(confirmed, na.rm=TRUE),
                                                          curr_recovered=sum(recovered, na.rm=TRUE),
                                                          curr_deaths   =sum(deaths   , na.rm=TRUE))
print(current_values)

###############
#### plots ####
###############

# set colors
attribute_colors <- c(  confirmed="#D95F02", # red
                     deaths   ="#7570B3", # purple
                     recovered="#1B9E77" # green
)
# brewer.pal(length(unique(dfplot$attribute)),"Dark2")
attribute_colors_scale <- scale_colour_manual(name = '',values = attribute_colors)

# plot NYC & other cities:
cities = c('New York', 'Hubei')#, 'France') 
ggplot(filter(dfplot, province %in% cities), aes(x=dates, color=attribute, y=value)) + 
  geom_point() + 
  geom_line() + 
  facet_grid(province~.) +
  scale_y_log10() +
  attribute_colors_scale


# plot NYC & other cities:
cities = c('_Japan', '_Italy')
#cities = c('_Japan')
ggplot(filter(dfplot, Location %in% cities), aes(x=dates, color=attribute, y=value)) + 
  geom_point() + 
  geom_line() + 
  facet_grid(Location~.) +
  scale_y_log10() +
  attribute_colors_scale

# plot all Locations in descending order of confirmed cases:
do_all_plotting=FALSE
if (do_all_plotting){
  dfp2 <- dfplot %>% 
    filter(attribute=='confirmed') %>%
    group_by(Location) %>%
    summarize(max_confirmed = max(value, na.rm=TRUE)) %>%
    ungroup
  
  locs_by_descending_confirmed_numbers = dfp2[order(dfp2$max_confirmed, decreasing=TRUE),]$Location
  for (loc in locs_by_descending_confirmed_numbers){
    print(loc)
    p1 <- ggplot(filter(dfplot, Location == loc), aes(x=dates, color=attribute, y=value)) + 
      geom_point() + 
      geom_line() + 
      facet_grid(Location~.) +
      scale_y_log10() +
      attribute_colors_scale
    
    show(p1)
    readline(prompt="Press [enter] to continue")  
  }
              
}

# Special considerations:
# for determining unadulterated growth rate, all regions in china aside from (including?) 
# Hubei need to be considered only before the date cutoff when they put in restrictions. 
# Check when that date is vs. what we see in the country curves

# China plots:
dfplot_China <- dfplot %>% filter(country == 'China', attribute=='confirmed') 
dfplot_Hubei <- dfplot %>% filter(province == 'Hubei', attribute=='confirmed') 
dfplot_Others <- dfplot %>% filter(country %in% c('China', 'Italy', 'France', 'Korea, South', 'Iran'), attribute=='confirmed') 
dfplot_Others$province <- if_else(dfplot_Others$country=='Italy'       , 'Italy'       , dfplot_Others$province)
dfplot_Others$province <- if_else(dfplot_Others$country=='Korea, South', 'Korea, South', dfplot_Others$province)
dfplot_Others$province <- if_else(dfplot_Others$country=='Iran'        , 'Iran'        , dfplot_Others$province)

#dfplot_China <- dfplot %>% filter(country == 'China', attribute=='confirmed') 
ggplot(dfplot_China, aes(x=dates, color=province, y=value)) + 
  theme_classic() +
  geom_line(size=1) + 
  scale_y_log10() +
  ylab('Confirmed cases') +
  geom_vline(xintercept=shutdown_date_other, color='red') +
  geom_vline(xintercept=shutdown_date_Wuhan, color='blue') +
  
  geom_vline(xintercept=last_day_cases_other, color='red') +
  geom_vline(xintercept=last_day_cases_Wuhan, color='blue') +
  
  annotate(geom='text', x=shutdown_date_other+0.25, y=20000, color='red' , label='China shuts 15 more cities', hjust=0) +
  annotate(geom='text', x=shutdown_date_Wuhan+0.25, y=40000, color='blue', label='China shuts Wuhan'         , hjust=0) +

  annotate(geom='text', x=last_day_cases_other+0.25, y=20000, color='red' , label='Cases stop in rest of China', hjust=0) +
  annotate(geom='text', x=last_day_cases_Wuhan+0.25, y=40000, color='blue', label='Cases stop in Wuhan', hjust=0)

last_day_cases_other - shutdown_date_other
last_day_cases_Wuhan - shutdown_date_Wuhan 
#### Plot Hubei data, confirmed and actual
# in https://jamanetwork.com/journals/jama/fullarticle/2762130: 
# cases lag confirmed by X days. 
# the lag will tell us 
# China data starts at 01

dfhplot <- melt(dfh, 'dates')
ggplot(dfhplot %>% filter(variable %in% c('actual_daily', 'confirmed_daily')), 
       aes(x=dates, y=value, fill=variable)) +
  geom_bar(stat='identity', position='dodge') + 
  ylab('new cases per day')
  
ggplot(dfhplot %>% filter(variable %in% c('actual', 'confirmed')), 
       aes(x=dates, y=value, fill=variable)) +
  geom_bar(stat='identity', position='dodge') + 
  ylab('Cumulative cases per day')

ggplot(dfhplot %>% filter(variable %in% c('actual', 'confirmed')), 
       aes(x=dates, y=value, fill=variable)) +
  geom_bar(stat='identity', position='dodge') + 
  ylab('Cumulative cases per day')

dfhplot_logversion <- dfhplot %>% filter(variable %in% c('actual', 'confirmed'), value>0)
ggplot(dfhplot_logversion,
       aes(x=dates, y=value, color=variable)) +
  geom_line() +
  geom_point() +
  ylab('Cumulative cases') +
  scale_y_log10() +
  geom_vline(xintercept=shutdown_date_Wuhan, color='blue') +
  stat_smooth(data=dfhplot_logversion %>% filter(dates<=shutdown_date_Wuhan), method='lm') +
  annotate(geom='text', x=shutdown_date_Wuhan+0.25, y=100, color='blue', label='China shuts Wuhan' , hjust=0)








######## plot NYC & other cities:
cities = c('_Italy', 
           '_Iran' ,
           #'_Korea, South',
           '_Spain',
           '_Germany',
           'France_France',
           '_Japan', 
           'New York_US', 
           'Washington_US'
           )
cities_slow = c('_Japan',
                'Washington_US'
)
cities_fast = c('_Italy', 
                '_Iran' ,
                '_Spain',
                '_Germany',
                'France_France',
                'New York_US' 
)
cities_w_italy = c('_Italy', 'New York_US')
  
cities = cities_slow
cities = cities_fast
cities = cities_w_italy

# elim counts of <100
dfplot_log_cities <- dfplot %>% filter(Location %in% cities, attribute=='confirmed', value>100)
dfplot_log_cities <- remove_redundant_timepoints_per_Location(dfplot_log_cities, remove_col='value')
ggplot(dfplot_log_cities, 
       aes(x=dates, y=value, color=Location)) + 
  geom_line() +
  geom_point() +
  ylab('Cumulative cases') +
  scale_y_log10() +
  stat_smooth(method='lm') +
  geom_vline(xintercept=italy_school_close_date, color='blue') +
  geom_vline(xintercept=italy_lockdown_date,  color='red') +
  annotate(geom='text', x=italy_school_close_date   +0.25, y=1500, color='blue', label='Italy closes all schools', hjust=0) +
  annotate(geom='text', x=italy_lockdown_date       +0.25, y=4000, color='red' , label='Italy lockdown', hjust=0) 




  # geom_vline(xintercept=shutdown_date_Wuhan, color='blue') +
  # annotate(geom='text', x=shutdown_date_Wuhan+0.25, y=100, color='blue', label='China shuts Wuhan' , hjust=0)


## fit all cities
dfplot_log_nonChina <- dfplot %>% filter(country!='China',
                                       #  country!='Korea, South',
                                         province!='Diamond Princess',
                                         attribute=='confirmed', 
                                         value>150)
dfplot_log_nonChina <- remove_redundant_timepoints_per_Location(dfplot_log_nonChina, remove_col='value')
ggplot(dfplot_log_nonChina, 
       aes(x=dates, y=value, color=Location)) + 
  geom_line() +
  geom_point() +
  ylab('Cumulative cases') +
  scale_y_log10() +
  stat_smooth(method='lm') 


# locs_by_descending_confirmed_numbers


#############################
####  fit all the cities ####
#############################
form_conf <- log_confirmed ~ dates
form_deaths <- log_deaths ~ dates

dfcurr_log_nonChina <- dfcurr %>% filter(country!='China',
                                         province!='Diamond Princess',
                                         confirmed>150)

dfcurr_log_nonChina <- remove_redundant_timepoints_per_Location(dfcurr_log_nonChina, remove_col='confirmed')

###
df <- dfcurr_log_nonChina
form = form_conf

# to record outputs
df_fits <- data.frame()
models  <- list()

all_Locations <- unique(df$Location)
for (i in seq(1, length(all_Locations))){
  loc = all_Locations[i]
  
  # grab Location and order by date
  df_curr <- df %>% 
    dplyr::filter(Location==loc)
  df_curr <- df_curr[order(df_curr$dates),]
  
  # stop if repeated dates, which shouldn't happen. this should be fixed 
  if (length(df_curr$dates)!=length(unique(df_curr$dates))){
    stop(paste0('you have repeated dates in the dataframe for this Location:\n\n', 
                loc, 
                '\n\nThis should be fixed upon data import.'))
  }
  
  # run a model on df_curr 
  mod_curr <- glm(formula=form, data=df_curr, family = gaussian(link="identity"))
  df_out_curr <- get_params_from_loglinear_model(mod_curr) %>% as.data.frame
  df_out_curr$Location <- loc
  df_out_curr$current_log_confirmed <- max(df_curr$log_confirmed)
  df_out_curr$current_log_deaths    <- max(df_curr$log_deaths)
 
  # record output
#  models[loc] = mod_curr
  df_fits <- rbind(df_fits, df_out_curr)
}

# remove failed fits
df_fits <- df_fits %>% dplyr::filter(!is.na(r)) 
df_fits <- df_fits[order(df_fits$r),]

# plot results of doubling times
actual_tD_Hubei = 3.8
confirmed_tD_Hubei = 0.9
ggplot(df_fits, aes(x=tD, y=t0_as_date, label=Location)) + 
  geom_point() + 
#  geom_text(check_overlap=TRUE, hjust=0, nudge_x=0.05) +
  geom_text( hjust=0, nudge_x=0.05) +
  xlab('doubling time in days') +
  ylab('date of epidemic start (predicted 1st patient)') +
  geom_vline(xintercept=actual_tD_Hubei   , color='blue') +
  geom_vline(xintercept=confirmed_tD_Hubei, color='red') +
  annotate(geom='text', x=actual_tD_Hubei   +0.25, y=as.Date('2020-01-10'), color='blue', label='Hubei, Actual', hjust=0) +
  annotate(geom='text', x=confirmed_tD_Hubei+0.25, y=as.Date('2020-01-15'), color='red' , label='Hubei, Confirmed', hjust=0) 
  
  
##### time shift the cities and overlay plots!
df_shifts <- df_fits %>% dplyr::filter(tD < 4, tD > 2)

ggplot(df_shifts, aes(x=tD, y=t0_as_date, label=Location)) + 
  geom_point() + 
  #  geom_text(check_overlap=TRUE, hjust=0, nudge_x=0.05) +
  geom_text( hjust=0, nudge_x=0.05) +
  xlab('doubling time in days') +
  ylab('date of epidemic start (predicted 1st patient)') 

cities = unique(df_shifts$Location)
dfplot_log_cities_shift <- dfplot %>% filter(Location %in% cities, attribute=='confirmed', value>150)
dfplot_log_cities_shift <- remove_redundant_timepoints_per_Location(dfplot_log_cities_shift, remove_col='value')

# subtract t0's from each city!:
for (i in 1:nrow(dfplot_log_cities_shift)){
  loc <- dfplot_log_cities_shift[i,'Location']
  shift <- (df_shifts %>% filter(Location==loc))$t0
  dfplot_log_cities_shift[i,'dates_shifted'] = dfplot_log_cities_shift[i,'dates'] - shift
  
}
dfplot_log_cities_shift$dates_shifted

# get avg points to fit line
date_shifted_mean <- mean(dfplot_log_cities_shift$dates_shifted)
value_mean        <- mean(dfplot_log_cities_shift$value)

# plot all of them 
ggplot(dfplot_log_cities_shift, 
       aes(x=dates_shifted, y=value, color=Location)) + 
  geom_line() +
  geom_point() +
  ylab('Cumulative cases') +
  scale_y_log10(limits=c(1, max(dfplot_log_cities_shift$value))) +
  xlim(0,  max(dfplot_log_cities_shift$dates_shifted)) +
  stat_smooth(method='lm') 
  
# 
# saved on time_shifts tab of the analysis excel
df_shifts$time_shift_from_NewYork = df_shifts$t0 - df_shifts[df_shifts$Location=='New York_US', 't0']
write.csv(df_shifts, 'G:/CISS-West/Value Institute/Users/Matthew/p_coronavirus/dfshifts.csv')  

# geom_vline(xintercept=shutdown_date_Wuhan, color='blue') +

# exponential fit to NYP data

## death rates




















###############################
#### analyze growth curves ####
###############################

form <- log_confirmed ~ dates

#### New York curves
dfny <- dfcurr %>% filter(Location=='New York_US', confirmed>0) 
m_ny_confirmed <- glm(formula=form, data=dfny, family = gaussian(link="identity"))
outNYconf <- get_params_from_loglinear_model(m_ny_confirmed)
outNYconf





  
    
  
#### Hubei curves, actual vs. confirmed - cut off at Wuhan shutdown
form_Hubei_confirmed <- log_confirmed ~ dates
form_Hubei_actual    <- log_actual    ~ dates

dfh_fit_confirmed <- dfh %>% filter(dates <= shutdown_date_Wuhan, log_confirmed>0)
dfh_fit_actual    <- dfh %>% filter(dates <= shutdown_date_Wuhan, log_actual   >0)

m_Hubei_confirmed <- glm(formula=form_Hubei_confirmed, data=dfh_fit_confirmed, family = gaussian(link="identity"))
m_Hubei_actual    <- glm(formula=form_Hubei_actual   , data=dfh_fit_actual   , family = gaussian(link="identity"))

dfh_fit_confirmed$log_confirmed_pred <- predict(m_Hubei_confirmed)
dfh_fit_actual$log_actual_pred       <- predict(m_Hubei_actual   )

out1_confirmed <- get_params_from_loglinear_model(m_Hubei_confirmed)
out1_actual    <- get_params_from_loglinear_model(m_Hubei_actual   )


# doublecheck vs other Hubei analysis (from later, not relevant time period)
dfh_fit_confirmed_0122_to_0205 <- dfh %>% filter(dates >= as.Date('2020-01-22') & dates <= as.Date('2020-02-05'), log_confirmed>0)
m_Hubei_confirmed_0122_to_0205 <- glm(formula=form_Hubei_confirmed, data=dfh_fit_confirmed_0122_to_0205, family = gaussian(link="identity"))
out1_confirmed_0122_to_0205 <- get_params_from_loglinear_model(m_Hubei_confirmed_0122_to_0205)
out1_confirmed_0122_to_0205
# matches Hopkins curve from same time period. good. this is not the exponential phase!


#### Hopkins curves

# model what our likely 'hopkins curve' looks like
# take demonstrative cities and fit epi curves
# where is start of curve?
# what is the characteristic curve structure?

### fit hubei as example BAD EXAMPLE - 
# SEE ABOVE section " Hubei curves, actual vs. confirmed"
# THIS IS ALREADY PAST THE EXPONENTIAL PHASE
df1 <- dfcurr %>% filter(province=='Hubei') %>% filter(confirmed < exp(10))
min(df1$dates)
max(df1$dates)

df1$log_confirmed = log(df1$confirmed)
df1$log_recovered = log(df1$recovered)
df1$log_deaths    = log(df1$deaths   )

m1 <- glm(formula=log_confirmed ~ dates, data=df1, family = gaussian(link="identity"))
df1$log_confirmed_pred <- predict(m1) # , newdata=list(dates=seq(min(df$dates), max(df$dates), by=1))

# plot and model summary
ggplot(df1, aes(x=dates, y=log_confirmed)) + geom_point() + geom_line(aes(y=log_confirmed_pred), color='red')
summary(m1)
out <- get_params_from_loglinear_model(m1)
print(out)


# fit deaths
m2 <- glm(formula=log_deaths ~ dates, data=df1, family = gaussian(link="identity"))
df1$log_deaths_pred <- predict(m2) # , newdata=list(dates=seq(min(df$dates), max(df$dates), by=1))

# plot and model summary
ggplot(df1, aes(x=dates, y=log_deaths)) + geom_point() + geom_line(aes(y=log_deaths_pred), color='blue')
summary(m2)

# how determine where to cut off growth curve?

# what doubling time was reported for Hubei? how compares to 2.3 days?
# how r0 ...
# 

# t0 is the time when there was exactly 1 patient.
# r_not: how many people will 1 person infect? can't get from this..

### italy:
df1 <- dfcurr %>% filter(country=='Italy') %>% filter(confirmed < exp(10)) %>% filter(confirmed>0) #! need to filter to confirmed>0
df1$log_confirmed = log(df1$confirmed)
df1$log_recovered = log(df1$recovered)
df1$log_deaths    = log(df1$deaths   )

m1 <- glm(formula=log_confirmed ~ dates, data=df1, family = gaussian(link="identity"))
df1$log_confirmed_pred <- predict(m1) # , newdata=list(dates=seq(min(df$dates), max(df$dates), by=1))

# plot and model summary
ggplot(df1, aes(x=dates, y=log_confirmed)) + geom_point() + geom_line(aes(y=log_confirmed_pred), color='red')
summary(m1)



###############################
#### Wuhan confirmed curve ####
###############################

str(dfw)
sum(dfw$confirmed_daily)
sum((dfw %>% filter(dates<as.Date('2020-03-01')))$confirmed_daily)

ggplot(dfw, aes(x=dates, y=confirmed_daily)) + geom_bar(stat='identity')


