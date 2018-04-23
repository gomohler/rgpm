library(sp)
library(maptools)
library(dplyr)
library(tidyr)
library(randomForest)
library(glmnet)


# make_grid()
#----------------------------------------------------------------------------
# make 2D rectangular grids
# inputs:
#   h: grid height
#   theta: rotation angle
#   offset: offset of grid
#   xrng, yrng: range for grid
#   area: area of a grid cell
# outputs:
#   SpatialPolygons object of grid cells
# notes:
#   requires(sp); requires(maptools)
#----------------------------------------------------------------------------
make_grid<-function(h, theta, offset, xrng, yrng, area=62500){
  fx=0
  fy=0
  w = area/h
  centerx=mean(xrng)
  centery=mean(yrng)
  Mx=ceiling((xrng[2]-xrng[1])/h)
  My=ceiling((yrng[2]-yrng[1])/w)
  Dx=ceiling(Mx*fx)*h
  Dy=ceiling(My*fy)*w
  grd <- GridTopology(cellcentre.offset=c(xrng[1]+offset-Dx, yrng[1]-Dy), 
                      cellsize=c(h,w), cells.dim=c(Mx+2*ceiling(Mx*fx),My+2*ceiling(My*fy)))
  grdp <- as.SpatialPolygons.GridTopology(grd)
  grid_rotate=elide(grdp, rotate = theta, center=c(centerx,centery))
return(grid_rotate)
}


# make_rf_features()
#----------------------------------------------------------------------------
# makes predictive features (predictors) from event history for RF model
# 1. # events [0-2] months before forecasting window
# 2. # events [3-5] months before forecasting window
# 3. # events [6-14] months before forecasting window
# 4. # events >14 months before forecasting window
# 5. # events in same 3 month period as forecasting window (e.g. March-May)
# - for crimetypes: burglary, mvt, street, other
# inputs:
#   events: event data with date, type, g_ind (grid index)
#   pred_date: date for making predictions. Only events<pred_date used.
# outputs:
#   data frame of predictors for every grid that has at least one event.
#   cols: g_ind, count1_burglary, ..., count5_other  
#----------------------------------------------------------------------------
make_rf_features <- function(events, pred_date){
  pred_date = as.Date(pred_date) 
  target_month = as.integer(format(pred_date, '%m'))
  dim = 30.4      # number of days in month
  e = events %>% 
    filter(!is.na(g_ind), date < pred_date) %>%
    mutate(d = as.numeric(pred_date - date)) %>%  # number of days before start
    mutate(period = case_when(
      d <= dim*2 ~ '1', 
      d > dim*2 & d <= dim*5 ~ '2',
      d > dim*5 & d <= dim*14 ~ '3',
      d > dim*14 ~ '4' 
    ))
  tmp1 = e %>% count(g_ind, period, type)  # period 1-4 counts
  tmp2 = e %>%                             # get period 5 counts
    mutate(month = as.integer(format(date, '%m'))) %>% 
    filter(month %in% (target_month + (0:2))) %>% 
    count(g_ind, type) %>% 
    mutate(period='5')
  tmp3 = bind_rows(tmp1, tmp2) %>% 
    group_by(period, type) %>% 
    mutate(x = ifelse(period %in% c(4,5), n/max(n), n)) %>% # rescale big counts
    ungroup  
  X =  tmp3 %>% 
    mutate(name = paste0('count', period, '_', type)) %>% 
    select(g_ind, name, x) %>% 
    spread(name, x, fill=0)
return(X)
}

# fit_grid_rf()
#----------------------------------------------------------------------------
# fits random forest model for hotspot modeling 
# inputs:
#   events: event data with x, y, date, type, g_ind (grid index)
#   crimetype: "burglary", "mvt", "street", "all"  
#   pred_start: start of prediction period
#   duration: length of training period
#   ntree: number of trees for random forest
# outputs:
#   pred: scores from function for [pred_start, pred_start+duration]
#   VI: variable importance 
#----------------------------------------------------------------------------
fit_grid_rf <- function(events, crimetype, pred_start = "2017-03-01", duration = 91,
                        ntree = 2000){
  if(crimetype[1] == "all") crimetype = c("burglary", "mvt", "street", "other")
  grid_only = na.omit(unique(events$g_ind)) # cells with at least one event
  pred_start = as.Date(pred_start)  
  train_start = pred_start - duration
  train_end = pred_start - 1 
  y = events %>% filter(!is.na(g_ind)) %>% 
    filter(between(date, train_start, train_end), type %in% crimetype) %>% 
    count(g_ind) %>% 
    mutate(y = log(n + 1)) %>% 
    tidyr::complete(g_ind = grid_only, fill=list(n=0, y=log(1)))
  x = make_rf_features(events, train_start)
  data = full_join(y, x, by="g_ind")   # expand x so it matches y
  data[is.na(data)] = 0    
  xmat = select(data, -g_ind, -n, -y)
  set.seed(1)
  rf = randomForest::randomForest(x=xmat, y=data$y, ntree=ntree)
  xpred = make_rf_features(events, pred_start) 
  pred = predict(rf, 
                 newdata = xpred %>% select(-g_ind)) %>% 
    as_tibble %>% rename(score=value) %>% 
    add_column(g_ind=xpred$g_ind, .before=1) %>% 
    arrange(-score)
  list(pred=pred, VI=rf$importance)
}  


# make_hawkes_features()
#----------------------------------------------------------------------------
# makes predictive 'hawkes' features (predictors) from event history in a grid
#  these use (shifted) geometric decay 
# inputs:
#   events: event data with date, type, g_ind (grid index)
#   pred_date: date for making predictions. Only events<pred_date used.
#   omega: vector of decay parameters
# outputs:
#   data frame of predictors for every grid that has at least one event.
#   cols: g_ind, burglary_w1, ...,other_w1,...,total_wp  
#----------------------------------------------------------------------------
make_hawkes_features <- function(events, pred_date, omega=c(.001, .005, .01, .02, .1)){
  pred_date = as.Date(pred_date) 
  tmp1 = events %>% select(g_ind, type, date) %>% 
    filter(!is.na(g_ind), date < pred_date) %>% 
    mutate(d = as.numeric(pred_date - date)) %>% 
    select(-date)
  dmax = min(max(tmp1$d), 365*3)  
  d = 1:dmax
  Decay = lapply(omega, function(p) dgeom(d-1, prob=p))  
  names(Decay) = omega
  Decay = as.tibble(Decay) %>% add_column(d, .before=1)
  tmp2 = tmp1 %>% filter(d <= dmax) %>% 
    left_join(Decay, by='d') %>%    
    select(-d) %>% 
    group_by(g_ind, type) %>% summarize_all(funs(sum)) 
  X = tmp2 %>% gather(parameter, weight, -g_ind, -type) %>% 
    unite(decay, type, parameter, sep='_') %>% 
    spread(decay, weight, fill=0L) %>% ungroup
  TYPES = unique(events$type)
  for(i in seq_along(omega)){
    p = omega[i]
    cnames = paste(TYPES, p, sep='_')
    new_cname = paste('total', p, sep='_')
    X[new_cname] = rowSums(X[cnames])
  }
  return(X)
}



# fit_grid_glmnet()
#----------------------------------------------------------------------------
# fits glmnet for hotspot modeling 
# inputs:
#   events: event data with x, y, date, type, g_ind (grid index)
#   crimetype: "burglary", "mvt", "street", "all"
#   pred_start: start of prediction period
#   duration: length of training period
#   thres: threshold for determining the event count that makes a grid a hotspot
#   omega: vector of decay parameters
#   K: number of cross-validation folds
#   alpha: alpha parameter for glmnet
#   nperiods: number of periods to use for training
#   nkeep: number of cold grids to keep for model fitting
# outputs:
#   pred: scores from function for [pred_start, pred_start+duration]
#   VI: variable importance
#----------------------------------------------------------------------------
fit_grid_glmnet <- function(events, crimetype, pred_start = "2017-03-01", 
                            duration = 91, thres=1, omega=c(.001, .005, .01, .02, .1), 
                            K = 10, alpha = 0.8, nperiods = 1, nkeep = 60000){
  if(crimetype[1] == "all") crimetype = c("burglary", "mvt", "street", "other")
  grid_only = na.omit(unique(events$g_ind)) # cells with at least one event
  pred_start = as.Date(pred_start)   
  data = tibble()
  for(k in 1:nperiods){
    train_start = pred_start - duration*k
    train_end = pred_start - 1 - duration*(k-1)
    y = events %>% filter(!is.na(g_ind)) %>% 
      filter(between(date, train_start, train_end), type %in% crimetype) %>% 
      count(g_ind) %>% 
      mutate(y = 1L*(n>=thres)) %>% 
      tidyr::complete(g_ind = grid_only, fill=list(n=0, y=0))
    x = make_hawkes_features(events, train_start, omega)
    dd = full_join(y, x, by="g_ind")   # expand x so it matches y
    dd[is.na(dd)] = 0    
    data = bind_rows(data, dd)
  }
  set.seed(112233)                      
  hot = which(data$y == 1)
  n.hot = length(hot)
  cold = which(data$y != 1)
  if(length(cold) > nkeep) cold = sample(cold, size=nkeep)
  data = data[c(hot, cold), ]
  set.seed(112233)           
  fold_id = data %>% select(y) %>% mutate(ind=row_number()) %>% 
    group_by(y) %>% mutate(k = sample(rep(seq(K),length=n()))) %>% 
    arrange(ind) %>% pull(k)
  xmat = select(data, -g_ind, -n, -y) %>% as.matrix() 
  m.cv = glmnet::cv.glmnet(x=xmat, y=data$y, family="binomial",alpha=alpha, foldid=fold_id)
  lam.opt = m.cv$lambda.min  
  coef = coef(m.cv, s=lam.opt)
  std = apply(xmat, 2, sd)
  VI = tibble(predictor = colnames(xmat), coef= coef[-1,1], sd = std, imp = coef*sd,
              imp2 = imp/sum(imp))
  xpred = make_hawkes_features(events, pred_start, omega) 
  pred = predict(m.cv, newx = xpred %>% select(-g_ind) %>% as.matrix, 
                 s = lam.opt) %>% 
    as_tibble %>% rename(score=`1`) %>% 
    add_column(g_ind=xpred$g_ind, .before=1) %>% 
    arrange(-score)
list(pred=pred, VI=VI)
}  