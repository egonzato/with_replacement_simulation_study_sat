weighted_variance_estimator=function(df=NULL,time=NULL,event=NULL,treatment=NULL,formula=NULL,cluster_id=NULL,pat_id=NULL,weight=NULL){
  # extract variables from df
  time_quo=enquo(time)
  event_quo=enquo(event)
  treatment_quo=enquo(treatment); treatment_name=as_label(treatment_quo)
  cluster_quo=enquo(cluster_id)
  pat_quo=enquo(pat_id)
  weight_quo=enquo(weight)
  # create unique id for each row
  df=df%>%mutate(unique_id=row_number())
  # build formula
  full_formula=reformulate(termlabels=attr(terms(formula),"term.labels"),
                           response = paste0("Surv(", as_label(time_quo), ",", as_label(event_quo), ")"))
  # extract clusters for models
  cluster_vec=df %>% pull(!!cluster_quo)
  pat_vec=df %>% pull(!!pat_quo)
  unique_vec=df$unique_id
  # function if weights argument is null
  if (rlang::quo_is_null(weight_quo)) {
    model_cluster=coxph(full_formula,df,cluster=cluster_vec,robust=TRUE)
    model_id=coxph(full_formula,df,cluster=pat_vec,robust=TRUE)
    model_unique=coxph(full_formula,df,cluster=unique_vec,robust=TRUE)
  } else {
    weight_quo=enquo(weight)
    weights_vec=df %>% pull(!!weight_quo)
    model_cluster=coxph(full_formula,df,cluster=cluster_vec,robust=TRUE,weights=weights_vec)
    model_id=coxph(full_formula,df,cluster=pat_vec,robust=TRUE,weights=weights_vec)
    model_unique=coxph(full_formula,df,cluster=unique_vec,robust=TRUE,weights=weights_vec)
  }
  # get variance from each model
  V1=vcov(model_cluster)
  V2=vcov(model_id)
  V3=vcov(model_unique)
  # calculate k1, k2, k3
  k1=length(unique(cluster_vec))
  k2=length(unique(pat_vec))
  k3=length(unique(unique_vec))
  # calculate variance
  var=((k1/(k1-1))*V1)+((k2/(k2-1))*V2)+((k3/(k3-1))*V3)
  # extract standard error
  ses=sqrt(diag(var))
  # get coefficient for beta
  betas=model_cluster$coefficients
  # calculate hr with lower and upper bound of confidence interval
  hr=exp(betas[[treatment_name]])
  lower=exp(betas[[treatment_name]]-1.96*ses[[treatment_name]]) 
  upper=exp(betas[[treatment_name]]+1.96*ses[[treatment_name]])
  # calculate z statistic and p-value
  z=betas[[treatment_name]]/ses[[treatment_name]]
  pval=2*(1-pnorm(abs(z)))
  # return
  return(list(dataset=df,
              models=list(cluster=model_cluster,
                          patient=model_id,
                          unique=model_unique),
              treatment=list(HR=hr,
                             lower=lower,
                             upper=upper),
              coefficients=betas,
              variance=var,
              stderr=ses,
              z=z,
              pvalue=pval))
}